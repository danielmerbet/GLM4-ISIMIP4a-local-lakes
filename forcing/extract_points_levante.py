#!/usr/bin/env python3
"""Extract ISIMIP4a GSWP3-ERA5 daily forcing at lake points -- on Levante (default --src) or on
Hydra from the group mirror (--src). Originally written for the local lakes ON LEVANTE,
so only small per-lake CSVs have to be transferred (the global files are ~301 GB).

    module load python3/2025.01-gcc-13.3.0
    python3 extract_points_levante.py --lakes lakes_points.csv --out out [--workers 6]

Adapted from ISIMIP4_LocalLakes_LERWQ/scripts/extract_isimip4a_points.py (same cell
choice and read strategy); changes: standalone (no lerwq config), ISIMIP4 file names
(..._global_30arcmin_daily_<y0>_<y1>.nc), all 11 variables, lake list from a CSV.

Cell choice: nearest 0.5 deg cell centre; if that cell is masked (GSWP3-ERA5 is
land-only, so large-lake/coastal cells can be missing), the nearest valid cell within
2 cells is used. points.csv records the cell used per lake and whether it was a fallback.

Read strategy: the files store one full 360x720 field per day, so reading one point
costs as much as reading the whole field. Each file is therefore read once, in yearly
blocks, and every lake is picked from each block. One process per variable.

Output, per lake: gswp3-era5_obsclim_<lake>_daily_1901_2025.csv with columns
time,tas,tasmin,tasmax,hurs,huss,pr,prsn,rsds,rlds,ps,sfcwind -- the same layout and
units (K, kg m-2 s-1, Pa, W m-2, m s-1, %, kg kg-1) as the ISIMIP3 local-lake files.
"""
import argparse
import glob
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd

SRC = "/work/bb0820/ISIMIP/ISIMIP4a/InputData/climate/atmosphere/obsclim/global/daily/historical/GSWP3-ERA5"
VARS = ["tas", "tasmin", "tasmax", "hurs", "huss", "pr", "prsn", "rsds", "rlds", "ps", "sfcwind"]
NFILES = 13   # 1901-1910 ... 2011-2020, 2021-2025


def choose_cells(sample_file, lakes):
    d = nc.Dataset(sample_file)
    lat, lon = d["lat"][:], d["lon"][:]
    v = [k for k in d.variables if k not in ("time", "lat", "lon")][0]
    field = d[v][0]
    mask = np.ma.getmaskarray(field) | ~np.isfinite(np.ma.filled(field, np.nan))
    d.close()
    rows = []
    for _, r in lakes.iterrows():
        i0, j0 = int(np.abs(lat - r.lat).argmin()), int(np.abs(lon - r.lon).argmin())
        best = (i0, j0, False)
        if mask[i0, j0]:
            cands = [(i, j % len(lon)) for i in range(i0 - 2, i0 + 3) for j in range(j0 - 2, j0 + 3)
                     if 0 <= i < len(lat) and not mask[i, j % len(lon)]]
            if not cands:
                raise SystemExit(f"{r.lake}: no valid land cell within 2 cells of ({r.lat}, {r.lon})")
            i, j = min(cands, key=lambda c: (lat[c[0]] - r.lat) ** 2 + (lon[c[1]] - r.lon) ** 2)
            best = (i, j, True)
        i, j, fb = best
        rows.append({"lake": r.lake, "lat": r.lat, "lon": r.lon, "ilat": i, "ilon": j,
                     "cell_lat": float(lat[i]), "cell_lon": float(lon[j]), "fallback_cell": fb})
    return pd.DataFrame(rows)


def extract_var(args):
    var, files, ilat, ilon = args
    out = []
    for f in files:
        d = nc.Dataset(f)
        t = nc.num2date(d["time"][:], d["time"].units, calendar=getattr(d["time"], "calendar", "standard"),
                        only_use_cftime_datetimes=False, only_use_python_datetimes=True)
        x = d[var]
        vals = np.empty((len(t), len(ilat)), dtype="f8")
        for s in range(0, len(t), 366):
            block = x[s:s + 366]                      # (days, 360, 720): whole fields
            vals[s:s + len(block)] = np.ma.filled(block[:, ilat, ilon], np.nan)
        d.close()
        out.append(pd.DataFrame(vals, index=pd.DatetimeIndex(t).normalize()))
        print(f"  {var} {Path(f).name}", flush=True)
    return var, pd.concat(out).sort_index()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lakes", required=True, help="CSV with columns lake,lat,lon")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--src", default=SRC, help="directory with the global GSWP3-ERA5 files "
                    "(default: DKRZ Levante; on Hydra: the group mirror, see forcing/README.md)")
    a = ap.parse_args()
    lakes = pd.read_csv(a.lakes)
    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)

    files = {v: sorted(glob.glob(f"{a.src}/gswp3-era5_obsclim_{v}_global_30arcmin_daily_*.nc")) for v in VARS}
    for v, fs in files.items():
        if len(fs) != NFILES:
            raise SystemExit(f"{v}: expected {NFILES} decade files in {a.src}, found {len(fs)}")
    cells = choose_cells(files["tas"][0], lakes)
    cells.to_csv(out / "points.csv", index=False)
    print(cells.to_string(index=False), flush=True)

    ilat, ilon = cells.ilat.values, cells.ilon.values
    with ProcessPoolExecutor(a.workers) as ex:
        res = dict(ex.map(extract_var, [(v, files[v], ilat, ilon) for v in VARS]))

    problems = []
    for k, lake in enumerate(cells.lake):
        df = pd.DataFrame({v: res[v][k] for v in VARS})
        df.index.name = "time"
        bad = df.isna().sum()
        if bad.any():
            problems.append(f"{lake}: missing values {bad[bad > 0].to_dict()}")
        y0, y1 = df.index.min().year, df.index.max().year
        df.to_csv(out / f"gswp3-era5_obsclim_{lake.lower()}_daily_{y0}_{y1}.csv", float_format="%.6g")
    print(f"-> {out} ({len(cells)} lakes)")
    if problems:
        print("PROBLEMS:\n" + "\n".join(problems))
        raise SystemExit(1)
    print("DONE: no missing values")


if __name__ == "__main__":
    main()
