#!/usr/bin/env python3
"""Checks of the ISIMIP4a local-lake files before upload.

    python3 local_lakes/qc_isimip.py <isimip dir> [out.csv]

Per lake folder (<isimip dir>/<Lake>/): 13 files with the expected names; 45,656 days
(1901-01-01..2025-12-31, days since 1901-01-01, step 1); required global attributes; SHA256SUMS
matching; values finite except where allowed (thermodepth without a thermocline); plausible ranges;
levlak increasing from 0 to the lake depth with contiguous bounds. One row per lake, problems listed.
"""
import hashlib
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
VARS = ["watertemp", "surftemp", "bottemp", "ice", "icethick", "strat", "thermodepth", "sensheatf",
        "latentheatf", "swup", "lwup", "lakeheatf", "momf"]
RANGE = {"watertemp": (270, 315), "surftemp": (270, 315), "bottemp": (270, 315), "ice": (0, 1), "icethick": (0, 5),
         "strat": (0, 1), "thermodepth": (0, 600), "sensheatf": (-500, 500), "latentheatf": (-200, 800),
         "swup": (0, 400), "lwup": (100, 600), "lakeheatf": (-1500, 1500), "momf": (0, 5)}
NAN_OK = {"thermodepth"}
ATTRS = ["contact", "institution", "author", "references", "acknowledgement", "license", "model",
         "simulation_round", "sector", "lake_site", "comment"]
ND = 45656


def check(folder, site):
    prob, info = [], {}
    sums = {}
    if (folder / "SHA256SUMS").exists():
        for line in (folder / "SHA256SUMS").read_text().splitlines():
            h, n = line.split(None, 1); sums[n.strip()] = h
    else:
        prob.append("no SHA256SUMS")
    for v in VARS:
        f = folder / f"glm4_gswp3-era5_obsclim_histsoc_default_{v}_{site}_daily_1901_2025.nc"
        if not f.exists():
            prob.append(f"missing {v}"); continue
        if sums and hashlib.sha256(f.read_bytes()).hexdigest() != sums.get(f.name):
            prob.append(f"{v}: checksum mismatch")
        d = nc.Dataset(f); d.set_auto_mask(False)
        t = d["time"][:]
        if len(t) != ND or t[0] != 0 or np.any(np.diff(t) != 1):
            prob.append(f"{v}: time axis ({len(t)} steps)")
        if d["time"].units != "days since 1901-01-01 00:00:00" or d["time"].calendar != "proleptic_gregorian":
            prob.append(f"{v}: time units/calendar")
        miss = [a for a in ATTRS if not str(getattr(d, a, "")).strip()]
        if miss:
            prob.append(f"{v}: attributes missing {miss}")
        x = d[v][:].astype("f8"); x[x > 1e19] = np.nan
        if v == "watertemp":
            # allowed: empty levels below the lake bottom on days the level is below the crest
            # (ice, evaporation), i.e. missing values only as one run at the bottom of a profile
            m = np.isnan(x)
            tail_ok = np.all(m[:, :-1] <= m[:, 1:]) and not m[:, 0].any()
            if not tail_ok:
                prob.append("watertemp: missing values above valid ones")
            info["days_bottom_level_empty_pct"] = round(100 * float(m.any(1).mean()), 1)
            info["max_empty_levels"] = int(m.sum(1).max())
        elif v not in NAN_OK and np.isnan(x).any():
            prob.append(f"{v}: {int(np.isnan(x).sum())} missing values")
        lo, hi = RANGE[v]
        if np.nanmin(x) < lo or np.nanmax(x) > hi:
            prob.append(f"{v}: range {np.nanmin(x):.4g}..{np.nanmax(x):.4g} outside {lo}..{hi}")
        if v == "watertemp":
            z, b = d["levlak"][:], d["levlak_bnds"][:]
            if z[0] != 0 or np.any(np.diff(z) <= 0) or np.any(b[1:, 0] != b[:-1, 1]) or b[0, 0] != 0:
                prob.append("levlak/levlak_bnds")
            info.update(levels=len(z), depth_m=float(b[-1, 1]))
        if v == "surftemp":
            info["surftemp_mean_C"] = round(float(np.nanmean(x)) - 273.15, 2)
        if v == "ice":
            info["ice_days_per_year"] = round(float(np.nansum(x)) / 125, 1)
        if v == "strat":
            info["strat_days_per_year"] = round(float(np.nansum(x)) / 125, 1)
            if "summer" not in d["strat"].comment:
                prob.append("strat: not the summer definition")
        if v == "lakeheatf":
            info["lakeheatf_mean_Wm2"] = round(float(np.nanmean(x)), 2)
        d.close()
    return prob, info


def main(isimip, out=None):
    sites = pd.read_csv(ROOT / "local_lakes" / "lake_sites.csv").set_index("lake")
    rows = []
    for lake in sites.index:
        folder = Path(isimip) / lake
        if not folder.exists():
            rows.append({"lake": lake, "ok": False, "problems": "no folder"}); continue
        prob, info = check(folder, sites.loc[lake, "lake_site"])
        rows.append({"lake": lake, "ok": not prob, "problems": "; ".join(prob), **info})
    r = pd.DataFrame(rows)
    if out:
        r.to_csv(out, index=False)
    bad = r[~r.ok]
    print(f"{r.ok.sum()} of {len(r)} lakes pass")
    for _, x in bad.iterrows():
        print(f"  {x.lake}: {x.problems}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
