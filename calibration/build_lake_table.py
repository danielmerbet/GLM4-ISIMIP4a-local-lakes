#!/usr/bin/env python3
"""Build the lake list for the GLM4 hydrodynamic calibration, with harmonised
temperature observations and hypsographs.

    python3 calibration/build_lake_table.py

Lakes = the 73 ISIMIP3 local lakes + the lakes new in the ISIMIP4 delivery
(Edward, Grossglienickersee, InreHarrsjon, MellerstaHarrsjon, Villasjon).
Carioca appears only in the ISIMIP4 water-quality tab and has no data files
here, so it is not included.

Sources, in order of preference per lake:
  * temperature + hypsometry: the ISIMIP4 delivery when the lake has one (via
    the ISIMIP4_LocalLakes_LERWQ harmoniser and its per-file catalog, which
    handles the mixed date formats/depth signs), else the ISIMIP3 database.
  * metadata (lat, lon, elevation, depth, area, Secchi, Kw): the ISIMIP4
    local-lakes Google Sheet, "MetaData Hydrothermal Corrected" tab for the 73
    ISIMIP3 lakes and "MetaData Hydrothermal" for the new ones (downloaded to
    reference/isimip4_locallakes_sheet/). Coordinates in ISIMIP4_LocalLakes_LERWQ
    config/lakes.csv override (Bosumtwi's protocol coordinates are ~300 km off).

Writes calibration/lakes.csv, calibration/obs/<lake>_temp.csv
(datetime,depth,value) and calibration/hypso/<lake>.csv (depth,area).
"""
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
LERWQ = Path(os.environ.get("LERWQ_DIR", "../ISIMIP4_LocalLakes_LERWQ"))        # LakeEnsemblR.WQ harmoniser
I3 = Path(os.environ.get("ISIMIP3_LOCAL_LAKES", "data/ISIMIP3_LocalLakes"))   # ISIMIP3 local-lake observations
I4 = Path(os.environ.get("ISIMIP4_LOCAL_LAKES", "data/ISIMIP4_LocalLakes"))   # ISIMIP4 local-lake delivery
SHEET = HERE.parent / "reference/isimip4_locallakes_sheet"
NEW = ["Edward", "Grossglienickersee", "InreHarrsjon", "MellerstaHarrsjon", "Villasjon"]
# Harp: the LERWQ catalog uses only the ISIMIP3 file (147 profile days, 2004-2014); the
# ISIMIP4 delivery adds 1,599 daily logger days (2011-2015). Combined for the calibration.
EXTRA_TEMP = {"Harp": [(I4 / "Harp" / "harp_watertemp_daily.csv", "WATERTEMP")]}

sys.path.insert(0, str(LERWQ))
from lerwq.config import load_config, load_lakes  # noqa: E402
from lerwq.obs import harmonize_lake, LIMITS      # noqa: E402

TMIN, TMAX = LIMITS["Temp_degreeCelcius"]

# Lake-specific removals of values that are physically impossible for that lake (checked by hand).
# Bosumtwi (ISIMIP4 delivery, bosumtwi_watertemp_monthly.csv): 10 values where the temperature
# equals the depth (0.20 at 0.2 m, 7.11 at 7 m, ..., 36.68 at 36 m) -- the depth was entered in
# the temperature column; observed range otherwise 26.0-32.7 C (deep water 27.2 C). Reported for the data contributor.
OBS_QC = {"Bosumtwi": {"min_value": 24.0, "max_value": 33.0,
                        "reason": "value = depth entries (column mix-up), outside 24-33 C"}}
# Hypsometry area factors (user decision 2026-10-23). Trout Lake's ISIMIP3 hypsometry gives 160,790,000 m2 at
# the surface, ten times the lake (NTL-LTER: 1,607.9 ha; sheet 15.65 km2): a decimal shift. All its areas / 10.
HYPSO_AREA_FACTOR = {"Trout": (0.1, "areas ten times too large in the source file (decimal shift); divided by 10")}


def num(x):
    """Sheet numbers use decimal commas in some tabs."""
    return pd.to_numeric(pd.Series(x).astype(str).str.replace(",", ".", regex=False), errors="coerce")


def isimip3_temp(lake):
    f = I3 / lake / f"{lake}_temp_daily.csv"
    d = pd.read_csv(f)
    d.columns = [c.strip() for c in d.columns]
    t = pd.to_datetime(d.TIMESTAMP.astype(str).str.replace(r"\.0$", "", regex=True), format="%Y%m%d", errors="coerce")
    o = pd.DataFrame({"datetime": t, "depth": pd.to_numeric(d.DEPTH, errors="coerce").abs().round(1),
                      "value": pd.to_numeric(d.WTEMP, errors="coerce")}).dropna()
    o = o[(o.value >= TMIN) & (o.value <= TMAX)]
    return o.groupby(["datetime", "depth"], as_index=False).value.mean(), f.name


def read_hypso(f):
    d = pd.read_csv(f)
    d.columns = [c.strip() for c in d.columns]      # Arendsee/Scharmutzel headers end in a space
    h = pd.DataFrame({"depth": pd.to_numeric(d.DEPTH, errors="coerce").abs(),
                      "area": pd.to_numeric(d.BATHYMETRY_AREA, errors="coerce")}).dropna()
    return h.groupby("depth", as_index=False).area.max().sort_values("depth")


def main():
    cfg = load_config()
    l4 = load_lakes(cfg).set_index("lake")
    meta3 = pd.read_csv(SHEET / "metadata_hydrothermal_corrected.csv")
    meta_all = pd.read_csv(SHEET / "metadata_hydrothermal.csv")
    meta_all = meta_all[meta_all["Lake Short Name"].notna()]
    norm = lambda s: s.str.normalize("NFKD").str.encode("ascii", "ignore").str.decode("ascii").str.replace(" ", "").str.lower()
    meta_all["key"] = norm(meta_all["Lake Short Name"])

    # "Zurich" and "LowerLakeZurich" are the same lake in the ISIMIP3 database
    # (byte-identical temperature and hypsometry files) -- keep one, or it counts
    # twice in the cross-lake median.
    lakes = [l for l in meta3["Lake.Short.Name"] if l != "Zurich"] + NEW
    (HERE / "obs").mkdir(exist_ok=True)
    (HERE / "hypso").mkdir(exist_ok=True)
    rows = []
    for lake in lakes:
        r = {"lake": lake}
        # --- metadata
        m3 = meta3[meta3["Lake.Short.Name"] == lake]
        if len(m3):
            m = m3.iloc[0]
            r.update(lat=m["latitude.dec.deg"], lon=m["longitude.dec.deg"], elev_m=m["elevation.m"],
                     mean_depth_m=m["mean.depth.m"], max_depth_m=m["max.depth.m"], area_km2=m["lake.area.sqkm"],
                     secchi_m=m["Average.Secchi.disk.depth.m"], kw_measured=m["Light.extinction.coefficient.m"],
                     meta_source="sheet:corrected")
        else:
            ma = meta_all[meta_all.key == lake.lower()]
            if not len(ma):
                ma = meta_all[meta_all.key.str.startswith(lake.lower()[:8])]
            m = ma.iloc[0]
            r.update(lat=num(m["Latitude (dec deg)"])[0], lon=num(m["Longitude (dec deg)"])[0],
                     elev_m=num(m["Elevation (m)"])[0], mean_depth_m=num(m["Mean depth (m)"])[0],
                     max_depth_m=num(m["Max depth (m)"])[0], area_km2=num(m["lake area (km²)"])[0],
                     secchi_m=num(m["Average Secchi disk depth [m]"])[0],
                     kw_measured=num(m["Light extinction coefficient [m-1]"])[0],
                     public=str(m.iloc[-1]), meta_source="sheet:hydrothermal")
        # LERWQ coordinates override (Bosumtwi correction, approximate Villasjon)
        if lake in l4.index:
            r.update(lat=l4.loc[lake, "lat"], lon=l4.loc[lake, "lon"], elev_m=l4.loc[lake, "elev_m"],
                     coords_note=str(l4.loc[lake, "notes"]))
        # --- temperature + hypsometry
        in4 = lake in l4.index
        if in4:
            obs, _ = harmonize_lake(cfg, lake)
            t = obs[obs.variable_global_name == "Temp_degreeCelcius"]
            from lerwq.obs import load_obs_catalog
            from lerwq.config import resolve_source
            tfiles = [str(resolve_source(cfg, e["file"])) for e in load_obs_catalog().get(lake, [])
                      if e["var"] == "Temp_degreeCelcius"]
            t = t[["datetime", "depth", "value"]]
            if t.empty:
                t, f = isimip3_temp(lake); tfiles = [str(I3 / lake / f)]
            hs = l4.loc[lake, "hypso_source"]
            prefix = l4.loc[lake, "prefix"]
            hf = I4 / lake / f"{prefix}_hypsometry.csv" if hs == "wq" else I3 / lake / f"{lake}_hypsometry.csv"
        else:
            t, f = isimip3_temp(lake); tfiles = [str(I3 / lake / f)]
            hf = I3 / lake / f"{lake}_hypsometry.csv"
        # Extra temperature files added on top of the catalog/ISIMIP3 choice (user decision,
        # 2026-10-21): averaged with the base data per day and 0.1 m depth.
        for f, col in EXTRA_TEMP.get(lake, []):
            d = pd.read_csv(f)
            d.columns = [c.strip() for c in d.columns]
            x = pd.DataFrame({"datetime": pd.to_datetime(d.TIMESTAMP.astype(str).str[:10], format="%Y-%m-%d", errors="coerce"),
                              "depth": pd.to_numeric(d.DEPTH, errors="coerce").abs().round(1),
                              "value": pd.to_numeric(d[col], errors="coerce")}).dropna()
            x = x[(x.value >= TMIN) & (x.value <= TMAX)]
            t = (pd.concat([t, x]).groupby(["datetime", "depth"], as_index=False).value.mean())
            tfiles.append(str(f))
            print(f"  {lake}: added {x.datetime.nunique()} days from {Path(f).name}")
        if not hf.exists():                           # some ISIMIP3 lakes keep it in ../hypso
            hf = I3.parent / "hypso" / f"{lake}_hypsometry.csv"
        h = read_hypso(hf)
        if lake in HYPSO_AREA_FACTOR:
            fac, why = HYPSO_AREA_FACTOR[lake]
            h = h.assign(area=h.area * fac)
            print(f"  {lake}: hypsometry areas x {fac} ({why})")
        # A fully inverted depth axis (area growing with depth at every level) is a height-
        # above-bottom column mislabelled DEPTH: Bosumtwi's ISIMIP4 file (0.06 km2 at "0 m",
        # 48.8 km2 at 73.4 m; the ISIMIP3 file has the same numbers the right way up).
        hyp_note = ""
        if len(h) > 2 and (h.area.diff().dropna() > 0).all():
            h = h.assign(depth=(h.depth.max() - h.depth).round(3)).sort_values("depth")
            hyp_note = "depth axis was inverted in source file; flipped"
            print(f"  {lake}: {hyp_note}")
        elif (h.area.diff().dropna() > 0).any():
            hyp_note = "area increases with depth at some level"
            print(f"  WARNING {lake}: {hyp_note}")
        qc = OBS_QC.get(lake)
        n_qc = 0
        if qc:
            keep = (t.value >= qc["min_value"]) & (t.value <= qc["max_value"])
            n_qc = int((~keep).sum()); t = t[keep]
            print(f"  {lake}: removed {n_qc} values ({qc['reason']})")
        t = t.sort_values(["datetime", "depth"])
        t.assign(datetime=t.datetime.dt.strftime("%Y-%m-%d")).to_csv(HERE / "obs" / f"{lake}_temp.csv", index=False)
        h.to_csv(HERE / "hypso" / f"{lake}.csv", index=False)
        r.update(obs_qc_removed=n_qc,
                 temp_files=";".join(tfiles), hypso_file=str(hf), hypso_note=hyp_note, hypso_max_depth=h.depth.max(), hypso_levels=len(h),
                 obs_first=t.datetime.min().date(), obs_last=t.datetime.max().date(),
                 obs_days=t.datetime.nunique(), obs_values=len(t), obs_depths=t.depth.nunique())
        rows.append(r)
        print(f"{lake:20s} obs {r['obs_first']}..{r['obs_last']} days={r['obs_days']:6d} "
              f"depths={r['obs_depths']:4d} hypso_levels={len(h)}")
    out = pd.DataFrame(rows)
    out.to_csv(HERE / "lakes.csv", index=False)
    print(f"\n{len(out)} lakes -> {HERE/'lakes.csv'}")
    print("obs span overall:", out.obs_first.min(), "->", out.obs_last.max())


if __name__ == "__main__":
    main()
