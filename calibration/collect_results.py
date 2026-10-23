#!/usr/bin/env python3
"""Collect per-lake GLM4 calibration results and compute the cross-lake medians used for
the global GLM4-AED3 runs.

    python3 calibration/collect_results.py [--exclude Taihu,...] [--min-days 20]

Reads calibration/results/<Lake>/best.csv (calibrate_lake.R). Writes
  calibration/results_all.csv        one row per lake, with lake metadata and flags
  calibration/median_parameters.csv  median, IQR and n per parameter, for the lakes kept
Lakes are dropped from the median only when explicitly excluded or when they have fewer
than --min-days observation days in the calibration half; every dropped lake is listed
with its reason. Nothing is dropped silently.
"""
import argparse
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
PARAMS = ["wind_factor", "sw_factor", "kw_scaling", "coef_mix_hyp", "coef_mix_conv", "coef_mix_turb"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exclude", default="", help="comma-separated lakes to leave out of the median")
    ap.add_argument("--min-days", type=int, default=0, help="min observation days in calibration half")
    a = ap.parse_args()
    excl = {x.strip() for x in a.exclude.split(",") if x.strip()}

    files = sorted((HERE / "results").glob("*/best.csv"))
    res = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    lakes = pd.read_csv(HERE / "lakes.csv")[["lake", "lat", "lon", "elev_m", "max_depth_m", "area_km2", "hypso_note"]]
    win = pd.read_csv(HERE / "windows.csv")[["lake", "days_obs_cal", "days_obs_val"]] if (HERE / "windows.csv").exists() else None
    res = res.merge(lakes, on="lake", how="left")
    if win is not None:
        res = res.merge(win, on="lake", how="left")
    res["improvement_cal"] = res.rmse_cal_default - res.rmse_cal
    res["improvement_val"] = res.rmse_val_default - res.rmse_val

    reason = pd.Series("", index=res.index)
    reason[res.lake.isin(excl)] = "excluded by user"
    if a.min_days and "days_obs_cal" in res:
        reason[(reason == "") & (res.days_obs_cal < a.min_days)] = f"< {a.min_days} calibration obs days"
    reason[(reason == "") & (res.rmse_cal >= 999)] = "model failed"
    res["dropped_from_median"] = reason
    res.to_csv(HERE / "results_all.csv", index=False)

    keep = res[reason == ""]
    rows = []
    for p in PARAMS + ["Kw", "rmse_cal", "rmse_val", "rmse_cal_default", "rmse_val_default"]:
        s = keep[p]
        rows.append({"parameter": p, "median": s.median(), "q25": s.quantile(0.25), "q75": s.quantile(0.75),
                     "min": s.min(), "max": s.max(), "n_lakes": s.notna().sum(),
                     "n_at_bound": keep.at_bound.fillna("").str.contains(p).sum() if p in PARAMS else None})
    med = pd.DataFrame(rows)
    med.to_csv(HERE / "median_parameters.csv", index=False)

    print(f"{len(res)} lakes calibrated; {len(keep)} used for the median")
    if (reason != "").any():
        print("dropped:\n" + res.loc[reason != "", ["lake", "dropped_from_median"]].to_string(index=False))
    print(med.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
