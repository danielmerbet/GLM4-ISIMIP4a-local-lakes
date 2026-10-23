#!/usr/bin/env python3
"""Refill the generated sections of README.md from the result tables.

    python3 local_lakes/fill_readme.py

Reads results/scores.csv (plot_local_lakes.py), results/qc.csv (qc_isimip.py), calibration/results_all.csv
(collect_results.py) and local_lakes/lake_sites.csv; rewrites the blocks between
<!-- BEGIN:name --> and <!-- END:name --> markers in README.md (summary, calibration, lakeheatf, evaluation,
lake_table). Nothing else in the README changes.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def blocks():
    sc = pd.read_csv(ROOT / "results" / "scores.csv")
    ra = pd.read_csv(ROOT / "calibration" / "results_all.csv")
    qc = pd.read_csv(ROOT / "results" / "qc.csv")
    sites = pd.read_csv(ROOT / "local_lakes" / "lake_sites.csv").set_index("lake")
    n, ntot = len(sc), len(sites)
    med = lambda c: sc[c].median()
    oth = sc.dropna(subset=["rmse_other_years"])
    missing = sorted(set(sites.index) - set(sc.lake))
    out = {}
    out["summary"] = (
        f"Results at a glance: {n} of {ntot} lakes delivered"
        + (f" ({', '.join(missing)} to follow)" if missing else "")
        + f". Median water-temperature RMSE of the delivered files: **{med('rmse_calibration'):.2f} °C** in the calibration "
        f"years, **{med('rmse_validation'):.2f} °C** in the validation years and **{oth.rmse_other_years.median():.2f} °C** in "
        f"years the calibration never used ({len(oth)} lakes); GLM with its default parameters: "
        f"{ra.rmse_val_default.median():.2f} °C in the validation years.")
    better = int((ra.rmse_val < ra.rmse_val_default).sum())
    pars = ["wind_factor", "sw_factor", "kw_scaling", "coef_mix_hyp", "coef_mix_conv", "coef_mix_turb"]
    out["calibration"] = (
        f"Over the {len(ra)} calibrated lakes (median, interquartile range in brackets): calibration RMSE "
        f"**{ra.rmse_cal.median():.2f} °C** ({ra.rmse_cal.quantile(.25):.2f}-{ra.rmse_cal.quantile(.75):.2f}), validation "
        f"**{ra.rmse_val.median():.2f} °C** ({ra.rmse_val.quantile(.25):.2f}-{ra.rmse_val.quantile(.75):.2f}); with GLM's default "
        f"parameters {ra.rmse_cal_default.median():.2f} and {ra.rmse_val_default.median():.2f} °C. The calibration improved the "
        f"validation RMSE in {better} of {len(ra)} lakes. Calibrated parameter medians: "
        + ", ".join(f"`{p}` {ra[p].median():.3g}" for p in pars)
        + f"; {int((ra.at_bound.fillna('') != '').sum())} lakes have at least one parameter at a bound (`at_bound` in `results_all.csv`).")
    h = qc.set_index("lake").lakeheatf_mean_Wm2.dropna()
    out["lakeheatf"] = (
        "**Note on `lakeheatf`.** Averaged over 1901-2025, the net surface heat flux is not zero in every lake although the "
        f"lakes do not warm at that rate: median {h.median():+.1f} W m⁻², from {h.min():+.0f} ({h.idxmin()}) to "
        f"{h.max():+.0f} W m⁻² ({h.idxmax()}); warm reservoirs are positive, cold ice-covered lakes negative. The flux is "
        "computed exactly as defined (surface fluxes of GLM, positive into the lake); the difference must leave the lake by "
        "other routes in GLM (e.g. the evaporation-compensating inflow and the overflow, or energy taken up by the ice). "
        "This is being investigated; use `lakeheatf` with care for long-term energy budgets.")
    ev = ["Median RMSE (°C) of the delivered files per lake, and the bias (model − observation):", "",
          "| Observations | Lakes | Median RMSE | Median bias |", "|---|---|---|---|"]
    for k, lab in [("calibration", "calibration half"), ("validation", "validation half"), ("other_years", "all other years 1901-2025")]:
        ev.append(f"| {lab} | {sc[f'rmse_{k}'].notna().sum()} | {sc[f'rmse_{k}'].median():.2f} | {sc[f'bias_{k}'].median():+.2f} |")
    ev += ["", "The delivered run differs from the calibration runs (125 years with a repeated-forcing spin-up instead of a "
           "20-year spin-up and an observed initial profile), so its calibration/validation scores differ slightly from `best.csv`.",
           "", "File checks (`local_lakes/qc_isimip.py`, `results/qc.csv`): every lake has the 13 files, 45,656 days, the required "
           "attributes and matching checksums; all values in plausible ranges and finite, except `thermodepth` on days without a "
           "thermocline. `watertemp` has missing values only for levels below the lake bottom on days when the water level is below "
           f"the spill level (ice, evaporation): median {qc.days_bottom_level_empty_pct.median():.1f} % of days, at most "
           f"{int(qc.max_empty_levels.max())} levels."]
    out["evaluation"] = "\n".join(ev)
    rows = ["| Lake | Depth (m) | RMSE cal / val / other (°C) | Figures |", "|---|---|---|---|"]
    for _, r in sc.sort_values("lake").iterrows():
        o = "–" if np.isnan(r.rmse_other_years) else f"{r.rmse_other_years:.2f}"
        rows.append(f"| {sites.loc[r.lake, 'title']} | {r.max_depth_m:g} | {r.rmse_calibration:.2f} / {r.rmse_validation:.2f} / {o} | "
                    f"[temperature](figures/lakes/{r.lake}_temperature.png), [all variables](figures/lakes/{r.lake}_variables.png) |")
    out["lake_table"] = "\n".join(rows)
    return out


def main():
    p = ROOT / "README.md"
    t = p.read_text()
    for name, text in blocks().items():
        pat = re.compile(rf"(<!-- BEGIN:{name} -->\n).*?(\n<!-- END:{name} -->)", re.S)
        if not pat.search(t):
            raise SystemExit(f"README.md: no markers for {name}")
        t = pat.sub(lambda m: m.group(1) + text + m.group(2), t)
    p.write_text(t)
    print("README.md: generated sections refilled")


if __name__ == "__main__":
    main()
