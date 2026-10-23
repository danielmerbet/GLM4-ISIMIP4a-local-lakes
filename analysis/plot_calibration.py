#!/usr/bin/env python3
"""Figures for the local-lake hydrodynamic calibration and the parameter-transfer tests.

    python3 analysis/plot_calibration.py <results_dir> <outdir> [label]

<results_dir>: calibration/results_v1_unpatched (or results/ once the recalibration is done);
transfer-test tables calibration/eval_*.csv are read when present.
"""
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plotlib as P
plt = P.plt
CAL = Path(__file__).resolve().parent.parent / "calibration"
PARS = [("wind_factor", 0.25, 1.5, 1.0), ("sw_factor", 0.7, 1.3, 1.0), ("kw_scaling", 0.5, 3.0, 1.0),
        ("coef_mix_hyp", 0.1, 2.0, 0.5), ("coef_mix_conv", 0.1, 0.3, 0.2), ("coef_mix_turb", 0.35, 0.65, 0.51)]


def main(results, outdir, label):
    b = pd.concat([pd.read_csv(f) for f in glob.glob(f"{results}/*/best.csv")]).set_index("lake")
    lk = pd.read_csv(CAL / "lakes.csv").set_index("lake")
    b = b.join(lk[["lat", "max_depth_m"]])
    out = Path(outdir); out.mkdir(parents=True, exist_ok=True)

    # 1) per-lake RMSE, default vs calibrated, calibration and validation halves
    d = b.sort_values("rmse_val")
    fig, axs = plt.subplots(1, 2, figsize=(12, 12), sharey=True, constrained_layout=True)
    y = np.arange(len(d))
    for ax, h in zip(axs, ["cal", "val"]):
        ax.hlines(y, d[f"rmse_{h}"], d[f"rmse_{h}_default"], color=P.AXIS, lw=1)
        ax.scatter(d[f"rmse_{h}_default"], y, s=18, color=P.CAT[1], label="GLM defaults", zorder=3)
        ax.scatter(d[f"rmse_{h}"], y, s=18, color=P.CAT[0], label="calibrated", zorder=3)
        ax.set_title(f"{'calibration' if h == 'cal' else 'validation'} half — RMSE (°C)")
        ax.set_xlim(0, max(6, d[f"rmse_{h}_default"].max() * 1.05)); ax.grid(axis="y", visible=False)
    axs[0].set_yticks(y, d.index, fontsize=7); axs[0].legend(fontsize=8, loc="lower right")
    fig.suptitle(f"Hydrodynamic calibration, {len(d)} local lakes ({label})", x=0.01, ha="left", fontweight="bold")
    fig.savefig(out / "cal_rmse_by_lake.png", dpi=140); plt.close(fig)

    # 2) parameter spread within the bounds
    fig, axs = plt.subplots(2, 3, figsize=(12, 6), constrained_layout=True)
    rng = np.random.default_rng(1)
    for ax, (p, lo, hi, dflt) in zip(axs.ravel(), PARS):
        v = b[p].values
        ax.scatter(v, rng.uniform(-0.3, 0.3, len(v)), s=14, color=P.CAT[0], alpha=0.7, lw=0)
        ax.axvline(np.median(v), color=P.INK, lw=1.6); ax.axvline(dflt, color=P.CAT[1], lw=1.6, ls="--")
        ax.text(np.median(v), 0.42, f"median {np.median(v):.2f}", ha="center", fontsize=8, color=P.INK)
        ax.text(dflt, -0.48, f"default {dflt}", ha="center", fontsize=8, color=P.CAT[1])
        ax.set_xlim(lo, hi); ax.set_ylim(-0.55, 0.55); ax.set_yticks([]); ax.set_title(f"{p}  (bounds {lo}–{hi})")
    fig.suptitle(f"Calibrated parameters across lakes ({label}); x-axis = search bounds", x=0.01, ha="left", fontweight="bold")
    fig.savefig(out / "cal_parameter_spread.png", dpi=140); plt.close(fig)

    # 3) parameter vs lake depth (the attribute with the strongest links)
    fig, axs = plt.subplots(2, 3, figsize=(12, 6), constrained_layout=True)
    for ax, (p, lo, hi, dflt) in zip(axs.ravel(), PARS):
        ax.scatter(b.max_depth_m, b[p], s=14, color=P.CAT[0], alpha=0.8, lw=0)
        rho = pd.Series(b[p]).corr(b.max_depth_m, method="spearman")
        ax.set_xscale("log"); ax.set_ylim(lo, hi); ax.set_title(f"{p} vs depth (Spearman {rho:+.2f})")
        ax.set_xlabel("max depth (m)")
    fig.savefig(out / "cal_parameters_vs_depth.png", dpi=140); plt.close(fig)

    # 4) transfer tests: one shared set applied to every lake
    sets = [("GLM defaults", None), ("median", "eval_median_75.csv"), ("median, leave-one-out", "eval_loo_median.csv"),
            ("depth/area rule, leave-one-out", "eval_loo_rule.csv"), ("median + GLM fixes", "eval_median_75_swfix.csv")]
    data, names = [b.rmse_val_default.values], ["GLM defaults"]
    for n, f in sets[1:]:
        if (CAL / f).exists():
            e = pd.read_csv(CAL / f).set_index("lake").reindex(b.index)
            data.append(e.rmse_val_shared.dropna().values); names.append(n)
    data.append(b.rmse_val.values); names.append("own calibration (not possible globally)")
    fig, ax = plt.subplots(figsize=(10, 4.2), constrained_layout=True)
    bp = ax.boxplot(data, vert=False, widths=0.5, patch_artist=True, showfliers=True,
                    medianprops=dict(color=P.INK, lw=1.6), flierprops=dict(marker="o", ms=3, mfc=P.MUTED, mec="none"))
    for k, patch in enumerate(bp["boxes"]):
        patch.set_facecolor(P.CAT[0] if k < len(data) - 1 else P.CAT[2]); patch.set_alpha(0.35); patch.set_edgecolor(P.AXIS)
    for k, dd in enumerate(data):
        ax.text(np.median(dd), k + 1.32, f"{np.median(dd):.2f}", ha="center", fontsize=8, color=P.INK)
    ax.set_yticks(range(1, len(names) + 1), names); ax.set_xlabel("validation RMSE (°C), one value per lake")
    ax.set_title("Transfer test — one parameter set for every local lake"); ax.grid(axis="y", visible=False)
    fig.savefig(out / "cal_transfer_test.png", dpi=140); plt.close(fig)

    # 5) effect of the GLM fixes per lake vs latitude (median parameters)
    a, c = CAL / "eval_median_75.csv", CAL / "eval_median_75_swfix.csv"
    if a.exists() and c.exists():
        e = pd.read_csv(a).set_index("lake").join(pd.read_csv(c).set_index("lake"), rsuffix="_fix").join(lk[["lat"]])
        e["d"] = e.rmse_val_shared_fix - e.rmse_val_shared
        fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
        ax.axhline(0, color=P.AXIS, lw=1)
        ax.scatter(e.lat, e.d, s=22, color=np.where(e.d < 0, P.CAT[0], P.CAT[1]), lw=0)
        for nm, r in e[(e.d.abs() > 0.6)].iterrows():
            ax.annotate(nm, (r.lat, r.d), xytext=(4, 2), textcoords="offset points", fontsize=7, color=P.INK2)
        ax.set_xlabel("latitude (°)"); ax.set_ylabel("Δ validation RMSE (°C)")
        ax.set_title("Effect of the GLM fixes per lake (median parameters): below 0 = better")
        fig.savefig(out / "cal_fix_effect_by_latitude.png", dpi=140); plt.close(fig)
    print("written to", out)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "")
