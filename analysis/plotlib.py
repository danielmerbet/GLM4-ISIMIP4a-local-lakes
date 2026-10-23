"""Shared plotting style and per-lake figures for the GLM4-AED3 outputs (reduced.nc).

Style follows the dataviz reference palette: categorical hues in fixed order, single-hue
sequential ramps for heatmaps (temperature = orange, oxygen = blue, chlorophyll = aqua),
recessive grids, 2 px lines, direct labels where there are few series.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

SURFACE, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
BLUE = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]
SEQ_T = LinearSegmentedColormap.from_list("temp", ["#fde8dd", "#f4a27a", "#eb6834", "#b8461c", "#6e2508"])
SEQ_O = LinearSegmentedColormap.from_list("oxy", BLUE)
SEQ_C = LinearSegmentedColormap.from_list("chl", ["#d9f3e9", "#7fd4b2", "#1baf7a", "#11815a", "#0a4f37"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "text.color": INK, "axes.labelcolor": INK2, "axes.edgecolor": AXIS, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "font.size": 9, "axes.titlesize": 10,
    "axes.titleweight": "bold", "axes.titlelocation": "left", "lines.linewidth": 1.4,
    "legend.frameon": False, "font.family": "sans-serif", "axes.axisbelow": True})


def datefmt(ax):
    loc = mdates.AutoDateLocator(maxticks=6)
    ax.xaxis.set_major_locator(loc); ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(loc))

from reducedio import MOL, TN, TP, load  # noqa: E402,F401  (data loading lives in reducedio)


def _heat(ax, t, depth, z, cmap, label, vmin=None, vmax=None):
    """Daily values (time, depth) on the saved levels; each cell spans its level's depth range."""
    m = ax.pcolormesh(t, depth, np.ma.masked_invalid(z.T), cmap=cmap, vmin=vmin, vmax=vmax,
                      shading="nearest", rasterized=True)
    ax.set_ylim(depth[-1], 0); ax.set_ylabel("depth (m)"); ax.grid(False)
    cb = plt.colorbar(m, ax=ax, pad=0.01, fraction=0.03); cb.set_label(label, color=INK2); cb.outline.set_visible(False)
    return m


def lake_physics(path, out, title, obs=None, years=None):
    """Temperature heatmap + thermocline + ice; surface/bottom series (+obs); monthly profiles;
    stratified days per year by definition."""
    t, depth, x, s = load(path)
    a = t >= pd.Timestamp(s.get("analysis_start", "2001-01-01"))
    if years:
        a &= (t.year >= years[0]) & (t.year <= years[1])
    tt = t[a]
    fig = plt.figure(figsize=(12, 10), constrained_layout=True)
    gs = fig.add_gridspec(3, 3, height_ratios=[1.3, 1, 1])
    ax = fig.add_subplot(gs[0, :])
    tv = x["temp"][a]
    _heat(ax, tt, depth, tv, SEQ_T, "water temperature (°C)", vmin=max(0.0, np.nanpercentile(tv, 0.5)),
          vmax=np.nanpercentile(tv, 99.5))
    ax.plot(tt, x["thermodepth"][a], color=INK, lw=0.6, label="thermocline (max density gradient)")
    ice = x["blue_ice_thickness"][a] + x["white_ice_thickness"][a]
    if np.nanmax(ice) > 0:
        ax.fill_between(tt, 0, np.where(ice > 0, depth[-1] * 0.02, np.nan), color=BLUE[2], step="mid", label="ice cover")
    ax.legend(loc="lower left", fontsize=8)
    ax.set_title(f"{title} — temperature profile, daily")
    ax2 = fig.add_subplot(gs[1, :2])
    ax2.plot(tt, x["temp_surface_layer"][a], color=CAT[1], lw=0.8, label="surface layer")
    ax2.plot(tt, x["temp_bottom_layer"][a], color=CAT[0], lw=0.8, label="bottom layer")
    if obs is not None and len(obs):
        o = obs[(obs.datetime >= tt[0]) & (obs.datetime <= tt[-1])]
        su = o[o.depth <= 1].groupby("datetime").value.mean()
        bt = o[o.depth >= depth[-1] - 2].groupby("datetime").value.mean()
        ax2.scatter(su.index, su.values, s=6, color=CAT[1], edgecolor=INK, lw=0.3, label="observed ≤1 m", zorder=3)
        if len(bt):
            ax2.scatter(bt.index, bt.values, s=6, color=CAT[0], edgecolor=INK, lw=0.3, label="observed near bottom", zorder=3)
    ax2.set_ylabel("°C"); ax2.set_title("surface and bottom temperature"); ax2.legend(fontsize=8, ncol=2)
    ax3 = fig.add_subplot(gs[1, 2])
    if np.nanmax(ice) > 0:
        ax3.plot(tt, ice, color=BLUE[4], lw=0.8); datefmt(ax3)
        ax3.text(0.02, 0.95, f"{(ice > 0).sum() / (len(tt) / 365.25):.0f} ice days/yr", transform=ax3.transAxes, va="top", color=INK2, fontsize=8)
    else:
        ax3.text(0.5, 0.5, "no ice", ha="center", va="center", color=MUTED, transform=ax3.transAxes); ax3.set_axis_off()
    ax3.set_title("ice thickness (m)")
    ax4 = fig.add_subplot(gs[2, 0])
    T = pd.DataFrame(x["temp"][a], index=tt, columns=depth)
    mon = T.groupby(T.index.month).mean()
    for k, mth in enumerate([2, 5, 8, 11]):
        ax4.plot(mon.loc[mth], depth, color=CAT[k], lw=1.6, label=pd.Timestamp(2001, mth, 1).strftime("%b"))
    ax4.set_ylim(depth[-1], 0); ax4.set_xlabel("°C"); ax4.set_ylabel("depth (m)"); ax4.set_title("mean profile by month")
    ax4.legend(fontsize=8)
    ax5 = fig.add_subplot(gs[2, 1])
    if obs is not None and len(obs):
        o = obs[(obs.datetime >= tt[0]) & (obs.datetime <= tt[-1])].copy()
        idx = pd.Series(np.arange(len(tt)), index=tt)
        o = o[o.datetime.isin(tt)]
        ti = idx.reindex(o.datetime.values).values.astype(int)
        o["sim"] = [np.interp(z, depth, x["temp"][a][i]) for z, i in zip(o.depth.values, ti)]
        o = o.dropna(subset=["sim"])
        ax5.scatter(o.value, o.sim, s=3, color=CAT[0], alpha=0.35, lw=0)
        lim = [min(o.value.min(), o.sim.min()), max(o.value.max(), o.sim.max())]
        ax5.plot(lim, lim, color=MUTED, lw=1, ls="--")
        rmse = float(np.sqrt(((o.sim - o.value) ** 2).mean())); bias = float((o.sim - o.value).mean())
        ax5.set_title(f"sim vs obs, all depths (RMSE {rmse:.2f}, bias {bias:+.2f} °C)")
        ax5.set_xlabel("observed °C"); ax5.set_ylabel("simulated °C")
    else:
        ax5.text(0.5, 0.5, "no observations", ha="center", va="center", color=MUTED, transform=ax5.transAxes)
        ax5.set_axis_off()
    ax6 = fig.add_subplot(gs[2, 2])
    defs = [("strat_dens0p005", "Δρ≥0.005 (ISIMIP4T)"), ("strat_dens0p01", "Δρ>0.01 & Ts>Tb"),
            ("strat_isimip3", "ISIMIP3 (Δρ≥0.1)"), ("strat_dens0p1", "Δρ>0.1 & Ts>Tb"), ("strat_temp1", "ΔT>1 K")]
    vals = [np.nansum(x[k][a]) / (a.sum() / 365.25) for k, _ in defs]
    ax6.barh(range(len(defs)), vals, color=CAT[0], height=0.6)
    ax6.set_yticks(range(len(defs)), [l for _, l in defs]); ax6.invert_yaxis(); ax6.set_xlim(0, 366)
    for i, v in enumerate(vals):
        ax6.text(v + 4, i, f"{v:.0f}", va="center", color=INK2, fontsize=8)
    ax6.set_title("stratified days per year"); ax6.grid(axis="y", visible=False)
    fig.savefig(out, dpi=140); plt.close(fig)


def lake_wq(path, out, title, years=None):
    """Oxygen and chlorophyll heatmaps; monthly epilimnion/hypolimnion cycles of the WQ variables."""
    t, depth, x, s = load(path)
    a = t >= pd.Timestamp(s.get("analysis_start", "2001-01-01"))
    if years:
        a &= (t.year >= years[0]) & (t.year <= years[1])
    tt = t[a]
    fig = plt.figure(figsize=(12, 10.5), constrained_layout=True)
    gs = fig.add_gridspec(4, 3, height_ratios=[1, 1, 0.9, 0.9])
    ax = fig.add_subplot(gs[0, :])
    _heat(ax, tt, depth, x["OXY_oxy"][a] * MOL["O2"] / 1000, SEQ_O, "dissolved oxygen (mg/L)", vmin=0)
    ax.plot(tt, x["thermodepth"][a], color=INK, lw=0.5)
    ax.set_title(f"{title} — dissolved oxygen")
    ax = fig.add_subplot(gs[1, :])
    _heat(ax, tt, depth, x["PHY_tchla"][a], SEQ_C, "chlorophyll-a (µg/L)", vmin=0,
          vmax=np.nanpercentile(x["PHY_tchla"][a], 99))
    ax.set_title("chlorophyll-a")
    def clim(v):
        e = pd.Series(v[0], index=tt); h = pd.Series(v[1], index=tt)
        return e.groupby(e.index.month).mean(), h.groupby(h.index.month).mean()
    panels = [
        ("chlorophyll-a (µg/L)", (x["PHY_tchla_epi"][a], x["PHY_tchla_hypo"][a])),
        ("total P (µg/L)", (sum(x[v + "_epi"][a] for v in TP) * MOL["P"], sum(x[v + "_hypo"][a] for v in TP) * MOL["P"])),
        ("total N (mg/L)", (sum(x[v + "_epi"][a] for v in TN) * MOL["N"] / 1000, sum(x[v + "_hypo"][a] for v in TN) * MOL["N"] / 1000)),
        ("dissolved O2 (mg/L)", (x["OXY_oxy_epi"][a] * MOL["O2"] / 1000, x["OXY_oxy_hypo"][a] * MOL["O2"] / 1000)),
        ("pH", (x["CAR_pH_epi"][a], x["CAR_pH_hypo"][a])),
        ("zooplankton (mmol C/m³)", (x["ZOO_daphnia_epi"][a], x["ZOO_daphnia_hypo"][a])) if "ZOO_daphnia_epi" in x else None]
    for k, p in enumerate([p for p in panels if p]):
        ax = fig.add_subplot(gs[2 + k // 3, k % 3])
        e, h = clim(p[1])
        ax.plot(e.index, e.values, color=CAT[1], lw=2, marker="o", ms=4, label="epilimnion")
        ax.plot(h.index, h.values, color=CAT[0], lw=2, marker="o", ms=4, label="hypolimnion")
        ax.set_xticks(range(1, 13), list("JFMAMJJASOND")); ax.set_title(p[0])
        if k == 0:
            ax.legend(fontsize=8)
    fig.savefig(out, dpi=140); plt.close(fig)


def phyto_groups(path, out, title, years=None):
    """Monthly epilimnion biomass by group (phytoplankton + zooplankton) and limitation context."""
    t, depth, x, s = load(path)
    a = t >= pd.Timestamp(s.get("analysis_start", "2001-01-01"))
    tt = t[a]
    fig, axs = plt.subplots(1, 2, figsize=(12, 3.6), constrained_layout=True)
    groups = [k for k in ("PHY_diatom", "PHY_cyano", "PHY_green") if k + "_epi" in x]
    for i, g in enumerate(groups + (["ZOO_daphnia"] if "ZOO_daphnia_epi" in x else [])):
        sr = pd.Series(x[g + "_epi"][a], index=tt)
        axs[0].plot(sr.index, sr.values, color=CAT[i], lw=0.8, label=g.split("_", 1)[1])
        m = sr.groupby(sr.index.month).mean()
        axs[1].plot(m.index, m.values, color=CAT[i], lw=2, marker="o", ms=4, label=g.split("_", 1)[1])
    axs[0].set_title(f"{title} — epilimnion biomass (mmol C/m³), daily"); axs[0].legend(fontsize=8, ncol=4)
    axs[1].set_title("monthly mean"); axs[1].set_xticks(range(1, 13), list("JFMAMJJASOND"))
    fig.savefig(out, dpi=140); plt.close(fig)
