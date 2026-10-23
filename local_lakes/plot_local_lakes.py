#!/usr/bin/env python3
"""Figures and scores of the delivered ISIMIP4a local-lake files against the observations.

    python3 local_lakes/plot_local_lakes.py <isimip dir> <out dir> [Lake ...]

<isimip dir>/<Lake>/ holds the 13 NetCDF files of write_isimip.py. The scores use the files as
delivered (1901-2025 run with spin-up), not the calibration runs:
  * calibration and validation observations: calibration/runs/<Lake>/obs_wtemp_{cal,val}.csv
    (the same sets the calibration used; first and second half of the last <= 10 years);
  * all other observations of the lake in 1901-2025 (calibration/obs/<Lake>_temp.csv): years
    the calibration never saw;
  * observation day D is compared with the file's day D; the model is interpolated linearly in
    depth between the levlak levels (the observation depth is below the surface).
Outputs:
  <out>/lakes/<Lake>_temperature.png  simulated profile with the observations on top, surface and
                                      bottom series, obs vs model, error by depth
  <out>/lakes/<Lake>_variables.png    all 13 variables: 1991-2020 seasonal cycles, 1901-2025 annual series
  <out>/scores.csv                    RMSE, bias, n per lake and observation set
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import netCDF4 as nc  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
CAL = ROOT / "calibration"
SETS = {"calibration": "#2a78d6", "validation": "#eb6834", "other years": "#8a8f98"}
T0 = pd.Timestamp("1901-01-01")
SCALARS = [("surftemp", "Surface temperature", "°C", -273.15), ("bottemp", "Bottom temperature", "°C", -273.15),
           ("icethick", "Ice thickness", "m", 0), ("ice", "Ice cover", "fraction of days", 0),
           ("strat", "Stratified", "fraction of days", 0), ("thermodepth", "Thermocline depth", "m", 0),
           ("swup", "Upward shortwave", "W m⁻²", 0), ("lwup", "Upward longwave", "W m⁻²", 0),
           ("sensheatf", "Sensible heat flux (up)", "W m⁻²", 0), ("latentheatf", "Latent heat flux (up)", "W m⁻²", 0),
           ("lakeheatf", "Net heat flux into lake", "W m⁻²", 0), ("momf", "Momentum flux", "kg m⁻¹ s⁻²", 0)]
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titlesize": 8.5, "axes.titleweight": "bold", "figure.dpi": 110})


def read_var(folder, var):
    f = next(Path(folder).glob(f"*_{var}_*_daily_1901_2025.nc"))
    d = nc.Dataset(f); d.set_auto_mask(False)
    t = T0 + pd.to_timedelta(d["time"][:], unit="D")
    v = d[var][:].astype("f8")
    v[v > 1e19] = np.nan
    z = d["levlak"][:].astype("f8") if "levlak" in d.variables else None
    d.close()
    return t, (v[:, :, 0, 0] if v.ndim == 4 else v[:, 0, 0]), z


def obs_sets(lake):
    run = CAL / "runs" / lake
    out = []
    for s, f in (("calibration", "obs_wtemp_cal.csv"), ("validation", "obs_wtemp_val.csv")):
        o = pd.read_csv(run / f)
        out.append(pd.DataFrame({"date": pd.to_datetime(o.datetime).dt.normalize(), "depth": o.Depth_meter,
                                 "obs": o.Water_Temperature_celsius, "set": s}))
    used = pd.concat(out)
    a = pd.read_csv(CAL / "obs" / f"{lake}_temp.csv")
    a = pd.DataFrame({"date": pd.to_datetime(a.datetime).dt.normalize(), "depth": a.depth, "obs": a.value,
                      "set": "other years"})
    key = lambda d: d.date.astype(str) + "_" + d.depth.round(3).astype(str)
    a = a[~key(a).isin(set(key(used)))]
    a = a[(a.date >= T0) & (a.date <= pd.Timestamp("2025-12-31"))]
    return pd.concat([used, a], ignore_index=True)


def match(obs, t, wt, z):
    """Model at each observation (same day, linear in depth between the levels)."""
    idx = pd.Series(np.arange(len(t)), index=t)
    obs = obs[obs.date.isin(idx.index) & (obs.depth >= 0) & (obs.depth <= z[-1])].copy()
    rows = idx.loc[obs.date].values
    sim = np.empty(len(obs))
    for r in np.unique(rows):
        m = rows == r
        prof = wt[r]
        ok = np.isfinite(prof)
        sim[m] = np.interp(obs.depth.values[m], z[ok], prof[ok]) if ok.sum() > 1 else np.nan
    obs["sim"] = sim
    return obs.dropna(subset=["sim"])


def scores(m):
    e = m.sim - m.obs
    return {"rmse": float(np.sqrt((e ** 2).mean())) if len(e) else np.nan, "bias": float(e.mean()) if len(e) else np.nan,
            "n": int(len(e)), "days": int(m.date.nunique())}


def fig_temperature(lake, title, t, wt, z, m, win, best, out):
    c0, c1 = pd.Timestamp(win.cal_start), pd.Timestamp(win.end)
    sel = (t >= c0) & (t <= c1)
    fig = plt.figure(figsize=(11, 10.5), constrained_layout=True)
    gs = fig.add_gridspec(4, 3, height_ratios=[0.9, 0.9, 1, 1.1])
    tt, ww = t[sel], wt[sel]
    w = m[(m.date >= c0) & (m.date <= c1)]
    vmin, vmax = np.nanpercentile(ww, 1), np.nanpercentile(ww, 99)
    ax = fig.add_subplot(gs[0, :])
    pc = ax.pcolormesh(tt, z, ww.T, cmap="RdYlBu_r", vmin=vmin, vmax=vmax, shading="nearest", rasterized=True)
    ax.set_title("Model: simulated water temperature (delivered file, levlak levels)")
    ax2 = fig.add_subplot(gs[1, :], sharex=ax, sharey=ax)
    ax2.set_facecolor("#eef0f3")
    ax2.scatter(w.date, w.depth, c=w.obs, cmap="RdYlBu_r", vmin=vmin, vmax=vmax, marker="s", lw=0,
                s=max(1.5, min(12, 4000 / max(1, w.date.nunique()))), rasterized=True)
    ax2.set_title("Observations (same colour scale; grey = no observation)")
    for a in (ax, ax2):
        a.axvline(pd.Timestamp(win.val_start), color="k", lw=1, ls="--")
        a.set_ylabel("depth (m)")
        fig.colorbar(pc, ax=a, label="°C", pad=0.01)
    ax.set_ylim(z[-1], 0); ax.set_xlim(c0, c1)
    ax.text(c0, z[-1] * 0.04, " calibration", va="top", fontsize=7.5)
    ax.text(pd.Timestamp(win.val_start), z[-1] * 0.04, " validation", va="top", fontsize=7.5)
    # surface and bottom series at the most observed shallow / deep depth of the window
    dep = w.groupby("depth").size()
    shallow = dep[dep.index <= dep.index.min() + 1].idxmax() if len(dep) else None   # within 1 m of the shallowest
    deep = dep[dep.index >= dep.index.max() - max(1.0, 0.1 * dep.index.max())].idxmax() if len(dep) else None
    for k, (dz, lab) in enumerate([(shallow, "near surface"), (deep, "near bottom")]):
        a = fig.add_subplot(gs[2, :2] if k == 0 else gs[3, :2])
        if dz is None:
            a.set_axis_off(); continue
        sim = np.array([np.interp(dz, z[np.isfinite(p)], p[np.isfinite(p)]) if np.isfinite(p).sum() > 1 else np.nan
                        for p in wt[sel]])
        a.plot(tt, sim, color="#3b3f46", lw=0.8, label="model")
        o = w[np.isclose(w.depth, dz)]
        for s, col in SETS.items():
            q = o[o.set == s]
            if len(q):
                a.scatter(q.date, q.obs, s=6, color=col, label=f"obs ({s})", zorder=3)
        a.axvline(pd.Timestamp(win.val_start), color="k", lw=0.8, ls="--")
        a.set_title(f"{lab}: {dz:g} m"); a.set_ylabel("°C"); a.legend(fontsize=6.5, ncol=4, loc="upper left", frameon=False)
    a = fig.add_subplot(gs[2, 2])
    lim = [min(m.obs.min(), m.sim.min()) - 1, max(m.obs.max(), m.sim.max()) + 1]
    for s, col in list(SETS.items())[::-1]:
        q = m[m.set == s]
        if len(q):
            a.scatter(q.obs, q.sim, s=3, alpha=0.35, color=col, label=s, rasterized=True)
    a.plot(lim, lim, color="k", lw=0.8); a.set_xlim(lim); a.set_ylim(lim)
    a.set_xlabel("observed (°C)"); a.set_ylabel("model (°C)"); a.set_title("All observations 1901-2025")
    a.legend(fontsize=6.5, markerscale=3, frameon=False)
    a = fig.add_subplot(gs[3, 2]); a.axis("off")
    lines = ["RMSE (°C), delivered run:"]
    for s in SETS:
        q = m[m.set == s]
        if len(q):
            sc = scores(q)
            lines.append(f"  {s:12s} {sc['rmse']:.2f}  bias {sc['bias']:+.2f}  ({sc['days']} days)")
    lines += ["", "Calibration (best.csv):", f"  calibration  {best.rmse_cal:.2f}  (GLM defaults {best.rmse_cal_default:.2f})",
              f"  validation   {best.rmse_val:.2f}  (GLM defaults {best.rmse_val_default:.2f})", "",
              "Parameters:", f"  wind_factor {best.wind_factor:.2f}   sw_factor {best.sw_factor:.2f}",
              f"  Kw {best.Kw:.3f} m⁻¹ (= {best.kw_scaling:.2f} × Håkanson)",
              f"  coef_mix_hyp {best.coef_mix_hyp:.2f}  conv {best.coef_mix_conv:.2f}  turb {best.coef_mix_turb:.2f}"]
    a.text(0, 1, "\n".join(lines), va="top", family="monospace", fontsize=7.2)
    fig.suptitle(f"{title} — GLM 4, ISIMIP4a GSWP3-ERA5 obsclim", fontsize=11, fontweight="bold")
    fig.savefig(out / f"{lake}_temperature.png"); plt.close(fig)


def fig_variables(lake, title, folder, out):
    fig, axs = plt.subplots(4, 4, figsize=(13, 10), constrained_layout=True)
    clim = slice("1991-01-01", "2020-12-31")
    for ax, (v, lab, unit, off) in zip(axs.flat[:12], SCALARS):
        t, x, _ = read_var(folder, v)
        s = pd.Series(x + off, index=t)
        c = s[clim]
        g = c.groupby(c.index.dayofyear)
        q = g.quantile([0.05, 0.5, 0.95]).unstack()
        ax.fill_between(q.index, q[0.05], q[0.95], color="#9ec5f4", lw=0, label="5-95 %")
        ax.plot(q.index, g.mean(), color="#184f95", lw=1, label="mean")
        ax.set_title(lab); ax.set_ylabel(unit); ax.set_xlim(1, 366)
        ax.set_xticks([1, 91, 182, 274], ["Jan", "Apr", "Jul", "Oct"])
    axs.flat[0].legend(fontsize=6.5, frameon=False, loc="best")
    t, wt, z = read_var(folder, "watertemp")
    s = pd.DataFrame(wt[(t >= "1991-01-01") & (t <= "2020-12-31")] - 273.15,
                     index=t[(t >= "1991-01-01") & (t <= "2020-12-31")])
    c = s.groupby(s.index.dayofyear).mean()
    ax = axs.flat[12]
    pc = ax.pcolormesh(c.index, z, c.values.T, cmap="RdYlBu_r", shading="nearest")
    ax.set_ylim(z[-1], 0); ax.set_ylabel("depth (m)"); ax.set_title("Water temperature, mean year")
    ax.set_xticks([1, 91, 182, 274], ["Jan", "Apr", "Jul", "Oct"]); fig.colorbar(pc, ax=ax, label="°C")
    ann = lambda v, f: pd.Series(read_var(folder, v)[1], index=read_var(folder, v)[0]).resample("YS").agg(f)
    ax = axs.flat[13]
    for v, col, lab in (("surftemp", "#eb6834", "surface"), ("bottemp", "#2a78d6", "bottom")):
        a = ann(v, "mean") - 273.15
        ax.plot(a.index.year, a.values, color=col, lw=0.8, label=lab)
    ax.set_title("Annual mean temperature 1901-2025"); ax.set_ylabel("°C"); ax.legend(fontsize=6.5, frameon=False)
    ax = axs.flat[14]
    for v, col, lab in (("ice", "#2a78d6", "ice-covered days"), ("strat", "#eb6834", "stratified days")):
        a = ann(v, "sum")
        ax.plot(a.index.year, a.values, color=col, lw=0.8, label=lab)
    ax.set_title("Days per year 1901-2025"); ax.set_ylabel("days"); ax.legend(fontsize=6.5, frameon=False)
    ax = axs.flat[15]
    a = ann("lakeheatf", "mean")
    ax.plot(a.index.year, a.values, color="#3b3f46", lw=0.8); ax.axhline(0, color="k", lw=0.5)
    ax.set_title("Annual mean net heat flux into lake"); ax.set_ylabel("W m⁻²")
    fig.suptitle(f"{title} — all 13 ISIMIP variables (seasonal cycles: 1991-2020)", fontsize=11, fontweight="bold")
    fig.savefig(out / f"{lake}_variables.png"); plt.close(fig)


def main(isimip, outdir, lakes):
    isimip, out = Path(isimip), Path(outdir)
    (out / "lakes").mkdir(parents=True, exist_ok=True)
    sites = pd.read_csv(ROOT / "local_lakes" / "lake_sites.csv").set_index("lake")
    wins = pd.read_csv(CAL / "windows.csv").set_index("lake")
    lakes = lakes or sorted(p.name for p in isimip.iterdir() if (p / "SHA256SUMS").exists())
    rows, pairs = [], []
    for lake in lakes:
        folder = isimip / lake
        best = pd.read_csv(CAL / "results" / lake / "best.csv").iloc[0]
        t, wt, z = read_var(folder, "watertemp")
        wt = wt - 273.15                                   # K -> °C
        m = match(obs_sets(lake), t, wt, z)
        title = sites.loc[lake, "title"]
        fig_temperature(lake, title, t, wt, z, m, wins.loc[lake], best, out / "lakes")
        fig_variables(lake, title, folder, out / "lakes")
        r = {"lake": lake, "max_depth_m": z[-1]}
        for s in SETS:
            for k, v in scores(m[m.set == s]).items():
                r[f"{k}_{s.replace(' ', '_')}"] = v
        r.update({k: best[k] for k in ("rmse_cal", "rmse_val", "rmse_cal_default", "rmse_val_default")})
        rows.append(r)
        pairs.append(m.assign(lake=lake)[["lake", "set", "obs", "sim"]])
        print(f"{lake}: RMSE cal {r['rmse_calibration']:.2f}, val {r['rmse_validation']:.2f}, "
              f"other {r['rmse_other_years']:.2f} °C", flush=True)
    sc = pd.DataFrame(rows)
    f = out / "scores.csv"
    if f.exists() and lakes:
        old = pd.read_csv(f); sc = pd.concat([old[~old.lake.isin(sc.lake)], sc]).sort_values("lake")
    sc.to_csv(f, index=False, float_format="%.4g")
    if len(sc) > 1:
        summary(sc, pd.concat(pairs), sites, out)


def summary(sc, pairs, sites, out):
    """Cross-lake figures: RMSE per lake, map, all observations vs model."""
    d = sc.set_index("lake").sort_values("rmse_validation")
    fig, ax = plt.subplots(figsize=(8, max(6, 0.17 * len(d))), constrained_layout=True)
    y = np.arange(len(d))
    ax.hlines(y, d.rmse_val_default, d.rmse_validation, color="#c9ccd1", lw=1)
    ax.scatter(d.rmse_val_default, y, s=16, color="#8a8f98", label="GLM defaults (calibration window)", zorder=3)
    ax.scatter(d.rmse_validation, y, s=16, color=SETS["validation"], label="delivered run: validation", zorder=3)
    ax.scatter(d.rmse_calibration, y, s=16, color=SETS["calibration"], marker="D", label="delivered run: calibration", zorder=3)
    ax.scatter(d.rmse_other_years, y, s=14, facecolors="none", edgecolors="k", lw=0.7, label="delivered run: other years", zorder=3)
    ax.set_yticks(y, [sites.loc[k, "title"] for k in d.index], fontsize=6.5); ax.set_ylim(-1, len(d))
    ax.set_xlabel("RMSE of water temperature (°C)"); ax.grid(axis="x", color="#eceef1")
    ax.legend(fontsize=7, loc="lower right", frameon=True)
    ax.set_title(f"Water temperature error per lake ({len(d)} lakes); sorted by validation RMSE")
    fig.savefig(out / "summary_rmse_by_lake.png"); plt.close(fig)

    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    fig = plt.figure(figsize=(11, 5.6), constrained_layout=True)
    ax = fig.add_subplot(projection=ccrs.Robinson())
    ax.set_global(); ax.add_feature(cfeature.LAND, color="#eceef1"); ax.coastlines(lw=0.4, color="#8a8f98")
    s = sites.loc[d.index]
    lat = np.where(s.index == "Bosumtwi", s.calib_lat, s.site_lat); lon = np.where(s.index == "Bosumtwi", s.calib_lon, s.site_lon)
    sca = ax.scatter(lon, lat, c=d.rmse_validation, cmap="viridis_r", vmin=0, vmax=3, s=12 + 3 * np.sqrt(s.max_depth_m),
                     edgecolors="k", lw=0.4, transform=ccrs.PlateCarree(), zorder=3)
    fig.colorbar(sca, ax=ax, shrink=0.7, label="validation RMSE (°C)", extend="max")
    ax.set_title(f"{len(d)} ISIMIP local lakes: validation RMSE of the delivered run (dot size grows with depth)")
    fig.savefig(out / "summary_map.png"); plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(12, 4.2), constrained_layout=True, sharex=True, sharey=True)
    for ax, s in zip(axs, SETS):
        q = pairs[pairs.set == s]
        hb = ax.hexbin(q.obs, q.sim, gridsize=70, bins="log", cmap="Blues", mincnt=1, extent=(-2, 32, -2, 32))
        ax.plot([-2, 32], [-2, 32], color="k", lw=0.8)
        e = q.sim - q.obs
        ax.set_title(f"{s}: {q.lake.nunique()} lakes, {len(q):,} obs\nRMSE {np.sqrt((e ** 2).mean()):.2f} °C, bias {e.mean():+.2f} °C")
        ax.set_xlabel("observed (°C)")
    axs[0].set_ylabel("model (°C)"); fig.colorbar(hb, ax=axs, label="observations per cell", shrink=0.8)
    fig.savefig(out / "summary_obs_vs_model.png"); plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3:])
