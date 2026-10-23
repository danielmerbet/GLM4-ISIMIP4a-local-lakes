"""Output levels and stratification diagnostics for the GLM4-AED3 lake outputs.

level_grid      : the saved depth levels (fixed count per lake, fine at the surface).
thermo_depth,
meta_depths     : line-by-line ports of rLakeAnalyzer 1.11.4.1 thermo.depth / meta.depths
                  (the functions ISIMIP3 LER used for `thermodepth`, via ts.thermo.depth),
                  including their quirks (findPeaks' first element, approx() tie handling).
                  Checked against R on GLM profiles: validate_thermocline.py.
strat_flags     : the stratification definitions in use (isi_cal's three from the manuscript
                  discussion, ISIMIP4T's 0.005 kg m-3, ISIMIP3's) -- see its docstring.
"""
import numpy as np

N_LEVELS = 20      # saved levels per lake, surface (0) and bottom (lake depth) included
DZ_TOP = 0.5       # m, spacing of the first interval below the surface


def level_grid(depth, n=N_LEVELS, dz0=DZ_TOP):
    """n levels from 0 to `depth`. Spacing grows geometrically from dz0 at the surface (ratio r
    solved so the n-1 intervals add up to the depth): fine where mixed layers and thermoclines
    sit, coarse in the hypolimnion. Lakes shallower than dz0 * (n-1) get uniform spacing."""
    if depth <= dz0 * (n - 1):
        return np.linspace(0.0, depth, n)
    lo, hi = 1.0 + 1e-12, 10.0
    for _ in range(200):                         # bisection on the geometric sum
        r = (lo + hi) / 2
        if dz0 * (r ** (n - 1) - 1) / (r - 1) > depth:
            hi = r
        else:
            lo = r
    z = np.concatenate([[0.0], np.cumsum(dz0 * r ** np.arange(n - 1))])
    z[-1] = depth
    return z


def level_bounds(z):
    """ISIMIP depth_bnds: half-way between levels, half-thickness top and bottom layers."""
    mid = (z[1:] + z[:-1]) / 2
    return np.column_stack([np.concatenate([[z[0]], mid]), np.concatenate([mid, [z[-1]]])])


def glm_layer_depths(ztop, ns):
    """GLM layers -> (depth of each layer mid-point below the actual surface, increasing from
    surface to bottom; the matching layer indices; the lake level). ztop = layer-top heights
    above the bottom (GLM `z`), ns = layers in use (GLM `NS`)."""
    top = ztop[:ns]
    mid = (top + np.concatenate([[0.0], top[:-1]])) / 2
    level = top[-1]
    return (level - mid)[::-1], np.arange(ns)[::-1], level


def density_mm(t):
    """rLakeAnalyzer water.density for salinity 0 (Martin & McCutcheon)."""
    return 1000 * (1 - (t + 288.9414) * (t - 3.9863) ** 2 / (508929.2 * (t + 68.12963)))


def density_unesco(t):
    """ISIMIP3_LER temp_to_dens.R (Millero & Poisson 1981)."""
    return (999.842594 + 6.793952e-2 * t - 9.095290e-3 * t ** 2 + 1.001685e-4 * t ** 3
            - 1.120083e-6 * t ** 4 + 6.536336e-9 * t ** 5)


def density_chen(t):
    """marelac 2.1.11 sw_dens(S = 0, t, method = "Chen") at surface pressure (p = 0, so the
    compressibility term rho0 / (1 - p/K) drops out): Chen & Millero (1986), kg m-3."""
    return 1000 * (0.9998395 + 6.7914e-05 * t - 9.0894e-06 * t * t + 1.0171e-07 * t * t * t
                   - 1.2846e-09 * t ** 4 + 1.1592e-11 * t ** 5 - 5.0125e-14 * t ** 6)


def strat_flags(t_surface, t_bottom):
    """Daily stratification flags (1/0) from the surface and bottom layer temperatures.

    isi_cal definitions (codes/annual_strat_definitions.R, 2026-09-29; proposed by Martin and
    Fabian after the first manuscript draft; density = marelac Chen):
      dens0.1  : rho_bot - rho_surf > 0.1  kg m-3 and T_surf > T_bot
      dens0.01 : rho_bot - rho_surf > 0.01 kg m-3 and T_surf > T_bot
      temp1    : T_surf - T_bot > 1 K
    dens0.005  : |rho_bot - rho_surf| >= 0.005 kg m-3 -- ISIMIP4T_simulations' `strat`
                 (config/variables.yaml; Baerenbold et al. 2026, validated against deepwater
                 oxygen recharge in 21 Swiss lakes), here with Chen density.
    isimip3    : rho_bot - rho_surf >= 0.1 kg m-3, UNESCO density -- the ISIMIP3 LER `strat`
                 (Scripts/5_Create_ISIMIP_netcdfs.R), for comparison with ISIMIP3 outputs.
    dens0.005s : rho_bot - rho_surf >= 0.005 kg m-3 and T_surf > T_bot -- the ISIMIP4T threshold,
                 direct (summer) stratification only: inverse stratification under ice (cold surface
                 over ~4 C water, also denser at the bottom) does not count. Used for the ISIMIP4a
                 local-lake files (user decision 2026-10-23)."""
    ts, tb = np.asarray(t_surface, float), np.asarray(t_bottom, float)
    dd = density_chen(tb) - density_chen(ts)
    dt = ts - tb
    f = lambda c: np.where(np.isfinite(dd), c.astype(float), np.nan)
    return {"dens0.1": f((dd > 0.1) & (dt > 0)),
            "dens0.01": f((dd > 0.01) & (dt > 0)),
            "temp1": f(dt > 1),
            "dens0.005": f(np.abs(dd) >= 0.005),
            "isimip3": f(density_unesco(tb) - density_unesco(ts) >= 0.1),
            "dens0.005s": f((dd >= 0.005) & (dt > 0))}


def _find_peaks(d, thresh):
    """rLakeAnalyzer:::findPeaks. R loops i in 1..n-1 over d[(i-1):(i+1)]; for i = 1 the
    index 0 is dropped, so the window is d[1:2] (a quirk kept here). Returns 0-based indices."""
    n = len(d)
    locs, peaks = [], []
    for i in range(n - 1):                       # R i = 1 .. n-1  ->  python i = 0 .. n-2
        w = d[max(i - 1, 0):i + 2]
        pk = int(np.argmax(w))                   # which.max: first maximum
        if pk == 1:                              # R: pkI == 2 (2nd element of the window)
            locs.append(i); peaks.append(np.max(w))
    return [l for l, p in zip(locs, peaks) if p > thresh]


def _interp_s(dn, up, sdn, sup):
    with np.errstate(divide="ignore", invalid="ignore"):
        return dn * (sdn / (sdn + sup)) + up * (sup / (sdn + sup))


def thermo_depth(wtr, depths, Smin=0.1, seasonal=True, mixed_cutoff=1.0):
    wtr = np.asarray(wtr, float); depths = np.asarray(depths, float)
    if np.isnan(wtr).any():
        return np.nan
    if np.ptp(wtr) < mixed_cutoff:
        return np.nan
    if len(wtr) < 3:
        return np.nan
    rho = density_mm(wtr)
    n = len(depths)
    with np.errstate(divide="ignore", invalid="ignore"):
        d = (rho[1:] - rho[:-1]) / (depths[1:] - depths[:-1])
    ti = int(np.argmax(d))
    m = d[ti]
    tD = depths[ti:ti + 2].mean()
    with np.errstate(divide="ignore", invalid="ignore"):
        if 0 < ti < n - 2:                       # R: thermoInd > 1 && thermoInd < numDepths-1
            sdn = -(depths[ti + 1] - depths[ti]) / (d[ti + 1] - d[ti])
            sup = (depths[ti] - depths[ti - 1]) / (d[ti] - d[ti - 1])
            if not np.isinf(sup) and not np.isinf(sdn):
                tD = _interp_s(depths[ti + 1], depths[ti], sdn, sup)
        locs = _find_peaks(d, max(0.15 * m, Smin))
        if not locs:
            sD, si = tD, ti
        else:
            si = locs[-1]
            if si > ti + 1:
                sD = depths[si:si + 2].mean()
                if 0 < si < n - 2:
                    sdn = -(depths[si + 1] - depths[si]) / (d[si + 1] - d[si])
                    sup = (depths[si] - depths[si - 1]) / (d[si] - d[si - 1])
                    if not np.isinf(sup) and not np.isinf(sdn):
                        sD = _interp_s(depths[si + 1], depths[si], sdn, sup)
            else:
                sD, si = tD, ti
    if sD < tD:                                  # NaN compares False, as in R
        sD = tD
    return sD if seasonal else tD


def _approx(x, y, xout):
    """stats::approx(x, y, xout), rule = 1, ties = mean: NA pairs dropped, sorted by x."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = ~(np.isnan(x) | np.isnan(y))
    x, y = x[ok], y[ok]
    if x.size == 0:
        return np.nan
    ux, inv = np.unique(x, return_inverse=True)
    uy = np.bincount(inv, weights=y) / np.bincount(inv)
    if ux.size == 1:
        return uy[0] if xout == ux[0] else np.nan
    if xout < ux[0] or xout > ux[-1]:
        return np.nan
    return float(np.interp(xout, ux, uy))


def meta_depths(wtr, depths, slope=0.1, seasonal=True, mixed_cutoff=1.0):
    wtr = np.asarray(wtr, float); depths = np.asarray(depths, float)
    if np.isnan(wtr).any():
        return np.nan, np.nan
    if len(wtr) < 3:
        return depths.max(), depths.max()
    o = np.argsort(depths, kind="stable")
    wtr, depths = wtr[o], depths[o]
    tD = thermo_depth(wtr, depths, seasonal=seasonal, mixed_cutoff=mixed_cutoff)
    if np.isnan(tD):
        return np.nan, np.nan
    rho = density_mm(wtr)
    n = len(depths)
    d = (rho[1:] - rho[:-1]) / (depths[1:] - depths[:-1])
    bot, top = depths[-1], depths[0]
    tdep = (depths[1:] + depths[:-1]) / 2
    sdep = np.unique(np.concatenate([tdep, [tD]]))
    nd = len(sdep)
    dz = np.array([_approx(tdep, d, s) for s in sdep])
    k = int(np.nonzero(sdep == tD)[0][0])
    i = k
    for i in range(k, nd):
        if not np.isnan(dz[i]) and dz[i] < slope:
            bot = sdep[i]
            break
    if i - k >= 1 and (not np.isnan(dz[k]) and dz[k] > slope):
        bot = _approx(dz[k:i + 1], sdep[k:i + 1], slope)
    if np.isnan(bot):
        bot = depths[-1]
    for i in range(k, -1, -1):
        if not np.isnan(dz[i]) and dz[i] < slope:
            top = sdep[i]
            break
    if k - i >= 1 and (not np.isnan(dz[k]) and dz[k] > slope):
        top = _approx(dz[i:k + 1], sdep[i:k + 1], slope)
    if np.isnan(top):
        top = depths[i]                          # R: depths[i] (original depths, sic)
    return top, bot
