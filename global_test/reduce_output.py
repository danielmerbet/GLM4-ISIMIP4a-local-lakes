#!/usr/bin/env python3
"""Reduce one GLM4-AED3 output.nc (~3.6 GB: 226 variables x max_layers, daily) to a small
NetCDF -- the input for the ISIMIP files and the checks.

    python3 global_test/reduce_output.py <GLM/output/output.nc> <setup.csv> <out.nc>

ISIMIP protocol (isimip-protocol-3 adcd97f / -4 9195c32, lake variables identical):
watertemp = "full profile of simulated water temperature layer averages" (levlak);
water-quality variables = "either full profile, or mean epi and mean hypo"; thermodepth =
"depth corresponding [to] the maximum water density gradient"; strat = 1/0, no definition given.
What is written here:

Profiles (*all* PROFILES), dims (time, depth), on N_LEVELS = 20 levels per lake
(lake_physics.level_grid: surface and bottom included, 0.5 m first spacing growing
geometrically; uniform in lakes < 9.5 m). Each value is the VOLUME-WEIGHTED LAYER AVERAGE over
the level's depth_bnds (half-way between levels; half-thickness top and bottom), integrated
exactly over GLM's Lagrangian layers using GLM's own layer volumes `V` (partial layers by
height; area taken as uniform within a GLM layer). Parts of a level below the actual bottom
(level dropped) are left out; a level entirely below it is NaN. The bottom level always extends
to the actual bottom (when rain lifts the level above the nominal depth), so the levels cover
the whole water column every day.

Epilimnion / hypolimnion means (<var>_epi, <var>_hypo), dims (time): volume-weighted averages
over 0..metatop and metabot..bottom, same integration; on days without a thermocline
(rLakeAnalyzer: top-bottom < 1 C) both are the whole-lake mean.

Thermocline, from the FULL GLM profile (every layer mid-point), daily:
  thermodepth          : rLakeAnalyzer thermo.depth, seasonal = FALSE -- the protocol's
                         definition (depth of the maximum density gradient, interpolated)
  thermodepth_seasonal : rLakeAnalyzer thermo.depth, seasonal = TRUE -- what ISIMIP3 LER used
  metatop / metabot    : rLakeAnalyzer meta.depths around `thermodepth` (seasonal = FALSE)
  NaN = no thermocline (ISIMIP3 LER wrote 0; apply such a convention when writing ISIMIP files).
Stratification flags strat_* (lake_physics.strat_flags), surface vs bottom GLM layer:
  strat_dens0p1, strat_dens0p01, strat_temp1 (isi_cal / manuscript), strat_dens0p005 (ISIMIP4T),
  strat_isimip3.
Scalars: GLM's own surface/bottom-layer temperatures (surftemp/bottemp = "upper/lowest layer"),
ice, level, fluxes (SCALARS).
"""
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from lake_physics import (N_LEVELS, DZ_TOP, glm_layer_depths, level_bounds, level_grid,  # noqa: E402
                          meta_depths, strat_flags, thermo_depth)

# raw GLM/AED states (units as in GLM: degC, mmol m-3, ug/L); ISIMIP variables are combinations
# of these, see ISIMIP_MAP
PROFILES = ["temp", "OXY_oxy", "PHY_tchla", "PHY_cyano", "PHY_diatom", "ZOO_daphnia", "PHY_in", "PHY_ip",
            "NIT_amm", "NIT_nit", "PHS_frp", "OGM_doc", "OGM_docr", "OGM_poc", "OGM_don", "OGM_donr",
            "OGM_pon", "OGM_dop", "OGM_dopr", "OGM_pop", "SIL_rsi", "CAR_dic", "CAR_pH", "extc"]
# AED's own TOT_tn / TOT_tp are NOT used: in the Mendota test they differ from the sum of their
# components in some layer on 15,953 of 16,435 days (validate_reduction.py, check 4).
ISIMIP_MAP = {   # documentation for the ISIMIP writer (mmol -> mol: x 1e-3; ug/L -> g m-3: x 1e-3)
    "watertemp": "temp + 273.15",
    "chl": "PHY_tchla * 1e-3",
    "phytobio": "PHY_diatom, PHY_cyano (one per functional group) * 1e-3",
    "tp": "(PHS_frp + OGM_dop + OGM_dopr + OGM_pop + PHY_ip) * 1e-3",
    "pp": "(OGM_pop + PHY_ip) * 1e-3",
    "tpd": "(PHS_frp + OGM_dop + OGM_dopr) * 1e-3",
    "tn": "(NIT_nit + NIT_amm + OGM_don + OGM_donr + OGM_pon + PHY_in) * 1e-3",
    "pn": "(OGM_pon + PHY_in) * 1e-3",
    "tdn": "(NIT_nit + NIT_amm + OGM_don + OGM_donr) * 1e-3",
    "do": "OXY_oxy * 1e-3",
    "doc": "(OGM_doc + OGM_docr) * 1e-3",
    "si": "SIL_rsi * 1e-3",
    "zoobio": "ZOO_daphnia * 1e-3"}
SCALARS = ["surface_temp", "lake_level", "lake_volume", "blue_ice_thickness", "white_ice_thickness",
           "snow_thickness", "evaporation", "rain", "overflow_vol", "tot_inflow_vol", "albedo",
           "daily_qsw", "daily_qlw", "daily_qe", "daily_qh"]
STRAT_NAMES = {"dens0.1": "strat_dens0p1", "dens0.01": "strat_dens0p01", "temp1": "strat_temp1",
               "dens0.005": "strat_dens0p005", "isimip3": "strat_isimip3"}
STRAT_DOC = {
    "dens0.1": "isi_cal: rho_bot - rho_surf > 0.1 kg m-3 and T_surf > T_bot (Chen density)",
    "dens0.01": "isi_cal: rho_bot - rho_surf > 0.01 kg m-3 and T_surf > T_bot (Chen density)",
    "temp1": "isi_cal: T_surf - T_bot > 1 K",
    "dens0.005": "ISIMIP4T: |rho_bot - rho_surf| >= 0.005 kg m-3 (Baerenbold et al. 2026; Chen density)",
    "isimip3": "ISIMIP3 LER strat: rho_bot - rho_surf >= 0.1 kg m-3 (UNESCO density)"}
F4 = dict(zlib=True, complevel=5, fill_value=np.float32(1e20))
BLOCK = 500     # days per block (memory ~ BLOCK x GLM layers x 20 levels)


class Integrator:
    """Exact volume-weighted averages over depth ranges [d0, d1] (m below the actual surface)
    for every time step, from GLM's Lagrangian layers. Works in height above the bottom:
    GLM layer i spans z[i-1]..z[i] (z[-1] = 0) and holds V[i]; within a layer the volume is
    taken as linear in height. Cumulative volume/content at a height h = sum of the layers
    below + the fraction of the layer containing h."""

    def __init__(self, z, V, ns):
        nt, nl = z.shape
        self.used = np.arange(nl)[None, :] < ns[:, None]
        self.z = np.where(self.used, z, np.inf)                 # unused layers never "below" h
        self.zb = np.concatenate([np.zeros((nt, 1)), np.where(self.used, z, 0)[:, :-1]], axis=1)
        self.V = np.where(self.used, V, 0.0)
        self.level = z[np.arange(nt), ns - 1]
        self.ns = ns
        self.cumV = np.concatenate([np.zeros((nt, 1)), np.cumsum(self.V, axis=1)], axis=1)

    def _locate(self, depths):
        """depths (nt, k) -> heights clipped to the water column, layer index and fraction."""
        h = np.clip(self.level[:, None] - depths, 0.0, self.level[:, None])
        i = (self.z[:, :, None] <= h[:, None, :]).sum(axis=1)  # first layer with top > h
        i = np.minimum(i, self.ns[:, None] - 1)
        rows = np.arange(h.shape[0])[:, None]
        zb, zt = self.zb[rows, i], self.z[rows, i]
        f = np.clip((h - zb) / np.where(zt > zb, zt - zb, 1.0), 0.0, 1.0)
        return i, f

    def prepare(self, d0, d1):
        """Precompute for ranges d0 < d1 (arrays (nt, k)); returns a handle for average()."""
        return self._locate(d0), self._locate(d1)

    def _cum(self, cum, x_v, loc):
        i, f = loc
        rows = np.arange(i.shape[0])[:, None]
        return cum[rows, i] + f * x_v[rows, i]

    def average(self, x, handle):
        """x (nt, nl) per-layer values -> volume-weighted averages (nt, k)."""
        xv = np.where(self.used, x, 0.0) * self.V
        cumx = np.concatenate([np.zeros((x.shape[0], 1)), np.cumsum(xv, axis=1)], axis=1)
        top, bot = handle                                      # d0 = shallower = higher
        num = self._cum(cumx, xv, top) - self._cum(cumx, xv, bot)
        den = self._cum(self.cumV, self.V, top) - self._cum(self.cumV, self.V, bot)
        with np.errstate(invalid="ignore", divide="ignore"):
            return np.where(den > 0, num / den, np.nan)


def main(src, setup_csv, out):
    s = pd.read_csv(setup_csv).iloc[0]
    d = nc.Dataset(src)
    d.set_auto_mask(False)
    ns = d["NS"][:].ravel().astype(int)
    nt = ns.size
    grid = level_grid(float(s.max_depth))
    bnds = level_bounds(grid)

    # outputs, filled block by block (memory bounded by BLOCK x layers, whatever the lake depth;
    # every variable of a block is read before the next block, so the file is read in one pass)
    prof = {v: np.full((nt, grid.size), np.nan, "f4") for v in PROFILES}
    epi_m = {v: np.full(nt, np.nan, "f4") for v in PROFILES}
    hyp_m = {v: np.full(nt, np.nan, "f4") for v in PROFILES}
    thermo = np.full(nt, np.nan); thermo_s = np.full(nt, np.nan)
    mtop = np.full(nt, np.nan); mbot = np.full(nt, np.nan)
    t_surf = np.empty(nt); t_bot = np.empty(nt)
    for v in PROFILES:
        if v not in d.variables:
            raise SystemExit(f"{v} not in {src}")
    for t0 in range(0, nt, BLOCK):
        t1 = min(t0 + BLOCK, nt)
        blk = lambda v: d[v][t0:t1, :, 0, 0].astype("f8")
        z, V, temp = blk("z"), blk("V"), blk("temp")
        nb = ns[t0:t1]
        # thermocline, metalimnion, surface/bottom layer temperature: full GLM profile, daily
        for j in range(t1 - t0):
            t = t0 + j
            dd, idx, _ = glm_layer_depths(z[j], nb[j])
            wt = temp[j, idx]
            thermo[t] = thermo_depth(wt, dd, seasonal=False)
            thermo_s[t] = thermo_depth(wt, dd, seasonal=True)
            mtop[t], mbot[t] = meta_depths(wt, dd, seasonal=False)
            t_surf[t], t_bot[t] = wt[0], wt[-1]
        integ = Integrator(z, V, nb)
        n = t1 - t0
        # the bottom level always reaches the actual bottom (level above nominal after rain)
        lower = np.tile(bnds[:, 1], (n, 1))
        lower[:, -1] = np.maximum(lower[:, -1], integ.level)
        lev = integ.prepare(np.tile(bnds[:, 0], (n, 1)), lower)
        # epilimnion 0..metatop, hypolimnion metabot..bottom; no thermocline -> both whole lake
        has_tc = np.isfinite(mtop[t0:t1]) & np.isfinite(mbot[t0:t1])
        zero = np.zeros(n)
        epi = integ.prepare(zero[:, None], np.where(has_tc, mtop[t0:t1], integ.level)[:, None])
        hyp = integ.prepare(np.where(has_tc, mbot[t0:t1], zero)[:, None], integ.level[:, None])
        for v in PROFILES:
            x = temp if v == "temp" else blk(v)
            prof[v][t0:t1] = integ.average(x, lev)
            epi_m[v][t0:t1] = integ.average(x, epi)[:, 0]
            hyp_m[v][t0:t1] = integ.average(x, hyp)[:, 0]

    o = nc.Dataset(out, "w", format="NETCDF4_CLASSIC")
    o.createDimension("time", None); o.createDimension("depth", grid.size); o.createDimension("bnds", 2)
    tv = o.createVariable("time", "f8", ("time",)); tv.units = d["time"].units
    tv.calendar = "proleptic_gregorian"; tv[:] = d["time"][:]
    dv = o.createVariable("depth", "f4", ("depth",)); dv.units = "m"; dv.positive = "down"
    dv.bounds = "depth_bnds"; dv.standard_name = "depth_below_surface"
    dv.comment = (f"{N_LEVELS} levels, first spacing {DZ_TOP} m, geometric growth to the bottom; "
                  "values are volume-weighted averages over depth_bnds")
    dv[:] = grid
    o.createVariable("depth_bnds", "f4", ("depth", "bnds"))[:] = bnds

    def put(name, dims, val, units, long_name):
        ov = o.createVariable(name, "f4", dims, **F4)
        ov.units = units; ov.long_name = long_name
        ov[:] = np.where(np.isfinite(val), val, 1e20)

    for v in PROFILES:
        u, ln = getattr(d[v], "units", ""), getattr(d[v], "long_name", v)
        put(v, ("time", "depth"), prof[v], u, f"{ln}, layer average over depth_bnds")
        put(f"{v}_epi", ("time",), epi_m[v], u, f"{ln}, epilimnion mean (0..metatop; whole lake if mixed)")
        put(f"{v}_hypo", ("time",), hyp_m[v], u, f"{ln}, hypolimnion mean (metabot..bottom; whole lake if mixed)")
    for v in SCALARS:
        ov = o.createVariable(v, "f4", ("time",), **F4)
        ov.units = getattr(d[v], "units", ""); ov[:] = d[v][:].ravel()
    put("temp_surface_layer", ("time",), t_surf, "celsius", "GLM top layer temperature (protocol surftemp)")
    put("temp_bottom_layer", ("time",), t_bot, "celsius", "GLM bottom layer temperature (protocol bottemp)")
    put("thermodepth", ("time",), thermo, "m", "thermocline depth: rLakeAnalyzer thermo.depth seasonal=FALSE "
        "(protocol: depth of the maximum density gradient), full GLM profile; NaN = no thermocline")
    put("thermodepth_seasonal", ("time",), thermo_s, "m", "rLakeAnalyzer thermo.depth seasonal=TRUE (ISIMIP3 LER), full GLM profile")
    put("metatop", ("time",), mtop, "m", "metalimnion top, rLakeAnalyzer meta.depths seasonal=FALSE, full GLM profile")
    put("metabot", ("time",), mbot, "m", "metalimnion bottom, rLakeAnalyzer meta.depths seasonal=FALSE, full GLM profile")
    flags = strat_flags(t_surf, t_bot)
    for k, name in STRAT_NAMES.items():
        put(name, ("time",), flags[k], "1", "stratified (surface vs bottom GLM layer), " + STRAT_DOC[k])
    for k, val in s.items():
        o.setncattr(f"setup_{k}", str(val))
    o.isimip_map = "; ".join(f"{k} = {v}" for k, v in ISIMIP_MAP.items())
    o.source = f"GLM 4.0.0 / AED 3.0.0, reduced from {src} by global_test/reduce_output.py"
    o.close(); d.close()


if __name__ == "__main__":
    main(*sys.argv[1:4])
