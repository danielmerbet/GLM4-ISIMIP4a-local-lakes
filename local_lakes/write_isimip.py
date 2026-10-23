#!/usr/bin/env python3
"""ISIMIP4a local-lake NetCDF files from one GLM run (local_lakes/setup_run.py).

    python3 local_lakes/write_isimip.py <run dir> <Lake> <out dir>

One file per variable, ISIMIP3a/4a lakes_local naming:
  glm4_gswp3-era5_obsclim_histsoc_default_<variable>_<lake-site>_daily_1901_2025.nc
Days 1901-01-01..2025-12-31 (the 1881-1900 spin-up and the 2026-01-01 pad day are dropped).
GLM writes each daily value with the time stamp of the END of its day (hours since start; checked:
daily shortwave correlates 0.991 with the forcing of the day before its stamp, 0.888 with the
stamp's own day), so the value stamped D 00:00 is assigned to day D - 1.

Variables (protocol definitions/variable/hydrothermal.yaml, isimip-protocol-4 9195c32):
  watertemp  K    layer averages on 1 m levels (0.5 m if the lake is < 10 m deep), surface to bottom;
                  levlak = level depth, levlak_bnds = half-way to the neighbours (top 0, bottom = lake
                  bottom); exact volume-weighted averages of GLM's layers (global_test/reduce_output.py
                  Integrator)
  surftemp   K    GLM top water layer; bottemp K: GLM bottom layer
  ice        1    ice (blue + white) present; icethick m: blue + white ice (snow not included)
  strat      1    STRAT_DEF of lake_physics.strat_flags on the top and bottom layer (default dens0.005s:
                  rho_bot - rho_surf >= 0.005 kg m-3 and T_surf > T_bot = the ISIMIP4T threshold of
                  Baerenbold et al. 2026, summer stratification only; user decision 2026-10-23)
  thermodepth m   rLakeAnalyzer thermo.depth seasonal=FALSE on the full GLM profile (protocol: depth of
                  the maximum density gradient); fill value when there is no thermocline
  sensheatf, latentheatf W m-2: -GLM daily mean Qh, -Qe (GLM positive into the water, ISIMIP upwards)
  swup  W m-2   swdown - GLM daily mean net shortwave (swdown = forcing ShortWave x sw_factor)
  lwup  W m-2   lwdown - GLM daily mean net longwave (lwdown = forcing LongWave x lw_factor)
  lakeheatf W m-2 swdown - swup + lwdown - lwup - sensheatf - latentheatf (protocol formula; net heat
                  into the lake)
  momf  kg m-1 s-2  rho_air * cd * U^2 with U = daily mean wind x wind_factor, cd from glm3.nml,
                  rho_air from the air temperature and the pressure at the lake elevation; magnitude
                  of the downward momentum flux (daily mean wind, so a lower bound of the mean stress)
"""
import os
import re
import sys
from pathlib import Path

import netCDF4 as nc
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "global_test"))
from lake_physics import glm_layer_depths, strat_flags, thermo_depth  # noqa: E402
from reduce_output import Integrator  # noqa: E402

MODEL, FORCING, CLIM, SOC, SENS = "glm4", "gswp3-era5", "obsclim", "histsoc", "default"
FIRST, LAST = pd.Timestamp("1901-01-01"), pd.Timestamp("2025-12-31")
STRAT_DEF = os.environ.get("STRAT_DEF", "dens0.005s")   # user decision 2026-10-23: summer stratification only
STRAT_TEXT = {"dens0.005s": "rho_bottom - rho_surface >= 0.005 kg m-3 and T_surface > T_bottom (ISIMIP4T threshold, direct/summer stratification only; inverse stratification under ice not counted; Chen density)",
              "dens0.005": "|rho_bottom - rho_surface| >= 0.005 kg m-3 (ISIMIP4T, Baerenbold et al. 2026; Chen density)",
              "dens0.1": "rho_bottom - rho_surface > 0.1 kg m-3 and T_surface > T_bottom (Chen density)",
              "dens0.01": "rho_bottom - rho_surface > 0.01 kg m-3 and T_surface > T_bottom (Chen density)",
              "temp1": "T_surface - T_bottom > 1 K", "isimip3": "rho_bottom - rho_surface >= 0.1 kg m-3 (UNESCO density)"}
FILL = np.float32(1e20)
BLOCK = 500
CONTACT = os.environ.get("ISIMIP_CONTACT", "Daniel Mercado-Bettín <daniel.mercado@csic.es>")
AUTHOR = "Daniel Mercado-Bettín"
INSTITUTION = os.environ.get("ISIMIP_INSTITUTION", (
    "Spanish National Research Council (CSIC), Centre for Advanced Studies of Blanes (CEAB), Spain; "
    "Vrije Universiteit Brussel (VUB), Belgium"))
REPOSITORY = os.environ.get("ISIMIP_REFERENCES", "https://github.com/danielmerbet/GLM4-ISIMIP4a-local-lakes")
ACKNOWLEDGEMENT = (
    "We thank Wim Thiery (Vrije Universiteit Brussel) for access to the VUB Hydra HPC system. "
    "The computational resources and services used in this work were provided by the VSC (Flemish "
    "Supercomputer Center), funded by the Research Foundation - Flanders (FWO) and the Flemish Government.")
LONG = {"watertemp": "Temperature of Lake Water", "surftemp": "Temperature of Lake Surface Water",
        "bottemp": "Temperature of Lake Bottom Water", "ice": "Lake Ice Cover", "icethick": "Ice Thickness",
        "strat": "Thermal Stratification", "thermodepth": "Depth of Thermocline",
        "sensheatf": "Sensible Heat Flux at Lake-Atmosphere Interface", "latentheatf": "Latent Heat Flux at Lake-Atmosphere Interface",
        "momf": "Momentum Flux at Lake-Atmosphere Interface", "swup": "Upward Short-Wave Radiation Flux at Lake-Atmosphere Interface",
        "lwup": "Upward Long-Wave Radiation Flux at Lake-Atmosphere Interface", "lakeheatf": "Downward Heat Flux at Lake-Atmosphere Interface"}
UNITS = {"watertemp": "K", "surftemp": "K", "bottemp": "K", "ice": "1", "icethick": "m", "strat": "1", "thermodepth": "m",
         "sensheatf": "W m-2", "latentheatf": "W m-2", "momf": "kg m-1 s-2", "swup": "W m-2", "lwup": "W m-2", "lakeheatf": "W m-2"}


def nml_value(txt, key):
    return float(re.search(rf"(?m)^\s*{key}\s*=\s*([-0-9.eE+]+)", txt).group(1))


def levels(depth):
    dz = 0.5 if depth < 10 else 1.0
    z = np.arange(0.0, depth, dz)
    z = np.append(z, depth) if depth - z[-1] > 1e-6 else z
    b = np.empty((z.size, 2))
    b[:, 0] = np.r_[0.0, (z[:-1] + z[1:]) / 2]
    b[:, 1] = np.r_[(z[:-1] + z[1:]) / 2, depth]
    return z, b


def main(run, lake, outdir):
    run, outdir = Path(run), Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    site = pd.read_csv(ROOT / "local_lakes" / "lake_sites.csv").set_index("lake").loc[lake]
    nml = (run / "glm3.nml").read_text()
    d = nc.Dataset(run / "output" / "output.nc"); d.set_auto_mask(False)
    stamp = pd.to_datetime(nc.num2date(d["time"][:], d["time"].units, only_use_cftime_datetimes=False,
                                       only_use_python_datetimes=True))
    day = stamp.normalize() - pd.Timedelta(days=1)
    keep = np.where((day >= FIRST) & (day <= LAST))[0]
    days = day[keep]
    exp = pd.date_range(FIRST, LAST, freq="D")
    if len(days) != len(exp) or not (days == exp).all():
        raise SystemExit(f"{lake}: output days {days[0]}..{days[-1]} ({len(days)}) do not cover {FIRST.date()}..{LAST.date()}")
    nt = len(keep)
    ns = d["NS"][:].ravel().astype(int)
    depth = float(site.max_depth_m)
    hy = pd.read_csv(ROOT / "calibration" / "runs" / lake / "hypsograph.csv")
    crest = nml_value(nml, "crest_elev") - float(re.search(r"(?m)^\s*H\s*=\s*([-0-9.eE+]+)", nml).group(1))
    # levels down to the deepest point below the spill level (MtBold: hypsograph 45.42 m, crest 44.5 m;
    # the water never reaches the top 0.92 m, so a grid on 45.42 m kept one level always empty)
    depth = min(float(hy.Depth_meter.max()), crest)
    z_lev, bnds = levels(depth)
    wt = np.full((nt, z_lev.size), np.nan, "f4")
    ts = np.full(nt, np.nan); tb = np.full(nt, np.nan); td = np.full(nt, np.nan)
    for b0 in range(0, nt, BLOCK):
        idx = keep[b0:b0 + BLOCK]; n = len(idx)
        blk = lambda v: d[v][idx[0]:idx[-1] + 1, :, 0, 0].astype("f8")
        z, V, temp = blk("z"), blk("V"), blk("temp"); nb = ns[idx]
        for j in range(n):
            dd, li, _ = glm_layer_depths(z[j], nb[j])
            w = temp[j, li]
            ts[b0 + j], tb[b0 + j] = w[0], w[-1]
            td[b0 + j] = thermo_depth(w, dd, seasonal=False)
        integ = Integrator(z, V, nb)
        lower = np.tile(bnds[:, 1], (n, 1)); lower[:, -1] = np.maximum(lower[:, -1], integ.level)
        wt[b0:b0 + n] = integ.average(temp, integ.prepare(np.tile(bnds[:, 0], (n, 1)), lower))
    one = lambda v: d[v][:].ravel()[keep].astype("f8")
    ice_m = one("blue_ice_thickness") + one("white_ice_thickness")
    met = pd.read_csv(run / "meteo_file.csv", parse_dates=["Date"], index_col="Date")
    met.index = met.index.normalize(); met = met.loc[days]
    swdown = met.ShortWave.values * nml_value(nml, "sw_factor")
    lwdown = met.LongWave.values * nml_value(nml, "lw_factor")
    qsw, qlw, qh, qe = one("daily_qsw"), one("daily_qlw"), one("daily_qh"), one("daily_qe")
    p = 101325.0 * (1 - 2.25577e-5 * max(float(site.elev_m), 0.0)) ** 5.25588
    rho_air = p / (287.05 * (met.AirTemp.values + 273.15))
    u = met.WindSpeed.values * nml_value(nml, "wind_factor")
    out = {"watertemp": wt + 273.15, "surftemp": ts + 273.15, "bottemp": tb + 273.15,
           "ice": (ice_m > 0).astype(float), "icethick": ice_m, "strat": strat_flags(ts, tb)[STRAT_DEF],
           "thermodepth": td, "sensheatf": -qh, "latentheatf": -qe, "swup": swdown - qsw, "lwup": lwdown - qlw,
           "lakeheatf": qsw + qlw + qh + qe, "momf": rho_air * nml_value(nml, "cd") * u ** 2}
    lat = float(site.calib_lat) if lake == "Bosumtwi" else float(site.site_lat)   # protocol entry has a sign error
    lon = float(site.calib_lon) if lake == "Bosumtwi" else float(site.site_lon)
    tvals = np.arange(nt, dtype="f8")
    written = []
    for var, val in out.items():
        f = outdir / f"{MODEL}_{FORCING}_{CLIM}_{SOC}_{SENS}_{var}_{site.lake_site}_daily_{FIRST.year}_{LAST.year}.nc"
        o = nc.Dataset(f, "w", format="NETCDF4_CLASSIC")
        o.createDimension("time", None)
        prof = val.ndim == 2
        if prof:
            o.createDimension("levlak", z_lev.size); o.createDimension("bnds", 2)
        o.createDimension("lat", 1); o.createDimension("lon", 1)
        # explicit chunks: with the unlimited time dimension netCDF would otherwise store one day per
        # chunk (45 MB for a 7 m lake, slow to read)
        tv = o.createVariable("time", "f8", ("time",), chunksizes=(nt,)); tv.standard_name = "time"; tv.long_name = "Time"
        tv.units = "days since 1901-01-01 00:00:00"; tv.calendar = "proleptic_gregorian"; tv.axis = "T"; tv[:] = tvals
        if prof:
            lv = o.createVariable("levlak", "f8", ("levlak",)); lv.long_name = "Depth of vertical layer below the lake surface"
            lv.units = "m"; lv.positive = "down"; lv.axis = "Z"; lv.bounds = "levlak_bnds"; lv[:] = z_lev
            o.createVariable("levlak_bnds", "f8", ("levlak", "bnds"))[:] = bnds
        la = o.createVariable("lat", "f8", ("lat",)); la.standard_name = "latitude"; la.long_name = "Latitude"
        la.units = "degrees_north"; la.axis = "Y"; la[:] = lat
        lo = o.createVariable("lon", "f8", ("lon",)); lo.standard_name = "longitude"; lo.long_name = "Longitude"
        lo.units = "degrees_east"; lo.axis = "X"; lo[:] = lon
        dims = ("time", "levlak", "lat", "lon") if prof else ("time", "lat", "lon")
        chunks = (min(nt, 3653), z_lev.size, 1, 1) if prof else (nt, 1, 1)
        vv = o.createVariable(var, "f4", dims, zlib=True, complevel=5, fill_value=FILL, chunksizes=chunks)
        vv.long_name = LONG[var]; vv.units = UNITS[var]; vv.missing_value = FILL
        if var == "strat":
            vv.comment = "1 if stratified, 0 if not: " + STRAT_TEXT[STRAT_DEF] + ", GLM top and bottom layer"
        if var == "watertemp":
            vv.comment = ("layer averages over levlak_bnds (exact volume-weighted averages of the GLM layers); "
                          f"{'0.5' if depth < 10 else '1'} m levels from the surface to the bottom ({depth:g} m)")
        if var == "momf":
            vv.comment = "rho_air * C_D * U^2 from the daily mean wind (x wind_factor); magnitude of the downward momentum flux"
        x = np.where(np.isfinite(val), val, FILL).astype("f4")
        vv[:] = x[:, :, None, None] if prof else x[:, None, None]
        o.contact = CONTACT; o.institution = INSTITUTION; o.author = AUTHOR
        o.references = REPOSITORY + " (model setup, calibration, code and documentation)"
        o.acknowledgement = ACKNOWLEDGEMENT
        o.license = "CC0 1.0 Universal (public domain dedication), https://creativecommons.org/publicdomain/zero/1.0/"
        o.model = "GLM 4.0.0 (glm-aed 83e20ed) with local fixes (patches/ in the repository above), physics only"
        o.simulation_round = "ISIMIP4a"; o.sector = "lakes_local"; o.lake_site = site.lake_site
        o.comment = (f"{lake}: own calibrated parameters (calibration/results/{lake}/best.csv); spin-up 1881-1900 = "
                     "forcing 1901-1920 repeated; closed lake with evaporation-compensating inflow")
        o.close(); written.append(f.name)
    d.close()
    print(f"{lake}: {len(written)} files, {nt} days, {z_lev.size} levels ({depth:g} m) -> {outdir}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
