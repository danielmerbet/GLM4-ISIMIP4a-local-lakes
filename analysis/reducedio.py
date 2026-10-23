"""Reading reduced.nc (global_test/reduce_output.py) without plotting dependencies, so the same
code runs on Hydra's Python module (no matplotlib) and in the plotting scripts (plotlib re-exports it)."""
import netCDF4 as nc
import numpy as np
import pandas as pd

MOL = {"P": 30.974, "N": 14.007, "O2": 32.0}
TP = ["PHS_frp", "OGM_dop", "OGM_dopr", "OGM_pop", "PHY_ip"]
TN = ["NIT_nit", "NIT_amm", "OGM_don", "OGM_donr", "OGM_pon", "PHY_in"]


def load(path):
    """reduced.nc -> (time index, depth grid, dict of arrays with NaN fills, setup attrs)."""
    d = nc.Dataset(path); d.set_auto_mask(False)
    t = pd.to_datetime(nc.num2date(d["time"][:], d["time"].units, only_use_cftime_datetimes=False,
                                   only_use_python_datetimes=True)).normalize()
    x = {v: np.where(d[v][:] > 1e19, np.nan, d[v][:]).astype("f8") for v in d.variables if v not in ("time",)}
    s = {k[6:]: d.getncattr(k) for k in d.ncattrs() if k.startswith("setup_")}
    depth = d["depth"][:].astype("f8")
    d.close()
    return t, depth, x, s
