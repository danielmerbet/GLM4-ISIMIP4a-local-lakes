"""Closed-lake GLM fixes, copied verbatim from ISIMIP4_LocalLakes_LERWQ (lerwq/patches.py) so this
repository is self-contained. Used by patch_glm_closed.py.
"""
import re
from pathlib import Path


def glm_drop_outflows(glm_dir):
    """Closed lakes: LER sets num_outlet = 0 but keeps the template's outflow_fl entries, whose
    outflow_1/2.csv are never written -> LER-WQ skips GLM-AED2. Excess water leaves over the crest."""
    f = Path(glm_dir) / "glm3.nml"
    nml = f.read_text()
    m = re.search(r"&outflow\b.*?\n/", nml, flags=re.S)
    fl = re.search(r"outflow_fl\s*=\s*(.+)", m.group(0)) if m else None
    if not fl:
        return False
    missing = [n for n in re.findall(r"'([^']+)'", fl.group(1)) if not (Path(glm_dir) / n).exists()]
    if not missing:
        return False
    f.write_text(nml.replace(m.group(0), "&outflow\n   num_outlet = 0\n/"))
    return True



def glm_evap_compensation(glm_dir, surface_area_m2, overshoot=1.25, ce=0.0013):
    """Constant level for GLM (no fixed-level switch; evap_file not supported in this build).

    Adds an inflow of solute-free water = overshoot x estimated evaporation, at the estimated
    surface temperature. GLM spills everything above the crest, so the level is pinned at the
    crest instead of drifting down in dry spells (Erken: 20.58-21.06 m without this).
    Solute-free water matches GOTM's fixed-level runs, where evaporation does not concentrate
    solutes either. Evaporation: Dalton bulk formula with GLM's own ce; water temperature ~
    15-day running mean air temperature, >= 0."""
    import numpy as np
    import pandas as pd
    d = Path(glm_dir)
    nml = (d / "glm3.nml").read_text()
    met = pd.read_csv(d / re.search(r"meteo_fl\s*=\s*'([^']+)'", nml).group(1))
    t = pd.to_datetime(met.iloc[:, 0])
    ta, rh, u = met["AirTemp"].values, met["RelHum"].values, met["WindSpeed"].values
    tw = pd.Series(ta).rolling(15, min_periods=1).mean().clip(lower=0).values

    def qsat(T):        # kg/kg at ~1013 hPa (Magnus)
        es = 6.112 * np.exp(17.67 * T / (T + 243.5))
        return 0.622 * es / 1013.25
    rho_a, rho_w = 1.2, 1000.0
    e = ce * rho_a * u * np.clip(qsat(tw) - rh / 100 * qsat(ta), 0, None) / rho_w      # m/s
    flow = overshoot * e * surface_area_m2                                             # m3/s
    m = re.search(r"&inflow\b.*?\n/", nml, flags=re.S)
    block = m.group(0) if m else "&inflow\n   num_inflows = 0\n/"
    n0 = int(re.search(r"num_inflows\s*=\s*(\d+)", block).group(1))
    vars_ = re.search(r"inflow_vars\s*=\s*(.+)", block)
    names = [v.strip().strip("'") for v in vars_.group(1).split(",")] if vars_ else ["FLOW", "SALT", "TEMP"]
    comp = pd.DataFrame({"Time": t.dt.strftime("%Y-%m-%d 00:00:00")})
    for v in names:
        comp[v] = 0.0
    comp["FLOW"], comp["TEMP"] = flow.round(6), tw.round(3)
    comp.to_csv(d / "inflow_evapcomp.csv", index=False)

    def lst(key, new, default):
        mm = re.search(rf"{key}\s*=\s*(.+)", block)
        vals = [x.strip() for x in mm.group(1).split(",")] if mm else []
        vals = (vals + [vals[-1] if vals else default] * n0)[:n0] if n0 else []
        return ", ".join(vals + [new])
    new = ("&inflow\n"
           f"   num_inflows = {n0 + 1}\n"
           f"   names_of_strms = {lst('names_of_strms', chr(39) + 'evapcomp' + chr(39), chr(39) + 'in' + chr(39))}\n"
           f"   subm_flag = {lst('subm_flag', '.false.', '.false.')}\n"
           f"   strm_hf_angle = {lst('strm_hf_angle', '65', '65')}\n"
           f"   strmbd_slope = {lst('strmbd_slope', '2', '2')}\n"
           f"   strmbd_drag = {lst('strmbd_drag', '0.016', '0.016')}\n"
           f"   inflow_factor = {lst('inflow_factor', '1', '1')}\n"
           f"   inflow_fl = {lst('inflow_fl', chr(39) + 'inflow_evapcomp.csv' + chr(39), '')}\n"
           f"   inflow_varnum = {len(names)}\n"
           f"   inflow_vars = {', '.join(repr(v) for v in names)}\n/")
    nml = nml[:m.start()] + new + nml[m.end():] if m else nml + "\n" + new + "\n"
    (d / "glm3.nml").write_text(nml)
    return float(flow.mean())

