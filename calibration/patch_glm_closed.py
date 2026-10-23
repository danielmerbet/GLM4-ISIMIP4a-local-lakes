#!/usr/bin/env python3
"""Post-export patches for the calibration lakes (closed lakes, no in/outflow data).

    python3 calibration/patch_glm_closed.py <Lake> [...]      # or: all

Calls ISIMIP4_LocalLakes_LERWQ's own patch functions (lerwq/patches.py), so the fix is
identical to the global pipeline's:
  * glm_drop_outflows: LER keeps the template's outflow file entries for closed lakes;
    replace with num_outlet = 0 (excess water leaves over the crest). MODEL_ISSUES #19.
  * glm_evap_compensation: GLM has no fixed-level switch, so a closed lake drains by
    evaporation (Allequash 8 m -> 6.4 m in 10 y in the first test). Adds a solute-free
    inflow of 1.25 x estimated evaporation; the surplus spills over the crest, pinning the
    level. MODEL_ISSUES #9. Surface area = top of the hypsograph.
Idempotent per freshly exported folder (setup_lake_glm.R re-exports from scratch).
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from closed_lake_patches import glm_drop_outflows, glm_evap_compensation  # noqa: E402  (copied from LERWQ)

HERE = Path(__file__).resolve().parent


TROPICS_LAT = 23.44     # same rule as global_test/patch_global.py
HIGH_ELEV_M = 1500.0


def set_stability(lake):
    """atm_stab = 1 (stability correction) for |lat| <= 23.44 OR elevation > 1500 m, 0 elsewhere
    (user decision 2026-10-22, analysis/stability/); written explicitly every time (idempotent), so
    the calibration uses the physics the global run will use."""
    import re
    r = pd.read_csv(HERE / "lakes.csv").set_index("lake").loc[lake]
    stab = 1 if (abs(float(r.lat)) <= TROPICS_LAT or float(r.elev_m) > HIGH_ELEV_M) else 0
    f = HERE / "runs" / lake / "GLM" / "glm3.nml"
    txt, n = re.subn(r"(?m)^(\s*atm_stab\s*=\s*)\d+", rf"\g<1>{stab}", f.read_text())
    if n != 1:
        raise SystemExit(f"{lake}: expected one atm_stab entry, found {n}")
    f.write_text(txt)
    return stab


def patch(lake):
    d = HERE / "runs" / lake / "GLM"
    stab = set_stability(lake)
    if (d / "inflow_evapcomp.csv").exists():
        return f"{lake}: already patched (atm_stab={stab})"
    area = pd.read_csv(HERE / "runs" / lake / "hypsograph.csv").sort_values("Depth_meter").Area_meterSquared.iloc[0]
    dropped = glm_drop_outflows(d)
    q = glm_evap_compensation(d, area)
    return f"{lake}: atm_stab={stab}, outflows dropped={dropped}, evap-compensation inflow mean {q:.4g} m3/s (area {area:.4g} m2)"


if __name__ == "__main__":
    lakes = sys.argv[1:]
    if lakes == ["all"]:
        # only lakes whose setup completed (window.csv is written last)
        lakes = sorted(p.name for p in (HERE / "runs").iterdir() if (p / "window.csv").exists())
    for lk in lakes:
        try:
            print(patch(lk))
        except Exception as e:
            print(f"{lk}: PATCH FAILED: {e}")
