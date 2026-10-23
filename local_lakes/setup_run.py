#!/usr/bin/env python3
"""GLM run folder for one ISIMIP4a local lake (obsclim, histsoc, 1901-2025) with its own calibrated
parameters.

    python3 local_lakes/setup_run.py <Lake> <run dir>

Source: calibration/runs/<Lake>/GLM (calibration setup: hypsograph, ISIMIP4a GSWP3-ERA5 obsclim
forcing 1901-2025 at the lake point, initial profile, closed lake with the evaporation-compensating
inflow, atm_stab rule) and calibration/results/<Lake>/best.csv (its calibrated parameters).
Changes:
  * parameters: wind_factor, sw_factor, Kw (= kw_base x kw_scaling, best.csv Kw), coef_mix_hyp,
    coef_mix_conv, coef_mix_turb from best.csv (as calibration/glm4_objective.R sets them);
  * spin-up (user decision 2026-10-23): the first 20 years of forcing (1901-1920) repeated in front
    of 1901, dated 1881-1900 (date - 20 years; leap years line up, except 1920-02-29, which has no
    1900 counterpart and is dropped from the copy). Same for the inflow file;
  * one padding day 2026-01-01 (copy of 2025-12-31) so that GLM completes 2025-12-31;
  * run 1881-01-01 .. 2026-01-01, daily output. write_isimip.py keeps 1901-01-01 .. 2025-12-31.
"""
import re
import shutil
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SPIN_FROM, SPIN_TO, SHIFT_Y = "1901-01-01", "1920-12-31", 20
START, END, PAD = "1881-01-01", "2025-12-31", "2026-01-01"
PARS = ["wind_factor", "sw_factor", "Kw", "coef_mix_hyp", "coef_mix_conv", "coef_mix_turb"]


def set_nml(txt, key, val, fmt="{:.10g}"):
    pat = rf"(?m)^(\s*{key}\s*=\s*)[^!,/\n]*"
    if len(re.findall(pat, txt)) != 1:
        raise SystemExit(f"expected exactly one '{key}' in glm3.nml")
    return re.sub(pat, lambda m: m.group(1) + (fmt.format(val) if not isinstance(val, str) else val), txt)


def extend(csv, time_col):
    """Prepend the 1901-1920 rows dated -20 years (without 1920-02-29) and append one pad day."""
    d = pd.read_csv(csv)
    t = pd.to_datetime(d[time_col])
    if t.iloc[0] > pd.Timestamp(SPIN_FROM) or t.iloc[-1] < pd.Timestamp(END):
        raise SystemExit(f"{csv}: covers {t.iloc[0].date()}..{t.iloc[-1].date()}, need {SPIN_FROM}..{END}")
    if (t.diff().dropna() != pd.Timedelta(days=1)).any():
        raise SystemExit(f"{csv}: not one row per day")
    main = d[(t >= pd.Timestamp(SPIN_FROM)) & (t <= pd.Timestamp(END))].copy()
    tm = pd.to_datetime(main[time_col])
    spin = main[(tm <= pd.Timestamp(SPIN_TO)) & ~((tm.dt.month == 2) & (tm.dt.day == 29) & (tm.dt.year == 1920))].copy()
    spin[time_col] = (pd.to_datetime(spin[time_col]) - pd.DateOffset(years=SHIFT_Y)).dt.strftime("%Y-%m-%d %H:%M:%S")
    pad = main.iloc[[-1]].copy(); pad[time_col] = f"{PAD} 00:00:00"
    out = pd.concat([spin, main, pad], ignore_index=True)
    tt = pd.to_datetime(out[time_col])
    assert (tt.diff().dropna() == pd.Timedelta(days=1)).all(), f"{csv}: gap after extension"
    assert tt.iloc[0] == pd.Timestamp(START) and tt.iloc[-1] == pd.Timestamp(PAD)
    out.to_csv(csv, index=False)
    return len(spin)


def main(lake, run):
    src = ROOT / "calibration" / "runs" / lake / "GLM"
    best = pd.read_csv(ROOT / "calibration" / "results" / lake / "best.csv").iloc[0]
    run = Path(run)
    shutil.rmtree(run, ignore_errors=True)
    shutil.copytree(src, run, ignore=shutil.ignore_patterns("output"))
    (run / "output").mkdir()
    n = extend(run / "meteo_file.csv", "Date")
    extend(run / "inflow_evapcomp.csv", "Time")
    t = (run / "glm3.nml").read_text()
    for k in PARS:
        t = set_nml(t, k, float(best[k]))
    t = set_nml(t, "start", f"'{START} 00:00:00'")
    t = set_nml(t, "stop", f"'{PAD} 00:00:00'")
    t = set_nml(t, "nsave", 24, "{:d}")
    (run / "glm3.nml").write_text(t)
    print(f"{lake}: {START}..{PAD} ({n} spin-up days), " + ", ".join(f"{k}={best[k]:.4g}" for k in PARS))


if __name__ == "__main__":
    main(*sys.argv[1:3])
