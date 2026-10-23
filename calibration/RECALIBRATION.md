# Recalibration with the patched GLM (2026-10-22)

The first calibration (results_v1_unpatched/, median_parameters_v1_unpatched.csv; 75 of 77
lakes, Kivu and Tahoe cancelled) used GLM 4.0.0 without our source fixes. Two bugs in GLM's
surface routine changed the physics of every lake (patches/, build_glm4.sh):
- GLM_daily_shortwave.patch: daily shortwave had the wrong-sign solar declination (winter day
  lengths in northern summers), no sun in polar day/night, and the half-sine centred on midnight;
- GLM_ice_turbulent_wind.patch: air-ice turbulent fluxes used the ~0 under-ice wind.
With the old medians on the patched model: validation RMSE (median, 77 lakes) 2.089 -> 2.051 C,
lakes > 66N 4.36 -> 2.37 C, but 25 lakes worse -- the medians had compensated for the bugs.
Hence a full recalibration (same setups, ranges, DEoptim settings and seed; only the binary
differs). Kivu and Tahoe get 48 h (Tahoe needed ~20 min per generation).

## Changes to the recalibration inputs (2026-10-22, after it was queued)
- **Bosumtwi observations**: 10 values with the depth entered as temperature removed
  (build_lake_table.py OBS_QC); its setup rebuilt (916 calibration / 906 validation values).
- **Atmospheric stability (user decisions)**: `atm_stab = 1` for |lat| <= 23.44 OR elevation > 1500 m
  (Argyle, Bosumtwi, Edward, Kivu, Laramie, Ngoring, Tahoe), `atm_stab = 0` for the other 70 lakes;
  set by patch_glm_closed.py (set_stability), the same rule as the global setup
  (analysis/stability/README.md). Evidence, GLM defaults, validation RMSE: tropical lakes
  2.57 -> 1.67 C (Edward 4.81 -> 1.16), but all lakes with the correction 2.06 -> 2.18 C
  (eval_default_stab0/1.csv).
- If the Hydra recalibration of these seven lakes (and Bosumtwi's new observations) started before these inputs were synced (Hydra
  login was refused from 2026-10-22 ~10:00), rerun those seven lakes.

## Provisional medians, 75 of 77 lakes (2026-10-23 morning; Tahoe and Kivu still running)
`collect_results.py` on calibration/results (fixed GLM, stability rule, Bosumtwi QC; Argyle and
Ngoring re-run with the new inputs; Kilpisjarvi re-run after a silent cluster-start failure):
wind_factor 1.057, sw_factor 0.970, kw_scaling 1.282, coef_mix_hyp 1.295, coef_mix_conv 0.170,
coef_mix_turb 0.501 (v1 unpatched: 1.059, 0.976, 1.266, 1.165, 0.203, 0.498). Median RMSE of each
lake's own calibration: 1.05 (cal) / 1.37 °C (val); GLM defaults 2.05 / 1.99 °C.
median_parameters.csv is overwritten when Tahoe and Kivu finish.
Transfer test with the provisional medians (evaluate_parameter_set.R, all 77 lakes, eval_median_recal75.csv;
comparison in transfer_test_recal75.csv): validation RMSE median 1.98 °C, GLM defaults with the stability
rule 1.99 °C (eval_default_stab0/1 combined by the rule), own calibration 1.37 °C; medians better in 26
lakes, worse in 43, same in 8.
