# GLM4 hydrodynamic calibration — ISIMIP local lakes, ISIMIP4a forcing

Goal: calibrate GLM 4.0.0 (built from the pinned glm-aed source, `../bin/glm-built`) on the
local lakes with ISIMIP4a forcing, then use the **median of the calibrated parameters** for
the global GLM4-AED3 runs (and the per-lake values for ISIMIP4a/4b local-lake runs).
Physics only here; water-quality (AED3) parameters are a separate, later step.

## Steps (in order)

| # | Script | Where | What |
|---|---|---|---|
| 1 | `build_lake_table.py` | local | lake list + harmonised temperature obs + hypsographs → `lakes.csv`, `obs/`, `hypso/` |
| 2 | `../forcing/extract_points_levante.py` | Levante login node | ISIMIP4a GSWP3-ERA5 at the 77 lake points (see `../forcing/README.md`) |
| 3 | `setup_lake_glm.R all` | local | per-lake LakeEnsemblR → GLM folder in `runs/<Lake>/` (window, spin-up, init profile, Kw base) |
| 4 | `patch_glm_closed.py all` | local | closed-lake fixes from ISIMIP4_LocalLakes_LERWQ (`lerwq/patches.py`): drop template outflows, evaporation-compensation inflow (constant level) |
| 5 | `sweep_defaults.R` | local | one run per lake with GLM defaults: checks all run, baseline RMSE, seconds per run |
| 6 | `../make_glm_launcher.sh` | each machine | `../bin/glm-run`: GLM with its build libraries pinned |
| 7 | `hydra_install_r_deps.sh` | Hydra login | R packages (data.table, DEoptim, ncdf4, glmtools) |
| 8 | `hydra_calibrate.slurm` | Hydra (sbatch array) | `calibrate_lake.R` per lake → `results/<Lake>/` |
| 9 | `collect_results.py` | local | `results_all.csv`, `median_parameters.csv` |

## Lakes (77)

The 73 ISIMIP3 local lakes plus the lakes new in ISIMIP4 (Edward, Grossglienickersee,
InreHarrsjon, MellerstaHarrsjon, Villasjon), minus "Zurich" — byte-identical duplicate of
LowerLakeZurich in the ISIMIP3 database. Carioca (ISIMIP4 water-quality tab only) has no data
files here. Lake list and metadata: ISIMIP4 local-lakes Google Sheet (saved in
`../reference/isimip4_locallakes_sheet/`).

Observations — `temp_files` / `hypso_file` in `lakes.csv` give the exact file per lake:
- `~/Documents/isi_cal/ISIMIP_Local_Lakes/LocalLakes/<Lake>/` (ISIMIP3 database, 65 lakes)
- `~/Documents/ISIMIP/WaterQuality/ISIMIP4_LocalLakes_WaterQuality/LocalLakes/<Lake>/`
  (ISIMIP4 delivery, 12 lakes; read with the LERWQ harmoniser and its `obs_catalog.yaml`)

## Setup choices (mirroring aemon-j/ISIMIP3_LER, commit 3ba8f23, in `../reference/`)

- **Window**: end = min(last obs, end of forcing); start = max(first obs, end − 10 y). First
  half calibrates (`obs_wtemp_cal.csv`), second half validates (`obs_wtemp_val.csv`).
  **Spin-up: 20 years** of the real ISIMIP4a forcing preceding the calibration start
  (ISIMIP3_LER used 1 recycled year). Initial profile: observed climatology ±1 month around
  the start day-of-year, ≤ 20 depths.
- **Depth**: from the observed hypsograph for every lake (ISIMIP3_LER took the larger of
  hypsograph and metadata; 28 lakes differed, by up to 22 m).
- **GLM namelist**: LakeEnsemblR `export_config(model = "GLM")` on ISIMIP3_LER's
  `LER_template/LakeEnsemblR.yaml` — all non-calibrated GLM settings as in ISIMIP3 GLM-LER.
- **Kw**: Håkanson (1995) 1.1925·zmax^−0.424 for every lake, times a calibrated scaling —
  the only rule that also exists for the 41,449 global lakes (no Secchi/Kw there).
- **Constant level**: closed lakes otherwise drain by evaporation (Allequash 8 → 6.4 m in
  10 y); see step 4.

## Calibration (`calibrate_lake.R`, `glm4_objective.R`)

| parameter | range | GLM default | as ISIMIP3_LER GLM? |
|---|---|---|---|
| `wind_factor` | 0.25–1.5 | 1 | yes (`wind_speed` scaling) |
| `sw_factor` | 0.7–1.3 | 1 | yes (`swr` scaling) |
| `kw_scaling` (× Håkanson Kw) | 0.5–3 | 1 | range widened (was ×0.7–1.3 of a per-lake Kw) |
| `coef_mix_hyp` | 0.1–2 | 0.5 | yes |
| `coef_mix_conv` | 0.1–0.3 | 0.2 | yes |
| `coef_mix_turb` | 0.35–0.65 | 0.51 | yes |

Objective: RMSE of water temperature at the observed depths (model interpolated to them),
calibration half. DEoptim, NP = 60 (10 × 6 parameters), strategy 2, ≤ 100 generations,
stops when the best RMSE improves < 0.1 % over 20 generations, seed 20261021. Failed runs
or runs matching < 50 % of observations score 999.

## Fixes found while setting up (all in the scripts, with comments)

| issue | effect | fix |
|---|---|---|
| ISIMIP3_LER `add_spin_up` parses dates in the local time zone, meteo is UTC | Pääjärvi (start 2012-03-26, after a DST switch): duplicated day → LER rejects meteo | run setup in UTC |
| some observations deeper than the hypsograph | LER rejects the initial profile (5 lakes) | initial profile ends 5 cm above the bottom |
| lakes observed at one depth (Taihu, Chao, Hulun) | **GLM4 hangs** in `initialise_lake` (no error) | isothermal 2-point initial profile — worth reporting upstream |
| Bosumtwi ISIMIP4 hypsometry: depth axis inverted | GLM: "H and A must be monotonically increasing" | detected and flipped automatically (ISIMIP3 file has the same numbers the right way up) |
| R module (2025a) and GLM build libraries (2024a) can't be loaded together on Hydra | `libnetcdf.so.19 not found` | `bin/glm-run` launcher |
| DEoptim 2.2-8 with a supplied cluster + `parallelType="parallel"` | closes the caller's cluster | `parallelType = "none"` + `cluster` |

## Open data questions

- **Taihu**: hypsograph reaches 13.31 m, volume 8.7 km³ vs metadata max 4.8 m, mean 2.1 m
  (~5.1 km³); observations at the surface only. **Decision: keep the observed hypsograph**,
  Taihu included like every lake.
- **Harp**: ISIMIP3 profiles (147 days, 2004–2014) **combined** with the ISIMIP4 daily logger
  file (1,599 days, 2011–2015) → 1,696 days (decision 2026-10-21).
- **Few calibration days** (Okauchee 12, Bosumtwi 13, Biel 20, Edward 24, …): **all lakes
  kept**, no minimum (decision 2026-10-21).

## Checks so far

- Uncalibrated (GLM defaults), 20-year spin-up: median RMSE 2.05 °C over 77 lakes; median
  25 s per run, Kivu 235 s (laptop, 14 parallel) (`sweep_defaults.csv`).
- Runtime tiers for the job limits (`lakes_tier{A,B,C}.txt`, worst case 101 × 60 runs / 16
  cores): A 55 lakes ≤ 6 h, B 21 lakes ≤ 18 h, C Kivu 48 h. ~4,400 core-h worst case.
- Memory (16 R workers each read the model output): tier A peaked 14–19 GB of 32 GB;
  Rappbode (tier B) 26.8 GB; **Kivu ran out of 32 GB** (job 13361138) -> resubmitted with 192 GB
  (job 13361966); pending tier-B tasks raised to 96 GB with `scontrol update`.
- Same seed, laptop vs Hydra: DEoptim generations 1–2 identical to all digits.
- Same Erken run on the laptop and on Hydra: RMSE 1.399034 °C on both.
