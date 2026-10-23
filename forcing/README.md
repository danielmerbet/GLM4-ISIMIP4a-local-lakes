# Local-lake forcing

## `isimip4a/gswp3-era5/` — ISIMIP4a obsclim GSWP3-ERA5, 1901–2025, 77 local lakes

- Source: DKRZ Levante
  `/work/bb0820/ISIMIP/ISIMIP4a/InputData/climate/atmosphere/obsclim/global/daily/historical/GSWP3-ERA5/`
  (143 files, 11 variables × 13 decade files, ~301 GB; never transferred).
- Extracted **on Levante** (login node, `nice -n 19`, 6 processes) on 2026-10-21 by
  `extract_points_levante.py`, started with `isimip4a/run.sh`; full log `isimip4a/extract.log`.
  Lake coordinates: `lakes_points.csv` (from `../calibration/lakes.csv`).
- Cell choice: nearest 0.5° cell centre; nearest valid cell within 2 cells if masked.
  `points.csv` records the cell per lake — **no lake needed a fallback cell**.
- One CSV per lake, `gswp3-era5_obsclim_<lake>_daily_1901_2025.csv`, 45,656 days each,
  columns `time,tas,tasmin,tasmax,hurs,huss,pr,prsn,rsds,rlds,ps,sfcwind` — same layout and
  units as the ISIMIP3 local-lake files (K, %, kg kg-1, kg m-2 s-1, W m-2, Pa, m s-1).
  No missing values.

### Check against ISIMIP3a GSWP3-W5E5 (1979–2019, 72 lakes present in both)

Median correlation of daily `tas` 0.999, median mean difference +0.02 K (IQR −0.30…+0.59).
Outliers, all explained:
- Bosumtwi (r 0.86, −0.9 K): ISIMIP3a forcing used protocol coordinates ~300 km off;
  ISIMIP4 extraction uses the corrected ones (ISIMIP4_LocalLakes_LERWQ config/lakes.csv).
- Kivu (r 0.78, −0.1 K): tropical, tiny daily variance; mean agrees.
- Kilpisjärvi +1.8 K, Ngoring −1.8 K, Thun −1.7 K: complex-terrain/arctic cells where the
  W5E5 and ERA5 references differ; r ≥ 0.99, so not a coordinate problem.

Not covered by ISIMIP3a (new lakes): Edward, Grossglienickersee, InreHarrsjon,
MellerstaHarrsjon, Villasjon.

## ISIMIP4a GSWP3-ERA5 global files on Hydra (for the global runs)

`/data/brussel/vo/000/bvo00012/data/dataset/ISIMIP/ISIMIP4a/InputData/climate/atmosphere/obsclim/global/daily/historical/GSWP3-ERA5/`
(group bvo00012 shared dataset folder, same layout as DKRZ). Mirrored 2026-10-21 directly
DKRZ -> Hydra with `hydra_mirror_isimip4a.sh` (rsync on the Hydra login1 node in tmux,
~66 min, run at nice 19 / lowest I/O priority). `verify_isimip4a_mirror.sh`: **143/143 files,
every byte size identical** to the DKRZ listing in `provenance/20261021T133707Z/`; 288 GB.
