# ISIMIP local-lake output files — what the ISIMIP3 LER runs delivered

**Authoritative reference: the published files on DKRZ**,
`levante:/work/bb0820/ISIMIP/ISIMIP3b/OutputData/lakes_local/{GLM,GOTM,Simstrat,FLake}-LER/<gcm>/{historical,future,pre-industrial}/`
(GLM-LER: 19,710 files, 60 GB = 73 lakes x 9 variables x 6 scenario periods x 5 GCMs, all
`2015soc_default`). ISIMIP3a has no LER submissions; ISIMIP4a `OutputData` has no lake files yet
(checked 2026-10-21). Where the local copies below differ, **the DKRZ version wins** — see
"Corrections from the DKRZ originals" at the end.

Reference for building the GLM4-AED3 local-lake NetCDFs. Sources checked 2026-10-21:
- **Published files**: 965 ISIMIP3 LER NetCDFs on disk, e.g.
  `~/Documents/intoDBP/CC_analysis/sau/ISIMIP3b_lakedata/` (QC'd by ISIMIP: CDI/`isimip_*`
  attributes added) and `~/Documents/Colaboraciones/HackingLimnology2024-ISIMIP/.../data/`
  (as written by the LER scripts, before ISIMIP's QC).
- **Code that wrote them**: aemon-j/ISIMIP3_LER `Scripts/5_Create_ISIMIP_netcdfs.R` and
  `Scripts/write_isimip_netcdf.R` (`reference/ISIMIP3_LER`, commit 3ba8f23).
- **Protocol variable table**: `~/Documents/ISIMIP/ISIMIP4T_simulations/config/variables.yaml`
  (ISIMIP3 lakes table carried forward — ISIMIP4's lake page had no sector table yet, 2026-08-24).

## File naming

```
<model>_<climate-forcing>_<bias-adjustment>_<climate-scenario>_<soc>_<sens>_<variable>_<lake>_daily_<y0>_<y1>.nc
glm-ler_mpi-esm1-2-hr_w5e5_ssp126_2015soc_default_surftemp_sau_daily_2015_2100.nc        (ISIMIP3b)
```
- `<model>` lower case, `-ler` suffix in ISIMIP3 (`glm-ler`, `gotm-ler`, `simstrat-ler`, `flake-ler`).
- `<lake>` = the "Lake Name in file name (reporting)" column of the local-lakes sheet
  (`lowerlakezurich`, `bigmuskellunge`, …), not the short name.
- One variable per file, daily. picontrol split into 1850–2014 and 2015–2100 files.
- ISIMIP4a: forcing `gswp3-era5`, scenario `obsclim`, soc `histsoc` (and `2021soc`), 1901–2025.
  ISIMIP4 has **not** published a lake-sector file pattern yet; ISIMIP4T decided to use the
  ISIMIP3a pattern. Whether ISIMIP4a keeps a bias-adjustment token (ISIMIP3a files have none
  for obsclim) must be checked against the protocol before submission.

## Variables the ISIMIP3 LER runs delivered (all four models)

| variable | units | shape | how LER derived it |
|---|---|---|---|
| `watertemp` | K | (time, levlak, lat, lon) | model profile at fixed output depths, °C + 273.15 |
| `surftemp` | K | (time, lat, lon) | first (shallowest) output depth |
| `bottemp` | K | (time, lat, lon) | deepest output depth with a value that day |
| `strat` | 1 | (time, lat, lon) | 1 if ρ(bottom) − ρ(surface) ≥ 0.1 kg m⁻³ (density from temperature only, `temp_to_dens.R`) |
| `thermodepth` | m | (time, lat, lon) | `rLakeAnalyzer::ts.thermo.depth` defaults; NA → 0 |
| `icethick` | m | (time, lat, lon) | LER `ice_height` |
| `ice` | 1 | (time, lat, lon) | 1 if icethick > 0 |
| `sensheatf-total` | W m⁻² | (time, lat, lon) | −1 × LER `q_sens` (positive upwards) |
| `latentheatf` | W m⁻² | (time, lat, lon) | −1 × LER `q_lat` |

Protocol variables **not** delivered by ISIMIP3 LER, but in the table: `icetemp` (K),
`lakealbedo` (1), `turbdiffheat` (m² s⁻¹, profile), `extcoeff` (m⁻¹), and **`do`
(dissolved oxygen, mole m⁻³, profile — "full profile, or mean epi and mean hypo")**, the one
water-quality variable in the table so far. Spin-up years are cut before writing.

## NetCDF structure

- NetCDF-4 (published files: `NETCDF4_CLASSIC` after ISIMIP's QC), zlib level 5, no shuffle,
  float32 data, `_FillValue` = `missing_value` = 1e20.
- Dimensions `time` (unlimited), `lat` (1), `lon` (1); profiles add `levlak` and `bnds` (2).
- `time`: float64, `days since 1601-01-01 00:00:00` (ISIMIP3b; **ISIMIP4a reference year is
  1901** per ISIMIP4T experiments.yaml), `calendar = proleptic_gregorian`, `axis = T`.
- `lon` in degrees_east **0–360** (negative longitudes + 360), `lat` degrees_north; both with
  `standard_name`, `axis`.
- Profiles: `levlak` = layer index 1…n (`units "-"`, `positive = down`, `axis = Z`,
  `standard_name = water_layer`); `depth(levlak, lat, lon)` = layer-centre depth (m, positive
  down, `bounds = depth_bnds`); `depth_bnds(bnds, levlak, lat, lon)` = layer top/bottom.
  Output depths: LER's fixed grid, resolution by lake depth (`get_output_resolution.R`:
  0.25 m ≤ 2 m, 0.5 m ≤ 20 m, 1 m ≤ 50 m, 2 m ≤ 100 m, 5 m ≤ 300 m, 10 m beyond).
- Data variable attributes: `standard_name` = variable name, `long_name`, `units`.
- Global attributes: `contact`, `institution`, `comment` ("Data prepared for ISIMIP3b");
  ISIMIP adds `isimip_id`, `isimip_protocol_version`, `isimip_qc_*` at publication.

## Corrections from the DKRZ originals (checked 2026-10-21, GLM-LER gfdl-esm4 Mendota)

- **Lake IDs in file names are hyphenated, umlauts transliterated**: `big-muskellunge`,
  `black-oak`, `burley-griffin`, `crystal-bog`, `crystal-lake`, `esthwaite-water`,
  `falling-creek`, `great`, `lower-zurich`, `mt-bold`, `mueggelsee`, `nohipalo-mustjaerv`,
  `nohipalo-valgejaerv`, `paaijarvi`, `trout-bog`, `trout-lake`, `two-sisters`, `vortsjaerv`, …
  The local copies in `intoDBP/` had been renamed (`fix_files_isimip.sh` strips the hyphens).
  The Google Sheet column "Lake Name in file name (reporting)" equals the DKRZ IDs for all
  ISIMIP3 lakes. DKRZ also has **both `zurich` and `lower-zurich`** (the duplicate lake was
  delivered twice in ISIMIP3).
- **New ISIMIP4 lakes have non-conforming IDs in the sheet**: `Grossglienickersee` (capital),
  `inre-harrsjön`, `mellersta-harrsjön`, `villasjön` (non-ASCII; ISIMIP3 practice would give
  `-harrsjoen`, `villasjoen`), `edward`. To be confirmed with the ISIMIP coordinators.
- `depth` attribute `standard_name = "depth_below_surface"` (not `"depth"`).
- `depth_bnds` dimensions are **`(levlak, bnds)`** (the LER script wrote `(bnds, levlak, lat, lon)`).
- Top and bottom layers are half-thickness: Mendota (25 m, 1 m resolution) centres
  0.25, 1, 2, …, 24, 24.75 m; bounds 0–0.5, 0.5–1.5, …, 24.5–25.
- All variables `NETCDF4_CLASSIC`, zlib 5, `time` chunk 512, data chunked one time step.
- `watertemp` files were reprocessed with NCO (QC 2024-07-09, protocol version 8587111c…);
  the scalar variables carry CDI attributes (QC 2023-09-27, protocol 95fa1f08…).
- Longitude 0–360 confirmed (Mendota 270.59, Alqueva 352.51).
