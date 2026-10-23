# Build one lake's GLM folder for the GLM4 hydrodynamic calibration.
#
#   Rscript calibration/setup_lake_glm.R <Lake> [<Lake> ...]      # or: all
#
# Mirrors aemon-j/ISIMIP3_LER (reference/ISIMIP3_LER, commit 3ba8f23), which built the
# ISIMIP3 GLM-LER setups, reusing its helper functions unchanged:
#   convert_ISIMIP_to_LER (meteo), create_init_profile2
#   (climatological profile around the start day-of-year, max 20 depths),
#   get_output_resolution, and LakeEnsemblR::export_config(model = "GLM") on its
#   LER_template/LakeEnsemblR.yaml.
# Choices, as in ISIMIP3_LER's calib_type = "cal_project":
#   * window: end = min(last obs, end of forcing); start = max(first obs, end - 10 y);
#     first half -> calibration (obs_wtemp_cal.csv), second half -> validation
#     (obs_wtemp_val.csv).
#   * spin-up: 20 years (user requirement, 2026-10-21; ISIMIP3_LER used 1 recycled year).
#     The real forcing of the 20 years before the calibration start is used -- ISIMIP4a
#     starts in 1901 and every window starts after 1990, so nothing needs recycling
#     (ISIMIP3_LER's add_spin_up recycles the years AFTER the start, which a 10-year
#     window cannot supply 20 of).
#   * Kw base = Hakanson (1995) 1.1925 * zmax^-0.424 for every lake -- the only rule
#     that also exists for the global lakes (no measured Kw/Secchi there), so the
#     calibrated Kw scaling transfers.
# Differences from ISIMIP3 GLM-LER: ISIMIP4a GSWP3-ERA5 forcing (forcing/isimip4a),
# observations/hypsographs from calibration/obs + calibration/hypso
# (build_lake_table.py), and the model is GLM4 (bin/glm-built), not GLM3r.
#
# Writes calibration/runs/<Lake>/ with LakeEnsemblR.yaml, meteo.csv, hypsograph.csv,
# obs_wtemp*.csv, init_temp_prof.csv, window.csv and GLM/glm3.nml (+ GLM/meteo_file.csv).

# UTC throughout (kept although add_spin_up is no longer used): it parsed the start date in the session time zone
# while meteo times are UTC; across a DST switch that duplicated a day (Paajarvi, window
# starting 2012-03-26) and LakeEnsemblR then rejects the meteo as irregular.
Sys.setenv(TZ = "UTC")
suppressMessages({
  library(data.table); library(lubridate); library(LakeEnsemblR); library(configr)
})

args <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))   # calibration/
root <- dirname(here)                                                                 # GLM4-AED3/
ler3 <- file.path(root, "reference", "ISIMIP3_LER")
for (f in c("convert_ISIMIP_to_LER.R", "create_init_profile2.R",
            "get_output_resolution.R")) source(file.path(ler3, "Scripts", f))

MAX_DURATION_RUN <- 10      # years, as ISIMIP3_LER Settings.R
SPIN_UP_YEARS    <- 20L     # user requirement (ISIMIP3_LER: 1)
MAX_DEPTHS_INIT  <- 20L     # as ISIMIP3_LER Settings.R

lakes_tab <- fread(file.path(here, "lakes.csv"))
sel <- commandArgs(TRUE)
if (length(sel) == 1 && sel == "all") sel <- lakes_tab$lake

setup_lake <- function(lk) {
  lake <- lk
  m <- lakes_tab[lakes_tab[["lake"]] %in% lk]
  if (nrow(m) != 1) stop("lake not in lakes.csv: ", lake)
  folder <- file.path(here, "runs", lake)
  unlink(folder, recursive = TRUE); dir.create(folder, recursive = TRUE)

  # --- inputs
  met_src <- file.path(root, "forcing", "isimip4a", "gswp3-era5",
                       sprintf("gswp3-era5_obsclim_%s_daily_1901_2025.csv", tolower(lake)))
  file.copy(met_src, folder)                              # convert_ISIMIP_to_LER finds *dddd.csv
  convert_ISIMIP_to_LER(folder, rounding_dec = 4L)
  file.remove(file.path(folder, basename(met_src)))
  hyps <- fread(file.path(here, "hypso", paste0(lake, ".csv")))
  setnames(hyps, c("Depth_meter", "Area_meterSquared"))
  fwrite(hyps, file.path(folder, "hypsograph.csv"))
  obs <- fread(file.path(here, "obs", paste0(lake, "_temp.csv")))
  obs <- obs[, .(datetime = as.POSIXct(paste(datetime, "00:00:00"), tz = "UTC"),
                 Depth_meter = depth, Water_Temperature_celsius = value)]
  setorder(obs, datetime, Depth_meter)

  # --- depth, output resolution, Kw
  # Depth from the observed hypsograph, not the metadata (user decision 2026-10-21; ISIMIP3_LER
  # used the larger of the two). GLM's geometry only spans the hypsograph anyway: a deeper
  # metadata value only spilled over the crest at the start and skewed Hakanson Kw
  # (28 lakes, up to 22 m, e.g. Windermere 64 vs 42 m).
  max_depth <- hyps[, max(Depth_meter)]
  output_res <- get_output_resolution(max_depth)
  kw_base <- round(1.1925 * max_depth^(-0.424), 3L)       # Hakanson (1995)

  # --- window (ISIMIP3_LER cal_project rule)
  met <- fread(file.path(folder, "meteo.csv"))
  end_val <- ceiling_date(min(max(met$datetime), max(obs$datetime)), "days")
  start_obs <- min(obs$datetime)
  start_cal <- if (as.numeric(difftime(end_val, start_obs, units = "days")) / 365 < MAX_DURATION_RUN)
    start_obs else end_val - years(MAX_DURATION_RUN)
  start_cal <- floor_date(start_cal, "days")
  start_val <- ceiling_date(start_cal + days(ceiling(as.numeric(difftime(end_val, start_cal, units = "days")) / 2)), "days")
  o_cal <- obs[datetime >= start_cal & datetime < start_val]
  o_val <- obs[datetime >= start_val & datetime <= end_val]
  fmt <- function(d) { d <- copy(d); d[, datetime := format(datetime, "%Y-%m-%d %H:%M:%S")]; d }
  fwrite(fmt(obs[datetime >= start_cal & datetime <= end_val]), file.path(folder, "obs_wtemp.csv"))
  fwrite(fmt(o_cal), file.path(folder, "obs_wtemp_cal.csv"))
  fwrite(fmt(o_val), file.path(folder, "obs_wtemp_val.csv"))

  # --- initial profile (climatology around the start day-of-year)
  init <- create_init_profile2(copy(obs), start_date = start_cal, max_depths = MAX_DEPTHS_INIT)
  # Some lakes have a few observations below the hypsograph's deepest point (e.g. 17 m in a
  # 16.46 m lake); LakeEnsemblR rejects an initial profile reaching the lake depth, so keep
  # it >= 5 cm above the bottom (same fix as ISIMIP4_LocalLakes_LERWQ MODEL_ISSUES #18).
  init <- init[Depth_meter <= max_depth - 0.05]
  # Lakes observed at a single depth (Taihu, Chao, Hulun) give a 1-point profile, and GLM4
  # then hangs in initialise_lake (no error, no output). Use an isothermal 2-point profile.
  if (nrow(init) < 2) {
    t0 <- mean(init$Water_Temperature_celsius)
    init <- data.table(Depth_meter = c(0, round(max_depth - 0.05, 2)), Water_Temperature_celsius = round(t0, 2))
  }
  fwrite(init, file.path(folder, "init_temp_prof.csv"))

  # --- LakeEnsemblR config (ISIMIP3_LER template; only per-lake values changed)
  cfg <- file.path(folder, "LakeEnsemblR.yaml")
  file.copy(file.path(ler3, "LER_template", "LakeEnsemblR.yaml"), cfg)
  set <- function(..., value) input_yaml_multiple(file = cfg, ..., value = value, verbose = FALSE)
  set(key1 = "location", key2 = "name", value = lake)
  set(key1 = "location", key2 = "latitude", value = m$lat)
  set(key1 = "location", key2 = "longitude", value = m$lon)
  set(key1 = "location", key2 = "elevation", value = m$elev_m)
  set(key1 = "location", key2 = "depth", value = max_depth)
  set(key1 = "location", key2 = "init_depth", value = max_depth)
  set(key1 = "time", key2 = "start", value = format(start_cal, "%Y-%m-%d %H:%M:%S"))
  set(key1 = "time", key2 = "stop", value = format(end_val, "%Y-%m-%d %H:%M:%S"))
  set(key1 = "output", key2 = "depths", value = output_res)
  set(key1 = "input", key2 = "light", key3 = "Kw", key4 = "GLM", value = kw_base)
  spinup_start <- start_cal - years(SPIN_UP_YEARS)
  if (spinup_start < min(met$datetime)) stop("forcing starts after the spin-up start")
  set(key1 = "time", key2 = "start", value = format(spinup_start, "%Y-%m-%d %H:%M:%S"))

  # --- GLM folder
  export_config(config_file = "LakeEnsemblR.yaml", model = "GLM", folder = folder)

  w <- data.table(lake = lake, max_depth = max_depth, output_res = output_res, kw_base = kw_base,
                  spinup_start = format(start_cal - years(SPIN_UP_YEARS), "%Y-%m-%d"),
                  cal_start = format(start_cal, "%Y-%m-%d"), val_start = format(start_val, "%Y-%m-%d"),
                  end = format(end_val, "%Y-%m-%d"), n_obs_cal = nrow(o_cal), n_obs_val = nrow(o_val),
                  days_obs_cal = uniqueN(o_cal$datetime), days_obs_val = uniqueN(o_val$datetime))
  fwrite(w, file.path(folder, "window.csv"))
  w
}

res <- rbindlist(lapply(sel, function(l) {
  r <- tryCatch(setup_lake(l), error = function(e) { message(l, ": FAILED: ", conditionMessage(e)); NULL })
  if (!is.null(r)) message(sprintf("%-20s cal %s..%s (%d days) val ..%s (%d days)  zmax=%.1f Kw0=%.3f",
                                   l, r$cal_start, r$val_start, r$days_obs_cal, r$end, r$days_obs_val,
                                   r$max_depth, r$kw_base))
  r
}))
if (length(sel) > 1) fwrite(res, file.path(here, "windows.csv"))
