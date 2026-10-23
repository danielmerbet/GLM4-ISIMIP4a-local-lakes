# Objective function for the GLM4 hydrodynamic calibration: set parameters in a copy of
# a lake's GLM folder, run GLM4, return the temperature RMSE against observations.
# Sourced by calibrate_lake.R; can be used on its own for checks:
#   source("calibration/glm4_objective.R"); glm4_rmse(...)

suppressMessages({ library(data.table); library(glmtools) })

# Calibrated parameters and ranges. Same six as aemon-j/ISIMIP3_LER's GLM calibration
# (LER_template/LakeEnsemblR.yaml) with the same ranges, except Kw: there a x0.7-1.3
# scaling of a per-lake base; here a scaling of the Hakanson depth rule for every lake
# (see setup_lake_glm.R), widened to x0.5-3 because Hakanson alone can be far from a
# humic or turbid lake's real extinction.
PARS <- data.table(
  name  = c("wind_factor", "sw_factor", "kw_scaling", "coef_mix_hyp", "coef_mix_conv", "coef_mix_turb"),
  lower = c(0.25,          0.7,         0.5,          0.1,            0.1,             0.35),
  upper = c(1.5,           1.3,         3.0,          2.0,            0.3,             0.65),
  glm_default = c(1,       1,           1,            0.5,            0.2,             0.51))

# Replace "name = value" in a namelist without touching anything else.
set_nml_values <- function(nml_file, values) {
  x <- readLines(nml_file)
  for (k in names(values)) {
    pat <- sprintf("^(\\s*%s\\s*=\\s*)[^!,/]*", k)
    hit <- grep(pat, x)
    if (length(hit) != 1) stop(nml_file, ": expected exactly one '", k, "' entry, found ", length(hit))
    x[hit] <- sub(pat, paste0("\\1", format(values[[k]], digits = 10)), x[hit])
  }
  writeLines(x, nml_file)
}

# Model-vs-obs pairs: obs at their own depths (below the surface), model interpolated to them.
match_obs <- function(nc_file, obs) {
  depths <- sort(unique(obs$Depth_meter))
  sim <- get_temp(nc_file, reference = "surface", z_out = depths)
  sim <- melt(as.data.table(sim), id.vars = "DateTime", variable.name = "var", value.name = "sim")
  sim[, Depth_meter := as.numeric(sub("temp_", "", var))][, var := NULL]
  sim[, date := as.Date(DateTime)][, DateTime := NULL]
  o <- copy(obs)[, date := as.Date(datetime)]
  merge(o, sim, by = c("date", "Depth_meter"))
}

# Run GLM4 in run_dir with parameter vector p (named as PARS$name); returns RMSE on obs.
# Failed runs return a large penalty so the optimiser moves away from them.
glm4_rmse <- function(p, run_dir, glm_bin, kw_base, obs, penalty = 999) {
  # resolve paths now: arguments are lazy, and we setwd() below
  run_dir <- normalizePath(run_dir, mustWork = TRUE)
  glm_bin <- normalizePath(glm_bin, mustWork = TRUE)
  names(p) <- PARS$name
  vals <- list(wind_factor = p[["wind_factor"]], sw_factor = p[["sw_factor"]],
               Kw = kw_base * p[["kw_scaling"]], coef_mix_hyp = p[["coef_mix_hyp"]],
               coef_mix_conv = p[["coef_mix_conv"]], coef_mix_turb = p[["coef_mix_turb"]])
  set_nml_values(file.path(run_dir, "glm3.nml"), vals)
  nc_file <- file.path(run_dir, "output", "output.nc")
  unlink(nc_file)
  old <- setwd(run_dir); on.exit(setwd(old))
  rc <- system2(glm_bin, c("--nml", "glm3.nml"), stdout = FALSE, stderr = FALSE)
  setwd(old)
  if (rc != 0 || !file.exists(nc_file)) return(penalty)
  pr <- tryCatch(match_obs(nc_file, obs), error = function(e) NULL)
  if (is.null(pr)) return(penalty)
  pr <- pr[is.finite(sim)]
  if (nrow(pr) < 0.5 * nrow(obs)) return(penalty)   # model failed or lost most obs depths
  r <- sqrt(mean((pr$sim - pr$Water_Temperature_celsius)^2))
  if (!is.finite(r)) penalty else r
}
