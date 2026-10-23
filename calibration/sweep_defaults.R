# One GLM4 run per lake with GLM's default parameters (before calibration): checks every
# lake runs, gives the baseline RMSE, and measures seconds per run (for job time limits).
#   Rscript calibration/sweep_defaults.R [cores]   -> calibration/sweep_defaults.csv
suppressMessages({ library(data.table); library(parallel) })
Sys.setenv(TZ = "UTC")
args <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
cores <- if (length(commandArgs(TRUE))) as.integer(commandArgs(TRUE)[1]) else max(1L, detectCores() - 1L)
lakes <- sort(basename(dirname(Sys.glob(file.path(here, "runs", "*", "window.csv")))))
glm_bin <- normalizePath(file.path(dirname(here), "bin", "glm-run"))   # launcher from make_glm_launcher.sh

one <- function(lake) {
  source(file.path(here, "glm4_objective.R"))
  d <- file.path(here, "runs", lake); w <- fread(file.path(d, "window.csv"))
  rd <- tempfile(paste0("sweep_", lake, "_")); dir.create(rd)
  file.copy(list.files(file.path(d, "GLM"), full.names = TRUE), rd, recursive = TRUE)
  rdo <- function(f) fread(file.path(d, f))[, datetime := as.POSIXct(datetime, tz = "UTC")]
  oc <- rdo("obs_wtemp_cal.csv"); ov <- rdo("obs_wtemp_val.csv")
  t0 <- Sys.time()
  rc <- glm4_rmse(PARS$glm_default, rd, glm_bin, w$kw_base, oc)
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  nc <- file.path(rd, "output", "output.nc")
  mv <- if (file.exists(nc)) tryCatch(match_obs(nc, ov)[is.finite(sim)], error = function(e) NULL) else NULL
  rv <- if (!is.null(mv) && nrow(mv)) sqrt(mean((mv$sim - mv$Water_Temperature_celsius)^2)) else NA
  unlink(rd, recursive = TRUE)
  data.table(lake = lake, rmse_cal_default = rc, rmse_val_default = rv,
             val_matched = if (is.null(mv)) 0L else nrow(mv), n_obs_val = nrow(ov),
             seconds_per_run = round(secs, 1), sim_years = round(as.numeric(as.Date(w$end) - as.Date(w$spinup_start)) / 365.25, 1),
             max_depth = w$max_depth)
}
cl <- makeCluster(cores); clusterExport(cl, c("here", "glm_bin"))
invisible(clusterEvalQ(cl, { suppressMessages(library(data.table)); Sys.setenv(TZ = "UTC") }))
res <- rbindlist(parLapplyLB(cl, lakes, one)); stopCluster(cl)
setorder(res, -seconds_per_run)
fwrite(res, file.path(here, "sweep_defaults.csv"))
cat("lakes:", nrow(res), " failed (RMSE 999):", res[rmse_cal_default >= 999, paste(lake, collapse = ", ")], "\n")
print(res[, .(median_rmse_cal = median(rmse_cal_default[rmse_cal_default < 999]),
              median_s = median(seconds_per_run), max_s = max(seconds_per_run))])
print(head(res, 10))
