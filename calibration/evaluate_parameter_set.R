# Run every local lake with ONE shared parameter set (e.g. the cross-lake medians) and score it
# on the calibration and validation halves -- the transfer test for the global run, where each
# lake gets the shared set and no calibration.
#   [NML_SET="atm_stab=1"] Rscript calibration/evaluate_parameter_set.R <params.csv> <out.csv> [cores]
# params.csv: either one shared set (columns parameter,median: median_parameters.csv from
# collect_results.py) or one row per lake (column lake + the 6 PARS names, e.g. a rule
# predicted from lake attributes); only the 6 calibrated parameters (glm4_objective.R PARS) are read.
suppressMessages({ library(data.table); library(parallel) })
Sys.setenv(TZ = "UTC")
args <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
a <- commandArgs(TRUE)
if (length(a) < 2) stop("usage: evaluate_parameter_set.R <params.csv> <out.csv> [cores]")
cores <- if (length(a) >= 3) as.integer(a[3]) else max(1L, detectCores() - 1L)
source(file.path(here, "glm4_objective.R"))
pset <- fread(a[1])
# two formats: one shared set (columns parameter, median) or one row per lake (lake + PARS$name)
per_lake <- "lake" %in% names(pset)
if (!per_lake) {
  p_shared <- setNames(pset$median, pset$parameter)[PARS$name]
  if (anyNA(p_shared)) stop("parameter file lacks: ", paste(PARS$name[is.na(p_shared)], collapse = ", "))
}
par_for <- function(lake) {
  if (!per_lake) return(p_shared)
  lk <- lake                                 # inside pset[...] the name `lake` is the column
  r <- pset[pset$lake == lk]
  if (nrow(r) != 1) stop("no parameter row for ", lake)
  unlist(r[, PARS$name, with = FALSE])
}
lakes <- sort(basename(dirname(Sys.glob(file.path(here, "runs", "*", "window.csv")))))
if (per_lake) lakes <- intersect(lakes, pset$lake)
glm_bin <- normalizePath(file.path(dirname(here), "bin", "glm-run"))

one <- function(lake) {
  source(file.path(here, "glm4_objective.R"))
  d <- file.path(here, "runs", lake); w <- fread(file.path(d, "window.csv"))
  rd <- tempfile(paste0("eval_", lake, "_")); dir.create(rd)
  file.copy(list.files(file.path(d, "GLM"), full.names = TRUE), rd, recursive = TRUE)
  # optional extra namelist settings, e.g. NML_SET="atm_stab=1" (comma-separated key=value)
  extra <- Sys.getenv("NML_SET", "")
  if (nzchar(extra)) {
    kv <- strsplit(strsplit(extra, ",")[[1]], "=")
    set_nml_values(file.path(rd, "glm3.nml"), setNames(lapply(kv, function(x) as.numeric(x[2])), sapply(kv, `[`, 1)))
  }
  rdo <- function(f) fread(file.path(d, f))[, datetime := as.POSIXct(datetime, tz = "UTC")]
  oc <- rdo("obs_wtemp_cal.csv"); ov <- rdo("obs_wtemp_val.csv")
  rc <- glm4_rmse(unname(par_for(lake)), rd, glm_bin, w$kw_base, oc)
  nc <- file.path(rd, "output", "output.nc")
  mv <- if (file.exists(nc)) tryCatch(match_obs(nc, ov)[is.finite(sim)], error = function(e) NULL) else NULL
  rv <- if (!is.null(mv) && nrow(mv)) sqrt(mean((mv$sim - mv$Water_Temperature_celsius)^2)) else NA
  unlink(rd, recursive = TRUE)
  data.table(lake = lake, rmse_cal_shared = rc, rmse_val_shared = rv)
}
res <- rbindlist(mclapply(lakes, one, mc.cores = cores, mc.preschedule = FALSE))
res[, parameter_file := normalizePath(a[1])]
fwrite(res, a[2])
cat(sprintf("%d lakes; shared set: median RMSE cal %.2f, val %.2f C\n", nrow(res),
            median(res$rmse_cal_shared, na.rm = TRUE), median(res$rmse_val_shared, na.rm = TRUE)))
