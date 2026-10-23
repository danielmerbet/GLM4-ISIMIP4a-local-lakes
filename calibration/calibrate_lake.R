# Calibrate GLM4 hydrodynamics for one lake with differential evolution (DEoptim).
#
#   Rscript calibration/calibrate_lake.R <Lake> [cores] [itermax]
#
# Parameters and ranges: PARS in glm4_objective.R (6 parameters, as ISIMIP3_LER's GLM
# calibration). Objective: RMSE of daily water temperature at the observed depths over the
# calibration half of the window (obs_wtemp_cal.csv); the validation half is only scored.
# DEoptim: population NP = 10 x 6 = 60 (the DEoptim manual's recommendation), strategy 2
# (local-to-best, the default), stops at itermax or when the best RMSE improves by
# < 0.1 % over 20 generations (reltol/steptol). Seed fixed -> reproducible.
# Each parallel worker runs GLM4 in its own copy of the lake's GLM folder.
#
# Writes calibration/results/<Lake>/: best.csv (parameters, Kw, RMSE cal/val, defaults'
# RMSE for reference), trace.csv (best RMSE per generation), deoptim.rds.

suppressMessages({ library(data.table); library(DEoptim); library(parallel) })
args <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)])))
source(file.path(here, "glm4_objective.R"))

a <- commandArgs(TRUE)
lake    <- a[1]
cores   <- if (length(a) >= 2) as.integer(a[2]) else max(1L, detectCores() - 1L)
itermax <- if (length(a) >= 3) as.integer(a[3]) else 100L
SEED <- 20261021

lake_dir <- file.path(here, "runs", lake)
glm_bin  <- normalizePath(file.path(dirname(here), "bin", "glm-run"))   # launcher from make_glm_launcher.sh
w   <- fread(file.path(lake_dir, "window.csv"))
rd  <- function(f) fread(file.path(lake_dir, f))[, datetime := as.POSIXct(datetime, tz = "UTC")]
obs_cal <- rd("obs_wtemp_cal.csv"); obs_val <- rd("obs_wtemp_val.csv")
out <- file.path(here, "results", lake); dir.create(out, recursive = TRUE, showWarnings = FALSE)
work <- file.path(Sys.getenv("TMPDIR", tempdir()), paste0("glm4cal_", lake, "_", Sys.getpid()))
dir.create(work, recursive = TRUE, showWarnings = FALSE)

# One private copy of the GLM folder per worker process.
objective <- function(p) {
  rd <- file.path(work, Sys.getpid())
  if (!dir.exists(rd)) {
    dir.create(rd)
    file.copy(list.files(file.path(lake_dir, "GLM"), full.names = TRUE), rd, recursive = TRUE)
  }
  glm4_rmse(p, rd, glm_bin, w$kw_base, obs_cal)
}

score <- function(p, o) {         # RMSE of a parameter set on an obs table (single run)
  rd <- file.path(work, "score"); unlink(rd, recursive = TRUE); dir.create(rd)
  file.copy(list.files(file.path(lake_dir, "GLM"), full.names = TRUE), rd, recursive = TRUE)
  glm4_rmse(p, rd, glm_bin, w$kw_base, o)
}

cat(sprintf("[%s] %s: %d cores, NP=60, itermax=%d, cal %s..%s (%d obs), val ..%s (%d obs)\n",
            format(Sys.time()), lake, cores, itermax, w$cal_start, w$val_start, nrow(obs_cal),
            w$end, nrow(obs_val)))
set.seed(SEED)
cl <- makeCluster(cores)
clusterExport(cl, c("here", "work", "lake_dir", "glm_bin", "w", "obs_cal"))
invisible(clusterEvalQ(cl, source(file.path(here, "glm4_objective.R"))))
t0 <- Sys.time()
res <- DEoptim(objective, lower = PARS$lower, upper = PARS$upper,
               # A supplied cluster is used whatever parallelType says; "parallel" would also
               # make DEoptim (2.2-8) call stopCluster(cl) on an undefined cl, which R
               # resolves to this script's global cl and closes it. Hence "none" + cluster.
               control = DEoptim.control(NP = 60, itermax = itermax, strategy = 2,
                                         reltol = 1e-3, steptol = 20, trace = 5,
                                         parallelType = "none", cluster = cl))
stopCluster(cl)
mins <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

best <- setNames(res$optim$bestmem, PARS$name)
row <- data.table(lake = lake, t(best), Kw = w$kw_base * best[["kw_scaling"]], kw_base = w$kw_base,
                  rmse_cal = res$optim$bestval, rmse_val = score(best, obs_val),
                  rmse_cal_default = score(PARS$glm_default, obs_cal),
                  rmse_val_default = score(PARS$glm_default, obs_val),
                  generations = res$optim$iter, evaluations = res$optim$nfeval, minutes = round(mins, 1),
                  n_obs_cal = nrow(obs_cal), n_obs_val = nrow(obs_val),
                  at_bound = paste(PARS$name[best <= PARS$lower * 1.01 | best >= PARS$upper * 0.99], collapse = ";"),
                  seed = SEED)
fwrite(row, file.path(out, "best.csv"))
fwrite(data.table(generation = seq_along(res$member$bestvalit), best_rmse = res$member$bestvalit),
       file.path(out, "trace.csv"))
saveRDS(res, file.path(out, "deoptim.rds"))
unlink(work, recursive = TRUE)
print(row)
