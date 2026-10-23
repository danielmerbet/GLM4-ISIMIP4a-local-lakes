#!/bin/bash
# R packages the GLM4 calibration needs on Hydra, installed into the same library as the
# LakeEnsemblR/LakeEnsemblR.WQ stack ($VSC_DATA/LERWQ/Rlibs, ~/LERWQ_setup on Hydra).
# Run once on the LOGIN node (needs internet). Only missing CRAN packages are installed;
# glmtools (GitHub) is already there from the LERWQ install.
#   bash calibration/hydra_install_r_deps.sh
set -euo pipefail
module purge
module load R/4.5.1-gfbf-2025a netCDF/4.9.3-gompi-2025a
export R_LIBS_USER="$VSC_DATA/LERWQ/Rlibs"
Rscript -e '
options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 4)
need <- c("data.table", "DEoptim", "ncdf4", "glmtools")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
cran <- setdiff(miss, "glmtools")
if (length(cran)) install.packages(cran)
if ("glmtools" %in% miss) stop("glmtools missing: install it with ~/LERWQ_setup/02_install_r_packages.sh")
for (p in need) cat(sprintf("%-11s %s\n", p, as.character(packageVersion(p))))'
