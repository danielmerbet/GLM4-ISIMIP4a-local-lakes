#!/bin/bash
# GLM4 hydrodynamic calibration on Hydra: one array task per lake (DEoptim, 16 cores).
#
#   cd $VSC_DATA/GLM4-AED3
#   sbatch --array=1-$(wc -l < calibration/lakes_to_calibrate.txt)%16 calibration/hydra_calibrate.slurm
# or per runtime tier (lakes_tier{A,B,C}.txt from sweep_defaults.csv, see README):
#   sbatch --export=ALL,LAKE_LIST=calibration/lakes_tierA.txt --array=1-55%10 --time=06:00:00 calibration/hydra_calibrate.slurm
#
# Lakes: calibration/lakes_to_calibrate.txt (one per line; line N = array task N).
# ITERMAX (env, default 100) caps DEoptim generations; a short test: --export=ALL,ITERMAX=2
# Skips lakes that already have results/<Lake>/best.csv, so resubmitting resumes.
# Time limit: the slowest lake (Kivu, 485 m) took 81 s per model run on the laptop under
# load; DEoptim stops at 100 generations x 60 = 6000 runs at most -> 6000 x 81 / 16 ~ 8.5 h.
#
#SBATCH --job-name=glm4-cal
#SBATCH --account=bvo00012
#SBATCH --partition=zen4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G             # tier A; deep lakes need more: tier B --mem=96G, Kivu --mem=192G
#SBATCH --time=12:00:00
#SBATCH --output=calibration/logs/%x-%A_%a.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
LAKE_LIST="${LAKE_LIST:-calibration/lakes_to_calibrate.txt}"
LAKE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LAKE_LIST")
[ -n "$LAKE" ] || { echo "no lake on line $SLURM_ARRAY_TASK_ID"; exit 1; }
if [ -f "calibration/results/$LAKE/best.csv" ]; then echo "$LAKE: already calibrated, skipping"; exit 0; fi

# R (+ ncdf4 built against netCDF 4.9.3); GLM itself runs through bin/glm-run, which pins
# its own build libraries (make_glm_launcher.sh), so the two toolchains don't mix.
module purge
module load R/4.5.1-gfbf-2025a netCDF/4.9.3-gompi-2025a
export R_LIBS_USER="$VSC_DATA/LERWQ/Rlibs"
export TZ=UTC
export TMPDIR="${TMPDIR:-/tmp}"          # node-local scratch for the per-worker GLM copies

echo "[$(date -u +%FT%TZ)] $LAKE on $(hostname), $SLURM_CPUS_PER_TASK cores"
Rscript calibration/calibrate_lake.R "$LAKE" "$SLURM_CPUS_PER_TASK" "${ITERMAX:-100}"
rc=$?            # keep R's status: a "$(date)" in the echo below would reset $? to 0
echo "[$(date -u +%FT%TZ)] $LAKE finished, exit $rc"
[ -f "calibration/results/$LAKE/best.csv" ] || { echo "$LAKE: no best.csv written"; rc=${rc/#0/1}; }
exit $rc
