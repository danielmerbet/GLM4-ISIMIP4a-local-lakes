#!/bin/bash
# ISIMIP4a local-lake files on Hydra: one array task per lake (line N+1 of local_lakes/lake_sites.csv).
#
#   cd $VSC_DATA/GLM4-AED3
#   sbatch --array=1-77 local_lakes/hydra_local_lakes.slurm          # or a subset, e.g. --array=5,9
#
# Per lake: setup_run.py (calibrated parameters, 1881-1900 spin-up = 1901-1920 repeated, run to
# 2026-01-01) -> GLM -> write_isimip.py (13 daily variables, 1901-2025).
#   run folder + GLM's full output: $VSC_SCRATCH_VO_USER/GLM4-AED3/local_lakes/runs/<Lake> (kept, so a
#       change of the strat definition or of the metadata only needs write_isimip.py again; ~160 GB
#       for all 77 lakes)
#   NetCDF files: $VSC_DATA_VO_USER/GLM4-AED3/local_lakes/isimip4a/<Lake>/ + SHA256SUMS
# Resumable: a lake with its SHA256SUMS already written is skipped (FORCE=1 redoes it).
# WRITE_ONLY=1 (with FORCE=1): rewrite the NetCDF files from the kept GLM output, no new GLM run.
# Lakes without calibration/results/<Lake>/best.csv stop with an error (calibration not finished).
#
#SBATCH --job-name=glm4-isimip
#SBATCH --account=bvo00012
#SBATCH --partition=zen4
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=local_lakes/logs/%x-%A_%a.out
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
LAKE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" local_lakes/lake_sites.csv | cut -d, -f1)
[ -n "$LAKE" ] || { echo "no lake on line $((SLURM_ARRAY_TASK_ID + 1))"; exit 1; }
RUN=${VSC_SCRATCH_VO_USER:?}/GLM4-AED3/local_lakes/runs/$LAKE
OUT=${VSC_DATA_VO_USER:?}/GLM4-AED3/local_lakes/isimip4a/$LAKE
if [ -f "$OUT/SHA256SUMS" ] && [ "${FORCE:-0}" != 1 ]; then echo "$LAKE: done already"; exit 0; fi
[ -f "calibration/results/$LAKE/best.csv" ] || { echo "$LAKE: no best.csv (calibration not finished)"; exit 1; }

module purge
module load netcdf4-python/1.7.2-foss-2025a
export TZ=UTC OMP_NUM_THREADS=1

echo "[$(date -u +%FT%TZ)] $LAKE on $(hostname)"
if [ "${WRITE_ONLY:-0}" = 1 ]; then
  [ -f "$RUN/output/output.nc" ] || { echo "$LAKE: WRITE_ONLY but no GLM output"; exit 1; }
else
  python3 local_lakes/setup_run.py "$LAKE" "$RUN" || exit 1
  t0=$(date +%s)
  (cd "$RUN" && "$SLURM_SUBMIT_DIR/bin/glm-run" --nml glm3.nml > glm.log 2>&1)
  rc=$?
  echo "[$(date -u +%FT%TZ)] GLM exit $rc after $(( $(date +%s) - t0 )) s; last lines of glm.log:"
  tail -3 "$RUN/glm.log"
  [ $rc -eq 0 ] && [ -f "$RUN/output/output.nc" ] || { echo "$LAKE: GLM failed"; exit 1; }
fi

rm -rf "$OUT"; mkdir -p "$OUT"
python3 local_lakes/write_isimip.py "$RUN" "$LAKE" "$OUT" || exit 1
cp "$RUN/glm3.nml" "$OUT/glm3.nml.used"
(cd "$OUT" && sha256sum *.nc glm3.nml.used > SHA256SUMS)
echo "[$(date -u +%FT%TZ)] $LAKE: $(ls "$OUT"/*.nc | wc -l) files, $(du -sh "$OUT" | cut -f1)"
