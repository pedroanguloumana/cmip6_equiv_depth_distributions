#!/bin/bash
# One independent dataset per array task; keep BLAS within the Slurm allocation.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
WORKFLOW_ROOT="${SLURM_SUBMIT_DIR:-$(cd -- "$SCRIPT_DIR/../.." && pwd -P)}"
PYTHON_BIN="${PYTHON_BIN:-$HOME/.conda/envs/cmip6_equiv_depth_distributions/bin/python}"
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 MKL_NUM_THREADS=2
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    cd "$WORKFLOW_ROOT"
    exec "$PYTHON_BIN" scripts/codex_olr_regression/run.py --task-index "$SLURM_ARRAY_TASK_ID" "$@"
fi
OUTPUT_DIR="$WORKFLOW_ROOT/data/olr_wave_regression/jobs"
mkdir -p "$OUTPUT_DIR"
sbatch --parsable --account=bb1153 --partition=interactive --nodes=1 --ntasks=1 \
    --cpus-per-task=2 --mem=24G --time=00:30:00 --array=0-2%3 \
    --job-name=olr_wave_maps --chdir="$WORKFLOW_ROOT" \
    --output="$OUTPUT_DIR/maps.%A_%a.log" \
    --export="ALL,PYTHON_BIN=$PYTHON_BIN" "$SCRIPT_DIR/submit.sh" "$@"
