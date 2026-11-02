#!/bin/bash
# One task per record: 0 NOAA, 1..59 CMIP6. Preserve all existing outputs.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
WORKFLOW_ROOT="${WORKFLOW_ROOT:-${SLURM_SUBMIT_DIR:-$(cd -- "$SCRIPT_DIR/../.." && pwd -P)}}"
PYTHON_BIN="${PYTHON_BIN:-$HOME/.conda/envs/cmip6_equiv_depth_distributions/bin/python}"
OUTPUT_DIR="$WORKFLOW_ROOT/data/wave_significance/codex_segment_heff"
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
    cd "$WORKFLOW_ROOT"
    exec "$PYTHON_BIN" scripts/codex_segment_heff/run.py --task-index "$SLURM_ARRAY_TASK_ID"
fi
case "${1:-}" in
    pilot) tasks="0-1%2" ;;
    full)
        [[ -f "$OUTPUT_DIR/NOAA_OLR.observations.not_applicable.nc" &&
           -f "$OUTPUT_DIR/ACCESS-CM2.historical.r1i1p1f1.nc" ]] || {
            echo "Run and inspect NOAA and ACCESS-CM2 pilots first" >&2; exit 1;
        }
        tasks="2-59%5" ;;
    *) echo "Usage: bash scripts/codex_segment_heff/submit.sh pilot|full" >&2; exit 2 ;;
esac
[[ -x "$PYTHON_BIN" ]] || { echo "Missing Python: $PYTHON_BIN" >&2; exit 1; }
mkdir -p "$OUTPUT_DIR/jobs"
sbatch --parsable --account=bb1153 --partition=interactive --nodes=1 --ntasks=1 \
    --cpus-per-task=2 --mem=16G --time=01:00:00 --array="$tasks" \
    --job-name=segment_heff --chdir="$WORKFLOW_ROOT" \
    --output="$OUTPUT_DIR/jobs/segments.%A_%a.log" \
    --export="ALL,WORKFLOW_ROOT=$WORKFLOW_ROOT,PYTHON_BIN=$PYTHON_BIN" "$SCRIPT_DIR/submit.sh"
