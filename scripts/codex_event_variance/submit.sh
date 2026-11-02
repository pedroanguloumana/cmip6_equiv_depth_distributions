#!/bin/bash
# Submit the 999 remaining held-out NOAA null records in disjoint 25-record blocks.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
WORKFLOW_ROOT="${WORKFLOW_ROOT:-$(cd -- "$SCRIPT_DIR/../.." && pwd -P)}"
PYTHON_BIN="${PYTHON_BIN:-$HOME/.conda/envs/cmip6_equiv_depth_distributions/bin/python}"
OUTPUT_DIR="$WORKFLOW_ROOT/data/wave_significance/gaussian_calibration_wk99_121/NOAA_OLR.observations.not_applicable.1990-2010.b898b1ae4307/codex_event_variance"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    start=$((1001 + 25 * SLURM_ARRAY_TASK_ID))
    stop=$((start + 25))
    if ((stop > 2000)); then stop=2000; fi
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    export PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
    cd "$WORKFLOW_ROOT"
    exec "$PYTHON_BIN" scripts/codex_event_variance/run.py --start "$start" --stop "$stop"
fi

[[ -x "$PYTHON_BIN" ]] || { echo "Missing Python: $PYTHON_BIN" >&2; exit 1; }
[[ -f "$OUTPUT_DIR/null.1000-1001.csv.gz" ]] || {
    echo "Run and validate null realization 1000 first" >&2; exit 1;
}
mkdir -p "$OUTPUT_DIR/jobs"
sbatch --parsable --account=bb1153 --partition=interactive --nodes=1 --ntasks=1 \
    --cpus-per-task=1 --mem=8G --time=00:25:00 --array=0-39%8 \
    --job-name=codex_event_variance --chdir="$WORKFLOW_ROOT" \
    --output="$OUTPUT_DIR/jobs/null.%A_%a.log" \
    --export="ALL,WORKFLOW_ROOT=$WORKFLOW_ROOT,PYTHON_BIN=$PYTHON_BIN" "$SCRIPT_DIR/submit.sh"
