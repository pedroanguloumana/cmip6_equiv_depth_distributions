#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
WORKFLOW_ROOT="${WORKFLOW_ROOT:-${SLURM_SUBMIT_DIR:-$(cd -- "$SCRIPT_DIR/../.." && pwd -P)}}"
PYTHON_BIN="${PYTHON_BIN:-$HOME/.conda/envs/cmip6_equiv_depth_distributions/bin/python}"
OUTPUT_DIR="$WORKFLOW_ROOT/data/wave_significance/codex_kelvin_komega/variance_bins_integer_centered"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
    cd "$WORKFLOW_ROOT"
    exec "$PYTHON_BIN" scripts/codex_kelvin_komega/variance.py --task-index "$SLURM_ARRAY_TASK_ID"
fi

[[ -x "$PYTHON_BIN" ]] || { echo "Missing Python: $PYTHON_BIN" >&2; exit 1; }
case "${1:-}" in
    pilot) tasks="0-1" ;;
    full)
        [[ -f "$OUTPUT_DIR/NOAA_OLR.observations.not_applicable.csv.gz" ]] || {
            echo "Run and audit the NOAA pilot first" >&2; exit 1;
        }
        tasks="2-59%5" ;;
    *) echo "Usage: bash scripts/codex_kelvin_komega/submit_variance.sh pilot|full" >&2; exit 2 ;;
esac
mkdir -p "$OUTPUT_DIR/jobs"
sbatch --parsable --account=bb1153 --partition=interactive --nodes=1 --ntasks=1 \
    --cpus-per-task=2 --mem=20G --time=01:00:00 --array="$tasks" \
    --job-name=kelvin_komega_variance --chdir="$WORKFLOW_ROOT" \
    --output="$OUTPUT_DIR/jobs/variance.%A_%a.log" \
    --export="ALL,WORKFLOW_ROOT=$WORKFLOW_ROOT,PYTHON_BIN=$PYTHON_BIN" "$SCRIPT_DIR/submit_variance.sh"
