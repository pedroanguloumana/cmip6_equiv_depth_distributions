#!/bin/bash
# Submit the 59 model/experiment records, with one audited pilot first.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
WORKFLOW_ROOT="${WORKFLOW_ROOT:-${SLURM_SUBMIT_DIR:-$(cd -- "$SCRIPT_DIR/../.." && pwd -P)}}"
PYTHON_BIN="${PYTHON_BIN:-$HOME/.conda/envs/cmip6_equiv_depth_distributions/bin/python}"
OUTPUT_DIR="$WORKFLOW_ROOT/data/wave_significance/codex_cmip_heff_comparison"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
    export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
    cd "$WORKFLOW_ROOT"
    exec "$PYTHON_BIN" scripts/codex_cmip_heff/run.py --task-index "$SLURM_ARRAY_TASK_ID"
fi

[[ -x "$PYTHON_BIN" ]] || { echo "Missing Python: $PYTHON_BIN" >&2; exit 1; }
case "${1:-}" in
    pilot) tasks="0" ;;
    full)
        [[ -f "$OUTPUT_DIR/ACCESS-CM2.historical.r1i1p1f1.event_energy.csv.gz" ]] || {
            echo "Run and audit the ACCESS-CM2 historical pilot first" >&2; exit 1;
        }
        [[ -f "$OUTPUT_DIR/ACCESS-CM2.historical.r1i1p1f1.spectral_power.csv.gz" ]] || {
            echo "Pilot spectral output is missing" >&2; exit 1;
        }
        tasks="1-58%5" ;;
    *) echo "Usage: bash scripts/codex_cmip_heff/submit.sh pilot|full" >&2; exit 2 ;;
esac
mkdir -p "$OUTPUT_DIR/jobs"
sbatch --parsable --account=bb1153 --partition=interactive --nodes=1 --ntasks=1 \
    --cpus-per-task=2 --mem=20G --time=02:00:00 --array="$tasks" \
    --job-name=codex_cmip_heff --chdir="$WORKFLOW_ROOT" \
    --output="$OUTPUT_DIR/jobs/model.%A_%a.log" \
    --export="ALL,WORKFLOW_ROOT=$WORKFLOW_ROOT,PYTHON_BIN=$PYTHON_BIN" "$SCRIPT_DIR/submit.sh"
