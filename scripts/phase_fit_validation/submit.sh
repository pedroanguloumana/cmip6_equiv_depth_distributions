#!/usr/bin/env bash
# Submit the 15 source/experiment/wave comparisons from the project root.
set -euo pipefail

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    # Slurm runs a spool copy of this script; its dirname is not the project.
    ROOT="$PWD"
    HERE="$ROOT/scripts/phase_fit_validation"
else
    HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
    ROOT="$(cd -- "$HERE/../.." && pwd -P)"
fi
PYTHON_BIN="${PYTHON_BIN:-$HOME/.conda/envs/cmip6_equiv_depth_distributions/bin/python}"
[[ -x "$PYTHON_BIN" ]] || { echo "Python not executable: $PYTHON_BIN" >&2; exit 1; }
[[ -f "$HERE/cases.tsv" ]] || { echo "Missing cases.tsv" >&2; exit 1; }

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    IFS=$'\t' read -r source experiment wave <<< "$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$HERE/cases.tsv")"
    [[ -n "$source" && -n "$experiment" && -n "$wave" ]] || { echo "Missing array case" >&2; exit 1; }
    cd "$ROOT"
    export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
    exec "$PYTHON_BIN" "$HERE/run.py" --source "$source" --experiment "$experiment" \
        --wave "$wave" --backgrounds 6
fi

[[ $# -eq 0 ]] || { echo "Usage: bash scripts/phase_fit_validation/submit.sh" >&2; exit 2; }
mkdir -p "$HERE/logs"
cd "$ROOT"
sbatch --parsable --job-name=phase_fit_validation --account=bb1153 --partition=interactive \
    --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=14G --time=02:00:00 \
    --array=0-14%5 --chdir="$ROOT" \
    --output="$HERE/logs/%A_%a.log" --error="$HERE/logs/%A_%a.log" \
    "$HERE/submit.sh"
