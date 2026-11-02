#!/bin/bash
# Submit disjoint held-out NOAA Gaussian-null catalog blocks from the project root.
# Use --dry-run to print the Slurm command without creating jobs or outputs.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
WORKFLOW_ROOT="${WORKFLOW_ROOT:-$(cd -- "$SCRIPT_DIR/.." && pwd -P)}"
PYTHON_BIN="${PYTHON_BIN:-$HOME/.conda/envs/cmip6_equiv_depth_distributions/bin/python}"
MANIFEST="${MANIFEST:-$WORKFLOW_ROOT/data/wave_significance/gaussian_calibration_wk99_121/NOAA_OLR.observations.not_applicable.1990-2010.b898b1ae4307/manifest.json}"
START_REALIZATION="${START_REALIZATION:-1000}"
STOP_REALIZATION="${STOP_REALIZATION:-2000}"
BLOCK_SIZE="${BLOCK_SIZE:-200}"
MAX_PARALLEL="${MAX_PARALLEL:-16}"
ACCOUNT="${ACCOUNT:-bb1153}"
PARTITION="${PARTITION:-interactive}"
JOB_CPUS="${JOB_CPUS:-24}"
JOB_MEM="${JOB_MEM:-40G}"
JOB_TIME="${JOB_TIME:-08:00:00}"
MAX_CONCURRENT="${MAX_CONCURRENT:-4}"

if [[ $# -gt 1 || ( $# -eq 1 && "$1" != "--dry-run" ) ]]; then
    echo "Usage: bash scripts/submit_noaa_background_fit.sh [--dry-run]" >&2
    exit 2
fi
[[ -x "$PYTHON_BIN" && -f "$MANIFEST" ]] || { echo "Missing Python or manifest" >&2; exit 1; }
for value in "$START_REALIZATION" "$STOP_REALIZATION" "$BLOCK_SIZE" "$MAX_PARALLEL" "$MAX_CONCURRENT"; do
    [[ "$value" =~ ^[1-9][0-9]*$ ]] || { echo "Indices, block size and concurrency must be positive integers" >&2; exit 2; }
done
[[ "$START_REALIZATION" -lt "$STOP_REALIZATION" ]] || { echo "Require START_REALIZATION < STOP_REALIZATION" >&2; exit 2; }
cd "$WORKFLOW_ROOT"

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    start=$((START_REALIZATION + SLURM_ARRAY_TASK_ID * BLOCK_SIZE))
    stop=$((start + BLOCK_SIZE))
    if [[ "$stop" -gt "$STOP_REALIZATION" ]]; then stop="$STOP_REALIZATION"; fi
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
    echo "host=$(hostname) job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID realizations=$start:$stop"
    exec "$PYTHON_BIN" scripts/run_noaa_background_fit_batch.py \
        --manifest "$MANIFEST" --start "$start" --stop "$stop" --max-parallel "$MAX_PARALLEL"
fi

blocks=$(((STOP_REALIZATION - START_REALIZATION + BLOCK_SIZE - 1) / BLOCK_SIZE))
log_dir="$(dirname "$MANIFEST")/jobs"
command=(
    sbatch --parsable --job-name=NOAA_null_component_fits
    "--account=$ACCOUNT" "--partition=$PARTITION"
    --nodes=1 --ntasks=1 "--cpus-per-task=$JOB_CPUS" "--mem=$JOB_MEM" "--time=$JOB_TIME"
    "--array=0-$((blocks-1))%$MAX_CONCURRENT" "--chdir=$WORKFLOW_ROOT"
    "--output=$log_dir/component_fits.%A_%a.log" "--error=$log_dir/component_fits.%A_%a.log"
    "--export=ALL,WORKFLOW_ROOT=$WORKFLOW_ROOT,PYTHON_BIN=$PYTHON_BIN,MANIFEST=$MANIFEST,START_REALIZATION=$START_REALIZATION,STOP_REALIZATION=$STOP_REALIZATION,BLOCK_SIZE=$BLOCK_SIZE,MAX_PARALLEL=$MAX_PARALLEL"
    "$SCRIPT_DIR/$(basename -- "${BASH_SOURCE[0]}")"
)
echo "Validation realizations: $START_REALIZATION:$STOP_REALIZATION; blocks=$blocks; five waves/block"
echo "Maximum parallel workers per block: $MAX_PARALLEL; concurrent blocks: $MAX_CONCURRENT"
echo "Resources/block: partition=$PARTITION cpus=$JOB_CPUS mem=$JOB_MEM time=$JOB_TIME"
if [[ "${1:-}" == "--dry-run" ]]; then
    printf '%q ' "${command[@]}"; printf '\n'
    exit 0
fi
mkdir -p "$log_dir"
job_id="$("${command[@]}")"
echo "Submitted array $job_id"
echo "Monitor: squeue -j ${job_id%%;*} -r"
echo "Logs: $log_dir/component_fits.${job_id%%;*}_<task>.log"
