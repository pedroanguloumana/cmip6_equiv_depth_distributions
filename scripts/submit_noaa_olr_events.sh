#!/bin/bash
# Submit one Slurm array task per wave, following the DYAMOND worker/launcher layout.
#
# Usage on Levante, from the project root (no active conda shell is required):
#   bash scripts/submit_noaa_olr_events.sh --dry-run
#   bash scripts/submit_noaa_olr_events.sh
#   WAVES="Kelvin" MAX_EVENTS=3 bash scripts/submit_noaa_olr_events.sh
#
# This script is saved by Slurm in its spool directory. WORKFLOW_ROOT is therefore
# passed explicitly to each task, instead of deriving the project from that copy.
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
WORKFLOW_ROOT="${WORKFLOW_ROOT:-$(cd -- "$SCRIPT_DIR/.." && pwd -P)}"
WORKER="$WORKFLOW_ROOT/scripts/noaa_olr_events.py"
PYTHON_BIN="${PYTHON_BIN:-$HOME/.conda/envs/cmip6_equiv_depth_distributions/bin/python}"

# Scientific parameters are explicit dictionaries/constants in the Python worker.
# Dates and a pilot limit are shared by every selected wave in this submission.
WAVES="${WAVES:-Kelvin n1ERW n0EIG MRG n1WIG}"
START_DATE="${START_DATE:-1990-01-01}"
END_DATE="${END_DATE:-2010-12-31}"
MAX_EVENTS="${MAX_EVENTS:-}"
DETECTION_SIGMA="${DETECTION_SIGMA:-1.5}"
FIT_SUPPORT="${FIT_SUPPORT:-window}"
case "$DETECTION_SIGMA" in 1|1.0|1.5|2|2.0) ;; *) echo "DETECTION_SIGMA must be 1, 1.5 or 2" >&2; exit 2 ;; esac
case "$FIT_SUPPORT" in window|component) ;; *) echo "FIT_SUPPORT must be window or component" >&2; exit 2 ;; esac
OVERWRITE="${OVERWRITE:-0}"
DRY_RUN="${DRY_RUN:-0}"

# Full-record, 25-event pilots (array 28552504) peaked at 10.1--10.2 GiB and
# projected about 6--32 minutes per complete wave catalog. Allow memory/runtime
# headroom: 14 GiB and one hour. The interactive partition permits 1940 MiB/CPU,
# so eight allocated CPUs accommodate this memory request (shared needs 16).
# BLAS remains single-threaded because these are many small, narrow least squares.
ACCOUNT="${ACCOUNT:-bb1153}"
PARTITION="${PARTITION:-interactive}"
JOB_MEM="${JOB_MEM:-14G}"
JOB_CPUS="${JOB_CPUS:-8}"
JOB_TIME="${JOB_TIME:-01:00:00}"
MAX_CONCURRENT="${MAX_CONCURRENT:-4}"

if [[ $# -gt 1 || ( $# -eq 1 && "$1" != "--dry-run" ) ]]; then
    echo "Usage: bash scripts/submit_noaa_olr_events.sh [--dry-run]" >&2
    exit 2
fi
if [[ "${1:-}" == "--dry-run" ]]; then DRY_RUN=1; fi

[[ -x "$PYTHON_BIN" ]] || { echo "Python not executable: $PYTHON_BIN" >&2; exit 1; }
[[ -f "$WORKER" ]] || { echo "Worker not found: $WORKER" >&2; exit 1; }
cd "$WORKFLOW_ROOT"

############################### WORKER MODE #################################
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    : "${TASKLIST:?The launcher must provide TASKLIST}"
    wave="$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$TASKLIST")"
    [[ -n "$wave" ]] || { echo "No wave at array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
    args=(--wave "$wave" --start "$START_DATE" --end "$END_DATE"
        --detection-sigma "$DETECTION_SIGMA" --fit-support "$FIT_SUPPORT")
    if [[ -n "$MAX_EVENTS" ]]; then args+=(--max-events "$MAX_EVENTS"); fi
    if [[ "$OVERWRITE" == 1 ]]; then args+=(--overwrite); fi

    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
    echo "host=$(hostname) wave=$wave job=${SLURM_JOB_ID:-} task=$SLURM_ARRAY_TASK_ID"
    if [[ "$DRY_RUN" == 1 ]]; then
        printf 'Would run: '; printf '%q ' "$PYTHON_BIN" "$WORKER" "${args[@]}"; printf '\n'
        exit 0
    fi
    exec "$PYTHON_BIN" "$WORKER" "${args[@]}"
fi

############################ ORCHESTRATE MODE ################################
read -r -a wave_list <<< "$WAVES"
[[ ${#wave_list[@]} -gt 0 ]] || { echo "WAVES is empty" >&2; exit 2; }
seen=" "
for wave in "${wave_list[@]}"; do
    case "$wave" in Kelvin|n1ERW|n0EIG|MRG|n1WIG) ;; *) echo "Unknown wave: $wave" >&2; exit 2 ;; esac
    [[ "$seen" != *" $wave "* ]] || { echo "Duplicate wave: $wave" >&2; exit 2; }
    seen+="$wave "
done
[[ "$MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]] || { echo "MAX_CONCURRENT must be positive" >&2; exit 2; }

# Get the output entry point from the installed project, not an archive path.
JOB_ROOT="$("$PYTHON_BIN" -c 'from src.paths import NOAA_OLR_JOB_DIR; print(NOAA_OLR_JOB_DIR)')"
LOGDIR="$JOB_ROOT/logs"
LISTDIR="$JOB_ROOT/tasklists"
TASKLIST="$LISTDIR/waves_$(date +%Y%m%d_%H%M%S)_$$.txt"
SELF="$SCRIPT_DIR/$(basename -- "${BASH_SOURCE[0]}")"

command=(
    sbatch --parsable --job-name=NOAA_OLR_events
    "--account=$ACCOUNT" "--partition=$PARTITION"
    --nodes=1 --ntasks=1 "--cpus-per-task=$JOB_CPUS" "--mem=$JOB_MEM" "--time=$JOB_TIME"
    "--array=0-$((${#wave_list[@]}-1))%$MAX_CONCURRENT"
    "--chdir=$WORKFLOW_ROOT"
    "--output=$LOGDIR/NOAA_OLR.%A_%a.log" "--error=$LOGDIR/NOAA_OLR.%A_%a.log"
    "--export=ALL,WORKFLOW_ROOT=$WORKFLOW_ROOT,PYTHON_BIN=$PYTHON_BIN,TASKLIST=$TASKLIST,START_DATE=$START_DATE,END_DATE=$END_DATE,DETECTION_SIGMA=$DETECTION_SIGMA,FIT_SUPPORT=$FIT_SUPPORT,MAX_EVENTS=$MAX_EVENTS,OVERWRITE=$OVERWRITE,DRY_RUN=0"
    "$SELF"
)

echo "Waves: $WAVES"
echo "Record: $START_DATE through $END_DATE; activity cutoff: ${DETECTION_SIGMA} sigma; fit support: $FIT_SUPPORT; pilot limit: ${MAX_EVENTS:-none} events/wave"
echo "Resources/task: partition=$PARTITION account=$ACCOUNT mem=$JOB_MEM cpus=$JOB_CPUS time=$JOB_TIME"
if [[ "$DRY_RUN" == 1 ]]; then
    echo "Dry run: no task list, logs, jobs, or catalogs are created."
    printf '%q ' "${command[@]}"; printf '\n'
    exit 0
fi

mkdir -p "$LOGDIR" "$LISTDIR"
printf '%s\n' "${wave_list[@]}" > "$TASKLIST"
job_id="$("${command[@]}")"
echo "Submitted array $job_id; task order: $WAVES"
echo "Monitor: squeue -j ${job_id%%;*} -r"
echo "Logs: $LOGDIR/NOAA_OLR.${job_id%%;*}_<task>.log"
