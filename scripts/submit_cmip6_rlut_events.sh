#!/bin/bash
# One Slurm task per model/experiment/member/wave. Run from the project root:
#   bash scripts/submit_cmip6_rlut_events.sh --dry-run
#   MODELS="CanESM5 HadGEM3-GC31-LL" MAX_EVENTS=3 bash scripts/submit_cmip6_rlut_events.sh
#   bash scripts/submit_cmip6_rlut_events.sh
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
WORKFLOW_ROOT="${WORKFLOW_ROOT:-$(cd -- "$SCRIPT_DIR/.." && pwd -P)}"
WORKER="$WORKFLOW_ROOT/scripts/cmip6_rlut_events.py"
PYTHON_BIN="${PYTHON_BIN:-$HOME/.conda/envs/cmip6_equiv_depth_distributions/bin/python}"

# Empty MODELS discovers all available pairs through data/cmip6_rlut_noaa1deg.
# Full native-calendar years avoid imposing Gregorian month lengths on models.
MODELS="${MODELS:-}"
EXPERIMENTS="${EXPERIMENTS:-historical amip}"
WAVES="${WAVES:-Kelvin n1ERW n0EIG MRG n1WIG}"
START_YEAR="${START_YEAR:-1990}"
END_YEAR="${END_YEAR:-2010}"
MAX_EVENTS="${MAX_EVENTS:-}"
DETECTION_SIGMA="${DETECTION_SIGMA:-1.5}"
FIT_SUPPORT="${FIT_SUPPORT:-window}"
case "$DETECTION_SIGMA" in 1|1.0|1.5|2|2.0) ;; *) echo "DETECTION_SIGMA must be 1, 1.5 or 2" >&2; exit 2 ;; esac
case "$FIT_SUPPORT" in window|component) ;; *) echo "FIT_SUPPORT must be window or component" >&2; exit 2 ;; esac
OVERWRITE="${OVERWRITE:-0}"
DRY_RUN="${DRY_RUN:-0}"

# The full 1990--2010, 1-degree NOAA pilots peaked at 10.2 GiB. The new CMIP
# inputs have the same selected time/latitude/longitude dimensions. Allow 14 GiB
# and eight CPUs on interactive (1940 MiB/CPU). The 1990--2010 CMIP pilot's ER
# extrapolation was 47 minutes; allow two hours for event-to-event variability.
# Library threads remain one; extra CPUs supply the required memory allowance.
ACCOUNT="${ACCOUNT:-bb1153}"
PARTITION="${PARTITION:-interactive}"
JOB_MEM="${JOB_MEM:-14G}"
JOB_CPUS="${JOB_CPUS:-8}"
JOB_TIME="${JOB_TIME:-02:00:00}"
MAX_CONCURRENT="${MAX_CONCURRENT:-4}"

if [[ $# -gt 1 || ( $# -eq 1 && "$1" != "--dry-run" ) ]]; then
    echo "Usage: bash scripts/submit_cmip6_rlut_events.sh [--dry-run]" >&2
    exit 2
fi
if [[ "${1:-}" == "--dry-run" ]]; then DRY_RUN=1; fi
[[ -x "$PYTHON_BIN" ]] || { echo "Python not executable: $PYTHON_BIN" >&2; exit 1; }
[[ -f "$WORKER" ]] || { echo "Worker not found: $WORKER" >&2; exit 1; }
cd "$WORKFLOW_ROOT"

############################### WORKER MODE #################################
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    : "${TASKLIST:?The launcher must provide TASKLIST}"
    IFS=$'\t' read -r model experiment member wave <<< "$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$TASKLIST")"
    [[ -n "$model" && -n "$experiment" && -n "$member" && -n "$wave" ]] || { echo "Missing simulation/wave at array index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
    args=(--model "$model" --experiment "$experiment" --member "$member" --wave "$wave"
        --start-year "$START_YEAR" --end-year "$END_YEAR"
        --detection-sigma "$DETECTION_SIGMA" --fit-support "$FIT_SUPPORT")
    if [[ -n "$MAX_EVENTS" ]]; then args+=(--max-events "$MAX_EVENTS"); fi
    if [[ "$OVERWRITE" == 1 ]]; then args+=(--overwrite); fi
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
    export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
    echo "host=$(hostname) model=$model experiment=$experiment member=$member wave=$wave job=${SLURM_JOB_ID:-} task=$SLURM_ARRAY_TASK_ID"
    if [[ "$DRY_RUN" == 1 ]]; then
        printf 'Would run: '; printf '%q ' "$PYTHON_BIN" "$WORKER" "${args[@]}"; printf '\n'
        exit 0
    fi
    exec "$PYTHON_BIN" "$WORKER" "${args[@]}"
fi

############################ ORCHESTRATE MODE ################################
# Discovery uses filenames only; every worker verifies the simulation against
# the NetCDF header. Available model/experiment pairs are not a Cartesian product.
available_simulations="$("$PYTHON_BIN" "$WORKER" --list-simulations)"
available_models="$(printf '%s\n' "$available_simulations" | cut -f1 | sort -u)"
available_models="${available_models//$'\n'/ }"
if [[ -z "$MODELS" ]]; then MODELS="$available_models"; fi
MODELS="${MODELS//$'\n'/ }"
read -r -a model_list <<< "$MODELS"
read -r -a experiment_list <<< "$EXPERIMENTS"
read -r -a wave_list <<< "$WAVES"
[[ ${#model_list[@]} -gt 0 && ${#experiment_list[@]} -gt 0 && ${#wave_list[@]} -gt 0 ]] || { echo "Empty model/experiment/wave selection" >&2; exit 2; }
for experiment in "${experiment_list[@]}"; do
    case "$experiment" in historical|amip) ;; *) echo "Unknown experiment: $experiment" >&2; exit 2 ;; esac
done
seen=" "
for model in "${model_list[@]}"; do
    [[ "$model" =~ ^[A-Za-z0-9_-]+$ ]] || { echo "Invalid model ID: $model" >&2; exit 2; }
    [[ " $available_models " == *" $model "* ]] || { echo "Model not found: $model" >&2; exit 2; }
    [[ "$seen" != *" $model "* ]] || { echo "Duplicate model: $model" >&2; exit 2; }
    seen+="$model "
done
selected_simulations=()
while IFS=$'\t' read -r model experiment member; do
    if [[ " $MODELS " == *" $model "* && " $EXPERIMENTS " == *" $experiment "* ]]; then
        selected_simulations+=("$model"$'\t'"$experiment"$'\t'"$member")
    fi
done <<< "$available_simulations"
[[ ${#selected_simulations[@]} -gt 0 ]] || { echo "No available model/experiment pairs match this selection" >&2; exit 2; }
seen=" "
for wave in "${wave_list[@]}"; do
    case "$wave" in Kelvin|n1ERW|n0EIG|MRG|n1WIG) ;; *) echo "Unknown wave: $wave" >&2; exit 2 ;; esac
    [[ "$seen" != *" $wave "* ]] || { echo "Duplicate wave: $wave" >&2; exit 2; }
    seen+="$wave "
done
[[ "$START_YEAR" =~ ^[1-9][0-9]{3}$ && "$END_YEAR" =~ ^[1-9][0-9]{3}$ && "$START_YEAR" -le "$END_YEAR" ]] || { echo "Require increasing four-digit years" >&2; exit 2; }
[[ "$MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]] || { echo "MAX_CONCURRENT must be positive" >&2; exit 2; }
if [[ -n "$MAX_EVENTS" && ! "$MAX_EVENTS" =~ ^[1-9][0-9]*$ ]]; then
    echo "MAX_EVENTS must be positive" >&2; exit 2
fi

JOB_ROOT="$("$PYTHON_BIN" -c 'from src.paths import CMIP6_EVENT_JOB_DIR; print(CMIP6_EVENT_JOB_DIR)')"
LOGDIR="$JOB_ROOT/logs"
LISTDIR="$JOB_ROOT/tasklists"
TASKLIST="$LISTDIR/simulations_waves_$(date +%Y%m%d_%H%M%S)_$$.tsv"
SELF="$SCRIPT_DIR/$(basename -- "${BASH_SOURCE[0]}")"
task_count=$((${#selected_simulations[@]} * ${#wave_list[@]}))
command=(
    sbatch --parsable --job-name=CMIP6_rlut_events
    "--account=$ACCOUNT" "--partition=$PARTITION"
    --nodes=1 --ntasks=1 "--cpus-per-task=$JOB_CPUS" "--mem=$JOB_MEM" "--time=$JOB_TIME"
    "--array=0-$((task_count-1))%$MAX_CONCURRENT"
    "--chdir=$WORKFLOW_ROOT"
    "--output=$LOGDIR/CMIP6_rlut.%A_%a.log" "--error=$LOGDIR/CMIP6_rlut.%A_%a.log"
    "--export=ALL,WORKFLOW_ROOT=$WORKFLOW_ROOT,PYTHON_BIN=$PYTHON_BIN,TASKLIST=$TASKLIST,START_YEAR=$START_YEAR,END_YEAR=$END_YEAR,DETECTION_SIGMA=$DETECTION_SIGMA,FIT_SUPPORT=$FIT_SUPPORT,MAX_EVENTS=$MAX_EVENTS,OVERWRITE=$OVERWRITE,DRY_RUN=0"
    "$SELF"
)

echo "Models selected: $MODELS"
echo "Experiments: $EXPERIMENTS; available simulation pairs selected=${#selected_simulations[@]}"
echo "Waves: $WAVES; tasks=$task_count; maximum concurrent=$MAX_CONCURRENT"
echo "Native-calendar years: $START_YEAR through $END_YEAR; activity cutoff: ${DETECTION_SIGMA} sigma; fit support: $FIT_SUPPORT; pilot limit: ${MAX_EVENTS:-none} events/task"
echo "Resources/task: partition=$PARTITION account=$ACCOUNT mem=$JOB_MEM cpus=$JOB_CPUS time=$JOB_TIME"
if [[ "$DRY_RUN" == 1 ]]; then
    echo "Dry run: no task list, logs, jobs, or catalogs are created."
    printf '%q ' "${command[@]}"; printf '\n'
    exit 0
fi
mkdir -p "$LOGDIR" "$LISTDIR"
for simulation in "${selected_simulations[@]}"; do
    for wave in "${wave_list[@]}"; do printf '%s\t%s\n' "$simulation" "$wave"; done
done > "$TASKLIST"
job_id="$("${command[@]}")"
echo "Submitted array $job_id; manifest: $TASKLIST"
echo "Monitor: squeue -j ${job_id%%;*} -r"
echo "Logs: $LOGDIR/CMIP6_rlut.${job_id%%;*}_<task>.log"
