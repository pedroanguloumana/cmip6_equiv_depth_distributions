#!/usr/bin/env bash
# Full-record Gaussian significance calibration, with dependent Slurm stages.
# Examples from project root:
#   bash scripts/submit_gaussian_event_calibration.sh --dry-run
#   bash scripts/submit_gaussian_event_calibration.sh
#   MODEL=ACCESS-CM2 EXPERIMENT=historical bash scripts/submit_gaussian_event_calibration.sh
set -euo pipefail

ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
PYTHON_BIN="${PYTHON_BIN:-$HOME/.conda/envs/cmip6_equiv_depth_distributions/bin/python}"
MODEL="${MODEL:-NOAA_OLR}"
EXPERIMENT="${EXPERIMENT:-observations}"
START_YEAR="${START_YEAR:-1990}"
END_YEAR="${END_YEAR:-2010}"
N_CALIBRATION="${N_CALIBRATION:-1000}"
N_VALIDATION="${N_VALIDATION:-1000}"
CHUNK_SIZE="${CHUNK_SIZE:-100}"
MAX_CONCURRENT="${MAX_CONCURRENT:-4}"
ACCOUNT="${ACCOUNT:-bb1153}"
PARTITION="${PARTITION:-interactive}"
BACKGROUND_METHOD="${BACKGROUND_METHOD:-gaussian}"
DETECTION_SIGMA="${DETECTION_SIGMA:-1.5}"
case "$DETECTION_SIGMA" in 1|1.0|1.5|2|2.0) ;; *) echo "DETECTION_SIGMA must be 1, 1.5 or 2" >&2; exit 2 ;; esac
JOB_PREFIX="OLR_null"
if [[ "$BACKGROUND_METHOD" == wk99_121 ]]; then JOB_PREFIX="OLR_wk99"; fi
WORKER="$ROOT/scripts/gaussian_event_calibration.py"
[[ $# -le 1 && ( $# -eq 0 || "$1" == --dry-run ) ]] || { echo "Usage: $0 [--dry-run]" >&2; exit 2; }
planning=("$PYTHON_BIN" "$WORKER" plan --model "$MODEL" --experiment "$EXPERIMENT"
    --start-year "$START_YEAR" --end-year "$END_YEAR" --n-calibration "$N_CALIBRATION"
    --n-validation "$N_VALIDATION" --chunk-size "$CHUNK_SIZE" --background-method "$BACKGROUND_METHOD"
    --detection-sigma "$DETECTION_SIGMA")
if [[ -n "${MEMBER:-}" ]]; then planning+=(--member "$MEMBER"); fi
if [[ "${1:-}" == --dry-run ]]; then
    printf 'Plan: '; printf '%q ' "${planning[@]}"; printf '\n'
    echo "Stages: prepare -> one pilot chunk per box profile -> complete simulation array -> threshold reduction -> two injection tasks -> final report"
    echo "Resources/task: 4 CPUs, 6 GiB; simulation concurrency $MAX_CONCURRENT; account $ACCOUNT / partition $PARTITION"
    exit 0
fi
cd "$ROOT"
manifest="$("${planning[@]}")"
jobdir="$(dirname -- "$manifest")/jobs"
mkdir -p "$jobdir"
[[ ! -e "$jobdir/submission.txt" ]] || { echo "Submission already recorded: $jobdir/submission.txt" >&2; exit 1; }
read -r n_tasks pilot_indices < <("$PYTHON_BIN" -c '
import json,sys
m=json.load(open(sys.argv[1])); first={}
for t in m["tasks"]: first.setdefault(t["wave"], t["task_id"])
print(len(m["tasks"]), ",".join(map(str,first.values())))' "$manifest")
common=(--parsable "--account=$ACCOUNT" "--partition=$PARTITION" --nodes=1 --ntasks=1
    --cpus-per-task=4 --mem=6G "--chdir=$ROOT")
base=(env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
    "$PYTHON_BIN" "$WORKER")
printf -v prepare_cmd '%q ' "${base[@]}" prepare --manifest "$manifest"
prepare_id="$(sbatch "${common[@]}" "--job-name=${JOB_PREFIX}_prepare" --time=00:20:00 \
    "--output=$jobdir/prepare.%j.log" --wrap "$prepare_cmd")"
prepare_id="${prepare_id%%;*}"
printf -v simulate_cmd '%q ' "${base[@]}" simulate --manifest "$manifest"
simulate_cmd+=" --task \${SLURM_ARRAY_TASK_ID}"
pilot_id="$(sbatch "${common[@]}" "--job-name=${JOB_PREFIX}_pilot" --time=01:00:00 \
    "--dependency=afterok:$prepare_id" "--array=$pilot_indices%$MAX_CONCURRENT" \
    "--output=$jobdir/pilot.%A_%a.log" --wrap "$simulate_cmd")"
pilot_id="${pilot_id%%;*}"
simulation_id="$(sbatch "${common[@]}" "--job-name=${JOB_PREFIX}_records" --time=01:00:00 \
    "--dependency=afterok:$pilot_id" "--array=0-$((n_tasks-1))%$MAX_CONCURRENT" \
    "--output=$jobdir/simulate.%A_%a.log" --wrap "$simulate_cmd")"
simulation_id="${simulation_id%%;*}"
printf -v summarize_cmd '%q ' "${base[@]}" summarize --manifest "$manifest"
summary_id="$(sbatch "${common[@]}" "--job-name=${JOB_PREFIX}_thresholds" --time=01:00:00 \
    "--dependency=afterok:$simulation_id" "--output=$jobdir/summarize.%j.log" --wrap "$summarize_cmd")"
summary_id="${summary_id%%;*}"
printf -v injection_cmd '%q ' "${base[@]}" inject --manifest "$manifest"
injection_cmd='if [ "$SLURM_ARRAY_TASK_ID" = 0 ]; then wave_null_symmetry=symmetric; else wave_null_symmetry=antisymmetric; fi; exec '"$injection_cmd"' --symmetry "$wave_null_symmetry"'
injection_id="$(sbatch "${common[@]}" "--job-name=${JOB_PREFIX}_injections" --time=01:00:00 \
    "--dependency=afterok:$summary_id" --array=0-1%2 "--output=$jobdir/inject.%A_%a.log" --wrap "$injection_cmd")"
injection_id="${injection_id%%;*}"
printf -v finalize_cmd '%q ' "${base[@]}" finalize --manifest "$manifest"
final_id="$(sbatch "${common[@]}" "--job-name=${JOB_PREFIX}_report" --time=00:10:00 \
    "--dependency=afterok:$injection_id" "--output=$jobdir/finalize.%j.log" --wrap "$finalize_cmd")"
final_id="${final_id%%;*}"
{
    echo "manifest=$manifest"
    echo "prepare=$prepare_id"
    echo "pilot=$pilot_id"
    echo "simulations=$simulation_id"
    echo "thresholds=$summary_id"
    echo "injections=$injection_id"
    echo "final_report=$final_id"
} | tee "$jobdir/submission.txt"
