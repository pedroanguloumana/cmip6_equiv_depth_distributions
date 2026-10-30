#!/bin/bash
# One array task per manifest entry (model/experiment/member). From project root:
#   CDO_BIN=cdo bash scripts/submit_cmip6_rlut_regridding.sh data/.../manifest.json --dry-run
#   CDO_BIN=cdo bash scripts/submit_cmip6_rlut_regridding.sh data/.../manifest.json
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
WORKFLOW_ROOT="${WORKFLOW_ROOT:-$(cd -- "$SCRIPT_DIR/.." && pwd -P)}"
PYTHON_BIN="${PYTHON_BIN:-$HOME/.conda/envs/cmip6_equiv_depth_distributions/bin/python}"
CDO_BIN="${CDO_BIN:-cdo}"
ACCOUNT="${ACCOUNT:-bb1153}"
PARTITION="${PARTITION:-interactive}"
JOB_MEM="${JOB_MEM:-8G}"
JOB_CPUS="${JOB_CPUS:-6}"
JOB_TIME="${JOB_TIME:-01:00:00}"
MAX_CONCURRENT="${MAX_CONCURRENT:-4}"

cd "$WORKFLOW_ROOT"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    : "${MANIFEST:?Launcher must provide MANIFEST}"
    echo "host=$(hostname) job=$SLURM_JOB_ID task=$SLURM_ARRAY_TASK_ID"
    exec "$PYTHON_BIN" "$WORKFLOW_ROOT/scripts/regrid_cmip6_rlut.py" \
        --manifest "$MANIFEST" --task "$SLURM_ARRAY_TASK_ID" --cdo "$CDO_BIN"
fi

if [[ $# -lt 1 || $# -gt 2 || ( $# -eq 2 && "$2" != --dry-run ) ]]; then
    echo "Usage: bash scripts/submit_cmip6_rlut_regridding.sh MANIFEST [--dry-run]" >&2
    exit 2
fi
MANIFEST="$1"
command -v "$CDO_BIN" >/dev/null || { echo "CDO not executable: $CDO_BIN" >&2; exit 1; }
[[ "$MAX_CONCURRENT" =~ ^[1-9][0-9]*$ ]] || { echo "MAX_CONCURRENT must be positive" >&2; exit 2; }
mapfile -t settings < <("$PYTHON_BIN" - "$MANIFEST" <<'PY'
import json
from pathlib import Path
import sys
from src.paths import CMIP6_NOAA_GRID_DIR, CMIP6_REGRID_JOB_DIR
path = Path(sys.argv[1]).absolute()
if not path.is_relative_to(CMIP6_REGRID_JOB_DIR):
    raise ValueError("Manifest must be under data/cmip6_regrid_jobs/.")
manifest = json.loads(path.read_text())
if manifest.get("schema_version") != 1 or not manifest["tasks"]:
    raise ValueError("Expected a nonempty version-1 regridding manifest.")
if not CMIP6_NOAA_GRID_DIR.is_symlink() or not CMIP6_NOAA_GRID_DIR.is_dir():
    raise FileNotFoundError("Configure data/cmip6_rlut_noaa1deg as an output archive symlink.")
print(path)
print(len(manifest["tasks"]))
print(CMIP6_REGRID_JOB_DIR / "logs")
PY
)
[[ ${#settings[@]} -eq 3 ]] || { echo "Manifest/output validation failed" >&2; exit 1; }
MANIFEST="${settings[0]}"
TASK_COUNT="${settings[1]}"
LOGDIR="${settings[2]}"
ARRAY_SPEC="${ARRAY_SPEC:-0-$((TASK_COUNT-1))%$MAX_CONCURRENT}"
command=(sbatch --parsable --job-name=CMIP6_rlut_regrid
    "--account=$ACCOUNT" "--partition=$PARTITION" --nodes=1 --ntasks=1
    "--cpus-per-task=$JOB_CPUS" "--mem=$JOB_MEM" "--time=$JOB_TIME"
    "--array=$ARRAY_SPEC" "--chdir=$WORKFLOW_ROOT"
    "--output=$LOGDIR/regrid.%A_%a.log" "--error=$LOGDIR/regrid.%A_%a.log"
    "--export=ALL,WORKFLOW_ROOT=$WORKFLOW_ROOT,PYTHON_BIN=$PYTHON_BIN,CDO_BIN=$CDO_BIN,MANIFEST=$MANIFEST"
    "$SCRIPT_DIR/$(basename -- "${BASH_SOURCE[0]}")")
echo "Manifest: $MANIFEST; tasks=$TASK_COUNT; array=$ARRAY_SPEC"
echo "Resources/task: $PARTITION $JOB_MEM $JOB_CPUS CPUs $JOB_TIME"
if [[ "${2:-}" == --dry-run ]]; then
    printf '%q ' "${command[@]}"; printf '\n'
    exit 0
fi
mkdir -p "$LOGDIR"
job_id="$("${command[@]}")"
echo "Submitted array $job_id"
echo "Logs: $LOGDIR/regrid.${job_id%%;*}_<task>.log"
