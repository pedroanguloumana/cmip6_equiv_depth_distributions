#!/bin/bash
# Submit with --dependency=afterok:<segment-array-id> and --chdir=<project-root>.
#SBATCH --account=bb1153
#SBATCH --partition=interactive
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=6G
#SBATCH --time=00:20:00
#SBATCH --job-name=segment_heff_summary
#SBATCH --output=data/wave_significance/codex_segment_heff/jobs/finalize.%j.log
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
WORKFLOW_ROOT="${WORKFLOW_ROOT:-${SLURM_SUBMIT_DIR:-$(cd -- "$SCRIPT_DIR/../.." && pwd -P)}}"
PYTHON_BIN="${PYTHON_BIN:-$HOME/.conda/envs/cmip6_equiv_depth_distributions/bin/python}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
cd "$WORKFLOW_ROOT"
"$PYTHON_BIN" scripts/codex_segment_heff/summarize.py
"$PYTHON_BIN" scripts/codex_segment_heff/execute_notebook.py
