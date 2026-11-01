#!/bin/bash
# Run the model-wide count/energy audit after the full Slurm array succeeds.
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)"
WORKFLOW_ROOT="${WORKFLOW_ROOT:-$(cd -- "$SCRIPT_DIR/../.." && pwd -P)}"
PYTHON_BIN="${PYTHON_BIN:-$HOME/.conda/envs/cmip6_equiv_depth_distributions/bin/python}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1
cd "$WORKFLOW_ROOT"
exec "$PYTHON_BIN" scripts/codex_cmip_heff/summarize.py
