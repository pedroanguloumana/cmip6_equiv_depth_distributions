#!/usr/bin/env python3
"""Regrid one manifest task; run inside a Slurm allocation.

    python scripts/regrid_cmip6_rlut.py --manifest data/.../manifest.json --task 0 --cdo cdo
"""

import argparse
import json
from pathlib import Path

from src.paths import CMIP6_NOAA_GRID_DIR, CMIP6_REGRID_JOB_DIR, PROJECT_ROOT
from src.regridding import regrid_rlut


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--task", required=True, type=int, help="Zero-based manifest index")
    parser.add_argument("--cdo", default="cdo", help="CDO executable available on the compute node")
    args = parser.parse_args()
    if PROJECT_ROOT != Path(__file__).resolve().parents[1]:
        raise RuntimeError("Install this checkout with pip -e .; a different src package is active.")
    if not args.manifest.absolute().is_relative_to(CMIP6_REGRID_JOB_DIR):
        parser.error("--manifest must be under data/cmip6_regrid_jobs/.")
    manifest = json.loads(args.manifest.read_text())
    if manifest.get("schema_version") != 1:
        raise ValueError("Unsupported regridding manifest schema.")
    if not 0 <= args.task < len(manifest["tasks"]):
        parser.error("--task is outside the manifest.")
    regrid_rlut(manifest["tasks"][args.task], CMIP6_NOAA_GRID_DIR, cdo=args.cdo)


if __name__ == "__main__":
    main()
