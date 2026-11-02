#!/usr/bin/env python3
"""Run a disjoint block of held-out NOAA null-fit catalogs in one Slurm job.

From the project root, set --start and --stop within the saved validation
realization indices (1000 through 1999). Each wave/realization worker has its
own log under <ensemble>/jobs/component_fits.<Slurm-job-id>/ and resumes any
incomplete catalog checkpoints.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import subprocess

from src.paths import DATA_DIR, PROJECT_ROOT
from src.wave_events import WAVE_PROFILES


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--stop", type=int, required=True, help="Exclusive realization index")
    parser.add_argument("--max-parallel", type=int, default=16)
    args = parser.parse_args()
    if args.start >= args.stop or args.max_parallel < 1:
        parser.error("Require start < stop and positive max-parallel")
    manifest = args.manifest.absolute()
    manifest.relative_to(DATA_DIR)
    saved = json.loads(manifest.read_text())
    first = saved["configuration"]["n_calibration"]
    last = first + saved["configuration"]["n_validation"]
    if saved["status"] != "complete" or not first <= args.start < args.stop <= last:
        parser.error(f"Require held-out realization indices within {first}:{last}")
    job_id = os.environ.get("SLURM_JOB_ID", "manual")
    log_dir = manifest.parent / "jobs" / f"component_fits.{job_id}"
    log_dir.mkdir(parents=True, exist_ok=True)
    worker = PROJECT_ROOT / "scripts" / "fit_noaa_background_events.py"
    python = Path(os.environ.get("PYTHON_BIN", str(Path.home() / ".conda/envs/cmip6_equiv_depth_distributions/bin/python")))
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1",
               PYTHONDONTWRITEBYTECODE="1", PYTHONUNBUFFERED="1")
    tasks = [(wave, realization) for realization in range(args.start, args.stop)
             for wave in WAVE_PROFILES]

    def run_task(task):
        wave, realization = task
        log = log_dir / f"{wave}.{realization:04d}.log"
        command = [str(python), str(worker), "--manifest", str(manifest),
                   "--wave", wave, "--realization", str(realization)]
        with log.open("w") as output:
            result = subprocess.run(command, cwd=PROJECT_ROOT, env=env,
                                    stdout=output, stderr=subprocess.STDOUT, check=False)
        return task, result.returncode, log

    print(f"{len(tasks)} null catalogs; realizations {args.start}:{args.stop}; "
          f"max_parallel={args.max_parallel}; logs={log_dir}", flush=True)
    failures = []
    with ThreadPoolExecutor(max_workers=args.max_parallel) as pool:
        for future in as_completed(pool.submit(run_task, task) for task in tasks):
            (wave, realization), returncode, log = future.result()
            print(f"{'OK' if returncode == 0 else 'FAILED'} {wave} {realization} "
                  f"exit={returncode} log={log}", flush=True)
            if returncode:
                failures.append((wave, realization))
    if failures:
        raise SystemExit(f"{len(failures)} null catalogs failed: {failures}")
    print(f"Complete: all {len(tasks)} null catalogs", flush=True)


if __name__ == "__main__":
    main()
