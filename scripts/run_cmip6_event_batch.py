#!/usr/bin/env python3
"""Run one experiment's CMIP6 catalog tasks concurrently inside a Slurm job.

From the project root, use a task list produced by submit_cmip6_rlut_events.sh.
Each worker writes its normal catalog/checkpoints; per-task stdout and stderr
go to data/cmip6_event_jobs/logs/CMIP6_pooled.<job-id>/.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import os
from pathlib import Path
import subprocess

from src.paths import CMIP6_EVENT_JOB_DIR, DATA_DIR, PROJECT_ROOT


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tasklist", type=Path, required=True)
    parser.add_argument("--experiment", choices=("historical", "amip"), required=True)
    parser.add_argument("--max-parallel", type=int, default=16)
    parser.add_argument("--detection-sigma", type=float, choices=(1.0, 1.5, 2.0), default=1.5)
    parser.add_argument("--fit-support", choices=("window", "component"), default="window")
    parser.add_argument("--start-year", type=int, default=1990)
    parser.add_argument("--end-year", type=int, default=2010)
    args = parser.parse_args()
    if args.max_parallel < 1:
        parser.error("--max-parallel must be positive")
    tasklist = args.tasklist.resolve()
    tasklist.relative_to(DATA_DIR)
    with tasklist.open(newline="") as source:
        tasks = [(index, model, experiment, member, wave)
                 for index, (model, experiment, member, wave) in
                 enumerate(csv.reader(source, delimiter="\t"))
                 if experiment == args.experiment]
    if not tasks:
        parser.error(f"No {args.experiment} tasks in {tasklist}")

    job_id = os.environ.get("SLURM_JOB_ID", "manual")
    log_dir = CMIP6_EVENT_JOB_DIR / "logs" / f"CMIP6_pooled.{job_id}"
    log_dir.mkdir(parents=True, exist_ok=True)
    worker = PROJECT_ROOT / "scripts" / "cmip6_rlut_events.py"
    python = Path(os.environ.get("PYTHON_BIN", str(Path.home() / ".conda/envs/cmip6_equiv_depth_distributions/bin/python")))
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1",
               PYTHONDONTWRITEBYTECODE="1", PYTHONUNBUFFERED="1")

    def run_task(task):
        index, model, experiment, member, wave = task
        command = [str(python), str(worker), "--model", model,
                   "--experiment", experiment, "--member", member, "--wave", wave,
                   "--start-year", str(args.start_year), "--end-year", str(args.end_year),
                   "--detection-sigma", str(args.detection_sigma),
                   "--fit-support", args.fit_support]
        log = log_dir / f"task_{index:03d}.{model}.{experiment}.{wave}.log"
        with log.open("w") as output:
            result = subprocess.run(command, cwd=PROJECT_ROOT, env=env,
                                    stdout=output, stderr=subprocess.STDOUT, check=False)
        return task, result.returncode, log

    print(f"{len(tasks)} {args.experiment} tasks; max_parallel={args.max_parallel}; "
          f"detection_sigma={args.detection_sigma}; fit_support={args.fit_support}; "
          f"logs={log_dir}", flush=True)
    failures = []
    with ThreadPoolExecutor(max_workers=args.max_parallel) as pool:
        for future in as_completed(pool.submit(run_task, task) for task in tasks):
            task, returncode, log = future.result()
            print(f"{'OK' if returncode == 0 else 'FAILED'} index={task[0]} "
                  f"model={task[1]} wave={task[4]} exit={returncode} log={log}", flush=True)
            if returncode:
                failures.append(task[0])
    if failures:
        raise SystemExit(f"{len(failures)} tasks failed: {failures}")
    print(f"Complete: all {len(tasks)} {args.experiment} catalogs", flush=True)


if __name__ == "__main__":
    main()
