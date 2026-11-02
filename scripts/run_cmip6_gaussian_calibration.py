#!/usr/bin/env python3
"""Run one CMIP6 WK99 Gaussian calibration within a Slurm array task.

From the project root on Levante, pass a data/ task list containing
model, experiment and member TSV rows. SLURM_ARRAY_TASK_ID selects one row.
Each existing calibration stage runs as a separate process. Simulation chunks
run concurrently and retain their normal checkpoints.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import json
import os
from pathlib import Path
import subprocess

from src.paths import DATA_DIR, PROJECT_ROOT


BACKGROUND_METHOD = "wk99_121"
DETECTION_SIGMA = 1.5
N_CALIBRATION = 1000
N_VALIDATION = 1000
CHUNK_SIZE = 100
SEED = 20261031


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tasklist", type=Path, required=True)
    parser.add_argument("--task-index", type=int,
                        default=int(os.environ.get("SLURM_ARRAY_TASK_ID", "-1")))
    parser.add_argument("--max-parallel", type=int, default=32)
    args = parser.parse_args()
    if args.max_parallel < 1 or args.task_index < 0:
        parser.error("Require a nonnegative task index and positive max-parallel")
    tasklist = args.tasklist.resolve()
    tasklist.relative_to(DATA_DIR)
    with tasklist.open(newline="") as source:
        simulations = list(csv.reader(source, delimiter="\t"))
    if args.task_index >= len(simulations):
        parser.error(f"Task index {args.task_index} exceeds {len(simulations)} simulations")
    model, experiment, member = simulations[args.task_index]
    if experiment not in ("historical", "amip"):
        parser.error(f"Unsupported experiment: {experiment}")

    python = Path(os.environ.get("PYTHON_BIN", str(Path.home() / ".conda/envs/cmip6_equiv_depth_distributions/bin/python")))
    worker = PROJECT_ROOT / "scripts" / "gaussian_event_calibration.py"
    base = [str(python), str(worker)]
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1",
               PYTHONDONTWRITEBYTECODE="1", PYTHONUNBUFFERED="1")
    print(f"task={args.task_index} model={model} experiment={experiment} member={member} "
          f"background={BACKGROUND_METHOD} detection_sigma={DETECTION_SIGMA} "
          f"calibration={N_CALIBRATION} validation={N_VALIDATION} "
          f"chunk_size={CHUNK_SIZE} max_parallel={args.max_parallel}", flush=True)

    plan = subprocess.run(base + ["plan", "--model", model, "--experiment", experiment,
        "--member", member, "--start-year", "1990", "--end-year", "2010",
        "--n-calibration", str(N_CALIBRATION), "--n-validation", str(N_VALIDATION),
        "--chunk-size", str(CHUNK_SIZE), "--seed", str(SEED),
        "--background-method", BACKGROUND_METHOD,
        "--detection-sigma", str(DETECTION_SIGMA)], cwd=PROJECT_ROOT, env=env,
        capture_output=True, text=True, check=True)
    manifest_path = Path(plan.stdout.strip().splitlines()[-1])
    print(f"manifest={manifest_path}", flush=True)
    root = manifest_path.parent
    log_dir = root / "jobs" / "pooled"
    log_dir.mkdir(parents=True, exist_ok=True)
    if (root / "report.json").is_file():
        report = json.loads((root / "report.json").read_text())
        if report.get("status") == "complete":
            print("Already complete", flush=True)
            return

    def run_stage(name, *options, log_name=None):
        log = log_dir / (log_name or f"{name}.log")
        with log.open("w") as output:
            result = subprocess.run(base + [name, "--manifest", str(manifest_path), *options],
                                    cwd=PROJECT_ROOT, env=env, stdout=output,
                                    stderr=subprocess.STDOUT, check=False)
        return result.returncode, log

    def require_stage(name, *options):
        returncode, log = run_stage(name, *options)
        print(f"{name}: exit={returncode} log={log}", flush=True)
        if returncode:
            raise SystemExit(f"{name} failed; inspect {log}")

    require_stage("prepare")
    manifest = json.loads(manifest_path.read_text())
    tasks = manifest["tasks"]
    pending = []
    for task in tasks:
        receipt = root / "chunks" / f"task_{task['task_id']:05d}.json"
        if receipt.is_file():
            status = json.loads(receipt.read_text())
            if status.get("status") == "complete" and status.get("ensemble_id") == manifest["ensemble_id"]:
                continue
        pending.append(task["task_id"])
    print(f"simulation chunks: {len(tasks) - len(pending)} complete; {len(pending)} pending", flush=True)
    failures = []
    with ThreadPoolExecutor(max_workers=args.max_parallel) as pool:
        futures = {pool.submit(run_stage, "simulate", "--task", str(index),
                               log_name=f"simulate_{index:05d}.log"): index
                   for index in pending}
        for future in as_completed(futures):
            index = futures[future]
            returncode, log = future.result()
            print(f"simulate {index}: exit={returncode} log={log}", flush=True)
            if returncode:
                failures.append(index)
    if failures:
        raise SystemExit(f"Simulation chunks failed: {failures}")

    manifest = json.loads(manifest_path.read_text())
    if manifest["status"] not in ("calibrated", "complete"):
        require_stage("summarize")
    injection_jobs = []
    for symmetry in ("symmetric", "antisymmetric"):
        receipt = root / f"injections.{symmetry}.json"
        if receipt.is_file():
            status = json.loads(receipt.read_text())
            if status.get("status") == "complete" and status.get("ensemble_id") == manifest["ensemble_id"]:
                continue
        injection_jobs.append(symmetry)
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = {pool.submit(run_stage, "inject", "--symmetry", symmetry,
                               log_name=f"inject_{symmetry}.log"): symmetry
                   for symmetry in injection_jobs}
        for future in as_completed(futures):
            symmetry = futures[future]
            returncode, log = future.result()
            print(f"inject {symmetry}: exit={returncode} log={log}", flush=True)
            if returncode:
                failures.append(symmetry)
    if failures:
        raise SystemExit(f"Injection stages failed: {failures}")
    require_stage("finalize")
    print(f"Complete: {model} {experiment} {member} {manifest_path}", flush=True)


if __name__ == "__main__":
    main()
