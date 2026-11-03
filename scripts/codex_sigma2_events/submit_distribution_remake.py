#!/usr/bin/env python3
"""Run the 60 record replays after successful NOAA/Gregorian/360-day pilots.

From the project root on Levante:
    python -m scripts.codex_sigma2_events.submit_distribution_remake
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
import shlex
import subprocess
import sys

from src.paths import DATA_DIR, PROJECT_ROOT
from .remake_distributions import OUTPUT, cohort


JOBS = DATA_DIR / "cmip6_event_jobs/sigma2_component/distribution_remake"
WORKERS = 12


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--shard-index", type=int, choices=(0, 1))
    args = parser.parse_args()
    _, records = cohort()
    JOBS.mkdir(parents=True, exist_ok=True)
    if args.shard_index is not None:
        def run_task(index):
            path = JOBS / f"task_{index:02d}.log"
            command = [sys.executable, "-u", "-m", "scripts.codex_sigma2_events.remake_distributions",
                       "--task-index", str(index)]
            env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
                       MKL_NUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
            with path.open("w") as output:
                result = subprocess.run(command, cwd=PROJECT_ROOT, env=env,
                                        stdout=output, stderr=subprocess.STDOUT)
            return index, result.returncode

        first = 30 * args.shard_index
        failures = []
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            for future in as_completed(pool.submit(run_task, i) for i in range(first, first + 30)):
                index, code = future.result()
                print(f"{'OK' if code == 0 else 'FAILED'} task {index}; exit {code}", flush=True)
                if code:
                    failures.append(index)
        if failures:
            raise SystemExit(f"Replay tasks failed: {failures}")
        return

    submission = JOBS / "submission.json"
    if submission.exists():
        raise FileExistsError(f"Preserving previous submissions: {submission}")
    worker = PROJECT_ROOT / "scripts/codex_sigma2_events/remake_distributions.py"
    worker_hash = hashlib.sha256(worker.read_bytes()).hexdigest()
    for index in (0, 2, 32):
        stem = ".".join(records.iloc[index][key] for key in ("model", "experiment", "member"))
        receipt = json.loads((OUTPUT / f"{stem}.json").read_text())
        if receipt["status"] != "complete" or receipt["worker_sha256"] != worker_hash:
            raise ValueError(f"Missing a successful current pilot: {stem}")
    base = ["sbatch", "--parsable", "--account=bb1153", "--partition=interactive",
            "--nodes=1", "--ntasks=1", f"--chdir={PROJECT_ROOT}"]
    wrap = shlex.join([sys.executable, "-u", "-m", "scripts.codex_sigma2_events.submit_distribution_remake",
                       "--shard-index"]) + ' "$SLURM_ARRAY_TASK_ID"'
    command = base + ["--cpus-per-task=48", "--mem=170G", "--time=01:00:00", "--array=0-1%2",
                      "--job-name=ds2_distribution_replay", f"--output={JOBS}/pooled.%A_%a.log", f"--wrap={wrap}"]
    job_id = subprocess.check_output(command, text=True).strip().split(";")[0]
    state = dict(status="submitting", model_array=job_id, workers_per_shard=WORKERS,
                 worker_sha256=worker_hash, pilot_array="28679407")
    submission.write_text(json.dumps(state, indent=2) + "\n")
    wrap = shlex.join([sys.executable, "-u", "-m", "scripts.codex_sigma2_events.remake_distributions", "--summarize"])
    command = base + ["--cpus-per-task=2", "--mem=6G", "--time=00:20:00",
                      "--job-name=ds2_distribution_audit", f"--dependency=afterok:{job_id}",
                      f"--output={JOBS}/audit.%j.log", f"--wrap={wrap}"]
    state["audit_job"] = subprocess.check_output(command, text=True).strip().split(";")[0]
    state["status"] = "submitted"
    submission.write_text(json.dumps(state, indent=2) + "\n")
    print(json.dumps(state, indent=2), flush=True)


if __name__ == "__main__":
    main()
