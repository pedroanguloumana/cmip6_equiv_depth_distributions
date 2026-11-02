#!/usr/bin/env python3
"""Run the existing 2-sigma component fitter for NOAA and the pinned CMIP cohort.

From the project root on Levante:
    python scripts/codex_sigma2_events/submit.py

Requires the 15 successful three-event NOAA/CanESM5 pilots. It submits NOAA,
four CMIP shards and an audit, preserving all earlier catalogs and checkpoints.
"""

import argparse
import csv
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys

from src.paths import CMIP6_EVENT_DIR, DATA_DIR, PROJECT_ROOT


TASKLIST = DATA_DIR / "cmip6_event_jobs/tasklists/simulations_waves_20261030_205901_2794174.tsv"
JOBS = DATA_DIR / "cmip6_event_jobs/sigma2_component"
OUTPUT = DATA_DIR / "event_threshold_comparison"
WAVES = ("Kelvin", "n1ERW", "MRG", "n0EIG", "n1WIG")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-index", type=int, choices=range(4))
    args = parser.parse_args()
    if args.task_index is not None:
        shard = json.loads((JOBS / "shards.json").read_text())[args.task_index]
        command = [sys.executable, str(PROJECT_ROOT / "scripts/run_cmip6_event_batch.py"),
                   "--tasklist", str(JOBS / shard["tasklist"]),
                   "--experiment", shard["experiment"], "--max-parallel", "12",
                   "--start-year", "1990", "--end-year", "2010",
                   "--detection-sigma", "2", "--fit-support", "component"]
        os.execv(sys.executable, command)

    submission = JOBS / "submission.json"
    if submission.exists():
        raise FileExistsError(f"Preserving recorded submissions: {submission}")
    with TASKLIST.open(newline="") as handle:
        tasks = list(csv.reader(handle, delimiter="\t"))
    if len(tasks) != 295 or len(set(map(tuple, tasks))) != 295:
        raise ValueError("Expected the pinned 295 model/experiment/member/wave tasks.")
    simulations = sorted(set(tuple(row[:3]) for row in tasks))
    if (len(simulations) != 59 or sum(s[1] == "historical" for s in simulations) != 34
            or sum(s[1] == "amip" for s in simulations) != 25):
        raise ValueError("Expected 34 Historical and 25 AMIP simulations.")
    for simulation in simulations:
        if {row[3] for row in tasks if tuple(row[:3]) == simulation} != set(WAVES):
            raise ValueError(f"Incomplete wave cohort: {simulation}")

    pilots = [path for directory in (DATA_DIR, CMIP6_EVENT_DIR)
              for path in directory.glob("*.ds2.component.*.pilot3.json")]
    required_pilots = {(model, experiment, wave)
                       for model, experiment in (("NOAA_CDR", "observations"),
                                                 ("CanESM5", "historical"), ("CanESM5", "amip"))
                       for wave in WAVES}
    complete = set()
    for path in pilots:
        metadata = json.loads(path.read_text())
        cfg = metadata["configuration"]
        identity = cfg["identity"]
        if (metadata["status"] == "complete" and metadata["catalog_candidates"] == 3
                and metadata["counts"]["fit_success"] == 3
                and cfg["detection"]["detection_sigma"] == 2
                and cfg["fit_support"] == "component"):
            complete.add((identity["source_id"], identity["experiment_id"], cfg["wave_type"]))
    if not required_pilots.issubset(complete):
        raise ValueError(f"Missing successful three-event pilots: {required_pilots - complete}")

    JOBS.mkdir(parents=True, exist_ok=True)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    shards = []
    for experiment in ("historical", "amip"):
        selected = [s for s in simulations if s[1] == experiment]
        midpoint = len(selected) // 2
        for number, group in enumerate((selected[:midpoint], selected[midpoint:])):
            rows = [row for row in tasks if tuple(row[:3]) in set(group)]
            path = JOBS / f"{experiment}.{number}.tsv"
            with path.open("w", newline="") as handle:
                csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)
            shards.append(dict(experiment=experiment, tasklist=path.name, tasks=len(rows)))
    (JOBS / "shards.json").write_text(json.dumps(shards, indent=2) + "\n")
    state = dict(status="submitting", detection_sigma=2, fit_support="component",
                 period="1990--2010, native calendars", source_tasklist=str(TASKLIST.relative_to(DATA_DIR)),
                 shards=shards, jobs={})
    submission.write_text(json.dumps(state, indent=2) + "\n")

    env = dict(os.environ, DETECTION_SIGMA="2", FIT_SUPPORT="component", MAX_EVENTS="",
               MAX_CONCURRENT="1", JOB_TIME="02:00:00", PYTHON_BIN=sys.executable)
    noaa = subprocess.run(["bash", "scripts/submit_noaa_olr_events.sh"], cwd=PROJECT_ROOT,
                          env=env, check=True, capture_output=True, text=True)
    print(noaa.stdout, flush=True)
    noaa_id = re.search(r"Submitted array (\d+)", noaa.stdout).group(1)
    state["jobs"]["noaa"] = noaa_id
    submission.write_text(json.dumps(state, indent=2) + "\n")
    wrap = shlex.join([sys.executable, "scripts/codex_sigma2_events/submit.py", "--task-index"])
    wrap += ' "$SLURM_ARRAY_TASK_ID"'
    command = ["sbatch", "--parsable", "--account=bb1153", "--partition=interactive",
               "--nodes=1", "--ntasks=1", "--cpus-per-task=96", "--mem=170G", "--time=04:00:00",
               "--array=0-3%4", "--job-name=CMIP6_sigma2_component", f"--chdir={PROJECT_ROOT}",
               f"--output={JOBS}/pooled.%A_%a.log", f"--export=ALL,PYTHON_BIN={sys.executable}",
               f"--wrap={wrap}"]
    model_id = subprocess.check_output(command, cwd=PROJECT_ROOT, text=True).strip().split(";")[0]
    state["jobs"]["cmip6"] = model_id
    submission.write_text(json.dumps(state, indent=2) + "\n")
    audit_wrap = shlex.join([sys.executable, "scripts/codex_sigma2_events/summarize.py"])
    command = ["sbatch", "--parsable", "--account=bb1153", "--partition=interactive",
               "--nodes=1", "--ntasks=1", "--cpus-per-task=2", "--mem=6G", "--time=00:30:00",
               "--job-name=sigma2_catalog_audit", f"--chdir={PROJECT_ROOT}",
               f"--dependency=afterok:{noaa_id}:{model_id}", f"--output={JOBS}/audit.%j.log",
               f"--wrap={audit_wrap}"]
    audit_id = subprocess.check_output(command, cwd=PROJECT_ROOT, text=True).strip().split(";")[0]
    state["jobs"]["audit"] = audit_id
    state["status"] = "submitted"
    submission.write_text(json.dumps(state, indent=2) + "\n")
    print(json.dumps(state, indent=2), flush=True)


if __name__ == "__main__":
    main()
