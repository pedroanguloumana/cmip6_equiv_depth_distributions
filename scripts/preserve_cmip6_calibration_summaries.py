#!/usr/bin/env python3
"""Preserve small, complete CMIP6 calibration products from scratch on home.

Run after the CMIP6 Gaussian calibration array. Full null chunks and spectra
remain at the scratch-backed data/ link; the copied tables are sufficient to
review thresholds, significance flags, injection summaries and provenance.
"""

import argparse
import json
import os
from pathlib import Path
import shutil

from src.paths import WAVE_SIGNIFICANCE_DIR


SOURCE = WAVE_SIGNIFICANCE_DIR / "gaussian_calibration_wk99_121_cmip6"
DESTINATION = WAVE_SIGNIFICANCE_DIR / "cmip6_calibration_summaries"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expected", type=int, default=59)
    args = parser.parse_args()
    if not SOURCE.is_symlink() or not SOURCE.resolve().is_dir():
        raise RuntimeError(f"Configure the CMIP calibration data symlink: {SOURCE}")
    DESTINATION.mkdir(parents=True, exist_ok=True)
    completed = []
    for run in sorted(SOURCE.iterdir()):
        if not run.is_dir() or not (run / "report.json").is_file():
            continue
        report = json.loads((run / "report.json").read_text())
        manifest = json.loads((run / "manifest.json").read_text())
        if report.get("status") != "complete" or manifest.get("status") != "complete":
            continue
        target = DESTINATION / run.name
        if target.exists():
            existing = json.loads((target / "report.json").read_text())
            if existing.get("ensemble_id") != report.get("ensemble_id"):
                raise RuntimeError(f"Conflicting preserved ensemble: {target}")
            completed.append(run.name)
            continue
        files = [
            "manifest.json", "report.json", "calibration_report.json",
            "thresholds.csv", "threshold_convergence.csv", "null_records.csv",
            "observed_distributions.csv", "background_diagnostics.json",
            "injection_recovery.csv", "injections.csv.gz",
        ]
        for wave in ("Kelvin", "n1ERW", "MRG", "n0EIG", "n1WIG"):
            files.extend((f"observed_significance.{wave}.csv",
                          f"threshold_curve.{wave}.csv",
                          f"record_max_ecdf.{wave}.csv",
                          f"spectral_summary.{wave}.csv"))
        missing = [name for name in files if not (run / name).is_file()]
        if missing:
            raise RuntimeError(f"{run.name} missing final products: {missing}")
        temporary = DESTINATION / f".{run.name}.tmp.{os.getpid()}"
        temporary.mkdir()
        for name in files:
            shutil.copy2(run / name, temporary / name)
        for folder in ("code_snapshot", "analysis_code_snapshot"):
            shutil.copytree(run / folder, temporary / folder)
        temporary.rename(target)
        completed.append(run.name)
        print(f"Preserved {run.name}", flush=True)
    index = dict(source=str(SOURCE), destination=str(DESTINATION),
                 complete_runs=len(completed), runs=completed,
                 scope="Reduced calibration tables, flags, injections and provenance; full null chunks remain on scratch")
    temporary_index = DESTINATION / f".index.tmp.{os.getpid()}"
    temporary_index.write_text(json.dumps(index, indent=2) + "\n")
    temporary_index.replace(DESTINATION / "index.json")
    print(f"Preserved {len(completed)} complete calibration summaries", flush=True)
    if len(completed) != args.expected:
        raise SystemExit(f"Expected {args.expected} complete calibrations; got {len(completed)}")


if __name__ == "__main__":
    main()
