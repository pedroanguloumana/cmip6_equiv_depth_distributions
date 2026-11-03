#!/usr/bin/env python3
"""Rebuild the old distribution inputs from the audited 2-sigma catalogs.

Run on Levante from the project root, with this checkout installed:
    python -m scripts.codex_sigma2_events.remake_distributions --task-index 0
    python -m scripts.codex_sigma2_events.remake_distributions --summarize

Indices 0..59 follow distinct records in the completed catalog_summary.csv.
Workers replay raw OLR through the exact event filter, measure the new active
footprint energies and bin the saved depths. Summarization requires all workers.
Existing threshold-independent observed spectral powers are reused; no null
ensemble, background estimation or catalog fitting is performed.
"""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from src.cmip import load_cmip_rlut
from src.paths import DATA_DIR, PROJECT_ROOT
from src import wave_events, wave_fitting
from src.wave_significance import project_olr
from scripts.noaa_olr_events import load_noaa_record
from scripts.codex_cmip_heff.run import event_energy_rows, preprocessed_anomaly
from scripts.codex_sigma2_events.summarize import read_catalog, catalog_statistics


ROOT = DATA_DIR / "event_threshold_comparison"
OUTPUT = ROOT / "distributions_sigma2"
PREVIOUS = DATA_DIR / "wave_significance/codex_cmip_heff_comparison"
WAVES = ("Kelvin", "MRG", "n0EIG", "n1ERW", "n1WIG")
EDGES_M = np.arange(1.0, 201.0)
LOG_EDGES_M = np.geomspace(0.5, 10_000.0, 51)
KEY = ["model", "experiment", "member"]


def cohort():
    """Return the completed source/wave audit and its 60 distinct records."""
    report = json.loads((ROOT / "audit.json").read_text())
    summary = pd.read_csv(ROOT / "catalog_summary.csv")
    selected = summary.loc[summary.detection_sigma.eq(2)].copy()
    records = selected[KEY].drop_duplicates().reset_index(drop=True)
    if report["status"] != "complete" or len(selected) != 300 or len(records) != 60:
        raise ValueError("Expected the completed 300-catalog 2-sigma cohort.")
    return selected, records


def write_table(frame, path):
    """Write a new compact table atomically, preserving previous results."""
    if path.exists():
        raise FileExistsError(path)
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(temporary, index=False, compression="gzip")
    temporary.replace(path)


def count_rows(frame, identity, wave, record_days):
    """Count saved depths on the old notebooks' linear/logarithmic bin grids."""
    valid = (frame.fit_success & frame.dispersion_valid
             & np.isfinite(frame.h_eff_m) & frame.h_eff_m.gt(0))
    rows = []
    for selection, keep in (("screened", frame.keep_for_histogram), ("relaxed", valid)):
        depths = frame.loc[keep, "h_eff_m"].to_numpy()
        row = dict(**identity, wave_type=wave, selection=selection, detection_sigma=2.0,
                   record_days=record_days, candidates=len(frame), selected_events=len(depths),
                   below_1m=int((depths < 1).sum()), above_200m=int((depths > 200).sum()),
                   below_0p5m=int((depths < .5).sum()), above_10000m=int((depths > 10_000).sum()))
        linear = np.histogram(depths, bins=EDGES_M)[0]
        logarithmic = np.histogram(depths, bins=LOG_EDGES_M)[0]
        if (linear.sum() + row["below_1m"] + row["above_200m"] != len(depths)
                or logarithmic.sum() + row["below_0p5m"] + row["above_10000m"] != len(depths)):
            raise ValueError("A displayed depth histogram loses selected events.")
        row.update({f"count_{int(left):03d}m": int(n) for left, n in zip(EDGES_M[:-1], linear)})
        row.update({f"log_count_{i:02d}": int(n) for i, n in enumerate(logarithmic)})
        rows.append(row)
    return rows


def run_record(task_index):
    """Replay one NOAA/model record and its five 2-sigma component catalogs."""
    summary, records = cohort()
    if not 0 <= task_index < len(records):
        raise ValueError("task-index must be in 0..59")
    identity = records.iloc[task_index].to_dict()
    model, experiment, member = (identity[key] for key in KEY)
    stem = f"{model}.{experiment}.{member}"
    OUTPUT.mkdir(parents=True, exist_ok=True)
    counts_path = OUTPUT / f"{stem}.counts.csv.gz"
    energy_path = OUTPUT / f"{stem}.energy.csv.gz"
    receipt_path = OUTPUT / f"{stem}.json"
    if receipt_path.exists() and counts_path.exists() and energy_path.exists():
        print(f"Already complete: {stem}", flush=True)
        return
    if any(path.exists() for path in (counts_path, energy_path, receipt_path)):
        raise RuntimeError(f"Partial output for {stem}; inspect before retrying.")

    catalogs, provenance = {}, []
    for wave in WAVES:
        frame, metadata, path = read_catalog(model, experiment, member, wave, 2)
        cfg = metadata["configuration"]
        for name, module in (("wave_events", wave_events), ("wave_fitting", wave_fitting)):
            if hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest() != cfg["code_sha256"][name]:
                raise ValueError(f"The catalog's scientific source changed: {name}")
        for item in cfg["input"]:
            stat = (PROJECT_ROOT / item["file"]).stat()
            if (stat.st_size, stat.st_mtime_ns) != (item["size_bytes"], item["mtime_ns"]):
                raise ValueError(f"The catalog input changed: {item['file']}")
        statistics, _ = catalog_statistics(frame, metadata)
        row = summary.loc[summary.model.eq(model) & summary.experiment.eq(experiment)
                          & summary.member.eq(member) & summary.wave.eq(wave)].squeeze()
        if (statistics["candidates"], statistics["valid_depths"], statistics["screened_events"]) != (
                row.candidates, row.valid_depths, row.screened_events):
            raise ValueError(f"Catalog counts changed since the completed audit: {stem}/{wave}")
        catalogs[wave] = frame
        provenance.append(dict(wave=wave, catalog=str(path.relative_to(DATA_DIR)),
                               sha256=hashlib.sha256(path.read_bytes()).hexdigest()))

    olr = (load_noaa_record("1990-01-01", "2010-12-31") if experiment == "observations"
           else load_cmip_rlut(model, experiment=experiment, member_id=member,
                              year_range=(1990, 2010), lat_range=wave_events.LAT_RANGE))
    components = project_olr(olr)
    denominator = float(sum(np.sum(preprocessed_anomaly(components[symmetry]).values ** 2)
                            for symmetry in ("symmetric", "antisymmetric")))
    del components
    if not np.isfinite(denominator) or denominator <= 0:
        raise ValueError(f"Nonpositive projected OLR energy: {stem}")
    counts, energies = [], []
    for wave in WAVES:
        field = wave_events.filter_olr_for_wave(olr, wave)
        frame = catalogs[wave]
        energies.extend(event_energy_rows(field, wave, frame.set_index("event_id"),
                                          denominator, identity, detection_sigma=2))
        counts.extend(count_rows(frame, identity, wave, olr.sizes["time"]))
        print(f"{stem}/{wave}: matched {len(frame)} 2-sigma events", flush=True)
        del field
    energy = pd.DataFrame(energies).assign(detection_sigma=2.0)
    write_table(pd.DataFrame(counts), counts_path)
    write_table(energy, energy_path)
    receipt_path.write_text(json.dumps(dict(
        status="complete", task_index=task_index, **identity, detection_sigma=2.0,
        catalogs=provenance, record_days=olr.sizes["time"],
        calculation="Raw-input replay of the exact catalog filter/detector; no catalog refitting or Gaussian calculations",
        worker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    ), indent=2) + "\n")
    print(f"Complete: {stem}", flush=True)


def summarize():
    """Audit all 60 replays and build the six tables read by the old notebooks."""
    report_path = OUTPUT / "audit.json"
    if report_path.exists():
        raise FileExistsError(report_path)
    _, records = cohort()
    count_frames, energy_frames = [], []
    for identity in records.to_dict("records"):
        stem = ".".join(identity[key] for key in KEY)
        receipt = json.loads((OUTPUT / f"{stem}.json").read_text())
        if receipt["status"] != "complete" or receipt["detection_sigma"] != 2:
            raise ValueError(f"Incomplete 2-sigma replay: {stem}")
        for kind, destination in (("counts", count_frames), ("energy", energy_frames)):
            frame = pd.read_csv(OUTPUT / f"{stem}.{kind}.csv.gz")
            if (len(frame) != 10 or set(frame.wave_type) != set(WAVES)
                    or frame.duplicated(["wave_type", "selection"]).any()
                    or not frame.detection_sigma.eq(2).all()
                    or not frame[KEY].eq(pd.Series(identity)).all().all()):
                raise ValueError(f"Misidentified or incomplete {kind}: {stem}")
            destination.append(frame)
    counts, energy = pd.concat(count_frames, ignore_index=True), pd.concat(energy_frames, ignore_index=True)
    selection_key = KEY + ["wave_type", "selection"]
    paired = counts.merge(energy, on=selection_key, validate="one_to_one", suffixes=("_count", "_energy"))
    if (not paired.candidates_count.eq(paired.candidates_energy).all()
            or not paired.selected_events_count.eq(paired.selected_events_energy).all()):
        raise ValueError("Replayed energies and saved depth counts disagree.")
    numeric = energy.select_dtypes("number")
    if (not np.isfinite(numeric).all().all() or (numeric < 0).any().any()
            or (energy.selected_fitted_energy > energy.selected_support_energy * (1+1e-8)).any()
            or (energy.selected_support_energy > energy.all_candidate_support_energy * (1+1e-8)).any()
            or (energy.all_candidate_support_energy > energy.filtered_energy * (1+1e-8)).any()):
        raise ValueError("An event energy invariant failed.")
    for metric in ("support", "fitted"):
        columns = [f"{metric}_{i:03d}m" for i in range(1, 200)]
        if (energy[columns].sum(axis=1) > energy[f"selected_{metric}_energy"] * (1+1e-8)).any():
            raise ValueError(f"Binned {metric} energy exceeds the full selected population.")
    for label, observed in (("all_model", False), ("noaa", True)):
        is_observed = counts.experiment.eq("observations")
        write_table(counts.loc[is_observed.eq(observed)], OUTPUT / f"{label}_event_counts.csv.gz")
        is_observed = energy.experiment.eq("observations")
        write_table(energy.loc[is_observed.eq(observed)], OUTPUT / f"{label}_event_energy.csv.gz")
        # Actual-data spectral power does not depend on the envelope threshold.
        spectral = pd.read_csv(PREVIOUS / f"{label}_spectral_power.csv.gz")
        spectral = spectral.loc[spectral.source.eq("observed")].copy()
        if len(spectral) != (5 if observed else 295) or spectral.duplicated(KEY + ["wave_type"]).any():
            raise ValueError(f"Incomplete threshold-independent observed spectra: {label}")
        write_table(spectral, OUTPUT / f"{label}_spectral_power.csv.gz")
    report_path.write_text(json.dumps(dict(
        status="complete", detection_sigma=2.0, source_records=60, catalogs=300,
        count_rows=len(counts), energy_rows=len(energy),
        checks="Source fingerprints and inputs; full raw-data detection replay IDs/peaks; saved fit/depth flags; selected counts; histogram closure; energy bounds",
        spectral_power="Reused existing observed-only power; independent of envelope threshold",
        background="No Gaussian-background calculation, null selection or comparison",
    ), indent=2) + "\n")
    print(f"Audited all 300 2-sigma catalogs and 60 raw-input replays: {OUTPUT}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--task-index", type=int)
    mode.add_argument("--summarize", action="store_true")
    args = parser.parse_args()
    summarize() if args.summarize else run_record(args.task_index)
