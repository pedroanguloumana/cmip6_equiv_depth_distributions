#!/usr/bin/env python3
"""Audit and combine all CMIP6 event-depth and spectral-power diagnostics.

Run on Levante from the project root after the 59 model tasks complete:
    python scripts/codex_cmip_heff/summarize.py
"""

import json

import numpy as np
import pandas as pd

from src.paths import DATA_DIR
from run import EDGES_M, OUTPUT, WAVES, event_catalog, simulation, spectral_energy_rows


LOG_EDGES_M = np.logspace(np.log10(0.5), np.log10(10_000), 51)
SELECTIONS = ("screened", "relaxed")
ENERGY_COLUMNS = {metric: [f"{metric}_{left:03d}m" for left in EDGES_M[:-1].astype(int)]
                  for metric in ("support", "fitted")}
POWER_COLUMNS = [f"power_{left:03d}m" for left in EDGES_M[:-1].astype(int)]
NOAA_RUN = (DATA_DIR / "wave_significance/gaussian_calibration_wk99_121"
            / "NOAA_OLR.observations.not_applicable.1990-2010.b898b1ae4307")


def write_csv(frame, path):
    if path.exists():
        raise FileExistsError(f"Preserving existing audit output: {path}")
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(temporary, index=False, compression="gzip")
    temporary.replace(path)


def count_rows(catalog, identity, wave, record_days):
    """Count screened and relaxed fitted depths on linear and logarithmic grids."""
    depth = catalog.h_eff_m.to_numpy(dtype=float)
    success = catalog.fit_success.fillna(False).to_numpy(dtype=bool)
    valid = catalog.dispersion_valid.fillna(False).to_numpy(dtype=bool)
    masks = dict(screened=catalog.keep_for_histogram.fillna(False).to_numpy(dtype=bool),
                 relaxed=success & valid & np.isfinite(depth) & (depth > 0))
    rows = []
    for selection, keep in masks.items():
        row = dict(**identity, record_days=record_days, wave_type=wave,
                   selection=selection, candidates=len(catalog),
                   selected_events=int(keep.sum()))
        linear = np.histogram(depth[keep], bins=EDGES_M)[0]
        logarithmic = np.histogram(depth[keep], bins=LOG_EDGES_M)[0]
        for left, count in zip(EDGES_M[:-1].astype(int), linear):
            row[f"count_{left:03d}m"] = int(count)
        for index, count in enumerate(logarithmic):
            row[f"log_count_{index:02d}"] = int(count)
        rows.append(row)
    return rows


def main():
    event_frames, spectral_frames, count_rows_all = [], [], []
    for task_index in range(59):
        model, experiment, member, run = simulation(task_index)
        stem = f"{model}.{experiment}.{member}"
        event_path = OUTPUT / f"{stem}.event_energy.csv.gz"
        spectral_path = OUTPUT / f"{stem}.spectral_power.csv.gz"
        if not event_path.is_file() or not spectral_path.is_file():
            raise FileNotFoundError(f"Missing diagnostic output pair: {stem}")
        event = pd.read_csv(event_path)
        spectral = pd.read_csv(spectral_path)
        expected = {(wave, selection) for wave in WAVES for selection in SELECTIONS}
        if (len(event) != 10 or set(zip(event.wave_type, event.selection)) != expected
                or event.duplicated(["wave_type", "selection"]).any()
                or not (event[["model", "experiment", "member"]] ==
                        [model, experiment, member]).all().all()):
            raise RuntimeError(f"Event rows incomplete or misidentified: {stem}")
        expected_spectral = {(wave, source) for wave in WAVES
                             for source in ("observed", "background")}
        if (len(spectral) != 10
                or set(zip(spectral.wave_type, spectral.source)) != expected_spectral
                or spectral.duplicated(["wave_type", "source"]).any()
                or not (spectral[["model", "experiment", "member"]] ==
                        [model, experiment, member]).all().all()):
            raise RuntimeError(f"Spectral rows incomplete or misidentified: {stem}")
        if (not np.isfinite(event.select_dtypes("number").to_numpy()).all()
                or not np.isfinite(spectral.select_dtypes("number").to_numpy()).all()
                or (event.denominator_energy <= 0).any()
                or (spectral.denominator_power <= 0).any()
                or (event.selected_fitted_energy > event.selected_support_energy * (1 + 1e-8)).any()
                or (event.selected_support_energy > event.all_candidate_support_energy * (1 + 1e-8)).any()
                or (event.all_candidate_support_energy > event.filtered_energy * (1 + 1e-8)).any()
                or (spectral.one_to_200_m_power > spectral.valid_depth_power * (1 + 1e-8)).any()
                or (spectral.valid_depth_power > spectral.full_box_power * (1 + 1e-8)).any()):
            raise RuntimeError(f"Energy invariant failed: {stem}")
        for metric, columns in ENERGY_COLUMNS.items():
            if not set(columns).issubset(event):
                raise RuntimeError(f"Missing event energy bins: {stem}/{metric}")
            if (event[columns].sum(axis=1) > event[f"selected_{metric}_energy"] * (1 + 1e-8)).any():
                raise RuntimeError(f"Event bins exceed selected energy: {stem}/{metric}")
        if not set(POWER_COLUMNS).issubset(spectral):
            raise RuntimeError(f"Missing spectral power bins: {stem}")
        if not np.allclose(spectral[POWER_COLUMNS].sum(axis=1),
                           spectral.one_to_200_m_power, rtol=1e-8, atol=1e-4):
            raise RuntimeError(f"Spectral bin accounting failed: {stem}")
        manifest = json.loads((run / "manifest.json").read_text())
        record_days = manifest["record_samples"]
        if not 7_500 <= record_days <= 7_700:
            raise RuntimeError(f"Unexpected record length: {stem}/{record_days}")
        event["record_days"] = record_days
        spectral["record_days"] = record_days
        event_frames.append(event)
        spectral_frames.append(spectral)

        for wave in WAVES:
            catalog = event_catalog(model, experiment, member, wave)
            for row in count_rows(catalog, dict(model=model, experiment=experiment,
                                                member=member), wave, record_days):
                saved = event.loc[(event.wave_type == wave)
                                  & (event.selection == row["selection"])].iloc[0]
                if row["selected_events"] != int(saved.selected_events) or len(catalog) != int(saved.candidates):
                    raise RuntimeError(f"Catalog count differs from energy worker: {stem}/{wave}/{row['selection']}")
                count_rows_all.append(row)
        print(f"Audited {task_index + 1}/59: {stem}", flush=True)

    event_all = pd.concat(event_frames, ignore_index=True)
    spectral_all = pd.concat(spectral_frames, ignore_index=True)
    counts_all = pd.DataFrame(count_rows_all)
    if len(event_all) != 590 or len(spectral_all) != 590 or len(counts_all) != 590:
        raise RuntimeError("Expected 590 rows in each model-wide summary")
    print("PASS: 59 simulations, 295 wave catalogs, 590 rows per model summary", flush=True)

    noaa_manifest = json.loads((NOAA_RUN / "manifest.json").read_text())
    if noaa_manifest["status"] != "complete":
        raise RuntimeError("NOAA calibration manifest is incomplete")
    noaa_days = noaa_manifest["record_samples"]
    noaa_identity = dict(model="NOAA_OLR", experiment="observations", member="not_applicable")
    noaa_event = pd.read_csv(NOAA_RUN / "codex_event_variance/observed.csv.gz")
    expected_noaa = {(wave, selection) for wave in WAVES for selection in SELECTIONS}
    if len(noaa_event) != 10 or set(zip(noaa_event.wave_type, noaa_event.selection)) != expected_noaa:
        raise RuntimeError("NOAA event-energy rows incomplete")
    noaa_event = noaa_event.assign(**noaa_identity, record_days=noaa_days)
    noaa_spectral = pd.DataFrame(spectral_energy_rows(NOAA_RUN, noaa_identity))
    noaa_spectral["record_days"] = noaa_days
    noaa_counts = []
    for wave in WAVES:
        matches = [path for path in DATA_DIR.glob(
            f"NOAA_OLR.{wave}.19900101-20101231.*.ds1p5.component.*.csv")
            if ".pilot" not in path.name]
        if len(matches) != 1:
            raise RuntimeError(f"Expected one full NOAA catalog for {wave}: {matches}")
        catalog = pd.read_csv(matches[0], usecols=[
            "event_id", "h_eff_m", "fit_success", "dispersion_valid",
            "keep_for_histogram",
        ]).set_index("event_id")
        for row in count_rows(catalog, noaa_identity, wave, noaa_days):
            saved = noaa_event.loc[(noaa_event.wave_type == wave)
                                   & (noaa_event.selection == row["selection"])].iloc[0]
            if row["selected_events"] != int(saved.selected_events):
                raise RuntimeError(f"NOAA selected count differs: {wave}/{row['selection']}")
            noaa_counts.append(row)
    write_csv(event_all, OUTPUT / "all_model_event_energy.csv.gz")
    write_csv(spectral_all, OUTPUT / "all_model_spectral_power.csv.gz")
    write_csv(counts_all, OUTPUT / "all_model_event_counts.csv.gz")
    write_csv(noaa_event, OUTPUT / "noaa_event_energy.csv.gz")
    write_csv(noaa_spectral, OUTPUT / "noaa_spectral_power.csv.gz")
    write_csv(pd.DataFrame(noaa_counts), OUTPUT / "noaa_event_counts.csv.gz")
    print("PASS: NOAA event, spectral, and depth-count comparisons", flush=True)


if __name__ == "__main__":
    main()
