#!/usr/bin/env python3
"""Bin saved Kelvin event-support variance and Fourier power by fitted K and f.

Run on Levante from the project root with --task-index 0..59; index 0 is NOAA.
The 59 model records follow kelvin_catalog_index.csv from extract.py.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from src import wave_events as events
from src.cmip import load_cmip_rlut
from src.paths import DATA_DIR


ROOT = DATA_DIR / "wave_significance/codex_kelvin_komega"
CMIP_RUNS = DATA_DIR / "wave_significance/gaussian_calibration_wk99_121_cmip6"
CMIP_INDEX = DATA_DIR / "wave_significance/cmip6_calibration_summaries/index.json"
NOAA_RUN = (DATA_DIR / "wave_significance/gaussian_calibration_wk99_121"
            / "NOAA_OLR.observations.not_applicable.1990-2010.b898b1ae4307")
CMIP_SUMMARY = DATA_DIR / "wave_significance/codex_cmip_heff_comparison"
K_EDGES = np.arange(0.5, 16.0, 1.0)
F_EDGES = np.arange(0.07, 0.4101, 0.02)


def calibration_run(record):
    if record.experiment == "observations":
        return NOAA_RUN
    index = json.loads(CMIP_INDEX.read_text())
    if index["complete_runs"] != 59:
        raise RuntimeError("Expected 59 completed CMIP6 calibrations")
    prefix = f"{record.model}.{record.experiment}.{record.member}.1990-2010."
    names = [name for name in index["runs"] if name.startswith(prefix)]
    if len(names) != 1:
        raise RuntimeError(f"Expected one calibration run for {prefix}: {names}")
    return CMIP_RUNS / names[0]


def summary_row(path, record, *, source=None):
    frame = pd.read_csv(path)
    match = (frame.model.eq(record.model) & frame.experiment.eq(record.experiment)
             & frame.member.eq(record.member) & frame.wave_type.eq("Kelvin"))
    if source is None:
        match &= frame.selection.eq("screened")
    else:
        match &= frame.source.eq(source)
    rows = frame.loc[match]
    if len(rows) != 1:
        raise RuntimeError(f"Expected one Kelvin summary row in {path}: {record.model}/{source}")
    return rows.iloc[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-index", type=int, required=True)
    args = parser.parse_args()
    index = pd.read_csv(ROOT / "kelvin_catalog_index.csv")
    if len(index) != 60 or not 0 <= args.task_index < 60:
        raise ValueError("Expected a task index in 0..59 for the audited Kelvin records")
    record = index.iloc[args.task_index]
    run = calibration_run(record)
    manifest = json.loads((run / "manifest.json").read_text())
    if (manifest["status"] != "complete"
            or manifest["configuration"]["detection"]["detection_sigma"] != 1.5
            or manifest["configuration"]["background_method"] != "wk99_121"):
        raise RuntimeError(f"Unexpected calibration configuration: {run}")
    catalog = pd.read_csv(DATA_DIR / record.catalog, usecols=[
        "event_id", "peak_amplitude_W_m2", "K_global", "frequency_cpd",
        "fit_success", "keep_for_histogram", "weighted_r2",
    ]).set_index("event_id")
    if len(catalog) != record.candidates or catalog.index.has_duplicates:
        raise RuntimeError(f"Catalog differs from extraction index: {record.catalog}")

    with xr.open_dataset(run / "observed_filtered.Kelvin.nc") as source:
        field = source.filtered_olr.load()
    detected, _, support = events.detect_events(field, detection_sigma=1.5, return_support=True)
    replayed_raw = False
    if {item["event_id"] for item in detected} != set(catalog.index):
        if record.experiment == "observations":
            raise RuntimeError(f"Saved filtered field does not reproduce NOAA event IDs: {run}")
        print(f"Saved filtered field changes event IDs; replaying catalog filter from raw OLR: {run}",
              flush=True)
        olr = load_cmip_rlut(record.model, experiment=record.experiment,
                             member_id=record.member, year_range=(1990, 2010),
                             lat_range=events.LAT_RANGE)
        field = events.filter_olr_for_wave(olr, "Kelvin")
        detected, _, support = events.detect_events(field, detection_sigma=1.5, return_support=True)
        replayed_raw = True
        if {item["event_id"] for item in detected} != set(catalog.index):
            raise RuntimeError(f"Raw-input replay still differs from fitted event IDs: {run}")
    peaks = pd.Series({item["event_id"]: item["peak_amplitude_W_m2"] for item in detected})
    if not np.allclose(peaks.sort_index(), catalog.peak_amplitude_W_m2.sort_index(),
                       rtol=1e-9 if replayed_raw else 1e-5,
                       atol=1e-8 if replayed_raw else 5e-5):
        raise RuntimeError(f"Saved filtered field does not reproduce event peaks: {run}")

    values = field.to_numpy()
    field_mean = float(values.mean())
    support_energy = np.empty(len(catalog))
    fitted_energy = np.empty(len(catalog))
    for j, (event_id, row) in enumerate(catalog.iterrows()):
        ti, xi = support[event_id]
        samples = values[ti, xi]
        support_energy[j] = np.sum((samples - field_mean) ** 2)
        local_sst = np.sum((samples - samples.mean()) ** 2)
        r2 = row.weighted_r2
        fitted_energy[j] = max(0.0, min(1.0, r2)) * local_sst if np.isfinite(r2) else 0.0
    if support_energy.sum() > np.sum((values - field_mean) ** 2) * (1 + 1e-8):
        raise RuntimeError("Event supports overlap or exceed filtered-field energy")

    summary_file = ("noaa_event_energy.csv.gz" if record.experiment == "observations"
                    else "all_model_event_energy.csv.gz")
    summary = summary_row(CMIP_SUMMARY / summary_file, record)
    denominator = float(summary.denominator_energy)
    if (not np.isclose(support_energy.sum(), summary.all_candidate_support_energy, rtol=1e-4)
            or not np.isclose(support_energy[catalog.keep_for_histogram.to_numpy()].sum(),
                               summary.selected_support_energy, rtol=1e-4)):
        raise RuntimeError(f"Saved event-support energy does not match prior audit: {run}")

    with xr.open_dataset(run / "spectrum.symmetric.nc") as source:
        frequency = source.fft_frequency.to_numpy()
        wavenumber = source.fft_wavenumber.to_numpy()
        fi = np.flatnonzero((frequency < 0) & (-frequency >= 1/14) & (-frequency <= 0.4))
        ki = np.flatnonzero((wavenumber >= 1) & (wavenumber <= 15))
        observed = 2 * source.observed_taper_corrected_power.isel(
            fft_frequency=fi, fft_wavenumber=ki).to_numpy()
        background = 2 * source.background_power.isel(
            fft_frequency=fi, fft_wavenumber=ki).to_numpy()
    K, f = np.meshgrid(wavenumber[ki], -frequency[fi])
    observed_bins = np.histogram2d(K.ravel(), f.ravel(), bins=(K_EDGES, F_EDGES),
                                   weights=observed.ravel())[0]
    background_bins = np.histogram2d(K.ravel(), f.ravel(), bins=(K_EDGES, F_EDGES),
                                     weights=background.ravel())[0]
    spectral_file = ("noaa_spectral_power.csv.gz" if record.experiment == "observations"
                     else "all_model_spectral_power.csv.gz")
    spectral_observed = summary_row(CMIP_SUMMARY / spectral_file, record, source="observed")
    spectral_background = summary_row(CMIP_SUMMARY / spectral_file, record, source="background")
    if (not np.isclose(observed_bins.sum(), spectral_observed.full_box_power, rtol=1e-8)
            or not np.isclose(background_bins.sum(), spectral_background.full_box_power, rtol=1e-8)):
        raise RuntimeError(f"Binned Fourier power disagrees with prior spectral audit: {run}")

    rows = []
    n_grid = values.size
    for selection, keep in (("all_successful", catalog.fit_success.to_numpy()),
                            ("screened", catalog.keep_for_histogram.to_numpy())):
        K_fit = catalog.K_global.to_numpy()[keep]
        f_fit = catalog.frequency_cpd.to_numpy()[keep]
        support_bins = np.histogram2d(K_fit, f_fit, bins=(K_EDGES, F_EDGES),
                                      weights=support_energy[keep])[0]
        fitted_bins = np.histogram2d(K_fit, f_fit, bins=(K_EDGES, F_EDGES),
                                     weights=fitted_energy[keep])[0]
        if not np.isclose(support_bins.sum(), support_energy[keep].sum(), rtol=1e-10):
            raise RuntimeError(f"Event energy falls outside Kelvin K-f bins: {run}/{selection}")
        for i, k_left in enumerate(K_EDGES[:-1]):
            for j, f_lower in enumerate(F_EDGES[:-1]):
                rows.append(dict(model=record.model, experiment=record.experiment,
                                 member=record.member, selection=selection,
                                 selected_events=int(keep.sum()),
                                 K_left=k_left, K_right=K_EDGES[i+1],
                                 f_lower_cpd=f_lower, f_upper_cpd=F_EDGES[j+1],
                                 support_variance_W2_m4=support_bins[i, j]/n_grid,
                                 fitted_variance_W2_m4=fitted_bins[i, j]/n_grid,
                                 support_pct_projected=100*support_bins[i, j]/denominator,
                                 fitted_pct_projected=100*fitted_bins[i, j]/denominator,
                                 spectral_variance_W2_m4=observed_bins[i, j]/n_grid,
                                 spectral_background_variance_W2_m4=background_bins[i, j]/n_grid,
                                 spectral_pct_projected=100*observed_bins[i, j]
                                                        / spectral_observed.denominator_power,
                                 spectral_observed_to_background=(
                                     observed_bins[i, j]/background_bins[i, j]
                                     if background_bins[i, j] > 0 else np.nan)))
    output_dir = ROOT / "variance_bins_integer_centered"
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / f"{record.model}.{record.experiment}.{record.member}.csv.gz"
    if path.exists():
        print(f"Already complete: {path}", flush=True)
        return
    temporary = path.with_name(path.name + ".tmp")
    pd.DataFrame(rows).to_csv(temporary, index=False, compression="gzip")
    temporary.replace(path)
    print(f"Wrote {path}: {len(catalog)} candidates, {int(catalog.keep_for_histogram.sum())} screened", flush=True)


if __name__ == "__main__":
    main()
