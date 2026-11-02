#!/usr/bin/env python3
"""Measure CMIP6 event-footprint energy and spectral power by equivalent depth.

Run from the project root on Levante after installing this checkout:
    python scripts/codex_cmip_heff/run.py --task-index 0

Each task is one model/experiment/member in the saved 59-simulation task list.
The 1.5-sigma, active-component event catalogs and WK99-style spectra must
already exist. Outputs are compact CSVs in data/wave_significance/.
"""

import argparse
import csv
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import xr_ccew as tw

from src.cmip import load_cmip_rlut
from src.paths import CMIP6_EVENT_DIR, DATA_DIR, PROJECT_ROOT
from src import wave_events as events
from src import wave_significance as significance


TASKLIST = DATA_DIR / "wave_significance/cmip6_calibration_jobs/simulations_1990-2010_wk99_ds1p5.tsv"
INDEX = DATA_DIR / "wave_significance/cmip6_calibration_summaries/index.json"
RUNS = DATA_DIR / "wave_significance/gaussian_calibration_wk99_121_cmip6"
OUTPUT = DATA_DIR / "wave_significance/codex_cmip_heff_comparison"
WAVES = ("Kelvin", "MRG", "n0EIG", "n1ERW", "n1WIG")
EDGES_M = np.arange(1.0, 201.0)
EARTH_RADIUS_M = 6_371_008.0
GRAVITY_M_S2 = 9.80665
BETA_S_INV_M_INV = 2 * 7.292e-5 / EARTH_RADIUS_M


def simulation(task_index):
    with TASKLIST.open(newline="") as source:
        simulations = list(csv.reader(source, delimiter="\t"))
    if len(simulations) != 59 or not 0 <= task_index < len(simulations):
        raise ValueError("Expected a task index for the 59 completed CMIP6 simulations")
    model, experiment, member = simulations[task_index]
    index = json.loads(INDEX.read_text())
    if index.get("complete_runs") != 59:
        raise RuntimeError("The 59 calibration summaries are not complete")
    names = [name for name in index["runs"]
             if name.startswith(f"{model}.{experiment}.{member}.1990-2010.")]
    if len(names) != 1:
        raise RuntimeError(f"Expected one calibration for {model}/{experiment}/{member}: {names}")
    run = RUNS / names[0]
    manifest = json.loads((run / "manifest.json").read_text())
    cfg = manifest["configuration"]
    identity = cfg["identity"]
    if (manifest["status"] != "complete" or cfg["background_method"] != "wk99_121"
            or cfg["detection"]["detection_sigma"] != 1.5
            or cfg["n_validation"] != 1000
            or (identity["source_id"], identity["experiment_id"], identity["member_id"])
            != (model, experiment, member)):
        raise RuntimeError(f"Unexpected calibration configuration: {run}")
    for wave in WAVES:
        if json.dumps(cfg["profiles"][wave], sort_keys=True) != json.dumps(
                events.WAVE_PROFILES[wave].as_dict(), sort_keys=True):
            raise RuntimeError(f"Profile differs from calibration: {wave}")
    return model, experiment, member, run


def event_catalog(model, experiment, member, wave):
    matches = [path for path in CMIP6_EVENT_DIR.glob(
        f"CMIP6.{model}.{experiment}.{member}.rlut.{wave}.*.ds1p5.component.*.csv")
        if ".pilot" not in path.name]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one full component catalog for {model}/{experiment}/{wave}: {matches}")
    path = matches[0]
    sidecar = json.loads(path.with_suffix(".json").read_text())
    cfg = sidecar["configuration"]
    if (sidecar["status"] != "complete" or cfg["fit_support"] != "component"
            or cfg["max_events"] is not None or cfg["detection"]["detection_sigma"] != 1.5
            or cfg["wave_type"] != wave
            or json.dumps(cfg["profile"], sort_keys=True) != json.dumps(
                events.WAVE_PROFILES[wave].as_dict(), sort_keys=True)):
        raise RuntimeError(f"Unexpected fitted catalog configuration: {path}")
    frame = pd.read_csv(path, usecols=[
        "event_id", "peak_amplitude_W_m2", "h_eff_m", "fit_success",
        "dispersion_valid", "keep_for_histogram", "weighted_r2",
    ])
    if frame.event_id.duplicated().any() or len(frame) != sidecar["catalog_candidates"]:
        raise RuntimeError(f"Fitted catalog count/IDs failed: {path}")
    return frame.set_index("event_id")


def preprocessed_anomaly(field):
    """Match the saved calibration's daily taper and detrending, W m-2."""
    return tw.apply_window(
        tw.remove_mean_and_linear_trend(
            tw.remove_harmonics_of_seasonal_cycle(
                field, num_harmonics=events.SEASONAL_HARMONICS)),
        window="tukey", pct=events.TAPER_FRACTION,
    )


def event_energy_rows(field, wave, catalog, denominator, identity):
    """Bin disjoint active-component energy by each saved fitted depth."""
    detected, _, support = events.detect_events(field, detection_sigma=1.5,
                                                 return_support=True)
    if {event["event_id"] for event in detected} != set(catalog.index):
        raise RuntimeError(f"Detected event IDs differ from fitted catalog: {identity}/{wave}")
    peaks = pd.Series({event["event_id"]: event["peak_amplitude_W_m2"]
                       for event in detected})
    if not np.allclose(peaks.sort_index().to_numpy(),
                       catalog.peak_amplitude_W_m2.sort_index().to_numpy(),
                       rtol=1e-9, atol=1e-8):
        raise RuntimeError(f"Detected peaks differ from fitted catalog: {identity}/{wave}")

    values = field.to_numpy()
    field_mean = float(values.mean())
    filtered_energy = float(np.sum((values - field_mean) ** 2))
    depth = catalog.h_eff_m.to_numpy(dtype=float)
    fit_success = catalog.fit_success.fillna(False).to_numpy(dtype=bool)
    dispersion_valid = catalog.dispersion_valid.fillna(False).to_numpy(dtype=bool)
    screened = catalog.keep_for_histogram.fillna(False).to_numpy(dtype=bool)
    relaxed = fit_success & dispersion_valid & np.isfinite(depth) & (depth > 0)
    if np.any(screened & ~relaxed):
        raise RuntimeError(f"Screened event lacks valid depth: {identity}/{wave}")

    support_energy = np.empty(len(catalog), dtype=float)
    fitted_energy = np.empty(len(catalog), dtype=float)
    for j, event_id in enumerate(catalog.index):
        ti, xi = support[event_id]
        samples = values[ti, xi]
        support_energy[j] = np.sum((samples - field_mean) ** 2)
        local_sst = np.sum((samples - samples.mean()) ** 2)
        r2 = catalog.weighted_r2.iloc[j]
        fitted_energy[j] = max(0.0, min(1.0, r2)) * local_sst if np.isfinite(r2) else 0.0
    if support_energy.sum() > filtered_energy * (1 + 1e-8):
        raise RuntimeError(f"Event footprints exceed filtered energy: {identity}/{wave}")
    if fitted_energy.sum() > support_energy.sum() * (1 + 1e-8):
        raise RuntimeError(f"Fitted energy exceeds event energy: {identity}/{wave}")

    rows = []
    for selection, keep in (("screened", screened), ("relaxed", relaxed)):
        row = dict(**identity, wave_type=wave, selection=selection,
                   candidates=len(catalog), selected_events=int(keep.sum()),
                   denominator_energy=denominator, filtered_energy=filtered_energy,
                   all_candidate_support_energy=float(support_energy.sum()),
                   selected_support_energy=float(support_energy[keep].sum()),
                   selected_fitted_energy=float(fitted_energy[keep].sum()))
        for label, weights in (("support", support_energy), ("fitted", fitted_energy)):
            histogram = np.histogram(depth[keep], bins=EDGES_M, weights=weights[keep])[0]
            for left, amount in zip(EDGES_M[:-1].astype(int), histogram):
                row[f"{label}_{left:03d}m"] = float(amount)
        rows.append(row)
    return rows


def spectral_equivalent_depth(wave, global_k, frequency_cpd):
    """Vectorized resting-branch dispersion inversion on positive physical f."""
    k = global_k / EARTH_RADIUS_M
    omega = 2 * np.pi * frequency_cpd / 86_400.0
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        if wave == "Kelvin":
            c0 = omega / k
        elif wave in ("MRG", "n0EIG"):
            denominator = BETA_S_INV_M_INV + k * omega
            c0 = np.where(denominator > 0, omega**2 / denominator, np.nan)
        else:
            coefficient = k**2 + BETA_S_INV_M_INV * k / omega
            discriminant = (3 * BETA_S_INV_M_INV)**2 + 4 * coefficient * omega**2
            if wave == "n1WIG":
                c0 = np.where(discriminant >= 0,
                              2 * omega**2 / (3 * BETA_S_INV_M_INV + np.sqrt(discriminant)),
                              np.nan)
            elif wave == "n1ERW":
                c0 = np.where((coefficient < 0) & (discriminant >= 0),
                              (-3 * BETA_S_INV_M_INV - np.sqrt(discriminant)) / (2 * coefficient),
                              np.nan)
            else:
                raise ValueError(wave)
        depth = c0**2 / GRAVITY_M_S2
    return np.where(np.isfinite(depth) & (depth > 0), depth, np.nan)


def spectral_energy_rows(run, identity):
    """Bin observed and smooth-background Fourier power by apparent depth."""
    spectra = {}
    for symmetry in ("symmetric", "antisymmetric"):
        with xr.open_dataset(run / f"spectrum.{symmetry}.nc") as dataset:
            spectra[symmetry] = dict(
                frequency=dataset.fft_frequency.to_numpy(),
                wavenumber=dataset.fft_wavenumber.to_numpy(),
                observed=dataset.observed_taper_corrected_power.to_numpy(),
                background=dataset.background_power.to_numpy(),
            )
        spectrum = spectra[symmetry]
        for source in ("observed", "background"):
            power = spectrum[source]
            if (power.shape != (len(spectrum["frequency"]), len(spectrum["wavenumber"]))
                    or not np.isfinite(power).all() or (power < 0).any()):
                raise RuntimeError(f"Invalid saved {symmetry}/{source} spectrum: {run}")
    if (not np.array_equal(spectra["symmetric"]["frequency"], spectra["antisymmetric"]["frequency"])
            or not np.array_equal(spectra["symmetric"]["wavenumber"], spectra["antisymmetric"]["wavenumber"])
            or not np.array_equal(spectra["symmetric"]["background"], spectra["antisymmetric"]["background"])):
        raise RuntimeError(f"Symmetric/antisymmetric spectral grids or common background differ: {run}")
    totals = {source: sum(spectra[symmetry][source].sum() for symmetry in spectra)
              for source in ("observed", "background")}
    if min(totals.values()) <= 0:
        raise RuntimeError(f"Nonpositive saved projected spectral power: {run}")

    rows = []
    for wave in WAVES:
        profile = events.WAVE_PROFILES[wave]
        spectrum = spectra[profile.symmetry]
        # Physical positive f maps to negative temporal NumPy FFT frequencies.
        frequency_indices = np.flatnonzero(
            (spectrum["frequency"] < 0)
            & (-spectrum["frequency"] >= profile.frequency_min)
            & (-spectrum["frequency"] <= profile.frequency_max))
        # Keep floating FFT coordinates: rounding can admit a filter-edge K.
        wavenumber_indices = np.flatnonzero(
            (spectrum["wavenumber"] >= profile.k_min)
            & (spectrum["wavenumber"] <= profile.k_max))
        if not len(frequency_indices) or not len(wavenumber_indices):
            raise RuntimeError(f"Empty spectral box: {identity}/{wave}")
        depth = spectral_equivalent_depth(
            wave, spectrum["wavenumber"][None, wavenumber_indices],
            -spectrum["frequency"][frequency_indices, None])
        valid = np.isfinite(depth)
        for source in ("observed", "background"):
            # Real-field conjugates carry the same power.
            power = 2 * spectrum[source][np.ix_(frequency_indices, wavenumber_indices)]
            histogram = np.histogram(depth[valid], bins=EDGES_M, weights=power[valid])[0]
            row = dict(**identity, wave_type=wave, source=source,
                       denominator_power=float(totals[source]),
                       full_box_power=float(power.sum()),
                       valid_depth_power=float(power[valid].sum()),
                       one_to_200_m_power=float(histogram.sum()),
                       twelve_to_50_m_power=float(power[
                           valid & (depth >= 12) & (depth < 50)].sum()))
            for left, amount in zip(EDGES_M[:-1].astype(int), histogram):
                row[f"power_{left:03d}m"] = float(amount)
            rows.append(row)
    return rows


def write_csv(rows, path):
    if path.exists():
        raise FileExistsError(f"Preserving existing diagnostic: {path}")
    temporary = path.with_name(path.name + ".tmp")
    pd.DataFrame(rows).to_csv(temporary, index=False, compression="gzip")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-index", type=int,
                        default=int(os.environ.get("SLURM_ARRAY_TASK_ID", "-1")))
    args = parser.parse_args()
    if PROJECT_ROOT != Path(__file__).resolve().parents[2]:
        raise RuntimeError("Install this checkout with pip -e . before running")
    model, experiment, member, run = simulation(args.task_index)
    identity = dict(model=model, experiment=experiment, member=member)
    stem = f"{model}.{experiment}.{member}"
    OUTPUT.mkdir(parents=True, exist_ok=True)
    event_path = OUTPUT / f"{stem}.event_energy.csv.gz"
    spectral_path = OUTPUT / f"{stem}.spectral_power.csv.gz"
    if event_path.exists() or spectral_path.exists():
        if event_path.exists() and spectral_path.exists():
            print(f"Already complete: {stem}", flush=True)
            return
        raise RuntimeError(f"Partial output for {stem}; inspect before retrying")

    olr = load_cmip_rlut(model, experiment=experiment, member_id=member,
                         year_range=(1990, 2010), lat_range=events.LAT_RANGE)
    components = significance.project_olr(olr)
    denominator = float(sum(np.sum(preprocessed_anomaly(components[symmetry]).to_numpy() ** 2)
                            for symmetry in ("symmetric", "antisymmetric")))
    if denominator <= 0:
        raise RuntimeError(f"Nonpositive projected OLR energy: {stem}")
    event_rows = []
    for wave in WAVES:
        filtered = events.filter_olr_for_wave(olr, wave)
        catalog = event_catalog(model, experiment, member, wave)
        event_rows.extend(event_energy_rows(filtered, wave, catalog, denominator, identity))
        print(f"{stem} {wave}: matched {len(catalog)} events", flush=True)
    spectral_rows = spectral_energy_rows(run, identity)
    write_csv(event_rows, event_path)
    write_csv(spectral_rows, spectral_path)
    print(f"Wrote {event_path} and {spectral_path}", flush=True)


if __name__ == "__main__":
    main()
