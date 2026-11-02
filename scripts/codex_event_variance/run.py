#!/usr/bin/env python3
"""Replay NOAA and held-out null events to measure variance on active supports.

Run from the project root on Levante:
    python scripts/codex_event_variance/run.py --observed
    python scripts/codex_event_variance/run.py --start 1000 --stop 1025

This reads fitted depths from the existing component catalogs. It does not
recalibrate the background or refit sinusoids. Output CSVs live under the
saved calibration's codex_event_variance/ directory.
"""

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
import xr_ccew as tw

from src import wave_events as events
from src import wave_significance as significance
from src.paths import DATA_DIR, NOAA_OLR_FILE, PROJECT_ROOT


RUN = (DATA_DIR / "wave_significance" / "gaussian_calibration_wk99_121"
       / "NOAA_OLR.observations.not_applicable.1990-2010.b898b1ae4307")
OUTPUT = RUN / "codex_event_variance"
WAVES = ("Kelvin", "MRG", "n0EIG", "n1ERW", "n1WIG")
EDGES_M = np.arange(1.0, 201.0, 1.0)


def preprocessed_anomaly(field):
    """Return the catalog's tapered (time, lon) projected OLR anomaly, W m-2."""
    return tw.apply_window(
        tw.remove_mean_and_linear_trend(
            tw.remove_harmonics_of_seasonal_cycle(
                field, num_harmonics=events.SEASONAL_HARMONICS)),
        window="tukey", pct=events.TAPER_FRACTION,
    )


def observed_catalog(wave):
    matches = [path for path in DATA_DIR.glob(
        f"NOAA_OLR.{wave}.19900101-20101231.*.ds1p5.component.*.csv")
        if ".pilot" not in path.name]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one full observed catalog for {wave}: {matches}")
    return matches[0]


def catalog_rows(path):
    frame = pd.read_csv(path, usecols=[
        "event_id", "peak_amplitude_W_m2", "h_eff_m", "fit_success",
        "dispersion_valid", "keep_for_histogram", "weighted_r2",
    ])
    if frame.event_id.duplicated().any():
        raise RuntimeError(f"Duplicate event IDs: {path}")
    return frame.set_index("event_id")


def variance_rows(field, wave, catalog, *, realization, denominator):
    """Bin disjoint event-support and fitted-carrier energy by fitted depth."""
    detected, _, support = events.detect_events(field, detection_sigma=1.5,
                                                 return_support=True)
    detected_ids = {event["event_id"] for event in detected}
    if detected_ids != set(catalog.index):
        raise RuntimeError(f"Detected IDs differ from catalog: {wave}/{realization}")
    peaks = pd.Series({event["event_id"]: event["peak_amplitude_W_m2"]
                       for event in detected})
    if not np.allclose(peaks.sort_index().to_numpy(),
                       catalog.peak_amplitude_W_m2.sort_index().to_numpy(),
                       rtol=1e-9, atol=1e-8):
        raise RuntimeError(f"Detected peaks differ from catalog: {wave}/{realization}")

    values = field.to_numpy()
    field_mean = float(values.mean())
    filtered_energy = float(np.sum((values - field_mean) ** 2))
    depth = catalog.h_eff_m.to_numpy(dtype=float)
    fit_success = catalog.fit_success.fillna(False).to_numpy(dtype=bool)
    dispersion_valid = catalog.dispersion_valid.fillna(False).to_numpy(dtype=bool)
    screened = catalog.keep_for_histogram.fillna(False).to_numpy(dtype=bool)
    relaxed = fit_success & dispersion_valid & np.isfinite(depth) & (depth > 0)
    if np.any(screened & ~relaxed):
        raise RuntimeError(f"Screened event lacks valid depth: {wave}/{realization}")

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
        raise RuntimeError(f"Event supports overlap or exceed field energy: {wave}/{realization}")
    if fitted_energy.sum() > support_energy.sum() * (1 + 1e-8):
        raise RuntimeError(f"Fitted energy exceeds support energy: {wave}/{realization}")

    rows = []
    for selection, keep in (("screened", screened), ("relaxed", relaxed)):
        row = dict(wave_type=wave, realization=realization, selection=selection,
                   candidates=len(catalog), selected_events=int(keep.sum()),
                   denominator_energy=denominator, filtered_energy=filtered_energy,
                   all_candidate_support_energy=float(support_energy.sum()),
                   selected_support_energy=float(support_energy[keep].sum()),
                   selected_fitted_energy=float(fitted_energy[keep].sum()))
        for label, weights in (("support", support_energy), ("fitted", fitted_energy)):
            histogram = np.histogram(depth[keep], bins=EDGES_M,
                                     weights=weights[keep])[0]
            for left, amount in zip(EDGES_M[:-1].astype(int), histogram):
                row[f"{label}_{left:03d}m"] = float(amount)
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--observed", action="store_true")
    mode.add_argument("--start", type=int)
    parser.add_argument("--stop", type=int)
    args = parser.parse_args()
    if PROJECT_ROOT != Path(__file__).resolve().parents[2]:
        raise RuntimeError("Install this checkout with pip -e . before running")
    manifest = json.loads((RUN / "manifest.json").read_text())
    cfg = manifest["configuration"]
    if (manifest["status"] != "complete" or cfg["background_method"] != "wk99_121"
            or cfg["detection"]["detection_sigma"] != 1.5
            or cfg["n_validation"] != 1000):
        raise RuntimeError("Expected the completed 1.5-sigma NOAA calibration")
    for wave in WAVES:
        if json.dumps(cfg["profiles"][wave], sort_keys=True) != json.dumps(
                events.WAVE_PROFILES[wave].as_dict(), sort_keys=True):
            raise RuntimeError(f"Profile changed since calibration: {wave}")
    if args.observed:
        if args.stop is not None:
            parser.error("--stop is only valid with --start")
        realizations = (-1,)
    else:
        if args.stop is None or not cfg["n_calibration"] <= args.start < args.stop <= cfg["n_calibration"] + 1000:
            parser.error("Use a held-out half-open realization range within [1000, 2000)")
        realizations = range(args.start, args.stop)

    index = pd.read_csv(RUN / "component_fit_catalogs" / "catalog_index.1000-2000.csv")
    if len(index) != 5000:
        raise RuntimeError("Expected 5,000 audited fitted null catalogs")
    catalog_paths = {(row.wave_type, row.realization): PROJECT_ROOT / row.catalog
                     for row in index.itertuples(index=False)}
    templates = {}
    power = {}
    for symmetry in ("symmetric", "antisymmetric"):
        with xr.open_dataset(RUN / f"template.{symmetry}.nc") as source:
            templates[symmetry] = source[symmetry].load()
        with xr.open_dataset(RUN / f"spectrum.{symmetry}.nc") as source:
            power[symmetry] = source.background_power.load()

    observed_olr = None
    projected_observed = None
    if args.observed:
        if not NOAA_OLR_FILE.parent.is_symlink() or not NOAA_OLR_FILE.is_file():
            raise FileNotFoundError("Configure data/noaa_olr for the observed replay")
        with xr.open_dataset(NOAA_OLR_FILE) as source:
            observed_olr = source.olr.sel(
                time=slice("1990-01-01", "2010-12-31"),
                lat=slice(*events.LAT_RANGE)).load()
        projected_observed = significance.project_olr(observed_olr)

    rows = []
    for realization in realizations:
        anomalies = {}
        for symmetry in templates:
            if realization == -1:
                source = projected_observed[symmetry]
            else:
                seed = np.random.SeedSequence([
                    cfg["seed"], 0 if symmetry == "symmetric" else 1, realization])
                source = significance.draw_background(
                    power[symmetry], templates[symmetry], np.random.default_rng(seed))
            anomalies[symmetry] = preprocessed_anomaly(source)
        denominator = float(sum(np.sum(anomaly.to_numpy() ** 2)
                                for anomaly in anomalies.values()))
        if denominator <= 0:
            raise RuntimeError(f"Nonpositive projected variance: {realization}")
        for wave in WAVES:
            symmetry = events.WAVE_PROFILES[wave].symmetry
            if realization == -1:
                # Reproduce the exact observed catalog path, including the
                # latitude projection/filter ordering, rather than reuse the
                # nearly identical calibration snapshot.
                filtered = events.filter_olr_for_wave(observed_olr, wave)
                path = observed_catalog(wave)
            else:
                filtered = tw.filter_field(anomalies[symmetry],
                                           events.WAVE_PROFILES[wave],
                                           apply_profile_symmetry=False)
                path = catalog_paths[wave, realization]
            catalog = catalog_rows(path)
            rows.extend(variance_rows(filtered, wave, catalog,
                                      realization=realization,
                                      denominator=denominator))
            print(f"{realization} {wave}: {len(catalog)} matched events", flush=True)

    OUTPUT.mkdir(parents=True, exist_ok=True)
    name = "observed.csv.gz" if args.observed else f"null.{args.start:04d}-{args.stop:04d}.csv.gz"
    output = OUTPUT / name
    if output.exists():
        raise FileExistsError(f"Preserve existing diagnostic output: {output}")
    temporary = output.with_name(output.name + ".tmp")
    pd.DataFrame(rows).to_csv(temporary, index=False, compression="gzip")
    temporary.replace(output)
    print(f"Wrote {len(rows)} record/selection rows to {output}", flush=True)


if __name__ == "__main__":
    main()
