#!/usr/bin/env python3
"""Audit full 1.5/2-sigma component catalogs and compare fits and implied depths.

Run on Levante from the project root after all 2-sigma catalogs finish:
    python scripts/codex_sigma2_events/summarize.py

This reads observed/model event catalogs only. It performs no background
calculation or significance selection and never matches events by their IDs
across thresholds, because their component boundaries and IDs have changed.
"""

import copy
import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd

from src.paths import CMIP6_EVENT_DIR, DATA_DIR


TASKLIST = DATA_DIR / "cmip6_event_jobs/tasklists/simulations_waves_20261030_205901_2794174.tsv"
OUTPUT = DATA_DIR / "event_threshold_comparison"
WAVES = ("Kelvin", "n1ERW", "MRG", "n0EIG", "n1WIG")
DEPTH_EDGES_M = np.arange(0.0, 200.1, 5.0)
FLAGS = ("fit_success", "at_search_edge", "dispersion_valid", "fit_quality_ok", "keep_for_histogram")


def read_catalog(model, experiment, member, wave, sigma):
    """Return one unambiguous full catalog, its completed metadata and path."""
    tag = "1p5" if sigma == 1.5 else "2"
    if experiment == "observations":
        directory = DATA_DIR
        pattern = f"NOAA_OLR.{wave}.19900101-20101231.*.ds{tag}.component.*.csv"
    else:
        directory = CMIP6_EVENT_DIR
        pattern = (f"CMIP6.{model}.{experiment}.{member}.rlut.{wave}."
                   f"19900101-2010123[01].*.ds{tag}.component.*.csv")
    matches = [path for path in directory.glob(pattern) if ".pilot" not in path.name]
    if len(matches) != 1:
        raise ValueError(f"Expected one full {sigma}-sigma catalog for {model}/{experiment}/{wave}: {matches}")
    path = matches[0]
    metadata = json.loads(path.with_suffix(".json").read_text())
    cfg = metadata["configuration"]
    identity = cfg["identity"]
    expected_source = "NOAA_CDR" if experiment == "observations" else model
    if (metadata["status"] != "complete" or cfg["max_events"] is not None
            or cfg["fit_support"] != "component" or cfg["wave_type"] != wave
            or cfg["detection"]["detection_sigma"] != sigma
            or cfg["detection"]["peak_sigma"] != sigma
            or (identity["source_id"], identity["experiment_id"], identity["member_id"])
            != (expected_source, experiment, member)):
        raise ValueError(f"Catalog recipe/identity mismatch: {path}")
    columns = [*identity, "wave_type", "event_id", "peak_amplitude_W_m2", "duration_days", "span_deg",
               "K_global", "frequency_cpd", "c0_m_s", "h_eff_m", "weighted_r2", "spatial_cycles",
               "temporal_cycles", "fit_status", *FLAGS]
    frame = pd.read_csv(path, usecols=columns, dtype={name: "boolean" for name in FLAGS})
    if len(frame) != metadata["catalog_candidates"] or frame.event_id.duplicated().any():
        raise ValueError(f"Catalog count or event IDs disagree with metadata: {path}")
    for name, value in {**identity, "wave_type": wave}.items():
        if not frame[name].eq(value).all():
            raise ValueError(f"CSV {name} disagrees with identity: {path}")
    if frame[list(FLAGS)].isna().any().any():
        raise ValueError(f"Missing fit flags: {path}")
    expected_counts = dict(fit_success=int(frame.fit_success.sum()),
                           retained=int(frame.keep_for_histogram.sum()),
                           fit_errors=int(frame.fit_status.eq("fit_error").sum()),
                           incomplete_windows=int(frame.fit_status.eq("incomplete_window").sum()))
    if any(metadata["counts"][name] != amount for name, amount in expected_counts.items()):
        raise ValueError(f"Completed sidecar counts disagree with CSV: {path}")
    return frame, metadata, path


def catalog_statistics(frame, metadata):
    """Describe all candidates and existing fit screens without new selection rules."""
    successful = frame.fit_success.to_numpy(dtype=bool)
    valid = (successful & frame.dispersion_valid.to_numpy(dtype=bool)
             & np.isfinite(frame.h_eff_m) & (frame.h_eff_m > 0)).to_numpy(dtype=bool)
    screened = frame.keep_for_histogram.to_numpy(dtype=bool)
    screens = metadata["configuration"]["fit_screens"]
    expected_quality = (successful & ~frame.at_search_edge
                        & (frame.weighted_r2 >= screens["min_weighted_r2"])
                        & (frame.spatial_cycles >= screens["min_spatial_cycles"])
                        & (frame.temporal_cycles >= screens["min_temporal_cycles"]))
    if (not np.array_equal(frame.fit_quality_ok, expected_quality)
            or not np.array_equal(screened, expected_quality & frame.dispersion_valid)
            or np.any(screened & ~valid)):
        raise ValueError("Saved quality flags disagree with the unchanged catalog screens.")
    if successful.any():
        if not np.isfinite(frame.loc[successful, ["K_global", "frequency_cpd", "weighted_r2"]]).all().all():
            raise ValueError("Successful fits contain nonfinite carrier coordinates or R2.")
    np.testing.assert_allclose(frame.loc[valid, "h_eff_m"],
                               frame.loc[valid, "c0_m_s"] ** 2 / 9.80665, rtol=1e-9, atol=1e-9)
    n = len(frame)
    stats = dict(candidates=n, successful_fits=int(successful.sum()),
                 valid_depths=int(valid.sum()), screened_events=int(screened.sum()),
                 valid_depth_fraction=float(valid.mean()) if n else np.nan,
                 screened_fraction=float(screened.mean()) if n else np.nan,
                 search_edge_fraction=float(frame.at_search_edge.mean()) if n else np.nan,
                 median_fit_r2=float(frame.loc[successful, "weighted_r2"].median()),
                 median_spatial_cycles=float(frame.loc[successful, "spatial_cycles"].median()),
                 median_temporal_cycles=float(frame.loc[successful, "temporal_cycles"].median()),
                 median_duration_days=float(frame.duration_days.median()),
                 median_span_deg=float(frame.span_deg.median()),
                 reference_sigma_W_m2=float(metadata["reference_sigma_W_m2"]),
                 median_peak_in_sigma=float(frame.peak_amplitude_W_m2.median()
                                            / metadata["reference_sigma_W_m2"]))
    for selection, keep in (("valid", valid), ("screened", screened)):
        depth = frame.loc[keep, "h_eff_m"].to_numpy()
        quantiles = np.quantile(depth, [.1, .25, .5, .75, .9]) if len(depth) else np.full(5, np.nan)
        stats.update({f"{selection}_h_{label}_m": float(value)
                      for label, value in zip(("p10", "p25", "median", "p75", "p90"), quantiles)})
        stats[f"{selection}_h_10_80_fraction"] = float(np.mean((depth >= 10) & (depth <= 80))) if len(depth) else np.nan
        stats[f"{selection}_h_above_200_fraction"] = float(np.mean(depth > 200)) if len(depth) else np.nan
    return stats, {"valid": valid, "screened": screened}


def main():
    summary_path = OUTPUT / "catalog_summary.csv"
    hist_path = OUTPUT / "depth_histograms.csv.gz"
    report_path = OUTPUT / "audit.json"
    if any(path.exists() for path in (summary_path, hist_path, report_path)):
        raise FileExistsError(f"Preserving previous threshold-comparison products in {OUTPUT}")
    with TASKLIST.open(newline="") as handle:
        tasks = list(csv.reader(handle, delimiter="\t"))
    simulations = sorted(set(tuple(row[:3]) for row in tasks))
    if len(tasks) != 295 or len(simulations) != 59:
        raise ValueError("Expected the pinned 59-simulation, five-wave cohort.")
    records = [("NOAA_OLR", "observations", "not_applicable"), *simulations]
    summary_rows, histogram_rows = [], []
    for model, experiment, member in records:
        print(f"Auditing {model}/{experiment}/{member}", flush=True)
        for wave in WAVES:
            baseline_cfg = None
            baseline_sigma = None
            for sigma in (1.5, 2.0):
                frame, metadata, path = read_catalog(model, experiment, member, wave, sigma)
                cfg = copy.deepcopy(metadata["configuration"])
                cfg["detection"].pop("detection_sigma")
                cfg["detection"].pop("peak_sigma")
                if sigma == 1.5:
                    baseline_cfg, baseline_sigma = cfg, metadata["reference_sigma_W_m2"]
                else:
                    if json.dumps(cfg, sort_keys=True) != json.dumps(baseline_cfg, sort_keys=True):
                        raise ValueError(f"More than the envelope threshold changed: {path}")
                    if not np.isclose(metadata["reference_sigma_W_m2"], baseline_sigma, rtol=1e-10, atol=1e-10):
                        raise ValueError(f"The filtered-field reference sigma changed: {path}")
                statistics, selections = catalog_statistics(frame, metadata)
                identity = dict(model=model, experiment=experiment, member=member,
                                calendar=cfg["identity"]["calendar"], wave=wave, detection_sigma=sigma)
                summary_rows.append(dict(**identity, **statistics,
                                         catalog=str(path.relative_to(DATA_DIR)),
                                         run_id=path.stem.rsplit(".", 1)[-1]))
                for selection, keep in selections.items():
                    depths = frame.loc[keep, "h_eff_m"].to_numpy()
                    counts = np.histogram(depths, bins=DEPTH_EDGES_M)[0]
                    if counts.sum() != np.count_nonzero((depths >= 0) & (depths <= 200)):
                        raise ValueError(f"Depth histogram loses in-range events: {path}")
                    for lower, upper, count in zip(DEPTH_EDGES_M[:-1], DEPTH_EDGES_M[1:], counts):
                        histogram_rows.append(dict(
                            **identity, selection=selection, lower_m=lower, upper_m=upper,
                            count=int(count), selected_events=len(depths),
                            above_200_events=int(np.count_nonzero(depths > 200)),
                            density_per_m=(count / len(depths) / (upper - lower) if len(depths) else np.nan)))
    summary = pd.DataFrame(summary_rows)
    if len(summary) != 600 or summary.duplicated(["model", "experiment", "member", "wave", "detection_sigma"]).any():
        raise ValueError("The expected 600 catalog records are incomplete or repeated.")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    for path, data in ((summary_path, summary), (hist_path, pd.DataFrame(histogram_rows))):
        temporary = path.with_name(path.name + ".tmp")
        data.to_csv(temporary, index=False, compression="gzip" if path == hist_path else None)
        temporary.replace(path)
    report = dict(
        status="complete", catalog_records=len(summary), source_records=len(records),
        detection_sigma=[1.5, 2.0], fit_support="component", years=[1990, 2010],
        checks="unique complete catalogs, identities, counts, unchanged recipes except thresholds, reference sigma, quality flags, c0/depth units, histogram closure",
        selection="All detected components retained in source catalogs; valid-depth and existing fit-quality subsets summarized separately.",
        caveat="Threshold changes event populations and support; IDs are not paired across thresholds. Fit statistics are descriptive, not background significance or evidence of true depth accuracy.",
    )
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(f"Audited {len(summary)} catalogs; wrote {summary_path} and {hist_path}", flush=True)


if __name__ == "__main__":
    main()
