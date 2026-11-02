#!/usr/bin/env python3
"""Extract fitted Kelvin event coordinates from the completed 1990–2010 catalogs.

Run from the project root on Levante after installing this checkout:
    python scripts/codex_kelvin_komega/extract.py

Inputs are the full 1.5σ active-component NOAA and CMIP6 Kelvin CSVs under
data/. Outputs are compact tables under data/wave_significance/codex_kelvin_komega/.
"""

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd

from src.paths import CMIP6_EVENT_DIR, DATA_DIR


TASKLIST = (DATA_DIR / "wave_significance/cmip6_calibration_jobs"
            / "simulations_1990-2010_wk99_ds1p5.tsv")
OUTPUT = DATA_DIR / "wave_significance/codex_kelvin_komega"
EVENT_COLUMNS = [
    "event_id", "K_global", "frequency_cpd", "k_rad_m", "omega_rad_s",
    "fit_success", "keep_for_histogram", "weighted_r2", "at_search_edge",
]
EARTH_RADIUS_M = 6_371_008.0


def full_catalog(model, experiment, member):
    if model == "NOAA_OLR":
        directory = DATA_DIR
        pattern = "NOAA_OLR.Kelvin.19900101-20101231.*.ds1p5.component.*.csv"
    else:
        directory = CMIP6_EVENT_DIR
        pattern = (f"CMIP6.{model}.{experiment}.{member}.rlut.Kelvin."
                   "19900101-2010123[01].*.ds1p5.component.*.csv")
    matches = [path for path in directory.glob(pattern) if ".pilot" not in path.name]
    if len(matches) != 1:
        raise RuntimeError(f"Expected one full Kelvin catalog for {model}/{experiment}/{member}: {matches}")
    path = matches[0]
    sidecar = json.loads(path.with_suffix(".json").read_text())
    configuration = sidecar["configuration"]
    if (sidecar["status"] != "complete" or configuration["wave_type"] != "Kelvin"
            or configuration["fit_support"] != "component"
            or configuration["max_events"] is not None
            or configuration["detection"]["detection_sigma"] != 1.5):
        raise RuntimeError(f"Unexpected catalog configuration: {path}")
    frame = pd.read_csv(path, usecols=EVENT_COLUMNS)
    if len(frame) != sidecar["catalog_candidates"] or frame.event_id.duplicated().any():
        raise RuntimeError(f"Catalog row count or event IDs differ from sidecar: {path}")
    good = frame.fit_success.fillna(False)
    if (frame.keep_for_histogram.fillna(False) & ~good).any():
        raise RuntimeError(f"Quality-screened fit is not successful: {path}")
    if not np.isfinite(frame.loc[good, ["K_global", "frequency_cpd", "k_rad_m", "omega_rad_s"]]).all().all():
        raise RuntimeError(f"Successful fit has nonfinite spectral coordinates: {path}")
    if not (np.allclose(frame.loc[good, "k_rad_m"], frame.loc[good, "K_global"] / EARTH_RADIUS_M,
                        rtol=1e-10, atol=1e-14)
            and np.allclose(frame.loc[good, "omega_rad_s"],
                            2 * np.pi * frame.loc[good, "frequency_cpd"] / 86400,
                            rtol=1e-10, atol=1e-12)):
        raise RuntimeError(f"Dimensional coordinates disagree with K and f: {path}")
    return frame, path


def main():
    with TASKLIST.open(newline="") as source:
        simulations = list(csv.reader(source, delimiter="\t"))
    if len(simulations) != 59 or len({tuple(row) for row in simulations}) != 59:
        raise RuntimeError("Expected 59 unique CMIP6 model/experiment/member records")

    catalog_rows = [("NOAA_OLR", "observations", "not_applicable"), *simulations]
    frames, index = [], []
    for model, experiment, member in catalog_rows:
        frame, path = full_catalog(model, experiment, member)
        frame.insert(0, "model", model)
        frame.insert(1, "experiment", experiment)
        frame.insert(2, "member", member)
        frames.append(frame)
        index.append({
            "model": model, "experiment": experiment, "member": member,
            "catalog": str(path.relative_to(DATA_DIR)), "candidates": len(frame),
            "successful_fits": int(frame.fit_success.sum()),
            "quality_screened": int(frame.keep_for_histogram.sum()),
        })

    OUTPUT.mkdir(parents=True, exist_ok=True)
    events_path = OUTPUT / "kelvin_event_fits.csv.gz"
    index_path = OUTPUT / "kelvin_catalog_index.csv"
    if events_path.exists() or index_path.exists():
        raise FileExistsError(f"Preserving existing extraction in {OUTPUT}")
    pd.concat(frames, ignore_index=True).to_csv(events_path, index=False, compression="gzip")
    pd.DataFrame(index).to_csv(index_path, index=False)
    print(f"Extracted {sum(len(frame) for frame in frames)} events from {len(index)} full catalogs into {OUTPUT}")


if __name__ == "__main__":
    main()
