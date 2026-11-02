#!/usr/bin/env python3
"""Fit one saved NOAA Gaussian-null realization with the observed catalog pipeline.

Run from the project root on Levante, for example:
    python scripts/fit_noaa_background_events.py \
        --manifest data/wave_significance/gaussian_calibration_wk99_121/NOAA_OLR.observations.not_applicable.1990-2010.b898b1ae4307/manifest.json \
        --wave Kelvin --realization 1000

The worker regenerates a held-out Gaussian record from its saved spectrum and
seed. It writes the observed catalog's CSV columns and JSON provenance under
<ensemble>/component_fit_catalogs/. Completed CSVs retain no per-event parts;
an interrupted catalog retains its parts for a resumable retry.
"""

import argparse
import json
from pathlib import Path
import shutil

import numpy as np
import xarray as xr

import gaussian_event_calibration as calibration
from src import wave_significance as significance
from src.paths import DATA_DIR, PROJECT_ROOT
from src.wave_events import (WAVE_PROFILES, catalog_path,
                             make_catalog_configuration, run_event_catalog)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--wave", choices=tuple(WAVE_PROFILES), required=True)
    parser.add_argument("--realization", type=int, required=True)
    parser.add_argument("--max-events", type=int, help="Pilot: fit only the first N detected events")
    args = parser.parse_args()
    if PROJECT_ROOT != Path(__file__).resolve().parents[1]:
        raise RuntimeError("The installed src package belongs to another checkout; install this project with pip -e .")
    manifest_path = args.manifest.absolute()
    manifest_path.relative_to(DATA_DIR)
    manifest = json.loads(manifest_path.read_text())
    cfg = manifest["configuration"]
    if (manifest["status"] != "complete" or cfg["model"] != "NOAA_OLR"
            or cfg["background_method"] != "wk99_121"
            or cfg["detection"]["detection_sigma"] != 1.5):
        raise ValueError("Expected the completed 1.5-sigma NOAA WK99 Gaussian ensemble.")
    if not cfg["n_calibration"] <= args.realization < cfg["n_calibration"] + cfg["n_validation"]:
        raise ValueError("The realization must belong to the held-out validation ensemble.")
    if args.max_events is not None and args.max_events < 1:
        parser.error("--max-events must be positive")
    saved_profile = json.dumps(cfg["profiles"][args.wave], sort_keys=True)
    current_profile = json.dumps(WAVE_PROFILES[args.wave].as_dict(), sort_keys=True)
    if saved_profile != current_profile:
        raise RuntimeError("The saved wave profile differs from the observed catalog profile.")

    # Saved templates and spectra replace the old preparation code. The draw,
    # filter and spectral-profile recipe must still match the saved ensemble.
    current_hashes = calibration.simulation_hashes("wk99_121", 1.5)
    for name in ("draw", "filter", "recipe", "xr_filtering", "xr_spectral", "xr_matsuno"):
        if current_hashes[name] != cfg["simulation_sha256"][name]:
            raise RuntimeError(f"Saved null ensemble and current {name} code differ.")

    symmetry = cfg["profiles"][args.wave]["symmetry"]
    template_path = manifest_path.parent / f"template.{symmetry}.nc"
    spectrum_path = manifest_path.parent / f"spectrum.{symmetry}.nc"
    template = calibration.open_template(manifest_path.parent, manifest, args.wave)
    start, end = str(template.time.values[0])[:10], str(template.time.values[-1])[:10]
    if (start != f"{cfg['start_year']}-01-01" or end != f"{cfg['end_year']}-12-31"):
        raise RuntimeError("The saved null template covers different dates.")
    identity = dict(dataset="NOAA_OLR_GAUSSIAN_NULL", source_id="NOAA_CDR_WK99_NULL",
                    experiment_id="validation", member_id=f"null_{args.realization:04d}",
                    source_grid_label="NOAA_1deg", analysis_grid="NOAA_1deg",
                    calendar=str(template.time.dt.calendar))
    configuration = make_catalog_configuration(
        f"NOAA_OLR_GAUSSIAN_NULL.validation.{args.realization:04d}", args.wave,
        start, end, identity=identity,
        input_paths=[manifest_path, template_path, spectrum_path],
        code_paths=dict(worker=Path(__file__), significance=Path(significance.__file__),
                        calibration=Path(calibration.__file__)),
        max_events=args.max_events, detection_sigma=1.5, fit_support="component",
    )
    configuration["null_realization"] = dict(
        ensemble_id=manifest["ensemble_id"], role="validation", index=args.realization,
        seed_components=[cfg["seed"], 0 if symmetry == "symmetric" else 1, args.realization],
        background_method=cfg["background_method"], symmetry=symmetry,
    )
    output_dir = manifest_path.parent / "component_fit_catalogs"
    output = catalog_path(configuration, output_dir)

    def load_filtered_record():
        with xr.open_dataset(spectrum_path) as source:
            power = source.background_power.load()
        seed = np.random.SeedSequence(configuration["null_realization"]["seed_components"])
        background = significance.draw_background(power, template, np.random.default_rng(seed))
        return significance.filter_projected_olr(
            background, args.wave, profile=WAVE_PROFILES[args.wave])

    run_event_catalog(configuration, load_filtered_record, output_dir=output_dir)
    metadata = json.loads(output.with_suffix(".json").read_text())
    same_configuration = (json.dumps(metadata["configuration"], sort_keys=True)
                          == json.dumps(configuration, sort_keys=True))
    if metadata["status"] != "complete" or not same_configuration or not output.is_file():
        raise RuntimeError(f"Catalog did not finish cleanly: {output}")
    parts = output.with_name(output.name + ".parts")
    if parts.exists():
        shutil.rmtree(parts)


if __name__ == "__main__":
    main()
