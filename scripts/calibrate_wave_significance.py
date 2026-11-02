#!/usr/bin/env python3
"""Calibrate candidate packet peaks using correlated Gaussian OLR backgrounds.

Run from the installed project, inside an allocation on Levante, for example:
  python scripts/calibrate_wave_significance.py --model NOAA_OLR --experiment observations --realizations 20
  python scripts/calibrate_wave_significance.py --model ACCESS-CM2 --experiment historical --realizations 20

Inputs enter through data/noaa_olr or data/cmip6_rlut_noaa1deg. Outputs under
data/wave_significance contain background spectra, null candidates, held-out
false-event counts, and observed candidate amplitudes. No event fits are changed.
"""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy
import xarray as xr

from src import wave_events as we, wave_significance as significance
from src.cmip import inspect_cmip_rlut, load_cmip_rlut
from src.paths import NOAA_OLR_FILE, PROJECT_ROOT, WAVE_SIGNIFICANCE_DIR


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--model", required=True, help="NOAA_OLR or exact CMIP6 source_id")
    parser.add_argument("--experiment", required=True, choices=("observations", "historical", "amip"))
    parser.add_argument("--member", help="Assert the CMIP6 member selected by the existing loader")
    parser.add_argument("--start-year", type=int, default=1990)
    parser.add_argument("--end-year", type=int, default=2010)
    parser.add_argument("--waves", nargs="+", choices=tuple(we.WAVE_PROFILES), default=list(we.WAVE_PROFILES))
    parser.add_argument("--realizations", type=int, default=20, help="Total; half calibrate, half validate. Four is a smoke test only.")
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--quantile", type=float, default=0.99, help="Pooled null-candidate peak quantile, not pixel/FDR significance")
    parser.add_argument("--frequency-smoothing-cpd", type=float, default=1/96)
    parser.add_argument("--wavenumber-smoothing", type=float, default=2.0)
    args = parser.parse_args()
    if PROJECT_ROOT != Path(__file__).resolve().parents[1]:
        raise RuntimeError("Install this checkout with pip -e . before running.")
    if (args.start_year > args.end_year or args.realizations < 4 or args.seed < 0
            or not 0 < args.quantile < 1):
        parser.error("Require ordered years, at least four realizations, nonnegative seed, and 0 < quantile < 1.")
    if args.model == "NOAA_OLR":
        if args.experiment != "observations" or args.member is not None:
            parser.error("NOAA_OLR requires observations and no ensemble member.")
        if not NOAA_OLR_FILE.parent.is_symlink() or not NOAA_OLR_FILE.is_file():
            raise FileNotFoundError("Configure the data/noaa_olr archive symlink.")
        path = NOAA_OLR_FILE
        with xr.open_dataset(path, chunks={"time": 365}) as dataset:
            olr = dataset.olr.sel(lat=slice(*we.LAT_RANGE), time=slice(str(args.start_year), str(args.end_year)))
            if (not olr.sizes["time"] or str(olr.time.values[0])[:10] != f"{args.start_year}-01-01"
                    or str(olr.time.values[-1])[:10] != f"{args.end_year}-12-31"):
                raise ValueError("NOAA record does not cover the requested full calendar years.")
            # Average lazily before loading: the calibration only needs two 2-D
            # projections, not a global stack of simulated latitude slices.
            components = significance.project_olr(olr)
        identity = dict(model="NOAA_OLR", experiment="observations", member="not_applicable")
    else:
        if args.experiment == "observations":
            parser.error("CMIP6 requires historical or amip.")
        info = inspect_cmip_rlut(args.model, experiment=args.experiment)
        path = info["path"]
        olr = load_cmip_rlut(args.model, experiment=args.experiment,
                            year_range=(args.start_year, args.end_year), member_id=args.member)
        components = significance.project_olr(olr)
        del olr
        identity = dict(model=args.model, experiment=args.experiment, member=info["member_id"])
    stat = path.stat()
    code_paths = dict(calibration=Path(significance.__file__), wave_events=Path(we.__file__), worker=Path(__file__))
    configuration = dict(identity=identity, start=str(components.time.values[0]), end=str(components.time.values[-1]),
        realizations=args.realizations, seed=args.seed, quantile=args.quantile,
        frequency_smoothing_cpd=args.frequency_smoothing_cpd, wavenumber_smoothing=args.wavenumber_smoothing,
        source=dict(path=str(path.relative_to(PROJECT_ROOT)), size_bytes=stat.st_size, mtime_ns=stat.st_mtime_ns),
        projection=dict(lat_range=we.LAT_RANGE, average_lat_range=we.AVERAGE_LAT_RANGE, latitude_weights="unweighted"),
        preprocessing=dict(seasonal_harmonics=we.SEASONAL_HARMONICS, remove_mean_and_trend=True, taper_fraction=we.TAPER_FRACTION),
        detection=dict(detection_sigma=we.DETECTION_SIGMA, peak_sigma=we.PEAK_SIGMA,
            min_duration_days=we.MIN_DURATION_DAYS, min_longitude_span_deg=we.MIN_LONGITUDE_SPAN_DEG,
            reference_exclude_days=we.REFERENCE_EXCLUDE_DAYS, periodic_longitude=True, connectivity=8),
        code_sha256={key: hashlib.sha256(value.read_bytes()).hexdigest() for key,value in code_paths.items()},
        versions=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, xarray=xr.__version__))
    WAVE_SIGNIFICANCE_DIR.mkdir(parents=True, exist_ok=True)
    for wave in args.waves:
        cfg = dict(configuration, wave_type=wave, profile=we.WAVE_PROFILES[wave].as_dict())
        run_id = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()[:12]
        name = f"{args.model}.{args.experiment}.{identity['member']}.{wave}.{args.start_year}-{args.end_year}.{run_id}"
        prefix = WAVE_SIGNIFICANCE_DIR/name
        report_path = prefix.with_name(prefix.name+".json")
        if report_path.exists():
            raise FileExistsError(f"Calibration already exists; preserve it and change settings for a new run: {report_path}")
        we.write_json(report_path, dict(status="running", configuration=cfg))
        started = time.monotonic()
        result = significance.calibrate_event_peaks(components[we.WAVE_PROFILES[wave].symmetry], wave,
            n_realizations=args.realizations, seed=args.seed, quantile=args.quantile,
            frequency_smoothing_cpd=args.frequency_smoothing_cpd, wavenumber_smoothing=args.wavenumber_smoothing)
        result["power"].to_netcdf(prefix.with_name(prefix.name+".background.nc"))
        for table in ("null_events", "realizations", "observed_events"):
            result[table].to_csv(prefix.with_name(prefix.name+f".{table}.csv"), index=False)
        report = dict(status="complete", configuration=cfg, **result["report"], elapsed_seconds=time.monotonic()-started)
        we.write_json(report_path, report)
        print(json.dumps(dict(wave=wave, report=str(report_path), threshold_W_m2=report["candidate_peak_threshold_W_m2"],
            heldout_false_candidates_per_record=report["heldout_false_candidates_per_record_mean"],
            observed_candidates=report["observed_candidates"], observed_above_threshold=report["observed_above_threshold"])), flush=True)


if __name__ == "__main__":
    main()
