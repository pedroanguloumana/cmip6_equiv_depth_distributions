#!/usr/bin/env python3
"""Fit one NOAA OLR wave catalog; use submit_noaa_olr_events.sh for Slurm jobs.

Example from the project root, with this project installed:
    python scripts/noaa_olr_events.py --wave Kelvin --start 2000-01-01 --end 2010-12-31

Shared profiles, thresholds and fit options are explicit in src/wave_events.py.
"""

import argparse
from datetime import date
from pathlib import Path

import xarray as xr

from src.paths import NOAA_OLR_FILE, PROJECT_ROOT
from src.wave_events import LAT_RANGE, WAVE_PROFILES, filter_olr_for_wave, make_catalog_configuration, run_event_catalog


def load_noaa_record(start, end):
    """Load the requested daily NOAA record through the local archive symlink."""
    with xr.open_dataset(NOAA_OLR_FILE) as dataset:
        field = dataset["olr"].sel(lat=slice(*LAT_RANGE), time=slice(start, end))
        if not field.sizes["time"]:
            raise ValueError("The requested NOAA period contains no data.")
        if str(field.time.values[0])[:10] != start or str(field.time.values[-1])[:10] != end:
            raise ValueError("The NOAA file does not cover both requested record endpoints.")
        return field.load()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--wave", required=True, choices=tuple(WAVE_PROFILES))
    parser.add_argument("--start", default="1990-01-01", help="Inclusive first date, YYYY-MM-DD")
    parser.add_argument("--end", default="2010-12-31", help="Inclusive last date, YYYY-MM-DD")
    parser.add_argument("--detection-sigma", type=float, choices=(1.0, 1.5, 2.0), default=1.5,
                        help="Envelope activity cutoff in signed-field standard deviations")
    parser.add_argument("--fit-support", choices=("window", "component"), default="window",
                        help="Fit the fixed weighted window or active envelope component pixels")
    parser.add_argument("--max-events", type=int, help="Fit only the first N candidates; filename includes pilotN")
    parser.add_argument("--overwrite", action="store_true", help="Recompute this exact catalog and its checkpoints")
    args = parser.parse_args()
    if date.fromisoformat(args.start) > date.fromisoformat(args.end):
        parser.error("--start must not follow --end")
    if args.max_events is not None and args.max_events <= 0:
        parser.error("--max-events must be positive")
    if PROJECT_ROOT != Path(__file__).resolve().parents[1]:
        raise RuntimeError("The installed src package belongs to another checkout; install this project with pip -e .")
    if not NOAA_OLR_FILE.parent.is_symlink() or not NOAA_OLR_FILE.is_file():
        raise FileNotFoundError(f"Configure data/noaa_olr to expose {NOAA_OLR_FILE.name}")
    with xr.open_dataset(NOAA_OLR_FILE) as dataset:
        calendar = str(dataset.time.dt.calendar)
    identity = dict(dataset="NOAA_OLR", source_id="NOAA_CDR", experiment_id="observations",
                    member_id="not_applicable", source_grid_label="NOAA_1deg",
                    analysis_grid="NOAA_1deg", calendar=calendar)
    configuration = make_catalog_configuration(
        "NOAA_OLR", args.wave, args.start, args.end, identity=identity,
        input_paths=[NOAA_OLR_FILE], code_paths={"worker": Path(__file__)}, max_events=args.max_events,
        detection_sigma=args.detection_sigma, fit_support=args.fit_support,
    )
    run_event_catalog(
        configuration,
        lambda: filter_olr_for_wave(load_noaa_record(args.start, args.end), args.wave),
        overwrite=args.overwrite,
    )


if __name__ == "__main__":
    main()
