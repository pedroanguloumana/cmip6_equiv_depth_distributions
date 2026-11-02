#!/usr/bin/env python3
"""Fit a CMIP6 model/wave event catalog, preserving its native calendar.

Examples from the project root with the project installed:
    python scripts/cmip6_rlut_events.py --list-models
    python scripts/cmip6_rlut_events.py --model CanESM5 --wave Kelvin --max-events 3

Input: data/cmip6_rlut_noaa1deg; output: DATA_DIR/cmip6_event_catalogs.
Full batches use submit_cmip6_rlut_events.sh, one task per model/experiment/wave.
"""

import argparse
from pathlib import Path

from src import cmip
from src.cmip import inspect_cmip_rlut, list_cmip_models, list_cmip_simulations, load_cmip_rlut
from src.paths import CMIP6_EVENT_DIR, PROJECT_ROOT
from src.wave_events import LAT_RANGE, WAVE_PROFILES, filter_olr_for_wave, make_catalog_configuration, run_event_catalog


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--list-models", action="store_true", help="List available model IDs and exit")
    parser.add_argument("--list-simulations", action="store_true", help="Print model/experiment/member TSV for batch discovery")
    parser.add_argument("--experiment", choices=("historical", "amip"), default="historical")
    parser.add_argument("--model", help="Exact CMIP6 source_id")
    parser.add_argument("--member", help="Assert the expected variant_label; otherwise use the unique file's member")
    parser.add_argument("--wave", choices=tuple(WAVE_PROFILES))
    parser.add_argument("--start-year", type=int, default=1990, help="Inclusive first complete native-calendar year")
    parser.add_argument("--end-year", type=int, default=2010, help="Inclusive last complete native-calendar year")
    parser.add_argument("--detection-sigma", type=float, choices=(1.0, 1.5, 2.0), default=1.5,
                        help="Envelope activity cutoff in signed-field standard deviations")
    parser.add_argument("--fit-support", choices=("window", "component"), default="window",
                        help="Fit the fixed weighted window or active envelope component pixels")
    parser.add_argument("--max-events", type=int, help="Pilot: fit only the first N candidates")
    parser.add_argument("--overwrite", action="store_true", help="Recompute this exact catalog and its checkpoints")
    args = parser.parse_args()
    if PROJECT_ROOT != Path(__file__).resolve().parents[1]:
        raise RuntimeError("The installed src package belongs to another checkout; install this project with pip -e .")
    if args.list_models:
        print("\n".join(list_cmip_models(args.experiment)))
        return
    if args.list_simulations:
        for item in list_cmip_simulations():
            print("\t".join(item[key] for key in ("source_id", "experiment_id", "member_id")))
        return
    if not args.model or not args.wave:
        parser.error("--model and --wave are required unless --list-models is used")
    if not 1 <= args.start_year <= args.end_year <= 9999:
        parser.error("Require 1 <= --start-year <= --end-year <= 9999")
    if args.max_events is not None and args.max_events <= 0:
        parser.error("--max-events must be positive")

    info = inspect_cmip_rlut(args.model, experiment=args.experiment)
    if args.member is not None and args.member != info["member_id"]:
        parser.error(f"{args.model} contains {info['member_id']}, not {args.member}")
    start = f"{args.start_year:04d}-01-01"
    # The analysis ends on the last actual day of the selected calendar year.
    end_day = 30 if info["calendar"] == "360_day" else 31
    end = f"{args.end_year:04d}-12-{end_day}"
    identity = {key: info[key] for key in (
        "source_id", "experiment_id", "member_id", "source_grid_label", "analysis_grid", "calendar",
    )}
    identity["dataset"] = "CMIP6"
    dataset_id = f"CMIP6.{info['source_id']}.{info['experiment_id']}.{info['member_id']}.rlut"
    provenance = info["path"].with_suffix(".json")
    if not provenance.is_file():
        raise FileNotFoundError(f"Missing regridding provenance: {provenance}")
    configuration = make_catalog_configuration(
        dataset_id, args.wave, start, end, identity=identity, input_paths=[info["path"], provenance],
        code_paths={"worker": Path(__file__), "cmip_loader": Path(cmip.__file__)},
        max_events=args.max_events, detection_sigma=args.detection_sigma,
        fit_support=args.fit_support,
    )
    print(f"Model={args.model} experiment={args.experiment} member={info['member_id']} calendar={info['calendar']} "
          f"source_grid={info['source_grid_label']} analysis_grid={info['analysis_grid']}", flush=True)
    run_event_catalog(
        configuration,
        lambda: filter_olr_for_wave(load_cmip_rlut(
            args.model, experiment=args.experiment, year_range=(args.start_year, args.end_year),
            lat_range=LAT_RANGE, member_id=info["member_id"],
        ), args.wave),
        output_dir=CMIP6_EVENT_DIR, overwrite=args.overwrite,
    )


if __name__ == "__main__":
    main()
