#!/usr/bin/env python3
"""Select daily CMIP6 historical/AMIP rlut with cmip_catalog; write a job manifest.

From the installed project root:
    python scripts/prepare_cmip6_rlut_regridding.py --start-year 1990 --end-year 2010
No climate fields are read or regridded by this planning step.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import re

import pandas as pd
from cmip_catalog import Catalog, MatchResult, build_filelists, load_catalog

from src.paths import CMIP6_ARCHIVE_DIR, CMIP6_CATALOG_FILE, CMIP6_REGRID_JOB_DIR, PROJECT_ROOT
from src.regridding import noaa_target_grid


def match_source_grids(catalog, start_year, end_year):
    """Prefer native gn, then other published grids in natural label order.

    This retains the earlier workflow's preference for GFDL-CM4 gr1 over gr2,
    rather than allowing the catalog's generic descending-label tie-breaker to
    choose its coarser published grid. Full-window coverage is still required.
    """
    matches, diagnostics = [], []
    for (institution, model, experiment, member), _ in catalog.df.groupby(
            ["institution_id", "source_id", "experiment_id", "member_id"], dropna=False):
        part = catalog.search(institution_id=institution, source_id=model,
                              experiment_id=experiment, member_id=member)
        grids = sorted(part.unique("grid_label"), key=lambda label: (
            label != "gn", [int(token) if token.isdigit() else token for token in re.split(r"(\d+)", label)]))
        for grid in grids:
            result = part.search(grid_label=grid).match(
                "rlut", frequency="daily", prefer_tables=("day",),
                time_window=(str(start_year), str(end_year)), coverage="full")
            diagnostics.append(result.diagnostics)
            if result.matches:
                matches.extend(result.matches)
                break
    return MatchResult(tuple(matches), pd.concat(diagnostics, ignore_index=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--start-year", type=int, default=1990)
    parser.add_argument("--end-year", type=int, default=2010)
    parser.add_argument("--models", nargs="+", help="Optional exact source_id selection")
    parser.add_argument("--experiments", nargs="+", choices=("historical", "amip"), default=["historical", "amip"])
    parser.add_argument("--manifest", type=Path, help="New manifest path under data/cmip6_regrid_jobs")
    parser.add_argument("--metadata", type=Path, help="Reuse explicitly saved inspected rlut metadata instead of probing again")
    args = parser.parse_args()
    if not 1 <= args.start_year <= args.end_year <= 9999:
        parser.error("Require increasing calendar years.")
    if PROJECT_ROOT != Path(__file__).resolve().parents[1]:
        raise RuntimeError("Install this checkout with pip -e .; a different src package is active.")
    if not CMIP6_ARCHIVE_DIR.is_symlink() or not CMIP6_ARCHIVE_DIR.is_dir():
        raise FileNotFoundError("Configure data/cmip6 as a symlink to the CMIP6 archive.")
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    manifest = args.manifest or CMIP6_REGRID_JOB_DIR / f"{args.start_year}-{args.end_year}.{stamp}" / "manifest.json"
    manifest = manifest.absolute()
    if not manifest.is_relative_to(CMIP6_REGRID_JOB_DIR) or manifest.exists():
        raise ValueError("Choose a new manifest path under data/cmip6_regrid_jobs/.")
    filters = dict(variable_id="rlut", experiment_id=args.experiments, frequency="daily", source_id=args.models)
    if args.metadata is None:
        candidates = load_catalog("CMIP6", source=CMIP6_CATALOG_FILE, filters=filters)
        # Derive the catalog URI prefix from the configured archive link. No
        # machine-specific archive root is built into the analysis workflow.
        archive_prefix = str(CMIP6_ARCHIVE_DIR.resolve()) + "/"
        if not candidates.df.uri.str.startswith(archive_prefix).all():
            raise ValueError("Catalog URIs do not match the archive exposed by data/cmip6.")
        candidates = candidates.remap_paths({archive_prefix: str(CMIP6_ARCHIVE_DIR) + "/"})
        groups = candidates.df.groupby(["institution_id", "source_id", "experiment_id"], dropna=False, sort=True)
        inspected, issues, diagnostics, matches = [], [], [], []
        # Try members in numeric order, stopping at the first complete match.
        # This avoids reopening hundreds of unused ensemble-member headers.
        for number, ((institution, model, experiment), group) in enumerate(groups, 1):
            members = sorted(group.member_id.unique(), key=lambda member: tuple(map(int, re.findall(r"\d+", member))))
            for member in members:
                print(f"{number}/{len(groups)}: inspecting {model}/{experiment}/{member}", flush=True)
                part = candidates.search(institution_id=institution, source_id=model,
                                         experiment_id=experiment, member_id=member)
                part = part.inspect_representatives(max_files=2*len(part.df), samples_per_dataset=2)
                matched = match_source_grids(part, args.start_year, args.end_year)
                inspected.append(part.df)
                issues.append(part.issues)
                diagnostics.append(matched.diagnostics)
                if matched.matches:
                    matches.extend(matched.matches)
                    break
        checked = Catalog(pd.concat(inspected, ignore_index=True), project="CMIP6")
        inspection_issues = pd.concat(issues, ignore_index=True)
        result = MatchResult(tuple(matches), pd.concat(diagnostics, ignore_index=True))
    else:
        if not args.metadata.absolute().is_relative_to(CMIP6_REGRID_JOB_DIR):
            raise ValueError("Saved study metadata must be under data/cmip6_regrid_jobs/.")
        # Restore the study table with pandas: an all-null Arrow conflict-flag
        # column cannot accept the matcher's fillna(False) in this environment.
        checked = Catalog(pd.read_parquet(args.metadata), project="CMIP6").search(**filters)
        if not checked.df.uri.str.startswith(str(CMIP6_ARCHIVE_DIR) + "/").all():
            raise ValueError("Saved metadata paths do not use this checkout's data/cmip6 entry point.")
        inspection_issues = checked.issues
        result = match_source_grids(checked, args.start_year, args.end_year)
    # One complete member per model/experiment, ranking the numeric (r,i,p,f)
    # indices. This prefers r1i1p1f1 when feasible, without requiring f1 for models
    # that publish f2/f3 or mixing members/versions to patch incomplete records.
    choices = {}
    for item in build_filelists(result):
        simulation = item.metadata["simulation"]
        key = tuple(simulation[name] for name in ("project", "institution_id", "source_id", "experiment_id"))
        rank = tuple(map(int, re.findall(r"\d+", simulation["member_id"])))
        if key not in choices or rank < choices[key][0]:
            choices[key] = (rank, item)
    if not choices:
        raise RuntimeError("No simulation covers all requested years; inspect the candidate metadata.")
    grid = noaa_target_grid()
    grid_hash = hashlib.sha256(b"".join(
        grid[name].values.astype("float64").tobytes()
        for name in ("lat", "lon", "lat_bnds", "lon_bnds")
    )).hexdigest()
    tasks = []
    for _, item in sorted(choices.values(), key=lambda pair: (
            pair[1].metadata["simulation"]["experiment_id"], pair[1].metadata["simulation"]["source_id"])):
        metadata = item.metadata
        files = []
        for filename in item.paths:
            path = Path(filename)
            stat = path.stat()
            files.append(dict(path=str(path.relative_to(PROJECT_ROOT)), size=stat.st_size, mtime_ns=stat.st_mtime_ns))
        tasks.append(dict(simulation=metadata["simulation"], table_id=metadata["table_id"],
                          source_grid_label=metadata["grid_label"], version=metadata["version"],
                          calendar=metadata["calendar"], year_range=[args.start_year, args.end_year],
                          target_grid_sha256=grid_hash, files=files, selection=metadata))
    manifest.parent.mkdir(parents=True, exist_ok=True)
    checked.df.to_parquet(manifest.parent / "inspected_metadata.parquet", index=False)
    inspection_issues.to_csv(manifest.parent / "inspection_issues.csv", index=False)
    result.diagnostics.to_csv(manifest.parent / "match_diagnostics.csv", index=False)
    report = dict(schema_version=1, created_utc=datetime.now(timezone.utc).isoformat(),
                  cmip_catalog_version=importlib.metadata.version("cmip-catalog-levante"),
                  member_selection="Lowest numeric r,i,p,f among complete simulations; one per model/experiment",
                  source_grid_selection="Native gn first, then other labels in ascending natural order; newest feasible version within each grid",
                  experiment_selection="All qualifying models independently; no cross-experiment intersection imposed",
                  catalog=str(CMIP6_CATALOG_FILE.relative_to(PROJECT_ROOT)),
                  catalog_mtime_ns=CMIP6_CATALOG_FILE.stat().st_mtime_ns,
                  target_grid="NOAA CDR OLR, global 360 x 180, exact centers and bounds",
                  tasks=tasks)
    manifest.write_text(json.dumps(report, indent=2) + "\n")
    for experiment in args.experiments:
        selected = [task for task in tasks if task["simulation"]["experiment_id"] == experiment]
        print(f"{experiment}: {len(selected)} model/member tasks", flush=True)
    print(f"Header inspection issues: {len(inspection_issues)}; see {manifest.parent}", flush=True)
    print(f"Manifest: {manifest}", flush=True)


if __name__ == "__main__":
    main()
