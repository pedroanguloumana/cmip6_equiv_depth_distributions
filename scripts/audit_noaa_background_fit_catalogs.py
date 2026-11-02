#!/usr/bin/env python3
"""Verify fitted NOAA null catalogs against saved calibration event counts.

From the project root on Levante:
    python scripts/audit_noaa_background_fit_catalogs.py --manifest data/wave_significance/gaussian_calibration_wk99_121/NOAA_OLR.observations.not_applicable.1990-2010.b898b1ae4307/manifest.json

Writes an index and five-wave summary only if every expected catalog passes.
Use --start/--stop to audit a disjoint subset while the full array is running.
"""

import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path

from src.paths import DATA_DIR, PROJECT_ROOT
from src.wave_events import CSV_COLUMNS, WAVE_PROFILES


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--start", type=int)
    parser.add_argument("--stop", type=int)
    args = parser.parse_args()
    manifest_path = args.manifest.absolute()
    manifest_path.relative_to(DATA_DIR)
    manifest = json.loads(manifest_path.read_text())
    cfg = manifest["configuration"]
    first = cfg["n_calibration"]
    last = first + cfg["n_validation"]
    start = first if args.start is None else args.start
    stop = last if args.stop is None else args.stop
    if not first <= start < stop <= last:
        parser.error(f"Require held-out realization range within {first}:{last}")
    source = manifest_path.parent / "null_records.csv"
    saved = {(row["wave"], int(row["realization"])): row
             for row in csv.DictReader(source.open(newline=""))}
    catalog_dir = manifest_path.parent / "component_fit_catalogs"
    catalog_paths = [path for path in catalog_dir.glob("NOAA_OLR_GAUSSIAN_NULL.*.csv")
                     if ".pilot" not in path.name]
    index = []
    seen = set()
    for path in sorted(catalog_paths):
        sidecar = path.with_suffix(".json")
        if not sidecar.is_file():
            raise FileNotFoundError(f"Missing provenance: {sidecar}")
        metadata = json.loads(sidecar.read_text())
        configuration = metadata["configuration"]
        wave = configuration["wave_type"]
        realization = configuration["null_realization"]["index"]
        if not start <= realization < stop:
            continue
        key = (wave, realization)
        if key in seen:
            raise RuntimeError(f"Duplicate null catalog: {key}")
        seen.add(key)
        if (metadata["status"] != "complete" or wave not in WAVE_PROFILES
                or configuration["fit_support"] != "component"
                or configuration["detection"]["detection_sigma"] != 1.5
                or configuration["null_realization"]["ensemble_id"] != manifest["ensemble_id"]
                or configuration["identity"]["member_id"] != f"null_{realization:04d}"):
            raise RuntimeError(f"Unexpected catalog configuration or status: {path}")
        n = fits = retained = 0
        peak = 0.0
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames != CSV_COLUMNS:
                raise RuntimeError(f"Catalog CSV columns differ from observed NOAA: {path}")
            for row in reader:
                n += 1
                peak = max(peak, float(row["peak_amplitude_W_m2"]))
                fits += row["fit_success"] == "True"
                retained += row["keep_for_histogram"] == "True"
                if (row["wave_type"] != wave or row["member_id"] != f"null_{realization:04d}"
                        or (row["keep_for_histogram"] == "True"
                            and not (math.isfinite(float(row["h_eff_m"])) and float(row["h_eff_m"]) > 0))):
                    raise RuntimeError(f"Invalid event identity or retained depth: {path}")
        reference = saved[key]
        if (reference["role"] != "validation" or n != metadata["catalog_candidates"]
                or n != int(reference["candidates"])
                or fits != metadata["counts"]["fit_success"]
                or retained != metadata["counts"]["retained"]
                or not math.isclose(peak, float(reference["maximum_candidate_peak_W_m2"]),
                                    rel_tol=0, abs_tol=1e-8)):
            raise RuntimeError(f"Catalog does not match saved null realization: {path}")
        index.append(dict(wave_type=wave, realization=realization,
                          catalog=str(path.relative_to(PROJECT_ROOT)), candidates=n,
                          fit_success=fits, keep_for_histogram=retained,
                          maximum_candidate_peak_W_m2=peak,
                          elapsed_seconds=metadata["elapsed_seconds"],
                          peak_rss_gib=metadata["peak_rss_gib"]))
    expected = {(wave, realization) for wave in WAVE_PROFILES
                for realization in range(start, stop)}
    if seen != expected:
        missing = sorted(expected - seen)
        raise RuntimeError(f"Expected {len(expected)} catalogs; found {len(seen)}. First missing: {missing[:10]}")

    summary = defaultdict(lambda: dict(records=0, candidates=0, fit_success=0,
                                       keep_for_histogram=0))
    for row in index:
        totals = summary[row["wave_type"]]
        totals["records"] += 1
        for name in ("candidates", "fit_success", "keep_for_histogram"):
            totals[name] += row[name]
    suffix = f"{start}-{stop}"
    for name, rows in ((f"catalog_index.{suffix}.csv", index),
                       (f"catalog_summary_by_wave.{suffix}.csv",
                        [dict(wave_type=wave, **summary[wave]) for wave in WAVE_PROFILES])):
        output = catalog_dir / name
        temporary = output.with_name(output.name + ".tmp")
        with temporary.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(output)
        print(output, flush=True)
    print(f"Verified {len(index)} fitted null catalogs, {sum(row['candidates'] for row in index)} events", flush=True)


if __name__ == "__main__":
    main()
