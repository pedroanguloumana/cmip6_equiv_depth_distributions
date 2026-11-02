#!/usr/bin/env python3
"""Audit all 60 records and save compact model/depth and segment tables."""

import json

import numpy as np
import pandas as pd
import xarray as xr

from src.paths import DATA_DIR


ROOT = DATA_DIR / "wave_significance/codex_segment_heff"


def main():
    destination = ROOT / "depth_summary.csv.gz"
    segments_destination = ROOT / "segment_summary.csv.gz"
    if destination.exists() or segments_destination.exists():
        raise FileExistsError("Preserving existing segment-analysis summary tables.")
    paths = sorted(ROOT.glob("*.nc"))
    if len(paths) != 60:
        raise RuntimeError(f"Expected NOAA + 59 CMIP6 records, found {len(paths)}.")
    rows, segment_rows, records = [], [], []
    reference_weight = None
    reference_sha = None
    for path in paths:
        with xr.open_dataset(path) as opened:
            ds = opened.load()
        if (ds.attrs["segment_days"] != 96 or ds.attrs["overlap_days"] != 60
                or ds.sizes["wave"] != 5 or ds.sizes["depth"] != 14
                or not np.array_equal(ds.depth_lower, np.arange(10, 80, 5))):
            raise ValueError(f"Unexpected spectral/depth recipe: {path}")
        if reference_weight is None:
            reference_weight = ds.frequency_cell_weight.values
            reference_sha = ds.attrs["analysis_sha256"]
        elif (ds.attrs["analysis_sha256"] != reference_sha
              or not np.array_equal(reference_weight, ds.frequency_cell_weight.values)):
            raise ValueError(f"Analysis source or integration geometry differs: {path}")
        identity = {name: ds.attrs[name] for name in ("model", "experiment", "member", "calendar")}
        records.append(dict(**identity, n_segments=ds.sizes["segment"], file=path.name))
        for wave in ds.wave.values:
            band = ds.sel(wave=wave)
            mean_positive = band.positive_excess.mean("segment").values
            pooled_positive = band.mean_spectrum_positive_excess.values
            np.testing.assert_allclose(band.observed_power - band.background_power,
                                       band.signed_excess, rtol=1e-9, atol=1e-9)
            if (not np.isfinite(band.observed_power).all()
                    or (band.positive_excess < -1e-10).any()
                    or (band.observed_power.sum("depth") > band.box_power + 1e-7).any()
                    or (band.positive_excess.sum("depth") > band.box_positive_excess + 1e-7).any()):
                raise ValueError(f"Invalid powers or wave-box closure: {path}/{wave}")
            for i, depth in enumerate(band.depth.values):
                item = band.isel(depth=i)
                scales = item.one_frequency_bin_depth_scale.values
                scales = scales[np.isfinite(scales)]
                cells = float(item.effective_frequency_cells)
                bg = float(item.background_power)
                signed = float(item.signed_excess.mean("segment"))
                rows.append(dict(
                    **identity, wave=str(wave), depth_m=float(depth),
                    lower_m=float(item.depth_lower), upper_m=float(item.depth_upper),
                    n_segments=ds.sizes["segment"],
                    mean_observed_power=float(item.observed_power.mean("segment")),
                    background_power=bg, mean_signed_excess=signed,
                    mean_positive_excess=float(mean_positive[i]),
                    mean_spectrum_positive_excess=float(pooled_positive[i]),
                    positive_excess_fraction=(float(mean_positive[i] / mean_positive.sum())
                                              if mean_positive.sum() > 0 else np.nan),
                    mean_spectrum_excess_fraction=(float(pooled_positive[i] / pooled_positive.sum())
                                                   if pooled_positive.sum() > 0 else np.nan),
                    signed_excess_over_background=signed / bg if bg > 0 else np.nan,
                    mean_positive_ratio_excess=float(item.positive_ratio_excess.mean("segment")),
                    mean_positive_ratio_per_cell=(float(item.positive_ratio_excess.mean("segment")) / cells
                                                  if cells > 0 else np.nan),
                    effective_frequency_cells=cells,
                    median_one_bin_depth_scale_m=float(np.median(scales)) if len(scales) else np.nan,
                    box_mean_power=float(band.box_power.mean("segment")),
                    box_mean_positive_excess=float(band.box_positive_excess.mean("segment")),
                ))
            amounts = band.positive_excess.values
            for j in range(ds.sizes["segment"]):
                total = float(amounts[j].sum())
                segment_rows.append(dict(
                    **identity, wave=str(wave), segment=j,
                    start=str(ds.segment_start.values[j]), end=str(ds.segment_end.values[j]),
                    positive_excess_10_80_m=total,
                    signed_excess_10_80_m=float(band.signed_excess.isel(segment=j).sum()),
                    largest_mass_bin_m=(float(band.depth.values[np.argmax(amounts[j])]) if total > 0 else np.nan),
                    **{f"excess_{int(lo):02d}_{int(hi):02d}m": float(amounts[j, i])
                       for i, (lo, hi) in enumerate(zip(band.depth_lower.values, band.depth_upper.values))},
                ))
    identities = [(r["model"], r["experiment"], r["member"]) for r in records]
    if len(set(identities)) != 60:
        raise ValueError("Repeated source identities.")
    counts = pd.DataFrame(records).experiment.value_counts().to_dict()
    if counts != {"historical": 34, "amip": 25, "observations": 1}:
        raise ValueError(f"Unexpected cohort: {counts}")
    for path, values in ((destination, rows), (segments_destination, segment_rows)):
        temporary = path.with_name(path.name + ".tmp")
        pd.DataFrame(values).to_csv(temporary, index=False, compression="gzip")
        temporary.replace(path)
    report = dict(records=records, experiment_counts=counts, analysis_sha256=reference_sha,
                  depth_range_m=[10, 80], depth_bin_width_m=5,
                  checks="recipe and geometry match, finite powers, signed algebra, band/box closure, unique cohort",
                  caveat="No significance/equivalence test; segment overlap, finite frequency resolution and clipping bias remain.")
    (ROOT / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"Audited {len(records)} records; {len(rows)} depth rows and {len(segment_rows)} segment rows.")


if __name__ == "__main__":
    main()
