#!/usr/bin/env python3
"""Bin the saved 96-day Kelvin segment spectrum and smooth WK-style background.

Run on Levante from the project root after extract.py:
    python scripts/codex_kelvin_komega/extract_wk.py
"""

import numpy as np
import pandas as pd
import xarray as xr

from variance import F_EDGES, K_EDGES, ROOT, calibration_run


def main():
    path = ROOT / "wk_segment_bins_integer_centered.csv.gz"
    if path.exists():
        raise FileExistsError(f"Preserving existing WK-style bins: {path}")
    index = pd.read_csv(ROOT / "kelvin_catalog_index.csv")
    if len(index) != 60:
        raise RuntimeError("Expected one NOAA and 59 CMIP6 records")
    rows = []
    for record in index.itertuples(index=False):
        run = calibration_run(record)
        with xr.open_dataset(run / "wk99_background.nc") as source:
            if source.attrs.get("background_method") != "wk99_121" or source.attrs.get("segment_days") != 96:
                raise RuntimeError(f"Unexpected spectral background method: {run}")
            observed = source.segment_power_symmetric.sel(
                frequency=slice(0.07, 0.41), zonal_wavenumber=slice(0.5, 15.5)).load()
            background = source.smoothed_positive_background.sel(
                frequency=slice(0.07, 0.41), zonal_wavenumber=slice(0.5, 15.5)).load()
        if (observed.shape != background.shape or observed.size == 0
                or not np.array_equal(observed.frequency, background.frequency)
                or not np.array_equal(observed.zonal_wavenumber, background.zonal_wavenumber)
                or not np.isfinite(observed).all() or not np.isfinite(background).all()
                or (observed < 0).any() or (background <= 0).any()):
            raise RuntimeError(f"Invalid segment spectrum or background: {run}")
        K, f = np.meshgrid(observed.zonal_wavenumber, observed.frequency)
        observed_bins = np.histogram2d(K.ravel(), f.ravel(), bins=(K_EDGES, F_EDGES),
                                       weights=observed.to_numpy().ravel())[0]
        background_bins = np.histogram2d(K.ravel(), f.ravel(), bins=(K_EDGES, F_EDGES),
                                         weights=background.to_numpy().ravel())[0]
        if (not np.isclose(observed_bins.sum(), observed.to_numpy().sum(), rtol=1e-10)
                or not np.isclose(background_bins.sum(), background.to_numpy().sum(), rtol=1e-10)
                or (background_bins <= 0).any()):
            raise RuntimeError(
                f"Segment cells lost in binning: {run}; observed {observed_bins.sum()}/"
                f"{observed.to_numpy().sum()}, background {background_bins.sum()}/"
                f"{background.to_numpy().sum()}, empty background bins "
                f"{int((background_bins <= 0).sum())}")
        for i, k_left in enumerate(K_EDGES[:-1]):
            for j, f_lower in enumerate(F_EDGES[:-1]):
                rows.append(dict(model=record.model, experiment=record.experiment,
                                 member=record.member, K_left=k_left, K_right=K_EDGES[i+1],
                                 f_lower_cpd=f_lower, f_upper_cpd=F_EDGES[j+1],
                                 segment_variance_W2_m4=observed_bins[i, j],
                                 smooth_background_variance_W2_m4=background_bins[i, j],
                                 segment_observed_to_background=(
                                     observed_bins[i, j] / background_bins[i, j])))
    temporary = path.with_name(path.name + ".tmp")
    pd.DataFrame(rows).to_csv(temporary, index=False, compression="gzip")
    temporary.replace(path)
    print(f"Wrote {len(rows)} bins for {len(index)} records to {path}", flush=True)


if __name__ == "__main__":
    main()
