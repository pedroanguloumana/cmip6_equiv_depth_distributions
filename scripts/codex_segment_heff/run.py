#!/usr/bin/env python3
"""Retain 96-day WK segments and integrate excess power along depth bands.

Run on Levante from the project root with --task-index 0..59 (0 is NOAA).
See README.md in this directory for the scientific definitions and limitations.
"""

import argparse
import csv
import hashlib
from pathlib import Path

import numpy as np
import xarray as xr
import xr_ccew as tw
from xr_ccew.matsuno import profile_curve_frequency

from src.cmip import inspect_cmip_rlut, load_cmip_rlut
from src.paths import DATA_DIR, NOAA_OLR_FILE, PROJECT_ROOT
from src.wave_profiles import BOX_PROFILES


OUTPUT = DATA_DIR / "wave_significance/codex_segment_heff"
TASKLIST = DATA_DIR / "wave_significance/cmip6_calibration_jobs/simulations_1990-2010_wk99_ds1p5.tsv"
DEPTH_EDGES_M = np.arange(10.0, 80.1, 5.0)
SEGMENT_DAYS = 96
OVERLAP_DAYS = 60
LAT_BOUND = 15.0
WAVES = ("Kelvin", "n1ERW", "MRG", "n0EIG", "n1WIG")


def segment_spectra(olr, *, segment_days=96, overlap_days=60):
    """Return individual symmetric/antisymmetric WK powers on (segment, f, K).

    Input is daily OLR [W m-2] on a regular global (time, lat, lon) grid with
    mirrored latitudes. Average *power*, equally over latitude, after symmetry
    projection. Each two-sided bin is variance [W2 m-4], corrected for taper
    loss. Native datetime64 and cftime dates are retained as segment labels.
    """
    olr = olr.transpose("time", "lat", "lon").sortby("lat")
    if olr.attrs.get("units") != "W m-2" or not np.isfinite(olr.values).all():
        raise ValueError("Expected finite daily OLR in W m-2.")
    nt, nlat, nx = olr.shape
    if (nt < segment_days or not 0 <= overlap_days < segment_days
            or nlat < 2 or not np.allclose(olr.lat, -olr.lat.values[::-1], atol=1e-9)
            or not np.all(np.diff(olr.lat) > 0)
            or nx < 4 or not np.allclose(np.diff(olr.lon), 360 / nx)):
        raise ValueError("Require complete segments, mirrored latitudes and increasing global longitude.")
    times = olr.time.values
    steps = (np.diff(times) / np.timedelta64(1, "D")
             if np.issubdtype(times.dtype, np.datetime64)
             else np.array([(b - a).total_seconds() / 86400 for a, b in zip(times[:-1], times[1:])]))
    if not np.allclose(steps, 1, rtol=0, atol=1e-9):
        raise ValueError("Require regular daily sampling in the native calendar.")

    # Integer slices ensure exactly 96 daily samples, including the final full
    # window when present. xr_ccew.segment_data uses inclusive endpoints.
    starts = np.arange(0, nt - segment_days + 1, segment_days - overlap_days)
    ones = xr.DataArray(np.ones(segment_days), dims="time")
    taper = tw.apply_window(ones, window="tukey", pct=0.1).values
    norm = (segment_days * nx) ** 2 * np.mean(taper ** 2)
    frequency = -np.fft.fftfreq(segment_days)
    wavenumber = np.rint(np.fft.fftfreq(nx) * nx).astype(int)
    fi, ki = np.argsort(frequency), np.argsort(wavenumber)
    spectra = {}
    segment_variance = {}
    for symmetry in ("symmetric", "antisymmetric"):
        print(f"  preparing {symmetry} latitude fields", flush=True)
        component = tw.symmetric_antisymmetric_component(olr, symmetry)
        anomaly = tw.remove_harmonics_of_seasonal_cycle(
            tw.remove_mean_and_linear_trend(component), num_harmonics=3).load()
        powers, variances = [], []
        for start in starts:
            segment = tw.remove_mean_and_linear_trend(
                anomaly.isel(time=slice(start, start + segment_days))).values
            tapered = segment * taper[:, None, None]
            ft = np.fft.fft2(tapered, axes=(0, 2))
            power = (abs(ft) ** 2).mean(axis=1) / norm
            variance = np.mean(tapered ** 2) / np.mean(taper ** 2)
            if not np.isclose(power.sum(), variance, rtol=1e-10, atol=1e-10):
                raise RuntimeError("Segment power fails Parseval variance closure.")
            powers.append(power[np.ix_(fi, ki)])
            variances.append(variance)
        spectra[symmetry] = (("segment", "frequency", "zonal_wavenumber"), np.stack(powers))
        segment_variance[f"{symmetry}_variance"] = ("segment", variances)
        del anomaly, component
    result = xr.Dataset(
        {**spectra, **segment_variance},
        coords=dict(segment=np.arange(len(starts)), frequency=frequency[fi],
                    zonal_wavenumber=wavenumber[ki],
                    segment_start=("segment", times[starts]),
                    segment_end=("segment", times[starts + segment_days - 1])),
        attrs=dict(segment_days=segment_days, overlap_days=overlap_days,
                   latitude_aggregation="unweighted mean of latitude powers after symmetry projection",
                   latitude_min=float(olr.lat.min()), latitude_max=float(olr.lat.max()),
                   seasonal_harmonics=3, detrend="full record before harmonics, then each segment",
                   taper="Tukey pct=0.1", taper_mean_square=float(np.mean(taper ** 2)),
                   normalization="abs(FFT)^2/(Nt*Nlon)^2/mean(taper^2); two-sided variance per bin",
                   calendar=str(olr.time.dt.calendar), units="W2 m-4"),
    )
    for variable in result.data_vars.values():
        variable.attrs["units"] = "W2 m-4"
    return result


def depth_band_weights(frequency, wavenumber, profile, depth_edges_m):
    """Fraction of each frequency cell between successive depth curves at K.

    K remains discrete. Fractions assume uniform power within each Fourier
    frequency cell; they interpolate bin allocation, not spectral information.
    Exclude zero frequency, K=0, and the directionally ambiguous daily Nyquist.
    Return weights (depth, frequency, K) and one-bin depth scales (depth, K).
    """
    f, k = np.asarray(frequency), np.asarray(wavenumber)
    edges = np.asarray(depth_edges_m, dtype=float)
    if (len(f) < 2 or not np.allclose(np.diff(f), np.diff(f)[0])
            or np.diff(f)[0] <= 0 or len(edges) < 2
            or not np.isfinite(edges).all() or edges[0] <= 0 or not np.all(np.diff(edges) > 0)):
        raise ValueError("Require increasing regular frequencies and positive depth edges.")
    df = float(np.diff(f)[0])
    curves = np.array([profile_curve_frequency(
        profile.curve_name, k, h, meridional_mode_number=profile.meridional_mode_number,
        background_u=0.0) for h in edges])
    k_ok = (k >= profile.k_min) & (k <= profile.k_max) & (k != 0)
    # The five selected positive-frequency branches are monotone in depth.
    if np.any(np.diff(curves[:, k_ok], axis=0) < -1e-12):
        raise ValueError("Depth curves must increase in frequency on the selected branch.")
    low = np.maximum(curves[:-1], profile.frequency_min)
    high = np.minimum(curves[1:], profile.frequency_max)
    overlap = np.maximum(
        0.0, np.minimum(high[:, None, :], f[None, :, None] + df / 2)
        - np.maximum(low[:, None, :], f[None, :, None] - df / 2))
    weights = np.where(k_ok[None, None, :], overlap / df, 0.0)
    weights[:, (f <= 0) | (f >= 0.5 - 1e-12), :] = 0
    if not np.isfinite(weights).all() or np.any(weights.sum(axis=0) > 1 + 1e-10):
        raise RuntimeError("Depth bins overlap or contain invalid curve weights.")
    centers = (edges[:-1] + edges[1:]) / 2
    center_f = np.array([profile_curve_frequency(
        profile.curve_name, k, h, meridional_mode_number=profile.meridional_mode_number,
        background_u=0.0) for h in centers])
    scales = df * np.diff(edges)[:, None] / np.diff(curves, axis=0)
    visible = (k_ok[None, :] & (center_f >= profile.frequency_min)
               & (center_f <= profile.frequency_max) & (center_f < 0.5 - df / 2))
    return weights, np.where(visible, scales, np.nan)


def integrate_depth_bands(spectra, *, depth_edges_m=DEPTH_EDGES_M):
    """Integrate signed/positive excess against one dataset's fixed background.

    Background: mean symmetric/antisymmetric powers pooled over all segments,
    smoothed with ten 1-2-1 passes each in K and frequency, matching the current
    WK plotting recipe. Positive excess is descriptive, not a significance test.
    """
    mean = spectra[["symmetric", "antisymmetric"]].mean("segment")
    background = tw.background_spectrum((mean.symmetric + mean.antisymmetric) / 2,
                                       passes_wavenumber=10, passes_frequency=10)
    if not np.isfinite(background).all() or (background < 0).any():
        raise ValueError("Invalid smoothed background.")
    positive = spectra.sel(frequency=slice(1e-10, 0.5 - 1e-10),
                           zonal_wavenumber=slice(-15, 15))
    bg = background.sel(frequency=positive.frequency, zonal_wavenumber=positive.zonal_wavenumber).values
    if (bg <= 0).any():
        raise ValueError("Positive background power is required in the analysis domain.")
    f, k = positive.frequency.values, positive.zonal_wavenumber.values
    edges = np.asarray(depth_edges_m)
    output = {}
    weights_all, scales_all, box_all, box_bg_all = [], [], [], []
    fields = {key: [] for key in ("observed_power", "signed_excess", "positive_excess",
                                  "positive_ratio_excess")}
    bg_all, mean_excess_all, box_positive_all = [], [], []
    df = float(np.diff(f)[0])
    for wave in WAVES:
        profile = BOX_PROFILES[wave]
        power = positive[profile.symmetry].values
        weights, scales = depth_band_weights(f, k, profile, edges)
        # Reunite each positive-frequency cell with its real-field conjugate.
        fields["observed_power"].append(2 * np.einsum("sfk,hfk->sh", power, weights))
        fields["signed_excess"].append(2 * np.einsum("sfk,hfk->sh", power - bg, weights))
        fields["positive_excess"].append(2 * np.einsum("sfk,hfk->sh", np.maximum(power - bg, 0), weights))
        fields["positive_ratio_excess"].append(np.einsum(
            "sfk,hfk->sh", np.maximum(power / bg - 1, 0), weights))
        bg_all.append(2 * np.einsum("fk,hfk->h", bg, weights))
        mean_excess_all.append(2 * np.einsum(
            "fk,hfk->h", np.maximum(power.mean(axis=0) - bg, 0), weights))
        box = (np.maximum(0, np.minimum(f + df / 2, profile.frequency_max)
                          - np.maximum(f - df / 2, profile.frequency_min))[:, None] / df
               * ((k >= profile.k_min) & (k <= profile.k_max))[None, :])
        box_all.append(2 * np.einsum("sfk,fk->s", power, box))
        box_positive_all.append(2 * np.einsum("sfk,fk->s", np.maximum(power - bg, 0), box))
        box_bg_all.append(2 * np.sum(bg * box))
        weights_all.append(weights)
        scales_all.append(scales)
    for key, values in fields.items():
        output[key] = (("wave", "segment", "depth"), np.stack(values))
    output.update(
        background_power=(("wave", "depth"), bg_all),
        mean_spectrum_positive_excess=(("wave", "depth"), mean_excess_all),
        box_power=(("wave", "segment"), box_all),
        box_positive_excess=(("wave", "segment"), box_positive_all),
        box_background_power=("wave", box_bg_all),
        frequency_cell_weight=(("wave", "depth", "frequency", "zonal_wavenumber"), np.stack(weights_all)),
        one_frequency_bin_depth_scale=(("wave", "depth", "zonal_wavenumber"), np.stack(scales_all)),
        symmetric_mean_power=(("all_frequency", "all_wavenumber"), mean.symmetric.values),
        antisymmetric_mean_power=(("all_frequency", "all_wavenumber"), mean.antisymmetric.values),
        smoothed_background=(("all_frequency", "all_wavenumber"), background.values),
    )
    result = xr.Dataset(output, coords=dict(
        wave=list(WAVES), segment=spectra.segment, depth=(edges[:-1] + edges[1:]) / 2,
        depth_lower=("depth", edges[:-1]), depth_upper=("depth", edges[1:]),
        segment_start=spectra.segment_start, segment_end=spectra.segment_end,
        frequency=f, zonal_wavenumber=k, all_frequency=mean.frequency.values,
        all_wavenumber=mean.zonal_wavenumber.values), attrs=dict(spectra.attrs))
    result["effective_frequency_cells"] = result.frequency_cell_weight.sum(("frequency", "zonal_wavenumber"))
    for variable in result.data_vars.values():
        variable.attrs["units"] = "W2 m-4"
    for name in ("frequency_cell_weight", "effective_frequency_cells", "positive_ratio_excess"):
        result[name].attrs["units"] = "1"
    result.one_frequency_bin_depth_scale.attrs.update(
        units="m", long_name="one Fourier-bin frequency spacing divided by the depth-band mean df/dh; not an uncertainty estimate")
    result.depth.attrs.update(units="m", long_name="apparent resting-background equivalent depth")
    result.attrs.update(
        background="shared S/A pooled full-record mean, 10 K and 10 frequency 1-2-1 passes",
        background_u_m_s=0.0, depth_assignment="fractional frequency-cell overlap at discrete K; no added resolution",
        positive_excess_definition="sum weights*max(P_segment-B,0), with conjugate power factor 2",
        signed_excess_definition="sum weights*(P_segment-B), with conjugate power factor 2",
        ratio_excess_definition="sum weights*max(P_segment/B-1,0); dimensionless cell sum, not variance",
        nyquist="exclude exact daily Nyquist because eastward/westward direction is ambiguous",
        range_normalization="distribution shares are conditional on 10--80 m within each wave box",
        statistical_status="descriptive; positive clipping has noise bias; overlapping segments and adjacent depth bins are dependent",
        wk_reference="https://doi.org/10.1175/1520-0469(1999)056<0374:CCEWAO>2.0.CO;2",
    )
    np.testing.assert_allclose(result.observed_power - result.background_power,
                               result.signed_excess, rtol=1e-9, atol=1e-9)
    if (result.observed_power.sum("depth") > result.box_power * (1 + 1e-9) + 1e-9).any():
        raise RuntimeError("Depth-band power exceeds its containing wave box.")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-index", type=int, required=True)
    args = parser.parse_args()
    with TASKLIST.open(newline="") as handle:
        records = list(csv.reader(handle, delimiter="\t"))
    if len(records) != 59 or not 0 <= args.task_index <= len(records):
        raise ValueError("Expected task index 0 (NOAA) or 1..59 (CMIP6).")
    if args.task_index == 0:
        model, experiment, member = "NOAA_OLR", "observations", "not_applicable"
        path = NOAA_OLR_FILE
    else:
        model, experiment, member = records[args.task_index - 1]
        path = inspect_cmip_rlut(model, experiment=experiment)["path"]
    destination = OUTPUT / f"{model}.{experiment}.{member}.nc"
    if destination.exists():
        raise FileExistsError(f"Preserving existing segment analysis: {destination}")
    print(f"Analyzing {model}/{experiment}/{member}", flush=True)
    if args.task_index == 0:
        if not path.parent.is_symlink() or not path.is_file():
            raise FileNotFoundError("Configure the data/noaa_olr symlink.")
        with xr.open_dataset(path) as source:
            olr = source.olr.sel(time=slice("1990", "2010"), lat=slice(-LAT_BOUND, LAT_BOUND)).load()
        if str(olr.time.values[0])[:10] != "1990-01-01" or str(olr.time.values[-1])[:10] != "2010-12-31":
            raise ValueError("NOAA input must cover complete 1990--2010 calendar years.")
    else:
        olr = load_cmip_rlut(model, experiment=experiment, member_id=member,
                             year_range=(1990, 2010), lat_range=(-LAT_BOUND, LAT_BOUND))
    result = integrate_depth_bands(segment_spectra(
        olr, segment_days=SEGMENT_DAYS, overlap_days=OVERLAP_DAYS))
    result.attrs.update(model=model, experiment=experiment, member=member,
                        time_range="1990--2010 inclusive, native calendar",
                        source_path=str(path.relative_to(PROJECT_ROOT)),
                        source_size_bytes=path.stat().st_size,
                        source_mtime_ns=str(path.stat().st_mtime_ns),
                        analysis_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                        profiles_sha256=hashlib.sha256((PROJECT_ROOT / "src/wave_profiles.py").read_bytes()).hexdigest(),
                        numpy_version=np.__version__, xarray_version=xr.__version__)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(".tmp.nc")
    result.to_netcdf(temporary, engine="netcdf4", encoding={
        name: dict(zlib=True, complevel=2) for name in result.data_vars})
    temporary.replace(destination)
    print(f"Wrote {result.sizes['segment']} segments, five waves, 14 depth bins: {destination}", flush=True)


if __name__ == "__main__":
    main()
