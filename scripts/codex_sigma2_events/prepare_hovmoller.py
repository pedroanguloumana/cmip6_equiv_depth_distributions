#!/usr/bin/env python3
"""Cache one NOAA Kelvin event for the threshold-comparison Hovmoller cells.

Run on Levante from the project root with the project environment active:
    python scripts/codex_sigma2_events/prepare_hovmoller.py

Replays the full 1990--2010 filter/detector, verifies the saved 2-sigma catalog,
and writes a small context slice to data/event_threshold_comparison/. It reuses
the saved carrier parameters; it neither refits the catalog nor writes figures.
"""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import hilbert
import xarray as xr

from src.paths import DATA_DIR, NOAA_OLR_FILE, PROJECT_ROOT
from src import wave_events, wave_fitting


OUTPUT = DATA_DIR / "event_threshold_comparison/noaa_kelvin_sigma2_example.nc"
TARGET_DATE = pd.Timestamp("2000-07-01")
TIME_PADDING_DAYS = 14
LONGITUDE_HALF_WIDTH_DEG = 90


def main():
    if OUTPUT.exists():
        raise FileExistsError(f"Preserving existing Hovmoller example: {OUTPUT}")
    if not NOAA_OLR_FILE.parent.is_symlink() or not NOAA_OLR_FILE.is_file():
        raise FileNotFoundError("Configure data/noaa_olr before preparing this example.")
    index = pd.read_csv(OUTPUT.parent / "catalog_summary.csv")
    record = index.loc[index.model.eq("NOAA_OLR") & index.wave.eq("Kelvin")
                       & index.detection_sigma.eq(2)]
    if len(record) != 1:
        raise ValueError("Expected one audited full 2-sigma NOAA Kelvin catalog.")
    catalog = wave_events.load_event_catalog("NOAA_OLR", "observations", "Kelvin",
                                             run_id=record.iloc[0].run_id)
    metadata = catalog.attrs["catalog_metadata"]
    cfg = metadata["configuration"]
    if cfg["fit_support"] != "component" or cfg["detection"]["detection_sigma"] != 2:
        raise ValueError("The example requires the audited 2-sigma component recipe.")
    for name, module in (("wave_events", wave_events), ("wave_fitting", wave_fitting)):
        if hashlib.sha256(Path(module.__file__).read_bytes()).hexdigest() != cfg["code_sha256"][name]:
            raise ValueError(f"Source code changed since catalog creation: {name}")
    for item in cfg["input"]:
        stat = (PROJECT_ROOT / item["file"]).stat()
        if (stat.st_size, stat.st_mtime_ns) != (item["size_bytes"], item["mtime_ns"]):
            raise ValueError("The NOAA input differs from the catalog provenance.")

    # Date proximity selects an illustration without selecting on depth or R2.
    eligible = catalog.loc[catalog.keep_for_histogram].copy()
    eligible["distance_to_target"] = (pd.to_datetime(eligible.peak_time) - TARGET_DATE).abs()
    chosen = eligible.sort_values(["distance_to_target", "event_id"]).iloc[0]
    print(f"Selected event {chosen.event_id}: {chosen.peak_time}", flush=True)
    with xr.open_dataset(NOAA_OLR_FILE) as source:
        olr = source.olr.sel(lat=slice(*wave_events.LAT_RANGE),
                             time=slice(cfg["start"], cfg["end"])).load()
    field = wave_events.filter_olr_for_wave(olr, "Kelvin")
    del olr
    events, sigma, supports = wave_events.detect_events(field, detection_sigma=2,
                                                        return_support=True)
    np.testing.assert_allclose(sigma, metadata["reference_sigma_W_m2"], rtol=1e-10)
    np.testing.assert_array_equal([e["event_id"] for e in events], catalog.event_id)
    np.testing.assert_array_equal(np.array([e["peak_time"] for e in events], dtype="datetime64[ns]"),
                                  pd.to_datetime(catalog.peak_time).to_numpy())
    np.testing.assert_allclose([e["peak_lon_deg"] for e in events], catalog.peak_lon_deg)
    np.testing.assert_allclose([e["peak_amplitude_W_m2"] for e in events],
                               catalog.peak_amplitude_W_m2, rtol=1e-10)

    ti, xi = supports[int(chosen.event_id)]
    peak_index = field.get_index("time").get_loc(pd.Timestamp(chosen.peak_time))
    relative_lon = (field.lon.values - chosen.peak_lon_deg + 180) % 360 - 180
    order = np.argsort(relative_lon)
    order = order[np.abs(relative_lon[order]) <= LONGITUDE_HALF_WIDTH_DEG]
    if not np.isin(xi, order).all():
        raise ValueError("The requested longitude context would truncate the chosen component.")
    first, stop = int(ti.min()) - TIME_PADDING_DAYS, int(ti.max()) + TIME_PADDING_DAYS + 1
    if first < 0 or stop > field.sizes["time"]:
        raise ValueError("The requested time context extends beyond the filtered record.")
    support = np.zeros(field.shape, dtype=np.uint8)
    support[ti, xi] = 1
    envelope = np.abs(hilbert(field.values, axis=1))
    patch = field.isel(time=slice(first, stop), lon=order).assign_coords(
        lon=chosen.peak_lon_deg + relative_lon[order],
        lag_day=("time", np.arange(first, stop) - peak_index))
    mask = support[first:stop, order].astype(bool)
    phase = (chosen.K_global * np.deg2rad(relative_lon[order])[None, :]
             - 2 * np.pi * chosen.frequency_cpd * patch.lag_day.values[:, None]
             - chosen.phase_rad)
    model = chosen.offset_W_m2 + chosen.amplitude_W_m2 * np.cos(phase)
    observed = patch.values[mask]
    r2 = 1 - np.sum((observed - model[mask]) ** 2) / np.sum((observed - observed.mean()) ** 2)
    np.testing.assert_allclose(r2, chosen.weighted_r2, rtol=1e-9, atol=1e-9)
    depth = wave_fitting.resting_wave_depth("Kelvin", chosen.k_rad_m, chosen.omega_rad_s)
    np.testing.assert_allclose(depth["h_eff_m"], chosen.h_eff_m, rtol=1e-9)
    if np.count_nonzero(mask) != len(ti) or not (envelope[ti, xi] > 2 * sigma).all():
        raise ValueError("The cached mask does not preserve the full active component.")

    example = xr.Dataset(dict(
        filtered_olr=patch,
        envelope=patch.copy(data=envelope[first:stop, order]),
        component_mask=patch.copy(data=mask.astype(np.uint8)),
        component_fit_olr=patch.copy(data=np.where(mask, model, np.nan)),
    ))
    for name in ("filtered_olr", "envelope", "component_fit_olr"):
        example[name].attrs = dict(units="W m-2")
    example.component_mask.attrs = dict(units="1", description="All active pixels of the chosen 2-sigma component")
    example.lon.attrs = dict(units="degrees_east", description="Longitudes unwrapped around the event peak")
    example.lag_day.attrs = dict(units="days", description="Elapsed days from the envelope peak")
    numeric = ("event_id", "peak_lon_deg", "duration_days", "span_deg", "K_global", "frequency_cpd",
               "phase_speed_m_s", "c0_m_s", "h_eff_m", "amplitude_W_m2", "offset_W_m2", "phase_rad",
               "weighted_r2", "spatial_cycles", "temporal_cycles")
    example.attrs = {name: float(chosen[name]) for name in numeric}
    example.attrs.update(
        wave_type="Kelvin", dataset="NOAA_OLR", detection_sigma=2.0,
        peak_time=str(chosen.peak_time), reference_sigma_W_m2=float(sigma),
        record_start=cfg["start"], record_end=cfg["end"],
        selection=f"Existing quality-screened event nearest {TARGET_DATE.date()}; ties use event_id; no depth selection",
        catalog=str(Path(catalog.attrs["catalog_path"]).relative_to(DATA_DIR)),
        catalog_run_id=str(record.iloc[0].run_id), configuration=json.dumps(cfg, sort_keys=True),
        preparation_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        verification="All replayed event IDs, peaks, longitudes and amplitudes match the saved catalog; complete support, R2 and branch depth verified",
    )
    temporary = OUTPUT.with_suffix(".nc.tmp")
    example.to_netcdf(temporary)
    temporary.replace(OUTPUT)
    print(f"Wrote {OUTPUT}: {dict(example.sizes)}; h_eff={chosen.h_eff_m:.3f} m, R2={r2:.6f}", flush=True)


if __name__ == "__main__":
    main()
