#!/usr/bin/env python3
"""Bin all five 2-sigma event families in WK coordinates and cache WK overlays.

Run from the project root on Levante, with this checkout installed:
    python scripts/codex_wave_event_komega/extract.py

Reads the completed 1990--2010 event catalogs and saved 96-day spectra only.
Writes compact plotting inputs under data/event_threshold_comparison/; it
does not refilter raw OLR, refit events or simulate any background realizations.
"""

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr
from xr_ccew.matsuno import profile_curve_frequency

from src.paths import DATA_DIR, PROJECT_ROOT
from src.wave_events import load_event_catalog
from src.wave_fitting import EARTH_RADIUS_M, resting_wave_depth
from src.wave_profiles import BOX_PROFILES


ROOT = DATA_DIR / "event_threshold_comparison"
OUTPUT = ROOT / "wave_event_komega_sigma2"
SPECTRA = DATA_DIR / "wave_significance/codex_segment_heff"
WAVES = ("Kelvin", "MRG", "n0EIG", "n1ERW", "n1WIG")
SELECTIONS = ("fit_success", "keep_for_histogram")
K_CENTERS = np.arange(-15.0, 16.0)
K_EDGES = np.arange(-15.5, 16.0)
F_CENTERS = np.arange(49) / 96
F_EDGES = np.r_[0.0, (np.arange(48) + .5) / 96, .5]
DEPTHS_M = (12.0, 25.0, 50.0)


def main():
    paths = [OUTPUT / name for name in ("histograms_and_wk.nc", "catalog_index.csv", "audit.json")]
    if any(path.exists() for path in paths):
        raise FileExistsError(f"Preserving previous event/WK extraction: {OUTPUT}")
    audit = json.loads((ROOT / "audit.json").read_text())
    summary = pd.read_csv(ROOT / "catalog_summary.csv")
    summary = summary.loc[summary.detection_sigma.eq(2)].copy()
    keys = ["model", "experiment", "member"]
    records = summary[keys].drop_duplicates().reset_index(drop=True)
    if (audit["status"] != "complete" or len(summary) != 300 or len(records) != 60
            or summary.duplicated(keys + ["wave"]).any()):
        raise ValueError("Expected the audited 300-catalog, 60-record 2-sigma cohort.")

    counts = np.zeros((60, 5, 2, 49, 31), dtype=np.int32)
    selected = np.zeros((60, 5, 2), dtype=np.int32)
    nyquist = np.zeros_like(selected)
    power = np.zeros((60, 2, 47, 31))
    background = np.zeros((60, 47, 31))
    profiles, fit_bounds, index_rows, source_rows = {}, {}, [], []
    for r, record in records.iterrows():
        raw_input = None
        for w, wave in enumerate(WAVES):
            row = summary.loc[summary.model.eq(record.model) & summary.experiment.eq(record.experiment)
                              & summary.member.eq(record.member) & summary.wave.eq(wave)].squeeze()
            frame = load_event_catalog(record.model, record.experiment, wave,
                                       member_id=record.member, run_id=row.run_id)
            metadata = frame.attrs["catalog_metadata"]
            cfg = metadata["configuration"]
            if (cfg["detection"]["detection_sigma"] != 2 or cfg["fit_support"] != "component"
                    or len(frame) != row.candidates
                    or int(frame.fit_success.sum()) != row.successful_fits
                    or int(frame.keep_for_histogram.sum()) != row.screened_events):
                raise ValueError(f"Catalog differs from the completed 2-sigma audit: {row.catalog}")
            if frame[list(SELECTIONS)].isna().any().any():
                raise ValueError(f"Missing selection flags: {row.catalog}")
            good = frame.fit_success.to_numpy(dtype=bool)
            if (frame.keep_for_histogram.to_numpy(dtype=bool) & ~good).any():
                raise ValueError("The quality subset contains an unsuccessful fit.")
            coordinates = frame.loc[good, ["K_global", "frequency_cpd", "k_rad_m", "omega_rad_s"]]
            if not np.isfinite(coordinates).all().all():
                raise ValueError("A successful fit has nonfinite spectral coordinates.")
            np.testing.assert_allclose(coordinates.k_rad_m * EARTH_RADIUS_M,
                                       coordinates.K_global, rtol=1e-10, atol=1e-12)
            np.testing.assert_allclose(coordinates.omega_rad_s * 86400 / (2*np.pi),
                                       coordinates.frequency_cpd, rtol=1e-10, atol=1e-12)
            for s, selection in enumerate(SELECTIONS):
                keep = frame[selection].to_numpy(dtype=bool)
                k = frame.loc[keep, "K_global"].to_numpy()
                f = frame.loc[keep, "frequency_cpd"].to_numpy()
                if np.any(f <= 0) or np.any(f > .5 + 1e-12):
                    raise ValueError("An event frequency lies outside the daily resolved range.")
                # Protect only the rounding of the mathematical upper bound.
                f = np.minimum(f, .5)
                histogram = np.histogram2d(f, k, bins=(F_EDGES, K_EDGES))[0]
                if int(histogram.sum()) != int(keep.sum()):
                    raise ValueError(f"Event coordinates fall outside the displayed WK grid: {row.catalog}")
                counts[r, w, s] = histogram.astype(np.int32)
                selected[r, w, s] = int(keep.sum())
                nyquist[r, w, s] = np.count_nonzero(np.isclose(f, .5, rtol=0, atol=1e-8))
            profile = cfg["profile"]
            options = cfg["fit_options"]
            bounds = dict(K_min=options["k_grid"][0], K_max=options["k_grid"][-1],
                          f_min=options["frequency_grid"][0], f_max=options["frequency_grid"][-1])
            if r == 0:
                profiles[wave], fit_bounds[wave] = profile, bounds
            elif profile != profiles[wave] or bounds != fit_bounds[wave]:
                raise ValueError(f"Different wave profiles or fit search domains: {wave}")
            if raw_input is None:
                raw_input = cfg["input"][0]
            elif raw_input != cfg["input"][0]:
                raise ValueError("Wave families do not share the same source record.")
            path = Path(frame.attrs["catalog_path"])
            index_rows.append(dict(**record.to_dict(), wave=wave, detection_sigma=2.0,
                                   catalog=str(path.relative_to(DATA_DIR)), run_id=row.run_id,
                                   candidates=len(frame), successful_fits=int(good.sum()),
                                   screened_events=int(frame.keep_for_histogram.sum()),
                                   sha256=hashlib.sha256(path.read_bytes()).hexdigest()))

        spectral_path = SPECTRA / f"{record.model}.{record.experiment}.{record.member}.nc"
        with xr.open_dataset(spectral_path) as spectra:
            attrs = spectra.attrs
            if ((attrs["model"], attrs["experiment"], attrs["member"])
                    != (record.model, record.experiment, record.member)
                    or attrs["segment_days"] != 96 or attrs["overlap_days"] != 60
                    or attrs.get("time_range") != "1990--2010 inclusive, native calendar"
                    or attrs["source_path"] != raw_input["file"]
                    or int(attrs["source_size_bytes"]) != raw_input["size_bytes"]
                    or int(attrs["source_mtime_ns"]) != raw_input["mtime_ns"]):
                raise ValueError(f"The WK spectrum does not match the event source: {spectral_path}")
            part = spectra.sel(all_frequency=slice(1e-10, .5-1e-10),
                               all_wavenumber=slice(-15, 15))
            np.testing.assert_allclose(part.all_frequency, np.arange(1, 48)/96, atol=1e-12)
            np.testing.assert_array_equal(part.all_wavenumber, K_CENTERS)
            for s, symmetry in enumerate(("symmetric", "antisymmetric")):
                power[r, s] = part[f"{symmetry}_mean_power"].values
            background[r] = part.smoothed_background.values
            source_rows.append(dict(**record.to_dict(), source=str(spectral_path.relative_to(DATA_DIR)),
                                    latitude_aggregation=attrs["latitude_aggregation"],
                                    latitude_min=float(attrs["latitude_min"]), latitude_max=float(attrs["latitude_max"]),
                                    background=attrs["background"]))
        print(f"Read five event families and WK spectrum: {record.model}/{record.experiment}", flush=True)
    if (not np.isfinite(power).all() or (power < 0).any()
            or not np.isfinite(background).all() or (background <= 0).any()):
        raise ValueError("Invalid saved WK power or background.")

    curve_k = np.linspace(-15, 15, 601)
    curves = np.full((5, 3, len(curve_k)), np.nan)
    for w, wave in enumerate(WAVES):
        profile = BOX_PROFILES[wave]
        for h, depth in enumerate(DEPTHS_M):
            curves[w, h] = profile_curve_frequency(
                profile.curve_name, curve_k, depth,
                meridional_mode_number=profile.meridional_mode_number, background_u=0.0)
            test_k = 5.0 if profile.direction == "eastward" else -5.0
            test_f = profile_curve_frequency(profile.curve_name, np.array([test_k]), depth,
                       meridional_mode_number=profile.meridional_mode_number, background_u=0.0)[0]
            restored = resting_wave_depth(wave, test_k/EARTH_RADIUS_M, 2*np.pi*test_f/86400)
            np.testing.assert_allclose(restored["h_eff_m"], depth, rtol=1e-8, atol=1e-8)

    result = xr.Dataset(dict(
        event_count=(("record", "wave", "selection", "frequency", "K"), counts),
        selected_events=(("record", "wave", "selection"), selected),
        nyquist_events=(("record", "wave", "selection"), nyquist),
        wk_mean_power=(("record", "symmetry", "wk_frequency", "K"), power),
        wk_background=(("record", "wk_frequency", "K"), background),
        curve_frequency=(("wave", "equivalent_depth", "curve_K"), curves),
    ), coords=dict(
        record=[f"{row.model}.{row.experiment}.{row.member}" for row in records.itertuples()],
        model=("record", records.model.to_numpy()), experiment=("record", records.experiment.to_numpy()),
        member=("record", records.member.to_numpy()), wave=list(WAVES), selection=list(SELECTIONS),
        symmetry=["symmetric", "antisymmetric"], frequency=F_CENTERS, wk_frequency=np.arange(1, 48)/96,
        K=K_CENTERS, equivalent_depth=list(DEPTHS_M), curve_K=curve_k,
        K_lower=("K", K_EDGES[:-1]), K_upper=("K", K_EDGES[1:]),
        f_lower=("frequency", F_EDGES[:-1]), f_upper=("frequency", F_EDGES[1:]),
        wave_symmetry=("wave", [BOX_PROFILES[wave].symmetry for wave in WAVES]),
        **{f"fit_{name}": ("wave", [fit_bounds[wave][name] for wave in WAVES])
           for name in ("K_min", "K_max", "f_min", "f_max")},
    ), attrs=dict(
        status="complete", detection_sigma=2.0, period="1990--2010, native calendars",
        earth_radius_m=EARTH_RADIUS_M, frequency_units="cycles/day", K_definition="K = k_rad_m * Earth radius",
        omega_definition="omega_rad_s = 2*pi*frequency_cpd/86400", histogram="One event, one count; no amplitude weighting",
        primary_selection="All successful finite-coordinate fits; no depth or dispersion-validity selection",
        spectral_method="Saved 96-day / 60-day overlap WK-style latitude-power spectra, 15 S--15 N; shared S/A background with ten K and ten frequency 1-2-1 passes",
        event_latitude="Symmetry projection over 20 S--20 N then northern 0--20 latitude mean",
        nyquist="All events retained; last half-width bin touches daily Nyquist and is directionally ambiguous; WK overlays exclude exact Nyquist",
        worker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    ))
    result.frequency.attrs["units"] = "cycles/day"
    result.wk_frequency.attrs["units"] = "cycles/day"
    result.curve_frequency.attrs["units"] = "cycles/day"
    for name in ("wk_mean_power", "wk_background"):
        result[name].attrs["units"] = "W2 m-4"
    OUTPUT.mkdir(parents=True, exist_ok=True)
    temporary = paths[0].with_suffix(".tmp.nc")
    result.to_netcdf(temporary, engine="netcdf4", encoding={name: dict(zlib=True, complevel=4)
                                                           for name in result.data_vars})
    temporary.replace(paths[0])
    pd.DataFrame(index_rows).to_csv(paths[1], index=False)
    paths[2].write_text(json.dumps(dict(
        status="complete", detection_sigma=2.0, source_records=60, catalogs=300,
        successful_fits=int(selected[:, :, 0].sum()), screened_events=int(selected[:, :, 1].sum()),
        checks="Catalog identities/provenance; audited candidate and subset counts; k/omega unit conversions; histogram closure; matched WK input identities; branch-curve inversion",
        spectra=source_rows,
    ), indent=2) + "\n")
    print(f"Wrote {paths[0]}: 300 catalogs, two selections, 60 matched WK spectra", flush=True)


if __name__ == "__main__":
    main()
