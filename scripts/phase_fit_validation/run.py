#!/usr/bin/env python3
"""Compare raw-value and phase-only carrier fits on known, modulated packets.

Run from the project root in the analysis environment. Inputs are completed
1.5-sigma WK99 calibration products under data/wave_significance/. One invocation
handles one source/experiment/wave; outputs stay beside this script as CSVs.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import hilbert
import xarray as xr
from xr_ccew.matsuno import profile_curve_frequency

from src.paths import PROJECT_ROOT, WAVE_SIGNIFICANCE_DIR
from src import wave_events as we, wave_significance as sg
from src.wave_fitting import fit_event_sinusoid
from src.wave_profiles import BOX_PROFILES


HERE = Path(__file__).resolve().parent
TRUTH_CASES = {
    "Kelvin": ((3, 25), (8, 8), (12, 8)),
    "n1ERW": ((-1, 25), (-3, 25), (-8, 8)),
    "n0EIG": ((1, 25), (3, 8), (8, 25)),
    "MRG": ((-1, 25), (-3, 25), (-8, 25)),
    "n1WIG": ((-1, 8), (-3, 8), (-8, 8)),
}
PEAK_RATIOS = (1.05, 1.5)


def calibration_root(source, experiment):
    parent = WAVE_SIGNIFICANCE_DIR / (
        "gaussian_calibration_wk99_121" if source == "NOAA_OLR"
        else "gaussian_calibration_wk99_121_cmip6"
    )
    matches = list(parent.glob(f"{source}.{experiment}.*.1990-2010.*/report.json"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one completed calibration for {source}/{experiment}; found {len(matches)}")
    root = matches[0].parent
    report = json.loads(matches[0].read_text())
    if (report["status"] != "complete" or report["configuration"]["background_method"] != "wk99_121"
            or report["configuration"]["detection"]["detection_sigma"] != 1.5):
        raise ValueError(f"Unexpected calibration configuration in {matches[0]}")
    return root


def read_threshold(root, wave):
    rows = [row for row in csv.DictReader((root / "thresholds.csv").open())
            if row["wave_type"] == wave and float(row["alpha"]) == 0.01]
    if len(rows) != 1:
        raise ValueError(f"Expected one 1% cutoff for {wave} in {root}")
    return float(rows[0]["threshold_W_m2"])


def open_array(path, use_cftime):
    coder = xr.coders.CFDatetimeCoder(use_cftime=use_cftime)
    with xr.open_dataarray(path, decode_times=coder) as data:
        return data.load()


def empirical_envelopes(observed, threshold, center_index, center_lon_index):
    """Return ordinary and rough observed amplitudes with a soft outer taper."""
    events, _, support = we.detect_events(observed, detection_sigma=1.5, return_support=True)
    analytic_envelope = abs(hilbert(observed.values, axis=1))
    candidates = [event for event in events if event["peak_amplitude_W_m2"] > threshold]
    selection = "significant"
    if len(candidates) < 2:
        candidates = sorted(events, key=lambda e: e["peak_amplitude_W_m2"], reverse=True)[:30]
        selection = "largest_candidates"
    if len(candidates) < 2:
        raise ValueError("Need at least two observed envelope components.")
    ranked = []
    for event in candidates:
        ti, xi = support[event["event_id"]]
        by_day = np.zeros(int(ti.max() - ti.min() + 1))
        np.maximum.at(by_day, ti - ti.min(), analytic_envelope[ti, xi])
        roughness = float(np.mean(abs(np.diff(by_day))) / max(np.mean(by_day), 1e-12))
        ranked.append((roughness, event, ti, xi))
    ranked.sort(key=lambda item: item[0])
    typical_index = len(ranked) // 2
    if typical_index == len(ranked) - 1:
        typical_index -= 1  # Two eligible events must yield two different shapes.
    choices = ("empirical_typical", ranked[typical_index]), ("empirical_rough", ranked[-1])
    result = {}
    for name, (roughness, event, ti, xi) in choices:
        peak_t = observed.get_index("time").get_loc(event["peak_time"])
        peak_x = int(np.argmin(abs(observed.lon.values - event["peak_lon_deg"])))
        # The active-pixel contour is a threshold boundary, not a physical
        # discontinuity in packet amplitude. Retain the actual envelope around
        # it and taper only outside the component's time/longitude extent.
        time_lag = np.arange(observed.sizes["time"]) - peak_t
        lon_lag = (observed.lon.values - event["peak_lon_deg"] + 180) % 360 - 180
        active_lon_lag = (observed.lon.values[xi] - event["peak_lon_deg"] + 180) % 360 - 180
        time_extent = max(abs(ti - peak_t).max() + max(5, .1 * event["duration_days"]), 5)
        lon_extent = max(abs(active_lon_lag).max() + max(10, .1 * event["span_deg"]), 10)
        taper_time = np.exp(-0.5 * (time_lag / time_extent) ** 4)
        taper_lon = np.exp(-0.5 * (lon_lag / lon_extent) ** 4)
        envelope = analytic_envelope * taper_time[:, None] * taper_lon[None, :]
        envelope = np.roll(envelope, (center_index - peak_t, center_lon_index - peak_x), axis=(0, 1))
        envelope /= envelope.max()
        result[name] = (envelope, dict(template_event_id=event["event_id"],
            template_selection=selection, template_duration_days=event["duration_days"],
            template_span_deg=event["span_deg"], template_roughness=roughness,
            template_peak_W_m2=event["peak_amplitude_W_m2"]))
    return result


def carrier_and_envelopes(template, wave, K, depth, center_index, center_lon_index, empirical):
    """Build positive envelopes around one exactly known carrier phase."""
    profile = BOX_PROFILES[wave]
    options = dict(meridional_mode_number=profile.meridional_mode_number)
    f = float(profile_curve_frequency(profile.curve_name, K, depth, **options))
    if not profile.frequency_min <= f <= profile.frequency_max:
        raise ValueError(f"Out-of-box truth: {wave} K={K}, depth={depth}, f={f}")
    step = 1e-3
    group_speed = 360 * (float(profile_curve_frequency(profile.curve_name, K + step, depth, **options))
                         - float(profile_curve_frequency(profile.curve_name, K - step, depth, **options))) / (2 * step)
    period = 1 / f
    lag = np.arange(template.sizes["time"]) - center_index
    center_lon = float(template.lon.values[center_lon_index])
    x = (template.lon.values - center_lon + 180) % 360 - 180
    moving_x = (x[None, :] - group_speed * lag[:, None] + 180) % 360 - 180
    spatial = np.exp(-0.5 * (moving_x / (0.5 * 360 / abs(K))) ** 2)
    temporal = np.exp(-0.5 * (lag / period) ** 2)
    envelopes = {
        "smooth": (temporal[:, None] * spatial, {}),
        "skewed": (temporal[:, None] * spatial * (1 + 0.7 * np.tanh(lag[:, None] / period)), {}),
        "bursty": (temporal[:, None] * spatial *
                   (0.35 + 0.65 * (0.5 + 0.5 * np.cos(2 * np.pi * lag[:, None] / (0.7 * period)))), {}),
        "double_lobe": ((np.exp(-0.5 * ((lag + period) / (0.65 * period)) ** 2)
                         + np.exp(-0.5 * ((lag - period) / (0.65 * period)) ** 2))[:, None] * spatial, {}),
        **empirical,
    }
    carrier = np.cos(K * np.deg2rad(x)[None, :] - 2 * np.pi * f * lag[:, None] - 0.4)
    for name, (envelope, metadata) in envelopes.items():
        envelope = envelope / envelope.max()
        packet = template.copy(data=envelope * carrier).rename("known_carrier_packet")
        yield name, packet, metadata, dict(K_global=K, frequency_cpd=f,
            period_days=period, h_eff_m=depth, center_index=center_index,
            center_lon_deg=center_lon, group_speed_deg_day=group_speed)


def phase_fit(field, event, support, k_grid, frequency_grid):
    """Fit a constant phase plane to the unit spatial analytic signal."""
    ti, xi = support
    peak_t = field.get_index("time").get_loc(event["peak_time"])
    lag = (ti - peak_t).astype(float)
    span = float(np.ptp(lag))
    if span <= 0:
        raise ValueError("Phase fit needs distinct time samples.")
    rel_lon = (field.lon.values[xi] - event["peak_lon_deg"] + 180) % 360 - 180
    x = np.deg2rad(rel_lon)
    unique_t, inverse = np.unique(ti, return_inverse=True)
    analytic_rows = hilbert(field.values[unique_t, :], axis=1)
    z = analytic_rows[inverse, xi]
    u = z / abs(z)
    direction = np.sign(k_grid[0])
    scaled_t = lag / span

    def model(K, cycles, phi):
        # Spatial Hilbert phases always have positive zonal wavenumber.
        return np.exp(1j * (direction * (K * x - 2 * np.pi * cycles * scaled_t) + phi))

    def residual(parameters):
        difference = u - model(*parameters)
        return np.r_[difference.real, difference.imag] / np.sqrt(2 * len(u))

    starts = []
    for K in k_grid:
        for f in frequency_grid:
            cycles = f * span
            base = model(K, cycles, 0)
            correlation = np.mean(u * np.conj(base))
            starts.append((abs(correlation), K, cycles, np.angle(correlation)))
    starts.sort(reverse=True)
    lower = np.array([k_grid[0], frequency_grid[0] * span, -np.inf])
    upper = np.array([k_grid[-1], frequency_grid[-1] * span, np.inf])
    solutions = [least_squares(residual, [K, cycles, phi], bounds=(lower, upper),
                               max_nfev=200, ftol=1e-9, xtol=1e-9, gtol=1e-9)
                 for _, K, cycles, phi in starts[:3]]
    solution = min(solutions, key=lambda fit: np.sum(fit.fun ** 2))
    K, cycles, phi = solution.x
    coherence = abs(np.mean(u * np.conj(model(K, cycles, 0))))
    return dict(K_global=float(K), frequency_cpd=float(cycles / span),
        fit_success=bool(solution.success), at_search_edge=bool(
            np.isclose(solution.x[:2], lower[:2], atol=1e-4, rtol=0).any()
            or np.isclose(solution.x[:2], upper[:2], atol=1e-4, rtol=0).any()),
        phase_coherence=float(coherence), phase_offset_rad=float(phi),
        spatial_cycles=float(abs(K) * np.ptp(rel_lon) / 360),
        temporal_cycles=float(cycles))


def matching_event(field, events, support, packet_mask):
    """Choose the strongest detected peak inside the known injected footprint."""
    matches = []
    time_index = field.get_index("time")
    longitude = np.asarray(field.lon.values)
    active_times = np.flatnonzero(packet_mask.any(axis=1))
    if not active_times.size:
        return None
    first, last = field.time.values[active_times[[0, -1]]]
    for event in events:
        if not first <= event["peak_time"] <= last:
            continue
        ti = time_index.get_loc(event["peak_time"])
        xi = int(np.argmin(abs(longitude - event["peak_lon_deg"])))
        if packet_mask[ti, xi]:
            matches.append(event)
    return max(matches, key=lambda item: item["peak_amplitude_W_m2"]) if matches else None


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"No rows for {path}")
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(dict.fromkeys(
            key for row in rows for key in row)))
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", choices=("NOAA_OLR", "ACCESS-CM2"), required=True)
    parser.add_argument("--experiment", choices=("observations", "historical", "amip"), required=True)
    parser.add_argument("--wave", choices=tuple(BOX_PROFILES), required=True)
    parser.add_argument("--backgrounds", type=int, default=6)
    parser.add_argument("--pilot", action="store_true", help="One truth, two shapes, one strength, one background")
    parser.add_argument("--er-k1-control", action="store_true",
                        help="ER K=-1 at 90 m depth (about 47 days), away from the 90-day box limit")
    args = parser.parse_args()
    if (args.source == "NOAA_OLR") != (args.experiment == "observations") or args.backgrounds < 1:
        parser.error("NOAA requires observations; CMIP requires historical/amip; backgrounds must be positive.")
    if args.er_k1_control and (args.wave != "n1ERW" or args.pilot):
        parser.error("--er-k1-control requires n1ERW without --pilot.")
    if PROJECT_ROOT != HERE.parents[1]:
        raise RuntimeError("Run this project's installed package from its own checkout.")
    root = calibration_root(args.source, args.experiment)
    threshold = read_threshold(root, args.wave)
    manifest = json.loads((root / "manifest.json").read_text())
    use_cftime = bool(manifest["use_cftime"])
    profile = BOX_PROFILES[args.wave]
    template = open_array(root / f"template.{profile.symmetry}.nc", use_cftime)
    observed = open_array(root / f"observed_filtered.{args.wave}.nc", use_cftime)
    with xr.open_dataset(root / f"spectrum.{profile.symmetry}.nc") as dataset:
        power = dataset.background_power.load()
    if template.shape != observed.shape or template.shape != power.shape:
        raise ValueError("Template, observed record and null spectrum must match.")
    suffix = f"{args.source}.{args.experiment}.{args.wave}" + (".pilot" if args.pilot else "")
    output_dir = HERE / "er_k1_control" if args.er_k1_control else HERE
    output_dir.mkdir(exist_ok=True)
    trial_path = output_dir / f"trials.{suffix}.csv"
    shape_path = output_dir / f"shapes.{suffix}.csv"
    if trial_path.exists() or shape_path.exists():
        raise FileExistsError(f"Preserving existing validation CSV for {suffix}")
    centers = (template.sizes["time"] // 2,)
    lon_centers = (180, 359, 60)
    empirical = empirical_envelopes(observed, threshold, centers[0], lon_centers[0])
    wave_index = tuple(BOX_PROFILES).index(args.wave)
    source_code = int.from_bytes(hashlib.sha256(f"{args.source}.{args.experiment}".encode()).digest()[:4], "little")
    backgrounds = []
    for replicate in range(1 if args.pilot else args.backgrounds):
        seed = [20261101, source_code, wave_index, replicate]
        random = np.random.default_rng(np.random.SeedSequence(seed))
        null = sg.draw_background(power, template, random)
        filtered = sg.filter_projected_olr(null, args.wave, profile=profile)
        backgrounds.append((replicate, filtered))
    reference_sigma = float(backgrounds[0][1].isel(time=slice(we.REFERENCE_EXCLUDE_DAYS,
        -we.REFERENCE_EXCLUDE_DAYS)).std())
    trials, shapes = [], []
    # Null maxima are diagnostic controls for fit scores, not a new significance
    # calibration. The six backgrounds per source cannot resolve a 1% false rate.
    for replicate, background in backgrounds:
        events, _, support = we.detect_events(background, detection_sigma=1.5, return_support=True)
        if not events:
            continue
        strongest = max(events, key=lambda event: event["peak_amplitude_W_m2"])
        for method in ("raw", "phase"):
            if method == "raw":
                fit, _ = fit_event_sinusoid(background, strongest["peak_time"],
                    strongest["peak_lon_deg"], k_grid=we.FIT_OPTIONS[args.wave]["k_grid"],
                    frequency_grid=we.FIT_OPTIONS[args.wave]["frequency_grid"],
                    support_indices=support[strongest["event_id"]])
            else:
                fit = phase_fit(background, strongest, support[strongest["event_id"]],
                    we.FIT_OPTIONS[args.wave]["k_grid"],
                    we.FIT_OPTIONS[args.wave]["frequency_grid"])
            trials.append(dict(source=args.source, experiment=args.experiment, wave=args.wave,
                calibration_root=str(root.relative_to(PROJECT_ROOT)), case="", morphology="",
                noise_type="null_maximum", background_replicate=replicate, peak_ratio_target="",
                threshold_W_m2=threshold, true_K_global="", true_frequency_cpd="",
                true_period_days="", true_h_eff_m="", template_event_id="", matched="",
                significant=strongest["peak_amplitude_W_m2"] > threshold,
                peak_amplitude_W_m2=strongest["peak_amplitude_W_m2"],
                support_pixels=len(support[strongest["event_id"]][0]),
                event_duration_days=strongest["duration_days"],
                event_span_deg=strongest["span_deg"], method=method,
                fit_success=fit["fit_success"], at_search_edge=fit["at_search_edge"],
                estimated_K_global=fit["K_global"],
                estimated_frequency_cpd=fit["frequency_cpd"], abs_K_error="",
                relative_frequency_error="", weighted_r2=fit.get("weighted_r2", ""),
                phase_coherence=fit.get("phase_coherence", ""), error=""))
    if args.er_k1_control:
        cases = ((-1, 90),)
    elif args.pilot:
        cases = TRUTH_CASES[args.wave][:1]
    else:
        cases = TRUTH_CASES[args.wave]
    for case_number, (K, depth) in enumerate(cases):
        center_lon_index = lon_centers[case_number]
        if case_number:
            empirical_case = empirical_envelopes(observed, threshold, centers[0], center_lon_index)
        else:
            empirical_case = empirical
        for morphology, packet, info, truth in carrier_and_envelopes(
                template, args.wave, K, depth, centers[0], center_lon_index, empirical_case):
            if args.pilot and morphology not in ("smooth", "empirical_rough"):
                continue
            filtered_packet = sg.filter_projected_olr(packet, args.wave, profile=profile)
            filtered_envelope = abs(hilbert(filtered_packet.values, axis=1))
            packet_peak = float(filtered_envelope.max())
            if packet_peak <= 0:
                raise ValueError("Filtered packet has no amplitude.")
            packet_mask = filtered_envelope >= 0.3 * packet_peak
            energy_fraction = float(np.sum(filtered_packet.values**2) / np.sum(packet.values**2))
            shapes.append(dict(source=args.source, experiment=args.experiment, wave=args.wave,
                case=case_number, morphology=morphology, true_K_global=K,
                true_frequency_cpd=truth["frequency_cpd"], true_period_days=truth["period_days"],
                center_lon_deg=truth["center_lon_deg"], retained_energy_fraction=energy_fraction,
                filtered_packet_peak_W_m2=packet_peak, **info))
            for ratio in PEAK_RATIOS[:1] if args.pilot else PEAK_RATIOS:
                scale = ratio * threshold / packet_peak
                clean = filtered_packet * scale
                clean_sigma = float(clean.isel(time=slice(we.REFERENCE_EXCLUDE_DAYS,
                    -we.REFERENCE_EXCLUDE_DAYS)).std())
                clean_events, _, clean_support = we.detect_events(clean,
                    detection_sigma=1.5 * reference_sigma / clean_sigma, return_support=True)
                all_fields = [("clean", -1, clean, clean_events, clean_support)]
                for replicate, background in backgrounds:
                    field = background + clean
                    events, _, support = we.detect_events(field, detection_sigma=1.5, return_support=True)
                    all_fields.append(("gaussian_background", replicate, field, events, support))
                for noise_type, replicate, field, events, support in all_fields:
                    matched = matching_event(field, events, support, packet_mask)
                    base = dict(source=args.source, experiment=args.experiment, wave=args.wave,
                        calibration_root=str(root.relative_to(PROJECT_ROOT)),
                        case=case_number, morphology=morphology, noise_type=noise_type,
                        background_replicate=replicate, peak_ratio_target=ratio,
                        threshold_W_m2=threshold, true_K_global=K,
                        true_frequency_cpd=truth["frequency_cpd"], true_period_days=truth["period_days"],
                        true_h_eff_m=depth, template_event_id=info.get("template_event_id", ""),
                        matched=matched is not None,
                        significant=bool(matched and matched["peak_amplitude_W_m2"] > threshold),
                        peak_amplitude_W_m2=matched["peak_amplitude_W_m2"] if matched else "",
                        support_pixels=len(support[matched["event_id"]][0]) if matched else "",
                        event_duration_days=matched["duration_days"] if matched else "",
                        event_span_deg=matched["span_deg"] if matched else "")
                    for method in ("raw", "phase"):
                        row = dict(base, method=method, fit_success="", at_search_edge="",
                            estimated_K_global="", estimated_frequency_cpd="", abs_K_error="",
                            relative_frequency_error="", weighted_r2="", phase_coherence="", error="")
                        if matched:
                            try:
                                indices = support[matched["event_id"]]
                                if method == "raw":
                                    fit, _ = fit_event_sinusoid(field, matched["peak_time"],
                                        matched["peak_lon_deg"], k_grid=we.FIT_OPTIONS[args.wave]["k_grid"],
                                        frequency_grid=we.FIT_OPTIONS[args.wave]["frequency_grid"],
                                        support_indices=indices)
                                    row["weighted_r2"] = fit["weighted_r2"]
                                else:
                                    fit = phase_fit(field, matched, indices,
                                        we.FIT_OPTIONS[args.wave]["k_grid"],
                                        we.FIT_OPTIONS[args.wave]["frequency_grid"])
                                    row["phase_coherence"] = fit["phase_coherence"]
                                row.update(fit_success=fit["fit_success"], at_search_edge=fit["at_search_edge"],
                                    estimated_K_global=fit["K_global"],
                                    estimated_frequency_cpd=fit["frequency_cpd"],
                                    abs_K_error=abs(fit["K_global"] - K),
                                    relative_frequency_error=abs(fit["frequency_cpd"] - truth["frequency_cpd"])
                                    / truth["frequency_cpd"])
                            except (ValueError, np.linalg.LinAlgError) as error:
                                row["error"] = f"{type(error).__name__}: {error}"
                        trials.append(row)
                print(args.source, args.experiment, args.wave, case_number, morphology,
                    f"target={ratio:g}x threshold", flush=True)
    write_csv(shape_path, shapes)
    write_csv(trial_path, trials)
    print("Wrote", shape_path, trial_path, flush=True)


if __name__ == "__main__":
    main()
