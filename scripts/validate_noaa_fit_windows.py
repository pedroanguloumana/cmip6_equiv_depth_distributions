#!/usr/bin/env python3
"""Compare original and alternative NOAA fits on known-depth 8-sigma injections.

Run from the project root on Levante after the WK99-style NOAA calibration:
    python scripts/validate_noaa_fit_windows.py --candidate component

Reads the saved Gaussian spectrum, seeds, templates and injection cases. Writes
trial fits and a provenance summary to data/wave_significance/. No figures or
event catalogs are changed. Null controls fit the strongest candidate in each
of the eight reused backgrounds; they do not estimate a 1% false-fit rate.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import xarray as xr

from src.paths import WAVE_SIGNIFICANCE_DIR
from src import wave_events as we, wave_significance as sg
from src.wave_fitting import fit_event_sinusoid, resting_wave_depth
from src.wave_profiles import BOX_PROFILES


RUN = (WAVE_SIGNIFICANCE_DIR / "gaussian_calibration_wk99_121"
       / "NOAA_OLR.observations.not_applicable.1990-2010.b898b1ae4307")
WAVES = ("n1ERW", "n0EIG", "n1WIG")
SHORT_WINDOWS = {
    "n1ERW": dict(half_window_days=45, half_width_deg=60, time_sigma_days=22.5, lon_sigma_deg=30),
    "n0EIG": dict(half_window_days=5, half_width_deg=60, time_sigma_days=2.5, lon_sigma_deg=30),
    "n1WIG": dict(half_window_days=5, half_width_deg=60, time_sigma_days=2.5, lon_sigma_deg=30),
}


def matching_significant_event(events, field, truth, threshold):
    """Pick the largest above-threshold detected peak within one packet sigma."""
    if not events:
        return None
    indices = field.get_index("time").get_indexer([event["peak_time"] for event in events])
    lag = indices - truth["center_index"]
    lon = np.array([event["peak_lon_deg"] for event in events])
    moving_lon = (lon - truth["center_lon_deg"] - truth["group_speed_deg_day"] * lag + 180) % 360 - 180
    matches = [event for event, nearby in zip(events,
        (np.abs(lag) <= truth["time_sigma_days"])
        & (np.abs(moving_lon) <= truth["lon_sigma_deg"]))
        if nearby and event["peak_amplitude_W_m2"] > threshold]
    return max(matches, key=lambda event: event["peak_amplitude_W_m2"]) if matches else None


def fit_peak(field, event, wave, options, *, label, source, realization,
             case=None, truth=None, support_indices=None):
    if label == "component":
        options = dict(k_grid=options["k_grid"], frequency_grid=options["frequency_grid"],
                       support_indices=support_indices)
    fit, _ = fit_event_sinusoid(field, event["peak_time"], event["peak_lon_deg"], **options)
    fit.update(resting_wave_depth(wave, fit["k_rad_m"], fit["omega_rad_s"]))
    quality = (fit["fit_success"] and not fit["at_search_edge"]
               and fit["weighted_r2"] >= we.MIN_WEIGHTED_R2
               and fit["spatial_cycles"] >= we.MIN_SPATIAL_CYCLES
               and fit["temporal_cycles"] >= we.MIN_TEMPORAL_CYCLES
               and fit["dispersion_valid"])
    row = dict(wave_type=wave, source=source, window=label, realization=realization,
               case=case, event_id=event["event_id"], peak_time=str(event["peak_time"]),
               peak_lon_deg=event["peak_lon_deg"], peak_amplitude_W_m2=event["peak_amplitude_W_m2"],
               quality=bool(quality), **fit)
    if truth is not None:
        row.update(true_K_global=truth["K_global"], true_frequency_cpd=truth["frequency_cpd"],
                   true_h_eff_m=truth["h_eff_m"], width=truth["width"],
                   relative_depth_error=abs(fit["h_eff_m"] / truth["h_eff_m"] - 1),
                   absolute_K_error=abs(fit["K_global"] - truth["K_global"]))
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", choices=("short", "component"), default="short")
    args = parser.parse_args()
    output_name = ("noaa_fit_window_validation_b898b1ae4307" if args.candidate == "short"
                   else "noaa_component_fit_validation_b898b1ae4307")
    output = WAVE_SIGNIFICANCE_DIR / output_name
    if output.exists():
        raise FileExistsError(f"Preserving existing validation output: {output}")
    report = json.loads((RUN / "report.json").read_text())
    if report["status"] != "complete" or report["configuration"]["background_method"] != "wk99_121":
        raise RuntimeError("The expected completed WK99-style calibration is unavailable.")
    cfg = report["configuration"]
    injections = pd.read_csv(RUN / "injections.csv.gz")
    injections = injections.loc[(injections.amplitude_sigma == 8)
        & injections.carrier_inside_rectangular_bounds & injections.carrier_below_daily_nyquist
        & injections.wave_type.isin(WAVES)]
    rows = []
    counts = []
    for wave in WAVES:
        profile = BOX_PROFILES[wave]
        symmetry = profile.symmetry
        with xr.open_dataarray(RUN / f"template.{symmetry}.nc") as source:
            template = source.load()
        with xr.open_dataset(RUN / f"spectrum.{symmetry}.nc") as source:
            power = source.background_power.load()
        backgrounds = {}
        for realization in range(cfg["n_calibration"], cfg["n_calibration"] + cfg["injection_realizations"]):
            rng = np.random.default_rng(np.random.SeedSequence([
                cfg["seed"], 0 if symmetry == "symmetric" else 1, realization]))
            null = sg.draw_background(power, template, rng)
            filtered = sg.filter_projected_olr(null, wave, profile=profile)
            null_events, _, null_support = we.detect_events(filtered,
                detection_sigma=cfg["detection"]["detection_sigma"], return_support=True)
            backgrounds[realization] = filtered
            time_index = filtered.get_index("time")
            half_window = we.FIT_OPTIONS[wave]["half_window_days"]
            eligible = [event for event in null_events
                        if half_window <= time_index.get_loc(event["peak_time"])
                        < filtered.sizes["time"] - half_window]
            if eligible:
                strongest = max(eligible, key=lambda event: event["peak_amplitude_W_m2"])
                alternative = (dict(we.FIT_OPTIONS[wave], **SHORT_WINDOWS[wave])
                               if args.candidate == "short" else we.FIT_OPTIONS[wave])
                for label, options in (("original", we.FIT_OPTIONS[wave]),
                                       (args.candidate, alternative)):
                    rows.append(fit_peak(filtered, strongest, wave, options, label=label,
                                         source="null_record_maximum", realization=realization,
                                         support_indices=null_support[strongest["event_id"]]
                                         if label == "component" else None))
        wave_cases = injections.loc[injections.wave_type == wave]
        for case, group in wave_cases.groupby("case"):
            first = group.iloc[0]
            time_width, lon_width = cfg["injection"]["widths"][first.width]
            packet, truth = sg.synthetic_wave_packet(template, wave, first.K_global, first.h_eff_m,
                center_index=int(first.center_index), center_lon=float(first.center_lon_deg),
                time_sigma_periods=time_width, lon_sigma_wavelengths=lon_width, profile=profile)
            filtered_packet = sg.filter_projected_olr(packet, wave, profile=profile)
            truth["width"] = first.width
            for trial in group.itertuples():
                field = backgrounds[int(trial.realization)] + trial.input_amplitude_W_m2 * filtered_packet
                events, _, support = we.detect_events(field,
                    detection_sigma=cfg["detection"]["detection_sigma"], return_support=True)
                matched = matching_significant_event(events, field, truth, trial.threshold_W_m2)
                counts.append(dict(wave_type=wave, case=int(case), realization=int(trial.realization),
                                   h_eff_m=truth["h_eff_m"], K_global=truth["K_global"], width=first.width,
                                   matched_significant=matched is not None))
                if matched is not None:
                    alternative = (dict(we.FIT_OPTIONS[wave], **SHORT_WINDOWS[wave])
                                   if args.candidate == "short" else we.FIT_OPTIONS[wave])
                    for label, options in (("original", we.FIT_OPTIONS[wave]),
                                           (args.candidate, alternative)):
                        rows.append(fit_peak(field, matched, wave, options, label=label,
                                             source="injected_8sigma", realization=int(trial.realization),
                                             case=int(case), truth=truth,
                                             support_indices=support[matched["event_id"]]
                                             if label == "component" else None))
            print(f"{wave}: case {int(case)} ({len(group)} backgrounds)", flush=True)
        del backgrounds
    fits = pd.DataFrame(rows)
    recovery = pd.DataFrame(counts)
    summary = []
    for wave in WAVES:
        n = len(recovery.loc[recovery.wave_type == wave])
        recovered = int(recovery.loc[recovery.wave_type == wave, "matched_significant"].sum())
        for label in ("original", args.candidate):
            injected = fits.loc[(fits.wave_type == wave) & (fits.source == "injected_8sigma")
                                & (fits.window == label)]
            good = injected.loc[injected.quality]
            null = fits.loc[(fits.wave_type == wave) & (fits.source == "null_record_maximum")
                            & (fits.window == label)]
            summary.append(dict(wave_type=wave, window=label, in_band_trials=n,
                matched_significant=recovered, quality_injected=int(injected.quality.sum()),
                quality_injected_within_25pct=int((good.relative_depth_error <= .25).sum()),
                quality_injected_within_factor2=int(
                    ((good.h_eff_m / good.true_h_eff_m >= .5)
                     & (good.h_eff_m / good.true_h_eff_m <= 2)).sum()),
                median_relative_depth_error=float(good.relative_depth_error.median()) if len(good) else None,
                null_record_maxima=len(null), quality_null_record_maxima=int(null.quality.sum())))
    output.mkdir(parents=True)
    fits.to_csv(output / "fits.csv", index=False)
    recovery.to_csv(output / "recovery.csv", index=False)
    pd.DataFrame(summary).to_csv(output / "summary.csv", index=False)
    (output / "provenance.json").write_text(json.dumps(dict(
        completed_utc=datetime.now(timezone.utc).isoformat(), calibration_run=str(RUN),
        ensemble_id=report["ensemble_id"], amplitude_sigma=8,
        inclusion="carrier inside box and strictly below daily Nyquist",
        matched_event="largest above-threshold detected peak within one packet sigma",
        null_control="largest fit-window-eligible candidate in each of eight reused backgrounds; not a 1% false-fit estimate",
        candidate_fit=args.candidate,
        original_fit_options={wave: {key: value.tolist() if isinstance(value, np.ndarray) else value
                                     for key, value in we.FIT_OPTIONS[wave].items()} for wave in WAVES},
        short_windows=SHORT_WINDOWS,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        calibration_code_sha256=cfg["code_sha256"]), indent=2) + "\n")
    print(pd.DataFrame(summary).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
