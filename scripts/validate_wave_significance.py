#!/usr/bin/env python3
"""Validate NOAA packet significance under Gaussian and empirical-marginal nulls.

Run each symmetry in a Slurm allocation, e.g. from the installed project root:
  python scripts/validate_wave_significance.py --symmetry symmetric --run-label pilot1
  python scripts/validate_wave_significance.py --symmetry antisymmetric --run-label pilot1

Uses NOAA 1990--2010 via data/noaa_olr. Outputs tables and provenance under
data/wave_significance/validation/<run-label>/<symmetry>/. No production
catalog/threshold changes, fitted depths or figures are produced.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import time

import numpy as np
import pandas as pd
from scipy.signal import hilbert
import xarray as xr
import xr_ccew as tw

from src import wave_events as we, wave_significance as sg
from src.paths import NOAA_OLR_FILE, PROJECT_ROOT, WAVE_SIGNIFICANCE_DIR


def filtered_projection(field, wave):
    """Apply the production temporal preprocessing and mask after projection."""
    anomaly = tw.apply_window(tw.remove_mean_and_linear_trend(
        tw.remove_harmonics_of_seasonal_cycle(field, num_harmonics=we.SEASONAL_HARMONICS)),
        window="tukey", pct=we.TAPER_FRACTION)
    return tw.filter_field(anomaly, we.WAVE_PROFILES[wave], apply_profile_symmetry=False)


def nearby_peaks(events, field, truth, threshold):
    """Match catalog peaks within one envelope sigma along the moving packet."""
    if not events:
        return dict(nearby_candidates=0, nearby_significant=0, nearby_maximum_W_m2=0.0)
    indices = field.get_index("time").get_indexer([event["peak_time"] for event in events])
    lag = indices-truth["center_index"]
    lon = np.array([event["peak_lon_deg"] for event in events])
    x = (lon-truth["center_lon_deg"]-truth["group_speed_deg_day"]*lag+180) % 360-180
    matched = (abs(lag)<=truth["time_sigma_days"]) & (abs(x)<=truth["lon_sigma_deg"])
    peaks = np.array([event["peak_amplitude_W_m2"] for event in events])[matched]
    return dict(nearby_candidates=len(peaks), nearby_significant=int(np.count_nonzero(peaks>=threshold)),
                nearby_maximum_W_m2=float(peaks.max()) if len(peaks) else 0.0)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--symmetry", required=True, choices=("symmetric", "antisymmetric"))
    parser.add_argument("--run-label", required=True)
    parser.add_argument("--realizations", type=int, default=8)
    parser.add_argument("--injection-realizations", type=int, default=2)
    parser.add_argument("--start-year", type=int, default=1990)
    parser.add_argument("--end-year", type=int, default=2010)
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--max-iterations", type=int, default=80)
    parser.add_argument("--spectral-tolerance", type=float, default=.03)
    args = parser.parse_args()
    if (not re.fullmatch(r"[A-Za-z0-9_-]+", args.run_label) or args.realizations<4
            or not 1<=args.injection_realizations<=args.realizations-args.realizations//2):
        parser.error("Use a simple run label, >=4 realizations and injections within the validation half.")
    if PROJECT_ROOT != Path(__file__).resolve().parents[1]:
        raise RuntimeError("Install this checkout with pip -e .")
    if not NOAA_OLR_FILE.parent.is_symlink() or not NOAA_OLR_FILE.is_file():
        raise FileNotFoundError("Configure data/noaa_olr.")
    # This initial comparison keeps one smoothing recipe fixed, changing only
    # the marginal constraint. Smoothing sensitivity is a separate comparison.
    frequency_smoothing_cpd, wavenumber_smoothing = 1/96, 2.0
    quantile = .99
    bootstrap_replicates = 200
    depths = [8.0, 25.0, 90.0]
    wavenumbers = {"Kelvin": [4, 10], "n1ERW": [-3, -8], "MRG": [-3, -8],
                   "n0EIG": [3, 10], "n1WIG": [-3, -8]}
    widths = {"compact": (.5, .2), "extended": (1.0, .5)}  # sigma in periods/wavelengths
    amplitudes_sigma = [2.0, 4.0]
    waves = [name for name,profile in we.WAVE_PROFILES.items() if profile.symmetry==args.symmetry]
    output = WAVE_SIGNIFICANCE_DIR/"validation"/args.run_label/args.symmetry
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    stat = NOAA_OLR_FILE.stat()
    config = dict(arguments=vars(args), source=dict(path=str(NOAA_OLR_FILE.relative_to(PROJECT_ROOT)),
        size_bytes=stat.st_size, mtime_ns=stat.st_mtime_ns),
        frequency_smoothing_cpd=frequency_smoothing_cpd, wavenumber_smoothing=wavenumber_smoothing,
        quantile=quantile, bootstrap_replicates=bootstrap_replicates, depths_m=depths,
        wavenumbers=wavenumbers, widths_periods_wavelengths=widths, amplitudes_sigma=amplitudes_sigma,
        profiles={wave:we.WAVE_PROFILES[wave].as_dict() for wave in waves},
        preprocessing=dict(lat_range=we.LAT_RANGE, average_lat_range=we.AVERAGE_LAT_RANGE,
            seasonal_harmonics=we.SEASONAL_HARMONICS, taper_fraction=we.TAPER_FRACTION,
            reference_exclude_days=we.REFERENCE_EXCLUDE_DAYS),
        detection=dict(detection_sigma=we.DETECTION_SIGMA, peak_sigma=we.PEAK_SIGMA,
            min_duration_days=we.MIN_DURATION_DAYS, min_longitude_span_deg=we.MIN_LONGITUDE_SPAN_DEG),
        code_sha256={name:hashlib.sha256(Path(path).read_bytes()).hexdigest() for name,path in
            dict(worker=__file__, significance=sg.__file__, wave_events=we.__file__).items()})
    we.write_json(output/"report.json", dict(status="running", configuration=config))
    with xr.open_dataset(NOAA_OLR_FILE, chunks={"time":365}) as dataset:
        olr = dataset.olr.sel(time=slice(str(args.start_year),str(args.end_year)), lat=slice(*we.LAT_RANGE))
        if (str(olr.time.values[0])[:10]!=f"{args.start_year}-01-01"
                or str(olr.time.values[-1])[:10]!=f"{args.end_year}-12-31"):
            raise ValueError("NOAA source does not cover the requested complete years.")
        projected = sg.project_olr(olr)
    field = projected[args.symmetry]
    del projected
    anomaly = tw.remove_mean_and_linear_trend(tw.remove_harmonics_of_seasonal_cycle(field, num_harmonics=we.SEASONAL_HARMONICS))
    power = sg.estimate_background_power(field, frequency_smoothing_cpd=frequency_smoothing_cpd,
                                          wavenumber_smoothing=wavenumber_smoothing)
    # Match both nulls' expected global variance to the empirical marginal. This
    # separates a distribution-shape comparison from a total-variance difference.
    variance_factor = float(anomaly.var())/float(power.mean())
    power.data *= variance_factor
    power.attrs.update(variance_rescaling_for_marginal_comparison=variance_factor,
                       expected_variance_W2_m4=float(power.mean()))
    power.to_netcdf(output/"target_background.nc")
    distributions, spectra, generation, event_rows, counts = [], [], [], [], []
    observed, injection_backgrounds = {}, {}
    reference = slice(we.REFERENCE_EXCLUDE_DAYS,-we.REFERENCE_EXCLUDE_DAYS)
    for row in sg.distribution_diagnostics(anomaly.isel(time=reference)):
        distributions.append(dict(model="observed", stage="before_filter", realization=-1, **row))
    for wave in waves:
        filtered = filtered_projection(field,wave)
        observed[wave], _ = we.detect_events(filtered)
        for row in sg.distribution_diagnostics(filtered.isel(time=reference)):
            distributions.append(dict(model="observed", stage=wave, realization=-1, **row))
    n_training = args.realizations//2
    for realization in range(args.realizations):
        split = "calibration" if realization<n_training else "validation"
        for model in ("gaussian", "empirical_marginal"):
            seed = np.random.SeedSequence([args.seed, 0 if args.symmetry=="symmetric" else 1, realization])
            rng = np.random.default_rng(seed)
            if model=="gaussian":
                background = sg.draw_background(power,anomaly,rng)
                quality = dict(converged=True, iterations=0, spectral_amplitude_relative_error=None)
            else:
                background, quality = sg.draw_marginal_background(power,anomaly,rng,
                    max_iterations=args.max_iterations, spectral_tolerance=args.spectral_tolerance)
            generation.append(dict(model=model, realization=realization, split=split, **quality))
            for row in sg.spectrum_diagnostics(background,power):
                spectra.append(dict(model=model, realization=realization, **row))
            for row in sg.distribution_diagnostics(background.isel(time=reference)):
                distributions.append(dict(model=model, stage="before_filter", realization=realization, **row))
            print(f"{args.symmetry} {model} {realization+1}/{args.realizations}: {quality}",flush=True)
            # Never calibrate with a failed spectral/marginal constraint match.
            if not quality["converged"]:
                continue
            for wave in waves:
                filtered = filtered_projection(background,wave)
                events, sigma = we.detect_events(filtered)
                counts.append(dict(model=model, realization=realization, split=split, wave=wave,
                                   candidates=len(events), sigma_W_m2=sigma))
                for event in events:
                    event_rows.append(dict(model=model, realization=realization, split=split, wave=wave, **event))
                for row in sg.distribution_diagnostics(filtered.isel(time=reference)):
                    distributions.append(dict(model=model, stage=wave, realization=realization, **row))
                if n_training<=realization<n_training+args.injection_realizations:
                    injection_backgrounds[wave,model,realization] = (filtered,events)
    for name,rows in (("distributions",distributions),("spectra",spectra),("generation",generation),
                      ("null_events",event_rows),("null_record_counts",counts)):
        pd.DataFrame(rows).to_csv(output/f"{name}.csv",index=False)
    events_frame, count_frame = pd.DataFrame(event_rows), pd.DataFrame(counts)
    comparisons, thresholds = [], {}
    for wave in waves:
        for model in ("gaussian","empirical_marginal"):
            selected = events_frame.loc[(events_frame.wave==wave)&(events_frame.model==model)]
            calibration = selected.loc[selected.split=="calibration"]
            validation = count_frame.loc[(count_frame.wave==wave)&(count_frame.model==model)&(count_frame.split=="validation")]
            groups = [group.peak_amplitude_W_m2.to_numpy() for _,group in calibration.groupby("realization")]
            if len(groups)<2 or len(validation)<2:
                comparisons.append(dict(wave=wave,model=model,usable=False,reason="Insufficient converged independent records"))
                continue
            peak_threshold = float(np.quantile(np.concatenate(groups),quantile,method="higher"))
            thresholds[wave,model] = peak_threshold
            above = selected.loc[(selected.split=="validation")&(selected.peak_amplitude_W_m2>=peak_threshold)].groupby("realization").size()
            false_counts = validation.realization.map(above).fillna(0).to_numpy()
            boot_rng = np.random.default_rng(args.seed+100)
            boot_thresholds = [np.quantile(np.concatenate([groups[index] for index in boot_rng.integers(len(groups),size=len(groups))]),quantile,method="higher") for _ in range(bootstrap_replicates)]
            boot_false = [np.mean(boot_rng.choice(false_counts,len(false_counts),replace=True)) for _ in range(bootstrap_replicates)]
            row = dict(wave=wave,model=model,usable=True,threshold_W_m2=peak_threshold,
                threshold_MC_interval_low=float(np.quantile(boot_thresholds,.025)), threshold_MC_interval_high=float(np.quantile(boot_thresholds,.975)),
                heldout_false_candidates_per_record=float(false_counts.mean()),
                false_count_MC_interval_low=float(np.quantile(boot_false,.025)), false_count_MC_interval_high=float(np.quantile(boot_false,.975)),
                observed_candidates=len(observed[wave]), observed_above_threshold=sum(e["peak_amplitude_W_m2"]>=peak_threshold for e in observed[wave]),
                calibration_records=len(groups), validation_records=len(validation), calibration_candidates=len(calibration),
                calibration_mean_sigma_W_m2=float(count_frame.loc[(count_frame.wave==wave)&(count_frame.model==model)&(count_frame.split=="calibration"),"sigma_W_m2"].mean()))
            comparisons.append(row)
            print("COMPARISON",json.dumps(row),flush=True)
    pd.DataFrame(comparisons).to_csv(output/"threshold_comparison.csv",index=False)
    injections = []
    for wave in waves:
        gaussian_records = count_frame.loc[(count_frame.wave==wave)&(count_frame.model=="gaussian")&(count_frame.split=="calibration")]
        amplitude_reference = float(gaussian_records.sigma_W_m2.mean())
        case = 0
        for depth in depths:
            for K in wavenumbers[wave]:
                for width_name,(time_width,lon_width) in widths.items():
                    center_index = field.sizes["time"]//2+(case%5-2)*120
                    center_lon = [60.5,180.5,300.5,359.5][case%4]
                    packet, truth = sg.synthetic_wave_packet(field,wave,K,depth,center_index=center_index,
                        center_lon=center_lon,time_sigma_periods=time_width,lon_sigma_wavelengths=lon_width)
                    filtered_packet = filtered_projection(packet,wave)
                    retained_energy = float((filtered_packet**2).sum()/(packet**2).sum())
                    unit_envelope_peak = float(abs(hilbert(filtered_packet.values,axis=1)).max())
                    for amplitude_sigma in amplitudes_sigma:
                        amplitude = amplitude_sigma*amplitude_reference
                        for (background_wave,model,realization),(background,baseline_events) in injection_backgrounds.items():
                            if background_wave!=wave or (wave,model) not in thresholds:
                                continue
                            threshold = thresholds[wave,model]
                            injected_events, injected_sigma = we.detect_events(background+amplitude*filtered_packet)
                            before = nearby_peaks(baseline_events,field,truth,threshold)
                            after = nearby_peaks(injected_events,field,truth,threshold)
                            injections.append(dict(truth, model=model,realization=realization,case=case,width=width_name,
                                amplitude_sigma=amplitude_sigma,input_amplitude_W_m2=amplitude,
                                filtered_packet_peak_W_m2=amplitude*unit_envelope_peak,retained_energy_fraction=retained_energy,
                                threshold_W_m2=threshold,injected_reference_sigma_W_m2=injected_sigma,
                                **{"baseline_"+key:value for key,value in before.items()},
                                **{"injected_"+key:value for key,value in after.items()}))
                    case += 1
                    print(f"INJECTION {wave} case {case}/12",flush=True)
                    pd.DataFrame(injections).to_csv(output/"injections.csv",index=False)
    report = dict(status="complete",configuration=config, completed_utc=datetime.now(timezone.utc).isoformat(),
        elapsed_seconds=time.monotonic()-started, comparisons=comparisons, injection_trials=len(injections),
        nonconverged_surrogates=sum(row["model"]=="empirical_marginal" and not row["converged"] for row in generation),
        interpretation="Descriptive distributions; conditional null-candidate peak thresholds; MC intervals resample whole null records only, not background-model uncertainty or catalog FDR.",
        limitations="Empirical marginal includes real waves; IAAFT is a covariance/marginal surrogate, not a model of nonlinear convection. Both nulls are stationary. Injection matching uses peak proximity and reports pre-existing matches. Injection tests assess detection, not fitted k/omega or h_eff accuracy. First-screen replication counts are small.")
    we.write_json(output/"report.json",report)
    print("COMPLETE",output,flush=True)


if __name__=="__main__":
    main()
