#!/usr/bin/env python3
"""Summarize one completed two-symmetry NOAA significance-validation run.

Example: python scripts/summarize_wave_significance_checks.py --run-label pilot1
Writes descriptive CSV/JSON summaries in that run's data/wave_significance folder.
No hypothesis tests treating individual correlated pixels as independent are used.
"""

import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
import pandas as pd

from src.paths import WAVE_SIGNIFICANCE_DIR


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-label", required=True)
    args = parser.parse_args()
    if not re.fullmatch(r"[A-Za-z0-9_-]+",args.run_label):
        parser.error("Use a simple run label.")
    root = WAVE_SIGNIFICANCE_DIR/"validation"/args.run_label
    distributions, comparisons, injection_summaries, provenance = [], [], [], []
    for symmetry in ("symmetric","antisymmetric"):
        directory = root/symmetry
        report_path = directory/"report.json"
        report = json.loads(report_path.read_text())
        if report["status"]!="complete":
            raise RuntimeError(f"Validation task has not completed: {directory}")
        provenance.append(dict(symmetry=symmetry, report_sha256=hashlib.sha256(report_path.read_bytes()).hexdigest(),
                               nonconverged_surrogates=report["nonconverged_surrogates"]))
        diagnostic = pd.read_csv(directory/"distributions.csv")
        observed = diagnostic.loc[diagnostic.model=="observed"]
        for stage, group in observed.groupby("stage"):
            global_row = group.loc[(group.region=="global")&(group.season=="all")].iloc[0]
            local = group.loc[(group.region!="global")&(group.season!="all")]
            if int(local.samples.sum())!=int(global_row.samples):
                raise ValueError("Season/sector subsets do not partition the global samples.")
            weights = local.samples.to_numpy()/local.samples.sum()
            delta = local.mean_W_m2.to_numpy()-global_row.mean_W_m2
            variance = local.sigma_W_m2.to_numpy()**2
            mixture_variance = np.sum(weights*(variance+delta**2))
            np.testing.assert_allclose(mixture_variance,global_row.sigma_W_m2**2,rtol=1e-8)
            # A descriptive Gaussian-mixture calculation: what pooled kurtosis
            # would these means/variances produce if each sector/season were
            # internally Gaussian? It is not a fitted convective-background null.
            mixture_fourth = np.sum(weights*(3*variance**2+6*variance*delta**2+delta**4))
            gaussian = diagnostic.loc[(diagnostic.model=="gaussian")&(diagnostic.stage==stage)
                                     &(diagnostic.region=="global")&(diagnostic.season=="all")]
            marginal = diagnostic.loc[(diagnostic.model=="empirical_marginal")&(diagnostic.stage==stage)
                                     &(diagnostic.region=="global")&(diagnostic.season=="all")]
            distributions.append(dict(symmetry=symmetry,stage=stage,
                observed_global_excess_kurtosis=float(global_row.excess_kurtosis),
                gaussian_mixture_predicted_excess_kurtosis=float(mixture_fourth/mixture_variance**2-3),
                observed_global_3sigma_fraction=float(global_row.lower_3sigma_fraction+global_row.upper_3sigma_fraction),
                weighted_local_3sigma_fraction=float(np.sum(weights*(local.lower_3sigma_fraction+local.upper_3sigma_fraction))),
                local_sigma_max_over_min=float(local.sigma_W_m2.max()/local.sigma_W_m2.min()),
                gaussian_global_mean_excess_kurtosis=float(gaussian.excess_kurtosis.mean()),
                empirical_marginal_global_mean_excess_kurtosis=float(marginal.excess_kurtosis.mean())))
        thresholds = pd.read_csv(directory/"threshold_comparison.csv")
        events = pd.read_csv(directory/"null_events.csv")
        counts = pd.read_csv(directory/"null_record_counts.csv")
        events = events.merge(counts[["model","realization","wave","sigma_W_m2"]],
                              on=["model","realization","wave"],validate="many_to_one")
        for wave, group in thresholds.loc[thresholds.usable].groupby("wave"):
            if set(group.model)!={"gaussian","empirical_marginal"}:
                continue
            by_model = group.set_index("model")
            g = by_model.loc["gaussian"]
            m = by_model.loc["empirical_marginal"]
            row = dict(wave=wave,gaussian_threshold_W_m2=float(g.threshold_W_m2),
                empirical_marginal_threshold_W_m2=float(m.threshold_W_m2),
                threshold_change_percent=float(100*(m.threshold_W_m2/g.threshold_W_m2-1)),
                filtered_sigma_change_percent=float(100*(m.calibration_mean_sigma_W_m2/g.calibration_mean_sigma_W_m2-1)),
                gaussian_observed_above_threshold=int(g.observed_above_threshold),
                empirical_marginal_observed_above_threshold=int(m.observed_above_threshold),
                gaussian_false_candidates_per_record=float(g.heldout_false_candidates_per_record),
                empirical_marginal_false_candidates_per_record=float(m.heldout_false_candidates_per_record))
            for model in ("gaussian","empirical_marginal"):
                selected = events.loc[(events.wave==wave)&(events.model==model)&(events.split=="calibration")]
                row[model+"_threshold_sigma_units"] = float(np.quantile(selected.peak_amplitude_W_m2/selected.sigma_W_m2,.99,method="higher"))
            comparisons.append(row)
        injections = pd.read_csv(directory/"injections.csv")
        # Explicitly distinguish intentionally out-of-band/aliased controls from
        # carriers that this daily wave mask is intended to retain.
        injections["carrier_below_daily_nyquist"] = injections.frequency_cpd < .5
        injections["in_band_and_resolved"] = injections.carrier_inside_rectangular_bounds & injections.carrier_below_daily_nyquist
        for keys, group in injections.groupby(["wave_type","model","amplitude_sigma","in_band_and_resolved"]):
            wave,model,amplitude,inside = keys
            before = group.baseline_nearby_significant.gt(0)
            after = group.injected_nearby_significant.gt(0)
            injection_summaries.append(dict(wave=wave,model=model,amplitude_sigma=float(amplitude),
                in_band_and_resolved=bool(inside),trials=len(group),
                nearby_candidate_fraction=float(group.injected_nearby_candidates.gt(0).mean()),
                baseline_significant_fraction=float(before.mean()), injected_significant_fraction=float(after.mean()),
                paired_change_significant_fraction=float(after.astype(int).sub(before.astype(int)).mean()),
                mean_retained_packet_energy=float(group.retained_energy_fraction.mean())))
    pd.DataFrame(distributions).to_csv(root/"distribution_summary.csv",index=False)
    pd.DataFrame(comparisons).to_csv(root/"model_comparison.csv",index=False)
    pd.DataFrame(injection_summaries).to_csv(root/"injection_recovery_summary.csv",index=False)
    summary = dict(status="complete",run_label=args.run_label, provenance=provenance,
        summary_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        distributions=distributions,model_comparison=comparisons,injections=injection_summaries,
        limitations="Pilot replication; stationary nulls; empirical marginal includes waves; Gaussian-mixture moment calculation is descriptive; no fitted k/omega validation or production threshold change.")
    (root/"analysis_summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
