#!/usr/bin/env python3
"""Summarize completed paired carrier-validation CSVs in this directory."""

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent


def summarize(group):
    matched = group.loc[group.matched]
    significant = matched.loc[matched.significant]
    fitted = significant.loc[significant.fit_success & significant.abs_K_error.notna()]
    return pd.Series(dict(
        trials=len(group), matched=len(matched), significant=len(significant), fitted=len(fitted),
        boundary=int(fitted.at_search_edge.sum()),
        median_signed_K_error=(fitted.estimated_K_global - fitted.true_K_global).median(),
        median_abs_K_error=fitted.abs_K_error.median(),
        p90_abs_K_error=fitted.abs_K_error.quantile(.9),
        median_signed_relative_frequency_error=(
            (fitted.estimated_frequency_cpd - fitted.true_frequency_cpd)
            / fitted.true_frequency_cpd).median(),
        median_relative_frequency_error=fitted.relative_frequency_error.median(),
        p90_relative_frequency_error=fitted.relative_frequency_error.quantile(.9),
        large_joint_error_fraction=(
            (fitted.abs_K_error > .5) | (fitted.relative_frequency_error > .1)
        ).mean() if len(fitted) else np.nan,
    ))


def main():
    files = sorted(HERE.glob("trials.*.csv"))
    files = [path for path in files if ".pilot." not in path.name]
    if len(files) != 15:
        raise RuntimeError(f"Expected all 15 non-pilot trial CSVs, found {len(files)}")
    trials = pd.concat([pd.read_csv(path) for path in files], ignore_index=True)
    injected = trials.loc[trials.noise_type.eq("gaussian_background")].copy()
    for column in ("matched", "significant", "fit_success", "at_search_edge"):
        injected[column] = injected[column].map({True: True, False: False,
            "True": True, "False": False}).fillna(False).astype(bool)
    # Each trial is paired: both methods must have identical detection outcomes.
    keys = ["source", "experiment", "wave", "case", "morphology", "background_replicate", "peak_ratio_target"]
    sizes = injected.groupby(keys).size()
    if not sizes.eq(2).all():
        raise RuntimeError("Incomplete raw/phase trial pairs.")
    groups = ["source", "experiment", "wave", "morphology", "peak_ratio_target", "method"]
    detail = injected.groupby(groups, dropna=False).apply(summarize).reset_index()
    detail.to_csv(HERE / "summary_by_shape.csv", index=False)
    wave = injected.groupby(["wave", "method"], dropna=False).apply(summarize).reset_index()
    wave.to_csv(HERE / "summary_by_wave.csv", index=False)
    by_morphology = injected.groupby(["wave", "morphology", "method"],
        dropna=False).apply(summarize).reset_index()
    by_morphology.to_csv(HERE / "summary_by_morphology.csv", index=False)
    by_strength = injected.groupby(["wave", "peak_ratio_target", "method"],
        dropna=False).apply(summarize).reset_index()
    by_strength.to_csv(HERE / "summary_by_strength.csv", index=False)
    truth = injected[["wave", "case", "true_K_global", "true_frequency_cpd"]].drop_duplicates()
    by_case = injected.groupby(["wave", "case", "morphology", "method"],
        dropna=False).apply(summarize).reset_index().merge(truth, on=["wave", "case"])
    by_case.to_csv(HERE / "summary_by_case.csv", index=False)
    paired = injected.pivot(index=keys, columns="method", values=[
        "matched", "significant", "fit_success", "abs_K_error", "relative_frequency_error"])
    paired.columns = [f"{quantity}_{method}" for quantity, method in paired.columns]
    paired = paired.reset_index()
    if not paired.matched_raw.eq(paired.matched_phase).all() or not paired.significant_raw.eq(paired.significant_phase).all():
        raise RuntimeError("Paired methods do not share detection/significance outcomes.")
    both_fitted = (paired.fit_success_raw & paired.fit_success_phase
                   & paired.abs_K_error_raw.notna() & paired.abs_K_error_phase.notna())
    paired["both_fitted"] = both_fitted
    paired["phase_lower_K_error"] = np.where(both_fitted,
        paired.abs_K_error_phase < paired.abs_K_error_raw, np.nan)
    paired["phase_lower_frequency_error"] = np.where(both_fitted,
        paired.relative_frequency_error_phase < paired.relative_frequency_error_raw, np.nan)
    paired.to_csv(HERE / "paired_errors.csv", index=False)
    paired_significant = paired.loc[paired.significant_raw & paired.both_fitted]
    paired_summary = paired_significant.groupby("wave").agg(
        paired_significant=("both_fitted", "size"),
        phase_lower_K_error_fraction=("phase_lower_K_error", "mean"),
        phase_lower_frequency_error_fraction=("phase_lower_frequency_error", "mean"),
        raw_median_abs_K_error=("abs_K_error_raw", "median"),
        phase_median_abs_K_error=("abs_K_error_phase", "median"),
        raw_median_relative_frequency_error=("relative_frequency_error_raw", "median"),
        phase_median_relative_frequency_error=("relative_frequency_error_phase", "median"),
    ).reset_index()
    paired_summary.to_csv(HERE / "summary_paired_by_wave.csv", index=False)
    clean = trials.loc[trials.noise_type.eq("clean")].copy()
    clean.to_csv(HERE / "clean_checks.csv", index=False)
    for column in ("matched", "significant", "fit_success", "at_search_edge"):
        clean[column] = clean[column].map({True: True, False: False,
            "True": True, "False": False}).fillna(False).astype(bool)
    clean.groupby(["wave", "morphology", "method"], dropna=False).apply(
        summarize).reset_index().to_csv(HERE / "summary_clean.csv", index=False)
    clean.groupby(["wave", "case", "morphology", "method"],
        dropna=False).apply(summarize).reset_index().merge(truth, on=["wave", "case"]).to_csv(
            HERE / "summary_clean_by_case.csv", index=False)

    raw = injected.loc[injected.method.eq("raw") & injected.significant & injected.fit_success].copy()
    raw["r2_at_least_0p3"] = raw.weighted_r2 >= .3
    raw["large_joint_error"] = ((raw.abs_K_error > .5)
                                | (raw.relative_frequency_error > .1))
    raw.groupby(["wave", "r2_at_least_0p3"]).agg(
        fitted=("large_joint_error", "size"),
        large_joint_error_fraction=("large_joint_error", "mean"),
        median_abs_K_error=("abs_K_error", "median"),
        median_relative_frequency_error=("relative_frequency_error", "median"),
    ).reset_index().to_csv(HERE / "raw_r2_screen.csv", index=False)

    # Backgrounds are reused across carrier/envelope cases. Resample whole
    # source/experiment/background blocks, never individual trial pixels/rows.
    random = np.random.default_rng(20261101)
    bootstrap_rows = []
    for wave_name, group in paired_significant.groupby("wave"):
        group = group.copy()
        group["large_raw"] = ((group.abs_K_error_raw > .5)
                              | (group.relative_frequency_error_raw > .1)).astype(float)
        group["large_phase"] = ((group.abs_K_error_phase > .5)
                                | (group.relative_frequency_error_phase > .1)).astype(float)
        blocks = []
        for _, block in group.groupby(["source", "experiment", "background_replicate"]):
            blocks.append(np.column_stack([
                block.abs_K_error_phase - block.abs_K_error_raw,
                block.relative_frequency_error_phase - block.relative_frequency_error_raw,
                block.large_phase - block.large_raw,
            ]))
        if len(blocks) != 18:
            raise RuntimeError(f"Expected 18 source/background blocks for {wave_name}; found {len(blocks)}")
        point = np.concatenate(blocks).mean(axis=0)
        replicates = np.array([np.concatenate([blocks[index] for index in
            random.integers(len(blocks), size=len(blocks))]).mean(axis=0)
            for _ in range(2000)])
        for index, metric in enumerate(("mean_abs_K_error", "mean_relative_frequency_error",
                                        "large_joint_error_fraction")):
            bootstrap_rows.append(dict(wave=wave_name, metric=metric,
                difference="phase_minus_raw", paired_trials=len(group), background_blocks=len(blocks),
                point_difference=float(point[index]),
                lower_95=float(np.quantile(replicates[:, index], .025)),
                upper_95=float(np.quantile(replicates[:, index], .975))))
    pd.DataFrame(bootstrap_rows).to_csv(HERE / "paired_block_bootstrap.csv", index=False)
    controls = trials.loc[trials.noise_type.eq("null_maximum")]
    controls.to_csv(HERE / "null_controls.csv", index=False)
    script_hash = hashlib.sha256((HERE / "run.py").read_bytes()).hexdigest()
    cases_hash = hashlib.sha256((HERE / "cases.tsv").read_bytes()).hexdigest()
    provenance = trials[["source", "experiment", "wave", "calibration_root"]].drop_duplicates().copy()
    provenance["reproduction_script_sha256"] = script_hash
    provenance["cases_sha256"] = cases_hash
    provenance.to_csv(HERE / "provenance.csv", index=False)
    edge_files = sorted((HERE / "er_k1_control").glob("trials.*.csv"))
    if edge_files:
        if len(edge_files) != 3:
            raise RuntimeError(f"Expected all three ER K=-1 control CSVs, found {len(edge_files)}")
        edge = pd.concat([pd.read_csv(path) for path in edge_files], ignore_index=True)
        for column in ("matched", "significant", "fit_success", "at_search_edge"):
            edge[column] = edge[column].map({True: True, False: False,
                "True": True, "False": False}).fillna(False).astype(bool)
        for noise_type in ("clean", "gaussian_background"):
            sample = edge.loc[edge.noise_type.eq(noise_type)]
            sample.groupby(["source", "experiment", "morphology", "method"],
                dropna=False).apply(summarize).reset_index().assign(noise_type=noise_type).to_csv(
                    HERE / f"er_k1_control_{noise_type}.csv", index=False)
        original_k1 = trials.loc[trials.wave.eq("n1ERW") & trials.case.eq(0)].copy()
        comparison = pd.concat([original_k1, edge], ignore_index=True)
        comparison = comparison.loc[comparison.noise_type.isin(("clean", "gaussian_background"))].copy()
        for column in ("matched", "significant", "fit_success", "at_search_edge"):
            comparison[column] = comparison[column].map({True: True, False: False,
                "True": True, "False": False}).fillna(False).astype(bool)
        comparison["envelope_class"] = np.where(
            comparison.morphology.str.startswith("empirical"), "observed_shape", "controlled_shape")
        comparison.groupby(["true_period_days", "envelope_class", "noise_type", "method"],
            dropna=False).apply(summarize).reset_index().to_csv(
                HERE / "er_k1_period_comparison.csv", index=False)
    print("Wrote wave, shape, morphology, strength, clean, paired, bootstrap,")
    print("null-control and provenance CSV diagnostics.")


if __name__ == "__main__":
    main()
