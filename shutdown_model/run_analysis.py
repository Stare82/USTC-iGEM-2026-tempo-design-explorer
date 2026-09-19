"""Run feasibility, tunability, load, and growth audits for shutdown module."""

from __future__ import annotations

import json
from dataclasses import asdict, replace
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from extract_waveform import extract_waveform
from shutdown_core import (
    CounterConfig,
    ShutdownConfig,
    config_dict,
    cycle_switch_audit,
    observables,
    parameter_provenance,
    pulse_metrics,
    simulate_downstream,
    simulate_full,
    toggle_score,
)


HERE = Path(__file__).resolve().parent


def median_metric(metrics: pd.DataFrame, column: str) -> float:
    if metrics.empty or column not in metrics:
        return np.nan
    values = metrics[column].replace([np.inf, -np.inf], np.nan).dropna()
    return float(values.median()) if len(values) else np.nan


def save_baseline_figure(
    out_dir: Path,
    baseline: pd.DataFrame,
    no_shutdown: pd.DataFrame,
    waveform,
) -> None:
    t = baseline["time_h"].to_numpy(float)
    fig, axes = plt.subplots(6, 1, figsize=(13, 15), sharex=True)
    axes[0].plot(
        t,
        waveform.translation_flux(t),
        color="0.55",
        lw=1.0,
        label="raw C31 translation flux",
    )
    axes[0].plot(
        t,
        waveform.effective_source(t),
        color="tab:purple",
        lw=1.1,
        label=f"Int production (RBS x{waveform.rbs_scale:g})",
    )
    axes[0].set_ylabel("C31 flux\n(uM/h)")
    axes[0].legend(fontsize=8, loc="upper right")
    axes[1].plot(t, baseline["LR_fraction"], color="tab:red", label="LR")
    axes[1].plot(t, baseline["PB_fraction"], color="tab:blue", label="PB")
    axes[1].set_ylabel("DNA fraction")
    axes[1].set_ylim(-0.05, 1.05)
    axes[1].legend(ncol=2, loc="upper right")
    axes[2].plot(t, baseline["BM3R1_total_uM"], color="tab:cyan", label="BM3R1 total")
    axes[2].plot(t, baseline["BM3R1_free_uM"], color="tab:blue", ls="--", label="BM3R1 free")
    axes[2].plot(t, baseline["RDF_total_uM"], color="black", label="RDF")
    axes[2].set_ylabel("Regulators\n(uM)")
    axes[2].legend(ncol=3, fontsize=8)
    axes[3].plot(t, baseline["shutdown_promoter_activity"], color="tab:orange", label="promoter activity")
    axes[3].plot(t, baseline["shutdown_sRNA_uM"], color="tab:green", label="sRNA (uM)")
    axes[3].set_ylabel("Shutdown")
    axes[3].legend(ncol=2, fontsize=8)
    axes[4].plot(t, baseline["output_translation_uM_h"], color="tab:red", label="with shutdown")
    axes[4].plot(t, no_shutdown["output_translation_uM_h"], color="0.45", ls="--", label="no shutdown")
    axes[4].set_ylabel("Translation\n(uM/h)")
    axes[4].legend(ncol=2, fontsize=8)
    axes[5].plot(t, baseline["output_protein_uM"], color="tab:red", label="with shutdown")
    axes[5].plot(t, no_shutdown["output_protein_uM"], color="0.45", ls="--", label="no shutdown")
    axes[5].set_ylabel("Normalized\nprotein")
    axes[5].set_xlabel("Time (h)")
    axes[5].legend(ncol=2, fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.suptitle("Oscillator-driven BM3R1/RDF counter with delayed sRNA shutdown")
    fig.tight_layout()
    fig.savefig(out_dir / "01_baseline_timecourse.png", dpi=170)
    plt.close(fig)


def tunability_scan(
    out_dir: Path,
    driver: pd.DataFrame,
    params: dict,
    base: ShutdownConfig,
) -> tuple[pd.DataFrame, dict]:
    threshold_values = np.geomspace(0.5, 2.0, 13)
    srna_tx_values = np.geomspace(0.006, 0.6, 13)
    reference_peak = median_metric(pulse_metrics(driver), "production_peak")
    rows = []
    for threshold_ratio in threshold_values:
        for srna_tx in srna_tx_values:
            cfg = replace(
                base,
                threshold_ratio_to_rdf=float(threshold_ratio),
                srna_max_tx_uM_h=float(srna_tx),
            )
            downstream = simulate_downstream(driver, params, cfg)
            metrics = pulse_metrics(downstream)
            width = median_metric(metrics, "production_fwhm_h")
            protein_width = median_metric(metrics, "protein_fwhm_h")
            closure = median_metric(metrics, "production_closure_ratio")
            peak = median_metric(metrics, "production_peak")
            auc = median_metric(metrics, "production_auc")
            feasible = bool(
                np.isfinite(width)
                and width > 0.05
                and np.isfinite(closure)
                and closure < 0.20
                and np.isfinite(peak)
                and np.isfinite(reference_peak)
                and peak >= 0.50 * reference_peak
            )
            rows.append(
                {
                    "threshold_ratio_to_rdf": threshold_ratio,
                    "srna_max_tx_uM_h": srna_tx,
                    "production_fwhm_h": width,
                    "protein_fwhm_h": protein_width,
                    "production_closure_ratio": closure,
                    "production_peak_uM_h": peak,
                    "production_auc_uM": auc,
                    "feasible": feasible,
                }
            )
    scan = pd.DataFrame(rows)
    scan.to_csv(out_dir / "02_tunability_scan.csv", index=False)

    matrix = scan.pivot(
        index="threshold_ratio_to_rdf",
        columns="srna_max_tx_uM_h",
        values="production_fwhm_h",
    )
    feasible_matrix = scan.pivot(
        index="threshold_ratio_to_rdf",
        columns="srna_max_tx_uM_h",
        values="feasible",
    )
    values = matrix.to_numpy(float)
    values[~feasible_matrix.to_numpy(bool)] = np.nan
    fig, ax = plt.subplots(figsize=(9, 6.5))
    masked = np.ma.masked_invalid(values)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("0.88")
    im = ax.imshow(masked, origin="lower", aspect="auto", cmap=cmap)
    ax.set_xticks(np.arange(len(matrix.columns)))
    ax.set_xticklabels([f"{x:.2g}" for x in matrix.columns], rotation=45, ha="right")
    ax.set_yticks(np.arange(len(matrix.index)))
    ax.set_yticklabels([f"{x:.3g}" for x in matrix.index])
    ax.set_xlabel("Maximum sRNA transcription (uM/h)")
    ax.set_ylabel("Shutdown/RDF BM3R1 threshold ratio")
    ax.set_title("Feasible output-production pulse width\n(gray = fails closure/peak criteria)")
    fig.colorbar(im, ax=ax, label="Median production FWHM (h)")
    fig.tight_layout()
    fig.savefig(out_dir / "02_tunability_heatmap.png", dpi=170)
    plt.close(fig)

    nearest_srna_tx = float(
        min(srna_tx_values, key=lambda value: abs(value - base.srna_max_tx_uM_h))
    )
    nearest_threshold = float(
        min(
            threshold_values,
            key=lambda value: abs(value - base.threshold_ratio_to_rdf),
        )
    )
    threshold_slice = scan[
        np.isclose(scan["srna_max_tx_uM_h"], nearest_srna_tx)
    ].sort_values(
        "threshold_ratio_to_rdf"
    )
    srna_tx_slice = scan[
        np.isclose(scan["threshold_ratio_to_rdf"], nearest_threshold)
    ].sort_values(
        "srna_max_tx_uM_h"
    )
    threshold_rho = float(
        spearmanr(
            threshold_slice["threshold_ratio_to_rdf"],
            threshold_slice["production_fwhm_h"],
        ).statistic
    )
    srna_tx_rho = float(
        spearmanr(
            srna_tx_slice["srna_max_tx_uM_h"],
            srna_tx_slice["production_fwhm_h"],
        ).statistic
    )
    response = pd.concat(
        [
            threshold_slice.assign(
                slice_name=f"srna_max_tx_uM_h={nearest_srna_tx:.4g}"
            ),
            srna_tx_slice.assign(
                slice_name=f"threshold_ratio_to_rdf={nearest_threshold:.4g}"
            ),
        ],
        ignore_index=True,
    )
    response.to_csv(out_dir / "02_parameter_response.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    axes[0].semilogx(
        threshold_slice["threshold_ratio_to_rdf"],
        threshold_slice["production_fwhm_h"],
        "o-",
        color="tab:blue",
    )
    axes[0].scatter(
        threshold_slice.loc[~threshold_slice["feasible"], "threshold_ratio_to_rdf"],
        threshold_slice.loc[~threshold_slice["feasible"], "production_fwhm_h"],
        facecolors="none",
        edgecolors="tab:red",
        s=80,
        label="fails closure/peak",
    )
    axes[0].set_xlabel("Shutdown/RDF BM3R1 threshold ratio")
    axes[0].set_ylabel("Production FWHM (h)")
    axes[0].set_title(f"Threshold control (Spearman rho={threshold_rho:.2f})")
    axes[0].legend(fontsize=8)
    axes[1].semilogx(
        srna_tx_slice["srna_max_tx_uM_h"],
        srna_tx_slice["production_fwhm_h"],
        "o-",
        color="tab:green",
    )
    axes[1].scatter(
        srna_tx_slice.loc[~srna_tx_slice["feasible"], "srna_max_tx_uM_h"],
        srna_tx_slice.loc[~srna_tx_slice["feasible"], "production_fwhm_h"],
        facecolors="none",
        edgecolors="tab:red",
        s=80,
        label="fails closure/peak",
    )
    axes[1].set_xlabel("Maximum sRNA transcription (uM/h)")
    axes[1].set_ylabel("Production FWHM (h)")
    axes[1].set_title(f"sRNA transcription control (Spearman rho={srna_tx_rho:.2f})")
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.25)
    fig.suptitle("One-parameter control of output duration")
    fig.tight_layout()
    fig.savefig(out_dir / "02_parameter_response.png", dpi=170)
    plt.close(fig)

    feasible_widths = scan.loc[scan["feasible"], "production_fwhm_h"].dropna()
    summary = {
        "grid_points": int(len(scan)),
        "feasible_points": int(scan["feasible"].sum()),
        "feasible_fraction": float(scan["feasible"].mean()),
        "min_feasible_width_h": float(feasible_widths.min()) if len(feasible_widths) else np.nan,
        "max_feasible_width_h": float(feasible_widths.max()) if len(feasible_widths) else np.nan,
        "threshold_ratio_to_width_spearman_rho": threshold_rho,
        "srna_tx_to_width_spearman_rho": srna_tx_rho,
        "threshold_slice_srna_tx_uM_h": nearest_srna_tx,
        "srna_tx_slice_threshold_ratio": nearest_threshold,
    }
    return scan, summary


def transcription_balance_scan(
    out_dir: Path,
    driver: pd.DataFrame,
    params: dict,
    base: ShutdownConfig,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """Map feasibility over literature-scale, explicitly declared transcription ranges."""
    srna_tx_values = np.geomspace(0.006, 0.6, 15)
    output_tx_values = np.geomspace(0.006, 0.06, 15)
    reference_peaks = {}
    for output_tx in output_tx_values:
        control = simulate_downstream(
            driver,
            params,
            replace(
                base,
                srna_max_tx_uM_h=0.0,
                copy_ratio_to_counter=0.0,
                output_max_tx_uM_h=float(output_tx),
            ),
        )
        reference_peaks[float(output_tx)] = median_metric(
            pulse_metrics(control), "production_peak"
        )

    rows = []
    for srna_tx in srna_tx_values:
        for output_tx in output_tx_values:
            cfg = replace(
                base,
                srna_max_tx_uM_h=float(srna_tx),
                output_max_tx_uM_h=float(output_tx),
            )
            result = simulate_downstream(driver, params, cfg)
            metrics = pulse_metrics(result)
            width = median_metric(metrics, "production_fwhm_h")
            closure = median_metric(metrics, "production_closure_ratio")
            peak = median_metric(metrics, "production_peak")
            reference_peak = reference_peaks[float(output_tx)]
            peak_retention = peak / reference_peak if reference_peak > 0 else np.nan
            feasible = bool(
                np.isfinite(width)
                and width > 0.05
                and np.isfinite(closure)
                and closure < 0.20
                and np.isfinite(peak_retention)
                and peak_retention >= 0.50
            )
            rows.append(
                {
                    "srna_max_tx_uM_h": srna_tx,
                    "output_max_tx_uM_h": output_tx,
                    "srna_to_output_tx_ratio": srna_tx / output_tx,
                    "production_fwhm_h": width,
                    "production_closure_ratio": closure,
                    "production_peak_uM_h": peak,
                    "no_shutdown_peak_uM_h": reference_peak,
                    "peak_retention": peak_retention,
                    "feasible": feasible,
                }
            )
    scan = pd.DataFrame(rows)
    scan.to_csv(out_dir / "04_transcription_balance_scan.csv", index=False)

    minimum_rows = []
    for output_tx, group in scan.groupby("output_max_tx_uM_h"):
        passing = group[group["feasible"]].sort_values("srna_max_tx_uM_h")
        minimum_rows.append(
            {
                "output_max_tx_uM_h": output_tx,
                "minimum_feasible_srna_tx_uM_h": (
                    float(passing.iloc[0]["srna_max_tx_uM_h"])
                    if len(passing)
                    else np.nan
                ),
                "minimum_feasible_supply_ratio": (
                    float(passing.iloc[0]["srna_to_output_tx_ratio"])
                    if len(passing)
                    else np.nan
                ),
            }
        )
    minimum_rule = pd.DataFrame(minimum_rows)
    minimum_rule.to_csv(out_dir / "04_minimum_srna_design_rule.csv", index=False)

    closure = scan.pivot(
        index="srna_max_tx_uM_h",
        columns="output_max_tx_uM_h",
        values="production_closure_ratio",
    )
    width = scan.pivot(
        index="srna_max_tx_uM_h",
        columns="output_max_tx_uM_h",
        values="production_fwhm_h",
    )
    feasible = scan.pivot(
        index="srna_max_tx_uM_h",
        columns="output_max_tx_uM_h",
        values="feasible",
    )
    x = closure.columns.to_numpy(float)
    y = closure.index.to_numpy(float)
    xx, yy = np.meshgrid(x, y)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.7), sharex=True, sharey=True)
    closure_values = np.clip(closure.to_numpy(float), 0.0, 1.0)
    left = axes[0].pcolormesh(
        xx,
        yy,
        closure_values,
        shading="nearest",
        cmap="RdYlGn_r",
        vmin=0.0,
        vmax=1.0,
    )
    axes[0].contour(xx, yy, closure_values, levels=[0.20], colors="black", linewidths=1.3)
    fig.colorbar(left, ax=axes[0], label="Late-LR production / peak")
    axes[0].set_title("Shutdown completeness")

    width_values = width.to_numpy(float)
    width_values[~feasible.to_numpy(bool)] = np.nan
    width_cmap = plt.get_cmap("viridis").copy()
    width_cmap.set_bad("0.84")
    right = axes[1].pcolormesh(
        xx,
        yy,
        np.ma.masked_invalid(width_values),
        shading="nearest",
        cmap=width_cmap,
    )
    fig.colorbar(right, ax=axes[1], label="Production FWHM (h)")
    axes[1].set_title("Controllable duration\n(gray = fails closure/peak criteria)")
    for axis in axes:
        axis.plot(x, x, ls="--", color="white", lw=1.5, label=r"$\alpha_S=\alpha_M$")
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlabel(r"Target-mRNA maximum transcription, $\alpha_M$ (uM/h)")
        axis.grid(alpha=0.18)
        axis.legend(loc="upper left", fontsize=8)
    axes[0].set_ylabel(r"sRNA maximum transcription, $\alpha_S$ (uM/h)")
    fig.suptitle("Literature-scale transcription-supply phase map")
    fig.tight_layout()
    fig.savefig(out_dir / "04_transcription_balance_phase_map.png", dpi=190)
    plt.close(fig)

    feasible_rows = scan[scan["feasible"]]
    boundary_ratios = minimum_rule["minimum_feasible_supply_ratio"].dropna()
    summary = {
        "grid_points": int(len(scan)),
        "feasible_points": int(scan["feasible"].sum()),
        "feasible_fraction": float(scan["feasible"].mean()),
        "output_tx_levels_with_a_feasible_solution": int(
            minimum_rule["minimum_feasible_srna_tx_uM_h"].notna().sum()
        ),
        "min_feasible_supply_ratio": (
            float(feasible_rows["srna_to_output_tx_ratio"].min())
            if len(feasible_rows)
            else np.nan
        ),
        "median_minimum_supply_ratio_by_output_level": (
            float(boundary_ratios.median()) if len(boundary_ratios) else np.nan
        ),
        "max_minimum_supply_ratio_by_output_level": (
            float(boundary_ratios.max()) if len(boundary_ratios) else np.nan
        ),
        "fixed_threshold_ratio_to_rdf": base.threshold_ratio_to_rdf,
        "interpretation": (
            "The boundary is conditional on the nominal counter driver and the declared closure/peak criteria."
        ),
    }
    return scan, minimum_rule, summary


def load_audit(
    out_dir: Path,
    waveform_csv: Path,
    counter: CounterConfig,
    base: ShutdownConfig,
    core_score: float,
) -> pd.DataFrame:
    ratios = [0.0, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0]
    rows = []
    for ratio in ratios:
        cfg = replace(base, copy_ratio_to_counter=ratio)
        sol, params, waveform = simulate_full(waveform_csv, counter, cfg)
        obs = observables(sol, params, cfg)
        _, troughs = waveform.peaks_and_troughs()
        score = toggle_score(
            obs["time_h"].to_numpy(), obs["LR_fraction"].to_numpy(), troughs
        )
        pulses = pulse_metrics(obs)
        rows.append(
            {
                "shutdown_copy_ratio": ratio,
                "toggle_score": score["score"],
                "score_change_from_no_load": score["score"] - core_score,
                "production_fwhm_h": median_metric(pulses, "production_fwhm_h"),
                "production_closure_ratio": median_metric(pulses, "production_closure_ratio"),
                "min_free_to_total_BM3R1": float(
                    np.min(
                        obs["BM3R1_free_uM"].to_numpy()
                        / np.maximum(obs["BM3R1_total_uM"].to_numpy(), 1e-12)
                    )
                ),
            }
        )
    audit = pd.DataFrame(rows)
    audit.to_csv(out_dir / "03_operator_load_audit.csv", index=False)
    fig, ax1 = plt.subplots(figsize=(8, 5.5))
    ax1.plot(audit["shutdown_copy_ratio"], audit["toggle_score"], "o-", color="tab:blue")
    ax1.axhline(core_score, color="0.5", ls="--", lw=1, label="no-load score")
    ax1.set_xlabel("Shutdown cassette copy ratio to counter")
    ax1.set_ylabel("Counter toggle score", color="tab:blue")
    finite_scores = audit["toggle_score"].replace([np.inf, -np.inf], np.nan).dropna()
    if len(finite_scores):
        lower = max(0.0, float(finite_scores.min()) - 0.01)
        upper = min(1.01, float(finite_scores.max()) + 0.01)
        if upper - lower < 0.02:
            lower = max(0.0, upper - 0.02)
        ax1.set_ylim(lower, upper)
    ax1.tick_params(axis="y", labelcolor="tab:blue")
    ax2 = ax1.twinx()
    ax2.plot(audit["shutdown_copy_ratio"], audit["production_fwhm_h"], "s-", color="tab:red")
    ax2.set_ylabel("Production pulse FWHM (h)", color="tab:red")
    ax2.tick_params(axis="y", labelcolor="tab:red")
    ax1.grid(alpha=0.25)
    fig.suptitle("Retroactivity audit: added BM3R1 operator load")
    fig.tight_layout()
    fig.savefig(out_dir / "03_operator_load_audit.png", dpi=170)
    plt.close(fig)
    return audit


def growth_audit(
    out_dir: Path,
    waveform_csv: Path,
    base_counter: CounterConfig,
    base_shutdown: ShutdownConfig,
) -> pd.DataFrame:
    rows = []
    for doubling_time_min in [40.0, 50.0, 60.0]:
        k_dil = float(np.log(2.0) / (doubling_time_min / 60.0))
        counter = replace(base_counter, k_dil_h=k_dil)
        sol, params, waveform = simulate_full(waveform_csv, counter, base_shutdown)
        obs = observables(sol, params, base_shutdown)
        _, troughs = waveform.peaks_and_troughs()
        score = toggle_score(
            obs["time_h"].to_numpy(), obs["LR_fraction"].to_numpy(), troughs
        )
        pulses = pulse_metrics(obs)
        rows.append(
            {
                "doubling_time_min": doubling_time_min,
                "nominal_growth_dilution_h": k_dil,
                "production_fwhm_h": median_metric(pulses, "production_fwhm_h"),
                "protein_fwhm_h": median_metric(pulses, "protein_fwhm_h"),
                "production_closure_ratio": median_metric(pulses, "production_closure_ratio"),
                "counter_toggle_score": score["score"],
                "scope": "full 41-state rerun without refitting the nominal counter parameters",
            }
        )
    audit = pd.DataFrame(rows)
    audit.to_csv(out_dir / "04_growth_sensitivity.csv", index=False)
    return audit


def main() -> None:
    waveform_csv = extract_waveform()
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = HERE / "outputs" / f"shutdown_analysis_{stamp}"
    out_dir.mkdir(parents=True, exist_ok=True)

    counter = CounterConfig()
    shutdown = ShutdownConfig()
    sol, params, waveform = simulate_full(waveform_csv, counter, shutdown)
    baseline = observables(sol, params, shutdown)
    baseline.to_csv(out_dir / "01_baseline_timecourse.csv", index=False)

    no_shutdown_cfg = replace(
        shutdown,
        srna_max_tx_uM_h=0.0,
        copy_ratio_to_counter=0.0,
    )
    no_sol, no_params, _ = simulate_full(waveform_csv, counter, no_shutdown_cfg)
    no_shutdown = observables(no_sol, no_params, no_shutdown_cfg)
    no_shutdown.to_csv(out_dir / "01_no_shutdown_control.csv", index=False)
    save_baseline_figure(out_dir, baseline, no_shutdown, waveform)

    _, troughs = waveform.peaks_and_troughs()
    core_toggle = toggle_score(
        no_shutdown["time_h"].to_numpy(),
        no_shutdown["LR_fraction"].to_numpy(),
        troughs,
    )
    coupled_toggle = toggle_score(
        baseline["time_h"].to_numpy(),
        baseline["LR_fraction"].to_numpy(),
        troughs,
    )
    no_load_cycle_audit, no_load_cycle_summary = cycle_switch_audit(
        no_shutdown,
        waveform,
    )
    loaded_cycle_audit, loaded_cycle_summary = cycle_switch_audit(
        baseline,
        waveform,
    )
    no_load_cycle_audit.to_csv(
        out_dir / "01_counter_cycle_audit_no_shutdown.csv", index=False
    )
    loaded_cycle_audit.to_csv(
        out_dir / "01_counter_cycle_audit_with_shutdown.csv", index=False
    )
    baseline_pulses = pulse_metrics(baseline)
    no_shutdown_pulses = pulse_metrics(no_shutdown)
    baseline_pulses.to_csv(out_dir / "01_baseline_pulse_metrics.csv", index=False)
    no_shutdown_pulses.to_csv(out_dir / "01_no_shutdown_pulse_metrics.csv", index=False)
    pd.DataFrame(parameter_provenance(counter, shutdown)).to_csv(
        out_dir / "00_parameter_provenance.csv", index=False
    )

    scan, scan_summary = tunability_scan(out_dir, no_shutdown, params, shutdown)
    balance, minimum_rule, balance_summary = transcription_balance_scan(
        out_dir, no_shutdown, params, shutdown
    )
    load = load_audit(
        out_dir,
        waveform_csv,
        counter,
        shutdown,
        core_score=core_toggle["score"],
    )
    growth = growth_audit(out_dir, waveform_csv, counter, shutdown)

    default_width = median_metric(baseline_pulses, "production_fwhm_h")
    default_closure = median_metric(baseline_pulses, "production_closure_ratio")
    counter_score_drop = coupled_toggle["score"] - core_toggle["score"]
    default_feasible = bool(
        np.isfinite(default_width)
        and default_width > 0.05
        and np.isfinite(default_closure)
        and default_closure < 0.20
        and np.isfinite(counter_score_drop)
        and counter_score_drop > -0.05
        and no_load_cycle_summary["system_pass"]
        and loaded_cycle_summary["system_pass"]
    )
    tunable = bool(
        scan_summary["feasible_points"] >= 5
        and scan_summary["max_feasible_width_h"]
        - scan_summary["min_feasible_width_h"]
        >= 0.5
        and abs(scan_summary["threshold_ratio_to_width_spearman_rho"]) >= 0.8
        and abs(scan_summary["srna_tx_to_width_spearman_rho"]) >= 0.8
    )
    summary = {
        "model_feasible_at_default": default_feasible,
        "duration_tunable_in_scan": tunable,
        "default_production_fwhm_h": default_width,
        "default_production_closure_ratio": default_closure,
        "counter_toggle_no_load": core_toggle,
        "counter_toggle_with_default_shutdown_load": coupled_toggle,
        "counter_score_change": counter_score_drop,
        "counter_cycle_audit_no_shutdown": no_load_cycle_summary,
        "counter_cycle_audit_with_shutdown": loaded_cycle_summary,
        "tunability": scan_summary,
        "transcription_balance": balance_summary,
        "interpretation_boundary": (
            "Feasibility is model-level. RNA kinetic parameters use the Levine RyhB/sodB proxy; design assumptions remain unfitted."
        ),
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True), encoding="utf-8"
    )
    metadata = config_dict(counter, shutdown)
    metadata["waveform_csv"] = str(waveform_csv)
    metadata["outputs"] = str(out_dir)
    metadata["scan_rows"] = int(len(scan))
    metadata["transcription_balance_rows"] = int(len(balance))
    metadata["transcription_boundary_rows"] = int(len(minimum_rule))
    metadata["load_audit_rows"] = int(len(load))
    metadata["growth_audit_rows"] = int(len(growth))
    (out_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2), encoding="utf-8"
    )

    report = f"""# Shutdown model result summary

- Model-feasible at the default exploratory parameters: **{default_feasible}**
- Duration tunable in the parameter scan: **{tunable}**
- Default translation-production FWHM: **{default_width:.3g} h**
- Output remaining late in LR: **{default_closure:.1%} of peak**
- Counter toggle score without shutdown load: **{core_toggle['score']:.4f}**
- Counter toggle score with default shutdown load: **{coupled_toggle['score']:.4f}**
- Audited oscillator cycles without shutdown load: **{no_load_cycle_summary['audited_cycles']}**
- Exactly one complete reversal without shutdown load: **{no_load_cycle_summary['fraction_one_and_complete']:.1%}**
- Exactly one complete reversal with shutdown load: **{loaded_cycle_summary['fraction_one_and_complete']:.1%}**
- Multiple reversals with shutdown load: **{loaded_cycle_summary['fraction_multiple_crossings']:.1%}**
- Alternation fidelity with shutdown load: **{loaded_cycle_summary['alternation_fidelity']:.1%}**
- Feasible scan points: **{scan_summary['feasible_points']}/{scan_summary['grid_points']}**
- Feasible production-width range: **{scan_summary['min_feasible_width_h']:.3g}-{scan_summary['max_feasible_width_h']:.3g} h**
- Threshold ratio vs width Spearman rho: **{scan_summary['threshold_ratio_to_width_spearman_rho']:.3f}**
- sRNA transcription vs width Spearman rho: **{scan_summary['srna_tx_to_width_spearman_rho']:.3f}**
- Feasible transcription-supply combinations: **{balance_summary['feasible_points']}/{balance_summary['grid_points']}**
- Smallest feasible sRNA/target-mRNA transcription ratio: **{balance_summary['min_feasible_supply_ratio']:.3g}**
- Median boundary ratio across target-mRNA levels: **{balance_summary['median_minimum_supply_ratio_by_output_level']:.3g}**

These are model-level results. RNA kinetic parameters are anchored to the RyhB/sodB proxy in Levine et al. (2007). The BM3R1 concentration calibration and cassette load assumptions still require construct-specific measurements before biological predictions are made.
"""
    (out_dir / "RESULTS.md").write_text(report, encoding="utf-8")
    (HERE / "latest_run.txt").write_text(str(out_dir), encoding="utf-8")

    print(json.dumps(summary, indent=2, allow_nan=True))
    print(f"Outputs: {out_dir}")


if __name__ == "__main__":
    main()
