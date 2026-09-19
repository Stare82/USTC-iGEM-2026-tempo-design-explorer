"""Extended growth-consistency, sensitivity, and robustness analyses."""

from __future__ import annotations

import json
import math
import sys
from dataclasses import asdict, replace
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from scipy.stats import qmc, spearmanr

from extract_waveform import extract_waveform
from shutdown_core import (
    CounterConfig,
    ShutdownConfig,
    cycle_switch_audit,
    observables,
    pulse_metrics,
    simulate_downstream,
    simulate_full,
)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
OSCILLATOR_CODE = ROOT / "oscillator" / "code"
if str(OSCILLATOR_CODE) not in sys.path:
    sys.path.append(str(OSCILLATOR_CODE))
import Oscillator_Engineering_Sensitivity_Analysis as osc  # noqa: E402


MOLECULES_PER_UM = 602.214076
DOUBLING_TIMES_MIN = (40.0, 50.0, 60.0)
MC_SAMPLES = 512
MC_SEED = 20260821


def median_metric(metrics: pd.DataFrame, column: str) -> float:
    if metrics.empty or column not in metrics:
        return np.nan
    values = metrics[column].replace([np.inf, -np.inf], np.nan).dropna()
    return float(values.median()) if len(values) else np.nan


def shutdown_metrics(frame: pd.DataFrame) -> dict:
    metrics = pulse_metrics(frame)
    return {
        "production_fwhm_h": median_metric(metrics, "production_fwhm_h"),
        "protein_fwhm_h": median_metric(metrics, "protein_fwhm_h"),
        "production_closure_ratio": median_metric(
            metrics, "production_closure_ratio"
        ),
        "protein_closure_ratio": median_metric(metrics, "protein_closure_ratio"),
        "production_peak_uM_h": median_metric(metrics, "production_peak"),
        "production_auc_uM": median_metric(metrics, "production_auc"),
    }


def shutdown_pass(metrics: dict, reference_peak: float) -> bool:
    return bool(
        np.isfinite(metrics["production_fwhm_h"])
        and metrics["production_fwhm_h"] > 0.0
        and np.isfinite(metrics["production_closure_ratio"])
        and metrics["production_closure_ratio"] < 0.20
        and np.isfinite(metrics["production_peak_uM_h"])
        and metrics["production_peak_uM_h"] >= 0.50 * reference_peak
    )


def generate_growth_waveform(doubling_time_min: float, output: Path) -> tuple[Path, dict]:
    """Re-solve the oscillator at one growth rate and save its C31 flux."""
    base = osc.make_reference_parameters()
    params = replace(base, Td_min=float(doubling_time_min))
    time_min = np.arange(0.0, 6000.0 + 0.5, 1.0)
    sol = solve_ivp(
        lambda t, y: osc.rhs(t, y, params),
        (float(time_min[0]), float(time_min[-1])),
        osc.initial_state(params),
        t_eval=time_min,
        method="LSODA",
        rtol=1e-7,
        atol=1e-9,
    )
    if not sol.success:
        raise RuntimeError(sol.message)

    flux_copies_min = params.c31_translation_per_mrna_per_min * sol.y[6]
    flux_uM_h = 60.0 * flux_copies_min / MOLECULES_PER_UM
    c31_uM = sol.y[7] / MOLECULES_PER_UM
    frame = pd.DataFrame(
        {
            "time_h": time_min / 60.0,
            "c31_uM": c31_uM,
            "translation_flux_uM_h": flux_uM_h,
        }
    )
    path = output / f"growth_{doubling_time_min:g}min_c31_waveform.csv"
    frame.to_csv(path, index=False)

    keep = time_min >= 1500.0
    local_t = time_min[keep] / 60.0
    local_flux = flux_uM_h[keep]
    dt_h = float(np.median(np.diff(local_t)))
    peaks, _ = find_peaks(
        local_flux,
        prominence=max(1e-9, 0.05 * float(np.ptp(local_flux))),
        distance=max(1, int(round(2.0 * doubling_time_min / 1.0))),
    )
    period_h = float(np.mean(np.diff(local_t[peaks][-6:])))
    widths = []
    for peak in peaks[-6:]:
        half = 0.5 * local_flux[peak]
        left = peak
        right = peak
        while left > 0 and local_flux[left] >= half:
            left -= 1
        while right + 1 < len(local_flux) and local_flux[right] >= half:
            right += 1
        widths.append((right - left - 1) * dt_h)
    info = {
        "doubling_time_min": doubling_time_min,
        "oscillator_period_h": period_h,
        "raw_flux_peak_uM_h": float(np.median(local_flux[peaks[-6:]])),
        "raw_flux_fwhm_h": float(np.median(widths)),
    }
    return path, info


def self_consistent_growth_analysis(output: Path) -> tuple[pd.DataFrame, dict]:
    base_shutdown = ShutdownConfig()
    no_shutdown = replace(
        base_shutdown,
        srna_max_tx_uM_h=0.0,
        copy_ratio_to_counter=0.0,
    )
    rows = []
    saved = {}
    for doubling_time_min in DOUBLING_TIMES_MIN:
        waveform_path, oscillator_info = generate_growth_waveform(
            doubling_time_min, output
        )
        k_dil = math.log(2.0) / (doubling_time_min / 60.0)
        counter = replace(CounterConfig(), k_dil_h=k_dil)
        sol, params, waveform = simulate_full(
            waveform_path,
            counter=counter,
            shutdown=base_shutdown,
            t_end_h=100.0,
            dt_h=0.02,
        )
        obs = observables(sol, params, base_shutdown)
        obs.to_csv(
            output / f"growth_{doubling_time_min:g}min_system_timecourse.csv",
            index=False,
        )
        _, counter_summary = cycle_switch_audit(obs, waveform)
        pulse = shutdown_metrics(obs)
        control_sol, control_params, _ = simulate_full(
            waveform_path,
            counter=counter,
            shutdown=no_shutdown,
            t_end_h=100.0,
            dt_h=0.02,
        )
        control = observables(control_sol, control_params, no_shutdown)
        control_peak = shutdown_metrics(control)["production_peak_uM_h"]
        row = {
            **oscillator_info,
            "growth_dilution_h": k_dil,
            **{f"counter_{key}": value for key, value in counter_summary.items()},
            **pulse,
            "no_shutdown_peak_uM_h": control_peak,
        }
        rows.append(row)
        saved[doubling_time_min] = (obs, waveform)

    result = pd.DataFrame(rows)
    result["shutdown_pass"] = result.apply(
        lambda row: shutdown_pass(
            row.to_dict(), float(row["no_shutdown_peak_uM_h"])
        ),
        axis=1,
    )
    result["joint_pass"] = (
        result["counter_system_pass"].astype(bool)
        & result["shutdown_pass"].astype(bool)
    )
    result.to_csv(output / "05_self_consistent_growth.csv", index=False)

    fig, axes = plt.subplots(3, 3, figsize=(14, 11), sharex="col")
    for row_index, doubling_time_min in enumerate(DOUBLING_TIMES_MIN):
        obs, waveform = saved[doubling_time_min]
        t = obs["time_h"].to_numpy(float)
        view = t >= 25.0
        axes[row_index, 0].plot(
            t[view], waveform.effective_source(t[view]), color="tab:purple"
        )
        axes[row_index, 1].plot(t[view], obs.loc[view, "LR_fraction"], color="tab:red")
        axes[row_index, 1].plot(t[view], obs.loc[view, "PB_fraction"], color="tab:blue")
        axes[row_index, 2].plot(
            t[view], obs.loc[view, "output_translation_uM_h"], color="tab:green"
        )
        axes[row_index, 0].set_ylabel(f"Td={doubling_time_min:g} min\nInt flux")
        for axis in axes[row_index]:
            axis.grid(alpha=0.25)
    axes[0, 0].set_title("RBS-scaled Integrase production")
    axes[0, 1].set_title("Counter state")
    axes[0, 2].set_title("Shutdown-limited output")
    for axis in axes[-1]:
        axis.set_xlabel("Time (h)")
    fig.suptitle("Self-consistent growth-rate analysis")
    fig.tight_layout()
    fig.savefig(output / "05_self_consistent_growth.png", dpi=170)
    plt.close(fig)

    summary = {
        "conditions": int(len(result)),
        "joint_pass_conditions": int(result["joint_pass"].sum()),
        "all_growth_conditions_pass": bool(result["joint_pass"].all()),
    }
    return result, summary


def nominal_driver() -> tuple[pd.DataFrame, dict, float]:
    waveform_csv = extract_waveform()
    off = ShutdownConfig(srna_max_tx_uM_h=0.0, copy_ratio_to_counter=0.0)
    sol, params, _ = simulate_full(waveform_csv, shutdown=off)
    driver = observables(sol, params, off)
    no_shutdown_metrics = shutdown_metrics(driver)
    return driver, params, no_shutdown_metrics["production_peak_uM_h"]


def oat_sensitivity_analysis(
    output: Path,
    driver: pd.DataFrame,
    params: dict,
    reference_peak: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    base = ShutdownConfig()
    value_grid = {
        "threshold_ratio_to_rdf": [0.5, 0.75, 1.0, 1.5, 2.0],
        "n_shutdown": [2.9, 3.1, 3.4, 3.6, 3.8],
        "promoter_leak": [0.004 / 0.6, 0.005 / 0.6, 0.009, 0.0105, 0.006 / 0.5],
        "srna_max_tx_uM_h": [0.006, 0.03, 0.06, 0.18, 0.6],
        "srna_decay_h": [0.6, 0.9, 1.2, 1.8, 2.4],
        "srna_cocodegradation_fraction": [0.5, 0.625, 0.75, 0.875, 1.0],
        "pair_rate_uM_h": [600.0, 900.0, 1200.0, 1800.0, 2400.0],
        "output_max_tx_uM_h": [0.006, 0.015, 0.03, 0.045, 0.06],
        "output_mrna_decay_h": [3.0, 4.5, 6.0, 9.0, 12.0],
        "copy_ratio_to_counter": [0.5, 0.75, 1.0, 1.5, 2.0],
    }
    rows = []
    for name, values in value_grid.items():
        baseline = float(getattr(base, name))
        for value in values:
            factor = float(value) / baseline
            cfg = replace(base, **{name: value})
            result = simulate_downstream(driver, params, cfg)
            metrics = shutdown_metrics(result)
            sample_reference_peak = (
                reference_peak
                * cfg.output_max_tx_uM_h
                / base.output_max_tx_uM_h
            )
            rows.append(
                {
                    "parameter": name,
                    "baseline": baseline,
                    "factor": factor,
                    "value": value,
                    **metrics,
                    "feasible": shutdown_pass(metrics, sample_reference_peak),
                }
            )
    scan = pd.DataFrame(rows)
    scan.to_csv(output / "06_shutdown_oat_sensitivity.csv", index=False)

    ranks = []
    for name, group in scan.groupby("parameter", sort=False):
        valid = group.replace([np.inf, -np.inf], np.nan).dropna(
            subset=["production_fwhm_h", "production_closure_ratio"]
        )
        width_rho = (
            float(spearmanr(valid["factor"], valid["production_fwhm_h"]).statistic)
            if valid["production_fwhm_h"].nunique() > 1
            else np.nan
        )
        ranks.append(
            {
                "parameter": name,
                "width_range_h": float(
                    valid["production_fwhm_h"].max()
                    - valid["production_fwhm_h"].min()
                ),
                "closure_range": float(
                    valid["production_closure_ratio"].max()
                    - valid["production_closure_ratio"].min()
                ),
                "protein_width_range_h": float(
                    valid["protein_fwhm_h"].max()
                    - valid["protein_fwhm_h"].min()
                ),
                "protein_closure_range": float(
                    valid["protein_closure_ratio"].max()
                    - valid["protein_closure_ratio"].min()
                ),
                "factor_to_width_spearman": width_rho,
                "feasible_fraction": float(valid["feasible"].mean()),
            }
        )
    ranking = pd.DataFrame(ranks).sort_values("width_range_h", ascending=False)
    ranking.to_csv(output / "06_shutdown_oat_ranking.csv", index=False)

    fig, axes = plt.subplots(2, 1, figsize=(10, 9), sharex=True)
    for name, group in scan.groupby("parameter", sort=False):
        axes[0].plot(
            group["factor"], group["production_fwhm_h"], marker="o", label=name
        )
        axes[1].plot(
            group["factor"],
            group["production_closure_ratio"],
            marker="o",
            label=name,
        )
    axes[0].set_ylabel("Production FWHM (h)")
    axes[1].set_ylabel("Late-LR output / peak")
    axes[1].set_xlabel("Parameter / baseline")
    axes[1].axhline(0.20, color="black", ls="--", lw=1)
    axes[0].legend(ncol=2, fontsize=8)
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.suptitle("Shutdown one-at-a-time sensitivity")
    fig.tight_layout()
    fig.savefig(output / "06_shutdown_oat_sensitivity.png", dpi=170)
    plt.close(fig)
    return scan, ranking


def literature_range_analysis(
    output: Path,
    driver: pd.DataFrame,
    params: dict,
    reference_peak: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    ranges = {
        "threshold_ratio_to_rdf": (0.5, 2.0, "log", "operator-variant design envelope"),
        "n_shutdown": (2.9, 3.8, "linear", "Shin 2020 BM3R1 B1-B3"),
        "promoter_leak": (0.004 / 0.6, 0.006 / 0.5, "linear", "Shin 2020 BM3R1 B1-B3 ymin/ymax"),
        "copy_ratio_to_counter": (0.5, 2.0, "log", "one-half to two-fold cassette dosage scenario"),
        "srna_max_tx_uM_h": (0.006, 0.6, "log", "Levine 2007 estimated RyhB range"),
        "srna_decay_h": (0.6, 2.4, "log", "factor-two context range around Levine 1.2/h"),
        "srna_cocodegradation_fraction": (0.5, 1.0, "linear", "Levine partial-to-full co-degradation cases"),
        "pair_rate_uM_h": (600.0, 2400.0, "log", "factor-two context range around Levine 1200/(uM h)"),
        "output_max_tx_uM_h": (
            0.006,
            0.06,
            "log",
            "project envelope below Levine's ~1 nM/min target scale",
        ),
        "output_mrna_decay_h": (3.0, 12.0, "log", "factor-two context range around Levine 6/h"),
    }
    pd.DataFrame(
        [
            {
                "parameter": name,
                "lower": spec[0],
                "upper": spec[1],
                "sampling": spec[2],
                "basis": spec[3],
            }
            for name, spec in ranges.items()
        ]
    ).to_csv(output / "07_parameter_ranges.csv", index=False)
    parameter_names = list(ranges)
    unit_samples = qmc.LatinHypercube(
        d=len(parameter_names), seed=MC_SEED
    ).random(n=MC_SAMPLES)
    rows = []
    for sample in range(MC_SAMPLES):
        values = {}
        for parameter_index, (name, (lower, upper, scale, _)) in enumerate(
            ranges.items()
        ):
            draw = float(unit_samples[sample, parameter_index])
            if scale == "log":
                values[name] = float(np.exp(np.log(lower) + draw * np.log(upper / lower)))
            else:
                values[name] = float(lower + draw * (upper - lower))
        cfg = replace(ShutdownConfig(), **values)
        result = simulate_downstream(driver, params, cfg)
        metrics = shutdown_metrics(result)
        sample_reference_peak = (
            reference_peak
            * values["output_max_tx_uM_h"]
            / ShutdownConfig().output_max_tx_uM_h
        )
        rows.append(
            {
                "sample": sample,
                **values,
                "effective_srna_to_target_tx_ratio": (
                    values["srna_max_tx_uM_h"]
                    * values["copy_ratio_to_counter"]
                    / values["output_max_tx_uM_h"]
                ),
                **metrics,
                "no_shutdown_reference_peak_uM_h": sample_reference_peak,
                "feasible": shutdown_pass(metrics, sample_reference_peak),
            }
        )
    samples = pd.DataFrame(rows)
    samples["closure_pass"] = samples["production_closure_ratio"] < 0.20
    samples["peak_pass"] = (
        samples["production_peak_uM_h"]
        >= 0.50 * samples["no_shutdown_reference_peak_uM_h"]
    )
    samples.to_csv(output / "07_literature_range_samples.csv", index=False)

    supply_bins = [0.0, 1.0, 2.0, 4.0, 8.0, np.inf]
    supply_labels = ["<1", "1-2", "2-4", "4-8", ">=8"]
    samples["supply_ratio_bin"] = pd.cut(
        samples["effective_srna_to_target_tx_ratio"],
        bins=supply_bins,
        labels=supply_labels,
        right=False,
    )
    supply_pass = (
        samples.groupby("supply_ratio_bin", observed=False)
        .agg(
            samples=("feasible", "size"),
            feasible_samples=("feasible", "sum"),
            combined_pass_fraction=("feasible", "mean"),
            closure_pass_fraction=("closure_pass", "mean"),
            peak_pass_fraction=("peak_pass", "mean"),
            median_closure=("production_closure_ratio", "median"),
        )
        .reset_index()
    )
    supply_pass.to_csv(output / "07_supply_ratio_pass_rate.csv", index=False)

    correlations = []
    for name in ranges:
        correlations.append(
            {
                "parameter": name,
                "rho_with_width": float(
                    spearmanr(np.log(samples[name]), samples["production_fwhm_h"]).statistic
                ),
                "rho_with_closure": float(
                    spearmanr(
                        np.log(samples[name]), samples["production_closure_ratio"]
                    ).statistic
                ),
                "rho_with_protein_closure": float(
                    spearmanr(
                        np.log(samples[name]), samples["protein_closure_ratio"]
                    ).statistic
                ),
            }
        )
    correlation_frame = pd.DataFrame(correlations)
    correlation_frame["max_abs_rho"] = correlation_frame[
        ["rho_with_width", "rho_with_closure", "rho_with_protein_closure"]
    ].abs().max(axis=1)
    correlation_frame = correlation_frame.sort_values(
        "max_abs_rho", ascending=False
    )
    correlation_frame.to_csv(output / "07_literature_range_correlations.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    axes[0].hist(samples["production_fwhm_h"], bins=22, color="tab:blue", alpha=0.8)
    axes[0].set_xlabel("Production FWHM (h)")
    axes[0].set_ylabel("Samples")
    plot_corr = correlation_frame.sort_values("rho_with_width")
    axes[1].barh(
        plot_corr["parameter"], plot_corr["rho_with_width"], color="tab:green"
    )
    axes[1].axvline(0.0, color="black", lw=0.8)
    axes[1].set_xlabel("Spearman rho with pulse width")
    positions = np.arange(len(supply_pass))
    axes[2].bar(
        positions - 0.18,
        supply_pass["closure_pass_fraction"],
        width=0.36,
        color="tab:green",
        label="closure <20%",
    )
    axes[2].bar(
        positions + 0.18,
        supply_pass["combined_pass_fraction"],
        width=0.36,
        color="tab:orange",
        label="closure + peak",
    )
    axes[2].set_xticks(positions)
    axes[2].set_xticklabels(supply_pass["supply_ratio_bin"].astype(str))
    axes[2].set_xlabel(r"Maximum effective supply ratio, $\alpha_S C/\alpha_M$")
    axes[2].set_ylabel("Fraction meeting criteria")
    axes[2].set_ylim(0.0, 1.0)
    axes[2].legend(fontsize=8)
    for axis in axes:
        axis.grid(alpha=0.22)
    fig.suptitle("Literature-anchored and declared-assumption parameter envelope")
    fig.tight_layout()
    fig.savefig(output / "07_shutdown_joint_uncertainty.png", dpi=170)
    plt.close(fig)

    width = samples["production_fwhm_h"].dropna()
    closure = samples["production_closure_ratio"].dropna()
    summary = {
        "samples": int(len(samples)),
        "sampling_method": (
            "Latin hypercube; transformed log-uniform for positive rates and "
            "linear-uniform for bounded dimensionless parameters"
        ),
        "criterion_pass_fraction": float(samples["feasible"].mean()),
        "width_q05_h": float(width.quantile(0.05)),
        "width_median_h": float(width.median()),
        "width_q95_h": float(width.quantile(0.95)),
        "closure_q95": float(closure.quantile(0.95)),
        "failed_samples": int((~samples["feasible"]).sum()),
        "median_effective_supply_ratio_feasible": float(
            samples.loc[
                samples["feasible"], "effective_srna_to_target_tx_ratio"
            ].median()
        ),
        "median_effective_supply_ratio_failed": float(
            samples.loc[
                ~samples["feasible"], "effective_srna_to_target_tx_ratio"
            ].median()
        ),
        "supply_ratio_bin_pass_rates": supply_pass.to_dict(orient="records"),
        "strongest_width_parameter": str(
            correlation_frame.iloc[
                correlation_frame["rho_with_width"].abs().argmax()
            ]["parameter"]
        ),
        "strongest_closure_parameter": str(
            correlation_frame.iloc[
                correlation_frame["rho_with_closure"].abs().argmax()
            ]["parameter"]
        ),
    }
    return samples, correlation_frame, summary


def main() -> None:
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output = HERE / "outputs" / f"shutdown_extended_{stamp}"
    output.mkdir(parents=True, exist_ok=False)

    growth, growth_summary = self_consistent_growth_analysis(output)
    driver, params, reference_peak = nominal_driver()
    _, oat_ranking = oat_sensitivity_analysis(
        output, driver, params, reference_peak
    )
    _, correlations, uncertainty_summary = literature_range_analysis(
        output, driver, params, reference_peak
    )

    summary = {
        "counter_model_source": "rdfmodel_new/model/zhao_core.py",
        "counter_input_source": (
            "rdfmodel_new/docs/week4_unloaded/"
            "unloaded_C31_translation_trajectory.csv"
        ),
        "self_consistent_growth": growth_summary,
        "growth_conditions": growth.to_dict(orient="records"),
        "literature_and_assumption_range_analysis": uncertainty_summary,
        "top_oat_width_parameters": oat_ranking.head(5).to_dict(orient="records"),
        "top_joint_correlations": correlations.head(5).to_dict(orient="records"),
        "acceptance_criteria": {
            "counter": "one complete crossing per cycle and >=95% alternation",
            "shutdown": (
                "finite production pulse; production late-LR output <20% of peak; "
                "production peak >=50% of no-shutdown control. These are declared "
                "engineering criteria, not literature constants."
            ),
        },
        "scope_boundary": (
            "The model is a single PB/LR bit with a local shutdown readout. "
            "A target-count decoder is not included because no target count or bit topology is specified."
        ),
        "base_counter": asdict(CounterConfig()),
        "base_shutdown": asdict(ShutdownConfig()),
    }
    (output / "extended_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True), encoding="utf-8"
    )
    report = f"""# Extended shutdown model summary

## Self-consistent growth analysis

- Growth conditions tested: **{growth_summary['conditions']}**
- Conditions passing counter and shutdown criteria: **{growth_summary['joint_pass_conditions']}**
- All 40/50/60 min conditions pass: **{growth_summary['all_growth_conditions_pass']}**

## Literature and declared-assumption range analysis

- Samples: **{uncertainty_summary['samples']}**
- Samples meeting the declared engineering criteria: **{uncertainty_summary['criterion_pass_fraction']:.1%}**
- Pulse-width 5/50/95% quantiles: **{uncertainty_summary['width_q05_h']:.3g} / {uncertainty_summary['width_median_h']:.3g} / {uncertainty_summary['width_q95_h']:.3g} h**
- 95th-percentile late-LR output: **{uncertainty_summary['closure_q95']:.1%} of peak**
- Samples failing the combined shutdown criteria: **{uncertainty_summary['failed_samples']}**
- Median effective supply ratio among passing samples: **{uncertainty_summary['median_effective_supply_ratio_feasible']:.3g}**
- Median effective supply ratio among failing samples: **{uncertainty_summary['median_effective_supply_ratio_failed']:.3g}**
- Strongest width-associated parameter: **{uncertainty_summary['strongest_width_parameter']}**
- Strongest closure-associated parameter: **{uncertainty_summary['strongest_closure_parameter']}**

The range analysis is conditional on the validated nominal counter trajectory.
It samples published RNA ranges plus explicitly declared design/context envelopes;
the pass fraction is not a wet-lab success probability.
The current scope is one PB/LR bit; a target-count decoder is not included.
"""
    (output / "EXTENDED_RESULTS.md").write_text(report, encoding="utf-8")
    (HERE / "latest_extended_run.txt").write_text(str(output), encoding="utf-8")
    print(json.dumps(summary, indent=2, allow_nan=True))
    print(f"Outputs: {output}")


if __name__ == "__main__":
    main()
