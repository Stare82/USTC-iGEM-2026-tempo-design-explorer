"""Adapter from browser design parameters to the deterministic TEMPO models.

This module deliberately keeps the existing scientific source files unchanged.
It imports the oscillator, counter, and shutdown implementations and translates
their named observables into the stable JSON contract used by the browser.
"""

from __future__ import annotations

from dataclasses import replace
import math
from pathlib import Path
import sys
import tempfile
from typing import Any

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks


PROJECT_ROOT = Path(__file__).resolve().parents[2]
OSCILLATOR_CODE = PROJECT_ROOT / "oscillator" / "code"
SHUTDOWN_DIR = PROJECT_ROOT / "shutdown_model"

for path in (OSCILLATOR_CODE, SHUTDOWN_DIR):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import Mechanistic_ODE_Global_Sensitivity_Analysis as oscillator  # noqa: E402
import Shared_PLtetO1_Period_Knob_Design_Map as period_map  # noqa: E402
import shutdown_core  # noqa: E402


MODEL_VERSION = "tempo-certified-period-knob-ode-2026.09"
MAX_HOURS = 84.0
PERIOD_KNOBS = {
    1: {"id": "K1", "rbs": 0.60, "mrna_half_life_min": 1.0,
        "median_generations": 7.25458, "nominal_period_h": 6.0455, "passing_tags": ()},
    2: {"id": "K2", "rbs": 0.80, "mrna_half_life_min": 1.0,
        "median_generations": 8.31446, "nominal_period_h": 6.9287, "passing_tags": ()},
    3: {"id": "K3", "rbs": 1.20, "mrna_half_life_min": 1.0,
        "median_generations": 9.76918, "nominal_period_h": 8.1408, "passing_tags": ()},
    4: {"id": "K4", "rbs": 1.20, "mrna_half_life_min": 1.5,
        "median_generations": 11.2208, "nominal_period_h": 9.3507, "passing_tags": (8.0, 12.0)},
    5: {"id": "K5", "rbs": 0.60, "mrna_half_life_min": 4.0,
        "median_generations": 12.5416, "nominal_period_h": 10.4513, "passing_tags": (12.0,)},
}


def _clamp(value: float, lower: float, upper: float) -> float:
    return min(upper, max(lower, float(value)))


def _validated_design(raw: dict[str, Any]) -> dict[str, float]:
    """Validate and bound public inputs to ranges covered by the UI."""
    return {
        "hours": _clamp(raw.get("hours", 60.0), 36.0, MAX_HOURS),
        "period": float(round(_clamp(raw.get("period", 4.0), 1.0, 5.0))),
        "pulseWidth": _clamp(raw.get("pulseWidth", 3.6), 1.5, 5.5),
        "promoter": _clamp(raw.get("promoter", 1.0), 0.5, 1.5),
        "oscRbs": _clamp(raw.get("oscRbs", 1.0), 0.5, 1.5),
        "c31": _clamp(raw.get("c31", 0.45), 0.3, 0.8),
        "degradation": _clamp(raw.get("degradation", 8.0), 2.8, 28.0),
        "rdf": _clamp(raw.get("rdf", 1.0), 0.4, 1.7),
        "bm3r1": _clamp(raw.get("bm3r1", 1.0), 0.4, 1.8),
        "srna": _clamp(raw.get("srna", 0.06), 0.02, 0.1),
    }


def _oscillator_parameters(design: dict[str, float]):
    knob = PERIOD_KNOBS[int(design["period"])]
    # The existing pulse-width target is implemented through the report's
    # independent C31-mRNA stability knob.  At 3.6 h it preserves the 2 min
    # frozen reference half-life used by the A->B certification.
    c31_half_life_min = _clamp(2.0 * design["pulseWidth"] / 3.6, 0.5, 5.0)
    settings = period_map.DesignSettings(
        mode="full",
        c31_RBS_scale=design["c31"],
        c31_mRNA_total_half_life_min=c31_half_life_min,
        peak_load_fraction=0.0,
    )
    vector, physical = period_map.build_vector(
        {
            "oscillator_translation_scale": knob["rbs"],
            "oscillator_mrna_total_half_life_min": knob["mrna_half_life_min"],
        },
        settings,
    )
    base = oscillator.make_reference_parameters()
    params, reconstructed = oscillator.apply_parameter_vector(vector, base)
    # Advanced controls are explicit multipliers around a certified K setting.
    # Shared PLtetO1 transcription is preserved for oscillator and C31.
    params = replace(
        params,
        oscillator_tx_per_plasmid_per_min=(
            params.oscillator_tx_per_plasmid_per_min * design["promoter"]
        ),
        c31_tx_per_plasmid_per_min=(
            params.c31_tx_per_plasmid_per_min * design["promoter"]
        ),
        oscillator_translation_per_mrna_per_min=(
            params.oscillator_translation_per_mrna_per_min * design["oscRbs"]
        ),
    )
    return params, reconstructed, knob, c31_half_life_min


def _median_fwhm_h(t_h: np.ndarray, signal: np.ndarray, peaks: np.ndarray) -> float:
    widths: list[float] = []
    for peak in peaks[-6:]:
        level = 0.5 * float(signal[peak])
        left = int(peak)
        right = int(peak)
        while left > 0 and signal[left] >= level:
            left -= 1
        while right < len(signal) - 1 and signal[right] >= level:
            right += 1
        if left > 0 and right < len(signal) - 1:
            widths.append(float(t_h[right] - t_h[left]))
    return float(np.median(widths)) if widths else float("nan")


def _solve_oscillator(design: dict[str, float]) -> tuple[pd.DataFrame, dict[str, Any]]:
    params, physical, knob, c31_half_life_min = _oscillator_parameters(design)
    duration_h = max(100.0, design["hours"])
    time_min = np.arange(0.0, duration_h * 60.0 + 0.5, 1.0)
    solution = solve_ivp(
        lambda t, y: oscillator.rhs(t, y, params, 1.0),
        (float(time_min[0]), float(time_min[-1])),
        oscillator.initial_state(params),
        t_eval=time_min,
        method="LSODA",
        rtol=1e-7,
        atol=1e-9,
    )
    if not solution.success:
        raise RuntimeError(f"Oscillator solve failed: {solution.message}")

    time_h = solution.t / 60.0
    c31_copies = np.maximum(0.0, solution.y[7])
    resource = np.array(
        [oscillator.resource_availability(value, params, 1.0) for value in solution.y[6]]
    )
    flux_copies_min = resource * params.c31_translation_per_mrna_per_min * np.maximum(
        0.0, solution.y[6]
    )
    sim_settings = oscillator.MODE_SETTINGS["full"]
    flux_uM_h = (
        60.0
        * flux_copies_min
        * sim_settings.active_C31_fraction
        / (
            sim_settings.molecules_per_uM_per_fL
            * sim_settings.cell_volume_fL
        )
    )
    c31_uM = c31_copies / (
        sim_settings.molecules_per_uM_per_fL * sim_settings.cell_volume_fL
    )

    span = float(np.ptp(flux_uM_h))
    peaks, _ = find_peaks(
        flux_uM_h,
        prominence=max(0.05 * span, 1e-9),
        distance=max(1, int(round(0.55 * knob["nominal_period_h"] * 60.0))),
    )
    stable_peaks = peaks[time_h[peaks] >= 30.0]
    period_h = (
        float(np.mean(np.diff(time_h[stable_peaks])))
        if len(stable_peaks) >= 2
        else float("nan")
    )
    width_h = _median_fwhm_h(time_h, flux_uM_h, stable_peaks)

    waveform = pd.DataFrame(
        {
            "time_h": time_h,
            "c31_uM": c31_uM,
            "translation_flux_uM_h": flux_uM_h,
        }
    )
    base_conditions = all(
        math.isclose(design[key], reference, rel_tol=0.0, abs_tol=1e-9)
        for key, reference in {
            "pulseWidth": 3.6,
            "promoter": 1.0,
            "oscRbs": 1.0,
            "c31": 0.45,
            "rdf": 1.0,
            "bm3r1": 1.0,
        }.items()
    )
    tag_certified = any(
        math.isclose(design["degradation"], tag, rel_tol=0.0, abs_tol=1e-9)
        for tag in knob["passing_tags"]
    )
    metadata = {
        "period_setting_id": knob["id"],
        "target_period_h": knob["nominal_period_h"],
        "realized_period_h": period_h,
        "target_pulse_width_h": design["pulseWidth"],
        "realized_pulse_width_h": width_h,
        "doubling_time_min": params.Td_min,
        "oscillator_promoter_scale": design["promoter"],
        "oscillator_rbs_scale": knob["rbs"] * design["oscRbs"],
        "period_knob_rbs_scale": knob["rbs"],
        "period_knob_oscillator_mrna_half_life_min": knob["mrna_half_life_min"],
        "period_knob_median_generations": knob["median_generations"],
        "c31_mrna_total_half_life_min": c31_half_life_min,
        "period_knob_passing_tags_h_inv": list(knob["passing_tags"]),
        "period_knob_base_conditions": base_conditions,
        "period_knob_ab_certified_for_tag": bool(base_conditions and tag_certified),
        "stochastic_ab_certified": False,
        "period_map_physical_inputs": physical,
    }
    return waveform, metadata


def _solve_coupled(
    design: dict[str, float], waveform: pd.DataFrame
) -> tuple[pd.DataFrame, Any, Any, Any]:
    counter = shutdown_core.CounterConfig(
        krep_tsl_h=15.0 * design["bm3r1"],
        krdf_tsl_h=200.0 * design["rdf"],
        k_tag_int_h=design["degradation"],
        # The A waveform already contains the selected physical C31 RBS.
        c31_rbs_scale=1.0,
        k_dil_h=math.log(2.0) / (50.0 / 60.0),
    )
    shutdown = shutdown_core.ShutdownConfig(srna_max_tx_uM_h=design["srna"])

    with tempfile.TemporaryDirectory(prefix="tempo-waveform-") as temp_dir:
        waveform_path = Path(temp_dir) / "waveform.csv"
        waveform.to_csv(waveform_path, index=False)
        solution, params, waveform_object = shutdown_core.simulate_full(
            waveform_path,
            counter=counter,
            shutdown=shutdown,
            # The existing audit discards the first 25 h and needs several
            # post-transient cycles to evaluate alternation. Always solve the
            # audited 100 h horizon, then crop the browser payload below.
            t_end_h=max(100.0, design["hours"]),
            dt_h=0.04,
            counter_initial_state="PB",
        )

    frame = shutdown_core.observables(solution, params, shutdown)
    return frame, params, waveform_object, shutdown


def _upcrossings(time: np.ndarray, signal: np.ndarray, threshold: float = 0.5) -> np.ndarray:
    indices = np.flatnonzero((signal[:-1] < threshold) & (signal[1:] >= threshold)) + 1
    return time[indices]


def _safe_number(value: Any, fallback: float = 0.0) -> float:
    number = float(value)
    return number if np.isfinite(number) else fallback


def _build_evaluation(
    frame: pd.DataFrame,
    waveform: Any,
    shutdown: Any,
    oscillator_metadata: dict[str, Any],
) -> tuple[dict[str, Any], dict[str, int]]:
    cycle_rows, cycle_summary = shutdown_core.cycle_switch_audit(frame, waveform)
    pulse_rows = shutdown_core.pulse_metrics(frame)

    time = frame["time_h"].to_numpy(float)
    lr = frame["LR_fraction"].to_numpy(float)
    output = frame["output_protein_uM"].to_numpy(float)
    integrase = frame["integrase_total_uM"].to_numpy(float)
    peaks, troughs = waveform.peaks_and_troughs()
    peaks = peaks[peaks <= time[-1]]
    troughs = troughs[troughs <= time[-1]]

    one_complete = _safe_number(cycle_summary.get("fraction_one_and_complete"), 0.0)
    alternation = _safe_number(cycle_summary.get("alternation_fidelity"), 0.0)
    flip_success = 100.0 * one_complete
    stable = bool(cycle_summary.get("system_pass", False))

    if len(troughs) and float(np.max(integrase)) > 0:
        trough_integrase = np.interp(troughs, time, integrase)
        recovery_ratio = float(np.max(trough_integrase) / np.max(integrase))
    else:
        recovery_ratio = 1.0
    recovery = recovery_ratio < 0.2

    peak_output = float(np.max(output)) if len(output) else 0.0
    off_mask = lr < 0.1
    leakage_fraction = (
        float(np.mean(output[off_mask]) / peak_output)
        if peak_output > 0 and np.any(off_mask)
        else 0.0
    )
    leakage_percent = 100.0 * leakage_fraction
    leakage = "Low" if leakage_percent < 10 else "Medium" if leakage_percent < 20 else "High"

    starts = _upcrossings(time, lr)
    output_start = float(starts[0]) if len(starts) else 0.0
    if not pulse_rows.empty:
        duration = _safe_number(pulse_rows["protein_fwhm_h"].median(), 0.0)
        residual = 100.0 * _safe_number(
            pulse_rows["production_closure_ratio"].median(), 1.0
        )
    else:
        duration = 0.0
        residual = 100.0
    shutdown_time = output_start + duration

    score = int(round(
        100.0
        * (
            0.55 * one_complete
            + 0.2 * alternation
            + 0.12 * (1.0 if recovery else max(0.0, 1.0 - recovery_ratio))
            + 0.13 * max(0.0, 1.0 - min(residual, 100.0) / 100.0)
        )
    ))
    score = int(_clamp(score, 0, 100))

    insights: list[dict[str, str]] = []
    knob_id = str(oscillator_metadata.get("period_setting_id", "K?"))
    passing_tags = oscillator_metadata.get("period_knob_passing_tags_h_inv", [])
    if oscillator_metadata.get("period_knob_ab_certified_for_tag", False):
        insights.append(
            {
                "level": "ok",
                "text": (
                    f"{knob_id} and the selected tag match the deterministic A->B certification."
                ),
            }
        )
    elif not oscillator_metadata.get("period_knob_base_conditions", False):
        insights.append(
            {
                "level": "warn",
                "text": (
                    f"{knob_id} is being solved outside the frozen certification settings."
                ),
            }
        )
    elif passing_tags:
        tag_text = ", ".join(f"{value:g}" for value in passing_tags)
        insights.append(
            {
                "level": "warn",
                "text": (
                    f"{knob_id} needs Integrase tag {tag_text} h^-1 for A->B certification."
                ),
            }
        )
    else:
        insights.append(
            {
                "level": "warn",
                "text": (
                    f"{knob_id} is oscillator-robust but not B-interface certified."
                ),
            }
        )
    if one_complete < 0.95:
        insights.append(
            {
                "level": "danger" if one_complete < 0.7 else "warn",
                "text": f"One-flip fidelity is {one_complete:.0%} across audited cycles.",
            }
        )
    if alternation < 0.95:
        insights.append(
            {
                "level": "warn",
                "text": f"Alternation fidelity is {alternation:.0%}; consecutive flips need attention.",
            }
        )
    if not recovery:
        insights.append(
            {
                "level": "warn",
                "text": f"Integrase carry-over is {recovery_ratio:.0%} at the inter-pulse trough.",
            }
        )
    if residual > 20:
        insights.append(
            {
                "level": "warn",
                "text": f"Residual output is {residual:.0f}%; the target is 20% or less.",
            }
        )
    if not insights:
        insights.append(
            {
                "level": "ok",
                "text": "The coupled ODE satisfies the current one-flip, alternation, recovery, and shutdown checks.",
            }
        )

    verdict = "success" if stable and recovery and residual <= 20 else "warning" if score >= 65 else "danger"
    title = {
        "success": "Within deterministic window",
        "warning": "Narrow deterministic margin",
        "danger": "Outside deterministic window",
    }[verdict]
    message = {
        "success": "The coupled ODE predicts complete alternating flips and a closing output window.",
        "warning": "The ODE solution completes part of the chain, but at least one audit has limited margin.",
        "danger": "The ODE solution fails at least one core counting or shutdown criterion.",
    }[verdict]

    evaluation = {
        "score": score,
        "verdict": verdict,
        "title": title,
        "message": message,
        "flipSuccess": flip_success,
        "stable": stable,
        "recovery": recovery,
        "recoveryMargin": 1.0 - recovery_ratio,
        "leakagePercent": leakage_percent,
        "leakage": leakage,
        "outputStart": output_start,
        "outputDuration": duration,
        "shutdownTime": shutdown_time,
        "residual": residual,
        "normalizedDose": 0.0,
        "insights": insights,
    }
    summary = {
        "pulses": int(len(peaks)),
        "flips": int(cycle_rows["one_and_complete"].sum()) if not cycle_rows.empty else 0,
        "outputs": int(len(starts)),
    }
    oscillator_metadata["audited_cycles"] = int(cycle_summary.get("audited_cycles", 0))
    oscillator_metadata["alternation_fidelity"] = alternation
    oscillator_metadata["recovery_ratio"] = recovery_ratio
    oscillator_metadata["shutdown_parameter_source"] = "shutdown_model.shutdown_core"
    return evaluation, summary


def _downsample_indices(length: int, limit: int = 900) -> np.ndarray:
    if length <= limit:
        return np.arange(length, dtype=int)
    return np.unique(np.linspace(0, length - 1, limit).round().astype(int))


def simulate_design(raw_design: dict[str, Any]) -> dict[str, Any]:
    design = _validated_design(raw_design)
    waveform_frame, oscillator_metadata = _solve_oscillator(design)
    coupled_frame, _, waveform, shutdown = _solve_coupled(design, waveform_frame)
    evaluation, summary = _build_evaluation(
        coupled_frame, waveform, shutdown, oscillator_metadata
    )

    visible_frame = coupled_frame[coupled_frame["time_h"] <= design["hours"]].copy()
    time = visible_frame["time_h"].to_numpy(float)
    indices = _downsample_indices(len(time))
    effective_input = waveform.effective_source(time)
    samples = [
        {
            "t": float(time[index]),
            "c31": float(effective_input[index]),
            "integrase": float(visible_frame["integrase_total_uM"].iloc[index]),
            "rdf": float(visible_frame["RDF_total_uM"].iloc[index]),
            "lr": float(visible_frame["LR_fraction"].iloc[index]),
            "pb": float(visible_frame["PB_fraction"].iloc[index]),
            "bm3r1": float(visible_frame["BM3R1_total_uM"].iloc[index]),
            "srna": float(visible_frame["shutdown_sRNA_uM"].iloc[index]),
            "output": float(visible_frame["output_protein_uM"].iloc[index]),
        }
        for index in indices
    ]

    peaks, _ = waveform.peaks_and_troughs()
    centers = [float(value) for value in peaks if value <= design["hours"]]
    visible_lr = visible_frame["LR_fraction"].to_numpy(float)
    visible_crossings = int(np.sum((visible_lr[:-1] - 0.5) * (visible_lr[1:] - 0.5) < 0))
    visible_outputs = int(len(_upcrossings(time, visible_lr)))
    summary = {
        "pulses": len(centers),
        "flips": visible_crossings,
        "outputs": visible_outputs,
    }
    return {
        "samples": samples,
        "centers": centers,
        "evaluation": evaluation,
        "summary": summary,
        "metadata": {
            **oscillator_metadata,
            "source": "deterministic_ode",
            "model_version": MODEL_VERSION,
            "solver": "SciPy LSODA",
            "states": {"oscillator": 8, "counter": 38, "shutdown": 3},
            "design": design,
        },
    }


__all__ = ["MODEL_VERSION", "simulate_design"]
