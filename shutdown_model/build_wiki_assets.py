"""Build web-oriented mechanism and literature-benchmark figures."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from shutdown_core import ShutdownConfig


HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "wiki_report" / "assets"


def paired_steady_state(alpha_m: float, cfg: ShutdownConfig) -> tuple[float, float, float]:
    """Return steady-state sRNA, mRNA, and pairing flux for constant production."""
    alpha_s = cfg.srna_max_tx_uM_h
    beta_s = cfg.srna_decay_h
    beta_m = cfg.output_mrna_decay_h
    p = cfg.srna_cocodegradation_fraction
    k = cfg.pair_rate_uM_h
    upper = min(alpha_m, alpha_s / p) if p > 0.0 else alpha_m

    def residual(pair_flux: float) -> float:
        srna = (alpha_s - p * pair_flux) / beta_s
        mrna = (alpha_m - pair_flux) / beta_m
        return pair_flux - k * srna * mrna

    pair_flux = brentq(residual, 0.0, upper)
    srna = (alpha_s - p * pair_flux) / beta_s
    mrna = (alpha_m - pair_flux) / beta_m
    return srna, mrna, pair_flux


def integrate_to_steady_state(alpha_m: float, cfg: ShutdownConfig) -> tuple[float, float]:
    """Integrate the reduced two-RNA system under constant full induction."""
    alpha_s = cfg.srna_max_tx_uM_h
    beta_s = cfg.srna_decay_h
    beta_m = cfg.output_mrna_decay_h
    p = cfg.srna_cocodegradation_fraction
    k = cfg.pair_rate_uM_h

    def rhs(_time: float, state: np.ndarray) -> np.ndarray:
        srna, mrna = np.maximum(state, 0.0)
        pair_loss = k * srna * mrna
        return np.array(
            [alpha_s - beta_s * srna - p * pair_loss, alpha_m - beta_m * mrna - pair_loss]
        )

    solution = solve_ivp(
        rhs,
        (0.0, 40.0),
        np.zeros(2),
        method="BDF",
        rtol=1e-9,
        atol=1e-12,
    )
    if not solution.success:
        raise RuntimeError(solution.message)
    return tuple(np.maximum(solution.y[:, -1], 0.0))


def draw_box(ax, x: float, y: float, width: float, height: float, text: str, color: str) -> None:
    box = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.015,rounding_size=0.02",
        linewidth=1.2,
        edgecolor=color,
        facecolor="white",
    )
    ax.add_patch(box)
    ax.text(x + width / 2, y + height / 2, text, ha="center", va="center", fontsize=9)


def draw_arrow(ax, start: tuple[float, float], end: tuple[float, float], color: str, inhibitory: bool = False) -> None:
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-[" if inhibitory else "-|>",
        mutation_scale=11,
        linewidth=1.5,
        color=color,
    )
    ax.add_patch(arrow)


def make_mechanism_figure() -> None:
    teal = "#0F6B67"
    coral = "#D4664A"
    ink = "#19313B"
    green = "#3A8D5D"
    gray = "#7A8588"

    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.9))
    stages = [
        ("PB state", "BM3R1 high", "Shutdown promoter OFF", "Target mRNA translated"),
        ("PB -> LR transition", "BM3R1 decays", "Repression is released", "RyhB accumulates"),
        ("Late LR state", "RyhB high", "Target mRNA is paired\nand co-degraded", "New translation OFF"),
    ]
    for index, (title, first, second, third) in enumerate(stages):
        ax = axes[index]
        ax.set_xlim(0, 1.12)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.set_title(title, fontsize=12, fontweight="bold", color=ink, pad=10)
        draw_box(ax, 0.19, 0.75, 0.62, 0.13, first, teal if index < 2 else green)
        draw_box(ax, 0.19, 0.48, 0.62, 0.13, second, coral if index == 0 else teal)
        draw_box(ax, 0.19, 0.20, 0.62, 0.15, third, green if index == 0 else coral)
        if index == 0:
            draw_arrow(ax, (0.5, 0.75), (0.5, 0.62), coral, inhibitory=True)
            draw_arrow(ax, (0.5, 0.48), (0.5, 0.36), gray)
        else:
            draw_arrow(ax, (0.5, 0.75), (0.5, 0.62), teal)
            draw_arrow(ax, (0.5, 0.48), (0.5, 0.36), coral)
        if index < 2:
            draw_arrow(ax, (0.87, 0.5), (1.08, 0.5), ink)

    fig.suptitle(
        "BM3R1-gated sRNA shutdown: the counter state starts expression, RNA pairing ends it",
        fontsize=13,
        fontweight="bold",
        color=ink,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(OUTPUT_DIR / "00_shutdown_mechanism.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def make_levine_benchmark() -> dict:
    cfg = ShutdownConfig()
    alpha_m_values = np.linspace(0.003, 0.18, 120)
    rows = []
    for alpha_m in alpha_m_values:
        srna_exact, mrna_exact, pair_flux = paired_steady_state(alpha_m, cfg)
        srna_ode, mrna_ode = integrate_to_steady_state(alpha_m, cfg)
        rows.append(
            {
                "alpha_m_uM_h": alpha_m,
                "alpha_m_nM_min": alpha_m * 1000.0 / 60.0,
                "srna_exact_uM": srna_exact,
                "mrna_exact_uM": mrna_exact,
                "srna_ode_uM": srna_ode,
                "mrna_ode_uM": mrna_ode,
                "pair_flux_uM_h": pair_flux,
                "mrna_no_srna_uM": alpha_m / cfg.output_mrna_decay_h,
                "mrna_ideal_threshold_uM": max(
                    0.0,
                    (alpha_m - cfg.srna_max_tx_uM_h / cfg.srna_cocodegradation_fraction)
                    / cfg.output_mrna_decay_h,
                ),
            }
        )
    frame = pd.DataFrame(rows)
    frame.to_csv(OUTPUT_DIR / "01_levine_threshold_validation.csv", index=False)

    scale = 1000.0
    max_abs_error = float(np.max(np.abs(frame["mrna_ode_uM"] - frame["mrna_exact_uM"])))
    max_scaled_error = float(
        np.max(
            np.abs(frame["mrna_ode_uM"] - frame["mrna_exact_uM"])
            / np.maximum(frame["mrna_exact_uM"], 1e-9)
        )
    )

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2))
    axes[0].plot(
        frame["alpha_m_nM_min"],
        scale * frame["mrna_no_srna_uM"],
        color="#7A8588",
        linestyle="--",
        linewidth=1.8,
        label="No sRNA",
    )
    axes[0].plot(
        frame["alpha_m_nM_min"],
        scale * frame["mrna_ideal_threshold_uM"],
        color="#D4664A",
        linestyle=":",
        linewidth=2.0,
        label="Ideal one-to-one threshold",
    )
    axes[0].plot(
        frame["alpha_m_nM_min"],
        scale * frame["mrna_exact_uM"],
        color="#0F6B67",
        linewidth=2.4,
        label="Finite-rate Levine model",
    )
    axes[0].axvline(1.0, color="#19313B", alpha=0.45, linewidth=1.0)
    axes[0].set_xlabel(r"Target transcription, $\alpha_M$ (nM min$^{-1}$)")
    axes[0].set_ylabel("Steady-state target mRNA (nM)")
    axes[0].set_title("Published threshold behavior is reproduced")
    axes[0].legend(frameon=False, fontsize=8)

    exact_nm = scale * frame["mrna_exact_uM"]
    ode_nm = scale * frame["mrna_ode_uM"]
    axes[1].scatter(exact_nm, ode_nm, s=13, color="#3A8D5D", alpha=0.8)
    upper = float(max(exact_nm.max(), ode_nm.max()))
    axes[1].plot([0, upper], [0, upper], color="#D4664A", linestyle="--", linewidth=1.5)
    axes[1].set_xlabel("Analytical steady state (nM)")
    axes[1].set_ylabel("Integrated ODE at 40 h (nM)")
    axes[1].set_title("Numerical implementation matches the steady state")
    axes[1].text(
        0.04,
        0.95,
        f"max absolute error = {max_abs_error * scale:.2e} nM",
        transform=axes[1].transAxes,
        va="top",
        fontsize=8,
        color="#19313B",
    )
    for ax in axes:
        ax.grid(alpha=0.2)
    fig.suptitle(
        r"RyhB/sodB benchmark ($\alpha_S=1$ nM min$^{-1}$; Levine et al., 2007)",
        fontsize=12.5,
        fontweight="bold",
        color="#19313B",
    )
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig(OUTPUT_DIR / "01_levine_threshold_validation.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    summary = {
        "benchmark": "Levine 2007 finite-rate RyhB/sodB steady-state response",
        "alpha_s_nM_min": cfg.srna_max_tx_uM_h * 1000.0 / 60.0,
        "alpha_m_min_nM_min": float(frame["alpha_m_nM_min"].min()),
        "alpha_m_max_nM_min": float(frame["alpha_m_nM_min"].max()),
        "max_absolute_mrna_error_nM": max_abs_error * scale,
        "max_relative_mrna_error_with_1e-6_uM_floor": max_scaled_error,
        "interpretation": (
            "Implementation benchmark only: agreement with the analytical steady state "
            "does not constitute independent wet-lab validation."
        ),
    }
    (OUTPUT_DIR / "01_levine_threshold_validation_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return summary


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    make_mechanism_figure()
    summary = make_levine_benchmark()
    print(json.dumps(summary, indent=2))
    print(f"Outputs: {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
