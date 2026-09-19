"""Audit the rdfmodel_new C31-flux/counter interface from both memory states."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from extract_waveform import extract_waveform
from shutdown_core import (
    CounterConfig,
    ShutdownConfig,
    cycle_switch_audit,
    observables,
    simulate_full,
    toggle_score,
)


HERE = Path(__file__).resolve().parent
RBS_VALUES = [0.30, 0.45, 0.60]
TAG_VALUES = [3.0, 6.0, 8.0, 12.0, 18.0, 22.0, 26.0]


def main() -> None:
    waveform_csv = extract_waveform()
    off = ShutdownConfig(srna_max_tx_uM_h=0.0, copy_ratio_to_counter=0.0)
    rows = []
    for rbs_scale in RBS_VALUES:
        for tag_rate in TAG_VALUES:
            for initial_state in ("PB", "LR"):
                counter = replace(
                    CounterConfig(),
                    c31_rbs_scale=rbs_scale,
                    k_tag_int_h=tag_rate,
                )
                sol, params, waveform = simulate_full(
                    waveform_csv,
                    counter=counter,
                    shutdown=off,
                    t_end_h=100.0,
                    dt_h=0.05,
                    counter_initial_state=initial_state,
                )
                obs = observables(sol, params, off)
                _, troughs = waveform.peaks_and_troughs()
                score = toggle_score(
                    obs["time_h"].to_numpy(),
                    obs["LR_fraction"].to_numpy(),
                    troughs,
                )
                _, cycle = cycle_switch_audit(obs, waveform)
                rows.append(
                    {
                        "C31_RBS_scale": rbs_scale,
                        "Int_tag_rate_h_inv": tag_rate,
                        "initial_state": initial_state,
                        "P_score": score["score"],
                        "P_high": score["high"],
                        "P_low": score["low"],
                        "P_fidelity": score["fidelity"],
                        **cycle,
                    }
                )
                print(rbs_scale, tag_rate, initial_state, score, cycle["system_pass"])

    detail = pd.DataFrame(rows)
    summary = (
        detail.groupby(["C31_RBS_scale", "Int_tag_rate_h_inv"], as_index=False)
        .agg(
            min_dual_initial_P_score=("P_score", "min"),
            min_one_and_complete=("fraction_one_and_complete", "min"),
            min_alternation_fidelity=("alternation_fidelity", "min"),
            both_initial_states_cycle_pass=("system_pass", "all"),
        )
    )
    summary["new_counter_pass"] = (
        (summary["min_dual_initial_P_score"] >= 0.95)
        & summary["both_initial_states_cycle_pass"]
    )

    data_dir = HERE / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    detail.to_csv(data_dir / "interface_scan.csv", index=False)
    summary.to_csv(data_dir / "interface_scan_summary.csv", index=False)

    matrix = summary.pivot(
        index="C31_RBS_scale",
        columns="Int_tag_rate_h_inv",
        values="min_dual_initial_P_score",
    )
    fig, ax = plt.subplots(figsize=(9.2, 4.2))
    im = ax.imshow(
        matrix.to_numpy(), vmin=0.0, vmax=1.0, cmap="viridis", aspect="auto"
    )
    ax.set_xticks(np.arange(len(matrix.columns)), [f"{x:g}" for x in matrix.columns])
    ax.set_yticks(np.arange(len(matrix.index)), [f"{x:.2f}" for x in matrix.index])
    ax.set_xlabel("Integrase tag loss rate (h$^{-1}$)")
    ax.set_ylabel("C31 RBS scale")
    ax.set_title("rdfmodel_new dual-initial-state counter score")
    for i in range(len(matrix.index)):
        for j in range(len(matrix.columns)):
            value = matrix.iloc[i, j]
            ax.text(
                j,
                i,
                f"{value:.2f}",
                ha="center",
                va="center",
                color="white" if value < 0.55 else "black",
            )
    fig.colorbar(im, ax=ax, label="minimum PB/LR-initial P score")
    fig.tight_layout()
    fig.savefig(data_dir / "interface_scan_heatmap.png", dpi=180)
    plt.close(fig)

    print(summary.to_string(index=False))
    print(data_dir / "interface_scan_summary.csv")


if __name__ == "__main__":
    main()
