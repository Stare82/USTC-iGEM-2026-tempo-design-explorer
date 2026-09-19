"""Fast structural checks for the coupled shutdown model."""

from __future__ import annotations

from dataclasses import replace

import numpy as np

from extract_waveform import extract_waveform
from shutdown_core import (
    ShutdownConfig,
    CounterConfig,
    counter_params,
    cycle_switch_audit,
    observables,
    pulse_metrics,
    simulate_full,
    shutdown_rhs,
)


def main() -> None:
    waveform_csv = extract_waveform()
    shutdown = ShutdownConfig()
    sol, params, waveform = simulate_full(
        waveform_csv,
        shutdown=shutdown,
        t_end_h=70.0,
        dt_h=0.05,
    )
    obs = observables(sol, params, shutdown)
    assert sol.success
    assert np.isfinite(sol.y).all()
    assert sol.y.min() > -1e-6
    dna_total = obs["LR_fraction"] + obs["PB_fraction"]
    assert float(np.max(np.abs(dna_total - 1.0))) < 5e-3
    assert float(obs["shutdown_sRNA_uM"].max()) > 0.0
    assert float(obs["output_protein_uM"].max()) > 0.0
    peaks, troughs = waveform.peaks_and_troughs()
    assert len(peaks) >= 2 and len(troughs) >= 1
    metrics = pulse_metrics(obs, skip_first=0)
    assert not metrics.empty
    _, counter_audit = cycle_switch_audit(obs, waveform)
    assert counter_audit["system_pass"]

    # With production and first-order turnover disabled, the same pair-loss
    # flux must remove one mRNA and p sRNA molecules.
    pair_cfg = replace(
        shutdown,
        srna_max_tx_uM_h=0.0,
        output_max_tx_uM_h=0.0,
        srna_decay_h=0.0,
        output_mrna_decay_h=0.0,
        output_translation_h=0.0,
        output_protein_decay_h=0.0,
    )
    pair_state = np.array([0.001, 0.002, 0.0])
    pair_rates, _ = shutdown_rhs(
        0.0,
        0.0,
        pair_state,
        counter_params(CounterConfig()),
        pair_cfg,
    )
    expected_pair_loss = pair_cfg.pair_rate_uM_h * pair_state[0] * pair_state[1]
    assert np.isclose(pair_rates[0], -pair_cfg.srna_cocodegradation_fraction * expected_pair_loss)
    assert np.isclose(pair_rates[1], -expected_pair_loss)
    print("shutdown model verification passed")


if __name__ == "__main__":
    main()
