"""Coupled oscillator -> BM3R1/RDF counter -> sRNA shutdown model.

The existing 38-state Zhao recombination model is used unchanged. Three
downstream states are appended:
    y[38] shutdown sRNA
    y[39] output mRNA
    y[40] output protein

The shutdown branch reads free BM3R1 through a repressed promoter. Its sRNA
and target mRNA are co-degraded with the reduced mass-action model of Levine
et al. (PLoS Biology, 2007), converting an LR state into a finite expression
pulse. Literature-derived and design-assumption parameters are distinguished
explicitly by parameter_provenance().
"""

from __future__ import annotations

import sys
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq
from scipy.signal import find_peaks


ROOT = Path(__file__).resolve().parents[1]
RDF_DIR = ROOT / "rdfmodel_new" / "model"
if str(RDF_DIR) not in sys.path:
    sys.path.insert(0, str(RDF_DIR))
import zhao_core as zc  # noqa: E402


@dataclass(frozen=True)
class CounterConfig:
    """Validated rdfmodel_new operating point for the A/B interface."""

    K_rep_uM: float = 0.0186
    n_rep: float = 3.4
    krep_tsl_h: float = 15.0
    krdf_tsl_h: float = 200.0
    k_tag_int_h: float = 8.0
    c31_rbs_scale: float = 0.45
    k_dil_h: float = 0.8317766


@dataclass(frozen=True)
class ShutdownConfig:
    """Literature-anchored RyhB/sodB proxy for the proposed shutdown branch."""

    # A ratio of one reuses the concentration threshold calibrated for the
    # counter's BM3R1-controlled RDF promoter.
    threshold_ratio_to_rdf: float = 1.0
    n_shutdown: float = 3.1
    promoter_leak: float = 0.005 / 0.6
    copy_ratio_to_counter: float = 1.0
    operator_load_stoich: float = 2.0
    srna_max_tx_uM_h: float = 0.06
    srna_decay_h: float = 1.2
    srna_cocodegradation_fraction: float = 1.0
    pair_rate_uM_h: float = 1200.0
    output_max_tx_uM_h: float = 0.03
    output_mrna_decay_h: float = 6.0
    output_translation_h: float = 1.0
    output_protein_decay_h: float = 0.0


def parameter_provenance(
    counter: CounterConfig | None = None,
    shutdown: ShutdownConfig | None = None,
) -> list[dict]:
    """Return a machine-readable audit trail for every exposed parameter."""
    counter = counter or CounterConfig()
    shutdown = shutdown or ShutdownConfig()
    levine = "Levine et al. 2007, doi:10.1371/journal.pbio.0050229"
    masse = "Masse et al. 2003, doi:10.1101/gad.1127103"
    moll = "Moll et al. 2003, doi:10.1261/rna.5850703"
    zhang = "Zhang et al. 2005, doi:10.1128/JB.187.3.980-990.2005"
    selinger = "Selinger et al. 2000, doi:10.1038/82367"
    shin = "Shin et al. 2020, doi:10.15252/msb.20199401"
    zhao = "Zhao et al. 2019, doi:10.1093/nar/gkz245"
    return [
        {
            "parameter": "threshold_ratio_to_rdf",
            "value": shutdown.threshold_ratio_to_rdf,
            "unit": "dimensionless",
            "evidence_class": "design assumption",
            "source": "none",
            "rationale": (
                "Same BM3R1 operator as the RDF promoter; avoids an unsupported "
                "RPU-to-uM conversion. Scanned from 0.5 to 2.0."
            ),
        },
        {
            "parameter": "n_shutdown",
            "value": shutdown.n_shutdown,
            "unit": "dimensionless",
            "evidence_class": "direct literature value",
            "source": shin,
            "rationale": "Hill coefficient reported for the B3-BM3R1 gate.",
        },
        {
            "parameter": "promoter_leak",
            "value": shutdown.promoter_leak,
            "unit": "fraction of maximum",
            "evidence_class": "derived from literature",
            "source": shin,
            "rationale": "B3-BM3R1 ymin/ymax = 0.005/0.6.",
        },
        {
            "parameter": "copy_ratio_to_counter",
            "value": shutdown.copy_ratio_to_counter,
            "unit": "dimensionless",
            "evidence_class": "design assumption",
            "source": "none",
            "rationale": "One shutdown cassette per counter plasmid at the reference design point.",
        },
        {
            "parameter": "operator_load_stoich",
            "value": shutdown.operator_load_stoich,
            "unit": "BM3R1 monomers/operator",
            "evidence_class": "mechanistic assumption",
            "source": "TetR-family homodimer architecture",
            "rationale": "Two monomer equivalents are sequestered per occupied operator.",
        },
        {
            "parameter": "srna_max_tx_uM_h",
            "value": shutdown.srna_max_tx_uM_h,
            "unit": "uM/h",
            "evidence_class": "within literature range",
            "source": levine,
            "rationale": (
                "1 nM/min * 60 min/h / 1000 nM/uM = 0.06 uM/h. This is the "
                "reference value used in Levine Figure 1C and lies inside the "
                "estimated 0.1-10 nM/min range. "
                "It also exceeds the 0.5 nM/min target reference rate."
            ),
        },
        {
            "parameter": "srna_decay_h",
            "value": shutdown.srna_decay_h,
            "unit": "1/h",
            "evidence_class": "converted literature value",
            "source": f"{levine}; {masse}; {moll}",
            "rationale": (
                "Levine's parameter-estimation section approximates beta_s as "
                "1/50 per min (~0.02/min) from an approximately 30 min RyhB "
                "half-life reported through Masse and Moll, then 0.02/min * "
                "60 min/h = 1.2/h. The Levine effective-rate approximation is "
                "retained rather than recomputing ln(2)/half-life."
            ),
        },
        {
            "parameter": "srna_cocodegradation_fraction",
            "value": shutdown.srna_cocodegradation_fraction,
            "unit": "dimensionless",
            "evidence_class": "declared literature-model scenario",
            "source": levine,
            "rationale": (
                "p=1 is Levine's idealized one-to-one limit, not an experimentally "
                "identified RyhB constant. Partial co-degradation p=0.5-1 is "
                "included in sensitivity and joint-envelope analyses."
            ),
        },
        {
            "parameter": "pair_rate_uM_h",
            "value": shutdown.pair_rate_uM_h,
            "unit": "1/(uM h)",
            "evidence_class": "converted literature value",
            "source": f"{levine}; {masse}",
            "rationale": (
                "Levine inferred 1/50 = 0.02/(nM min) from Masse's observation "
                "that RyhB disappeared within about 3 min in the presence of an "
                "estimated 20 nM target. Conversion: 0.02 * 1000 nM/uM * "
                "60 min/h = 1200/(uM h); this is not a direct binding-rate fit."
            ),
        },
        {
            "parameter": "output_max_tx_uM_h",
            "value": shutdown.output_max_tx_uM_h,
            "unit": "uM/h",
            "evidence_class": "literature-scale design point",
            "source": f"{levine}; {zhang}; {selinger}",
            "rationale": (
                "Project baseline 0.5 nM/min * 60 min/h / 1000 nM/uM = "
                "0.03 uM/h. Levine estimated alpha_m near 1 nM/min from "
                "10-20 copies/cell transcriptome measurements by Zhang and "
                "Selinger, and noted roughly ten-fold physiological variation, "
                "but did not specify the one-sided 0.1-1 nM/min interval used "
                "here as a conservative project envelope."
            ),
        },
        {
            "parameter": "output_mrna_decay_h",
            "value": shutdown.output_mrna_decay_h,
            "unit": "1/h",
            "evidence_class": "converted literature value",
            "source": f"{levine}; {masse}",
            "rationale": (
                "Levine approximated beta_m as 1/10 per min (~0.1/min) from "
                "Masse's approximately 6 min sodB mRNA half-life; 0.1/min * "
                "60 min/h = 6/h. Used as total effective RNA loss."
            ),
        },
        {
            "parameter": "output_translation_h",
            "value": shutdown.output_translation_h,
            "unit": "1/h",
            "evidence_class": "normalization convention",
            "source": "none",
            "rationale": "Unit scaling only; normalized pulse timing is invariant to this factor.",
        },
        {
            "parameter": "output_protein_decay_h",
            "value": shutdown.output_protein_decay_h,
            "unit": "1/h",
            "evidence_class": "reporter assumption",
            "source": "stable GFP readout convention",
            "rationale": "No active degradation; protein is removed by growth dilution only.",
        },
        {
            "parameter": "counter_K_rep_uM",
            "value": counter.K_rep_uM,
            "unit": "uM",
            "evidence_class": "model-selected engineering point",
            "source": "rdfmodel_new scan; Zhao network",
            "rationale": "Inherited counter calibration; not a direct Zhao or Cello constant.",
        },
        {
            "parameter": "counter_n_rep",
            "value": counter.n_rep,
            "unit": "dimensionless",
            "evidence_class": "model-selected engineering point",
            "source": "rdfmodel_new scan; BM3R1 gate-family context",
            "rationale": "Inherited counter working point; not refitted in this module.",
        },
        {
            "parameter": "counter_krep_tsl_h",
            "value": counter.krep_tsl_h,
            "unit": "1/h",
            "evidence_class": "model-selected engineering point",
            "source": "rdfmodel_new scan; Zhao network",
            "rationale": "Selected for one input/one reversal; Zhao default is 0.3/h.",
        },
        {
            "parameter": "counter_krdf_tsl_h",
            "value": counter.krdf_tsl_h,
            "unit": "1/h",
            "evidence_class": "model-selected engineering point",
            "source": "rdfmodel_new scan; Zhao network",
            "rationale": "Selected for one input/one reversal; Zhao default is 4/h.",
        },
        {
            "parameter": "counter_k_tag_int_h",
            "value": counter.k_tag_int_h,
            "unit": "1/h",
            "evidence_class": "model-selected engineering point",
            "source": "rdfmodel_new tag scan",
            "rationale": "Selected inside the counter's one-peak/one-reversal region.",
        },
        {
            "parameter": "counter_c31_rbs_scale",
            "value": counter.c31_rbs_scale,
            "unit": "dimensionless",
            "evidence_class": "model-selected interface point",
            "source": "rdfmodel_new interface scan",
            "rationale": "Scales oscillator C31 translation flux at the A/B interface.",
        },
        {
            "parameter": "counter_k_dil_h",
            "value": counter.k_dil_h,
            "unit": "1/h",
            "evidence_class": "derived from project growth condition",
            "source": "oscillator project nominal culture condition",
            "rationale": (
                "ln(2)/(50/60 h) for the nominal 50 min doubling time; this is "
                "not a Zhao kinetic constant."
            ),
        },
    ]


class C31Waveform:
    """Physical A-to-B interface based on C31 translation-production flux."""

    def __init__(
        self,
        csv_path: Path,
        k_dil_h: float,
        rbs_scale: float = 0.45,
    ) -> None:
        df = pd.read_csv(csv_path)
        required = {"time_h", "c31_uM", "translation_flux_uM_h"}
        missing = required.difference(df.columns)
        if missing:
            raise ValueError(f"Extracted oscillator waveform is missing: {sorted(missing)}")
        self.t_h = df["time_h"].to_numpy(float)
        self.c31_uM = np.maximum(0.0, df["c31_uM"].to_numpy(float))
        self.raw_flux_uM_h = np.maximum(
            0.0, df["translation_flux_uM_h"].to_numpy(float)
        )
        self.scaled_flux_uM_h = float(rbs_scale) * self.raw_flux_uM_h
        self.k_dil_h = float(k_dil_h)
        self.rbs_scale = float(rbs_scale)
        self._concentration_spl = PchipInterpolator(
            self.t_h, self.c31_uM, extrapolate=False
        )
        self._raw_flux_spl = PchipInterpolator(
            self.t_h, self.raw_flux_uM_h, extrapolate=False
        )
        self._source_spl = PchipInterpolator(
            self.t_h, self.scaled_flux_uM_h, extrapolate=False
        )
        self._peak_times_h, self._trough_times_h = self._detect_peaks_and_troughs()

    @property
    def end_h(self) -> float:
        return float(self.t_h[-1])

    def concentration(self, t: float | np.ndarray) -> float | np.ndarray:
        tc = np.clip(t, self.t_h[0], self.t_h[-1])
        return np.maximum(0.0, self._concentration_spl(tc))

    def translation_flux(self, t: float | np.ndarray) -> float | np.ndarray:
        tc = np.clip(t, self.t_h[0], self.t_h[-1])
        return np.maximum(0.0, self._raw_flux_spl(tc))

    def source(self, t: float) -> float:
        tc = float(np.clip(t, self.t_h[0], self.t_h[-1]))
        return max(0.0, float(self._source_spl(tc)))

    def _detect_peaks_and_troughs(self) -> tuple[np.ndarray, np.ndarray]:
        span = float(np.ptp(self.scaled_flux_uM_h))
        prominence = max(1e-9, 0.05 * span)
        dt = float(np.median(np.diff(self.t_h)))
        min_distance = max(1, int(round(2.0 / dt)))
        peak_idx, _ = find_peaks(
            self.scaled_flux_uM_h,
            prominence=prominence,
            distance=min_distance,
        )
        trough_idx = []
        for left, right in zip(peak_idx[:-1], peak_idx[1:]):
            local = left + int(np.argmin(self.scaled_flux_uM_h[left : right + 1]))
            trough_idx.append(local)
        return self.t_h[peak_idx], self.t_h[np.asarray(trough_idx, dtype=int)]

    def peaks_and_troughs(self) -> tuple[np.ndarray, np.ndarray]:
        return self._peak_times_h.copy(), self._trough_times_h.copy()

    def effective_source(self, t: np.ndarray) -> np.ndarray:
        return np.asarray([self.source(float(value)) for value in t])


def counter_params(config: CounterConfig) -> dict:
    params = zc.default_params()
    params.update(
        K_rep=config.K_rep_uM,
        n_rep=config.n_rep,
        krep_tsl=config.krep_tsl_h,
        krdf_tsl=config.krdf_tsl_h,
        k_tag_int=config.k_tag_int_h,
        k_dil=config.k_dil_h,
    )
    return params


def repression_off_fraction(rep_free: float, K: float, n: float, leak: float) -> float:
    rep_free = max(0.0, float(rep_free))
    return leak + (1.0 - leak) / (1.0 + (rep_free / K) ** n)


def operator_occupancy(rep_free: float, K: float, n: float) -> float:
    ratio = (max(0.0, rep_free) / K) ** n
    return ratio / (1.0 + ratio)


def free_repressor(
    rep_total: float,
    dna_load_uM: float,
    cfg: ShutdownConfig,
    shutdown_K_uM: float,
) -> float:
    """Rapid-equilibrium mass balance for the added BM3R1 operator load."""
    rep_total = max(0.0, float(rep_total))
    if rep_total == 0.0 or dna_load_uM <= 0.0:
        return rep_total

    def balance(rep_free: float) -> float:
        bound = (
            cfg.operator_load_stoich
            * dna_load_uM
            * operator_occupancy(
                rep_free, shutdown_K_uM, cfg.n_shutdown
            )
        )
        return rep_free + bound - rep_total

    return float(brentq(balance, 0.0, rep_total, xtol=1e-12, rtol=1e-10))


def initial_state(
    params: dict,
    shutdown: ShutdownConfig,
    counter_initial_state: str = "PB",
) -> np.ndarray:
    """Build a full state from either counter memory state."""
    rep_mrna_ss = params["k_tscr"] * params["Dtot"] / params["k_rna"]
    rep_ss = params["krep_tsl"] * rep_mrna_ss / params["k_dil"]
    y = np.zeros(41)
    if counter_initial_state == "PB":
        y[:38] = zc.y0_PB(params, rep_mrna=rep_mrna_ss, rep=rep_ss)
    elif counter_initial_state == "LR":
        y[0] = params["Dtot"]
        y[22] = 1.0
        y[32] = 0.5
        y[36] = 1e-4
        y[37] = 1e-4
    else:
        raise ValueError("counter_initial_state must be 'PB' or 'LR'")

    rep_total = float(y[37])
    shutdown_K = params["K_rep"] * shutdown.threshold_ratio_to_rdf
    rep_free = free_repressor(
        rep_total,
        shutdown.copy_ratio_to_counter * params["Dtot"],
        shutdown,
        shutdown_K,
    )
    activity = repression_off_fraction(
        rep_free,
        shutdown_K,
        shutdown.n_shutdown,
        shutdown.promoter_leak,
    )
    srna_prod = (
        shutdown.srna_max_tx_uM_h
        * shutdown.copy_ratio_to_counter
        * activity
    )
    y[38] = srna_prod / shutdown.srna_decay_h
    return y


def shutdown_rhs(
    lr_fraction: float,
    rep_free_uM: float,
    state: np.ndarray,
    params: dict,
    shutdown: ShutdownConfig,
) -> tuple[np.ndarray, float]:
    """Evaluate the three-state Levine-style shutdown subsystem."""
    srna, output_mrna, output_protein = np.maximum(state, 0.0)
    shutdown_K = params["K_rep"] * shutdown.threshold_ratio_to_rdf
    activity = repression_off_fraction(
        rep_free_uM,
        shutdown_K,
        shutdown.n_shutdown,
        shutdown.promoter_leak,
    )
    srna_production = (
        shutdown.srna_max_tx_uM_h
        * shutdown.copy_ratio_to_counter
        * activity
    )
    pair_loss = shutdown.pair_rate_uM_h * srna * output_mrna
    rates = np.array(
        [
            srna_production
            - shutdown.srna_decay_h * srna
            - shutdown.srna_cocodegradation_fraction * pair_loss,
            shutdown.output_max_tx_uM_h * max(0.0, float(lr_fraction))
            - shutdown.output_mrna_decay_h * output_mrna
            - pair_loss,
            shutdown.output_translation_h * output_mrna
            - (shutdown.output_protein_decay_h + params["k_dil"])
            * output_protein,
        ]
    )
    return rates, activity


def coupled_rhs(
    t: float,
    y: np.ndarray,
    params: dict,
    constants: dict,
    int_source: Callable[[float], float],
    shutdown: ShutdownConfig,
) -> np.ndarray:
    core_y = y[:38]
    dy = np.zeros(41)
    dy[:38] = zc.rhs(t, core_y, params, constants, int_source)

    lr_total = float(zc.LR_total(core_y[np.newaxis, :])[0])
    lr_fraction = lr_total / params["Dtot"]
    rep_total = max(0.0, float(core_y[37]))
    dna_load = shutdown.copy_ratio_to_counter * params["Dtot"]
    shutdown_K = params["K_rep"] * shutdown.threshold_ratio_to_rdf
    rep_free = free_repressor(rep_total, dna_load, shutdown, shutdown_K)

    # Replace the original RDF transcription term with free-BM3R1 regulation.
    rdf_activity = repression_off_fraction(
        rep_free,
        params["K_rep"],
        params["n_rep"],
        params["leak_rep"],
    )
    dy[32] = params["k_tscr"] * lr_total * rdf_activity - params["k_rna"] * core_y[32]

    dy[38:41], _ = shutdown_rhs(
        lr_fraction,
        rep_free,
        y[38:41],
        params,
        shutdown,
    )
    return dy


def simulate_full(
    waveform_csv: Path,
    counter: CounterConfig | None = None,
    shutdown: ShutdownConfig | None = None,
    t_end_h: float = 100.0,
    dt_h: float = 0.02,
    counter_initial_state: str = "PB",
):
    counter = counter or CounterConfig()
    shutdown = shutdown or ShutdownConfig()
    params = counter_params(counter)
    waveform = C31Waveform(
        waveform_csv,
        params["k_dil"],
        rbs_scale=counter.c31_rbs_scale,
    )
    t_end_h = min(float(t_end_h), waveform.end_h)
    t_eval = np.arange(0.0, t_end_h + 0.5 * dt_h, dt_h)
    constants = zc._rate_constants(params)
    y0 = initial_state(params, shutdown, counter_initial_state)
    sol = solve_ivp(
        lambda t, y: coupled_rhs(
            t, y, params, constants, waveform.source, shutdown
        ),
        (0.0, t_end_h),
        y0,
        method="LSODA",
        rtol=1e-6,
        atol=1e-11,
        max_step=0.05,
        t_eval=t_eval,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol, params, waveform


def observables(sol, params: dict, shutdown: ShutdownConfig) -> pd.DataFrame:
    y = sol.y.T
    lr_frac = zc.LR_total(y[:, :38]) / params["Dtot"]
    pb_frac = zc.PB_total(y[:, :38]) / params["Dtot"]
    rep_total = y[:, 37]
    dna_load = shutdown.copy_ratio_to_counter * params["Dtot"]
    shutdown_K = params["K_rep"] * shutdown.threshold_ratio_to_rdf
    rep_free = np.array(
        [
            free_repressor(value, dna_load, shutdown, shutdown_K)
            for value in rep_total
        ]
    )
    shutdown_activity = np.array(
        [
            repression_off_fraction(
                value,
                shutdown_K,
                shutdown.n_shutdown,
                shutdown.promoter_leak,
            )
            for value in rep_free
        ]
    )
    return pd.DataFrame(
        {
            "time_h": sol.t,
            "LR_fraction": lr_frac,
            "PB_fraction": pb_frac,
            "integrase_total_uM": zc.int_total(y[:, :38]),
            "RDF_total_uM": zc.rdf_total(y[:, :38]),
            "BM3R1_total_uM": rep_total,
            "BM3R1_free_uM": rep_free,
            "shutdown_promoter_activity": shutdown_activity,
            "shutdown_sRNA_uM": y[:, 38],
            "output_mRNA_uM": y[:, 39],
            "output_translation_uM_h": shutdown.output_translation_h * y[:, 39],
            "output_protein_uM": y[:, 40],
        }
    )


def simulate_downstream(
    driver: pd.DataFrame,
    params: dict,
    shutdown: ShutdownConfig,
) -> pd.DataFrame:
    """Fast readout-only simulation for tunability scans.

    This intentionally does not feed operator load back into the counter. Full
    coupled load effects are assessed separately with simulate_full().
    """
    t = driver["time_h"].to_numpy(float)
    lr_spl = PchipInterpolator(t, driver["LR_fraction"].to_numpy(float))
    rep_spl = PchipInterpolator(t, driver["BM3R1_total_uM"].to_numpy(float))
    dna_load = shutdown.copy_ratio_to_counter * params["Dtot"]

    def rhs(time: float, state: np.ndarray) -> np.ndarray:
        lr_total = params["Dtot"] * float(np.clip(lr_spl(time), 0.0, 1.0))
        lr_fraction = lr_total / params["Dtot"]
        shutdown_K = params["K_rep"] * shutdown.threshold_ratio_to_rdf
        rep_free = free_repressor(
            float(rep_spl(time)), dna_load, shutdown, shutdown_K
        )
        rates, _ = shutdown_rhs(
            lr_fraction,
            rep_free,
            state,
            params,
            shutdown,
        )
        return rates

    sol = solve_ivp(
        rhs,
        (float(t[0]), float(t[-1])),
        np.zeros(3),
        method="LSODA",
        rtol=2e-6,
        atol=1e-10,
        max_step=0.05,
        t_eval=t,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    out = driver[["time_h", "LR_fraction", "BM3R1_total_uM"]].copy()
    out["shutdown_sRNA_uM"] = sol.y[0]
    out["output_mRNA_uM"] = sol.y[1]
    out["output_translation_uM_h"] = shutdown.output_translation_h * sol.y[1]
    out["output_protein_uM"] = sol.y[2]
    return out


def toggle_score(
    t: np.ndarray,
    lr_fraction: np.ndarray,
    troughs_h: np.ndarray,
    skip: int = 2,
    transient_h: float = 25.0,
) -> dict:
    """Banded P score using the rdfmodel_new steady-state sampling rule."""
    lower = max(float(t[0]), float(transient_h))
    troughs = troughs_h[(troughs_h >= lower) & (troughs_h <= t[-1])]
    values = np.interp(troughs, t, lr_fraction)
    if len(values) <= skip + 2:
        return {"score": np.nan, "high": np.nan, "low": np.nan, "fidelity": np.nan}
    values = values[skip:]
    high_mask = values > 0.6
    low_mask = values < 0.4
    valid = high_mask | low_mask
    transitions = valid[:-1] & valid[1:]
    fidelity = (
        np.mean(high_mask[:-1][transitions] != high_mask[1:][transitions])
        if transitions.any()
        else 0.0
    )
    high = float(values[high_mask].mean()) if high_mask.any() else np.nan
    low = float(values[low_mask].mean()) if low_mask.any() else np.nan
    core = min(np.nan_to_num(high, nan=0.0), 1.0 - np.nan_to_num(low, nan=1.0))
    return {
        "score": float(core * fidelity),
        "high": high,
        "low": low,
        "fidelity": float(fidelity),
    }


def cycle_switch_audit(
    df: pd.DataFrame,
    waveform: C31Waveform,
    transient_h: float = 25.0,
    midpoint: float = 0.5,
    complete_threshold: float = 0.95,
) -> tuple[pd.DataFrame, dict]:
    """Count zero, single, and multiple DNA-state crossings per input cycle."""
    t = df["time_h"].to_numpy(float)
    lr = df["LR_fraction"].to_numpy(float)
    pb = df["PB_fraction"].to_numpy(float)
    peaks, _ = waveform.peaks_and_troughs()
    peaks = peaks[(peaks >= transient_h) & (peaks <= t[-1])]
    boundaries = 0.5 * (peaks[:-1] + peaks[1:])
    rows = []
    for cycle_index, (left, right) in enumerate(
        zip(boundaries[:-1], boundaries[1:]), start=1
    ):
        local = np.flatnonzero((t >= left) & (t < right))
        if len(local) < 3:
            continue
        events = []
        for i0, i1 in zip(local[:-1], local[1:]):
            a = lr[i0] - midpoint
            b = lr[i1] - midpoint
            if a * b < 0.0:
                fraction = abs(a) / (abs(a) + abs(b))
                event_time = t[i0] + fraction * (t[i1] - t[i0])
                direction = "PB_to_LR" if b > a else "LR_to_PB"
                events.append((event_time, direction))
        end = local[-1]
        complete_fraction = max(float(lr[end]), float(pb[end]))
        rows.append(
            {
                "cycle_index": cycle_index,
                "cycle_start_h": float(left),
                "cycle_end_h": float(right),
                "total_crossings": len(events),
                "single_crossing_direction": events[0][1] if len(events) == 1 else "",
                "endpoint_complete_fraction": complete_fraction,
                "endpoint_complete": bool(complete_fraction >= complete_threshold),
                "one_and_complete": bool(
                    len(events) == 1 and complete_fraction >= complete_threshold
                ),
                "input_peak_uM_h": float(
                    np.max(waveform.effective_source(t[local]))
                ),
            }
        )
    audit = pd.DataFrame(rows)
    if audit.empty:
        summary = {
            "audited_cycles": 0,
            "fraction_zero_crossings": np.nan,
            "fraction_exactly_one_crossing": np.nan,
            "fraction_multiple_crossings": np.nan,
            "fraction_one_and_complete": np.nan,
            "alternation_fidelity": np.nan,
            "system_pass": False,
        }
        return audit, summary

    counts = audit["total_crossings"].to_numpy(int)
    single = audit["single_crossing_direction"].to_numpy(str)
    eligible = (counts[:-1] == 1) & (counts[1:] == 1)
    alternation = (
        float(np.mean(single[:-1][eligible] != single[1:][eligible]))
        if eligible.any()
        else np.nan
    )
    one_complete = float(audit["one_and_complete"].mean())
    summary = {
        "audited_cycles": int(len(audit)),
        "fraction_zero_crossings": float(np.mean(counts == 0)),
        "fraction_exactly_one_crossing": float(np.mean(counts == 1)),
        "fraction_multiple_crossings": float(np.mean(counts > 1)),
        "fraction_complete_endpoint": float(audit["endpoint_complete"].mean()),
        "fraction_one_and_complete": one_complete,
        "alternation_fidelity": alternation,
        "system_pass": bool(
            one_complete >= 0.95
            and np.isfinite(alternation)
            and alternation >= 0.95
        ),
    }
    return audit, summary


def _upcrossings(t: np.ndarray, signal: np.ndarray, threshold: float = 0.5) -> np.ndarray:
    idx = np.flatnonzero((signal[:-1] < threshold) & (signal[1:] >= threshold)) + 1
    return t[idx]


def _downcrossings(t: np.ndarray, signal: np.ndarray, threshold: float = 0.5) -> np.ndarray:
    idx = np.flatnonzero((signal[:-1] >= threshold) & (signal[1:] < threshold)) + 1
    return t[idx]


def pulse_metrics(df: pd.DataFrame, skip_first: int = 1) -> pd.DataFrame:
    """Measure production and protein pulse widths between PB->LR events."""
    t = df["time_h"].to_numpy(float)
    lr = df["LR_fraction"].to_numpy(float)
    production = df["output_translation_uM_h"].to_numpy(float)
    protein = df["output_protein_uM"].to_numpy(float)
    starts = _upcrossings(t, lr)
    stops = _downcrossings(t, lr)
    rows = []
    for pulse_index, start in enumerate(starts):
        if pulse_index < skip_first:
            continue
        end = starts[pulse_index + 1] if pulse_index + 1 < len(starts) else t[-1]
        lr_stops = stops[(stops > start) & (stops < end)]
        lr_end = float(lr_stops[0]) if len(lr_stops) else float(end)
        mask = (t >= start) & (t < end)
        if mask.sum() < 5:
            continue
        row = {
            "pulse_index": pulse_index,
            "start_h": float(start),
            "LR_end_h": lr_end,
            "window_end_h": float(end),
        }
        for name, signal in (("production", production), ("protein", protein)):
            local_t = t[mask]
            local = signal[mask]
            peak = float(np.max(local))
            area = float(np.trapz(local, local_t))
            row[f"{name}_peak"] = peak
            row[f"{name}_auc"] = area
            if peak <= 0.0:
                row[f"{name}_fwhm_h"] = np.nan
                row[f"{name}_closure_ratio"] = np.nan
                continue
            above = np.flatnonzero(local >= 0.5 * peak)
            row[f"{name}_fwhm_h"] = (
                float(local_t[above[-1]] - local_t[above[0]]) if len(above) else np.nan
            )
            # Evaluate shutdown while LR is still active. Sampling after the
            # LR->PB transition would falsely credit state inversion for sRNA
            # closure.
            late_lr_start = start + 0.8 * (lr_end - start)
            late_lr = (local_t >= late_lr_start) & (local_t < lr_end)
            row[f"{name}_closure_ratio"] = (
                float(np.mean(local[late_lr]) / peak) if late_lr.any() else np.nan
            )
        rows.append(row)
    return pd.DataFrame(rows)


def config_dict(counter: CounterConfig, shutdown: ShutdownConfig) -> dict:
    return {
        "counter_model": {
            "source": "rdfmodel_new/model/zhao_core.py",
            "states": 38,
            "nominal_growth_basis": "50 min doubling time",
            "dual_initial_state_audit": True,
        },
        "counter": asdict(counter),
        "shutdown": asdict(shutdown),
        "parameter_provenance": parameter_provenance(counter, shutdown),
    }


__all__ = [
    "C31Waveform",
    "CounterConfig",
    "ShutdownConfig",
    "config_dict",
    "counter_params",
    "cycle_switch_audit",
    "observables",
    "parameter_provenance",
    "pulse_metrics",
    "replace",
    "simulate_downstream",
    "simulate_full",
    "shutdown_rhs",
    "toggle_score",
]
