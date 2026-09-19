export const PERIOD_KNOBS = Object.freeze({
  1: Object.freeze({ id: "K1", rbs: 0.60, mrnaHalfLifeMin: 1.0, medianGenerations: 7.25458, nominalHours: 6.0455, passingTags: [] }),
  2: Object.freeze({ id: "K2", rbs: 0.80, mrnaHalfLifeMin: 1.0, medianGenerations: 8.31446, nominalHours: 6.9287, passingTags: [] }),
  3: Object.freeze({ id: "K3", rbs: 1.20, mrnaHalfLifeMin: 1.0, medianGenerations: 9.76918, nominalHours: 8.1408, passingTags: [] }),
  4: Object.freeze({ id: "K4", rbs: 1.20, mrnaHalfLifeMin: 1.5, medianGenerations: 11.2208, nominalHours: 9.3507, passingTags: [8, 12] }),
  5: Object.freeze({ id: "K5", rbs: 0.60, mrnaHalfLifeMin: 4.0, medianGenerations: 12.5416, nominalHours: 10.4513, passingTags: [12] }),
});

export function getPeriodKnob(value) {
  const index = Math.min(5, Math.max(1, Math.round(Number(value) || 4)));
  return PERIOD_KNOBS[index];
}

export const DEFAULTS = Object.freeze({
  hours: 60,
  period: 4,
  pulseWidth: 3.6,
  promoter: 1,
  oscRbs: 1,
  c31: 0.45,
  degradation: 8,
  rdf: 1,
  bm3r1: 1,
  srna: 0.06,
});

export const PRESETS = Object.freeze({
  tempo: { ...DEFAULTS },
  stable: {
    ...DEFAULTS,
    period: 5,
    pulseWidth: 3.3,
    c31: 0.42,
    degradation: 12,
    rdf: 1.05,
    srna: 0.065,
  },
  fast: {
    ...DEFAULTS,
    period: 3,
    pulseWidth: 2.5,
    c31: 0.52,
    degradation: 14,
    rdf: 1.15,
    bm3r1: 0.85,
    srna: 0.072,
  },
  long: {
    ...DEFAULTS,
    bm3r1: 1.35,
    srna: 0.04,
  },
  short: {
    ...DEFAULTS,
    bm3r1: 0.72,
    srna: 0.085,
  },
});

export const CONTROL_IDS = Object.freeze([
  "period",
  "pulse-width",
  "promoter",
  "osc-rbs",
  "c31",
  "degradation",
  "rdf",
  "bm3r1",
  "srna",
]);

export const PARAMETER_MAP = Object.freeze([
  { design: "Certified period setting K1–K5", model: "oscillator_translation_scale + oscillator_mrna_total_half_life_min", reference: "Robust five-setting library at Td = 50 min" },
  { design: "C31 pulse-width target", model: "c31_mrna_total_half_life_min", reference: "2.0 min at 3.6 h target" },
  { design: "Oscillator promoter strength", model: "oscillator_tx_per_plasmid_per_min + c31_tx_per_plasmid_per_min", reference: "shared PLtetO1 multiplier" },
  { design: "Oscillator RBS multiplier", model: "oscillator_translation_per_mrna_per_min", reference: "multiplier around selected K setting" },
  { design: "C31 RBS strength", model: "c31_translation_scale", reference: "0.45× candidate, applied once in A waveform" },
  { design: "Integrase degradation tag", model: "k_tag_int_h", reference: "8 h⁻¹ candidate" },
  { design: "RDF expression", model: "krdf_tsl_h", reference: "200 h⁻¹ at 1.00×" },
  { design: "BM3R1 expression", model: "krep_tsl_h", reference: "15 h⁻¹ at 1.00×" },
  { design: "sRNA promoter strength", model: "srna_max_tx_uM_h", reference: "0.06 μM/h" },
]);
