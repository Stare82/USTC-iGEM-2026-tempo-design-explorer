import { getPeriodKnob } from "./parameters.js";

const clamp = (value, min, max) => Math.min(max, Math.max(min, value));
const sigmoid = (value) => 1 / (1 + Math.exp(-value));

function proximityScore(value, low, high, softness) {
  if (value >= low && value <= high) return 1;
  const distance = value < low ? low - value : value - high;
  return Math.exp(-Math.pow(distance / softness, 2));
}

export function evaluateDesign(params) {
  const expressionGain = Math.sqrt(params.promoter * params.oscRbs);
  const normalizedDose = (params.c31 / 0.45) * (params.pulseWidth / 3.6) * expressionGain;
  const doseScore = proximityScore(normalizedDose, 0.68, 1.48, 0.36);
  const tagScore = proximityScore(params.degradation, 6, 18, 5.2);
  const rdfScore = proximityScore(params.rdf, 0.62, 1.38, 0.34);
  const periodScore = proximityScore(params.period, 7.2, 14.5, 2.4);

  const recoveryTime = params.pulseWidth + 1.15 / params.rdf + 7.2 / params.degradation;
  const recoveryMargin = params.period - recoveryTime;
  const recoveryScore = clamp((recoveryMargin + 0.3) / 2.8, 0, 1);
  const flipSuccess = clamp(
    100 * (0.42 * doseScore + 0.26 * tagScore + 0.18 * rdfScore + 0.14 * recoveryScore),
    0,
    100,
  );

  const stable = flipSuccess >= 90 && recoveryMargin >= 0.7 && params.rdf <= 1.5;
  const risky = flipSuccess >= 70 && recoveryMargin >= -0.2;
  const leakagePercent = clamp(
    3.5 + Math.max(0, params.promoter - 1) * 9 + Math.max(0, 0.38 - params.c31) * 22 + Math.max(0, params.rdf - 1.2) * 14,
    2,
    34,
  );

  const firstPulse = params.period * 0.58;
  const outputStart = firstPulse + 0.62 + 0.3 / Math.max(0.4, params.c31 / 0.45);
  const outputDuration = clamp(
    4.25 * Math.pow(params.bm3r1, 0.58) * Math.pow(0.06 / params.srna, 0.72),
    1.2,
    Math.max(1.3, params.period - 1.1),
  );
  const shutdownTime = outputStart + outputDuration;
  const residual = clamp(
    5 + (0.06 / params.srna - 1) * 11 + (params.bm3r1 - 1) * 6,
    2,
    38,
  );
  const closureScore = clamp(1 - residual / 45, 0, 1);

  let score = Math.round(
    100 * (0.47 * flipSuccess / 100 + 0.2 * recoveryScore + 0.12 * (1 - leakagePercent / 40) + 0.21 * closureScore),
  );
  score = clamp(score, 0, 100);

  const insights = [];
  if (recoveryMargin < 0.7) {
    insights.push({ level: recoveryMargin < 0 ? "danger" : "warn", text: "The RDF / Integrase reset is predicted to overlap the next oscillator pulse." });
  }
  if (normalizedDose < 0.68) {
    insights.push({ level: "warn", text: "C31 dose is below the calibrated flip window; increase RBS strength or pulse width." });
  } else if (normalizedDose > 1.48) {
    insights.push({ level: "warn", text: "C31 exposure is above the single-flip window and may permit premature reversal." });
  }
  if (params.degradation > 24) {
    insights.push({ level: "danger", text: "Integrase clearance is too fast for reliable PB → LR completion at this input level." });
  } else if (params.degradation < 5) {
    insights.push({ level: "warn", text: "Slow Integrase clearance increases carry-over between pulses." });
  }
  if (params.rdf > 1.45) {
    insights.push({ level: "warn", text: "RDF remains high relative to the next Integrase pulse; reversal timing becomes fragile." });
  }
  if (residual > 18) {
    insights.push({ level: "warn", text: "The sRNA branch leaves substantial residual output; raise sRNA production or reduce BM3R1 expression." });
  }
  if (!insights.length) {
    insights.push({ level: "ok", text: "RBS 0.45× and the 8 h⁻¹ tag sit inside the current robust counter window." });
  }

  const verdict = stable && residual <= 18 ? "success" : risky ? "warning" : "danger";
  return {
    score,
    verdict,
    title: verdict === "success" ? "Within design window" : verdict === "warning" ? "Narrow design margin" : "Outside design window",
    message: verdict === "success"
      ? "This design is predicted to separate consecutive flips and close the output window."
      : verdict === "warning"
        ? "The chain may function, but at least one interface has little timing margin."
        : "At least one stage is predicted to fail before a complete timed output is produced.",
    flipSuccess,
    stable,
    recovery: recoveryMargin >= 0.7,
    recoveryMargin,
    leakagePercent,
    leakage: leakagePercent < 10 ? "Low" : leakagePercent < 20 ? "Medium" : "High",
    outputStart,
    outputDuration,
    shutdownTime,
    residual,
    normalizedDose,
    insights,
  };
}

function pulseCenters(params) {
  const centers = [];
  for (let center = params.period * 0.58; center <= params.hours + params.period; center += params.period) {
    centers.push(center);
  }
  return centers;
}

function gaussianPulse(t, center, fwhm) {
  const sigma = fwhm / 2.355;
  return Math.exp(-0.5 * Math.pow((t - center) / sigma, 2));
}

function dnaState(t, centers, quality) {
  let state = 0.015;
  centers.forEach((center, index) => {
    const progress = sigmoid((t - center - 0.28) / 0.2);
    const target = index % 2 === 0 ? 0.985 : 0.015;
    const completeness = clamp(0.38 + 0.62 * quality, 0.2, 1);
    state += (target - state) * progress * completeness;
  });
  return clamp(state, 0, 1);
}

export function simulateTempo(params) {
  const knob = getPeriodKnob(params.period);
  const modelParams = { ...params, period: knob.nominalHours };
  const evaluation = evaluateDesign(modelParams);
  const dt = 0.08;
  const centers = pulseCenters(modelParams);
  const samples = [];
  let integrase = 0;
  let rdf = 0;
  let bm3r1 = 0;
  let srna = 0;
  let outputMrna = 0;
  let outputProtein = 0;
  let previousLR = 0.015;
  let lastLRRise = -999;

  const expressionGain = Math.sqrt(modelParams.promoter * modelParams.oscRbs);
  const inputAmplitude = 6.639 * (modelParams.c31 / 0.45) * expressionGain;

  for (let t = 0; t <= modelParams.hours + dt / 2; t += dt) {
    const c31 = centers.reduce((sum, center) => sum + gaussianPulse(t, center, modelParams.pulseWidth), 0) * inputAmplitude;
    const delayed = centers.reduce((sum, center) => sum + gaussianPulse(t, center + modelParams.pulseWidth * 0.58, modelParams.pulseWidth * 0.82), 0);

    const integraseSource = 0.72 * c31;
    integrase += dt * (integraseSource - (0.7 + 0.07 * modelParams.degradation) * integrase);
    rdf += dt * (2.1 * modelParams.rdf * delayed - (0.72 + 0.055 * modelParams.degradation) * rdf);

    const lr = dnaState(t, centers, evaluation.flipSuccess / 100);
    const pb = 1 - lr;
    if (previousLR < 0.5 && lr >= 0.5) lastLRRise = t;
    previousLR = lr;

    const lrAge = Math.max(0, t - lastLRRise);
    const bmTarget = lr > 0.5 ? modelParams.bm3r1 * Math.exp(-lrAge / (2.55 * modelParams.bm3r1)) : 0.08 * modelParams.bm3r1;
    bm3r1 += dt * (1.25 * (bmTarget - bm3r1));
    const srnaGate = 1 / (1 + Math.pow(Math.max(0, bm3r1) / 0.24, 3.1));
    srna += dt * (modelParams.srna * srnaGate - 1.2 * srna);

    const pairLoss = 22 * srna * outputMrna;
    outputMrna += dt * (0.9 * lr - 1.55 * outputMrna - pairLoss);
    outputProtein += dt * (1.2 * outputMrna - 0.42 * outputProtein);

    samples.push({
      t,
      c31: Math.max(0, c31),
      integrase: Math.max(0, integrase),
      rdf: Math.max(0, rdf),
      lr,
      pb,
      bm3r1: Math.max(0, bm3r1),
      srna: Math.max(0, srna),
      output: Math.max(0, outputProtein),
    });
  }

  const visibleCenters = centers.filter((center) => center <= modelParams.hours);
  const successfulFlips = Math.max(0, Math.round(visibleCenters.length * evaluation.flipSuccess / 100));
  return {
    samples,
    centers: visibleCenters,
    evaluation,
    summary: {
      pulses: visibleCenters.length,
      flips: successfulFlips,
      outputs: Math.ceil(successfulFlips / 2),
    },
    metadata: {
      source: "browser_surrogate",
      period_setting_id: knob.id,
      target_period_h: knob.nominalHours,
      realized_period_h: knob.nominalHours,
    },
  };
}
