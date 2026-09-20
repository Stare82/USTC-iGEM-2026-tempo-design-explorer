import { CONTROL_IDS, DEFAULTS, PARAMETER_MAP, PRESETS, getPeriodKnob } from "./model/parameters.js";
import { simulateTempo } from "./model/tempo-model.js";
import { renderChart, renderStateRail } from "./charts.js";

const CONTROL_KEYS = {
  period: "period",
  "pulse-width": "pulseWidth",
  promoter: "promoter",
  "osc-rbs": "oscRbs",
  c31: "c31",
  degradation: "degradation",
  rdf: "rdf",
  bm3r1: "bm3r1",
  srna: "srna",
};

const FORMATTERS = {
  period: (value) => {
    const knob = getPeriodKnob(value);
    return `${knob.id} · ${knob.nominalHours.toFixed(2)} h`;
  },
  "pulse-width": (value) => `${value.toFixed(1)} h`,
  promoter: (value) => `${value.toFixed(2)}×`,
  "osc-rbs": (value) => `${value.toFixed(2)}×`,
  c31: (value) => `${value.toFixed(2)}×`,
  degradation: (value) => `${value.toFixed(1)} h⁻¹`,
  rdf: (value) => `${value.toFixed(2)}×`,
  bm3r1: (value) => `${value.toFixed(2)}×`,
  srna: (value) => `${value.toFixed(3)} μM/h`,
};

const OUTPUT_IDS = {
  period: "period-output",
  "pulse-width": "pulse-width-output",
  promoter: "promoter-output",
  "osc-rbs": "osc-rbs-output",
  c31: "c31-output",
  degradation: "degradation-output",
  rdf: "rdf-output",
  bm3r1: "bm3r1-output",
  srna: "srna-output",
};

const COLORS = {
  integrase: "#19788a",
  rdf: "#87549a",
  input: "#bd8b2e",
  pb: "#356f91",
  lr: "#d26b35",
  output: "#b13f68",
  srna: "#5b8a45",
};

const CHART_FORMATTERS = {
  concentration: (value) => `${value.toFixed(3)} μM`,
  rate: (value) => `${value.toFixed(3)} μM/h`,
  fraction: (value) => `${(value * 100).toFixed(1)}%`,
};

const PRESET_DESCRIPTIONS = {
  custom: "Current controls differ from a named preset.",
  tempo: "Certified K4 baseline with the reference counter and shutdown settings.",
  stable: "K5 with tag 12 h⁻¹, tuned for slower pulses and robust state recovery.",
  fast: "K3 with stronger C31 input and faster clearance for shorter switching cycles.",
  long: "Reference counting with weaker shutdown to extend the output window.",
  short: "Reference counting with stronger shutdown to close the output window sooner.",
};

const SCORE_EXPLANATIONS = {
  deterministic: "Reference only — not calibrated against experimental outcomes. Heuristic ODE score: 55% complete one-flip rate + 20% alternation fidelity + 12% inter-pulse Integrase recovery + 13% output closure. The pass/fail verdict still requires every core check to pass.",
  surrogate: "Reference only — browser approximation, not an experimental prediction. Preview score: 47% predicted flip fidelity + 20% recovery margin + 12% leakage control + 21% output closure. Start the local server for the deterministic ODE score.",
};

let params = { ...DEFAULTS };
let lastResult;
let framePending = false;
let serverAvailable = false;
let odeTimer;
let requestRevision = 0;
let odeInFlight = false;
let odeQueued = false;
let designDirty = false;

const byId = (id) => document.getElementById(id);

function setRangeFill(input) {
  const min = Number(input.min);
  const max = Number(input.max);
  const percentage = ((Number(input.value) - min) / (max - min)) * 100;
  input.style.setProperty("--pct", `${percentage}%`);
}

function syncControls() {
  CONTROL_IDS.forEach((id) => {
    const input = byId(id);
    const value = params[CONTROL_KEYS[id]];
    input.value = String(value);
    byId(OUTPUT_IDS[id]).textContent = FORMATTERS[id](value);
    setRangeFill(input);
  });
  document.querySelectorAll("[data-period-setting]").forEach((button) => {
    const selected = Number(button.dataset.periodSetting) === Math.round(params.period);
    button.classList.toggle("active", selected);
    button.setAttribute("aria-pressed", String(selected));
  });
}

function setPresetLabel(name) {
  byId("preset-select").value = name;
  byId("preset-description").textContent = PRESET_DESCRIPTIONS[name];
}

function controlContainer(target) {
  if (target === "period-library") return document.querySelector(".period-picker");
  return byId(target)?.closest(".range-field, .select-field") ?? byId(target);
}

function buildTuningSuggestions(evaluation, metadata = {}) {
  const suggestions = [];
  const add = (target, title, copy) => {
    if (!suggestions.some((item) => item.target === target)) suggestions.push({ target, title, copy });
  };

  if (metadata.period_knob_base_conditions === false) {
    add("preset-select", "Restore the certified baseline", "Choose TEMPO reference to restore K4, C31 0.45× and the frozen expression settings before fine-tuning.");
  } else if (metadata.period_knob_ab_certified_for_tag === false) {
    const passingTags = metadata.period_knob_passing_tags_h_inv ?? [];
    if (passingTags.length) {
      add("degradation", "Match the certified Integrase tag", `Move degradation toward ${passingTags.join(" or ")} h⁻¹ for this K setting.`);
    } else {
      add("period-library", "Use a certified period setting", "Start from K4 with tag 8 or 12 h⁻¹, or K5 with tag 12 h⁻¹.");
    }
  }

  if (!evaluation.recovery) {
    add("degradation", "Clear Integrase before the next pulse", "Increase degradation in small steps; if flip fidelity drops, slightly reduce C31 RBS.");
  }

  if (evaluation.flipSuccess < 95) {
    if (params.degradation > 18) {
      add("degradation", "Slow an over-fast degradation tag", "Lower degradation so Integrase has enough time to complete each DNA flip.");
    } else {
      add("c31", "Recover one-flip fidelity", "Adjust C31 RBS in small steps around 0.42–0.52× and recheck complete flips.");
    }
  }

  if (evaluation.residual > 20) {
    add("srna", "Close the output more completely", "Increase sRNA production gradually; this directly targets excessive residual output.");
  }

  if (evaluation.outputDuration < 1) {
    add("srna", "Preserve a measurable output window", "Reduce sRNA production slightly so shutdown does not suppress the output too early.");
  }

  if (!suggestions.length) {
    add("preset-select", "Return to a known starting point", "Load TEMPO reference, then change one control at a time and wait for the ODE result.");
  }
  return suggestions.slice(0, 3);
}

function renderTuningGuide(evaluation, metadata) {
  const guide = byId("tuning-guide");
  const list = byId("tuning-suggestions");
  const toggle = byId("tuning-guide-toggle");
  document.querySelectorAll(".recommended-control").forEach((item) => item.classList.remove("recommended-control"));
  if (evaluation.verdict === "success") {
    guide.hidden = true;
    guide.classList.remove("compact", "expanded");
    list.replaceChildren();
    return;
  }

  const suggestions = buildTuningSuggestions(evaluation, metadata);
  const nodes = suggestions.map((suggestion, index) => {
    controlContainer(suggestion.target)?.classList.add("recommended-control");
    const button = document.createElement("button");
    button.type = "button";
    button.className = "tuning-tip";
    button.dataset.target = suggestion.target;
    button.innerHTML = `<span>${index + 1}</span><div><b>${suggestion.title}</b><small>${suggestion.copy}</small></div>`;
    return button;
  });
  list.replaceChildren(...nodes);
  const compact = evaluation.verdict === "warning";
  guide.classList.toggle("compact", compact);
  guide.classList.remove("expanded");
  toggle.hidden = !compact;
  toggle.textContent = "Show details";
  toggle.setAttribute("aria-expanded", "false");
  byId("tuning-guide-title").textContent = compact ? "Suggested adjustment" : "Suggested adjustments";
  guide.hidden = false;
}

function renderEvaluation(evaluation, metadata, source) {
  byId("verdict-title").textContent = {
    success: "Design passes",
    warning: "Needs attention",
    danger: "Design fails",
  }[evaluation.verdict];
  const countingPasses = evaluation.flipSuccess >= 95 && evaluation.recovery;
  const warningCopy = countingPasses && evaluation.residual > 20
    ? "Counting passes; residual output remains above the shutdown target."
    : "The model found a limited counting, timing or shutdown margin.";
  byId("verdict-copy").textContent = evaluation.verdict === "success"
    ? "Counting, recovery and shutdown checks pass."
    : evaluation.verdict === "warning"
      ? warningCopy
      : "A core counting or shutdown check fails.";
  byId("score-value").textContent = String(evaluation.score);
  byId("score-ring").style.setProperty("--score", evaluation.score);
  byId("score-ring").setAttribute("aria-label", `Reference-only heuristic score ${evaluation.score} out of 100`);
  byId("score-summary").className = `score-summary ${evaluation.verdict}`;
  byId("score-tooltip").textContent = SCORE_EXPLANATIONS[source] ?? SCORE_EXPLANATIONS.surrogate;
  byId("score-rule").hidden = evaluation.verdict === "success";
  byId("score-rule").textContent = evaluation.verdict === "warning"
    ? "The reference-only heuristic score can remain high even when one required core check misses its target."
    : "A failed core check determines the verdict even when the other heuristic metrics perform well.";

  const banner = byId("verdict-banner");
  banner.className = `verdict-banner ${evaluation.verdict}`;
  banner.querySelector(".verdict-icon").textContent = evaluation.verdict === "success" ? "✓" : evaluation.verdict === "warning" ? "!" : "×";

  byId("flip-success").textContent = `${Math.round(evaluation.flipSuccess)}%`;

  const recovery = byId("recovery-status");
  recovery.textContent = evaluation.recovery ? "Yes" : "No";
  recovery.className = `status ${evaluation.recovery ? "success" : "danger"}`;

  byId("duration-time").textContent = `${evaluation.outputDuration.toFixed(1)} h`;
  byId("residual-value").textContent = `${Math.round(evaluation.residual)}%`;

  const priorityInsights = evaluation.insights.slice(0, 2);
  byId("diagnostic-count").textContent = `${priorityInsights.length} ${priorityInsights.length === 1 ? "note" : "notes"}`;
  byId("diagnostic-items").replaceChildren(...priorityInsights.map((insight) => {
    const item = document.createElement("p");
    item.className = `diagnostic-item ${insight.level}`;
    item.textContent = insight.text;
    return item;
  }));
  renderTuningGuide(evaluation, metadata);
}

function setModelStatus(mode, metadata = null) {
  const status = byId("model-status");
  const facts = byId("model-facts");
  const scope = byId("model-scope-copy");
  const certificationCard = byId("certification-card");
  const certificationStatus = byId("certification-status");
  const indicator = byId("update-indicator");
  const runButton = byId("run-model-button");
  const updating = mode === "loading" || mode === "queued";
  document.body.classList.toggle("model-updating", updating);
  byId("workspace").setAttribute("aria-busy", String(updating));
  byId("score-summary").setAttribute("aria-busy", String(updating));
  indicator.hidden = !["pending", "queued", "loading", "error"].includes(mode);
  runButton.hidden = !["pending", "error"].includes(mode);
  if (mode === "pending") {
    status.textContent = "Parameters changed · result not updated";
    indicator.textContent = odeInFlight
      ? "Current solve is finishing · latest changes are pending"
      : "Parameters changed · release the control or run the latest design";
    scope.innerHTML = "<b>Preview only.</b> The controls have changed, but these values are not yet a completed Python ODE result.";
    return;
  }
  if (mode === "queued") {
    const knob = getPeriodKnob(params.period);
    status.textContent = `${knob.id} queued · waiting for current solve`;
    indicator.textContent = `Latest ${knob.id} design queued · it will run next`;
    scope.innerHTML = "<b>Queued.</b> The current solve must finish before the latest controls can enter the SciPy model.";
    return;
  }
  if (mode === "loading") {
    const knob = getPeriodKnob(params.period);
    status.textContent = `Solving ${knob.id} coupled ODE…`;
    indicator.textContent = `Solving ${knob.id} locally with Python/SciPy · the full coupled ODE may take 1–2 min`;
    certificationStatus.textContent = "Checking current design…";
    certificationCard.className = "certification-card";
    scope.innerHTML = "<b>Computing.</b> The current controls are being solved by the coupled Python model.";
    return;
  }
  if (mode === "deterministic") {
    status.textContent = `Deterministic ODE · ${metadata?.period_setting_id ?? "latest design"} updated`;
    facts.hidden = false;
    byId("realized-period").textContent = Number.isFinite(metadata?.realized_period_h)
      ? `${metadata?.period_setting_id ?? ""} · ${metadata.realized_period_h.toFixed(2)} h`
      : "Unresolved";
    byId("realized-width").textContent = Number.isFinite(metadata?.realized_pulse_width_h)
      ? `${metadata.realized_pulse_width_h.toFixed(2)} h`
      : "Unresolved";
    if (metadata?.period_knob_ab_certified_for_tag) {
      certificationStatus.textContent = "Certified for this tag";
      certificationCard.className = "certification-card success";
    } else if (metadata?.period_knob_base_conditions === false) {
      certificationStatus.textContent = "Outside frozen conditions";
      certificationCard.className = "certification-card danger";
    } else if (metadata?.period_knob_passing_tags_h_inv?.length) {
      certificationStatus.textContent = `Use tag ${metadata.period_knob_passing_tags_h_inv.join(" / ")} h⁻¹`;
      certificationCard.className = "certification-card";
    } else {
      certificationStatus.textContent = "Oscillator only · B not certified";
      certificationCard.className = "certification-card danger";
    }
    scope.innerHTML = "<b>Deterministic model.</b> These curves were solved from the existing 8-state oscillator, 38-state counter, and 3-state shutdown equations. The diagnostic states separately whether the selected K setting lies inside the certified A→B window.";
    return;
  }
  if (mode === "error") {
    status.textContent = "Model update failed · retry available";
    indicator.textContent = "The latest controls were not solved · run the latest design again";
    scope.innerHTML = "<b>Update failed.</b> The visible curves are a browser preview until the deterministic request succeeds.";
    return;
  }
  status.textContent = "Browser surrogate · fallback";
  facts.hidden = true;
  scope.innerHTML = "<b>Fallback preview.</b> Start the local TEMPO server to replace these approximate curves with the coupled Python ODE solution.";
}

function renderResult(result, source) {
  lastResult = result;
  const { samples, centers, evaluation, summary } = result;

  renderChart(byId("protein-chart"), samples, [
    { label: "Integrase", accessor: (d) => d.integrase, color: COLORS.integrase, area: true, format: CHART_FORMATTERS.concentration },
    { label: "RDF", accessor: (d) => d.rdf, color: COLORS.rdf, format: CHART_FORMATTERS.concentration },
    { label: "C31 input", accessor: (d) => d.c31, color: COLORS.input, dash: "4 4", opacity: 0.9, format: CHART_FORMATTERS.rate },
  ], centers);
  renderChart(byId("counter-chart"), samples, [
    { label: "PB fraction", accessor: (d) => d.pb, color: COLORS.pb, area: true, format: CHART_FORMATTERS.fraction },
    { label: "LR fraction", accessor: (d) => d.lr, color: COLORS.lr, format: CHART_FORMATTERS.fraction },
  ], centers);
  renderChart(byId("output-chart"), samples, [
    { label: "Output", accessor: (d) => d.output, color: COLORS.output, area: true, format: CHART_FORMATTERS.concentration },
    { label: "sRNA", accessor: (d) => d.srna, color: COLORS.srna, dash: "6 3", format: CHART_FORMATTERS.concentration },
  ], centers);
  renderStateRail(byId("state-rail"), samples, centers, params.hours);

  byId("pulse-count").textContent = String(summary.pulses);
  byId("flip-count").textContent = String(summary.flips);
  byId("output-count").textContent = String(summary.outputs);
  renderEvaluation(evaluation, result.metadata, source);
  setModelStatus(source, result.metadata);
}

function renderSurrogate() {
  framePending = false;
  renderResult(simulateTempo(params), "surrogate");
  if (serverAvailable && designDirty) setModelStatus("pending");
}

function schedulePreview() {
  if (framePending) return;
  framePending = true;
  requestAnimationFrame(renderSurrogate);
}

async function requestOde() {
  if (!serverAvailable) return;
  if (odeInFlight) {
    odeQueued = true;
    setModelStatus("queued");
    return;
  }
  const revision = requestRevision;
  const requestedParams = { ...params };
  let requestFailed = false;
  odeInFlight = true;
  setModelStatus("loading");
  try {
    const response = await fetch("/api/simulate", {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      cache: "no-store",
      body: JSON.stringify(requestedParams),
    });
    if (!response.ok) throw new Error(`ODE request failed: ${response.status}`);
    const result = await response.json();
    if (revision !== requestRevision) return;
    designDirty = false;
    renderResult(result, "deterministic");
  } catch (error) {
    if (revision !== requestRevision) return;
    console.warn(error);
    requestFailed = true;
    designDirty = true;
    setModelStatus("error");
  } finally {
    odeInFlight = false;
    if (odeQueued) {
      odeQueued = false;
      requestOde();
    } else if (designDirty && !requestFailed) {
      setModelStatus("pending");
    }
  }
}

function scheduleOde(delay = 180) {
  clearTimeout(odeTimer);
  if (!serverAvailable) return;
  if (odeInFlight) {
    odeQueued = true;
    setModelStatus("queued");
    return;
  }
  odeTimer = setTimeout(requestOde, delay);
}

function markDesignChanged() {
  requestRevision += 1;
  designDirty = true;
  clearTimeout(odeTimer);
  schedulePreview();
  if (serverAvailable) setModelStatus("pending");
}

function scheduleRender() {
  markDesignChanged();
  scheduleOde();
}

async function connectModelServer() {
  try {
    const response = await fetch("/api/health", { cache: "no-store" });
    if (!response.ok) throw new Error("Model server unavailable");
    const health = await response.json();
    serverAvailable = health.source === "deterministic_ode";
  } catch {
    serverAvailable = false;
  }
  if (serverAvailable) requestOde();
  else setModelStatus("surrogate");
}

function setPreset(name) {
  params = { ...PRESETS[name], hours: params.hours };
  setPresetLabel(name);
  syncControls();
  scheduleRender();
}

function registerWebMcpTools() {
  const context = document.modelContext;
  if (!context?.registerTool) return;
  const lifecycle = new AbortController();
  const numericProperties = {
    periodSetting: { key: "period", minimum: 1, maximum: 5, integer: true },
    pulseWidthHours: { key: "pulseWidth", minimum: 1.5, maximum: 5.5 },
    promoterScale: { key: "promoter", minimum: 0.5, maximum: 1.5 },
    oscillatorRbsMultiplier: { key: "oscRbs", minimum: 0.5, maximum: 1.5 },
    c31RbsScale: { key: "c31", minimum: 0.3, maximum: 0.8 },
    integraseTagPerHour: { key: "degradation", minimum: 2.8, maximum: 28 },
    rdfScale: { key: "rdf", minimum: 0.4, maximum: 1.7 },
    bm3r1Scale: { key: "bm3r1", minimum: 0.4, maximum: 1.8 },
    srnaProduction: { key: "srna", minimum: 0.02, maximum: 0.1 },
    durationHours: { key: "hours", minimum: 36, maximum: 84 },
  };

  const register = (tool) => {
    try {
      void Promise.resolve(context.registerTool(tool, { signal: lifecycle.signal }))
        .catch((error) => console.warn("WebMCP registration failed", error));
    } catch (error) {
      console.warn("WebMCP registration failed", error);
    }
  };

  register({
    name: "stage_tempo_design",
    title: "Stage TEMPO design",
    description: "Set one or more visible TEMPO controls and start the same deterministic ODE workflow used by the page. Returns the staged design immediately; use read_tempo_result to inspect completion.",
    inputSchema: {
      type: "object",
      properties: Object.fromEntries(Object.entries(numericProperties).map(([name, spec]) => [name, {
        type: spec.integer ? "integer" : "number",
        minimum: spec.minimum,
        maximum: spec.maximum,
      }])),
      additionalProperties: false,
    },
    annotations: { readOnlyHint: false, untrustedContentHint: false },
    execute(input) {
      if (!input || typeof input !== "object" || Array.isArray(input)) throw new TypeError("Input must be an object.");
      for (const [name, value] of Object.entries(input)) {
        const spec = numericProperties[name];
        if (!spec) throw new RangeError(`Unknown parameter: ${name}`);
        if (!Number.isFinite(value) || value < spec.minimum || value > spec.maximum || (spec.integer && !Number.isInteger(value))) {
          throw new RangeError(`${name} must be within ${spec.minimum}–${spec.maximum}.`);
        }
        params[spec.key] = value;
      }
      setPresetLabel("custom");
      syncControls();
      renderSurrogate();
      scheduleOde();
      const knob = getPeriodKnob(params.period);
      return { status: "staged", periodSetting: knob.id, nominalPeriodHours: knob.nominalHours, design: { ...params } };
    },
  });

  register({
    name: "read_tempo_result",
    title: "Read TEMPO result",
    description: "Read the current visible TEMPO design, model provenance, deterministic evaluation, and pulse/flip/output counts.",
    inputSchema: { type: "object", properties: {}, additionalProperties: false },
    annotations: { readOnlyHint: true, untrustedContentHint: false },
    execute(input) {
      if (!input || typeof input !== "object" || Array.isArray(input) || Object.keys(input).length) {
        throw new TypeError("Input must be an empty object.");
      }
      return {
        design: { ...params },
        provenance: lastResult?.metadata ?? { source: "pending" },
        evaluation: lastResult?.evaluation ?? null,
        summary: lastResult?.summary ?? null,
      };
    },
  });
}

CONTROL_IDS.forEach((id) => {
  const input = byId(id);
  input.addEventListener("input", () => {
    const value = Number(input.value);
    params[CONTROL_KEYS[id]] = value;
    byId(OUTPUT_IDS[id]).textContent = FORMATTERS[id](value);
    setRangeFill(input);
    setPresetLabel("custom");
    markDesignChanged();
  });
  input.addEventListener("change", () => scheduleOde());
});

document.querySelectorAll(".segmented button").forEach((button) => {
  button.addEventListener("click", () => {
    document.querySelectorAll(".segmented button").forEach((item) => item.classList.toggle("active", item === button));
    const advanced = button.dataset.mode === "advanced";
    document.querySelectorAll(".advanced-only").forEach((section) => { section.hidden = !advanced; });
    byId("mode-description").textContent = advanced
      ? "Adds pulse-shape and expression multipliers. Changing them moves beyond the frozen certification conditions."
      : "Four core decisions: period, C31 input, Integrase reset and shutdown strength.";
  });
});

document.querySelectorAll("[data-period-setting]").forEach((button) => {
  button.addEventListener("click", () => {
    params.period = Number(button.dataset.periodSetting);
    setPresetLabel("custom");
    syncControls();
    scheduleRender();
  });
});

document.querySelectorAll(".time-window button").forEach((button) => {
  button.addEventListener("click", () => {
    document.querySelectorAll(".time-window button").forEach((item) => item.classList.toggle("active", item === button));
    params.hours = Number(button.dataset.hours);
    scheduleRender();
  });
});

byId("tuning-suggestions").addEventListener("click", (event) => {
  const tip = event.target.closest(".tuning-tip");
  if (!tip) return;
  const target = controlContainer(tip.dataset.target);
  const focusTarget = tip.dataset.target === "period-library"
    ? document.querySelector("[data-period-setting].active")
    : byId(tip.dataset.target);
  target?.scrollIntoView({ behavior: "smooth", block: "center" });
  focusTarget?.focus({ preventScroll: true });
  target?.classList.add("guidance-focus");
  window.setTimeout(() => target?.classList.remove("guidance-focus"), 1600);
});

byId("tuning-guide-toggle").addEventListener("click", () => {
  const guide = byId("tuning-guide");
  const expanded = guide.classList.toggle("expanded");
  byId("tuning-guide-toggle").textContent = expanded ? "Hide details" : "Show details";
  byId("tuning-guide-toggle").setAttribute("aria-expanded", String(expanded));
});

byId("preset-select").addEventListener("change", (event) => setPreset(event.target.value));
byId("reset-button").addEventListener("click", () => {
  setPreset("tempo");
});
byId("run-model-button").addEventListener("click", () => {
  if (serverAvailable) scheduleOde(0);
  else connectModelServer();
});

const aboutDialog = byId("about-dialog");
byId("about-button").addEventListener("click", () => aboutDialog.showModal());

byId("export-button").addEventListener("click", () => {
  const isDeterministic = lastResult?.metadata?.source === "deterministic_ode";
  const payload = {
    exported_at: new Date().toISOString(),
    result_source: isDeterministic ? "deterministic_ode" : "browser_surrogate",
    model_metadata: lastResult?.metadata ?? null,
    scientific_notice: isDeterministic
      ? "Solved with the existing TEMPO Python ODE pipeline. Deterministic A→B certification applies only where the exported metadata explicitly marks it true."
      : "Interactive calibrated surrogate; start the local model server for reportable ODE output.",
    design_parameters: { ...params },
    model_parameter_mapping: PARAMETER_MAP,
    evaluation: lastResult.evaluation,
    summary: lastResult.summary,
  };
  const blob = new Blob([JSON.stringify(payload, null, 2)], { type: "application/json" });
  const url = URL.createObjectURL(blob);
  const anchor = document.createElement("a");
  anchor.href = url;
  anchor.download = "tempo-design-summary.json";
  anchor.click();
  URL.revokeObjectURL(url);
});

syncControls();
setPresetLabel("tempo");
renderSurrogate();
connectModelServer();
registerWebMcpTools();
