const SVG_NS = "http://www.w3.org/2000/svg";

function svgElement(name, attributes = {}) {
  const node = document.createElementNS(SVG_NS, name);
  Object.entries(attributes).forEach(([key, value]) => node.setAttribute(key, value));
  return node;
}

function pathFor(samples, accessor, maxValue, width, height, pad) {
  return samples.map((sample, index) => {
    const x = pad.left + (sample.t / samples.at(-1).t) * (width - pad.left - pad.right);
    const y = pad.top + (1 - accessor(sample) / maxValue) * (height - pad.top - pad.bottom);
    return `${index ? "L" : "M"}${x.toFixed(2)},${y.toFixed(2)}`;
  }).join(" ");
}

function areaFor(samples, accessor, maxValue, width, height, pad) {
  const line = pathFor(samples, accessor, maxValue, width, height, pad);
  const xEnd = width - pad.right;
  const baseline = height - pad.bottom;
  return `${line} L${xEnd},${baseline} L${pad.left},${baseline} Z`;
}

function nearestSample(samples, targetTime) {
  let low = 0;
  let high = samples.length - 1;
  while (low < high) {
    const middle = Math.floor((low + high) / 2);
    if (samples[middle].t < targetTime) low = middle + 1;
    else high = middle;
  }
  if (low === 0) return samples[0];
  const previous = samples[low - 1];
  return Math.abs(previous.t - targetTime) <= Math.abs(samples[low].t - targetTime) ? previous : samples[low];
}

export function renderChart(container, samples, series, centers = []) {
  const width = 760;
  const height = 176;
  const pad = { top: 10, right: 12, bottom: 25, left: 36 };
  const svg = svgElement("svg", { viewBox: `0 0 ${width} ${height}`, preserveAspectRatio: "none", "aria-hidden": "true" });
  const maxValue = Math.max(1e-6, ...series.flatMap((item) => samples.map(item.accessor)));
  const plotHeight = height - pad.top - pad.bottom;

  [0, 0.5, 1].forEach((fraction) => {
    const y = pad.top + (1 - fraction) * plotHeight;
    svg.append(svgElement("line", { x1: pad.left, y1: y, x2: width - pad.right, y2: y, class: "grid-line" }));
    const label = svgElement("text", { x: pad.left - 7, y: y + 3, class: "axis-label", "text-anchor": "end" });
    label.textContent = fraction === 1 ? "max" : fraction.toFixed(1);
    svg.append(label);
  });

  const hours = samples.at(-1).t;
  const tickCount = hours <= 36 ? 6 : hours <= 60 ? 5 : 7;
  for (let i = 0; i <= tickCount; i += 1) {
    const fraction = i / tickCount;
    const x = pad.left + fraction * (width - pad.left - pad.right);
    const label = svgElement("text", { x, y: height - 6, class: "axis-label", "text-anchor": i === 0 ? "start" : i === tickCount ? "end" : "middle" });
    label.textContent = `${Math.round(hours * fraction)} h`;
    svg.append(label);
  }

  centers.forEach((center) => {
    const x = pad.left + (center / hours) * (width - pad.left - pad.right);
    svg.append(svgElement("line", { x1: x, y1: pad.top, x2: x, y2: height - pad.bottom, class: "pulse-line" }));
  });

  series.forEach((item) => {
    if (item.area) {
      svg.append(svgElement("path", {
        d: areaFor(samples, item.accessor, maxValue, width, height, pad),
        fill: item.color,
        class: "chart-area",
      }));
    }
    svg.append(svgElement("path", {
      d: pathFor(samples, item.accessor, maxValue, width, height, pad),
      stroke: item.color,
      class: "chart-path",
      opacity: item.opacity ?? 1,
      "stroke-dasharray": item.dash ?? "",
    }));
  });

  const hoverLine = svgElement("line", {
    x1: pad.left,
    y1: pad.top,
    x2: pad.left,
    y2: height - pad.bottom,
    class: "chart-hover-line",
    visibility: "hidden",
  });
  svg.append(hoverLine);

  const markers = series.map((item) => {
    const marker = svgElement("circle", {
      cx: pad.left,
      cy: height - pad.bottom,
      r: 3.2,
      fill: item.color,
      class: "chart-hover-marker",
      visibility: "hidden",
    });
    svg.append(marker);
    return marker;
  });

  const hitArea = svgElement("rect", {
    x: pad.left,
    y: pad.top,
    width: width - pad.left - pad.right,
    height: height - pad.top - pad.bottom,
    class: "chart-hit-area",
  });
  svg.append(hitArea);

  const tooltip = document.createElement("div");
  tooltip.className = "chart-tooltip";
  tooltip.hidden = true;

  const hideTooltip = () => {
    tooltip.hidden = true;
    hoverLine.setAttribute("visibility", "hidden");
    markers.forEach((marker) => marker.setAttribute("visibility", "hidden"));
  };

  hitArea.addEventListener("pointermove", (event) => {
    const svgRect = svg.getBoundingClientRect();
    const containerRect = container.getBoundingClientRect();
    const svgX = Math.min(width - pad.right, Math.max(pad.left, ((event.clientX - svgRect.left) / svgRect.width) * width));
    const fraction = (svgX - pad.left) / (width - pad.left - pad.right);
    const sample = nearestSample(samples, fraction * hours);
    const sampleX = pad.left + (sample.t / hours) * (width - pad.left - pad.right);

    hoverLine.setAttribute("x1", sampleX);
    hoverLine.setAttribute("x2", sampleX);
    hoverLine.setAttribute("visibility", "visible");
    markers.forEach((marker, index) => {
      const value = series[index].accessor(sample);
      const y = pad.top + (1 - value / maxValue) * plotHeight;
      marker.setAttribute("cx", sampleX);
      marker.setAttribute("cy", y);
      marker.setAttribute("visibility", "visible");
    });

    const heading = document.createElement("strong");
    heading.textContent = `Time ${sample.t.toFixed(2)} h`;
    const rows = series.map((item) => {
      const row = document.createElement("div");
      row.className = "chart-tooltip-row";
      const dot = document.createElement("i");
      dot.style.background = item.color;
      const label = document.createElement("b");
      label.textContent = item.label ?? "Series";
      const value = document.createElement("span");
      const rawValue = item.accessor(sample);
      value.textContent = item.format ? item.format(rawValue) : rawValue.toFixed(3);
      row.append(dot, label, value);
      return row;
    });
    tooltip.replaceChildren(heading, ...rows);
    tooltip.hidden = false;

    const screenX = svgRect.left + (sampleX / width) * svgRect.width - containerRect.left;
    const tooltipWidth = Math.min(210, Math.max(168, tooltip.offsetWidth));
    tooltip.style.left = `${Math.max(8, Math.min(containerRect.width - tooltipWidth - 8, screenX + 12))}px`;
    tooltip.style.top = "8px";
  });
  svg.addEventListener("pointerleave", hideTooltip);

  container.replaceChildren(svg, tooltip);
}

export function renderStateRail(container, samples, centers, hours) {
  const boundaries = [0, ...centers.map((center) => Math.min(hours, center + 0.5)), hours]
    .filter((value, index, values) => index === 0 || value > values[index - 1]);
  const fragments = [];
  for (let index = 0; index < boundaries.length - 1; index += 1) {
    const left = boundaries[index];
    const right = boundaries[index + 1];
    if (right - left < 0.1) continue;
    const midpoint = (left + right) / 2;
    const sample = samples.reduce((best, current) => Math.abs(current.t - midpoint) < Math.abs(best.t - midpoint) ? current : best);
    const state = sample.lr > 0.7 ? "LR" : sample.pb > 0.7 ? "PB" : "mixed";
    const span = document.createElement("span");
    span.className = state.toLowerCase();
    span.style.flexGrow = right - left;
    span.textContent = state === "mixed" ? "MIX" : state;
    fragments.push(span);
  }
  container.replaceChildren(...fragments);
}
