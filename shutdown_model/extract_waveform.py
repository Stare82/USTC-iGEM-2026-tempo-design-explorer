"""Extract the oscillator C31 trajectory used to drive the shutdown model."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SOURCE = (
    ROOT
    / "rdfmodel_new"
    / "docs"
    / "week4_unloaded"
    / "unloaded_C31_translation_trajectory.csv"
)
DATA_DIR = Path(__file__).resolve().parent / "data"
OUTPUT = DATA_DIR / "c31_waveform.csv"
METADATA = DATA_DIR / "c31_waveform_metadata.json"
COPIES_PER_UM = 602.0


def extract_waveform() -> Path:
    """Create a compact C31 waveform in hours and uM."""
    if not SOURCE.exists():
        raise FileNotFoundError(f"Oscillator trajectory not found: {SOURCE}")

    df = pd.read_csv(SOURCE)
    required = {"time_min", "C31_protein"}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Oscillator trajectory is missing columns: {sorted(missing)}")

    out = pd.DataFrame(
        {
            "time_h": df["time_min"].astype(float) / 60.0,
            "c31_uM": df["C31_protein"].astype(float) / COPIES_PER_UM,
        }
    )
    flux_column = "actual_C31_translation_flux_uM_h"
    if flux_column not in df.columns:
        raise ValueError(
            "Oscillator trajectory is missing the physical A-to-B interface: "
            f"{flux_column}"
        )
    out["translation_flux_uM_h"] = df[flux_column].astype(float)

    if not out["time_h"].is_monotonic_increasing:
        raise ValueError("Oscillator time points must be monotonically increasing")
    if out["time_h"].duplicated().any():
        raise ValueError("Oscillator trajectory contains duplicate time points")

    DATA_DIR.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUTPUT, index=False)
    metadata = {
        "source": str(SOURCE),
        "copies_per_uM": COPIES_PER_UM,
        "rows": int(len(out)),
        "start_h": float(out["time_h"].iloc[0]),
        "end_h": float(out["time_h"].iloc[-1]),
        "c31_min_uM": float(out["c31_uM"].min()),
        "c31_max_uM": float(out["c31_uM"].max()),
        "translation_flux_min_uM_h": float(out["translation_flux_uM_h"].min()),
        "translation_flux_max_uM_h": float(out["translation_flux_uM_h"].max()),
        "note": (
            "Frozen input distributed with rdfmodel_new. C31 protein concentration "
            "is retained for diagnostics only; the counter is driven by the "
            "translation-production flux."
        ),
    }
    METADATA.write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return OUTPUT


if __name__ == "__main__":
    path = extract_waveform()
    print(path)
