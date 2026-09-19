"""Small deterministic deployment check using the public reference design."""

from __future__ import annotations

import json
from pathlib import Path
import sys


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
APP_ROOT = REPOSITORY_ROOT / "tempo-design-explorer"
sys.path.insert(0, str(APP_ROOT))

from backend.model_adapter import MODEL_VERSION, simulate_design  # noqa: E402


REFERENCE_DESIGN = {
    "hours": 36,
    "period": 4,
    "pulseWidth": 3.6,
    "promoter": 1.0,
    "oscRbs": 1.0,
    "c31": 0.45,
    "degradation": 8.0,
    "rdf": 1.0,
    "bm3r1": 1.0,
    "srna": 0.06,
}


def main() -> None:
    result = simulate_design(REFERENCE_DESIGN)
    metadata = result["metadata"]
    evaluation = result["evaluation"]

    assert metadata["source"] == "deterministic_ode"
    assert metadata["model_version"] == MODEL_VERSION
    assert metadata["states"] == {"oscillator": 8, "counter": 38, "shutdown": 3}
    assert len(result["samples"]) > 100
    assert 0 <= evaluation["score"] <= 100
    assert evaluation["verdict"] in {"success", "warning", "danger"}
    assert result["summary"]["pulses"] >= 1

    print(json.dumps({
        "status": "ok",
        "model_version": MODEL_VERSION,
        "score": evaluation["score"],
        "verdict": evaluation["verdict"],
        "samples": len(result["samples"]),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()

