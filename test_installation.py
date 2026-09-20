"""Verify dependencies, package layout, and the coupled reference simulation."""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parent
APP_ROOT = ROOT / "tempo-design-explorer"
EXAMPLE = ROOT / "examples" / "baseline.json"
os.environ.setdefault("MPLBACKEND", "Agg")


def run_installation_test() -> dict[str, object]:
    if sys.version_info < (3, 11):
        raise RuntimeError("Python 3.11 or newer is required; Python 3.13.9 is recommended.")

    try:
        import matplotlib
        import numpy
        import pandas
        import scipy
    except ModuleNotFoundError as error:
        raise RuntimeError(
            f"Missing dependency {error.name!r}. Run: python -m pip install -r requirements.txt"
        ) from error

    required = (
        APP_ROOT / "server.py",
        APP_ROOT / "index.html",
        APP_ROOT / "styles.css",
        APP_ROOT / "src" / "app.js",
        APP_ROOT / "backend" / "model_adapter.py",
        APP_ROOT / "public" / "tempo-logo.png",
        ROOT / "oscillator" / "code" / "Mechanistic_ODE_Global_Sensitivity_Analysis.py",
        ROOT / "oscillator" / "code" / "Shared_PLtetO1_Period_Knob_Design_Map.py",
        ROOT / "rdfmodel_new" / "model" / "zhao_core.py",
        ROOT / "shutdown_model" / "shutdown_core.py",
        EXAMPLE,
    )
    missing = [str(path.relative_to(ROOT)) for path in required if not path.is_file()]
    if missing:
        raise RuntimeError(f"Package files are missing: {', '.join(missing)}")

    sys.path.insert(0, str(APP_ROOT))
    from backend.model_adapter import MODEL_VERSION, simulate_design

    design = json.loads(EXAMPLE.read_text(encoding="utf-8"))
    result = simulate_design(design)
    metadata = result["metadata"]
    evaluation = result["evaluation"]

    assert metadata["source"] == "deterministic_ode"
    assert metadata["solver"] == "SciPy LSODA"
    assert metadata["model_version"] == MODEL_VERSION
    assert metadata["states"] == {"oscillator": 8, "counter": 38, "shutdown": 3}
    assert len(result["samples"]) == 900
    assert result["summary"]["pulses"] >= 1
    assert 0 <= evaluation["score"] <= 100

    return {
        "python": ".".join(map(str, sys.version_info[:3])),
        "numpy": numpy.__version__,
        "scipy": scipy.__version__,
        "pandas": pandas.__version__,
        "matplotlib": matplotlib.__version__,
        "model_version": MODEL_VERSION,
        "solver": metadata["solver"],
        "samples": len(result["samples"]),
        "score": evaluation["score"],
        "verdict": evaluation["verdict"],
    }


def main() -> None:
    try:
        details = run_installation_test()
    except Exception as error:
        print(f"TEMPO Design Explorer installation test failed: {error}", file=sys.stderr)
        raise SystemExit(1) from error

    print(json.dumps(details, indent=2, ensure_ascii=False))
    print("TEMPO Design Explorer installation test passed.")


if __name__ == "__main__":
    main()
