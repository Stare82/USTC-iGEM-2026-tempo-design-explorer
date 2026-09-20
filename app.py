"""User-facing launcher for the reproducible TEMPO Design Explorer package."""

from __future__ import annotations

from pathlib import Path
import sys


MINIMUM_PYTHON = (3, 11)
ROOT = Path(__file__).resolve().parent
WEB_APP = ROOT / "tempo-design-explorer"


def _check_package_layout() -> None:
    required = (
        WEB_APP / "server.py",
        WEB_APP / "index.html",
        WEB_APP / "styles.css",
        WEB_APP / "src" / "app.js",
        WEB_APP / "backend" / "model_adapter.py",
        WEB_APP / "public" / "tempo-logo.png",
        WEB_APP / "docs" / "PLATFORM_ARCHITECTURE.md",
        ROOT / "oscillator" / "code" / "Mechanistic_ODE_Global_Sensitivity_Analysis.py",
        ROOT / "oscillator" / "code" / "Shared_PLtetO1_Period_Knob_Design_Map.py",
        ROOT / "rdfmodel_new" / "model" / "zhao_core.py",
        ROOT / "rdfmodel_new" / "docs" / "week4_unloaded" / "unloaded_C31_translation_trajectory.csv",
        ROOT / "shutdown_model" / "shutdown_core.py",
    )
    missing = [path.relative_to(ROOT) for path in required if not path.is_file()]
    if missing:
        paths = "\n".join(f"  - {path}" for path in missing)
        raise SystemExit(
            "TEMPO Design Explorer is incomplete. Re-download and extract the full ZIP.\n"
            f"Missing files:\n{paths}"
        )


def main() -> None:
    if sys.version_info < MINIMUM_PYTHON:
        found = ".".join(map(str, sys.version_info[:3]))
        raise SystemExit(
            f"Python {MINIMUM_PYTHON[0]}.{MINIMUM_PYTHON[1]} or newer is required "
            f"(found {found}). Python 3.13.9 is the tested release version."
        )

    _check_package_layout()
    sys.path.insert(0, str(WEB_APP))
    try:
        from server import main as run_server
    except ModuleNotFoundError as error:
        if error.name in {"numpy", "scipy", "pandas", "matplotlib"}:
            raise SystemExit(
                f"Missing Python dependency: {error.name}\n"
                "Please run:\n  python -m pip install -r requirements.txt"
            ) from error
        raise

    run_server(default_open=True)


if __name__ == "__main__":
    main()
