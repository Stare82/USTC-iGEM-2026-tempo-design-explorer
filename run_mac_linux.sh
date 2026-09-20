#!/usr/bin/env sh
set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
cd "$SCRIPT_DIR"

if ! command -v python3 >/dev/null 2>&1; then
  echo "Python 3 was not found. Install Python 3.13.9 and try again." >&2
  exit 1
fi

if ! python3 -c 'import sys; raise SystemExit(0 if sys.version_info >= (3, 11) else 1)' >/dev/null 2>&1; then
  echo "Python 3.11 or newer is required. Python 3.13.9 is recommended." >&2
  python3 --version >&2
  exit 1
fi

if [ ! -x ".venv/bin/python" ]; then
  echo "Creating the TEMPO virtual environment..."
  python3 -m venv .venv
fi

if ! .venv/bin/python -c "import numpy, scipy, pandas, matplotlib" >/dev/null 2>&1; then
  echo "Installing the tested TEMPO dependencies..."
  .venv/bin/python -m pip install -r requirements.txt
fi

echo "Starting TEMPO Design Explorer..."
exec .venv/bin/python app.py
