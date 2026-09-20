# TEMPO Design Explorer

TEMPO Design Explorer is the reproducible local interface for the USTC iGEM
2026 TEMPO modeling workflow. It connects the existing deterministic models
into one design-to-behavior pipeline:

```text
8-state oscillator
        ↓ C31(t)
38-state recombinase counter
        ↓ DNA state
3-state timed-shutdown module
        ↓
system output
```

The browser provides parameter controls, synchronized plots, design evaluation,
and an exportable Design Summary. The calculations are performed locally by the
existing Python/SciPy reference implementation. The package does not replace the
models with a JavaScript approximation.

## What users can explore

Users can adjust:

- the characterized K1–K5 oscillator period setting;
- pulse width, shared promoter strength, and oscillator RBS multiplier;
- C31 expression and integrase degradation;
- RDF and BM3R1 expression;
- shutdown sRNA production;
- the displayed simulation window.

The interface reports:

- C31, Integrase, and RDF dynamics;
- PB/LR counter-state transitions;
- shutdown sRNA and output dynamics;
- pulse, flip, and output-window counts;
- flip fidelity, recovery, leakage, residual output, and output duration;
- the deterministic A→B certification status where applicable;
- a reference-only heuristic score and design-adjustment suggestions.

## Requirements

- **Tested and recommended Python:** 3.13.9
- **Minimum supported Python:** 3.11
- **Operating systems:** Windows, macOS, or Linux
- A modern browser

For the most reproducible environment, use Python 3.13.9 and the exact package
versions in [`requirements.txt`](requirements.txt).

## Download

Download the repository ZIP from the official USTC iGEM GitLab project, then
extract the complete archive. Do not copy only the `tempo-design-explorer`
folder: the application also needs the scientific source files in `oscillator`,
`rdfmodel_new`, and `shutdown_model`.

## Quick start

### Windows one-click launcher

Double-click:

```text
run_windows.bat
```

The script creates `.venv` when needed, installs the tested dependencies, starts
the application, and opens the local browser page.

### macOS or Linux launcher

From a terminal in the extracted folder:

```bash
chmod +x run_mac_linux.sh
./run_mac_linux.sh
```

You can also run it without changing permissions:

```bash
sh run_mac_linux.sh
```

## Manual installation

Create an isolated environment:

```bash
python -m venv .venv
```

On Windows, activate it with:

```powershell
.venv\Scripts\activate
```

On macOS or Linux:

```bash
source .venv/bin/activate
```

Install only the required runtime packages:

```bash
python -m pip install -r requirements.txt
```

Start the application:

```bash
python app.py
```

The terminal prints a local address similar to:

```text
TEMPO Design Explorer is running at:

http://127.0.0.1:4173/
```

The default browser opens automatically. Keep the terminal window open while
using the application. Press `Ctrl+C` in that terminal to stop it.

To start without automatically opening a browser:

```bash
python app.py --no-open
```

To choose a different local port:

```bash
python app.py --port 5000
```

If the default port is already occupied and no explicit port was requested,
TEMPO automatically selects a free local port and prints the new URL.

## Verify an installation

Run the deterministic installation test before using the package for a report or
demonstration. If you used the one-click launcher, run the test with its isolated
environment:

Windows:

```powershell
.venv\Scripts\python.exe test_installation.py
```

macOS or Linux:

```bash
.venv/bin/python test_installation.py
```

If the virtual environment is already activated, the shorter command is:

```bash
python test_installation.py
```

The test checks:

- the Python version;
- all four pinned dependencies;
- the required frontend and scientific source files;
- import of the integrated adapter;
- one complete oscillator → counter → shutdown reference simulation;
- the 8 + 38 + 3 state contract, SciPy LSODA provenance, and output payload.

A successful installation ends with:

```text
TEMPO Design Explorer installation test passed.
```

The reference simulation is intentionally a real ODE solve and can take tens of
seconds on a typical laptop.

## Example workflow

1. Start from **TEMPO reference** or select an oscillator setting K1–K5.
2. Adjust C31 expression and the Integrase degradation tag.
3. Adjust RDF/BM3R1 controls if advanced tuning is required.
4. Adjust sRNA production to change the shutdown behavior.
5. Release the control and wait for the deterministic ODE update to finish.
6. Inspect the synchronized input, counter, and output plots.
7. Read the certification, system metrics, and suggested next adjustment.
8. Export **Design Summary** to preserve the parameters and result metadata.

Representative parameter sets are stored in [`examples/`](examples/README.md):

- `baseline.json`;
- `fast_oscillator.json`;
- `slow_oscillator.json`;
- `strong_shutdown.json`;
- `weak_shutdown.json`.

These examples contain interface-level parameters. They do not replace or
duplicate the scientific model source.

## Model architecture and source provenance

The interactive adapter is located at
`tempo-design-explorer/backend/model_adapter.py`. It maps the visible controls to
the existing model implementations and returns named observables for the UI.

The local application imports these reference sources directly:

| Module | Runtime source |
|---|---|
| Oscillator | `oscillator/code/Mechanistic_ODE_Global_Sensitivity_Analysis.py` |
| Period-knob mapping | `oscillator/code/Shared_PLtetO1_Period_Knob_Design_Map.py` |
| Recombinase counter | `rdfmodel_new/model/zhao_core.py` |
| Coupled counter/shutdown | `shutdown_model/shutdown_core.py` |

The package retains the existing governing equations, parameters, model
connections, `scipy.integrate.solve_ivp`, `method="LSODA"`, tolerances, and
evaluation logic. Packaging changes must not silently change these scientific
settings.

## Repository layout

```text
TEMPO-Design-Explorer/
├── app.py                         # single local entry point
├── requirements.txt               # minimal pinned runtime dependencies
├── run_windows.bat                # Windows environment + launcher
├── run_mac_linux.sh               # macOS/Linux environment + launcher
├── test_installation.py           # clean-install deterministic check
├── examples/                      # reproducible interface parameter sets
├── tempo-design-explorer/         # browser UI, HTTP server, model adapter
├── oscillator/                    # oscillator reference implementation
├── rdfmodel_new/                  # recombinase-counter reference implementation
├── shutdown_model/                # shutdown and coupled-system implementation
├── LICENSE
├── NOTICE
└── THIRD_PARTY_NOTICES.md
```

The broader model directories also contain research scripts and reports retained
for provenance. The runtime imports only the four source files listed above.

## Reproducibility and scientific scope

- Deterministic results are solved locally from the TEMPO Python ODE pipeline.
- The UI displays the model version, solver, state counts, and result source.
- The browser surrogate is only a responsive fallback/preview and must not be
  reported as a deterministic scientific result.
- A deterministic prediction is not experimental validation.
- The heuristic score is manually weighted and provided for design comparison;
  it is not calibrated against experimental outcomes.
- **Upstream oscillation does not imply reliable downstream counting.** K1–K5 are
  oscillator-robust settings, while the current deterministic counter interface
  certifies only K4/K5 under their specified Integrase degradation tags.
- A green **Simulation checks pass** result and the heuristic score are separate
  from A→B certification; always read the certification card and exported metadata.

For a frozen result, retain the exported Design Summary together with the Git
commit, Python version, dependency versions, and model version.

## Troubleshooting

### A dependency is missing

Run from the repository root:

```bash
python -m pip install -r requirements.txt
```

Confirm that the activated Python belongs to `.venv`:

```bash
python -c "import sys; print(sys.executable)"
```

### `venv` cannot be created on Linux

Install the operating system's Python venv package (often named
`python3-venv`), then rerun the launcher.

### The default port is occupied

The normal launcher automatically chooses a free port. Alternatively specify
one explicitly:

```bash
python app.py --port 5000
```

### The page was opened with `file://`

Do not double-click `index.html`. Start `app.py` and use the printed
`http://127.0.0.1:...` URL. Direct file access cannot call the local Python ODE
API reliably.

### A simulation takes time

The complete 49-state workflow includes a stiff 38-state recombinase model.
Parameter changes display a fast browser preview while the deterministic Python
result is being solved. Wait until the status reports **Deterministic ODE** before
recording a result.

### A simulation fails

Reset to **TEMPO reference** and retry. Public controls are bounded to the model's
supported ranges. If failure persists, run `python test_installation.py` and
include its output when reporting the problem.

## Online demo and Wiki

The Render deployment is a convenient online preview. The reproducible package
in the iGEM GitLab repository is the long-term scientific deliverable. The Wiki
should explain the architecture and worked examples with screenshots or video,
then point users to this package for complete local simulation.

## Team and citation

Developed by **USTC iGEM 2026** for the **TEMPO** project.

When reusing the toolkit, cite the USTC iGEM 2026 TEMPO Wiki Model page and
record the repository commit used for the simulation. Final Wiki and attribution
URLs should be added before the competition freeze.

## License

Team-authored code and documentation are released under the
[Apache License 2.0](LICENSE). Third-party dependencies retain their own
licenses; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). Team member roles
should be completed in [AUTHORS.md](AUTHORS.md) and on the official iGEM
attribution page.
