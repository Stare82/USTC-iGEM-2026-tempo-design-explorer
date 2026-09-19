# TEMPO Modeling Platform / Design Explorer

TEMPO Design Explorer connects the project’s oscillator, recombinase counter, and timed-shutdown models in one interactive interface. The public application serves both the static frontend and a deterministic Python/SciPy ODE API from the same origin.

## Repository layout

- `tempo-design-explorer/`: frontend, API adapter, documentation, and tests.
- `oscillator/code/`: certified oscillator and period-knob runtime sources.
- `rdfmodel_new/model/zhao_core.py`: 38-state recombinase counter source.
- `shutdown_model/shutdown_core.py`: coupled shutdown source.
- `Dockerfile`: minimal production image containing only the runtime sources above.
- `render.yaml`: one-service Render Blueprint.
- `.github/workflows/model-smoke.yml`: deterministic smoke test for GitHub.
- `.gitlab-ci.yml`: equivalent pipeline retained for a future iGEM GitLab mirror.

Generated model outputs are intentionally excluded from version control and the deployment image. They can be reproduced from the documented source models.

## Local development

```powershell
python -m pip install -r tempo-design-explorer/requirements.txt
python tempo-design-explorer/server.py --open
```

Open `http://127.0.0.1:4173/`. Do not open `index.html` through `file://`; that bypasses the Python ODE service.

## Verification

```powershell
python -m py_compile tempo-design-explorer/server.py tempo-design-explorer/backend/model_adapter.py
python tempo-design-explorer/tests/test_deployment_smoke.py
```

The smoke test solves the public K4 reference design and confirms the 8 + 38 + 3 state deterministic response contract.

## Cloud deployment

The repository is prepared as a single Docker web service. See [DEPLOYMENT.md](DEPLOYMENT.md) for the exact Render and iGEM GitLab workflow.

## Scientific boundary

The deterministic ODE results are model predictions, not experimental validation. The interface clearly labels its browser fallback and its manually weighted heuristic score. Pass/fail certification boundaries are kept separate from that reference-only score.

## License

Team-authored source code and repository documentation are released under the
[Apache License 2.0](LICENSE), with copyright attributed to the USTC iGEM Team
contributors. The project name and logo may not be used to imply endorsement.
Third-party dependencies retain their own licenses; see
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

Downloaded literature used during model research is not redistributed. Team
members should record individual roles in [AUTHORS.md](AUTHORS.md) and on the
official iGEM attribution page before the competition freeze.
