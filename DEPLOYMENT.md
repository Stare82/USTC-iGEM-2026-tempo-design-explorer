# Deployment guide

## Recommended publication model

1. Use the dedicated iGEM GitLab software repository as the canonical source.
2. Deploy this repository as one Docker web service so the frontend and `/api/*` remain same-origin.
3. Link the resulting HTTPS URL from the team wiki’s Model and Software pages.
4. Keep GitHub only as an optional mirror.

GitHub Pages and GitLab Pages are static hosts and cannot run this project’s Python/SciPy backend.

## Files included in the runtime image

The Docker build deliberately copies only:

- `tempo-design-explorer/`;
- `oscillator/code/Mechanistic_ODE_Global_Sensitivity_Analysis.py`;
- `oscillator/code/Shared_PLtetO1_Period_Knob_Design_Map.py`;
- `shutdown_model/shutdown_core.py`;
- `rdfmodel_new/model/zhao_core.py`.

Do not change these paths without also updating `backend/model_adapter.py` or the `Dockerfile`.

## Local container check

After installing Docker Desktop:

```powershell
docker build -t tempo-design-explorer .
docker run --rm -p 10000:10000 -e PORT=10000 tempo-design-explorer
```

Then verify:

- `http://127.0.0.1:10000/`
- `http://127.0.0.1:10000/api/health`

## Render deployment

The root `render.yaml` defines a single Docker service with `/api/health` as its health check.

1. Push the prepared repository to
   `https://github.com/Stare82/USTC-iGEM-2026-tempo-design-explorer`.
2. In Render, create a Blueprint or Web Service from that GitHub repository.
3. Keep the official USTC iGEM GitLab project for the team wiki. Link the live
   platform and this source repository from the wiki's Model page.
4. Keep the free plan for testing only. Upgrade before judging or public demonstrations to avoid cold starts.
5. After deployment, open the site and verify the upper-right badge reads `Deterministic ODE · 49 states`.
6. Change at least one parameter and confirm that `/api/simulate` returns a new deterministic result.

`autoDeployTrigger` is intentionally `off` in `render.yaml`. This prevents an unreviewed commit from replacing the public judging build. Trigger a manual deployment after CI passes.

## iGEM release checklist

- [x] Confirm the public source repository URL.
- [x] Add the OSI-approved Apache-2.0 license and dependency notices.
- [ ] Add the official team name, authors, and acknowledgements.
- [ ] Push the source and confirm the GitHub Actions smoke-test workflow passes.
- [ ] Deploy the exact passing commit to the cloud service.
- [ ] Record the deployed commit hash and model version on the wiki.
- [ ] Test desktop and mobile layouts from an external network.
- [ ] Test first-load cold-start behavior before judging.
- [ ] Add the live URL and source URL to the wiki’s Model and Software pages.
- [ ] Preserve the deterministic/fallback and heuristic-score disclaimers.
