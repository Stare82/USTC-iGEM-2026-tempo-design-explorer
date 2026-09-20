# TEMPO Design Explorer

> A reproducible local modeling interface that connects our oscillator, recombinase counter, and timed-shutdown models into one design workflow.

<!--
页面定位：Modeling contribution 的 interactive integration layer，不申报 Best Software。
正式长期交付：iGEM GitLab 中的可下载本地 Python 包。
Wiki：负责解释、截图、视频、worked example 和安装指引，不承担实时 SciPy 计算。
不要把 Render 或 GitHub 写成正式交付入口；如保留在线演示，也只作为 optional preview。
所有图片上传至 static.igem.wiki，视频使用当年 iGEM 允许的官方视频托管方式。
模型机制、ODE 推导、参数来源和 sensitivity analysis 仍放在 Model 页面。
-->

![TEMPO Design Explorer interface](INTERNAL_PATH/software-hero-overview.webp)

*Figure 1. TEMPO Design Explorer places biologically meaningful controls, three synchronized model responses, deterministic certification, and explainable design evaluation in one interface. Each chart is normalized independently, so timing and curve shape—not absolute height between charts—should be compared.*

<!--
图 1：16:9 软件总览图，不要使用很长的 full-page screenshot。
保留 K4 参数区、三张同步曲线、certification card 和 evaluation inset。
截图时必须显示 “Deterministic ODE · K4 updated”，不能处于 loading、waiting 或 fallback 状态。
建议同时露出新增的 Interface boundary 提示。
-->

## Why We Built TEMPO

Our timing system contains three computational modules: an oscillator, a recombinase counter, and a timed-shutdown circuit. Each model can be studied independently, but a parameter change in one module also changes the inputs and behavior of the downstream modules. Separate scripts and static plots therefore made system-level design difficult, especially for team members without modeling experience.

We developed **TEMPO Design Explorer** to connect the complete computational workflow in one interface. Users can select an oscillator setting, adjust biologically meaningful parameters, run the coupled Python models, inspect the resulting dynamics, and read why a candidate design passes a simulation check or needs further adjustment.

> **Central modeling lesson:** robust upstream oscillation does not automatically produce reliable downstream counting.

TEMPO is not a new biological model. It is the interactive integration layer for our existing oscillator, counter, and shutdown models.

## What TEMPO Does

TEMPO supports five tasks:

1. **Design input** — select a K1–K5 oscillator setting and adjust C31, Integrase-reset, RDF/BM3R1, and shutdown parameters.
2. **Coupled simulation** — run the existing 8-state oscillator, 38-state counter, and 3-state shutdown models through one parameter interface.
3. **Dynamic visualization** — compare input, counter-state, and timed-output responses on synchronized time axes.
4. **Design evaluation** — report switching fidelity, recovery, output duration, residual output, deterministic certification, and diagnostic suggestions.
5. **Reproducible reuse** — download the complete source package, create an isolated Python environment, and reproduce the simulation locally.

The biological mechanisms, governing equations, parameter sources, sensitivity analyses, and certification experiments are documented on our [Model page](INTERNAL_MODEL_PAGE).

## Model-to-Interface Architecture

```text
K1–K5 oscillator design
        ↓ C31(t)
8-state oscillator model
        ↓ Integrase input
38-state recombinase counter
        ↓ PB/LR DNA state
3-state timed-shutdown model
        ↓
system output + diagnostics
```

The interface does not replace the scientific models with JavaScript. Reportable results are calculated locally by the existing Python/SciPy implementation using `solve_ivp` with `LSODA`. A browser approximation is displayed only as a clearly labelled transient fallback and is not a scientific result.

## How to Use TEMPO

<!--
视频不链接 Render。将录屏上传到 iGEM 官方允许的媒体托管位置后嵌入。
推荐视频：35–50 秒、无配音、带英文字幕、1440×900、Basic mode、60 h。
视频流程：启动本地程序 → K4 baseline → 查看 warning → sRNA 0.09 → 等待 ODE → checks pass。
-->

<video controls muted playsinline poster="INTERNAL_PATH/tempo-walkthrough-poster.webp">
  <source src="INTERNAL_VIDEO_PATH/tempo-walkthrough.webm" type="video/webm">
  <source src="INTERNAL_VIDEO_PATH/tempo-walkthrough.mp4" type="video/mp4">
  Your browser does not support the video element.
</video>

*Video 1. A short walkthrough from launching the local package to selecting the K4 reference design, reading the diagnostics, and strengthening the shutdown branch.*

### In-App Workflow

1. Select a named preset or a K1–K5 oscillator setting.
2. Check whether the selected setting is inside the current deterministic counter-certification region.
3. Adjust C31 expression and Integrase degradation.
4. Adjust RDF/BM3R1 parameters if advanced tuning is required.
5. Adjust sRNA production to shape output closure.
6. Wait until the status reads **Deterministic ODE**.
7. Read simulation checks, certification, key outcomes, and suggested adjustments separately.
8. Export the Design Summary to preserve parameters, metrics, and model metadata.

### Screenshot Alternative

<!-- 四张图可以做成一张横向 step strip，避免重复占据页面。 -->

| Step | Screenshot | What the reader should notice |
|---|---|---|
| 1. Launch | `INTERNAL_PATH/quickstart-01-launch.webp` | Localhost URL and local Python server |
| 2. Select | `INTERNAL_PATH/quickstart-02-select.webp` | K4, biological controls, and interface-boundary note |
| 3. Evaluate | `INTERNAL_PATH/quickstart-03-result.webp` | Simulation checks and certification shown separately |
| 4. Adjust | `INTERNAL_PATH/quickstart-04-adjust.webp` | Stronger shutdown and the updated deterministic result |

## Run the Full Explorer Locally

The Wiki explains the platform and shows worked results. The complete interactive simulator is distributed as a reproducible local Python package in the official USTC iGEM GitLab repository.

- [Browse the TEMPO source package on iGEM GitLab](https://gitlab.igem.org/2026/ustc/-/tree/main/model/tempo-design-explorer)
- [Download the official USTC repository ZIP](https://gitlab.igem.org/2026/ustc/-/archive/main/ustc-main.zip)

After extracting the ZIP, open:

```text
<extracted USTC repository>/
└── model/
    └── tempo-design-explorer/
```

### Requirements

- Python 3.11 or newer;
- Python 3.13.9 recommended and tested;
- Windows, macOS, or Linux;
- an internet connection during the first dependency installation;
- a modern browser.

### Windows

Double-click:

```text
run_windows.bat
```

The launcher checks Python, creates `.venv`, installs the pinned dependencies, starts the local server, and opens the browser automatically.

### macOS or Linux

From a terminal in the package directory:

```bash
sh run_mac_linux.sh
```

### Manual Installation

```bash
python -m venv .venv
```

Windows:

```powershell
.venv\Scripts\activate
```

macOS/Linux:

```bash
source .venv/bin/activate
```

Then install and run:

```bash
python -m pip install -r requirements.txt
python app.py
```

The terminal displays a local address such as:

```text
TEMPO Design Explorer is running at:

http://127.0.0.1:4173/
```

Keep the terminal open while using the platform. Press `Ctrl+C` to stop it.

### Verify the Installation

Windows:

```powershell
.venv\Scripts\python.exe test_installation.py
```

macOS/Linux:

```bash
.venv/bin/python test_installation.py
```

A successful check ends with:

```text
TEMPO Design Explorer installation test passed.
```

The test verifies the pinned dependencies, required source and reference files, the complete oscillator → counter → shutdown solve, the 8 + 38 + 3 state contract, and SciPy LSODA provenance.

### Local Environment and Privacy

The launcher creates `.venv` inside the extracted TEMPO folder. It does not replace the user's system Python packages, modify other projects, require administrator privileges, or expose the application publicly. The server binds to `127.0.0.1`, so the model and parameters remain on the user's computer. Removing the extracted folder also removes its isolated environment.

## Worked Example: From Diagnosis to Design Adjustment

<!--
以下数字已经用当前冻结代码重新求解：
K4, tag 8, sRNA 0.060, 60 h → score 97, warning, residual 22.36%, 7/7 flips。
K4, tag 8, sRNA 0.090, 60 h → score 98, success, residual 13.84%, 7/7 flips。
截图仍应在最终 release tag 后重新生成，避免 UI 与数字版本不一致。
-->

### Step 1 — Start from the K4 Reference Design

| Design parameter | Reference value |
|---|---:|
| Oscillator period setting | K4 |
| Nominal / realized period | 9.35 / 9.44 h |
| C31 RBS strength | 0.45× |
| Integrase degradation | 8 h⁻¹ |
| sRNA production | 0.060 μM h⁻¹ |
| Simulation window | 60 h |

### Step 2 — Read the Result

The reference design produces seven detected input pulses and seven complete counter flips in the displayed 60-hour window. The interface reports 100% flip fidelity, recovery before the next pulse, and deterministic A→B certification for K4 with the selected `8 h⁻¹` Integrase tag.

However, the overall simulation still needs attention. Residual output is `22.36%`, slightly above the `20%` shutdown target. The reference-only heuristic score is `97/100`, demonstrating why a high score or an upstream interface certification must not be interpreted as proof that every downstream requirement is satisfied.

![K4 diagnostic before adjustment](INTERNAL_PATH/k4-diagnosis-before.webp)

*Figure 2. The K4 counter interface is certified at the selected tag, but residual output remains above the shutdown target. Certification, heuristic score, and required system checks are intentionally reported separately.*

### Step 3 — Follow the Suggested Adjustment

TEMPO recommends strengthening sRNA production. Increasing the setting from `0.060` to `0.090 μM h⁻¹` preserves seven complete flips and counter recovery while reducing residual output to `13.84%`. The heuristic score becomes `98/100`, and the interface reports **Simulation checks pass**.

![K4 result after adjustment](INTERNAL_PATH/k4-diagnosis-after.webp)

*Figure 3. Stronger sRNA production closes the output more completely while retaining the certified K4 counter interface.*

### Design Lesson

```text
Run the reference design
        ↓
Identify the failed requirement
        ↓
Connect it to a controllable biological parameter
        ↓
Rerun the coupled Python model
        ↓
Compare the new engineering outcome
```

## Upstream Oscillation Does Not Guarantee Downstream Counting

The K1–K5 settings form a robust upstream oscillator library. They do not all belong to the certified oscillator-to-counter operating region. The current counter has an approximately `8.75 h` deterministic recovery boundary under the tested `12 h⁻¹` Integrase tag condition.

| Setting | Nominal period | Oscillator status | Current deterministic A→B certification |
|---|---:|---|---|
| K1 | 6.05 h | Oscillator robust | Not counter-certified |
| K2 | 6.93 h | Oscillator robust | Not counter-certified |
| K3 | 8.14 h | Oscillator robust | Not counter-certified |
| K4 | 9.35 h | Oscillator robust | Certified at tags 8 and 12 h⁻¹ |
| K5 | 10.45 h | Oscillator robust | Certified at tag 12 h⁻¹ |

K1–K3 remain useful for exploring the upstream clock, but their shorter pulse spacing lies below or near the present downstream recovery boundary. K3 can satisfy the displayed simulation checks under some advanced shutdown settings, but that does not retroactively certify its A→B interface. Accordingly, TEMPO uses two separate statements:

- **Simulation checks pass** — the current deterministic run meets the displayed flip, recovery, and residual-output checks;
- **A→B certified** — the selected K setting, tag, and frozen interface conditions belong to the previously tested deterministic certification region.

These statements must not be merged into a single claim.

![Deterministic counter recovery boundary](INTERNAL_PATH/counter-recovery-boundary.webp)

*Figure 4. The upstream oscillator library spans a broader period range than the current downstream counter can reliably accept.*

## From Interface Controls to Engineering Decisions

![Mapping interface controls to design decisions](INTERNAL_PATH/controls-to-decisions.svg)

*Figure 5. TEMPO exposes design-level controls rather than internal ODE state indices.*

| Interface control | Question answered by TEMPO |
|---|---|
| K1–K5 period setting | Is the oscillator setting compatible with the current downstream recovery time? |
| C31 RBS strength | Is the Integrase pulse sufficient for one complete DNA flip? |
| Integrase degradation | Can the counter recover before the next pulse? |
| RDF/BM3R1 expression | Does the counter switch and reset with sufficient timing margin? |
| sRNA production | Does output close within the intended time window? |
| Simulation duration | Is the selected observation window long enough to reveal repeated behavior? |

## How TEMPO Informed Our Project

### Selecting the Oscillator–Counter Interface

- **Question:** Which robust oscillator settings also satisfy the current downstream counter-interface criteria?
- **TEMPO result:** K1–K3 are robust oscillator settings but fall outside the current deterministic counter certification. K4 passes with tags `8/12 h⁻¹`, while K5 requires `12 h⁻¹`.
- **Modeling recommendation:** use K4 as the reference coupled design and retain K5 as a slower alternative; do not present K1–K3 as certified end-to-end settings.
- **Team decision:** `[TEAM TO CONFIRM: insert the construct or condition prioritized for experiments.]`
- **Experimental follow-up:** measure pulse spacing, flip completion, and inter-pulse recovery in the selected construct.

### Adjusting the Shutdown Window

- **Question:** Is residual expression below the target after the timed output window?
- **TEMPO result:** the K4 reference predicts `22.36%` residual output at `0.060 μM h⁻¹` sRNA production; `0.090 μM h⁻¹` reduces it to `13.84%` while preserving the counter checks.
- **Modeling recommendation:** prioritize the stronger sRNA setting as a candidate for experimental comparison.
- **Team decision:** `[TEAM TO CONFIRM: insert the sRNA construct selected for testing.]`
- **Experimental follow-up:** collect time-series output and sRNA measurements to test the predicted closure behavior.

## Scientific Scope

TEMPO deliberately separates four levels of evidence:

1. **Browser approximation** — responsive preview only; never used as a reportable model result.
2. **Deterministic ODE result** — calculated locally from the existing coupled Python models.
3. **Deterministic A→B certification** — displayed only for frozen combinations inside the tested counter-interface boundary.
4. **Stochastic or experimental validation** — not implied by a deterministic result or heuristic score.

The heuristic score is manually weighted for design comparison and is not calibrated against experimental outcomes. A green **Simulation checks pass** result is therefore not equivalent to A→B certification, stochastic robustness, or experimental validation.

For biological mechanisms, mathematical derivations, source parameters, and validation work, see the [Model page](INTERNAL_MODEL_PAGE).

## Implementation and Reproducibility

```text
Browser controls on localhost
        ↓ JSON request
Python model adapter
        ↓
Existing deterministic source models
        ↓
Named time series + metadata + evaluation metrics
        ↓
Charts, certification card, diagnostics, and export
```

The reproducible package includes:

- a single `app.py` entry point;
- isolated Windows and macOS/Linux launchers;
- pinned versions of NumPy, SciPy, pandas, and Matplotlib;
- the unchanged reference model sources;
- baseline and alternative example parameter sets;
- `test_installation.py` for a complete deterministic installation check;
- an exportable Design Summary linking parameters to outcomes and provenance.

The scientific reference implementation retains the existing equations, coupling, parameter sets, `solve_ivp(method="LSODA")`, and solver tolerances.

## Limitations

- TEMPO currently presents deterministic ODE predictions, not experimental validation.
- K1–K3 are upstream oscillator settings but are not certified for the current downstream counter interface.
- K4/K5 deterministic certification does not imply stochastic success; the current strict stochastic paired-cycle success rates are `50.0%` and `60.7%`, respectively, and must be reported separately on the Model page.
- The browser approximation must not be used for scientific parameter selection or Wiki result figures.
- The heuristic score is reference-only and should never replace the individual checks or certification card.
- Several design scales still require mapping to final physical part IDs and uncertainty ranges.
- Experimental time-series data are required for quantitative validation.

<!--
发布前检查清单：
1. 将所有 INTERNAL_PATH / INTERNAL_VIDEO_PATH / INTERNAL_MODEL_PAGE 替换为正式站内路径。
2. 用最终 GitLab main 或 release tag 重新运行 K4 baseline 和 strong_shutdown。
3. 截图必须来自 Deterministic ODE 状态，不能使用 fallback。
4. 核对数值：baseline residual 22.36%，strong shutdown residual 13.84%。
5. 补充并确认两个 TEAM TO CONFIRM 项，不能把 modeling recommendation 写成已经完成的实验决定。
6. 确认 GitLab 下载链接、Windows launcher、macOS/Linux launcher 和 test_installation.py。
7. Wiki 只链接官方 iGEM GitLab，不把 GitHub 或 Render 作为正式交付入口。
-->
