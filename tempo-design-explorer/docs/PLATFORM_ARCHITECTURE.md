# TEMPO Modeling Platform：信息架构与 MVP 技术方案

## 1. 产品定位

平台的首要任务不是“展示曲线”，而是帮助合成生物学设计者完成一个判断：

> 给定一组可设计的生物元件，这套 oscillator–counter–shutdown 时间控制链是否仍处于可行范围？

因此产品采用 **working surface**，页面打开后直接看到参数、动态轨迹与设计评价。Wiki 的机制说明、推导与文献证据仍保留在叙事页面；Explorer 负责参数试验、因果反馈与设计决策。

## 2. 信息架构

### 2.1 单页面 MVP

```text
TEMPO Design Explorer
├── Design input
│   ├── Basic mode
│   │   ├── Certified period setting K1–K5 / C31 pulse-width target
│   │   ├── C31 RBS / Integrase degradation / RDF expression
│   │   └── BM3R1 expression / sRNA production
│   ├── Advanced mode
│   │   └── Oscillator promoter / RBS scales（第一版）
│   └── Presets
├── Dynamic response
│   ├── Input & protein dynamics
│   ├── Recombinase counter state
│   └── Timed output & shutdown
├── Design evaluation
│   ├── Counter quality
│   ├── Output timing
│   └── Explainable diagnostics
└── Export
    └── Parameters + mapping + evaluation JSON
```

### 2.2 后续扩展，不纳入当前 MVP

- Parameter scan：二维 heatmap、feasible region 与 failure boundary；
- Cascade：Bit 1 → carry → Bit 2 / Bit 3；
- Design recommendation：从目标脉冲数与 shutdown 时间反推参数范围；
- 模型版本对比与实验数据 overlay；
- 保存、分享和团队协作。

这些功能应在真实 ODE 计算接口稳定后再进入，避免第一版出现多个不一致的计算入口。

## 3. 页面布局

### 3.1 桌面端

三栏工作台：

| 区域 | 宽度建议 | 作用 |
|---|---:|---|
| 左侧 Parameter Control | 280–320 px | 参数分组、Basic / Advanced、presets |
| 中间 Dynamic Visualization | 弹性主区域 | 三张同步时间轴、pulse marker、PB/LR state rail |
| 右侧 Design Evaluation | 300–340 px | score、判据、时间指标、失败原因、导出 |

顶部只保留 TEMPO 标识、三阶段因果链和模型状态。没有 marketing hero，用户进入页面后可立即操作。

### 3.2 响应式

- 1180 px 以下：evaluation 移到曲线下方并改为三列摘要；
- 820 px 以下：参数区、曲线区、评价区纵向排列；
- 560 px 以下：图例、标题、指标卡改为单列，避免横向滚动；
- 触摸设备上的 slider 保留足够命中区域，主要文字不小于 14 px。

## 4. MVP 功能边界

### 已实现

1. 七个 Basic 参数实时调节；
2. 两个 Advanced 元件参数；
3. 五组预设；
4. 三组同步动态轨迹；
5. one-pulse-one-flip、flip success、counter recovery、leakage；
6. output start、duration、shutdown time、residual expression；
7. 解释性诊断；
8. JSON 设计摘要导出；
9. browser surrogate / deterministic ODE / deterministic certification 边界提示；
10. K1–K5 稳健周期库与 A→B tag-specific 认证提示；
11. `stage_tempo_design` 与 `read_tempo_result` WebMCP 工具。

### 明确不做

- 不把确定性 ODE 结果称为随机模型或实验验证；
- 不将 surrogate 结果用于实验参数定稿；
- 不在第一版加入账号、数据库、多人协作或批量扫描；
- 不把未来功能做成无效按钮。

## 5. Design parameter → Model parameter

前端只暴露设计者能理解的术语，计算层负责转换：

| 设计参数 | 真实模型参数 | 当前参考点 |
|---|---|---:|
| Certified period setting K1–K5 | `oscillator_translation_scale` + `oscillator_mrna_total_half_life_min` | K4: 1.20× + 1.5 min |
| Oscillator promoter strength | `oscillator_tx_per_plasmid_per_min` | 0.9609 min⁻¹ at 1× |
| Oscillator RBS strength | `oscillator_translation_per_mrna_per_min` | 0.5 min⁻¹ at 1× |
| C31 RBS strength | `c31_rbs_scale` | 0.45× |
| Integrase degradation tag | `k_tag_int_h` | 8 h⁻¹ |
| RDF expression | `krdf_tsl_h` | 200 h⁻¹ at 1× |
| BM3R1 expression | `krep_tsl_h` | 15 h⁻¹ at 1× |
| sRNA promoter strength | `srna_max_tx_uM_h` | 0.06 μM/h |

K1–K5 来自 `ranked_five_setting_libraries.csv`。K4 在 8/12 h⁻¹ tag 下通过当前确定性 A→B 认证，K5 仅在 12 h⁻¹ 下通过；K1–K3 仅是稳健振荡器档位，不是已通过 B 接口的档位。定性选项 `Weak / Medium / Strong` 后续仍应映射到团队实际选定的 part ID、测量条件与不确定度。

## 6. 前端与模型计算的连接

### 6.1 稳定数据契约

前端只依赖如下计算请求：

```json
{
  "model_version": "tempo-certified-period-knob-ode-2026.09",
  "simulation_hours": 60,
  "initial_counter_state": "PB",
  "design": {
    "period_setting": 4,
    "oscillator_promoter_scale": 1.0,
    "oscillator_rbs_scale": 1.0,
    "c31_rbs_scale": 0.45,
    "integrase_tag_rate_h": 8.0,
    "rdf_expression_scale": 1.0,
    "bm3r1_expression_scale": 1.0,
    "srna_max_tx_uM_h": 0.06
  }
}
```

计算层返回：

```json
{
  "time_h": [0.0, 0.02],
  "series": {
    "c31_flux_uM_h": [],
    "integrase_total_uM": [],
    "rdf_total_uM": [],
    "pb_fraction": [],
    "lr_fraction": [],
    "bm3r1_total_uM": [],
    "shutdown_srna_uM": [],
    "output_protein_uM": []
  },
  "evaluation": {
    "fraction_one_and_complete": 1.0,
    "alternation_fidelity": 1.0,
    "recovery_before_next_pulse": true,
    "output_start_h": 0.0,
    "output_duration_h": 0.0,
    "shutdown_h": 0.0,
    "residual_fraction": 0.0,
    "warnings": []
  },
  "provenance": {
    "solver": "LSODA",
    "source_commit": "...",
    "parameter_set": "..."
  }
}
```

### 6.2 推荐实现路径

**阶段 A：浏览器 fallback。** `src/model/tempo-model.js` 即时计算 surrogate，保证模型服务不可用时界面仍可操作，并明确标记结果来源。

**阶段 B：预计算网格。** 用现有 Python 模型在可信参数范围内批量运行，导出规则网格或 Latin hypercube 样本。浏览器用插值即时反馈，同时清楚显示“interpolated”。该路径最适合 iGEM Wiki 的纯静态限制。

**阶段 C：真实求解服务（当前已实现本地版）。** `server.py` 调用 `oscillator → zhao_core → shutdown_core`，返回统一 contract。前端对 slider 输入做 debounce，并使用 request revision 丢弃过期响应。服务端缓存规范化参数哈希，限制参数范围和计算时长。

如果 Wiki 不允许外部计算服务，优先使用阶段 B，而不是把 SciPy 41-state 求解器匆忙翻译到 JavaScript。

### 6.3 模型适配层

现有代码可以按三层封装：

1. `OscillatorAdapter`：返回真实 `translation_flux_uM_h(t)`，不把 C31 concentration 重复积分；
2. `CounterAdapter`：以 C31 production flux 驱动 `zhao_core.py`，返回 PB / LR、Integrase、RDF、BM3R1；
3. `ShutdownAdapter`：读取 LR fraction 与 free BM3R1，返回 sRNA、output mRNA / protein 与 pulse metrics。

前端不得依赖 ODE state index（例如 y[38]），只依赖命名后的 observable。

## 7. 仍需团队提供或确认的数据

### 高优先级：决定模型是否可用于设计

1. 最终认定的 oscillator 代码版本与默认参数集；
2. promoter / RBS 实物 part 与 scale 的映射表，而不只是 0.5× / 1×；
3. 最终 `T_d` 条件：宿主、培养基、温度、抗生素、质粒负载；
4. PLtetO1 leak、Hill coefficient、binding threshold 的采用值及不确定度；
5. C31 active fraction、细胞体积和 copies ↔ μM 换算口径；
6. C31/RDF/BM3R1 的实测降解或半衰期；
7. counter 的初态定义、目标 state、完整 flip 阈值；
8. one-pulse-one-flip、recovery、leakage 的最终 pass / fail 阈值。

### 中优先级：决定 shutdown 指标能否定量解释

9. BM3R1 promoter 的 transfer function；
10. RyhB/sodB proxy 是否为最终构建；若不是，需新 sRNA–target 配对速率；
11. 输出蛋白的翻译率、成熟时间与主动降解；
12. “shutdown time”的实验定义：production rate、mRNA 还是 protein 低于哪个阈值；
13. 单 bit 与 cascade 中 output activation 的目标逻辑；
14. 至少一组实验时间轨迹，用于校准而非仅验证曲线形状。

### 展示与复现元数据

15. 每个参数的 evidence level、来源、单位、有效范围；
16. 模型 commit / version、求解器、容差、采样间隔；
17. 参数扫描的完整结果表，而非只提供图片；
18. Wiki 页面计划使用的语言、品牌色与嵌入尺寸。

## 8. iGEM Wiki 部署方案

### 方案 1：纯静态 bundle（MVP 推荐）

- HTML、CSS、原生 ES modules，无外部 CDN；
- 构建后使用相对路径，放入 Wiki 允许的静态资源目录；
- 通过 Wiki 页面中的容器或 iframe 嵌入；
- 模型采用 surrogate 或预计算网格；
- 优点：最稳、最容易审核、无服务端依赖；
- 限制：无法对任意参数运行完整 SciPy ODE。

当前原型即按此方案实现。部署前应根据当年 iGEM Wiki 的 CSP、脚本与 iframe 规则做一次实际验证。

### 方案 2：静态前端 + 外部计算 API

- 前端仍由 Wiki 托管；
- API 单独部署，Python 环境复用现有 SciPy 代码；
- 必须配置 CORS、HTTPS、限流、缓存、超时与版本 provenance；
- 前端在 API 不可用时保留 reference preset 或预计算 fallback；
- 优点：真实 ODE、任意参数；
- 风险：竞赛展示时依赖外部服务可用性。

### 方案 3：离线 WebAssembly

只有在确定必须完全离线、且阶段 B 的插值精度不足时再评估。把 SciPy/LSODA 直接搬到浏览器成本较高，包体和兼容性也不适合第一版。

## 9. 验收标准

MVP 可以进入 Wiki 联调的最低标准：

- 默认参数打开后立即显示完整因果链；
- 任一参数变化后 300 ms 内先刷新明确标注的 surrogate；确定性 ODE 在后台完成后替换；
- PB + LR 在任意时间点近似守恒为 1；
- 状态评价与曲线不会明显矛盾；
- 所有单位、默认值和模型版本可追溯；
- 键盘可操作 slider、select、toggle 与导出；
- 320 px 到 1440 px 无横向溢出；
- browser surrogate、deterministic ODE、deterministic A→B certification 三种状态不会混淆；
- 真实 ODE 接入后，默认点与现有 Python 输出的关键指标在预先约定容差内一致。

## 10. 与现有仓库的对应关系

| 平台层 | 当前来源 |
|---|---|
| Oscillator | `oscillator/code/Mechanistic_ODE_Global_Sensitivity_Analysis.py` |
| Period knob | `Shared_PLtetO1_Period_Knob_Design_Map.py`、`ranked_five_setting_libraries.csv` |
| A → B interface | `Robust_Continuous_Period_Knob_AB_Certification.py`、`AB_tag_specific_certification.csv` |
| Counter | `rdfmodel_new/model/zhao_core.py` 与扫描脚本 |
| Shutdown | `shutdown_model/shutdown_core.py` |
| 当前默认与窗口 | 三部分 README、CSV summary 与 parameter provenance |

平台没有修改三套科学模型；`backend/model_adapter.py` 只负责参数映射、求解编排与稳定 JSON 契约。
