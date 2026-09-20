# TEMPO Design Explorer — local application

这是 TEMPO Modeling Platform 的本地交互应用。它把现有的三部分模型组织为一条连续设计链：

`Oscillator → C31 / Integrase pulse → Recombinase counter → Timed output / Shutdown`

## 当前版本包含

- Basic / Advanced 两层参数控制；
- 五组设计预设与一键恢复默认值；
- Input / protein、PB/LR counter state、output / sRNA 三组同步曲线；
- one-pulse-one-flip、恢复时间、leakage、输出时间窗等自动评价；
- 失败原因解释，而不只显示状态颜色；
- JSON 设计摘要导出；
- 桌面、平板和手机响应式布局；
- 前端无框架；本地科学计算使用仓库现有 Python/SciPy/Pandas 依赖。

## 模型模式与科学边界

平台现在有两个明确区分的计算模式：

- **Deterministic ODE**：通过本地服务调用仓库中现有的 8-state oscillator、38-state counter 和 3-state shutdown 模型，使用 SciPy LSODA 求解；
- **Browser surrogate · fallback**：模型服务不可用时才显示，用于保持界面可操作，不作为科研预测。

页面右上角和导出的 JSON 都会写明当前结果来源，不会把 fallback 曲线冒充为 ODE 输出。周期旋钮现在直接使用 `oscillator` 报告中的五档稳健构建设计，而不是按目标周期拉伸曲线：

- K1: RBS `0.60×`, oscillator mRNA half-life `1.0 min`, nominal `6.05 h`；
- K2: RBS `0.80×`, half-life `1.0 min`, nominal `6.93 h`；
- K3: RBS `1.20×`, half-life `1.0 min`, nominal `8.14 h`；
- K4: RBS `1.20×`, half-life `1.5 min`, nominal `9.35 h`；
- K5: RBS `0.60×`, half-life `4.0 min`, nominal `10.45 h`；
- C31 flux peak: `6.639 μM/h`；
- C31 RBS candidate: `0.45×`；
- Integrase tag candidate: `8 h⁻¹`；
- sRNA production reference: `0.06 μM/h`。

确定性 A→B 认证边界也会显示在诊断中：K4 在 `8/12 h⁻¹` tag 下通过，K5 仅在 `12 h⁻¹` 下通过，K1–K3 未通过当前 B 恢复认证。改变冻结接口参数后，平台仍会求解真实 ODE，但不会继续显示“已认证”。随机模型认证与实验验证不能由确定性 ODE 替代。

## 普通浏览器启动

推荐从仓库根目录双击：

```text
run_windows.bat
```

本目录中的 `start_platform.bat` 也会转交给同一个根目录启动器。启动器会创建独立 `.venv`、安装根目录 `requirements.txt` 中锁定的依赖，并打开本地页面。不要直接双击 `index.html`：普通浏览器会限制 `file://` 页面加载 ES modules，而且这种方式无法调用 Python ODE。

也可以在仓库根目录手动运行：

```powershell
python app.py
```

首次完整 ODE 求解在普通笔记本上可能需要 1–2 分钟，具体取决于处理器；同一组参数会由服务端缓存。

## 目录

```text
tempo-design-explorer/
├── index.html
├── styles.css
├── server.py
├── start_platform.bat
├── backend/
│   └── model_adapter.py
├── public/favicon.svg
├── src/
│   ├── app.js
│   ├── charts.js
│   └── model/
│       ├── parameters.js
│       └── tempo-model.js
└── docs/
    └── PLATFORM_ARCHITECTURE.md
```

完整的信息架构、真实 ODE 接入约定、待补数据与部署方案见 [`docs/PLATFORM_ARCHITECTURE.md`](docs/PLATFORM_ARCHITECTURE.md)。
