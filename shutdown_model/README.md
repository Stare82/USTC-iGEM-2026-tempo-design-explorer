# BM3R1-driven sRNA shutdown model

本目录连接振荡器、`rdfmodel_new`计数器与BM3R1驱动的sRNA Shutdown模块。版本3.3采用RyhB/sodB文献代理，补充RyhB序列审查、参数底层来源和单/多bit解释，并为参数包络和生长分析使用逐条件无Shutdown对照。

## 方程

```text
dS/dt = alpha_S * H(BM3R1_free) - beta_S * S - p * k_pair * S * M
dM/dt = alpha_M * LR_fraction - beta_M * M - k_pair * S * M
dO/dt = beta_translation * M - (delta_O + dilution) * O
```

旧版遗漏了第一条方程中的`-p*k_pair*S*M`，已经修正。`beta_S`和`beta_M`采用文献中的有效RNA周转率，不再额外加入生长稀释；蛋白读出仍包含稀释。

## 基准参数

|参数|值|来源或理由|
|---|---:|---|
|`srna_max_tx_uM_h`|0.06 µM/h|Levine计算使用的1 nM/min参考量级；并非文献上限|
|`output_max_tx_uM_h`|0.03 µM/h|相当于0.5 nM/min；Levine给出约1 nM/min量级，本项目采用0.1-1 nM/min保守扫描区间|
|`srna_decay_h`|1.2 h⁻¹|Levine参数估计小节取约0.02 min⁻¹，底层半衰期依据Massé 2003与Moll 2003；`0.02×60=1.2`|
|`output_mrna_decay_h`|6 h⁻¹|Levine由Massé 2003的`sodB`半衰期估计约0.1 min⁻¹；`0.1×60=6`|
|`pair_rate_uM_h`|1200 µM⁻¹h⁻¹|Levine根据Massé 2003推得约0.02 nM⁻¹min⁻¹；`0.02×1000×60=1200`|
|`srna_cocodegradation_fraction`|1|Levine理想一对一模型情形；不是实验辨识常数|
|`n_shutdown`|3.1|Shin 2020 B3-BM3R1 gate|
|`promoter_leak`|0.005/0.6|Shin 2020 B3-BM3R1 `ymin/ymax`|
|`threshold_ratio_to_rdf`|1|设计假设：复用RDF的BM3R1 operator|
|`copy_ratio_to_counter`|1|设计假设：每个counter plasmid一个cassette|
|`operator_load_stoich`|2|TetR家族同源二聚体机制假设|
|`output_translation_h`|1 h⁻¹|归一化比例，不作绝对蛋白预测|
|`output_protein_decay_h`|0 h⁻¹|稳定GFP读出假设|

每次核心运行会生成`00_parameter_provenance.csv`，逐项记录证据等级、来源和换算。`metadata.json`中也保存同一份审计信息。

## 上游接口

模型读取：

`rdfmodel_new/docs/week4_unloaded/unloaded_C31_translation_trajectory.csv`

计数器输入为：

```text
Int production(t) = 0.45 * actual_C31_translation_flux_uM_h(t)
```

计数器反应网络来自Zhao等人（2019）；当前BM3R1翻译、RDF翻译和Integrase标签参数是`rdfmodel_new`筛选出的工程工作点，不是Zhao原文常数。

## 运行

```powershell
python verify_model.py
python scan_interface.py
python run_analysis.py
python run_extended_analysis.py
python build_wiki_assets.py
```

`verify_model.py`除检查数值稳定性、DNA守恒和一次输入一次反转外，还显式检查每个`k_pair*S*M`通量会同时消耗一份mRNA和`p`份sRNA。

核心结果写入`outputs/shutdown_analysis_YYYYMMDD_HHMMSS/`，扩展结果写入`outputs/shutdown_extended_YYYYMMDD_HHMMSS/`。`latest_run.txt`与`latest_extended_run.txt`指向当前有效结果。

## 结果解释

文献锚定基准下，新生输出FWHM约6.98 h，LR后期产生约为峰值17.0%，并维持一次输入一次完整反转。修正RyhB范围后，阈值-强度扫描有91/169组通过；转录供给相图有124/225组通过，最低可行$\alpha_S/\alpha_M$约1.39、边界中位数约2.68；512点文献量级与明示假设包络中有54.1%通过。通过比例是所选参数包络的覆盖结果，不是wet-lab成功概率。

基准中的无Shutdown FWHM 10.61 h来自当前单bit交替轨迹，只覆盖约一个LR驻留间隔，不是多bit计数器的表达上限。多bit目标状态若保持多个输入周期，无Shutdown表达可能远长于10.61 h；准确时长须在实验组确定计数位、状态译码和复位逻辑后再计算。

当前sRNA参考为*E. coli* K-12 MG1655 RyhB（NCBI Gene 2847761，当前注释95 nt）。Levine实验使用`+1..+96`克隆片段；这一边界差异和完整序列审查记录在主报告中。ODE不读取碱基序列，固定的配对速率只对RyhB与相容`sodB`调控片段构成的代理系统有依据。

sRNA直接控制mRNA和新的翻译通量，不会主动清除已有稳定蛋白。若指定的最终蛋白需要快速下降，必须为该蛋白加入有来源的主动降解参数。

## Wiki版交付

`wiki_report/BM3R1驱动的sRNA_Shutdown模块_Wiki版报告.md`按Wiki叙事重新组织了问题、设计、方程、验证、结果、设计决策与局限性，并配套可直接上传的`wiki_report/assets/`图片。`build_wiki_assets.py`生成机制图和Levine阈值行为实现基准；该基准检验数值ODE与稳态解析解的一致性，不等同于wet-lab验证。

当前模型验证的是单个PB/LR bit局部状态变化触发的Shutdown动力学。完整计数器在第几个计数状态启用该模块、采用何种译码或接线方式，保留为与实验组共同确认的接口问题；确认前不把现有结果表述为完整的count-to-N输出系统。

## 主要来源

- Zhao et al. 2019: https://doi.org/10.1093/nar/gkz245
- Levine et al. 2007: https://doi.org/10.1371/journal.pbio.0050229
- Zhang et al. 2005: https://doi.org/10.1128/JB.187.3.980-990.2005
- Selinger et al. 2000: https://doi.org/10.1038/82367
- Moll et al. 2003: https://doi.org/10.1261/rna.5850703
- Shin et al. 2020: https://doi.org/10.15252/msb.20199401
- Hao et al. 2011: https://doi.org/10.1073/pnas.1100432108
- NCBI RyhB Gene 2847761: https://www.ncbi.nlm.nih.gov/gene/2847761
