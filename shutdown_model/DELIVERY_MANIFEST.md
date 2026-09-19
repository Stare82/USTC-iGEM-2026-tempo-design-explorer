# Shutdown modeling delivery manifest

交付版本：`3.3`  
交付日期：`2026-09-09`

## 正式文件

- `BM3R1驱动的sRNA Shutdown模块建模报告.pdf`
- `BM3R1驱动的sRNA Shutdown模块建模报告.md.zip`
- `report_source/BM3R1驱动的sRNA Shutdown模块建模报告.md`
- `report_source/Shutdown_report.tex`
- `BM3R1驱动的sRNA_Shutdown模块_Wiki版报告.pdf`
- `BM3R1驱动的sRNA_Shutdown模块_Wiki版报告.md.zip`
- `wiki_report/BM3R1驱动的sRNA_Shutdown模块_Wiki版报告.md`
- `wiki_report/Wiki_report.tex`
- `wiki_report/assets/`
- `build_wiki_assets.py`
- `README.md`
- `requirements.txt`

工作区根目录的`shutdown_model.zip`是可复现总包：除本目录外，还包含运行所需的最小上游依赖`rdfmodel_new/model/zhao_core.py`、冻结C31轨迹，以及`oscillator/code/Oscillator_Engineering_Sensitivity_Analysis.py`。解压后须保留`shutdown_model`、`rdfmodel_new`和`oscillator`三个同级目录。

## 当前有效结果

- 核心分析：`outputs/shutdown_analysis_20260909_150254/`
- 扩展分析：`outputs/shutdown_extended_20260909_133421/`
- 参数溯源：核心分析目录中的`00_parameter_provenance.csv`
- 参数包络：扩展目录中的`07_parameter_ranges.csv`和`07_literature_range_samples.csv`

`latest_run.txt`与`latest_extended_run.txt`指向上述结果。旧版6.26 h、99.2%结果、等转录率失败基准以及使用固定无Shutdown峰值对照的统计均不再作为有效交付。

## 关键修正

1. sRNA方程加入`-p*k_pair*S*M`，与mRNA共同消耗；
2. RNA参数改为Levine 2007 RyhB/sodB文献代理，并明确区分文献直接值、文献量级与项目扫描包络；
3. BM3R1 Hill/leak采用Shin 2020 B3 gate数据，阈值使用相对RDF的无量纲比例；
4. 每个参数均标注证据等级、来源与理由；
5. 新增$\alpha_S$-$\alpha_M$转录供给相图，直接给出可行供给边界；
6. 联合分析改为512点拉丁超立方文献量级与明示假设包络，并为每个靶mRNA转录率使用匹配的无Shutdown对照；
7. 采用NCBI当前95 nt RyhB参考序列作为模型注释，并明确Levine `+1..+96`构建边界仍需原质粒序列确认。
8. 新增Levine阈值响应实现基准，ODE稳态与同一质量作用方程的解析稳态最大绝对误差约`1.04e-10 nM`；这验证数值实现，而非实验可行性。
9. 新增Wiki版叙事、机制图、精确定义的输出指标与“模型结果到设计决策”表，并将完整计数逻辑的接入位置明确列为待实验组确认项。
10. 将所有换算参数补充为“论文原值—换算式—论文位置—底层实验来源”四层证据，并明确10.61 h仅为当前单bit无Shutdown对照，不是多bit系统的表达上限。

## 当前结论

- 文献锚定基准的新生输出FWHM：`6.98 h`
- LR后期新生输出：峰值的`17.0%`
- 计数器一次输入一次完整反转：`5/5`周期
- 二维扫描通过：`91/169`，可行脉宽`5.42-7.36 h`
- 转录供给相图通过：`124/225`
- 最低可行供给比/边界中位数：`1.39/2.68`
- 参数包络达到工程标准：`277/512 (54.1%)`

这些是RyhB/sodB代理、当前BM3R1计数器工作点与明示工程判据下的条件性结果，不是wet-lab成功率或固定生物学常数。

## 尚待实验组确认

完整计数器应在第几个状态启用Shutdown，以及采用状态特异启动子、逻辑译码器或其他接线方式，尚未冻结。现阶段可交付的是单bit局部Shutdown模块的机制可行性、参数调节规律和鲁棒性分析；在接口确定前，不主张已经完成完整count-to-N系统建模。

## 验证命令

```powershell
python -m py_compile shutdown_core.py run_analysis.py run_extended_analysis.py scan_interface.py verify_model.py
python verify_model.py
```
