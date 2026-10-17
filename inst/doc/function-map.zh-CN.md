# 原函数与模块接口对照

原分析脚本中的 34 个顶层函数名称均保留并导出。新接口按分析任务组织：通过 `panfam_setup()` 一次设置路径、基因家族和癌种，再把同一个 `project` 传给所需模块。旧函数名与模块名不同，不表示相应分析被删除。

下面的对照说明分析用途及主要输出。部分旧绘图变体仍通过原函数提供，不能据此认定两套接口的参数、统计处理及所有图形完全相同。新用户建议从模块接口开始；已有脚本直接调用旧函数时，应查看该函数当前的帮助与参数。

## 1. 先配置项目，再选择模块

```r
library(PanFamR)

project <- panfam_setup(
  work_dir = "D:/PanFamR_analysis",
  output_dir = "results",
  genes = c("COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1"),
  family_name = "MyFamily",
  cancers = c("BRCA", "LUAD")
)

# 可选：先汇总输入样本数量，不需要评分。
sample_result <- panfam_sample_summary(project)

# 计算评分，随后按需要逐项运行。
scores <- panfam_ssgsea(project)
```

只分析一个癌种时，在设置中写 `cancers = "LUAD"`。模块函数使用项目中的癌种筛选，无需逐次填写 `cancer` 参数。要改变癌种或基因列表，重新创建项目。

癌种在评分前筛选，改变所选癌种也可能改变家族评分的归一化范围；不同癌种选择项目的分值应结合各自评分范围解释。

## 2. 免疫和通路面板：八个原函数的对应关系

下表中的相关性结果保存在 `结果目录/家族名/模块名/`，对应 PDF 和 PNG 在该目录下的 `plots/` 中。

| 原函数 | 推荐模块调用 | 模块目录 | 主要 CSV 表 |
| --- | --- | --- | --- |
| `ssGsea_Score_pan_cancer_chemokine(df, tpm)` | `panfam_chemokines(project)` | `chemokines` | `chemokines_correlations.csv`、`chemokines_gene_coverage.csv` |
| `ssGsea_Score_pan_cancer_chemokine_receptor(df, tpm)` | `panfam_chemokine_receptors(project)` | `chemokine_receptors` | `chemokine_receptors_correlations.csv`、`chemokine_receptors_gene_coverage.csv` |
| `ssGsea_Score_pan_cancer_immune_stimulator(df, tpm)` | `panfam_immune_stimulators(project)` | `immune_stimulators` | `immune_stimulators_correlations.csv`、`immune_stimulators_gene_coverage.csv` |
| `ssGsea_Score_pan_cancer_immune_inhibitor(df, tpm)` | `panfam_immune_inhibitors(project)` | `immune_inhibitors` | `immune_inhibitors_correlations.csv`、`immune_inhibitors_gene_coverage.csv` |
| `ssGsea_Score_pan_cancer_immunecellratio(df)` | `panfam_immune_cells(project)` | `immune_cells` | `immune_cells_correlations.csv` |
| `ssGsea_Score_pan_cancer_immunescore(df)` | `panfam_immune_scores(project)` | `immune_scores` | `immune_scores_correlations.csv` |
| `ssGsea_Score_pan_cancer_TME_GSVA(df)` | `panfam_tme(project)` | `tme` | `tme_pathways_correlations.csv` |
| `ssGsea_Score_pan_cancer_HALLMARK_GSVA(df)` | `panfam_hallmark(project)` | `hallmark` | `hallmark_pathways_correlations.csv` |

模块保留相应的标志基因列表及内置辅助表，自动取得家族评分。相关性表包含 `Panel`、`Cancer`、`Feature`、`n`、`correlation`、`p_value`、`reason` 和 `fdr` 等字段；标志基因面板另有基因覆盖表。图形使用相关表的同名文件前缀。

默认相关方法为 Pearson。需要 Spearman 时，在项目设置中指定 `correlation_method = "spearman"`。模块会按肿瘤患者去重、匹配样本，并保留 BH 校正结果；这些处理应结合 README 的统计规则解读，不能把全部数值视为与原脚本每次运行完全一致。

## 3. 四个单癌函数的调用方式

| 原函数 | 推荐模块调用 | 对应内容 |
| --- | --- | --- |
| `ssGsea_Score_Tumor_Normal(df, cancer, tpm = NULL)` | `panfam_expression(project)` | 所选癌种中 Normal 与 Tumor 的家族评分比较，独立图为 `plots/score_tumor_normal_<Cancer>.pdf` 及 PNG，同时保留多癌种图和成员基因表达结果。这里的分组是正常/肿瘤，不是 Stage、Grade 或性别分析。 |
| `ssGsea_Score_one_paired_boxplot(df, cancer)` | `panfam_paired(project)` | 患者内肿瘤/正常配对评分比较；包括 `paired_samples.csv`、`paired_score_differences.csv`，以及 `plots/paired_scores_<Cancer>.pdf` 和 PNG。 |
| `ssGsea_Score_ROC(df, cancer, tpm = NULL)` | `panfam_roc(project)` | 所选癌种的肿瘤/正常评分 ROC；包括 `roc_summary.csv`、`roc_curve.csv` 和 `plots/tumor_normal_roc_<Cancer>.pdf` 及 PNG。 |
| `ssGsea_Score_KM(df, cancer)` | `panfam_survival(project, endpoints = "OS")` | 所选癌种的家族评分 OS 分析；包括 `survival_cox.csv`、`survival_km.csv`、`survival_samples.csv`、`plots/survival_KM/OS_<Cancer>.pdf` 及 PNG，并另存带 `_with_risk_table` 后缀的风险表图。 |

例如，只运行 LUAD 的四项分析：

```r
project <- panfam_setup(
  output_dir = "LUAD_results",
  genes = c("COX4I1", "COX5A", "COX5B"),
  family_name = "MyFamily",
  cancers = "LUAD"
)

expression_result <- panfam_expression(project)
paired_result <- panfam_paired(project)
roc_result <- panfam_roc(project)
survival_result <- panfam_survival(project, endpoints = "OS")
```

需要评分的模块会自动计算并复用评分，不必为每项手工准备 `df` 或 `tpm`。生存模块还支持 DSS、PFI、DFI。KM 主图和伴随风险表图使用相同拟合及分组，不重新估计模型。单癌评分图要求 Normal/Tumor 两组均有有限评分，缺少其中一组时不生成该图；图中复用 `score_differences` 表中的统计结果，样本不足或无法估计 P 值时不添加显著性标记。ROC 模块采用固定“较高评分预测肿瘤”的方向和未平滑经验曲线，保留低于 0.5 的 AUC；这与旧 ROC 的自动方向和平滑设置不同，也不代表外部验证。

直接使用原 `ssGsea_Score_Tumor_Normal(df, cancer)` 或 `ssGsea_Score_ROC(df, cancer)` 时，第三个参数可省略，此时分组直接取自评分表自身的 `Cancer`、`Group`，不会额外加载完整表达矩阵。显式传入 `tpm` 时，按样本 ID 精确匹配元数据。

## 4. 其余原函数的对应用途

| 原函数 | 模块接口或保留用法 | 主要内容 |
| --- | --- | --- |
| `ssGsea_Score(tpm, cox_genes, ...)` | `panfam_ssgsea(project)`；也可直接使用 `panfam_score()` | 家族 ssGSEA 评分。项目接口自动把评分及覆盖表保存到 `ssgsea/`。旧接口需要指定独立评分文件时，可显式传入 `output_file`。 |
| `plot_diff_heatmap(rt)` | `panfam_expression(project)` | 成员基因的肿瘤/正常表达比较；包括 `gene_differences.csv` 及文件名前缀为 `gene_expression_differences` 的均值差热图。 |
| `plot_diff_heatmap1(...)` | `panfam_tcga_gtex(project)` | 内置 TCGA/GTEx 表达比较；当前原接口需要基因向量及整合表达表，并非仅接收一个目录名。 |
| `plot_OS_KM_expression(rt, ...)` | `panfam_gene_survival(project, endpoints = "OS")` | 成员基因 OS 分析，包括 Cox、KM 统计表和基因覆盖信息。 |
| `plot_cox_analysis_forest(rt_cox)`、`plot_cox_analysis_forest_picture(folder_path)` | `panfam_gene_survival(project, endpoints = "OS")`；原表输入/森林图变体继续通过旧接口使用 | 新模块汇总成员基因的 Cox 统计与分类热图，不要求中途手工生成、读取每个基因的 `.cox.txt` 文件。 |
| `plot_cox_heatmap1(folder_path)`、`plot_cox_heatmap2(folder_path)` | `panfam_gene_survival(project, endpoints = "OS")` | 成员基因 Cox 结果的 Risky/Protective 分类及汇总热图；旧目录输入绘图变体仍保留。 |
| `ssGsea_Score_diff(df)`、`ssGsea_Score_diff1(df)` | `panfam_expression(project)` | 多癌种的家族评分肿瘤/正常差异与分布图。 |
| `ssGsea_Score_diff_Tumor(df)` | `panfam_tumor_scores(project)` | 仅肿瘤样本的家族评分汇总和分布，模块按患者去重。 |
| `ssGsea_Score_pan_cancer_TMB(df)` | `panfam_tmb(project)` | 家族评分与 TMB 的关联及雷达图。 |
| `ssGsea_Score_pan_cancer_MSI(df)` | `panfam_msi(project)` | 家族评分与 MSI 的关联及雷达图。 |
| `ssGsea_Score_pan_cancer_ICG(df, tpm)` | `panfam_checkpoints(project)` | 家族评分与免疫检查点基因表达的关联。 |
| `ssGsea_Score_pan_cancer_fores(df, Target)` | `panfam_survival(project, endpoints = "OS")` 或选择其他终点 | 家族评分生存关联、森林图和 KM 曲线；主图使用真实中位数组的 High/Low 模型。 |
| `ssGsea_Score_paired_boxplot(df)` | `panfam_paired(project)` | 多癌种与逐癌种的患者内配对评分比较。 |
| `Alterations(folder_path)` | `panfam_alterations(project)` | 自定义家族在现有变异表覆盖范围内的事件、频率、覆盖表与计数热图。 |
| `CNV_Methylation_mRNA(analysis_target)`、`CNV_Methylation_mRNA1(analysis_target)` | `panfam_cnv(project)` 或 `panfam_methylation(project)` | 分别检索现有 CNV/甲基化相关表，输出相关性、源 FDR、覆盖表及气泡图；原绘图变体仍保留。 |
| `plot_km(file_path)`、`plot_km1(file_path)`、`plot_km2(file_path)` | `panfam_immunotherapy(project)` 提供家族评分的同类分析；其他旧变体继续通过原接口使用 | 项目接口完成四队列家族评分、可用生存分析、伴随风险表图与应答组成图。旧函数中基因表达 KM 等变体仍需各自要求的兼容输入，不能等同于家族评分分析。 |

补充的 `run_ssGSEA(tpm, gene_sets, output_dir)` 为直接批量评分接口；项目工作流使用 `panfam_ssgsea(project)`。样本数量统计另见[样本数量汇总教程](sample-summary.zh-CN.md)：基础接口为 `get_all_tpmsample(tpm)`，自动读取并保存项目结果的接口为 `panfam_sample_summary(project)`。

模块对样本 ID 和患者进行匹配及适用的去重处理，新计算的检验在对应结果表中保留 BH 校正值；未知响应或事件标签保留为缺失，并按各分析所需信息筛选。真实中位数分组不会通过手工移动样本或反转标签调整。这些规则与保留旧函数的部分处理方式存在差别，因此不能保证新旧结果逐值相同。

## 5. 查看结果位置

家族级模块返回 `panfam_step`，通过以下字段查看对应统计表和目录：

```r
expression_result$tables[["MyFamily"]]
expression_result$output_dirs[["MyFamily"]]
expression_result$run_log
```

样本汇总是项目级模块，其数量表位于 `sample_result$tables$Project$sample_summary`。已运行步骤可用 `panfam_status(project)` 查看，也可用 `panfam_report(project)` 生成 HTML 结果索引。
