# PanFamR 中文分步骤入门

安装好 PanFamR 后，先设置工作目录、结果目录和自定义基因家族，再逐一调用需要的分析函数。每个函数自动读取内置数据，生成这一部分的结果表和图片。

## 1. 安装一次

需要 R 4.3 或更新版本：

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install(c("GSVA", "BiocParallel", "ComplexHeatmap"), ask = FALSE)
options(timeout = 3600)
remotes::install_github("JianGao0852/PanFamR", upgrade = "never")
```

包内置输入约 1.29 GiB，建议预留约 6 GB 磁盘空间用于下载、安装和缓存。TCGA/GTEx 表达矩阵、配对表达数据及 Braun 免疫治疗表达 RDS 等大文件以不超过 20 MiB 的分片随包保存，所有分片数据在首次需要时自动重组到本地缓存目录。安装完成后，分析使用本地数据，无需再下载分析数据。

## 2. 设置路径和基因家族

修改以下工作路径和基因符号，运行这一段：

```r
library(PanFamR)

work_dir <- "D:/PanFamR_analysis"
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

project <- panfam_setup(
  work_dir = getwd(),
  output_dir = "results",
  genes = c("COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1"),
  family_name = "MyFamily",
  correlation_method = "pearson"
)
```

`output_dir = "results"` 表示结果保存在工作目录下的 `results` 文件夹，也可以写完整地址，例如 `"D:/MyResults"`。`family_name` 是你给基因家族起的名称，可以使用中文。`project` 保存这一轮分析的设置，后续函数都传入同一个 `project`。

设置函数只准备目录、记录参数，不执行分析，也不修改当前工作目录。基因无需具有相同前缀，也无需属于原来预设的家族。基因符号区分大小写；未匹配的基因会列入覆盖表，整个家族都匹配不到时会提示错误。默认至少匹配一个基因，可用 `min_size` 调整。

相关分析默认使用原始分析脚本的 Pearson 方法。若研究设计需要 Spearman 方法，在设置时改为 `correlation_method = "spearman"`。

## 3. 逐一运行需要的部分

在 RStudio 中选中需要运行的那一行或那一段，点击 **Run**，或按 **Ctrl+Enter**。

```r
# 第一步：计算基因家族的 ssGSEA 评分。
scores <- panfam_ssgsea(project)

# 第二步：得到 TCGA 基因表达差异和家族评分差异。
expression_result <- panfam_expression(project)

# 第三步：生存分析，这里只运行 OS。
survival_result <- panfam_survival(project, endpoints = "OS")

# 第四步：TMB 关联。
tmb_result <- panfam_tmb(project)

# 第五步：MSI 关联。
msi_result <- panfam_msi(project)
```

每运行一个分析函数，便保存这一部分的结果。需要评分的分析会自动检查并计算评分，之后在这个 `project` 中复用。基因、路径和原始数据不需要在每一步重复填写或手工加载。

### 各部分对应的函数

| 需要的分析 | 运行这一行 |
| --- | --- |
| ssGSEA 基因家族评分 | `panfam_ssgsea(project)` |
| TCGA 基因表达差异和家族评分差异 | `panfam_expression(project)` |
| 仅肿瘤样本的评分分布 | `panfam_tumor_scores(project)` |
| 肿瘤/正常评分探索性 ROC | `panfam_roc(project)` |
| TCGA/GTEx 基因表达比较 | `panfam_tcga_gtex(project)` |
| 同一患者的肿瘤/正常配对评分比较 | `panfam_paired(project)` |
| OS、DSS、PFI、DFI 生存分析 | `panfam_survival(project, endpoints = c("OS", "DSS", "PFI", "DFI"))` |
| 家族中各个基因的生存分析 | `panfam_gene_survival(project, endpoints = "OS")` |
| TMB 关联 | `panfam_tmb(project)` |
| MSI 关联 | `panfam_msi(project)` |
| 免疫细胞关联 | `panfam_immune_cells(project)` |
| 免疫评分关联 | `panfam_immune_scores(project)` |
| 肿瘤微环境通路关联 | `panfam_tme(project)` |
| Hallmark 通路关联 | `panfam_hallmark(project)` |
| 免疫检查点基因关联 | `panfam_checkpoints(project)` |
| 趋化因子关联 | `panfam_chemokines(project)` |
| 趋化因子受体关联 | `panfam_chemokine_receptors(project)` |
| 免疫刺激因子关联 | `panfam_immune_stimulators(project)` |
| 免疫抑制因子关联 | `panfam_immune_inhibitors(project)` |
| 四个免疫治疗队列分析 | `panfam_immunotherapy(project)` |
| 已内置基因变异表的分析 | `panfam_alterations(project)` |
| 已内置 CNV 关联表的分析 | `panfam_cnv(project)` |
| 已内置甲基化关联表的分析 | `panfam_methylation(project)` |

生存分析可以只选一个或几个结局，例如 `endpoints = c("OS", "PFI")`。选择自己需要的函数逐步运行即可。

## 4. 查看每一步的结果

输出按“家族 → 分析部分”组织：

```text
results/
  MyFamily/
    ssgsea/
    expression/
    survival_OS/
    tmb/
    ...
```

各部分文件夹中保存可用的 CSV 结果表、PDF 图片和 300 dpi 白底 PNG 图片。重复运行同一部分会生成带编号的新文件夹，保留先前结果；在已有分析设置的结果目录再次创建项目时，会生成带时间戳的相邻目录。实际结果根目录可以通过 `project$output_dir` 查看。

也可以直接在 R 中查看：

```r
# 查看评分差异表。
expression_result$tables[["MyFamily"]][["score_differences"]]

# 查看这一部分的状态和保存地址。
expression_result$run_log
expression_result$output_dirs

# 查看已经运行过哪些部分。
panfam_status(project)

# 汇总已运行步骤，生成并打开 HTML 结果索引。
report_file <- panfam_report(project)
utils::browseURL(report_file)
```

生成索引不会启动其他分析。`panfam_ssgsea()` 返回各家族的样本评分表；其他分析返回的结果对象包含 `tables`、`run_log`、`output_dirs` 和 `module`。

## 同时定义多个家族、选择癌种

这些参数仍然只在设置时填写，后续函数的调用方式相同：

```r
families <- list(
  FamilyA = c("COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1"),
  MySignature = c("TP53", "MDM2", "CDKN1A", "BAX", "BCL2")
)
project <- panfam_setup(
  output_dir = "two_families",
  gene_sets = families,
  cancers = c("BRCA", "LUAD")
)
scores <- panfam_ssgsea(project)
panfam_survival(project, endpoints = c("OS", "PFI"))
```

单个家族用 `genes`，多个家族用带名称的 `gene_sets`，二者选择其一。家族名可以使用中文，程序自动生成不同的安全目录名。要更换家族、癌种或其他参数时，重新调用 `panfam_setup()` 创建项目。完整矩阵已在 32 GB 内存的 Windows 电脑上实跑；先选少量癌种可缩小计算规模。

## 数据覆盖和分析规则

### 0.3.1 恢复原始绘图风格和主要估计量

自定义基因集合、分模块调用和内置数据接口继续保留。各模块恢复以下原始绘图形式：

| 分析模块 | 结果图及其含义 |
| --- | --- |
| 基因表达与 TCGA/GTEx 比较 | 六档热图，以来源表达尺度的 `mean(tumor) - mean(normal)` 着色并显示数值，分档阈值为 −2、−1、0、1、2；原始 P >= 0.05 的格子为白色。原脚本称该数值为 `logFC`，不能直接称为 log2 fold change。 |
| TMB、MSI | 红/青配色雷达图，展示所选相关系数。默认 Pearson，可在设置时选择 Spearman。 |
| 免疫细胞、免疫评分、通路及标志基因面板 | 蓝–白–红 ComplexHeatmap 热图，星号依据原始 P 值；CSV 中另保留 BH 校正值。 |
| 家族生存分析 | High 对 Low 中位数组的表格式森林图，以及黄/蓝 `ggsurvplot` 曲线，标注 HR、置信区间和 P 值。已知以天为单位的随访在图中按 365 天/年显示，来源单位不明时保留来源标签。 |
| 成员基因生存分析 | 原始 Risky/Protective 分类热图及计数差条形注释，依据基因表达 Cox HR 的方向和原始 P < 0.05 分类。 |
| 基因变异 | 计数热图及行、列变异类型堆叠注释。基因/癌种格子统计变异事件数，同一样本的不同变异类型可分别计数；Total 行统计自定义输入家族中至少有一项变异的唯一样本数。频率表与覆盖表继续保留。 |
| CNV、甲基化 | 气泡图展示源表相关系数，边框编码源表 FDR；不重新计算源表 FDR。 |
| 免疫治疗 | 在各队列内按评分真实中位数分为 High/Low，绘制来源可用的 Yes/No/SD 应答组成堆叠图及配套 KM 曲线。 |

家族 Cox 主图采用 `hr_high_vs_low`、`group_ci_lower`、`group_ci_upper`、`group_p_value` 和 `group_fdr`，表示 **Low：评分 <= 中位数；High：评分 > 中位数** 的组间关联。附加连续模型仍保留 `hr_per_sd`、`ci_lower`、`ci_upper`、`p_value` 和 `fdr`，其中 HR 对应评分增加一个标准差。成员基因 Cox 表另提供原始表达单位的 `hr_per_unit`、`unit_ci_lower`、`unit_ci_upper` 和 `unit_p_value`。解读时应先确认图注和字段对应的模型。

新计算的统计结果在每张结果表的有效检验内进行 BH 多重检验校正，分组模型与连续模型分别校正。图中沿用原始 P 值的着色或星号会明确标注，不能替代校正值。来源 CNV/甲基化表保留原始 FDR。样本、事件或数据不足时会在结果和状态中注明。

恢复原图不表示所有数值与继承脚本的每一次输出完全相同：程序不手工移动中位数组中的样本，不反转 High/Low 标签。ROC 继续保留固定方向的经验曲线，不进行平滑，与原脚本的平滑 ROC 有差别。比较使用同一 TCGA 数据中的肿瘤和正常样本，固定以高评分预测肿瘤，因此 AUC 可以低于 0.5。这是探索性比较，未经过外部队列验证。

现有变异、CNV 和甲基化输入只覆盖原项目的部分基因。自定义家族按基因符号检查覆盖范围，无法支持的结果会显示 `skipped` 或 `partial`，这些状态不能解释为“没有生物学关联”。TCGA/GTEx 比较使用来源标记的正常样本，未校正队列或批次差异。

免疫治疗分析优先使用已标注的治疗前样本，治疗时点缺失的样本予以保留并记录，排除已知非免疫治疗对照，每位患者保留一个样本。应答组成图保留来源可用类别，与另行输出的应答者/非应答者比较表含义不同。GSE135222 使用 PFS，其余三个队列使用 OS。GEO 未提供 GSE135222 的 PFS 时间单位，因此保留来源单位。未知结局保留为缺失；非有限表达值的基因会被排除并记录。

## 其他接口和发布准备

`panfam_data_files()` 可以列出内置文件。需要自带表达数据时，在设置中传入 `expression = my_expression`：每行一个样本，行名为唯一样本 ID，包含 `Cancer`、`Group` 和以基因符号命名的数值表达列；`Group` 取值为 `Tumor` 或 `Normal`。内置临床等辅助数据只在样本 ID 能匹配时可用。

完整批量运行仍可使用 `panfam_analyze()`；分步骤调用不需要它。可直接编辑 `inst/scripts/custom-gene-family-analysis.R`，按其中标注逐步运行。

维护者为 Serena，联系邮箱为 <g512827252@outlook.com>。使用问题可在 [GitHub Issues](https://github.com/JianGao0852/PanFamR/issues) 反馈。`inst/extdata/README-data.md` 中的数据来源、引用和再分发条款说明仍在补充。MIT 许可证适用于代码，不替代各数据集自身的条款。MSigDB GMT 文件未内置，自定义基因向量可以直接使用。
