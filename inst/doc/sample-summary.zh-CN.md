# 样本数量汇总

PanFamR 0.3.2 新增样本数量汇总接口。项目设置完成后，即可按癌种查看输入数据中的正常与肿瘤样本数，无需先计算 ssGSEA，也无需自行加载表达矩阵。

## 1. 设置项目

安装方法见项目的 [中文入门教程](https://github.com/JianGao0852/PanFamR/blob/main/README.zh-CN.md)。在 RStudio 中修改目录和基因列表，选中需要的代码后点击 **Run**：

```r
library(PanFamR)

project <- panfam_setup(
  work_dir = "D:/PanFamR_analysis",
  output_dir = "results",
  genes = c("TP53", "MDM2", "CDKN1A"),
  family_name = "MyFamily",
  cancers = c("BRCA", "LUAD")
)
```

`cancers` 决定汇总哪些癌种。设为 `NULL` 时保留所选输入中的全部癌种。样本汇总不依赖家族成员的匹配结果；基因列表供后续基因家族分析使用。

## 2. 选择输入并运行

### TCGA 输入，或项目的自带表达表

```r
sample_result <- panfam_sample_summary(project)
```

默认使用 `data_source = "tcga"`：未在设置时提供 `expression`，则自动读取内置 TCGA `all_tpm`；提供了 `expression = my_expression`，则汇总该表。两种情况都会遵循项目中的癌种筛选。

### 内置 TCGA/GTEx 整合输入

```r
tcga_gtex_result <- panfam_sample_summary(project, data_source = "tcga_gtex")
```

这个选项始终读取内置 TCGA/GTEx 整合表，与项目设置中的自带表达表独立。源表的 `project` 和 `sample_type` 分别提供癌种与分组信息；`TCGA_tumor` 计入 `Tumor`，`TCGA_normal` 和 `GTEx_normal` 合并计入 `Normal`。

这两个调用各汇总所选输入一次。即使项目定义了多个基因家族，也不会为每个家族重复统计。

## 3. 查看和保存结果

运行时会直接打印不带行号的汇总表；同时返回可赋值保存的模块结果：

```r
# 数量表：Cancer、Normal、Tumor 三列，末行为 Total。
sample_result$tables$Project$sample_summary

# 数据来源、计数单位、所选癌种等记录。
sample_result$tables$Project$sample_summary_metadata

# 保存地址与运行记录。
sample_result$output_dirs[["Project"]]
sample_result$run_log
panfam_status(project)
```

样本汇总保存在项目级目录，独立于各基因家族的分析文件夹：

```text
results/
  sample_summary_tcga/
    sample_summary.csv
    sample_summary_metadata.csv
    result.rds
  sample_summary_tcga_gtex/
    sample_summary.csv
    sample_summary_metadata.csv
    result.rds
  MyFamily/
    ...
```

这里只生成统计表，不生成图片。再次调用会新建 `sample_summary_tcga_002` 等带编号的目录，保留先前结果。运行记录中的家族字段为 `Project`，模块名为 `sample_summary_tcga` 或 `sample_summary_tcga_gtex`。`panfam_report(project)` 生成的 HTML 索引也会包含已运行的汇总步骤。

## 4. 已有 tpm 对象时的兼容接口

如果已经在 R 中读入含 `Cancer` 和 `Group` 两列的数据框，可以直接使用原函数名：

```r
tpm <- data.frame(
  Cancer = c("BRCA", "BRCA", "LUAD", "LUAD"),
  Group = c("Tumor", "Tumor", "Normal", "Tumor")
)

counts <- get_all_tpmsample(tpm)
```

上面的示例数据对应以下布局：

```text
Cancer Normal Tumor
  BRCA      0     2
  LUAD      1     1
 Total      1     3
```

`get_all_tpmsample()` 会打印数量表，并不可见返回该数据框，便于赋值后继续处理。只使用 `Cancer`、`Group` 两列，其他列不参与统计；某一分组不存在时补零。这个基础接口不自动加载数据、不保存文件，也不记录项目历史。初次使用者可以直接使用前面的 `panfam_sample_summary()`。

## 5. 计数口径

- 统计单位是输入数据的样本行，不是独立患者；不会按患者 ID 去重。
- `Normal`、`Tumor` 始终作为独立列，`Total` 是末尾的总计行。
- 分组标签区分大小写，只接受 `Normal`、`Tumor`；癌种或分组缺失、空白，或存在未知分组时不会静默漏计。独立接口会报错，项目接口会发出警告、将运行状态记为 `error` 并保存 `error.txt`，该次不生成数量表。`Cancer = "Total"` 为总计行保留，不能作为输入癌种。
- 此步骤不进行基因匹配、评分、生存结局匹配或各下游分析的有效性筛选。各模块最终使用多少样本，应查看对应结果表。

完成汇总后，可以继续使用同一个项目逐步运行需要的分析：

```r
scores <- panfam_ssgsea(project)
expression_result <- panfam_expression(project)
```
