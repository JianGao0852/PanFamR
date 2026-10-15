# PanFamR

PanFamR is an R package for step-by-step pan-cancer analysis of custom gene families and signatures. Set your working directory, result directory, and gene symbols once. Then run the function for each analysis you want. The package loads its bundled data and saves the corresponding tables and figures automatically.

中文分步骤入门见 [README.zh-CN.md](README.zh-CN.md). A ready-to-edit script is provided in [inst/scripts/custom-gene-family-analysis.R](inst/scripts/custom-gene-family-analysis.R).

## Install once

Requires R 4.3 or newer:

```r
install.packages(c("BiocManager", "remotes"))
BiocManager::install(c("GSVA", "BiocParallel", "ComplexHeatmap"), ask = FALSE)
options(timeout = 3600)
remotes::install_github("JianGao0852/PanFamR", upgrade = "never")
```

Bundled inputs occupy approximately 1.29 GiB. Allow about 6 GB of free disk space for download, installation, and cached data files. Large inputs, including the TCGA/GTEx expression matrices, paired expression data, and Braun immunotherapy expression RDS, are distributed in numbered parts of at most 20 MiB. All split inputs are reconstructed automatically in the local cache when first needed. Once installed, analyses use local data without downloading datasets. The full-matrix workflow has been run on a 32 GB RAM Windows computer; selecting fewer cancer types reduces calculation size.

## 1. Set paths and genes once

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

Replace the paths and gene symbols with your own. An absolute result path is also accepted. Relative result paths are resolved against `work_dir`. Setup records your settings and prepares the output directory; it does not run analyses or change the working directory. Keep using the same `project` object for subsequent steps.

Genes do not need a shared prefix or a predefined family name. Matching is exact and case sensitive. Missing symbols are recorded in coverage tables; a wholly unmatched family produces an error. The default `min_size = 1` permits small families; increase it to require more matched genes.

Pearson correlation is the default, matching the original analysis script. Use `correlation_method = "spearman"` during setup if your study requires Spearman correlation.

## 2. Run individual analysis functions

In RStudio, select the line or block you want and click **Run** or press **Ctrl+Enter**.

```r
# Score the specified family using the bundled TCGA expression data.
scores <- panfam_ssgsea(project)

# Compare gene expression and family scores between tumor and normal samples.
expression_result <- panfam_expression(project)

# Run survival analysis; choose the endpoints you need.
survival_result <- panfam_survival(project, endpoints = "OS")

# Each association is a separate step.
tmb_result <- panfam_tmb(project)
msi_result <- panfam_msi(project)
```

Each call runs the requested analysis and saves its results immediately. Steps that need ssGSEA scores obtain them automatically if necessary and reuse the scores in `project` afterward. You do not have to supply the genes, load input files, or calculate the same scores again for every step.

| Analysis | Function |
| --- | --- |
| ssGSEA family scoring | `panfam_ssgsea(project)` |
| TCGA expression and score differences | `panfam_expression(project)` |
| Tumor-only score distributions | `panfam_tumor_scores(project)` |
| Exploratory tumor/normal score ROC | `panfam_roc(project)` |
| TCGA/GTEx gene-expression comparison | `panfam_tcga_gtex(project)` |
| Matched tumor/normal score comparison | `panfam_paired(project)` |
| OS, DSS, PFI, or DFI survival | `panfam_survival(project, endpoints = c("OS", "DSS", "PFI", "DFI"))` |
| Individual-gene survival | `panfam_gene_survival(project, endpoints = "OS")` |
| TMB association | `panfam_tmb(project)` |
| MSI association | `panfam_msi(project)` |
| Immune-cell association | `panfam_immune_cells(project)` |
| Immune-score association | `panfam_immune_scores(project)` |
| Tumor microenvironment pathways | `panfam_tme(project)` |
| Hallmark pathways | `panfam_hallmark(project)` |
| Immune checkpoint genes | `panfam_checkpoints(project)` |
| Chemokine genes | `panfam_chemokines(project)` |
| Chemokine receptors | `panfam_chemokine_receptors(project)` |
| Immune stimulators | `panfam_immune_stimulators(project)` |
| Immune inhibitors | `panfam_immune_inhibitors(project)` |
| Four immunotherapy cohorts | `panfam_immunotherapy(project)` |
| Available gene-alteration tables | `panfam_alterations(project)` |
| Available CNV associations | `panfam_cnv(project)` |
| Available methylation associations | `panfam_methylation(project)` |

## 3. Find results and check progress

Results are organized by family and step:

```text
results/
  MyFamily/
    ssgsea/
    expression/
    survival_OS/
    tmb/
    ...
```

Step folders contain the available CSV tables, PDF figures, and PNG figures at 300 dpi with a white background. Repeating a step creates a new numbered folder, preserving its earlier output. Starting a new project in an already configured result directory uses a timestamped sibling directory. The actual root is available as `project$output_dir`.

```r
# Inspect a result directly in R.
expression_result$tables[["MyFamily"]][["score_differences"]]
expression_result$run_log
expression_result$output_dirs

# Check steps already run.
panfam_status(project)

# Create an HTML index of the steps you have run.
report_file <- panfam_report(project)
utils::browseURL(report_file)
```

Creating this index does not execute additional analyses. Module result objects contain `tables`, `run_log`, `output_dirs`, and `module`; `panfam_ssgsea()` returns a named list of sample-score data frames.

## Multiple families and selected cancer types

Set these options during setup, then use the same step functions:

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

Provide either `genes` for a single family or `gene_sets` for named families. Family labels can contain Chinese characters; distinct safe folder names are generated automatically. To change genes, cancer types, or other settings, create a new project with `panfam_setup()`.

For your own expression table, pass `expression = my_expression` to setup. It must have one sample per row, unique sample row names, `Cancer` and `Group` columns, and numeric gene-expression columns named by gene symbol. `Group` must be `Tumor` or `Normal`. Bundled clinical and association data are usable only where sample IDs match.

## Analysis rules and data coverage

### Source plotting styles restored in 0.3.1

Version 0.3.1 restores the plotting styles and principal estimates used by the original scripts while retaining custom gene sets, the modular interface, and bundled-data access:

| Module | Figure and interpretation |
| --- | --- |
| Gene expression and TCGA/GTEx | A six-bin heatmap using `mean(tumor) - mean(normal)` on the supplied expression scale, with numeric labels and the original cut points of -2, -1, 0, 1, and 2. Cells with raw P >= 0.05 are white. The source script calls this value `logFC`; it is not automatically a log2 fold change. |
| TMB and MSI | Red/cyan radar plots of the selected correlation coefficient. Pearson is the default; Spearman can be selected during setup. |
| Immune cells, scores, pathways, and marker panels | Blue/white/red ComplexHeatmap figures; significance stars use raw P values. BH-adjusted values remain available in CSV tables. |
| Family survival | Table-style forest plots for High versus Low median-score groups, and yellow/blue `ggsurvplot` curves with HR, confidence interval, and P annotations. Known day-based follow-up is displayed in years using 365 days per year; unspecified source units are retained. |
| Individual-gene survival | The original Risky/Protective categorical heatmap and count-difference annotation, using the direction of the gene-expression Cox HR and raw P < 0.05. |
| Gene alterations | A count heatmap with stacked alteration-type annotations along rows and columns. Gene/cancer cells count alteration events and can count multiple types in the same sample; the Total row counts unique samples with at least one alteration in the custom input family. Frequency and coverage tables are retained. |
| CNV and methylation | Source-correlation bubble plots with source FDR shown through border coding; no new source FDR is calculated. |
| Immunotherapy | Within-cohort median-score High/Low groups, Yes/No/SD response-composition stacked bars where source categories are available, and accompanying KM curves. |

Family Cox plots primarily report `hr_high_vs_low`, `group_ci_lower`, `group_ci_upper`, `group_p_value`, and `group_fdr`. Groups are fixed as **Low: score <= median; High: score > median**. The continuous model remains available separately through `hr_per_sd`, `ci_lower`, `ci_upper`, `p_value`, and `fdr`, where HR is per score standard deviation. Individual-gene Cox tables also include `hr_per_unit`, `unit_ci_lower`, `unit_ci_upper`, and `unit_p_value` for the original expression unit. These estimates must be interpreted according to the model named by their columns and figure labels.

Newly calculated statistical tables use Benjamini-Hochberg adjustment across valid tests within each result table; group-model and continuous-model P values are adjusted separately. Raw-P plotting conventions are explicitly labelled and do not replace the adjusted values. Source CNV/methylation FDR values are retained. Insufficient samples, missing outcomes, and unavailable data are recorded in results and run status.

The restored plots do not imply that every number must match every inherited script output. Samples are never manually moved between median groups, and High/Low labels are not reversed. The ROC step retains an empirical, unsmoothed curve with a fixed direction: higher scores predict tumor. It therefore differs from the original script's smoothed ROC, and AUC below 0.5 is retained. This is a comparison within the bundled TCGA data and is not an external validation.

Mutation, CNV, and methylation inputs cover the source project's gene panels rather than every human gene. Custom families are searched by gene symbol. Unsupported genes are reported explicitly, and unavailable results are marked `skipped` or `partial`; these statuses do not mean a negative biological finding. TCGA/GTEx comparisons pool source-labelled normal samples and do not correct study or batch effects.

Immunotherapy uses baseline samples where annotated, retains and records samples with an unannotated treatment time, excludes known non-immunotherapy controls, and retains one sample per patient. Response-composition figures preserve available source categories; they are distinct from the separate responder/non-responder comparison tables. GSE135222 uses its native expression TSV and bundled GEO metadata with **PFS**; the other cohorts use OS. The GSE135222 PFS unit is not specified in GEO, so source units are retained. Unknown outcomes remain missing. Non-finite expression genes are excluded and recorded.

## Data access and optional batch interface

```r
panfam_data_files()
panfam_load_data("all_tpm")
panfam_data_path("clinical/ALL_OS.Rdata")
```

Other dataset keys include `tcga_gtex`, `paired_tpm`, `ALL_OS`, `DSS_PFI_DFI`, `All_TMB`, `All_MSI`, `ssgsea_immucells`, `immuscore`, `TME_Pathaway_GSVA_Score`, and `HALLMARK_GSVA_Score`. Change the reconstructed-data cache with `options(PanFamR.cache_dir = "path/to/cache")` before loading data.

For users who prefer a batch run, `panfam_analyze()` remains available. It is optional; the modular workflow above does not require it. `panfam_score()` also remains available for direct scoring without a project.

## Data provenance and license

Maintainer: Serena (<g512827252@outlook.com>). Report problems through [GitHub Issues](https://github.com/JianGao0852/PanFamR/issues).

The [dataset provenance and reference record](inst/DATA_SOURCES.zh-CN.md) maps TCGA, GTEx, GEO, TIGER, MSigDB Hallmark, and GENCODE to the bundled inputs and provides cohort papers and 16 references. Some download portals, versions, algorithms, and precomputed-score generation records remain unconfirmed. Patient-level TMB/MSI and CNV/methylation tables are not attributed to MSigDB. The MIT license covers package code and does not replace dataset-specific terms; adding citations does not confirm all redistribution conditions. MSigDB GMT collections are not bundled.
