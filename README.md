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

The installation command above can be used after the GitHub repository is published. Bundled inputs occupy approximately 1.29 GiB. Allow about 6 GB of free disk space for download, installation, and cached data files. Large inputs, including the TCGA/GTEx expression matrices, paired expression data, and Braun immunotherapy expression RDS, are distributed in numbered parts of at most 20 MiB. All split inputs are reconstructed automatically in the local cache when first needed. Once installed, analyses use local data without downloading datasets. The full-matrix workflow has been run on a 32 GB RAM Windows computer; selecting fewer cancer types reduces calculation size.

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
  family_name = "MyFamily"
)
```

Replace the paths and gene symbols with your own. An absolute result path is also accepted. Relative result paths are resolved against `work_dir`. Setup records your settings and prepares the output directory; it does not run analyses or change the working directory. Keep using the same `project` object for subsequent steps.

Genes do not need a shared prefix or a predefined family name. Matching is exact and case sensitive. Missing symbols are recorded in coverage tables; a wholly unmatched family produces an error. The default `min_size = 1` permits small families; increase it to require more matched genes.

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

Step folders contain the available CSV tables and PDF/PNG plots. Repeating a step creates a new numbered folder, preserving its earlier output. Starting a new project in an already configured result directory uses a timestamped sibling directory. The actual root is available as `project$output_dir`.

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

Survival uses a fixed median split and univariate Cox models per score SD. Newly calculated statistical tables use Benjamini-Hochberg adjustment; source CNV/methylation FDR values are retained. Insufficient samples, missing outcomes, and unavailable data are recorded in results and run status.

The ROC step compares tumor and normal samples within the bundled TCGA data. Higher scores are always treated as predicting tumor, so AUC can be below 0.5; this comparison is exploratory and is not an external validation.

Mutation, CNV, and methylation inputs cover the source project's gene panels rather than every human gene. Custom families are searched by gene symbol. Unsupported genes are reported explicitly, and unavailable results are marked `skipped` or `partial`; these statuses do not mean a negative biological finding. TCGA/GTEx comparisons pool source-labelled normal samples and do not correct study or batch effects.

Immunotherapy uses baseline samples where annotated, excludes known non-immunotherapy controls, and retains one sample per patient. GSE135222 uses its native expression TSV and bundled GEO metadata with **PFS**; the other cohorts use OS. The GSE135222 PFS unit is not specified in GEO, so source units are retained. Unknown outcomes remain missing. Non-finite expression genes are excluded and recorded.

## Data access and optional batch interface

```r
panfam_data_files()
panfam_load_data("all_tpm")
panfam_data_path("clinical/ALL_OS.Rdata")
```

Other dataset keys include `tcga_gtex`, `paired_tpm`, `ALL_OS`, `DSS_PFI_DFI`, `All_TMB`, `All_MSI`, `ssgsea_immucells`, `immuscore`, `TME_Pathaway_GSVA_Score`, and `HALLMARK_GSVA_Score`. Change the reconstructed-data cache with `options(PanFamR.cache_dir = "path/to/cache")` before loading data.

For users who prefer a batch run, `panfam_analyze()` remains available. It is optional; the modular workflow above does not require it. `panfam_score()` also remains available for direct scoring without a project.

## Before publishing

Maintainer: Serena (<g512827252@outlook.com>). Report problems through [GitHub Issues](https://github.com/JianGao0852/PanFamR/issues). Record source URLs, citations, and redistribution terms for the author-supplied datasets in `inst/extdata/README-data.md`. The MIT license covers package code and does not replace dataset-specific terms. MSigDB GMT collections are not bundled; custom gene vectors work directly, while bundled pathway-score tables require their own provenance and terms to be documented.
