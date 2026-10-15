# PanFamR: run the setup block once, then run individual analysis steps.
# In RStudio, select a line or block and click Run (Ctrl+Enter).
# Edit the working directory, output directory, and gene family below.

library(PanFamR)

# ---- 1. Working directory and result directory ----
work_dir <- "D:/PanFamR_analysis"
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

# ---- 2. Define your family once ----
project <- panfam_setup(
  work_dir = getwd(),
  output_dir = "results",
  genes = c("COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1"),
  family_name = "MyFamily",
  cancers = NULL # Replace NULL with c("BRCA", "LUAD") to select cancer types.
)

# ---- 3. Calculate ssGSEA scores ----
scores <- panfam_ssgsea(project)

# ---- 4. Expression and score differences ----
expression_result <- panfam_expression(project)

# ---- 5. Choose further steps ----
# Remove the leading # from a step you need, select that line, and click Run.
# Use the same project object for every step.

# tcga_gtex_result <- panfam_tcga_gtex(project)
# paired_result <- panfam_paired(project)
# tumor_scores_result <- panfam_tumor_scores(project)
# roc_result <- panfam_roc(project)
# survival_result <- panfam_survival(project, endpoints = "OS")
# survival_result <- panfam_survival(project, endpoints = c("OS", "DSS", "PFI", "DFI"))
# gene_survival_result <- panfam_gene_survival(project, endpoints = "OS")
# tmb_result <- panfam_tmb(project)
# msi_result <- panfam_msi(project)
# immune_cells_result <- panfam_immune_cells(project)
# immune_scores_result <- panfam_immune_scores(project)
# tme_result <- panfam_tme(project)
# hallmark_result <- panfam_hallmark(project)
# checkpoints_result <- panfam_checkpoints(project)
# chemokines_result <- panfam_chemokines(project)
# chemokine_receptors_result <- panfam_chemokine_receptors(project)
# immune_stimulators_result <- panfam_immune_stimulators(project)
# immune_inhibitors_result <- panfam_immune_inhibitors(project)
# immunotherapy_result <- panfam_immunotherapy(project)
# alterations_result <- panfam_alterations(project)
# cnv_result <- panfam_cnv(project)
# methylation_result <- panfam_methylation(project)

# ---- 6. Inspect a result or summarize completed steps ----
# expression_result$tables[["MyFamily"]][["score_differences"]]
# expression_result$run_log
# expression_result$output_dirs
# panfam_status(project)
# report_file <- panfam_report(project)
# utils::browseURL(report_file)
