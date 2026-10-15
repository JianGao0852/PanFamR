#' PanFamR pan-cancer analysis functions
#'
#' Functions for ssGSEA scoring and downstream pan-cancer plots. Score inputs
#' should contain one sample per row, `Cancer` and `Group` metadata columns,
#' and gene symbols as the remaining column names. Functions that need
#' precomputed clinical or analysis tables use bundled files by default;
#' [panfam_set_paths()] can point to compatible user-supplied replacements.
#'
#' @param tpm Expression data frame used for ssGSEA and expression-based
#'   analyses.
#' @param cox_genes Character vector of genes to score.
#' @param df Sample-level ssGSEA score data frame, typically returned by
#'   [ssGsea_Score()].
#' @param rt Expression or survival input table used by the corresponding
#'   analysis function.
#' @param cancer Cancer type abbreviation.
#' @param Target Survival endpoint such as `OS`, `DSS`, `PFI`, or `DFI`.
#' @param folder_path Input directory for the corresponding analysis files.
#' @param file_path Input file path.
#' @param analysis_target Either `CNV` or `Methylation`.
#' @param dir_name Output name prefix for the TCGA-GTEx comparison.
#' @param tcga_gtex_mrna_pheno TCGA/GTEx expression and phenotype table.
#' @param output_file Optional `.Rdata` path in which to save the score object.
#' @param family_name Optional name for the gene set.
#' @param gene_sets Named list of gene vectors.
#' @param output_dir Directory for score files.
#' @param ... Additional arguments reserved for future extensions.
#'
#' @name PanFamR
#' @aliases PanFamR Alterations CNV_Methylation_mRNA CNV_Methylation_mRNA1
#'   plot_cox_analysis_forest plot_cox_analysis_forest_picture plot_cox_heatmap1
#'   plot_cox_heatmap2 plot_diff_heatmap plot_diff_heatmap1 plot_km plot_km1
#'   plot_km2 plot_OS_KM_expression run_ssGSEA ssGsea_Score ssGsea_Score_diff
#'   ssGsea_Score_diff1 ssGsea_Score_diff_Tumor ssGsea_Score_KM
#'   ssGsea_Score_one_paired_boxplot ssGsea_Score_paired_boxplot
#'   ssGsea_Score_pan_cancer_chemokine ssGsea_Score_pan_cancer_chemokine_receptor
#'   ssGsea_Score_pan_cancer_HALLMARK_GSVA ssGsea_Score_pan_cancer_ICG
#'   ssGsea_Score_pan_cancer_immune_inhibitor
#'   ssGsea_Score_pan_cancer_immunecellratio ssGsea_Score_pan_cancer_immunescore
#'   ssGsea_Score_pan_cancer_immune_stimulator ssGsea_Score_pan_cancer_MSI
#'   ssGsea_Score_pan_cancer_TMB ssGsea_Score_pan_cancer_TME_GSVA
#'   ssGsea_Score_pan_cancer_fores ssGsea_Score_ROC ssGsea_Score_Tumor_Normal
#'   panfam_set_paths
#' @details The package includes TCGA/GTEx expression, survival, mutation,
#'   immune, and supporting analysis tables. Review the data notes in
#'   `inst/extdata/README-data.md` before redistributing the package. The
#'   legacy plotting functions write figures to the current working directory.
#' @return `ssGsea_Score()` returns a data frame with `Score`, `Cancer`, and
#'   `Group` columns and sample IDs as row names. `run_ssGSEA()` returns a
#'   named list of such data frames. Plotting functions generally save a PDF
#'   and return `NULL` or a plot object, depending on the function.
#' @seealso [panfam_set_paths()], [run_ssGSEA()]
#' @keywords package
NULL

utils::globalVariables(c(
  "Cancer", "CancerType", "DFI", "DFI.time", "DSS", "DSS.time", "Gene",
  "GeneType", "Group", "ID", "MSI", "OS", "OS.time", "PFI", "PFI.time",
  "Score", "Status", "TMB", "Total", "Type", "altered", "cancer",
  "cancertype", "event", "fdr", "fill_bin", "freq", "gene", "group",
  "logFC", "n", "p.value", "pval", "response", "sample", "spm", "symbol",
  "time", "value", "variable", "x", "y", "Feature", "rank_biserial",
  "mark", "hr_per_sd", "ci_lower", "ci_upper", "Endpoint", "Response",
  "Association", "Patient", "correlation", ".", "Altered", "DSS_PFI_DFI", "HR",
  "HR (95% CI)", "HR_gt1", "HR_lower", "HR_lt1", "HR_upper", "HR_val",
  "P", "P value", "P_value", "RiskyScore", "TME_Pathaway", "Total(N)",
  "immucells", "immuscore", "logFC_text", "merged_data", "neg_log10_fdr",
  "paired", "paired_tpm", "percent", "point_size", "sig", "status",
  "total_neg_log10_fdr"
))

