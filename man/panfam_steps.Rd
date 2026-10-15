\name{panfam_steps}
\alias{panfam_steps}
\alias{panfam_expression}
\alias{panfam_tcga_gtex}
\alias{panfam_paired}
\alias{panfam_survival}
\alias{panfam_gene_survival}
\alias{panfam_tumor_scores}
\alias{panfam_roc}
\alias{panfam_tmb}
\alias{panfam_msi}
\alias{panfam_immune_cells}
\alias{panfam_immune_scores}
\alias{panfam_tme}
\alias{panfam_hallmark}
\alias{panfam_checkpoints}
\alias{panfam_chemokines}
\alias{panfam_chemokine_receptors}
\alias{panfam_immune_stimulators}
\alias{panfam_immune_inhibitors}
\alias{panfam_immunotherapy}
\alias{panfam_alterations}
\alias{panfam_cnv}
\alias{panfam_methylation}
\alias{print.panfam_step}
\title{Run Individual Pan-Cancer Analysis Steps}
\usage{
panfam_expression(project)
panfam_tcga_gtex(project)
panfam_paired(project)
panfam_survival(project, endpoints = c("OS", "DSS", "PFI", "DFI"))
panfam_gene_survival(project, endpoints = "OS")
panfam_tumor_scores(project)
panfam_roc(project)
panfam_tmb(project)
panfam_msi(project)
panfam_immune_cells(project)
panfam_immune_scores(project)
panfam_tme(project)
panfam_hallmark(project)
panfam_checkpoints(project)
panfam_chemokines(project)
panfam_chemokine_receptors(project)
panfam_immune_stimulators(project)
panfam_immune_inhibitors(project)
panfam_immunotherapy(project)
panfam_alterations(project)
panfam_cnv(project)
panfam_methylation(project)
}
\arguments{
\item{project}{Project returned by panfam_setup(). Paths, gene families and settings need not be repeated.}
\item{endpoints}{Non-empty subset of OS, DSS, PFI and DFI. Only selected endpoints are analyzed.}
}
\value{A panfam_step list with module, tables, run_log and output_dirs. Access a table with result$tables[[family_name]][[table_name]]. CSV tables and result.rds are saved under output_dir/family/module; supported figures are under its plots directory. Repeating a step creates a numbered directory, preserving its previous results.}
\description{Each function runs its named analysis only. Inputs are bundled and loaded automatically; TCGA scores are computed on demand and reused within the project. Users can execute any desired steps individually.}
\details{
\describe{
\item{panfam_expression}{Tumor-normal differences in TCGA family scores and individual gene expression, with boxplots and an expression effect heatmap.}
\item{panfam_tcga_gtex}{Individual-gene differences pooling TCGA/GTEx source-labelled normal samples. This comparison does not correct study/batch differences.}
\item{panfam_paired}{TCGA patient-matched tumor-normal score differences and paired plots.}
\item{panfam_survival}{Univariate Cox models, fixed-median KM/log-rank analysis and survival plots for selected endpoints.}
\item{panfam_gene_survival}{Separate gene-expression Cox and median log-rank summaries for selected endpoints, with an association heatmap. Defaults to OS.}
\item{panfam_tumor_scores}{Family-score summaries and distributions in tumor samples, retaining one tumor per patient.}
\item{panfam_roc}{Exploratory tumor-versus-normal score ROC curves in the same TCGA dataset. The direction is fixed: higher score predicts tumor. AUC below 0.5 is retained. This is an in-sample description, not external predictive validation.}
\item{panfam_tmb, panfam_msi}{Family-score associations with TMB or MSI.}
\item{panfam_immune_cells, panfam_immune_scores}{Associations with precomputed immune-cell or stromal/immune/ESTIMATE scores.}
\item{panfam_tme, panfam_hallmark}{Associations with precomputed TME or Hallmark pathway scores.}
\item{panfam_checkpoints, panfam_chemokines, panfam_chemokine_receptors, panfam_immune_stimulators, panfam_immune_inhibitors}{Separate associations with each marker-expression panel, with gene coverage recorded.}
\item{panfam_immunotherapy}{Within-cohort family scoring, response comparisons and available survival in the four bundled immunotherapy cohorts.}
\item{panfam_alterations}{Available mutation/fusion/amplification/deletion frequencies from bundled gene-panel tables.}
\item{panfam_cnv, panfam_methylation}{Separate source CNV-expression or methylation-expression correlations and source FDR, with coverage tables.}
}
TCGA correlation and survival analyses use one tumor sample per patient and unambiguous ID matching. Wilcoxon differences require at least three samples per group or three pairs; correlations require at least five matched samples and variable inputs. Cox models require at least ten complete samples and five events; median groups require three samples each. Median groups are never manually reassigned. Newly computed tests use Benjamini-Hochberg adjustment across valid tests in each returned result table for that family and step. Source CNV/methylation FDR is retained.

All functions support arbitrary family labels and gene vectors. Bundled alteration/CNV/methylation inputs have limited gene coverage; unavailable inputs produce explicit coverage/status records rather than invented results. Errors for individual families are saved in error.txt and run_log while remaining families can proceed. Immunotherapy scoring uses its own cohort data, not TCGA scores. Individual-gene survival, TCGA/GTEx, alterations, CNV and methylation do not compute TCGA scores.
}
\seealso{\code{\link{panfam_setup}}, \code{\link{panfam_ssgsea}}, \code{\link{panfam_status}}, \code{\link{panfam_report}}}
