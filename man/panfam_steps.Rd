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
\value{A panfam_step list with module, tables, run_log and output_dirs. Access a table with result$tables[[family_name]][[table_name]]. CSV tables and result.rds are saved under output_dir/family/module; supported PDF figures and 300 dpi white-background PNG figures are under its plots directory. Repeating a step creates a numbered directory, preserving its previous results.}
\description{Each function runs its named analysis only. Inputs are bundled and loaded automatically; TCGA scores are computed on demand and reused within the project. Users can execute any desired steps individually.}
\details{
\describe{
\item{panfam_expression}{Tumor-normal differences in TCGA family scores and individual gene expression, with pan-cancer boxplots, separate score_tumor_normal_<Cancer> plots and a six-bin expression heatmap. Single-cancer plots reuse the corresponding score_differences statistics. Numeric cell labels and colors use mean(tumor) minus mean(normal) on the supplied expression scale, with cut points of -2, -1, 0, 1 and 2; raw P >= 0.05 is shown in white. The source-script name logFC does not establish that the value is log2 fold change.}
\item{panfam_tcga_gtex}{Individual-gene differences pooling TCGA/GTEx source-labelled normal samples, with the same mean-difference heatmap. This comparison does not correct study/batch differences.}
\item{panfam_paired}{TCGA patient-matched tumor-normal score differences and paired plots.}
\item{panfam_survival}{Genuine fixed-median High-versus-Low Cox models, KM/log-rank analysis and table-style forest plots for selected endpoints. Primary plot estimates use hr_high_vs_low, group_ci_lower, group_ci_upper, group_p_value and group_fdr. Additional continuous score-SD estimates remain in hr_per_sd, ci_lower, ci_upper, p_value and fdr. Yellow/blue ggsurvplot curves show the real groups and HR/CI/P annotations. Each KM main figure is accompanied by a file ending _with_risk_table, using the same fitted curve and groups.}
\item{panfam_gene_survival}{Separate gene-expression Cox and median log-rank summaries for selected endpoints. Original expression-unit estimates are retained in hr_per_unit, unit_ci_lower, unit_ci_upper and unit_p_value, alongside score-SD estimates. The source Risky/Protective heatmap classifies HR direction using raw P < 0.05, with a count-difference annotation. Defaults to OS.}
\item{panfam_tumor_scores}{Family-score summaries and distributions in tumor samples, retaining one tumor per patient.}
\item{panfam_roc}{Exploratory tumor-versus-normal score ROC curves in the same TCGA dataset. The direction is fixed: higher score predicts tumor. AUC below 0.5 is retained. Curves are empirical and unsmoothed, differing from the inherited script's smoothed ROC. This is an in-sample description, not external predictive validation.}
\item{panfam_tmb, panfam_msi}{Family-score associations with TMB or MSI, presented as red/cyan radar plots.}
\item{panfam_immune_cells, panfam_immune_scores}{Associations with precomputed immune-cell or stromal/immune/ESTIMATE scores, with blue/white/red ComplexHeatmap figures and raw-P significance stars.}
\item{panfam_tme, panfam_hallmark}{Associations with precomputed TME or Hallmark pathway scores, with the same raw-P heatmap convention.}
\item{panfam_checkpoints, panfam_chemokines, panfam_chemokine_receptors, panfam_immune_stimulators, panfam_immune_inhibitors}{Separate associations with each marker-expression panel, with gene coverage recorded and raw-P heatmap significance stars.}
\item{panfam_immunotherapy}{Within-cohort family scoring, response comparisons and available survival in the four bundled immunotherapy cohorts. Real median-score High/Low groups are used for Yes/No/SD response-composition stacked bars where source categories are available, with accompanying KM curves and companion risk-table figures. Composition categories are distinct from the separate responder/non-responder comparison tables.}
\item{panfam_alterations}{Available mutation/fusion/amplification/deletion frequencies and coverage from bundled gene-panel tables, with a count heatmap and stacked alteration-type annotations along rows and columns. Gene/cancer cells count alteration events and may count multiple types in the same sample. The Total row counts unique samples with at least one alteration in the custom input family; individual cells should not be interpreted as unique-patient counts.}
\item{panfam_cnv, panfam_methylation}{Separate source CNV-expression or methylation-expression correlations and source FDR, with coverage tables and source-correlation bubble plots using FDR border coding.}
}
TCGA correlation and survival analyses use one tumor sample per patient and unambiguous ID matching. Pearson is the project default, matching the source script; Spearman can be selected with correlation_method in panfam_setup(). Wilcoxon differences require at least three samples per group or three pairs; correlations require at least five matched samples and variable inputs. Cox models require at least ten complete samples and five events; median groups require three samples each. Groups are Low for score <= median and High for score > median. Groups are never manually reassigned, and High/Low labels are not reversed. Known day-based follow-up is displayed in years using 365 days per year, while unspecified source units are retained.

Newly computed tests use Benjamini-Hochberg adjustment across valid tests in each returned result table for that family and step. Group-model and continuous-model P values are adjusted separately. Raw-P plotting conventions are labelled explicitly and do not replace the adjusted values in the tables. Source CNV/methylation FDR is retained. Restoring the source plotting styles does not imply equality with every inherited numeric result: manual sample changes and label swaps are excluded, and ROC direction and smoothing differ as stated above.

All functions support arbitrary family labels and gene vectors. Bundled alteration/CNV/methylation inputs have limited gene coverage; unavailable inputs produce explicit coverage/status records rather than invented results. Errors for individual families are saved in error.txt and run_log while remaining families can proceed. Immunotherapy scoring uses its own cohort data, not TCGA scores. Individual-gene survival, TCGA/GTEx, alterations, CNV and methylation do not compute TCGA scores.

Input sample counts are available through the separate \code{\link{panfam_sample_summary}} interface. This project-level step does not depend on gene matching or scores and produces count and metadata tables, without figures or repeated summaries for each family.

Cancer selection in \code{panfam_setup()} is applied before scoring and can change the score-normalization range. Interpret scores within the cancer selection used for that project.
}
\seealso{\code{\link{panfam_setup}}, \code{\link{panfam_sample_summary}}, \code{\link{panfam_ssgsea}}, \code{\link{panfam_status}}, \code{\link{panfam_report}}}
