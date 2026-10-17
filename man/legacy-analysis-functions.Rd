\name{legacy-analysis-functions}
\alias{legacy-analysis-functions}
\alias{Alterations}
\alias{CNV_Methylation_mRNA}
\alias{CNV_Methylation_mRNA1}
\alias{plot_cox_analysis_forest}
\alias{plot_cox_analysis_forest_picture}
\alias{plot_cox_heatmap1}
\alias{plot_cox_heatmap2}
\alias{plot_diff_heatmap}
\alias{plot_diff_heatmap1}
\alias{plot_km}
\alias{plot_km1}
\alias{plot_km2}
\alias{plot_OS_KM_expression}
\alias{ssGsea_Score_diff}
\alias{ssGsea_Score_diff1}
\alias{ssGsea_Score_diff_Tumor}
\alias{ssGsea_Score_KM}
\alias{ssGsea_Score_one_paired_boxplot}
\alias{ssGsea_Score_paired_boxplot}
\alias{ssGsea_Score_pan_cancer_chemokine}
\alias{ssGsea_Score_pan_cancer_chemokine_receptor}
\alias{ssGsea_Score_pan_cancer_HALLMARK_GSVA}
\alias{ssGsea_Score_pan_cancer_ICG}
\alias{ssGsea_Score_pan_cancer_immune_inhibitor}
\alias{ssGsea_Score_pan_cancer_immunecellratio}
\alias{ssGsea_Score_pan_cancer_immunescore}
\alias{ssGsea_Score_pan_cancer_immune_stimulator}
\alias{ssGsea_Score_pan_cancer_MSI}
\alias{ssGsea_Score_pan_cancer_TMB}
\alias{ssGsea_Score_pan_cancer_TME_GSVA}
\alias{ssGsea_Score_pan_cancer_fores}
\alias{ssGsea_Score_ROC}
\alias{ssGsea_Score_Tumor_Normal}
\title{Legacy Pan-Cancer Analysis Functions}
\usage{
ssGsea_Score_Tumor_Normal(df, cancer, tpm = NULL)
ssGsea_Score_ROC(df, cancer, tpm = NULL)
}
\arguments{
\item{df}{A sample-level score data frame containing \code{Score}, \code{Cancer} and \code{Group}, with sample IDs as row names.}
\item{cancer}{Cancer label to analyze.}
\item{tpm}{Optional expression metadata containing \code{Cancer} and \code{Group}. If \code{NULL}, labels are taken directly from \code{df}, preserving two-argument calls without loading an expression matrix. If supplied, labels are matched to score rows by exact sample IDs.}
}
\description{Additional analysis and plotting functions adapted from the original PanFamR scripts. These functions cover expression differences, survival, mutation, CNV, methylation, TMB/MSI, immune association, and immunotherapy cohorts.}
\details{Functions generally accept an expression or score table and may write figures or intermediate files to the current working directory. Inputs not bundled for a particular analysis must be supplied by the user in the format expected by the original workflow. The modular workflow starts with \code{\link{panfam_setup}} and the functions documented in \code{\link{panfam_steps}}; \code{panfam_analyze()} is an optional batch interface. See the README and \code{panfam_data_files()} for bundled data.

\code{ssGsea_Score_Tumor_Normal()} compares Normal versus Tumor scores within one cancer; it does not perform stage or grade analysis. The legacy ROC retains its smoothed curve and automatic direction, while \code{panfam_roc()} uses a fixed higher-score-predicts-tumor direction and an empirical unsmoothed curve. Retained names do not imply identical parameter requirements, filtering or statistical output between all legacy functions and the modular workflow.}
\value{Most plotting functions save a PDF and return a plot or \code{NULL}; consult the function source for the exact return value and arguments.}
