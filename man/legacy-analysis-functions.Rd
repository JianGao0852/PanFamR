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
\description{Additional analysis and plotting functions adapted from the original PanFamR scripts. These functions cover expression differences, survival, mutation, CNV, methylation, TMB/MSI, immune association, and immunotherapy cohorts.}
\details{Functions generally accept an expression or score table and may write figures or intermediate files to the current working directory. Inputs not bundled for a particular analysis must be supplied by the user in the format expected by the original workflow. The generic one-call workflow is \code{panfam_analyze()}. See the README and \code{panfam_data_files()} for bundled data.}
\value{Most plotting functions save a PDF and return a plot or \code{NULL}; consult the function source for the exact return value and arguments.}
