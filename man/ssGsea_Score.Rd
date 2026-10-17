\name{ssGsea_Score}
\alias{ssGsea_Score}
\title{Score a Gene Set with ssGSEA}
\usage{ssGsea_Score(tpm, cox_genes, output_file = NULL, family_name = NULL)}
\arguments{
\item{tpm}{A data frame with one sample per row, \code{Cancer} and \code{Group} metadata columns, and gene symbols as the remaining column names.}
\item{cox_genes}{Character vector of gene symbols. At least five must occur in the expression data.}
\item{output_file}{Optional path for saving the result as an RData file.}
\item{family_name}{Optional label for the gene set.}
}
\value{A data frame with sample-level \code{Score}, \code{Cancer}, and \code{Group} columns.}
\description{Calculates ssGSEA scores for one user-defined gene set.}
\details{The default \code{output_file = NULL} returns scores without automatically writing an RData file. Supply \code{output_file} explicitly to save one, or use \code{\link{panfam_ssgsea}} for automatic scoring and output within a configured project.}
