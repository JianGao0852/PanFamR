\name{run_ssGSEA}
\alias{run_ssGSEA}
\title{Score Custom Gene Sets with ssGSEA}
\usage{run_ssGSEA(tpm, gene_sets, output_dir = ".")}
\arguments{
\item{tpm}{A data frame with one sample per row, \code{Cancer} and \code{Group} metadata columns, and gene symbols as the remaining column names.}
\item{gene_sets}{A non-empty named list of character vectors. Each list element is a user-defined gene family or gene set.}
\item{output_dir}{Directory in which one score file per gene set is written.}
}
\value{A named list of sample-level score data frames.}
\description{Calculates ssGSEA scores for multiple user-defined gene sets.}
