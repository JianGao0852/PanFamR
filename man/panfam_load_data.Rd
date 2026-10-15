\name{panfam_load_data}
\alias{panfam_load_data}
\title{Load a Bundled PanFamR Dataset}
\usage{panfam_load_data(name = "all_tpm")}
\arguments{\item{name}{Dataset key, such as \code{all_tpm}, \code{tcga_gtex}, \code{ALL_OS}, or \code{All_TMB}.}}
\value{A named list containing the objects stored in the selected RData file.}
\description{Loads one of the data files bundled with PanFamR. Large split expression files are joined and cached on first use.}
