\name{panfam_set_paths}
\alias{panfam_set_paths}
\title{Configure Analysis Data Locations}
\usage{panfam_set_paths(analysis_dir = NULL, clinical_dir = NULL)}
\arguments{
\item{analysis_dir}{Optional directory containing compatible auxiliary analysis tables.}
\item{clinical_dir}{Optional directory containing compatible survival tables.}
}
\value{The configured paths, invisibly.}
\description{Overrides the bundled analysis and clinical data directories used by legacy plotting functions.}
