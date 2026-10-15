\name{panfam_ssgsea}
\alias{panfam_ssgsea}
\title{Compute and Save Project ssGSEA Scores}
\usage{panfam_ssgsea(project)}
\arguments{\item{project}{Project returned by panfam_setup().}}
\value{Named list of sample-level score data frames with Score, Cancer and Group columns. Attributes include gene_coverage and background_notes.}
\description{Loads expression data, scores each configured gene family, and writes scores and coverage under each family's ssgsea directory.}
\details{Scores use the original alpha of 0.25 and each family's own score range for normalization. Ranking is shared across families. The full expression background is used. Repeated calls reuse the scores without creating another score run. Score-based analysis functions call this function automatically if scores have not yet been computed.}
\seealso{\code{\link{panfam_setup}}, \code{\link{panfam_score}}, \code{\link{panfam_steps}}}
