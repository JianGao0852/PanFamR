\name{panfam_score}
\alias{panfam_score}
\title{Score Arbitrary Gene Families}
\usage{panfam_score(gene_sets = NULL, expression = NULL, genes = NULL,
  family_name = "CustomFamily", min_size = 1L, output_dir = NULL)}
\arguments{
\item{gene_sets}{A named list of gene vectors, or a single gene vector.}
\item{expression}{Optional sample-by-gene data frame with Cancer and Group columns.}
\item{genes}{A single gene vector, as an alternative to gene_sets.}
\item{family_name}{Label for a single gene vector.}
\item{min_size}{Minimum matched genes; defaults to one.}
\item{output_dir}{Optional score-output directory.}
}
\value{A named list of score data frames, with a gene_coverage attribute.}
\description{Calculates ssGSEA for one or more custom gene families, sharing the ranking calculation across families and normalizing each family by its own score range.}
