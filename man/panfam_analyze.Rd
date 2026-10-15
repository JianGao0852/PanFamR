\name{panfam_analyze}
\alias{panfam_analyze}
\title{Run an Automated Pan-Cancer Gene-Family Workflow}
\usage{
panfam_analyze(gene_sets = NULL, expression = NULL,
  output_dir = "PanFamR_results", genes = NULL,
  family_name = "CustomFamily",
  analyses = c("expression", "tcga_gtex", "paired", "survival", "tmb",
    "msi", "immune", "immunotherapy", "alterations", "cnv_methylation"),
  cancers = NULL, min_size = 1L,
  correlation_method = c("pearson", "spearman"))
}
\arguments{
\item{gene_sets}{A named list of gene vectors, or a single gene vector.}
\item{expression}{Optional sample-by-gene data frame with Cancer and Group columns. Uses the bundled TCGA expression table by default.}
\item{output_dir}{Result directory. An existing report causes a timestamped sibling directory to be used. NULL returns scores only.}
\item{genes}{A single gene vector, as an alternative to gene_sets.}
\item{family_name}{Label for a single gene vector.}
\item{analyses}{Modules to run. Defaults to all listed modules; character(0) runs scoring only.}
\item{cancers}{Optional cancer labels to select.}
\item{min_size}{Minimum matched genes for ssGSEA; positive integer, default one.}
\item{correlation_method}{Pearson by default, matching the source script; Spearman is also available.}
}
\value{A panfam_result containing scores, tables, gene coverage, run status, settings, and output directory. When output_dir is NULL, returns a named list of score data frames.}
\description{Scores arbitrary gene families with the bundled data and automatically writes analysis tables, plots, and an HTML report.}
\details{Family survival figures use fixed median-group High versus Low Cox estimates. Continuous per-standard-deviation estimates and their P values remain in separate columns; gene-level tables also contain per-expression-unit estimates. Source-style figure marks use raw P values as labelled, while newly computed result tables retain Benjamini-Hochberg adjustment. Pearson correlation is the source-script default. Modules with insufficient data are recorded in the run log. Source mutation and CNV/methylation inputs cover only the genes present in the original project.}
\examples{
\dontrun{
result <- panfam_analyze(genes = c("COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1"),
                        family_name = "MyFamily")
}
}
