\name{panfam_setup}
\alias{panfam_setup}
\alias{print.panfam_project}
\title{Configure a Modular Gene-Family Analysis Project}
\usage{
panfam_setup(work_dir = getwd(), output_dir = "PanFamR_results",
  genes = NULL, family_name = "CustomFamily", gene_sets = NULL,
  cancers = NULL, min_size = 1L,
  correlation_method = c("pearson", "spearman"), expression = NULL)
}
\arguments{
\item{work_dir}{Working directory used to resolve a relative output directory. Created if necessary. Does not change the R working directory.}
\item{output_dir}{Result directory. An existing project/report/log causes a timestamped sibling directory to be created to preserve previous results.}
\item{genes}{Character vector of gene symbols for one user-defined family. Supply either genes or gene_sets.}
\item{family_name}{Display label for genes.}
\item{gene_sets}{Alternative named list of gene vectors for multiple user-defined families.}
\item{cancers}{Optional cancer labels to select, such as BRCA or LUAD. NULL uses all cancer types.}
\item{min_size}{Minimum matched gene count for ssGSEA; a positive integer, default one.}
\item{correlation_method}{Pearson (source-script default) or Spearman correlation.}
\item{expression}{Optional sample-by-gene data frame with Cancer and Group columns; Group must be Tumor or Normal. NULL uses bundled expression.}
}
\value{A panfam_project with fixed configuration and an internal cache. Pass this object to each analysis function. The absolute result directory is project$output_dir.}
\description{Configures paths and arbitrary gene families once. This function does not run analyses or score genes. Subsequent functions load their required bundled inputs automatically.}
\details{Expression and scores are cached within this project. Configuration bindings are locked to avoid stale scores after changing genes or settings. Create another project when settings change. The saved project_config.rds contains lightweight settings only and is a record, not a resume checkpoint. Do not serialize the project itself: its cache can contain the large expression table. The project remains valid if getwd() changes later.}
\seealso{\code{\link{panfam_ssgsea}}, \code{\link{panfam_steps}}, \code{\link{panfam_status}}, \code{\link{panfam_report}}}
\examples{
\dontrun{
project <- panfam_setup(work_dir = "D:/research",
  output_dir = "results", genes = c("TP53", "MDM2", "CDKN1A"),
  family_name = "MyFamily")
scores <- panfam_ssgsea(project)
expression_result <- panfam_expression(project)
survival_result <- panfam_survival(project, endpoints = "OS")
}
}
