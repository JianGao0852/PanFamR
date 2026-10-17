\name{panfam_sample_summary}
\alias{panfam_sample_summary}
\alias{get_all_tpmsample}
\title{Summarize Input Sample Counts by Cancer and Group}
\usage{
panfam_sample_summary(project, data_source = c("tcga", "tcga_gtex"))

get_all_tpmsample(tpm)
}
\arguments{
\item{project}{Project returned by \code{\link{panfam_setup}}. The configured cancer selection is applied; gene families and ssGSEA scores are not used for counting.}
\item{data_source}{Input to summarize. \code{"tcga"}, the default, uses the project's custom expression table if supplied, otherwise bundled TCGA \code{all_tpm}. \code{"tcga_gtex"} always uses the bundled TCGA/GTEx integrated data, independently of the project's custom expression table.}
\item{tpm}{A data frame containing \code{Cancer} and \code{Group} columns. Extra columns are ignored. Groups must be exactly \code{"Normal"} or \code{"Tumor"}; missing or blank cancer/group labels and unsupported groups are errors. The cancer label \code{"Total"} is reserved for the summary row.}
}
\value{
\code{get_all_tpmsample()} prints the count table without row numbers and invisibly returns a data frame with \code{Cancer}, \code{Normal} and \code{Tumor} columns. A missing group is filled with zero, and the final \code{Total} row sums each group's counts across cancers.

\code{panfam_sample_summary()} prints the count table and invisibly returns a \code{panfam_step}. Its tables are \code{result$tables$Project$sample_summary} and \code{result$tables$Project$sample_summary_metadata}. The latter records the data source, counting unit and selected cancer types. The result also contains \code{run_log}, \code{output_dirs} and \code{module}; the saved directory is \code{result$output_dirs[["Project"]]}.
}
\description{
Counts input sample rows by cancer type and Normal/Tumor group. The project interface automatically loads the selected input and saves the summary. It can be run immediately after setup, before ssGSEA or any downstream analysis.
}
\details{
Counts are numbers of input rows, not unique patients. No patient deduplication, gene matching, ssGSEA scoring or downstream complete-case filtering is performed. Actual sample sizes used by each downstream analysis are reported by that analysis.

The project wrapper applies \code{project$cancers} to either input source. For bundled TCGA/GTEx data, the source \code{project} and \code{sample_type} fields supply the cancer and group labels. \code{TCGA_tumor} becomes \code{Tumor}; \code{TCGA_normal} and \code{GTEx_normal} become \code{Normal}. Missing, blank or unsupported labels are rejected rather than silently omitted from counts.

This is a project-level step: each call summarizes the selected input once, regardless of the number of configured gene families. Files are saved under \code{project$output_dir/sample_summary_tcga/} or \code{project$output_dir/sample_summary_tcga_gtex/}. Each directory contains \code{sample_summary.csv}, \code{sample_summary_metadata.csv} and \code{result.rds}. The step produces tables only, with no figure. Repeated calls use numbered directories such as \code{sample_summary_tcga_002}, preserving earlier results.

The run log uses \code{Family = "Project"} and \code{Module = "sample_summary_tcga"} or \code{"sample_summary_tcga_gtex"}. The step appears in \code{\link{panfam_status}} and the index created by \code{\link{panfam_report}}. The standalone \code{get_all_tpmsample()} helper does not load bundled data, record project history or save files.

Invalid input causes an error in the standalone helper. The project wrapper records a failed run as \code{Status = "error"}, writes \code{error.txt} and emits a warning; its result has no count table for that run.
}
\seealso{\code{\link{panfam_setup}}, \code{\link{panfam_ssgsea}}, \code{\link{panfam_status}}, \code{\link{panfam_report}}}
\examples{
tpm <- data.frame(
  Cancer = c("BRCA", "BRCA", "LUAD", "LUAD"),
  Group = c("Tumor", "Tumor", "Normal", "Tumor")
)
counts <- get_all_tpmsample(tpm)

\dontrun{
project <- panfam_setup(
  output_dir = "results", genes = c("TP53", "MDM2"),
  family_name = "MyFamily", cancers = c("BRCA", "LUAD")
)
sample_result <- panfam_sample_summary(project)
sample_result$tables$Project$sample_summary
tcga_gtex_samples <- panfam_sample_summary(project, data_source = "tcga_gtex")
}
}
