\name{panfam_report}
\alias{panfam_report}
\title{Build a Report for Previously Executed Steps}
\usage{panfam_report(project)}
\arguments{\item{project}{Project returned by panfam_setup().}}
\value{The absolute path to report.html, returned invisibly.}
\description{Creates or updates an HTML index of results from the steps already executed, with run status, table links and figure previews. It does not run any analysis or scoring.}
\examples{
\dontrun{
report_file <- panfam_report(project)
utils::browseURL(report_file)
}
}
