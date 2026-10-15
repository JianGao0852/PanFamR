\name{panfam_status}
\alias{panfam_status}
\title{Inspect Previously Executed Analysis Steps}
\usage{panfam_status(project)}
\arguments{\item{project}{Project returned by panfam_setup().}}
\value{A data frame with Family, Module, Run, Status, Rows, Message, Output_dir and Seconds. Empty before any analysis.}
\description{Returns the history of executed steps, including repeated runs, partial results and errors. Does not trigger any analysis. The same history is saved as run_log.csv.}
