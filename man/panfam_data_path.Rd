\name{panfam_data_path}
\alias{panfam_data_path}
\title{Get the Installed Path to a Bundled Data File}
\usage{panfam_data_path(file)}
\arguments{\item{file}{Relative bundled path, or the original logical name of a file distributed in numbered parts.}}
\value{The absolute path to the installed file, or to its reconstructed copy in the user's R cache directory.}
\description{Returns the absolute path to a bundled file. Files distributed in numbered \code{.partNN} pieces are joined automatically in the user's R cache directory when their original logical name is requested. Existing complete files are returned directly.}
