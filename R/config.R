.panfam_paths <- function() {
  bundled <- system.file("extdata", package = "PanFamR")
  getOption("PanFamR.paths", list(
    analysis = file.path(bundled, "auxiliary"),
    clinical = file.path(bundled, "clinical")
  ))
}

#' Configure locations for analysis data
#'
#' Plotting functions use bundled precomputed TCGA and clinical tables by
#' default. This helper optionally points them to user-supplied replacements.
#' The analysis directory should contain files such as `All_TMB.Rdata`,
#' `All_MSI.Rdata`, and `paired_tpm.Rdata`; the clinical directory should
#' contain `ALL_OS.Rdata` and `DSS_PFI_DFI.Rdata`.
#'
#' @param analysis_dir Directory containing precomputed analysis data.
#' @param clinical_dir Directory containing clinical survival data.
#' @return The configured paths, invisibly.
#' @export
#' @examples
#' \dontrun{
#' panfam_set_paths("data/analysis", "data/clinical")
#' }
panfam_set_paths <- function(analysis_dir = NULL, clinical_dir = NULL) {
  paths <- .panfam_paths()
  if (!is.null(analysis_dir)) paths$analysis <- normalizePath(analysis_dir, mustWork = TRUE)
  if (!is.null(clinical_dir)) paths$clinical <- normalizePath(clinical_dir, mustWork = TRUE)
  options(PanFamR.paths = paths)
  invisible(paths)
}

.panfam_path <- function(area, filename) {
  paths <- .panfam_paths()
  directory <- paths[[area]]
  if (is.null(directory) || !nzchar(directory)) {
    stop(
      sprintf("Configure the '%s' data directory with panfam_set_paths() before running this analysis.", area),
      call. = FALSE
    )
  }
  path <- file.path(directory, filename)
  if (!file.exists(path)) stop("Required data file not found: ", path, call. = FALSE)
  path
}

.panfam_prefix <- function() {
  basename(normalizePath(getwd(), winslash = "/", mustWork = FALSE))
}

#' Score one or more gene families with ssGSEA
#'
#' @param tpm A data frame with `Cancer` and `Group` columns, one sample per
#'   row, and gene symbols in the remaining column names.
#' @param gene_sets A named list of character vectors, one vector per gene set.
#' @param output_dir Directory in which score `.Rdata` files are written.
#' @return A named list of sample-level score data frames.
#' @export
#' @examples
#' \dontrun{
#' scores <- run_ssGSEA(tpm, list(COX = c("COX6A1", "COX4I1", "COX5A", "COX5B", "COX6B1")))
#' }
run_ssGSEA <- function(tpm, gene_sets, output_dir = ".") {
  if (!is.list(gene_sets) || !length(gene_sets) || is.null(names(gene_sets)) ||
      any(!nzchar(names(gene_sets)))) {
    stop("`gene_sets` must be a non-empty named list of gene vectors.", call. = FALSE)
  }
  if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
    stop("Could not create output directory: ", output_dir, call. = FALSE)
  }
  safe_names <- gsub("[^A-Za-z0-9_.-]+", "_", names(gene_sets))
  if (anyDuplicated(safe_names)) stop("Gene set names must remain unique after filename cleanup.", call. = FALSE)

  result <- lapply(seq_along(gene_sets), function(i) {
    path <- file.path(output_dir, paste0(safe_names[[i]], "_ssGSEA_scores.Rdata"))
    ssGsea_Score(tpm, gene_sets[[i]], output_file = path, family_name = names(gene_sets)[[i]])
  })
  names(result) <- names(gene_sets)
  result
}
