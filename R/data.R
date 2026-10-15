# Bundled data accessors. Large RData files are split into <100 MB pieces so
# the repository can be stored on GitHub without Git LFS.

.panfam_data_index <- function() {
  c(
    all_tpm = "expression/all_tpm.Rdata",
    tcga_gtex = "expression/TCGA_GTEx_pancancer_mrna_pheno.rdata",
    paired_tpm = "expression/paired_tpm.Rdata",
    ALL_OS = "clinical/ALL_OS.Rdata",
    DSS_PFI_DFI = "clinical/DSS_PFI_DFI.Rdata",
    All_TMB = "auxiliary/All_TMB.Rdata",
    All_MSI = "auxiliary/All_MSI.Rdata",
    ssgsea_immucells = "auxiliary/ssgsea_immucells.Rdata",
    immuscore = "auxiliary/immuscore.Rdata",
    TME_Pathaway_GSVA_Score = "auxiliary/TME_Pathaway_GSVA_Score.Rdata",
    HALLMARK_GSVA_Score = "auxiliary/HALLMARK_GSVA_Score.Rdata"
  )
}

#' List data bundled with PanFamR
#'
#' @return A data frame with the bundled file name, category, and size.
#' @export
panfam_data_files <- function() {
  root <- system.file("extdata", package = "PanFamR")
  files <- list.files(root, recursive = TRUE, full.names = TRUE)
  files <- files[file.info(files)$isdir %in% FALSE]
  data.frame(
    file = substring(files, nchar(root) + 2L),
    size_mb = round(file.info(files)$size / 1024^2, 2),
    row.names = NULL,
    check.names = FALSE
  )
}

#' Get the installed path to a bundled data file
#'
#' @param file Relative bundled path, or the original logical name of a file
#'   distributed in numbered parts.
#' @return The absolute path to the installed file, or to its reconstructed
#'   copy in the user's R cache directory.
#' @export
panfam_data_path <- function(file) {
  if (!is.character(file) || length(file) != 1L || is.na(file) ||
      grepl("(^|[/\\\\])\\.\\.([/\\\\]|$)", file)) {
    stop("`file` must be a single safe relative path from panfam_data_files().", call. = FALSE)
  }
  .panfam_materialize_file(file)
}

.panfam_materialize <- function(name) {
  index <- .panfam_data_index()
  if (!name %in% names(index)) {
    stop("Unknown bundled dataset: ", name, ". See panfam_data_files().", call. = FALSE)
  }
  .panfam_materialize_file(unname(index[[name]]))
}

.panfam_materialize_file <- function(rel) {
  stem <- basename(rel)
  extdata <- system.file("extdata", package = "PanFamR")
  original <- file.path(extdata, rel)
  if (file.exists(original)) return(original)

  candidates <- list.files(dirname(original), full.names = FALSE)
  prefix <- paste0(stem, ".part")
  suffix <- substring(candidates, nchar(prefix) + 1L)
  keep <- startsWith(candidates, prefix) & grepl("^[0-9]+$", suffix)
  candidates <- candidates[keep]
  part_numbers <- as.numeric(suffix[keep])
  parts <- file.path(dirname(original), candidates[order(part_numbers)])
  if (!length(parts)) stop("Bundled data file is missing: ", rel, call. = FALSE)

  cache <- file.path(getOption("PanFamR.cache_dir", tools::R_user_dir("PanFamR", "cache")), "data")
  if (!dir.exists(cache) && !dir.create(cache, recursive = TRUE)) {
    stop("Cannot create PanFamR data cache: ", cache, call. = FALSE)
  }
  target <- file.path(cache, stem)
  if (!file.exists(target) || file.info(target)$size != sum(file.info(parts)$size)) {
    tmp <- tempfile(pattern = paste0(stem, "-"), tmpdir = cache)
    on.exit(unlink(tmp), add = TRUE)
    out <- file(tmp, open = "wb")
    out_open <- TRUE
    on.exit(if (out_open) close(out), add = TRUE)
    for (part in parts) {
      input <- file(part, open = "rb")
      tryCatch({
        repeat {
          block <- readBin(input, what = "raw", n = 8L * 1024L^2)
          if (!length(block)) break
          writeBin(block, out)
        }
      }, finally = close(input))
    }
    close(out)
    out_open <- FALSE
    if (!file.rename(tmp, target)) {
      if (!file.copy(tmp, target, overwrite = TRUE)) stop("Could not cache bundled data.", call. = FALSE)
    }
  }
  target
}

#' Load a bundled PanFamR dataset
#'
#' Large expression matrices are joined from package-internal file parts on the
#' first call and cached in the user's R cache directory. The returned value is
#' a named list of objects stored in the original RData file.
#'
#' @param name Dataset key. See `panfam_data_files()` and package README.
#' @return A named list containing the objects in the selected RData file.
#' @export
panfam_load_data <- function(name = "all_tpm") {
  if (!is.character(name) || length(name) != 1L || is.na(name)) {
    stop("`name` must be one dataset key.", call. = FALSE)
  }
  path <- .panfam_materialize(name)
  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)
  stats::setNames(lapply(object_names, get, envir = env, inherits = FALSE), object_names)
}

