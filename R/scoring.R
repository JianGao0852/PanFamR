.panfam_safe_names <- function(labels) {
  safe <- substr(gsub("[^A-Za-z0-9_.-]+", "_", labels), 1L, 80L)
  safe <- sub("[.]+$", "", safe)
  invalid <- !grepl("[A-Za-z0-9]", safe)
  safe[invalid] <- paste0("Family_", sprintf("%03d", which(invalid)))
  reserved <- grepl("^(CON|PRN|AUX|NUL|COM[1-9]|LPT[1-9])([.]|$)", safe, ignore.case = TRUE)
  safe[reserved] <- paste0("Family_", safe[reserved])
  used <- character()
  for (i in seq_along(safe)) {
    candidate <- safe[[i]]
    suffix <- 1L
    while (tolower(candidate) %in% used) {
      suffix <- suffix + 1L
      candidate <- paste0(safe[[i]], "_", suffix)
    }
    safe[[i]] <- candidate
    used <- c(used, tolower(candidate))
  }
  stats::setNames(safe, labels)
}

.panfam_gene_sets <- function(gene_sets = NULL, genes = NULL, family_name = "CustomFamily") {
  if (!is.null(genes)) {
    if (!is.null(gene_sets)) stop("Supply either `genes` or `gene_sets`, not both.", call. = FALSE)
    if (!is.character(family_name) || length(family_name) != 1L || is.na(family_name) ||
        !nzchar(trimws(family_name))) stop("`family_name` must be a non-empty string.", call. = FALSE)
    gene_sets <- stats::setNames(list(genes), family_name)
  }
  if (is.character(gene_sets)) gene_sets <- stats::setNames(list(gene_sets), family_name)
  if (!is.list(gene_sets) || !length(gene_sets) || is.null(names(gene_sets)) ||
      anyNA(names(gene_sets)) || any(!nzchar(trimws(names(gene_sets)))) || anyDuplicated(names(gene_sets))) {
    stop("Provide `genes` or a non-empty named list in `gene_sets`.", call. = FALSE)
  }
  lapply(gene_sets, function(x) {
    if (!is.character(x) || !length(x) || anyNA(x) || any(!nzchar(trimws(x)))) {
      stop("Each gene set must be a non-empty character vector of gene symbols without missing values.", call. = FALSE)
    }
    unique(trimws(x))
  })
}

.panfam_expression <- function(expression = NULL) {
  if (is.null(expression)) {
    objects <- panfam_load_data("all_tpm")
    if (!is.null(objects$tpm)) {
      expression <- objects$tpm
    } else {
      candidates <- objects[vapply(objects, function(x) {
        is.data.frame(x) && all(c("Cancer", "Group") %in% names(x))
      }, logical(1))]
      if (length(candidates) != 1L) stop("Cannot identify the bundled expression table.", call. = FALSE)
      expression <- candidates[[1L]]
    }
  }
  if (!is.data.frame(expression) || !all(c("Cancer", "Group") %in% names(expression))) {
    stop("`expression` must be a data frame with Cancer and Group columns and one sample per row.", call. = FALSE)
  }
  if (!nrow(expression) || anyDuplicated(names(expression)) || anyDuplicated(rownames(expression)) ||
      anyNA(expression$Cancer) || anyNA(expression$Group)) {
    stop("Expression input must have samples, unique gene and sample names, and complete Cancer/Group labels.", call. = FALSE)
  }
  if (!all(as.character(expression$Group) %in% c("Tumor", "Normal"))) {
    stop("`Group` must use the labels 'Tumor' and 'Normal'.", call. = FALSE)
  }
  expression
}

.panfam_score_matrix <- function(matrix, gene_sets, min_size = 1L) {
  if (!is.numeric(min_size) || length(min_size) != 1L || !is.finite(min_size) ||
      min_size < 1 || min_size > .Machine$integer.max || min_size != as.integer(min_size)) {
    stop("`min_size` must be a positive integer.", call. = FALSE)
  }
  if (!is.matrix(matrix) || !is.numeric(matrix) || nrow(matrix) < 2L || ncol(matrix) < 2L ||
      is.null(rownames(matrix)) || anyDuplicated(rownames(matrix)) || anyDuplicated(colnames(matrix))) {
    stop("The score matrix needs at least two genes and two samples with unique names.", call. = FALSE)
  }
  if (anyNA(matrix) || any(!is.finite(matrix))) stop("Expression values must be finite numbers without NA.", call. = FALSE)
  matched <- lapply(gene_sets, intersect, y = rownames(matrix))
  too_small <- lengths(matched) < min_size
  if (any(too_small)) {
    stop("Too few matched genes for: ", paste(names(matched)[too_small], collapse = ", "),
         ". Required minimum: ", min_size, ". Check gene symbols and min_size.", call. = FALSE)
  }
  params <- GSVA::ssgseaParam(exprData = matrix, geneSets = matched, minSize = min_size,
                            maxSize = Inf, alpha = 0.25, normalize = FALSE)
  background_notes <- character()
  raw <- withCallingHandlers(GSVA::gsva(params, verbose = FALSE, BPPARAM = BiocParallel::SerialParam()),
    warning = function(w) {
      if (grepl("rows with constant values throughout the samples", conditionMessage(w), fixed = TRUE)) {
        background_notes <<- c(background_notes, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    })
  ranges <- apply(raw, 1L, function(x) diff(range(x)))
  if (any(!is.finite(ranges) | ranges <= 0)) stop("A gene set produced constant or non-finite ssGSEA scores.", call. = FALSE)
  # Normalize each family by its own score range, matching separate single-set
  # calls in the original source. Ranking is shared across all requested sets.
  scores <- sweep(raw, 1L, ranges, "/")
  attr(scores, "background_notes") <- background_notes
  scores
}

#' Score arbitrary gene families using bundled or custom expression data
#' @param gene_sets A named list of gene vectors, or a single gene vector.
#' @param expression Optional sample-by-gene data frame with Cancer and Group columns.
#' @param genes Optional single gene vector, as an alternative to gene_sets.
#' @param family_name Label used with a single gene vector.
#' @param min_size Minimum number of matched genes; defaults to one.
#' @param output_dir Optional directory for score files.
#' @return A named list of sample-level score data frames.
#' @export
panfam_score <- function(gene_sets = NULL, expression = NULL, genes = NULL,
                         family_name = "CustomFamily", min_size = 1L, output_dir = NULL) {
  sets <- .panfam_gene_sets(gene_sets, genes, family_name)
  expression <- .panfam_expression(expression)
  gene_columns <- setdiff(names(expression), c("Cancer", "Group"))
  matrix <- as.matrix(expression[, gene_columns, drop = FALSE])
  suppressWarnings(storage.mode(matrix) <- "double")
  matrix <- t(matrix)
  message("Scoring ", length(sets), " gene set(s) in ", ncol(matrix), " samples ...")
  values <- .panfam_score_matrix(matrix, sets, min_size)
  scores <- lapply(names(sets), function(name) {
    data.frame(Score = as.numeric(values[name, rownames(expression)]), Cancer = as.character(expression$Cancer),
               Group = as.character(expression$Group), row.names = rownames(expression),
               check.names = FALSE)
  })
  names(scores) <- names(sets)
  attr(scores, "background_notes") <- attr(values, "background_notes")
  attr(scores, "gene_coverage") <- lapply(sets, function(x) {
    data.frame(Gene = x, Matched = x %in% gene_columns, stringsAsFactors = FALSE)
  })
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(output_dir)) stop("Cannot create output directory.", call. = FALSE)
    for (name in names(scores)) {
      safe <- .panfam_safe_names(names(scores))[[name]]
      df <- scores[[name]]
      save(df, file = file.path(output_dir, paste0(safe, "_ssGSEA_scores.Rdata")))
      .panfam_write_csv(data.frame(Sample = rownames(df), df, row.names = NULL),
                       file.path(output_dir, paste0(safe, "_ssGSEA_scores.csv")), row.names = FALSE)
    }
  }
  scores
}
