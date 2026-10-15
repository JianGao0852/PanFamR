.panfam_absolute_path <- function(path, base = getwd()) {
  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(trimws(path))) {
    stop("Directory paths must be non-empty strings.", call. = FALSE)
  }
  path <- path.expand(path)
  if (!grepl("^([A-Za-z]:[/\\\\]|[/\\\\])", path)) path <- file.path(base, path)
  path
}

.panfam_create_directory <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(path)) stop("Cannot create directory: ", path, call. = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.panfam_unique_directory <- function(path) {
  candidate <- path
  counter <- 1L
  while (file.exists(candidate)) {
    counter <- counter + 1L
    candidate <- paste0(path, "_", sprintf("%03d", counter))
  }
  .panfam_create_directory(candidate)
}

.panfam_check_project <- function(project) {
  if (!inherits(project, "panfam_project") || !is.environment(project) ||
      !is.environment(project$cache)) {
    stop("Create a project with panfam_setup() first.", call. = FALSE)
  }
  invisible(project)
}

.panfam_empty_status <- function() {
  data.frame(Family = character(), Module = character(), Run = integer(),
             Status = character(), Rows = integer(), Message = character(),
             Output_dir = character(), Seconds = numeric(), stringsAsFactors = FALSE)
}

.panfam_project_progress <- function(project, text) {
  message(text)
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text, "\n",
      file = file.path(project$output_dir, "progress.log"), append = TRUE)
}

.panfam_project_record <- function(project, rows) {
  cache <- project$cache
  cache$history <- rbind(cache$history, rows)
  .panfam_write_csv(cache$history, file.path(project$output_dir, "run_log.csv"),
                   row.names = FALSE, na = "")
}

#' Configure a modular gene-family analysis project
#'
#' Configuration is supplied once. Each analysis function accepts this project,
#' loads its required bundled data and writes its own results. Configuration is
#' read-only; create another project to change genes or settings.
#' @param work_dir Working directory used to resolve a relative output_dir.
#' @param output_dir Directory for results. An existing project is preserved by
#'   creating a timestamped sibling directory.
#' @param genes Character vector of gene symbols for one custom family.
#' @param family_name Label for genes.
#' @param gene_sets Alternative named list of custom gene vectors.
#' @param cancers Optional vector of TCGA cancer labels; NULL selects all.
#' @param min_size Positive minimum matched gene count for ssGSEA.
#' @param correlation_method Pearson (source-script default) or Spearman correlation.
#' @param expression Optional custom sample-by-gene data frame.
#' @return A panfam_project passed to the individual analysis functions.
#' @export
panfam_setup <- function(work_dir = getwd(), output_dir = "PanFamR_results",
                          genes = NULL, family_name = "CustomFamily", gene_sets = NULL,
                          cancers = NULL, min_size = 1L,
                          correlation_method = c("pearson", "spearman"), expression = NULL) {
  sets <- .panfam_gene_sets(gene_sets, genes, family_name)
  method <- match.arg(correlation_method)
  if (!is.numeric(min_size) || length(min_size) != 1L || !is.finite(min_size) ||
      min_size < 1 || min_size > .Machine$integer.max || min_size != as.integer(min_size)) {
    stop("`min_size` must be a positive integer.", call. = FALSE)
  }
  if (!is.null(cancers) && (!is.character(cancers) || !length(cancers) || anyNA(cancers) ||
                           any(!nzchar(trimws(cancers))))) {
    stop("`cancers` must be NULL or a non-empty vector of cancer labels.", call. = FALSE)
  }
  if (!is.null(expression)) expression <- .panfam_expression(expression)
  work <- .panfam_create_directory(.panfam_absolute_path(work_dir))
  root <- .panfam_absolute_path(output_dir, work)
  if (any(file.exists(file.path(root, c("project_config.rds", "report.html", "run_log.csv"))))) {
    root <- .panfam_unique_directory(paste0(root, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  } else {
    root <- .panfam_create_directory(root)
  }
  project <- new.env(parent = emptyenv())
  project$work_dir <- work
  project$output_dir <- root
  project$gene_sets <- sets
  project$family_directories <- .panfam_safe_names(names(sets))
  project$cancers <- if (is.null(cancers)) NULL else unique(cancers)
  project$min_size <- as.integer(min_size)
  project$correlation_method <- method
  project$data_source <- if (is.null(expression)) "Bundled all_tpm (TCGA)" else "User-supplied expression"
  project$created <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "Asia/Shanghai")
  project$cache <- new.env(parent = emptyenv())
  project$cache$expression <- expression
  project$cache$expression_ready <- FALSE
  project$cache$scores <- NULL
  project$cache$score_runs <- 0L
  project$cache$results <- list()
  project$cache$history <- .panfam_empty_status()
  class(project) <- "panfam_project"
  config <- as.list.environment(project, all.names = TRUE)
  config$cache <- NULL
  config$package_version <- as.character(utils::packageVersion("PanFamR"))
  saveRDS(config, file.path(root, "project_config.rds"))
  .panfam_write_csv(.panfam_bind(lapply(names(sets), function(name) {
    data.frame(Family = name, Gene = sets[[name]], stringsAsFactors = FALSE)
  })), file.path(root, "gene_families.csv"), row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(root, "sessionInfo.txt"))
  lockEnvironment(project, bindings = TRUE)
  .panfam_project_progress(project, paste("Project configured:", length(sets), "family/families. Results:", root))
  project
}

.panfam_project_expression <- function(project) {
  .panfam_check_project(project)
  cache <- project$cache
  if (!cache$expression_ready) {
    .panfam_project_progress(project, "Loading expression data ...")
    expression <- .panfam_expression(cache$expression)
    if (!is.null(project$cancers)) {
      unknown <- setdiff(project$cancers, unique(as.character(expression$Cancer)))
      if (length(unknown)) stop("Cancer labels absent from expression: ", paste(unknown, collapse = ", "), call. = FALSE)
      expression <- expression[expression$Cancer %in% project$cancers, , drop = FALSE]
    }
    cache$expression <- expression
    cache$expression_ready <- TRUE
  }
  project$cache$expression
}

.panfam_project_cancers <- function(project) {
  sort(unique(as.character(.panfam_project_expression(project)$Cancer)))
}

#' Compute and save ssGSEA scores for a configured project
#' @param project Project returned by panfam_setup().
#' @return Named list of score data frames, one per family. Repeated calls reuse
#'   scores in this project. Other score-based functions compute them if needed.
#' @export
panfam_ssgsea <- function(project) {
  .panfam_check_project(project)
  cache <- project$cache
  if (!is.null(project$cache$scores)) {
    message("Reusing ssGSEA scores for this project.")
    return(project$cache$scores)
  }
  started <- proc.time()["elapsed"]
  expression <- .panfam_project_expression(project)
  scores <- panfam_score(project$gene_sets, expression = expression, min_size = project$min_size)
  coverage <- attr(scores, "gene_coverage")
  rows <- list()
  for (name in names(scores)) {
    directory <- .panfam_unique_directory(file.path(project$output_dir, project$family_directories[[name]], "ssgsea"))
    .panfam_write_csv(data.frame(Sample = rownames(scores[[name]]), scores[[name]], row.names = NULL),
                     file.path(directory, "scores.csv"), row.names = FALSE)
    .panfam_write_csv(coverage[[name]], file.path(directory, "gene_coverage.csv"), row.names = FALSE)
    df <- scores[[name]]
    save(df, file = file.path(directory, "ssGSEA_scores.Rdata"))
    saveRDS(list(scores = df, gene_coverage = coverage[[name]],
                 background_notes = attr(scores, "background_notes")), file.path(directory, "result.rds"))
    rows[[name]] <- data.frame(Family = name, Module = "ssgsea", Run = 1L,
                               Status = if (all(coverage[[name]]$Matched)) "complete" else "partial",
                               Rows = nrow(df), Message = paste(sum(coverage[[name]]$Matched), "of", nrow(coverage[[name]]), "genes matched"),
                               Output_dir = directory, Seconds = unname(proc.time()["elapsed"] - started), stringsAsFactors = FALSE)
  }
  cache$scores <- scores
  cache$score_runs <- cache$score_runs + 1L
  .panfam_project_record(project, .panfam_bind(rows))
  .panfam_project_progress(project, paste("ssGSEA finished. Results:", project$output_dir))
  scores
}

.panfam_run_step <- function(project, module, analyze, needs_scores = TRUE) {
  .panfam_check_project(project)
  cache <- project$cache
  if (needs_scores && is.null(project$cache$scores)) panfam_ssgsea(project)
  tables <- logs <- directories <- list()
  for (name in names(project$gene_sets)) {
    started <- proc.time()["elapsed"]
    previous <- project$cache$history
    run <- sum(previous$Family == name & previous$Module == module) + 1L
    directory <- .panfam_unique_directory(file.path(project$output_dir, project$family_directories[[name]], module))
    directories[[name]] <- directory
    plot_dir <- .panfam_create_directory(file.path(directory, "plots"))
    .panfam_project_progress(project, paste0("[", name, "] ", module, " ..."))
    value <- tryCatch(analyze(name, plot_dir), error = function(e) e)
    if (inherits(value, "error")) {
      summary <- list(status = "error", rows = 0L, message = conditionMessage(value))
      writeLines(conditionMessage(value), file.path(directory, "error.txt"))
      value <- list()
    } else {
      summary <- .panfam_module_summary(value)
      for (table in names(value)) {
        if (ncol(value[[table]])) .panfam_write_csv(value[[table]], file.path(directory, paste0(table, ".csv")),
                                                  row.names = FALSE, na = "")
      }
    }
    tables[[name]] <- value
    logs[[name]] <- data.frame(Family = name, Module = module, Run = run, Status = summary$status,
                               Rows = summary$rows, Message = summary$message, Output_dir = directory,
                               Seconds = unname(proc.time()["elapsed"] - started), stringsAsFactors = FALSE)
    saveRDS(list(family = name, module = module, tables = value, run_log = logs[[name]]), file.path(directory, "result.rds"))
    .panfam_project_progress(project, paste0("[", name, "] ", module, ": ", summary$status, ". Results: ", directory))
  }
  result <- structure(list(module = module, tables = tables, run_log = .panfam_bind(logs),
                           output_dirs = unlist(directories), output_dir = project$output_dir), class = "panfam_step")
  .panfam_project_record(project, result$run_log)
  cache$results[[module]] <- result
  if (any(result$run_log$Status == "error")) warning("Analysis failed for some families; see result$run_log and error.txt.", call. = FALSE)
  result
}

#' Run one configured analysis step
#' @param project Project returned by panfam_setup().
#' @param endpoints Survival endpoints to run: OS, DSS, PFI and/or DFI.
#' @return A panfam_step with tables per family, run_log and output_dirs.
#' @name panfam_steps
#' @export
panfam_expression <- function(project) {
  .panfam_run_step(project, "expression", function(name, directory) {
    .panfam_expression_analysis(project$cache$scores[[name]], .panfam_project_expression(project),
                                project$gene_sets[[name]], directory, name)
  })
}

#' @rdname panfam_steps
#' @export
panfam_tcga_gtex <- function(project) {
  .panfam_run_step(project, "tcga_gtex", function(name, directory) {
    .panfam_tcga_gtex_analysis(project$gene_sets[[name]], .panfam_project_cancers(project), directory, name)
  }, needs_scores = FALSE)
}

#' @rdname panfam_steps
#' @export
panfam_paired <- function(project) {
  .panfam_run_step(project, "paired", function(name, directory) .panfam_paired_analysis(project$cache$scores[[name]], directory, name))
}

#' @rdname panfam_steps
#' @export
panfam_survival <- function(project, endpoints = c("OS", "DSS", "PFI", "DFI")) {
  if (!is.character(endpoints) || !length(endpoints) || anyNA(endpoints) ||
      any(!endpoints %in% c("OS", "DSS", "PFI", "DFI"))) {
    stop("`endpoints` must contain OS, DSS, PFI and/or DFI.", call. = FALSE)
  }
  endpoints <- unique(endpoints)
  .panfam_run_step(project, paste0("survival_", paste(endpoints, collapse = "_")), function(name, directory) {
    .panfam_survival_analysis(project$cache$scores[[name]], directory, name, endpoints = endpoints)
  })
}

#' @rdname panfam_steps
#' @export
panfam_gene_survival <- function(project, endpoints = "OS") {
  if (!is.character(endpoints) || !length(endpoints) || anyNA(endpoints) ||
      any(!endpoints %in% c("OS", "DSS", "PFI", "DFI"))) {
    stop("`endpoints` must contain OS, DSS, PFI and/or DFI.", call. = FALSE)
  }
  endpoints <- unique(endpoints)
  .panfam_run_step(project, paste0("gene_survival_", paste(endpoints, collapse = "_")), function(name, directory) {
    expression <- .panfam_project_expression(project)
    .panfam_gene_survival_analysis(expression[, c("Cancer", "Group"), drop = FALSE], expression,
                                   project$gene_sets[[name]], directory, name, endpoints = endpoints)
  }, needs_scores = FALSE)
}

#' @rdname panfam_steps
#' @export
panfam_tumor_scores <- function(project) {
  .panfam_run_step(project, "tumor_scores", function(name, directory) {
    .panfam_tumor_scores_analysis(project$cache$scores[[name]], directory, name)
  })
}

#' @rdname panfam_steps
#' @export
panfam_roc <- function(project) {
  .panfam_run_step(project, "roc", function(name, directory) {
    .panfam_roc_analysis(project$cache$scores[[name]], directory, name)
  })
}

#' @rdname panfam_steps
#' @export
panfam_tmb <- function(project) {
  .panfam_run_step(project, "tmb", function(name, directory) {
    .panfam_biomarker_analysis(project$cache$scores[[name]], "All_TMB", "TMB", directory, project$correlation_method)
  })
}

#' @rdname panfam_steps
#' @export
panfam_msi <- function(project) {
  .panfam_run_step(project, "msi", function(name, directory) {
    .panfam_biomarker_analysis(project$cache$scores[[name]], "All_MSI", "MSI", directory, project$correlation_method)
  })
}

.panfam_immune_step <- function(project, module, panel, needs_expression = FALSE) {
  .panfam_run_step(project, module, function(name, directory) {
    expression <- if (needs_expression) .panfam_project_expression(project) else NULL
    .panfam_immune_analysis(project$cache$scores[[name]], expression, directory,
                            project$correlation_method, panels = panel)
  })
}

#' @rdname panfam_steps
#' @export
panfam_immune_cells <- function(project) .panfam_immune_step(project, "immune_cells", "Immune_cells")
#' @rdname panfam_steps
#' @export
panfam_immune_scores <- function(project) .panfam_immune_step(project, "immune_scores", "Immune_scores")
#' @rdname panfam_steps
#' @export
panfam_tme <- function(project) .panfam_immune_step(project, "tme", "TME_pathways")
#' @rdname panfam_steps
#' @export
panfam_hallmark <- function(project) .panfam_immune_step(project, "hallmark", "Hallmark_pathways")
#' @rdname panfam_steps
#' @export
panfam_checkpoints <- function(project) .panfam_immune_step(project, "checkpoints", "Immune checkpoints", TRUE)
#' @rdname panfam_steps
#' @export
panfam_chemokines <- function(project) .panfam_immune_step(project, "chemokines", "Chemokines", TRUE)
#' @rdname panfam_steps
#' @export
panfam_chemokine_receptors <- function(project) .panfam_immune_step(project, "chemokine_receptors", "Chemokine receptors", TRUE)
#' @rdname panfam_steps
#' @export
panfam_immune_stimulators <- function(project) .panfam_immune_step(project, "immune_stimulators", "Immune stimulators", TRUE)
#' @rdname panfam_steps
#' @export
panfam_immune_inhibitors <- function(project) .panfam_immune_step(project, "immune_inhibitors", "Immune inhibitors", TRUE)

#' @rdname panfam_steps
#' @export
panfam_immunotherapy <- function(project) {
  .panfam_run_step(project, "immunotherapy", function(name, directory) {
    .panfam_immunotherapy_analysis(project$gene_sets[[name]], project$min_size, directory, name)
  }, needs_scores = FALSE)
}

#' @rdname panfam_steps
#' @export
panfam_alterations <- function(project) {
  .panfam_run_step(project, "alterations", function(name, directory) {
    .panfam_alterations_analysis(project$gene_sets[[name]], .panfam_project_cancers(project), directory, name)
  }, needs_scores = FALSE)
}

#' @rdname panfam_steps
#' @export
panfam_cnv <- function(project) {
  .panfam_run_step(project, "cnv", function(name, directory) {
    .panfam_cnv_methylation_analysis(project$gene_sets[[name]], .panfam_project_cancers(project), directory, associations = "CNV")
  }, needs_scores = FALSE)
}

#' @rdname panfam_steps
#' @export
panfam_methylation <- function(project) {
  .panfam_run_step(project, "methylation", function(name, directory) {
    .panfam_cnv_methylation_analysis(project$gene_sets[[name]], .panfam_project_cancers(project), directory, associations = "Methylation")
  }, needs_scores = FALSE)
}

#' Inspect the status of completed analysis steps
#' @param project Project returned by panfam_setup().
#' @return Data frame of all executed steps, including repeats and failures.
#' @export
panfam_status <- function(project) {
  .panfam_check_project(project)
  project$cache$history
}

#' Write a report for previously run steps
#' @param project Project returned by panfam_setup().
#' @return Path to report.html. No analysis or scoring is triggered.
#' @export
panfam_report <- function(project) {
  .panfam_check_project(project)
  root <- project$output_dir
  history <- panfam_status(project)
  parts <- c('<!doctype html><html lang="en"><meta charset="utf-8">',
             '<meta name="viewport" content="width=device-width, initial-scale=1"><title>PanFamR step results</title>',
             '<style>body{font:16px/1.6 system-ui,sans-serif;max-width:1200px;margin:32px auto;padding:0 24px;color:#243344}table{border-collapse:collapse;font-size:13px;width:100%}td,th{border:1px solid #d9e1e8;padding:8px}th{background:#eff5f8}img{max-width:100%}details{margin:18px 0}a{color:#166ea0}</style>',
             '<h1>PanFamR: results by analysis step</h1>',
             paste0('<p>Project created: ', .panfam_html(project$created), '</p>'),
             paste0('<p>Data: ', .panfam_html(project$data_source), '. Families: ', .panfam_html(paste(names(project$gene_sets), collapse = ", ")), '.</p>'),
             '<p>This index includes only steps you have executed. Missing source genes/assays and insufficient samples are recorded in coverage/status tables. Bundled mutation, CNV and methylation inputs cover only the genes available in the original source project.</p>',
             '<h2>Run status</h2>', .panfam_html_table(history))
  relative <- function(path) substring(path, nchar(root) + 2L)
  for (i in seq_len(nrow(history))) {
    row <- history[i, , drop = FALSE]
    directory <- row$Output_dir
    files <- list.files(directory, pattern = "[.](csv|rds|Rdata|pdf|txt)$", recursive = TRUE, full.names = TRUE)
    parts <- c(parts, paste0('<h2>', .panfam_html(row$Family), ' / ', .panfam_html(row$Module), ' (run ', row$Run, ')</h2><ul>'),
               vapply(files, function(file) paste0('<li><a href="', .panfam_html(relative(file)), '">', .panfam_html(substring(file, nchar(directory) + 2L)), '</a></li>'), character(1)), '</ul>')
    for (file in list.files(directory, pattern = "[.]png$", recursive = TRUE, full.names = TRUE)) {
      parts <- c(parts, paste0('<details><summary>', .panfam_html(basename(file)), '</summary><img alt="',
                              .panfam_html(basename(file)), '" src="', .panfam_html(relative(file)), '"></details>'))
    }
  }
  parts <- c(parts, '<p><a href="gene_families.csv">Gene families</a> | <a href="run_log.csv">Run log</a> | <a href="sessionInfo.txt">R session</a></p></html>')
  path <- file.path(root, "report.html")
  writeLines(enc2utf8(parts), path, useBytes = TRUE)
  message("Report for completed steps: ", path)
  invisible(path)
}

#' @export
print.panfam_project <- function(x, ...) {
  cat("PanFamR project:", length(x$gene_sets), "family/families\n")
  cat("Families:", paste(names(x$gene_sets), collapse = ", "), "\n")
  cat("Results:", x$output_dir, "\n")
  cat("Executed family steps:", nrow(x$cache$history), "\n")
  invisible(x)
}

#' @export
print.panfam_step <- function(x, ...) {
  cat("PanFamR step:", x$module, "\n")
  print(x$run_log[, c("Family", "Status", "Rows", "Output_dir"), drop = FALSE], row.names = FALSE)
  invisible(x)
}
