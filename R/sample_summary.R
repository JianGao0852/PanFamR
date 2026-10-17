# Count expression-input rows before scoring or patient-level filtering.
.panfam_sample_counts <- function(tpm) {
  required <- c("Cancer", "Group")
  if (!is.data.frame(tpm) ||
      any(vapply(required, function(column) sum(names(tpm) == column) != 1L, logical(1)))) {
    stop("`tpm` must be a data frame with one Cancer column and one Group column.", call. = FALSE)
  }
  cancer <- as.character(tpm$Cancer)
  group <- as.character(tpm$Group)
  if (anyNA(cancer) || any(!nzchar(trimws(cancer)))) {
    stop("Cancer labels must be complete and non-empty.", call. = FALSE)
  }
  if (anyNA(group) || any(!group %in% c("Normal", "Tumor"))) {
    stop("Group must contain only 'Normal' and 'Tumor', without missing values.", call. = FALSE)
  }
  if ("Total" %in% cancer) {
    stop("'Total' is reserved for the final summary row; supply unaggregated sample rows.", call. = FALSE)
  }
  cancers <- sort(unique(cancer))
  counts <- table(factor(cancer, levels = cancers),
                  factor(group, levels = c("Normal", "Tumor")))
  result <- data.frame(Cancer = cancers,
                       Normal = as.integer(counts[, "Normal"]),
                       Tumor = as.integer(counts[, "Tumor"]),
                       stringsAsFactors = FALSE)
  total <- data.frame(Cancer = "Total", Normal = sum(result$Normal),
                      Tumor = sum(result$Tumor), stringsAsFactors = FALSE)
  result <- rbind(result, total)
  rownames(result) <- NULL
  result
}

#' Summarize normal and tumor sample counts by cancer
#'
#' Preserves the original function name and its printed table. Counts refer to
#' rows in the supplied input, without patient deduplication or gene filtering.
#' Both count columns are retained when one sample group is absent.
#' @param tpm A data frame with Cancer and Group columns, one sample per row.
#'   Group must contain Normal or Tumor. Gene-expression columns are optional.
#' @return Invisibly returns the printed Cancer, Normal and Tumor data frame,
#'   with a final Total row.
#' @export
get_all_tpmsample <- function(tpm) {
  result <- .panfam_sample_counts(tpm)
  print(result, row.names = FALSE)
  invisible(result)
}

.panfam_sample_summary_input <- function(project, data_source) {
  if (data_source == "tcga") {
    expression <- .panfam_project_expression(project)
    return(list(samples = expression[, c("Cancer", "Group"), drop = FALSE],
                source = project$data_source))
  }
  expression <- .panfam_aux_object("tcga_gtex")
  if (!is.data.frame(expression) || !all(c("project", "sample_type") %in% names(expression))) {
    stop("Bundled TCGA/GTEx data must contain project and sample_type columns.", call. = FALSE)
  }
  if (!is.null(project$cancers)) {
    unknown <- setdiff(project$cancers, unique(as.character(expression$project)))
    if (length(unknown)) {
      stop("Cancer labels absent from TCGA/GTEx input: ", paste(unknown, collapse = ", "), call. = FALSE)
    }
    expression <- expression[expression$project %in% project$cancers, , drop = FALSE]
  }
  group_map <- c(TCGA_normal = "Normal", GTEx_normal = "Normal", TCGA_tumor = "Tumor")
  sample_type <- as.character(expression$sample_type)
  if (anyNA(sample_type) || any(!sample_type %in% names(group_map))) {
    stop("TCGA/GTEx sample_type must use TCGA_normal, GTEx_normal or TCGA_tumor.", call. = FALSE)
  }
  list(samples = data.frame(Cancer = as.character(expression$project),
                            Group = unname(group_map[sample_type]), stringsAsFactors = FALSE),
       source = "Bundled TCGA/GTEx expression")
}

#' Save sample counts for a configured analysis project
#'
#' Summarizes the expression input once per project, independently of gene
#' families and ssGSEA scores. No scoring, patient deduplication or plotting is
#' performed. Repeated calls preserve earlier results in numbered directories.
#' @param project Project returned by panfam_setup().
#' @param data_source Either tcga (the project expression input, including custom
#'   expression when supplied) or tcga_gtex (the bundled combined input).
#' @return Invisibly returns a panfam_step. The printed table is also available
#'   at result$tables$Project$sample_summary. CSV files and result.rds are saved
#'   in a project-level sample_summary_tcga or sample_summary_tcga_gtex directory.
#' @export
panfam_sample_summary <- function(project, data_source = c("tcga", "tcga_gtex")) {
  .panfam_check_project(project)
  data_source <- match.arg(data_source)
  module <- paste0("sample_summary_", data_source)
  started <- proc.time()["elapsed"]
  history <- project$cache$history
  run <- sum(history$Family == "Project" & history$Module == module) + 1L
  directory <- .panfam_unique_directory(file.path(project$output_dir, module))
  .panfam_project_progress(project, paste0("[Project] ", module, " ..."))

  value <- tryCatch({
    input <- .panfam_sample_summary_input(project, data_source)
    counts <- .panfam_sample_counts(input$samples)
    metadata <- data.frame(
      Field = c("data_source", "source_label", "counting_unit", "patient_deduplication",
                "selected_cancers", "input_rows"),
      Value = c(data_source, input$source, "Expression-input sample rows", "None",
                paste(sort(unique(as.character(input$samples$Cancer))), collapse = ";"),
                as.character(nrow(input$samples))), stringsAsFactors = FALSE)
    list(sample_summary = counts, sample_summary_metadata = metadata)
  }, error = function(e) e)

  if (inherits(value, "error")) {
    status <- "error"
    message_text <- conditionMessage(value)
    rows <- 0L
    writeLines(message_text, file.path(directory, "error.txt"))
    value <- list()
  } else {
    status <- "complete"
    rows <- nrow(value$sample_summary)
    total <- value$sample_summary[rows, , drop = FALSE]
    message_text <- paste(sum(total$Normal, total$Tumor), "sample rows across",
                          rows - 1L, "cancer types; no patient deduplication")
    for (name in names(value)) {
      .panfam_write_csv(value[[name]], file.path(directory, paste0(name, ".csv")), row.names = FALSE)
    }
  }
  log <- data.frame(Family = "Project", Module = module, Run = run, Status = status,
                    Rows = rows, Message = message_text, Output_dir = directory,
                    Seconds = unname(proc.time()["elapsed"] - started), stringsAsFactors = FALSE)
  result <- structure(list(module = module, tables = list(Project = value), run_log = log,
                            output_dirs = c(Project = directory), output_dir = project$output_dir),
                       class = "panfam_step")
  saveRDS(list(scope = "Project", module = module, tables = value, run_log = log),
          file.path(directory, "result.rds"))
  .panfam_project_record(project, log)
  project$cache$results[[module]] <- result
  .panfam_project_progress(project, paste0("[Project] ", module, ": ", status, ". Results: ", directory))
  if (status == "error") {
    warning("Sample summary failed; see result$run_log and error.txt: ", message_text, call. = FALSE)
  } else {
    print(value$sample_summary, row.names = FALSE)
  }
  invisible(result)
}
