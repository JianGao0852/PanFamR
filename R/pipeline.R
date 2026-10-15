.panfam_html <- function(x) {
  x <- gsub("&", "&amp;", as.character(x), fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  gsub('"', "&quot;", x, fixed = TRUE)
}

.panfam_html_table <- function(x) {
  if (!nrow(x)) return("<p>No records.</p>")
  head <- paste0("<tr>", paste0("<th>", .panfam_html(names(x)), "</th>", collapse = ""), "</tr>")
  rows <- apply(x, 1L, function(row) paste0("<tr>", paste0("<td>", .panfam_html(row), "</td>", collapse = ""), "</tr>"))
  paste0("<table>", head, paste(rows, collapse = ""), "</table>")
}

.panfam_report <- function(result) {
  root <- result$output_dir
  parts <- c('<!doctype html><html lang="zh-CN"><meta charset="utf-8">',
             '<meta name="viewport" content="width=device-width, initial-scale=1"><title>PanFamR report</title>',
             '<style>body{font:16px/1.6 system-ui,sans-serif;max-width:1250px;margin:36px auto;padding:0 24px;color:#243344}h1,h2{color:#164762}table{border-collapse:collapse;width:100%;font-size:13px;margin:16px 0}th,td{border:1px solid #d9e1e8;padding:8px;text-align:left}th{background:#eff5f8}img{max-width:100%;height:auto}a{color:#166ea0}details{margin:18px 0}code{background:#eff5f8;padding:3px}li{margin:4px 0}.note{padding:16px;background:#eff5f8;border-left:4px solid #397fb6}</style>',
             '<h1>PanFamR \u6cdb\u764c\u5206\u6790\u62a5\u544a / Pan-cancer analysis</h1>',
             paste0('<p>Generated: ', .panfam_html(result$settings$created), ' | R ', .panfam_html(getRversion()), '</p>'),
             paste0('<p>Expression input: ', .panfam_html(result$settings$data_source), ' | Samples: ', nrow(result$scores[[1L]]),
                    ' | Cancer types: ', length(unique(result$scores[[1L]]$Cancer)), '</p>'),
             '<div class="note">\u6bcf\u4e2a\u7ed3\u679c\u76ee\u5f55\u5305\u542b CSV \u7ed3\u679c\u3001PDF/PNG \u56fe\u4ee5\u53ca\u57fa\u56e0\u8986\u76d6\u8868\u3002\u8bf7\u5148\u67e5\u770b\u4e0b\u9762\u7684\u8fd0\u884c\u72b6\u6001\uff1a\u7f3a\u5931\u6570\u636e\u548c\u672a\u80fd\u5b8c\u6210\u7684\u5206\u6790\u4f1a\u660e\u786e\u5217\u51fa\u3002CNV\u3001\u7532\u57fa\u5316\u548c\u7a81\u53d8\u8f93\u5165\u53ea\u8986\u76d6\u539f\u9879\u76ee\u5df2\u6709\u7684\u57fa\u56e0\uff0c\u8986\u76d6\u4e0d\u8db3\u65f6\u4e0d\u4f1a\u8865\u9020\u7ed3\u679c\u3002</div>',
             '<h2>\u8fd0\u884c\u72b6\u6001 / Run status</h2>', .panfam_html_table(result$run_log),
             '<h2>\u5206\u6790\u89c4\u5219 / Methods</h2>',
             '<p>ssGSEA uses alpha=0.25, with each family normalized by its own score range. Tumor-versus-normal tests use two-sided Wilcoxon rank-sum tests; paired tests use matched TCGA patient IDs. Effects are mean/median differences and rank-biserial effects, without assuming the input expression scale is logarithmic.</p>',
             paste0('<p>Associations use ', .panfam_html(result$settings$correlation_method), ' correlation. Survival uses a univariate Cox model per one score standard deviation and a fixed median split for log-rank/KM analysis. Low: score \u2264 median; high: score &gt; median. No samples are manually reassigned. Models require at least 10 complete samples and 5 events; groups require at least 3 samples. Each newly computed result table uses Benjamini-Hochberg adjustment across its valid tests. Source CNV/methylation FDR values are retained.</p>'),
             '<p>TCGA survival and correlation analyses retain one tumor sample per patient, with exact ID matches preferred and unambiguous TCGA prefixes used where necessary. Immunotherapy analyses use baseline samples and exclude known non-immunotherapy controls; unannotated treatment stages and unavailable OS data are recorded. Duplicate gene-symbol rows are averaged, and genes with non-finite expression in selected cohort samples are excluded with counts and family coverage recorded. Expression values are not imputed. Unknown response/event labels remain missing. Scores are interpreted within each cohort. TCGA/GTEx tumor-normal comparisons pool source-labelled normal samples and do not correct study/batch differences.</p>')
  for (name in names(result$scores)) {
    safe <- .panfam_safe_names(names(result$scores))[[name]]
    dir <- file.path(root, safe)
    coverage <- result$gene_coverage[[name]]
    parts <- c(parts, paste0('<h2>', .panfam_html(name), '</h2>'),
               paste0('<p>Matched ', sum(coverage$Matched), ' / ', nrow(coverage), ' requested genes.</p>'),
               '<details><summary>Gene coverage</summary>', .panfam_html_table(coverage), '</details>')
    plots <- list.files(file.path(dir, "plots"), pattern = "[.]png$", full.names = FALSE)
    for (plot in plots) {
      path <- paste(safe, "plots", plot, sep = "/")
      parts <- c(parts, paste0('<details', if (plot == "score_distribution.png") ' open' else '',
                              '><summary>', .panfam_html(sub("[.]png$", "", plot)), '</summary><img alt="',
                              .panfam_html(plot), '" src="', path, '"></details>'))
    }
    files <- list.files(dir, pattern = "[.](csv|pdf)$", recursive = TRUE)
    links <- vapply(files, function(file) paste0('<li><a href="', safe, '/', file, '">', .panfam_html(file), '</a></li>'), character(1))
    parts <- c(parts, '<details><summary>\u4e0b\u8f7d\u8868\u683c\u548c\u56fe / Tables and figures</summary><ul>', links, '</ul></details>')
  }
  parts <- c(parts, '<p>Run metadata: <a href="run_log.csv">run_log.csv</a>, <a href="sessionInfo.txt">sessionInfo.txt</a>, <a href="analysis_result.rds">analysis_result.rds</a>.</p>', '</html>')
  writeLines(enc2utf8(parts), file.path(root, "report.html"), useBytes = TRUE)
}

.panfam_module_summary <- function(tables) {
  rows <- sum(vapply(tables, nrow, integer(1)))
  p <- unlist(lapply(tables, function(x) if ("p_value" %in% names(x)) x$p_value else numeric()), use.names = FALSE)
  coverage <- unlist(lapply(tables, function(x) if ("Available" %in% names(x)) x$Available else logical()), use.names = FALSE)
  reasons <- unique(unlist(lapply(tables, function(x) if ("reason" %in% names(x)) as.character(x$reason) else character()), use.names = FALSE))
  reasons <- reasons[!is.na(reasons) & nzchar(reasons)]
  status <- if (length(coverage) && !any(coverage)) "skipped" else if (length(coverage) && any(!coverage)) "partial" else "complete"
  if (!rows) status <- "skipped"
  if (length(p) && !any(is.finite(p))) status <- "skipped"
  if (length(reasons) && status == "complete") status <- "partial"
  message <- paste0(rows, " rows; ", sum(is.finite(p)), " valid statistical tests")
  if (length(coverage) && any(!coverage)) message <- paste0(message, "; ", sum(!coverage), " unavailable gene/assay entries")
  if (length(reasons)) message <- paste0(message, "; ", paste(reasons, collapse = " | "))
  list(status = status, rows = rows, message = message)
}

#' Run a complete automated pan-cancer gene-family workflow
#' @param gene_sets Named list of gene vectors, or a single gene vector.
#' @param expression Optional expression data frame. The bundled TCGA table is the default.
#' @param output_dir Result directory. Existing reports cause a timestamped sibling directory to be used. NULL returns scores only.
#' @param genes Single gene vector, as an alternative to gene_sets.
#' @param family_name Label for a single gene vector.
#' @param analyses Modules to run. Defaults to all supported modules; character(0) runs scoring only.
#' @param cancers Optional cancer labels to analyze.
#' @param min_size Minimum number of matched genes for ssGSEA.
#' @param correlation_method Correlation method, either spearman or pearson.
#' @return A panfam_result containing scores, result tables, gene coverage, run status, and output paths.
#' @export
panfam_analyze <- function(gene_sets = NULL, expression = NULL, output_dir = "PanFamR_results",
                           genes = NULL, family_name = "CustomFamily",
                           analyses = c("expression", "tcga_gtex", "paired", "survival", "tmb", "msi", "immune", "immunotherapy", "alterations", "cnv_methylation"),
                           cancers = NULL, min_size = 1L, correlation_method = c("spearman", "pearson")) {
  sets <- .panfam_gene_sets(gene_sets, genes, family_name)
  correlation_method <- match.arg(correlation_method)
  supported <- c("expression", "tcga_gtex", "paired", "survival", "tmb", "msi", "immune", "immunotherapy", "alterations", "cnv_methylation")
  if (!is.character(analyses) || anyNA(analyses) || any(!analyses %in% supported)) {
    stop("Unknown analysis module. Supported: ", paste(supported, collapse = ", "), call. = FALSE)
  }
  analyses <- unique(analyses)
  bundled <- is.null(expression)
  if (is.null(output_dir)) return(panfam_score(sets, expression = expression, min_size = min_size))
  if (!is.character(output_dir) || length(output_dir) != 1L || is.na(output_dir) || !nzchar(output_dir)) {
    stop("`output_dir` must be one non-empty directory path.", call. = FALSE)
  }
  if (file.exists(file.path(output_dir, "report.html"))) {
    output_dir <- paste0(output_dir, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Cannot create result directory: ", output_dir, call. = FALSE)
  root <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  log_path <- file.path(root, "progress.log")
  progress <- function(text) {
    message(text)
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), text, "\n", file = log_path, append = TRUE)
  }
  started <- proc.time()["elapsed"]
  progress("Loading expression data ...")
  expression <- .panfam_expression(expression)
  if (!is.null(cancers)) {
    if (!is.character(cancers) || anyNA(cancers) || !length(cancers) || any(!cancers %in% unique(expression$Cancer))) {
      stop("Every requested cancer must occur in expression$Cancer.", call. = FALSE)
    }
    expression <- expression[expression$Cancer %in% cancers, , drop = FALSE]
  }
  progress(paste("Scoring", length(sets), "family/families across", nrow(expression), "samples ..."))
  scores <- panfam_score(sets, expression = expression, min_size = min_size)
  coverage <- attr(scores, "gene_coverage")
  logs <- list()
  tables <- list()
  for (name in names(sets)) {
    safe <- .panfam_safe_names(names(sets))[[name]]
    directory <- file.path(root, safe)
    plot_dir <- file.path(directory, "plots")
    dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
    .panfam_write_csv(data.frame(Sample = rownames(scores[[name]]), scores[[name]], row.names = NULL), file.path(directory, "scores.csv"), row.names = FALSE)
    .panfam_write_csv(coverage[[name]], file.path(directory, "gene_coverage.csv"), row.names = FALSE)
    df <- scores[[name]]
    save(df, file = file.path(directory, "ssGSEA_scores.Rdata"))
    logs[[paste(name, "score")]] <- data.frame(Family = name, Module = "score", Status = if (all(coverage[[name]]$Matched)) "complete" else "partial",
                                             Rows = nrow(df), Message = paste(sum(coverage[[name]]$Matched), "of", nrow(coverage[[name]]), "requested genes matched"))
    tables[[name]] <- list()
    for (module in analyses) {
      progress(paste0("[", name, "] ", module, " ..."))
      value <- tryCatch(switch(module,
        expression = .panfam_expression_analysis(df, expression, sets[[name]], plot_dir, name),
        tcga_gtex = .panfam_tcga_gtex_analysis(sets[[name]], sort(unique(df$Cancer)), plot_dir, name),
        paired = .panfam_paired_analysis(df),
        survival = .panfam_survival_analysis(df, plot_dir, name),
        tmb = .panfam_biomarker_analysis(df, "All_TMB", "TMB", plot_dir, correlation_method),
        msi = .panfam_biomarker_analysis(df, "All_MSI", "MSI", plot_dir, correlation_method),
        immune = .panfam_immune_analysis(df, expression, plot_dir, correlation_method),
        immunotherapy = .panfam_immunotherapy_analysis(sets[[name]], min_size, plot_dir, name),
        alterations = .panfam_alterations_analysis(sets[[name]], sort(unique(df$Cancer))),
        cnv_methylation = .panfam_cnv_methylation_analysis(sets[[name]], sort(unique(df$Cancer)), plot_dir)),
        error = function(e) e)
      key <- paste(name, module)
      if (inherits(value, "error")) {
        logs[[key]] <- data.frame(Family = name, Module = module, Status = "error", Rows = 0L, Message = conditionMessage(value))
        progress(paste0("[", name, "] ", module, " failed: ", conditionMessage(value)))
        next
      }
      summary <- .panfam_module_summary(value)
      logs[[key]] <- data.frame(Family = name, Module = module, Status = summary$status, Rows = summary$rows, Message = summary$message)
      for (table in names(value)) {
        if (ncol(value[[table]])) .panfam_write_csv(value[[table]], file.path(directory, paste0(table, ".csv")), row.names = FALSE, na = "")
      }
      tables[[name]] <- c(tables[[name]], value)
    }
  }
  result <- structure(list(scores = scores, tables = tables, gene_coverage = coverage,
                           run_log = .panfam_bind(logs), output_dir = root,
                           settings = list(created = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z", tz = "Asia/Shanghai"),
                                           data_source = if (bundled) "Bundled all_tpm (TCGA)" else "User-supplied expression",
                                           analyses = analyses, min_size = min_size, correlation_method = correlation_method,
                                           elapsed_seconds = unname(proc.time()["elapsed"] - started), gene_sets = sets,
                                           background_notes = attr(scores, "background_notes"),
                                           family_directories = .panfam_safe_names(names(sets)))), class = "panfam_result")
  .panfam_write_csv(result$run_log, file.path(root, "run_log.csv"), row.names = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(root, "sessionInfo.txt"))
  saveRDS(result, file.path(root, "analysis_result.rds"))
  .panfam_report(result)
  failures <- sum(result$run_log$Status == "error")
  progress(paste0("Finished with ", failures, " failed module(s). Report: ", file.path(root, "report.html")))
  if (failures) warning("Some requested modules failed. See result$run_log and report.html.", call. = FALSE)
  result
}

#' @export
print.panfam_result <- function(x, ...) {
  cat("PanFamR results:", length(x$scores), "family/families;", nrow(x$scores[[1L]]), "samples;",
      length(unique(x$scores[[1L]]$Cancer)), "cancer types\n")
  cat("Module status:", paste(names(table(x$run_log$Status)), as.integer(table(x$run_log$Status)), collapse = ", "), "\n")
  cat("Report:", file.path(x$output_dir, "report.html"), "\n")
  invisible(x)
}
