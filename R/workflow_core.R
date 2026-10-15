.panfam_bind <- function(rows) {
  rows <- Filter(function(x) is.data.frame(x) && nrow(x) > 0L, rows)
  if (!length(rows)) return(data.frame())
  out <- do.call(rbind, unname(rows))
  rownames(out) <- NULL
  out
}

# Write UTF-8 directly: write.table() can translate non-native labels to
# <U+XXXX> before connection encoding is applied on an ASCII R session.
.panfam_write_csv <- function(x, file, row.names = FALSE, na = "") {
  if (!is.data.frame(x) || isTRUE(row.names)) stop("CSV output requires a data frame without row names.", call. = FALSE)
  quote_field <- function(value) paste0('"', gsub('"', '""', enc2utf8(value), fixed = TRUE), '"')
  columns <- lapply(x, function(value) {
    missing <- is.na(value)
    value <- if (is.numeric(value) && !inherits(value, c("Date", "POSIXt"))) {
      format(value, digits = 17L, trim = TRUE)
    } else as.character(value)
    value[missing] <- na
    quote_field(value)
  })
  header <- paste(quote_field(names(x)), collapse = ",")
  rows <- if (nrow(x)) do.call(paste, c(unname(columns), sep = ",")) else character()
  connection <- base::file(file, open = "wb")
  on.exit(close(connection), add = TRUE)
  writeLines(enc2utf8(c(header, rows)), connection, useBytes = TRUE)
  invisible(file)
}

.panfam_fdr <- function(x) {
  if (nrow(x) && "p_value" %in% names(x)) x$fdr <- stats::p.adjust(x$p_value, method = "BH")
  if (nrow(x) && "group_p_value" %in% names(x)) x$group_fdr <- stats::p.adjust(x$group_p_value, method = "BH")
  x
}

.panfam_match <- function(ids, reference) {
  idx <- match(ids, reference)
  tcga <- grepl("^TCGA-", ids)
  for (width in c(16L, 15L, 12L)) {
    keys <- substr(reference, 1L, width)
    ambiguous <- duplicated(keys) | duplicated(keys, fromLast = TRUE)
    available <- which(!ambiguous & grepl("^TCGA-", reference))
    missing <- which(is.na(idx) & tcga)
    if (length(missing) && length(available)) {
      idx[missing] <- available[match(substr(ids[missing], 1L, width), keys[available])]
    }
  }
  idx
}

.panfam_tumors <- function(scores) {
  out <- scores[scores$Group == "Tumor", , drop = FALSE]
  out <- out[order(rownames(out)), , drop = FALSE]
  ids <- rownames(out)
  keys <- ifelse(grepl("^TCGA-", ids), substr(ids, 1L, 12L), ids)
  out[!duplicated(keys), , drop = FALSE]
}

.panfam_save_plot <- function(plot, directory, name, width = 11, height = 5) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(file.path(directory, paste0(name, ".pdf")), plot, width = width, height = height, device = pdf_device, bg = "white")
  ggplot2::ggsave(file.path(directory, paste0(name, ".png")), plot, width = width, height = height, dpi = 300, bg = "white")
}

.panfam_difference <- function(tumor, normal) {
  tumor <- tumor[is.finite(tumor)]
  normal <- normal[is.finite(normal)]
  valid <- length(tumor) >= 3L && length(normal) >= 3L
  fit <- if (valid) suppressWarnings(tryCatch(stats::wilcox.test(tumor, normal), error = function(e) NULL)) else NULL
  data.frame(n_tumor = length(tumor), n_normal = length(normal),
             mean_difference = if (length(tumor) && length(normal)) mean(tumor) - mean(normal) else NA_real_,
             median_difference = if (length(tumor) && length(normal)) stats::median(tumor) - stats::median(normal) else NA_real_,
             rank_biserial = if (!is.null(fit)) 2 * as.numeric(fit$statistic) / (length(tumor) * length(normal)) - 1 else NA_real_,
             p_value = if (!is.null(fit)) fit$p.value else NA_real_,
             reason = if (valid) "" else "Fewer than three samples in a group", stringsAsFactors = FALSE)
}

.panfam_expression_analysis <- function(scores, expression, genes, directory, label) {
  cancers <- sort(unique(scores$Cancer))
  score_table <- .panfam_fdr(.panfam_bind(lapply(cancers, function(cancer) {
    tmp <- scores[scores$Cancer == cancer, , drop = FALSE]
    cbind(Cancer = cancer, .panfam_difference(tmp$Score[tmp$Group == "Tumor"], tmp$Score[tmp$Group == "Normal"]))
  })))
  genes <- intersect(genes, names(expression))
  gene_table <- .panfam_fdr(.panfam_bind(lapply(genes, function(gene) {
    .panfam_bind(lapply(cancers, function(cancer) {
      keep <- expression$Cancer == cancer
      values <- suppressWarnings(as.numeric(expression[[gene]][keep]))
      groups <- expression$Group[keep]
      cbind(Gene = gene, Cancer = cancer, .panfam_difference(values[groups == "Tumor"], values[groups == "Normal"]))
    }))
  })))
  .panfam_original_score_plot(scores, score_table, directory, label)
  if (nrow(gene_table)) {
    .panfam_original_expression_heatmap(gene_table, directory, "gene_expression_differences", label)
  }
  list(score_differences = score_table, gene_differences = gene_table)
}

.panfam_paired_analysis <- function(scores, directory = NULL, label = "Gene family") {
  ids <- rownames(scores)
  keep <- grepl("^TCGA-", ids)
  tmp <- scores[keep, , drop = FALSE]
  tmp$Patient <- substr(rownames(tmp), 1L, 12L)
  tumor <- tmp[tmp$Group == "Tumor", , drop = FALSE]
  normal <- tmp[tmp$Group == "Normal", , drop = FALSE]
  tumor <- tumor[order(rownames(tumor)), , drop = FALSE]
  normal <- normal[order(rownames(normal)), , drop = FALSE]
  tumor <- tumor[!duplicated(tumor$Patient), , drop = FALSE]
  normal <- normal[!duplicated(normal$Patient), , drop = FALSE]
  pairs <- merge(tumor[, c("Patient", "Cancer", "Score")], normal[, c("Patient", "Cancer", "Score")],
                 by = c("Patient", "Cancer"), suffixes = c("_tumor", "_normal"))
  summary <- .panfam_fdr(.panfam_bind(lapply(sort(unique(scores$Cancer)), function(cancer) {
    x <- pairs[pairs$Cancer == cancer, , drop = FALSE]
    p <- if (nrow(x) >= 3) suppressWarnings(tryCatch(stats::wilcox.test(x$Score_tumor, x$Score_normal, paired = TRUE)$p.value,
                                                     error = function(e) NA_real_)) else NA_real_
    data.frame(Cancer = cancer, n_pairs = nrow(x),
               median_difference = if (nrow(x)) stats::median(x$Score_tumor - x$Score_normal) else NA_real_, p_value = p,
               reason = if (nrow(x) < 3) "Fewer than three matched pairs" else if (!is.finite(p)) "Paired test could not be estimated" else "")
  })))
  if (!is.null(directory) && nrow(pairs)) {
    .panfam_original_paired_plot(pairs, directory, label, summary)
  }
  list(paired_samples = pairs, paired_score_differences = summary)
}

.panfam_correlation <- function(x, y, method = "pearson") {
  valid <- is.finite(x) & is.finite(y)
  x <- x[valid]
  y <- y[valid]
  enough <- length(x) >= 3L && stats::sd(x) > 0 && stats::sd(y) > 0
  fit <- if (enough) suppressWarnings(tryCatch(stats::cor.test(x, y, method = method, exact = FALSE), error = function(e) NULL)) else NULL
  data.frame(n = length(x), correlation = if (!is.null(fit)) unname(fit$estimate) else NA_real_,
             p_value = if (!is.null(fit)) fit$p.value else NA_real_,
             reason = if (!is.null(fit)) "" else "Insufficient matched samples or constant values", stringsAsFactors = FALSE)
}

.panfam_correlations <- function(scores, features, panel, method) {
  scores <- .panfam_tumors(scores)
  features <- as.data.frame(features)
  idx <- .panfam_match(rownames(scores), rownames(features))
  .panfam_fdr(.panfam_bind(lapply(sort(unique(scores$Cancer)), function(cancer) {
    rows <- which(scores$Cancer == cancer)
    .panfam_bind(lapply(names(features), function(feature) {
      values <- suppressWarnings(as.numeric(features[[feature]][idx[rows]]))
      cbind(Panel = panel, Cancer = cancer, Feature = feature,
            .panfam_correlation(scores$Score[rows], values, method))
    }))
  })))
}

.panfam_correlation_plot <- function(table, directory, name, title) {
  .panfam_original_correlation_plot(table, directory, name, title)
}

.panfam_aux_object <- function(name) {
  values <- panfam_load_data(name)
  if (length(values) != 1L) stop("Bundled dataset must contain exactly one table: ", name, call. = FALSE)
  values[[1L]]
}

.panfam_biomarker_analysis <- function(scores, dataset, panel, directory, method) {
  table <- .panfam_correlations(scores, .panfam_aux_object(dataset), panel, method)
  attr(table, "correlation_method") <- method
  .panfam_correlation_plot(table, directory, paste0(tolower(panel), "_correlations"), paste(panel, "association"))
  stats::setNames(list(table), paste0(tolower(panel), "_correlations"))
}

.panfam_survival_fit <- function(data) {
  data <- data[is.finite(data$Time) & data$Time > 0 & data$Event %in% c(0, 1) & is.finite(data$Score), , drop = FALSE]
  events <- sum(data$Event)
  row <- data.frame(n = nrow(data), events = events, hr_per_sd = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                    p_value = NA_real_, ph_p_value = NA_real_,
                    hr_per_unit = NA_real_, unit_ci_lower = NA_real_, unit_ci_upper = NA_real_, unit_p_value = NA_real_,
                    hr_high_vs_low = NA_real_, group_ci_lower = NA_real_, group_ci_upper = NA_real_,
                    group_p_value = NA_real_, group_reason = "", reason = "", stringsAsFactors = FALSE)
  km <- data.frame(n = nrow(data), events = events, n_low = 0L, n_high = 0L, cutoff = NA_real_, p_value = NA_real_, reason = "")
  if (nrow(data) < 10L || events < 5L || !is.finite(stats::sd(data$Score)) || stats::sd(data$Score) <= 0) {
    row$reason <- km$reason <- "Fewer than 10 complete samples, fewer than 5 events, or constant score"
    return(list(cox = row, km = km, data = data, fit = NULL))
  }
  data$ScaledScore <- as.numeric(scale(data$Score))
  problems <- character()
  fit <- tryCatch(withCallingHandlers(survival::coxph(survival::Surv(Time, Event) ~ ScaledScore, data = data, ties = "efron"),
                                      warning = function(w) { problems <<- c(problems, conditionMessage(w)); invokeRestart("muffleWarning") }), error = function(e) NULL)
  if (!is.null(fit) && !length(problems)) {
    sm <- summary(fit)
    ci <- sm$conf.int[1, ]
    if (all(is.finite(ci))) {
      row$hr_per_sd <- ci["exp(coef)"]
      row$ci_lower <- ci["lower .95"]
      row$ci_upper <- ci["upper .95"]
      row$p_value <- sm$coefficients[1, "Pr(>|z|)"]
      row$ph_p_value <- tryCatch(survival::cox.zph(fit)$table[1, "p"], error = function(e) NA_real_)
      score_sd <- stats::sd(data$Score)
      row$hr_per_unit <- exp(log(row$hr_per_sd) / score_sd)
      row$unit_ci_lower <- exp(log(row$ci_lower) / score_sd)
      row$unit_ci_upper <- exp(log(row$ci_upper) / score_sd)
      row$unit_p_value <- row$p_value
    }
  } else row$reason <- if (length(problems)) paste(problems, collapse = "; ") else "Cox model could not be fitted"
  km$cutoff <- stats::median(data$Score)
  data$ScoreGroup <- factor(ifelse(data$Score > km$cutoff, "High", "Low"), levels = c("Low", "High"))
  km$n_low <- sum(data$ScoreGroup == "Low")
  km$n_high <- sum(data$ScoreGroup == "High")
  kmfit <- NULL
  if (min(km$n_low, km$n_high) >= 3L) {
    diff <- survival::survdiff(survival::Surv(Time, Event) ~ ScoreGroup, data = data)
    km$p_value <- stats::pchisq(diff$chisq, df = 1, lower.tail = FALSE)
    kmfit <- survival::survfit(survival::Surv(Time, Event) ~ ScoreGroup, data = data)
    group_problems <- character()
    group_fit <- tryCatch(withCallingHandlers(
      survival::coxph(survival::Surv(Time, Event) ~ ScoreGroup, data = data, ties = "efron"),
      warning = function(w) { group_problems <<- c(group_problems, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) { group_problems <<- c(group_problems, conditionMessage(e)); NULL })
    if (!is.null(group_fit) && !length(group_problems)) {
      group_summary <- summary(group_fit)
      group_ci <- group_summary$conf.int[1, ]
      if (all(is.finite(group_ci))) {
        row$hr_high_vs_low <- group_ci["exp(coef)"]
        row$group_ci_lower <- group_ci["lower .95"]
        row$group_ci_upper <- group_ci["upper .95"]
        row$group_p_value <- group_summary$coefficients[1, "Pr(>|z|)"]
      }
    } else row$group_reason <- if (length(group_problems)) paste(group_problems, collapse = "; ") else "Grouped Cox model could not be fitted"
  } else km$reason <- "Fewer than three samples in a median-split group"
  list(cox = row, km = km, data = data, fit = kmfit)
}

.panfam_km_plot <- function(result, path, title, time_label = "Time (days)") {
  .panfam_original_km_plot(result, path, title, time_label)
}

.panfam_survival_analysis <- function(scores, directory, label, endpoints = c("OS", "DSS", "PFI", "DFI")) {
  endpoints <- unique(match.arg(endpoints, c("OS", "DSS", "PFI", "DFI"), several.ok = TRUE))
  tumors <- .panfam_tumors(scores)
  os <- if ("OS" %in% endpoints) .panfam_aux_object("ALL_OS") else NULL
  other_endpoints <- if (any(endpoints != "OS")) .panfam_aux_object("DSS_PFI_DFI") else NULL
  rows <- kms <- samples <- list()
  for (endpoint in endpoints) {
    clinical <- if (endpoint == "OS") os else other_endpoints
    idx <- .panfam_match(rownames(tumors), rownames(clinical))
    data <- data.frame(Sample = rownames(tumors), Cancer = tumors$Cancer,
                       Endpoint = endpoint, Score = tumors$Score,
                       Time = as.numeric(clinical[[paste0(endpoint, ".time")]][idx]),
                       Event = as.numeric(clinical[[endpoint]][idx]), stringsAsFactors = FALSE)
    for (cancer in sort(unique(data$Cancer))) {
      fit <- .panfam_survival_fit(data[data$Cancer == cancer, , drop = FALSE])
      key <- paste(endpoint, cancer, sep = "_")
      rows[[key]] <- cbind(Endpoint = endpoint, Cancer = cancer, fit$cox)
      kms[[key]] <- cbind(Endpoint = endpoint, Cancer = cancer, fit$km)
      if (nrow(fit$data)) samples[[key]] <- fit$data[, intersect(c("Sample", "Cancer", "Endpoint", "Score", "Time", "Event", "ScoreGroup"), names(fit$data)), drop = FALSE]
      .panfam_km_plot(fit, file.path(directory, "survival_KM", paste0(key, ".pdf")), paste(label, endpoint, cancer))
    }
  }
  cox <- .panfam_fdr(.panfam_bind(rows))
  km <- .panfam_fdr(.panfam_bind(kms))
  .panfam_original_survival_forest(cox, directory, label)
  # Some small cancers may lack a valid score group; keep a uniform schema.
  samples <- lapply(samples, function(x) { if (!"ScoreGroup" %in% names(x)) x$ScoreGroup <- NA_character_; x })
  list(survival_cox = cox, survival_km = km, survival_samples = .panfam_bind(samples))
}
