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
  ggplot2::ggsave(file.path(directory, paste0(name, ".pdf")), plot, width = width, height = height, device = pdf_device)
  ggplot2::ggsave(file.path(directory, paste0(name, ".png")), plot, width = width, height = height, dpi = 120)
}

.panfam_difference <- function(tumor, normal) {
  tumor <- tumor[is.finite(tumor)]
  normal <- normal[is.finite(normal)]
  valid <- length(tumor) >= 3L && length(normal) >= 3L
  fit <- if (valid) tryCatch(stats::wilcox.test(tumor, normal, exact = FALSE), error = function(e) NULL) else NULL
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
  plot <- ggplot2::ggplot(scores, ggplot2::aes(x = Cancer, y = Score, fill = Group)) +
    ggplot2::geom_boxplot(outlier.size = 0.3) + ggplot2::theme_bw() +
    ggplot2::scale_fill_manual(values = c(Normal = "#397FB6", Tumor = "#D65A45")) +
    ggplot2::labs(title = paste(label, "ssGSEA scores"), x = NULL, y = "ssGSEA score") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  .panfam_save_plot(plot, directory, "score_distribution")
  if (nrow(gene_table)) {
    plot <- ggplot2::ggplot(gene_table, ggplot2::aes(x = Cancer, y = Gene, fill = rank_biserial)) +
      ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "#397FB6", mid = "white", high = "#D65A45", limits = c(-1, 1), na.value = "grey90") +
      ggplot2::theme_minimal() + ggplot2::labs(title = paste(label, "tumor versus normal"), fill = "Rank effect", x = NULL, y = NULL) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    .panfam_save_plot(plot, directory, "gene_expression_differences", height = max(4, min(14, length(genes) * 0.22 + 2)))
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
    p <- if (nrow(x) >= 3) suppressWarnings(tryCatch(stats::wilcox.test(x$Score_tumor, x$Score_normal, paired = TRUE, exact = FALSE)$p.value,
                                                     error = function(e) NA_real_)) else NA_real_
    data.frame(Cancer = cancer, n_pairs = nrow(x),
               median_difference = if (nrow(x)) stats::median(x$Score_tumor - x$Score_normal) else NA_real_, p_value = p,
               reason = if (nrow(x) < 3) "Fewer than three matched pairs" else if (!is.finite(p)) "Paired test could not be estimated" else "")
  })))
  if (!is.null(directory) && nrow(pairs)) {
    normal_plot <- data.frame(Patient = pairs$Patient, Cancer = pairs$Cancer, Group = "Normal", Score = pairs$Score_normal)
    tumor_plot <- data.frame(Patient = pairs$Patient, Cancer = pairs$Cancer, Group = "Tumor", Score = pairs$Score_tumor)
    values <- rbind(normal_plot, tumor_plot)
    values$Group <- factor(values$Group, levels = c("Normal", "Tumor"))
    plot <- ggplot2::ggplot(values, ggplot2::aes(x = Group, y = Score, group = Patient)) +
      ggplot2::geom_line(colour = "grey75", alpha = 0.45) +
      ggplot2::geom_point(ggplot2::aes(colour = Group), size = 0.6, alpha = 0.7) +
      ggplot2::facet_wrap(~Cancer, scales = "free_y", ncol = 6) +
      ggplot2::scale_colour_manual(values = c(Normal = "#397FB6", Tumor = "#D65A45")) +
      ggplot2::theme_bw() + ggplot2::labs(title = paste(label, "matched tumor-normal scores"), x = NULL, y = "ssGSEA score")
    .panfam_save_plot(plot, directory, "paired_scores", height = max(4, ceiling(length(unique(pairs$Cancer)) / 6) * 2))
  }
  list(paired_samples = pairs, paired_score_differences = summary)
}

.panfam_correlation <- function(x, y, method = "spearman") {
  valid <- is.finite(x) & is.finite(y)
  x <- x[valid]
  y <- y[valid]
  enough <- length(x) >= 5L && stats::sd(x) > 0 && stats::sd(y) > 0
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
  if (!nrow(table) || !any(is.finite(table$correlation))) return(invisible(NULL))
  table$mark <- ifelse(!is.na(table$fdr) & table$fdr < 0.05, "*", "")
  plot <- ggplot2::ggplot(table, ggplot2::aes(x = Cancer, y = Feature, fill = correlation)) +
    ggplot2::geom_tile() + ggplot2::geom_text(ggplot2::aes(label = mark), size = 3) +
    ggplot2::scale_fill_gradient2(low = "#397FB6", mid = "white", high = "#D65A45", limits = c(-1, 1), na.value = "grey90") +
    ggplot2::theme_minimal() + ggplot2::labs(title = title, subtitle = "* BH-adjusted P < 0.05", fill = "Correlation", x = NULL, y = NULL) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  .panfam_save_plot(plot, directory, name, height = max(3, min(16, length(unique(table$Feature)) * 0.22 + 2)))
}

.panfam_aux_object <- function(name) {
  values <- panfam_load_data(name)
  if (length(values) != 1L) stop("Bundled dataset must contain exactly one table: ", name, call. = FALSE)
  values[[1L]]
}

.panfam_biomarker_analysis <- function(scores, dataset, panel, directory, method) {
  table <- .panfam_correlations(scores, .panfam_aux_object(dataset), panel, method)
  .panfam_correlation_plot(table, directory, paste0(tolower(panel), "_correlations"), paste(panel, "association"))
  stats::setNames(list(table), paste0(tolower(panel), "_correlations"))
}

.panfam_survival_fit <- function(data) {
  data <- data[is.finite(data$Time) & data$Time > 0 & data$Event %in% c(0, 1) & is.finite(data$Score), , drop = FALSE]
  events <- sum(data$Event)
  row <- data.frame(n = nrow(data), events = events, hr_per_sd = NA_real_, ci_lower = NA_real_, ci_upper = NA_real_,
                    p_value = NA_real_, ph_p_value = NA_real_, reason = "", stringsAsFactors = FALSE)
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
  } else km$reason <- "Fewer than three samples in a median-split group"
  list(cox = row, km = km, data = data, fit = kmfit)
}

.panfam_km_plot <- function(result, path, title, time_label = "Time (days)") {
  if (is.null(result$fit)) return(invisible(NULL))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (capabilities("cairo")) grDevices::cairo_pdf(path, width = 7, height = 6) else grDevices::pdf(path, width = 7, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::plot(result$fit, col = c("#397FB6", "#D65A45"), lwd = 2, xlab = time_label,
                 ylab = "Survival probability", main = title, conf.int = FALSE, mark.time = TRUE)
  graphics::legend("bottomleft", legend = paste0(c("Low", "High"), " (n=", c(result$km$n_low, result$km$n_high), ")"),
                   col = c("#397FB6", "#D65A45"), lty = 1, lwd = 2, bty = "n")
  graphics::mtext(paste0("Median split; log-rank P = ", format.pval(result$km$p_value, digits = 3)), side = 3, line = 0.2, cex = 0.8)
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
  valid <- cox[is.finite(cox$hr_per_sd) & is.finite(cox$ci_lower) & is.finite(cox$ci_upper), , drop = FALSE]
  if (nrow(valid)) {
    plot <- ggplot2::ggplot(valid, ggplot2::aes(x = hr_per_sd, y = Cancer)) +
      ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey60") +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
      ggplot2::geom_point(ggplot2::aes(colour = fdr < 0.05)) + ggplot2::scale_x_log10() +
      ggplot2::facet_wrap(~Endpoint, nrow = 1) + ggplot2::theme_bw() +
      ggplot2::labs(title = paste(label, "univariate survival association"), x = "Hazard ratio per score SD (95% CI)", y = NULL, colour = "FDR < 0.05")
    .panfam_save_plot(plot, directory, "survival_forest", width = 13, height = 11)
  }
  # Some small cancers may lack a valid score group; keep a uniform schema.
  samples <- lapply(samples, function(x) { if (!"ScoreGroup" %in% names(x)) x$ScoreGroup <- NA_character_; x })
  list(survival_cox = cox, survival_km = km, survival_samples = .panfam_bind(samples))
}
