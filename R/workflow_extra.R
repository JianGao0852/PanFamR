.panfam_roc_analysis <- function(scores, directory, label) {
  rows <- curves <- list()
  for (cancer in sort(unique(scores$Cancer))) {
    data <- scores[scores$Cancer == cancer & scores$Group %in% c("Normal", "Tumor") &
                     is.finite(scores$Score), , drop = FALSE]
    n_normal <- sum(data$Group == "Normal")
    n_tumor <- sum(data$Group == "Tumor")
    row <- data.frame(Cancer = cancer, n_normal = n_normal, n_tumor = n_tumor,
                      AUC = NA_real_, CIlower = NA_real_, CIupper = NA_real_,
                      reason = "", stringsAsFactors = FALSE)
    if (min(n_normal, n_tumor) < 3L) {
      row$reason <- "Fewer than three finite scores in a group"
    } else if (length(unique(data$Score)) < 2L) {
      row$reason <- "Constant score"
    } else {
      fit <- tryCatch(pROC::roc(data$Group, data$Score, levels = c("Normal", "Tumor"),
                               direction = "<", quiet = TRUE), error = function(e) e)
      if (inherits(fit, "error")) {
        row$reason <- conditionMessage(fit)
      } else {
        row$AUC <- as.numeric(pROC::auc(fit))
        ci <- tryCatch(suppressWarnings(pROC::ci.auc(fit, method = "delong")), error = function(e) e)
        if (inherits(ci, "error")) {
          row$reason <- paste("AUC confidence interval unavailable:", conditionMessage(ci))
        } else {
          row$CIlower <- as.numeric(ci[1L])
          row$CIupper <- as.numeric(ci[3L])
        }
        curves[[cancer]] <- data.frame(Cancer = cancer, specificity = fit$specificities,
                                       sensitivity = fit$sensitivities, threshold = fit$thresholds,
                                       stringsAsFactors = FALSE)
        .panfam_original_roc_plot(fit, row, directory, label, cancer)
      }
    }
    rows[[cancer]] <- row
  }
  summary <- .panfam_bind(rows)
  curve <- .panfam_bind(curves)
  if (nrow(curve)) {
    plot_data <- curve
    plot_data$x <- 1 - plot_data$specificity
    plot_data$y <- plot_data$sensitivity
    plot <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey60") +
      ggplot2::geom_path(colour = "#56B4E9", linewidth = 0.6) +
      ggplot2::facet_wrap(~Cancer, ncol = min(6L, length(unique(curve$Cancer)))) +
      ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + ggplot2::theme_bw() +
      ggplot2::labs(title = paste(label, "tumor versus normal ROC"),
                    subtitle = "In-sample exploratory ROC; fixed direction: higher score predicts tumor",
                    x = "False positive rate (1 - specificity)", y = "Sensitivity")
    .panfam_save_plot(plot, directory, "tumor_normal_roc", width = 12,
                      height = max(4, min(16, ceiling(length(unique(curve$Cancer)) / 6) * 2.5)))
  }
  list(roc_summary = summary, roc_curve = curve)
}

.panfam_tumor_scores_analysis <- function(scores, directory, label) {
  tumors <- .panfam_tumors(scores)
  samples <- data.frame(Sample = rownames(tumors), tumors, row.names = NULL,
                        stringsAsFactors = FALSE)
  summary <- .panfam_bind(lapply(sort(unique(tumors$Cancer)), function(cancer) {
    values <- tumors$Score[tumors$Cancer == cancer & is.finite(tumors$Score)]
    data.frame(Cancer = cancer, n = length(values),
               mean = if (length(values)) mean(values) else NA_real_,
               median = if (length(values)) stats::median(values) else NA_real_,
               sd = if (length(values) > 1L) stats::sd(values) else NA_real_)
  }))
  if (nrow(tumors) && any(is.finite(tumors$Score))) {
    .panfam_original_tumor_plot(tumors, directory, label)
  }
  list(tumor_score_samples = samples, tumor_score_summary = summary)
}

.panfam_gene_survival_analysis <- function(scores, expression, genes, directory, label,
                                          endpoints = "OS") {
  endpoints <- unique(match.arg(endpoints, c("OS", "DSS", "PFI", "DFI"), several.ok = TRUE))
  tumors <- .panfam_tumors(scores)
  available_genes <- intersect(genes, names(expression))
  coverage <- data.frame(Gene = genes, Available = genes %in% available_genes,
                         reason = ifelse(genes %in% available_genes, "", "Gene absent from expression input"),
                         stringsAsFactors = FALSE)
  if (!length(available_genes)) {
    return(list(gene_survival_cox = data.frame(), gene_survival_km = data.frame(),
                gene_survival_gene_coverage = coverage))
  }
  os <- if ("OS" %in% endpoints) .panfam_aux_object("ALL_OS") else NULL
  other_endpoints <- if (any(endpoints != "OS")) .panfam_aux_object("DSS_PFI_DFI") else NULL
  expression_idx <- .panfam_match(rownames(tumors), rownames(expression))
  rows <- kms <- list()
  for (endpoint in endpoints) {
    clinical <- if (endpoint == "OS") os else other_endpoints
    clinical_idx <- .panfam_match(rownames(tumors), rownames(clinical))
    data <- data.frame(Sample = rownames(tumors), Cancer = tumors$Cancer,
                       Time = as.numeric(clinical[[paste0(endpoint, ".time")]][clinical_idx]),
                       Event = as.numeric(clinical[[endpoint]][clinical_idx]), stringsAsFactors = FALSE)
    for (gene in available_genes) {
      data$Score <- suppressWarnings(as.numeric(expression[[gene]][expression_idx]))
      for (cancer in sort(unique(data$Cancer))) {
        fit <- .panfam_survival_fit(data[data$Cancer == cancer, , drop = FALSE])
        key <- paste(endpoint, gene, cancer, sep = "_")
        rows[[key]] <- cbind(Gene = gene, Endpoint = endpoint, Cancer = cancer, fit$cox)
        kms[[key]] <- cbind(Gene = gene, Endpoint = endpoint, Cancer = cancer, fit$km)
      }
    }
  }
  cox <- .panfam_fdr(.panfam_bind(rows))
  km <- .panfam_fdr(.panfam_bind(kms))
  .panfam_original_gene_survival_heatmap(cox, directory, label)
  list(gene_survival_cox = cox, gene_survival_km = km, gene_survival_gene_coverage = coverage)
}
