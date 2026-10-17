# Plotting adapters for the user's original expression, score and ROC figures.
# These helpers receive the analysis data/results; plotting does not change
# sample groups, effect estimates, P values or the direction of an ROC curve.

.panfam_original_significance <- function(p) {
  out <- rep("", length(p))
  keep <- is.finite(p)
  out[keep] <- as.character(cut(p[keep],
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "ns")))
  out
}

.panfam_original_expression_heatmap <- function(table, directory, name, label) {
  required <- c("Gene", "Cancer", "mean_difference", "p_value")
  if (!all(required %in% names(table))) {
    stop("Expression plot table requires Gene, Cancer, mean_difference and p_value.", call. = FALSE)
  }
  if (!nrow(table)) return(invisible(NULL))
  values <- table
  if ("n_normal" %in% names(values)) {
    eligible_cancers <- unique(as.character(values$Cancer[is.finite(values$n_normal) & values$n_normal >= 3L]))
    values <- values[as.character(values$Cancer) %in% eligible_cancers, , drop = FALSE]
    if (!nrow(values)) return(invisible(NULL))
  }
  values$Gene <- factor(as.character(values$Gene), levels = unique(as.character(values$Gene)))
  values$Cancer <- factor(as.character(values$Cancer), levels = sort(unique(as.character(values$Cancer))))
  values$effect_text <- ifelse(is.finite(values$mean_difference), sprintf("%.2f", values$mean_difference), "")
  values$fill_bin <- as.character(cut(values$mean_difference,
    breaks = c(-Inf, -2, -1, 0, 1, 2, Inf),
    labels = c("<-2", "-2~-1", "-1~0", "0~1", "1~2", ">2")))
  values$fill_bin[!is.finite(values$p_value) | values$p_value >= 0.05 |
                    !is.finite(values$mean_difference)] <- "P=>0.05"
  values$fill_bin <- factor(values$fill_bin,
    levels = c(">2", "1~2", "0~1", "-1~0", "-2~-1", "<-2", "P=>0.05"))
  color_map <- c("<-2" = "#2f5d91", "-2~-1" = "#669dbb", "-1~0" = "#ded3cf",
    "0~1" = "#fdb1a4", "1~2" = "#f45859", ">2" = "#99253a", "P=>0.05" = "white")
  # The source called this mean(tumor) - mean(normal) statistic "logFC".
  # The package also accepts unlogged expression, so the legend names the
  # quantity actually calculated rather than assuming a logarithmic unit.
  plot <- ggplot2::ggplot(values, ggplot2::aes(x = Cancer, y = Gene)) +
    ggplot2::geom_tile(ggplot2::aes(fill = fill_bin), color = "grey90") +
    ggplot2::scale_fill_manual(values = color_map, name = "Mean difference", drop = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = effect_text), size = 3.2) +
    ggplot2::coord_fixed(ratio = 0.8) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(x = "Pan-Cancer", y = paste(label, "gene family")) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 10), panel.grid = ggplot2::element_blank())
  .panfam_save_plot(plot, directory, name, width = 13, height = 9)
  invisible(plot)
}

.panfam_original_score_plot <- function(scores, score_table, directory, label) {
  values <- scores[is.finite(scores$Score) & scores$Group %in% c("Normal", "Tumor"), , drop = FALSE]
  if (!nrow(values)) return(invisible(NULL))
  values$Cancer <- factor(as.character(values$Cancer), levels = sort(unique(as.character(values$Cancer))))
  values$Group <- factor(as.character(values$Group), levels = c("Normal", "Tumor"))
  annotations <- score_table[is.finite(score_table$p_value), c("Cancer", "p_value"), drop = FALSE]
  annotations$Cancer <- factor(as.character(annotations$Cancer), levels = levels(values$Cancer))
  annotations$label <- .panfam_original_significance(annotations$p_value)
  score_range <- range(values$Score)
  score_span <- diff(score_range)
  if (!is.finite(score_span) || score_span == 0) score_span <- max(1, abs(score_range))
  annotations$y <- score_range[2] + score_span * 0.04
  y_label <- paste(label, "ssGSEA score")
  plot <- ggpubr::ggboxplot(values, x = "Cancer", y = "Score", fill = "Group",
    xlab = "", ylab = y_label, color = "black", palette = "jco") +
    ggpubr::rotate_x_text(angle = 90) + ggpubr::grids(linetype = "dashed") +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 12), axis.title.x = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(size = 12), axis.title.y = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 12), legend.position = "right") +
    ggpubr::border("black") +
    ggplot2::geom_text(data = annotations, ggplot2::aes(x = Cancer, y = y, label = label),
      inherit.aes = FALSE, size = 3.8, vjust = 0)
  .panfam_save_plot(plot, directory, "score_distribution", width = 12, height = 3.6)

  # The second boxplot follows ssGsea_Score_diff1, including its Set1 palette.
  jitter_plot <- ggplot2::ggplot(values, ggplot2::aes(x = Cancer, y = Score, fill = Group, color = Group)) +
    ggplot2::geom_boxplot(alpha = 0.5, width = 0.6,
      position = ggplot2::position_dodge2(width = 0.6, preserve = "single", padding = 0.2)) +
    ggplot2::geom_point(shape = 20, size = 0.1,
      position = ggplot2::position_jitterdodge(jitter.width = 0.2, dodge.width = 0.6, seed = 0), alpha = 0.5) +
    ggplot2::theme_classic() +
    ggplot2::scale_fill_brewer(palette = "Set1", direction = -1) +
    ggplot2::scale_color_brewer(palette = "Set1", direction = -1) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(colour = "black", size = 11, angle = 45, hjust = 0.5, vjust = 0.5),
      axis.text.y = ggplot2::element_text(size = 12), axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(size = 14), legend.title = ggplot2::element_text(size = 14),
      legend.text = ggplot2::element_text(size = 12), legend.position = "right") +
    ggplot2::labs(y = y_label) +
    ggplot2::geom_text(data = annotations, ggplot2::aes(x = Cancer, y = y, label = label),
      inherit.aes = FALSE, color = "black", size = 3.8, vjust = 0)
  .panfam_save_plot(jitter_plot, directory, "score_distribution_jitter", width = 14, height = 3)
  invisible(list(boxplot = plot, jitter = jitter_plot))
}

.panfam_original_single_score_plots <- function(scores, score_table, directory, label) {
  for (cancer in sort(unique(as.character(scores$Cancer)))) {
    values <- scores[scores$Cancer == cancer & scores$Group %in% c("Normal", "Tumor") &
                       is.finite(scores$Score), , drop = FALSE]
    if (!all(c("Normal", "Tumor") %in% as.character(values$Group))) next
    values$Group <- factor(as.character(values$Group), levels = c("Normal", "Tumor"))
    plot <- ggpubr::ggboxplot(values, x = "Group", y = "Score", color = "Group",
                              add = "jitter", palette = "jco", xlab = cancer,
                              ylab = paste(label, "ssGSEA score")) +
      ggpubr::theme_classic2(base_size = 12) +
      ggplot2::theme(legend.position = "none", axis.title.x = ggplot2::element_text(size = 14),
                     axis.text.x = ggplot2::element_text(size = 12),
                     axis.title.y = ggplot2::element_text(size = 14),
                     axis.text.y = ggplot2::element_text(size = 12),
                     panel.border = ggplot2::element_rect(color = "black", fill = NA))
    row <- score_table[score_table$Cancer == cancer, , drop = FALSE]
    if (nrow(row) == 1L && is.finite(row$p_value)) {
      score_range <- range(values$Score)
      offset <- max(diff(score_range) * 0.08, abs(score_range[2L]) * 0.05, .Machine$double.eps)
      annotation <- data.frame(group1 = "Normal", group2 = "Tumor",
                                label = .panfam_original_significance(row$p_value),
                                y.position = score_range[2L] + offset)
      plot <- plot + ggpubr::stat_pvalue_manual(annotation, label = "label", tip.length = 0.01, size = 4)
    }
    .panfam_save_plot(plot, directory,
                      paste0("score_tumor_normal_", gsub("[^A-Za-z0-9_.-]", "_", cancer)),
                      width = 5, height = 4)
  }
  invisible(NULL)
}

.panfam_original_paired_plot <- function(pairs, directory, label, paired_table = NULL) {
  required <- c("Patient", "Cancer", "Score_tumor", "Score_normal")
  if (!all(required %in% names(pairs))) stop("Paired plot requires matched patient scores.", call. = FALSE)
  pairs <- pairs[is.finite(pairs$Score_tumor) & is.finite(pairs$Score_normal), , drop = FALSE]
  if (!nrow(pairs)) return(invisible(NULL))
  cancers <- sort(unique(as.character(pairs$Cancer)))
  normal <- data.frame(Patient = pairs$Patient, Cancer = pairs$Cancer, Group = "Normal", Score = pairs$Score_normal)
  tumor <- data.frame(Patient = pairs$Patient, Cancer = pairs$Cancer, Group = "Tumor", Score = pairs$Score_tumor)
  values <- rbind(normal, tumor)
  values$Cancer <- factor(as.character(values$Cancer), levels = cancers)
  values$Group <- factor(values$Group, levels = c("Normal", "Tumor"))
  # Include cancer in the pair ID so shared/custom patient IDs cannot connect
  # samples from different cancers. No groups or score values are reassigned.
  values$pair_id <- interaction(values$Cancer, values$Patient, drop = TRUE)
  if (is.null(paired_table)) {
    paired_table <- .panfam_bind(lapply(cancers, function(cancer) {
      subset <- pairs[pairs$Cancer == cancer, , drop = FALSE]
      p <- if (nrow(subset) >= 3L) suppressWarnings(tryCatch(
        stats::wilcox.test(subset$Score_tumor, subset$Score_normal, paired = TRUE, exact = FALSE)$p.value,
        error = function(e) NA_real_)) else NA_real_
      data.frame(Cancer = cancer, p_value = p)
    }))
  }
  annotations <- paired_table[is.finite(paired_table$p_value), c("Cancer", "p_value"), drop = FALSE]
  annotations$Cancer <- factor(as.character(annotations$Cancer), levels = cancers)
  annotations$group1 <- "Normal"
  annotations$group2 <- "Tumor"
  annotations$label <- .panfam_original_significance(annotations$p_value)
  maxima <- tapply(values$Score, values$Cancer, max)
  annotations$y.position <- as.numeric(maxima[as.character(annotations$Cancer)]) + 0.05
  palette <- c(Normal = "#4DBBD5", Tumor = "#E64B35")
  y_label <- paste(label, "ssGSEA score")
  plot <- ggpubr::ggpaired(values, x = "Group", y = "Score", id = "pair_id",
      color = "black", fill = "Group", line.size = 0.5, line.color = "grey80") +
    ggplot2::facet_wrap(~ Cancer, scales = "free_x", nrow = 2) +
    ggplot2::scale_fill_manual(values = palette) + ggplot2::scale_color_manual(values = palette) +
    ggpubr::theme_classic2(base_size = 9) + ggplot2::labs(x = "Pan-Cancer", y = y_label)
  if (nrow(annotations)) {
    plot <- plot + ggpubr::stat_pvalue_manual(annotations, label = "label", tip.length = 0.01, size = 2)
  }
  .panfam_save_plot(plot, directory, "paired_scores", width = 8, height = 5)
  single_plots <- vector("list", length(cancers))
  names(single_plots) <- cancers
  for (cancer in cancers) {
    subset <- values[values$Cancer == cancer, , drop = FALSE]
    plot_one <- ggpubr::ggpaired(subset, x = "Group", y = "Score", id = "pair_id",
        color = "black", fill = "Group", palette = palette, line.size = 0.5,
        line.color = "grey80", add = "jitter") +
      ggpubr::theme_classic2(base_size = 9) +
      ggplot2::theme(legend.position = "none",
        panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = ggplot2::element_line(color = "black"), axis.title.x = ggplot2::element_text(size = 14),
        axis.text.x = ggplot2::element_text(size = 12, colour = "black"),
        axis.title.y = ggplot2::element_text(size = 14), axis.text.y = ggplot2::element_text(size = 12, colour = "black")) +
      ggplot2::labs(x = cancer, y = y_label)
    annotation_one <- annotations[as.character(annotations$Cancer) == cancer, , drop = FALSE]
    if (nrow(annotation_one)) {
      plot_one <- plot_one + ggpubr::stat_pvalue_manual(annotation_one, label = "label", tip.length = 0.01, size = 3)
    }
    .panfam_save_plot(plot_one, directory, paste0("paired_scores_", gsub("[^A-Za-z0-9_.-]", "_", cancer)), width = 5, height = 4)
    single_plots[[cancer]] <- plot_one
  }
  invisible(list(combined = plot, by_cancer = single_plots))
}

.panfam_original_tumor_plot <- function(scores, directory, label) {
  values <- scores[scores$Group == "Tumor" & is.finite(scores$Score), , drop = FALSE]
  if (!nrow(values)) return(invisible(NULL))
  order <- names(sort(tapply(values$Score, values$Cancer, stats::median)))
  values$Cancer <- factor(as.character(values$Cancer), levels = order)
  plot <- ggplot2::ggplot(values, ggplot2::aes(x = Cancer, y = Score, fill = Cancer)) +
    ggplot2::geom_boxplot() + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1), legend.position = "none") +
    ggplot2::labs(y = paste(label, "ssGSEA score"), x = "Pan-Cancer")
  .panfam_save_plot(plot, directory, "tumor_score_distribution", width = 12, height = 3.6)
  invisible(plot)
}

.panfam_original_roc_plot <- function(fit, summary, directory, label, cancer) {
  # Use the ROC estimate already returned by analysis, retaining its direction,
  # smoothing choice and confidence interval rather than fitting a new curve.
  if (is.null(fit)) return(invisible(NULL))
  auc <- if ("AUC" %in% names(summary)) summary$AUC[1] else as.numeric(pROC::auc(fit))
  lower <- if ("CIlower" %in% names(summary)) summary$CIlower[1] else NA_real_
  upper <- if ("CIupper" %in% names(summary)) summary$CIupper[1] else NA_real_
  auc_label <- if (is.finite(auc)) sprintf("AUC: %.3f", auc) else "AUC: unavailable"
  ci_label <- if (is.finite(lower) && is.finite(upper)) sprintf("CI: %.3f - %.3f", lower, upper) else "CI: unavailable"
  plot <- pROC::ggroc(fit, color = "#56B4E9", linewidth = 1, legacy.axes = TRUE) +
    ggplot2::annotate("segment", x = 0, xend = 1, y = 0, yend = 1,
      color = "#505050", linetype = "dashed", linewidth = 0.75) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.title = ggplot2::element_blank(), axis.title.x = ggplot2::element_text(size = 14),
      axis.text.x = ggplot2::element_text(size = 12), axis.title.y = ggplot2::element_text(size = 14),
      axis.text.y = ggplot2::element_text(size = 12), panel.grid.major = ggplot2::element_line(colour = "lightgrey"),
      panel.grid.minor = ggplot2::element_line(colour = "lightgrey"), panel.background = ggplot2::element_rect(fill = "white")) +
    ggplot2::annotate("text", x = 0.98, y = 0.05, label = auc_label, size = 4, hjust = 1, fontface = "plain") +
    ggplot2::annotate("text", x = 0.98, y = 0, label = ci_label, size = 4, hjust = 1, fontface = "plain") +
    ggplot2::annotate("text", x = 0.98, y = 0.1,
      label = paste(strwrap(paste(label, "ssGSEA score"), width = 30), collapse = "\n"),
      size = 4, hjust = 1, vjust = 0, fontface = "plain") +
    ggplot2::labs(x = "False positive rate (1 - specificity)", y = "Sensitivity")
  .panfam_save_plot(plot, directory, paste0("tumor_normal_roc_", gsub("[^A-Za-z0-9_.-]", "_", cancer)), width = 6, height = 5)
  invisible(plot)
}
