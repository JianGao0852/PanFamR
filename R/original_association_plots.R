# Plot adapters for the original analysis_functions.R visual designs.
# These adapters consume the modular result tables without changing their
# estimates, raw P values, BH FDR values, or unavailable-data records.

.panfam_original_assoc_save <- function(draw, directory, name, width, height) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  paths <- file.path(directory, paste0(name, c(".pdf", ".png")))
  for (i in seq_along(paths)) {
    if (i == 1L) {
      if (capabilities("cairo")) {
        grDevices::cairo_pdf(paths[[i]], width = width, height = height, bg = "white")
      } else {
        grDevices::pdf(paths[[i]], width = width, height = height, bg = "white")
      }
    } else {
      arguments <- list(filename = paths[[i]], width = width, height = height,
                        units = "in", res = 300, bg = "white")
      if (capabilities("cairo")) arguments$type <- "cairo"
      do.call(grDevices::png, arguments)
    }
    tryCatch(draw(), finally = grDevices::dev.off())
  }
  invisible(paths)
}

.panfam_original_assoc_stars <- function(p, radar = FALSE) {
  stars <- rep("", length(p))
  valid <- is.finite(p) & p >= 0 & p <= 1
  stars[valid & p < 0.05] <- "*"
  stars[valid & p < 0.01] <- "**"
  if (!radar) stars[valid & p < 0.001] <- "***"
  stars
}

.panfam_original_assoc_method <- function(table) {
  method <- attr(table, "correlation_method", exact = TRUE)
  if (is.null(method) && "method" %in% names(table)) method <- unique(table$method)
  method <- unique(tolower(as.character(method)))
  method <- method[!is.na(method) & nzchar(method)]
  if (length(method) == 1L && method %in% c("pearson", "spearman", "kendall")) {
    paste0(toupper(substr(method, 1L, 1L)), substring(method, 2L))
  } else "Correlation"
}

.panfam_original_radar_plot <- function(table, directory, name, title) {
  # Sources: ssGsea_Score_pan_cancer_TMB / ssGsea_Score_pan_cancer_MSI.
  # Original significance labels use unadjusted P, not the retained FDR column.
  panel <- toupper(as.character(table$Panel[[1L]]))
  colour <- if (panel == "MSI") "#00AFBB" else "#E64B35"
  valid <- is.finite(table$correlation) & !is.na(table$Cancer)
  available <- table[valid & !duplicated(table$Cancer), , drop = FALSE]
  if (!nrow(available)) return(invisible(NULL))
  omitted <- sum(!valid)
  note <- paste0("Raw P: * < 0.05; ** < 0.01. BH FDR retained in CSV.",
                 if (omitted) paste0(" Unavailable coefficients omitted: ", omitted, ".") else "")
  draw <- function() {
    old <- graphics::par(mar = c(3.5, 3, 3, 3))
    on.exit(graphics::par(old), add = TRUE)
    if (nrow(available) < 3L) {
      # fmsb cannot draw a radar with fewer than three axes. Do not invent or
      # repeat cancer types to make up a polygon.
      graphics::plot(available$correlation, seq_len(nrow(available)), pch = 19,
                     col = colour, xlim = c(-1, 1), yaxt = "n", ylab = "",
                     xlab = .panfam_original_assoc_method(table), main = title,
                     ylim = c(0.5, nrow(available) + 0.5))
      graphics::abline(v = 0, col = "grey80", lty = 2)
      graphics::axis(2, at = seq_len(nrow(available)), las = 1,
                     labels = paste0(available$Cancer,
                                     .panfam_original_assoc_stars(available$p_value, TRUE)))
      graphics::mtext("Radar requires at least three evaluable cancer types.",
                      side = 3, line = 0.15, cex = 0.75)
    } else {
      rho <- available$correlation
      limits <- range(rho)
      if (diff(limits) < 1e-8) {
        limits <- c(max(-1, limits[[1L]] - 0.05), min(1, limits[[2L]] + 0.05))
      }
      matrix <- rbind(Max = rep(limits[[2L]], length(rho)),
                      Min = rep(limits[[1L]], length(rho)), rho = rho)
      colnames(matrix) <- paste0(available$Cancer,
                                 .panfam_original_assoc_stars(available$p_value, TRUE))
      # Preserve the original clockwise axis order.
      matrix <- matrix[, c(1L, seq.int(ncol(matrix), 2L)), drop = FALSE]
      radar <- as.data.frame(matrix, check.names = FALSE)
      fmsb::radarchart(
        radar, axistype = 1, pcol = colour,
        pfcol = grDevices::adjustcolor(colour, alpha.f = 0.1),
        plwd = 2, plty = 1, cglcol = "grey", cglty = 1, cglwd = 0.8,
        axislabcol = "black", title = title,
        vlcex = if (ncol(radar) > 33L) 0.8 else 1,
        vlabels = colnames(radar),
        caxislabels = format(round(seq(limits[[1L]], limits[[2L]], length.out = 5L), 2),
                             trim = TRUE)
      )
    }
    graphics::mtext(paste(strwrap(note, width = 100), collapse = "\n"),
                    side = 1, line = 2, cex = 0.62)
  }
  .panfam_original_assoc_save(draw, directory, name, width = 7.5, height = 7.5)
}

.panfam_original_heatmap_plot <- function(table, directory, name, title) {
  # Sources: original ICG, chemokine/receptor, immune stimulator/inhibitor,
  # immunecellratio, immunescore, TME_GSVA and HALLMARK_GSVA functions.
  cancers <- unique(as.character(table$Cancer[!is.na(table$Cancer)]))
  features <- unique(as.character(table$Feature[!is.na(table$Feature)]))
  if (!length(cancers) || !length(features)) return(invisible(NULL))
  values <- matrix(NA_real_, nrow = length(features), ncol = length(cancers),
                   dimnames = list(features, cancers))
  marks <- matrix("", nrow = length(features), ncol = length(cancers),
                  dimnames = list(features, cancers))
  rows <- match(as.character(table$Feature), features)
  columns <- match(as.character(table$Cancer), cancers)
  keep <- !is.na(rows) & !is.na(columns)
  positions <- cbind(rows[keep], columns[keep])
  values[positions] <- table$correlation[keep]
  marks[positions] <- .panfam_original_assoc_stars(table$p_value[keep])
  marks[!is.finite(values)] <- ""
  finite <- values[is.finite(values)]
  if (!length(finite)) return(invisible(NULL))
  # Keep the original observed-range blue/white/red scale, with a strictly
  # increasing range when all observed values fall on one side of zero.
  low <- min(c(finite, -1e-6))
  high <- max(c(finite, 1e-6))
  col_fun <- circlize::colorRamp2(c(low, 0, high), c("blue", "white", "red"))
  panel <- as.character(table$Panel[[1L]])
  is_hallmark <- grepl("hallmark", panel, ignore.case = TRUE)
  is_tme <- grepl("tme", panel, ignore.case = TRUE)
  baseline_height <- if (is_hallmark) 10 else if (is_tme) 4.5 else 2
  width <- max(if (is_hallmark || is_tme) 10 else 8,
               min(20, length(cancers) * 0.22 + max(nchar(features)) * 0.055 + 1.8))
  height <- max(baseline_height, length(features) * 0.16 + 1.5)
  display_names <- if (is_hallmark) sub("^HALLMARK_", "", features) else features
  draw <- function() {
    heatmap <- ComplexHeatmap::Heatmap(
      values, name = .panfam_original_assoc_method(table), col = col_fun,
      na_col = "grey90", cluster_rows = FALSE, cluster_columns = FALSE,
      column_title = paste0(title, "\nRaw P: * < 0.05; ** < 0.01; *** < 0.001. BH FDR in CSV."),
      column_title_gp = grid::gpar(fontsize = 10),
      row_labels = display_names, row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = if (is_hallmark) 8.4 else 10),
      row_names_max_width = grid::unit(max(6, max(nchar(display_names)) * 0.15), "cm"),
      column_names_side = "bottom", column_names_rot = 90,
      column_names_gp = grid::gpar(fontsize = if (is_hallmark) 10 else 12),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid::grid.text(marks[i, j], x, y, gp = grid::gpar(fontsize = 8))
      },
      heatmap_legend_param = list(title = .panfam_original_assoc_method(table),
                                  title_gp = grid::gpar(fontsize = 8),
                                  labels_gp = grid::gpar(fontsize = 8),
                                  legend_height = grid::unit(2, "cm")),
      rect_gp = grid::gpar(col = "white")
    )
    ComplexHeatmap::draw(heatmap)
  }
  .panfam_original_assoc_save(draw, directory, name, width, height)
}

.panfam_original_correlation_plot <- function(table, directory, name, title) {
  if (!nrow(table) || !any(is.finite(table$correlation))) return(invisible(NULL))
  panels <- unique(as.character(table$Panel[!is.na(table$Panel)]))
  paths <- list()
  for (panel in panels) {
    data <- table[!is.na(table$Panel) & table$Panel == panel, , drop = FALSE]
    attr(data, "correlation_method") <- attr(table, "correlation_method", exact = TRUE)
    output_name <- if (length(panels) == 1L) name else {
      paste0(name, "_", gsub("[^[:alnum:]_-]", "_", tolower(panel)))
    }
    heading <- if (length(panels) == 1L) title else paste(title, panel)
    if (toupper(panel) %in% c("TMB", "MSI")) {
      paths[[panel]] <- .panfam_original_radar_plot(data, directory, output_name, heading)
    } else {
      paths[[panel]] <- .panfam_original_heatmap_plot(data, directory, output_name, heading)
    }
  }
  invisible(paths)
}

.panfam_original_cnv_methylation_plot <- function(table, directory, name, label) {
  if (!nrow(table) || !any(is.finite(table$correlation))) return(invisible(NULL))
  associations <- unique(as.character(table$Association[!is.na(table$Association)]))
  paths <- list()
  for (association in associations) {
    data <- table[!is.na(table$Association) & table$Association == association, , drop = FALSE]
    data <- data[is.finite(data$correlation) & is.finite(data$fdr) &
                   data$fdr >= 0 & data$fdr <= 1, , drop = FALSE]
    if (!nrow(data)) next
    data$Cancer <- factor(data$Cancer, levels = sort(unique(as.character(data$Cancer))))
    data$Gene <- factor(data$Gene, levels = sort(unique(as.character(data$Gene))))
    data$sig <- factor(ifelse(data$fdr <= 0.05, "<=0.05", ">0.05"),
                        levels = c("<=0.05", ">0.05"))
    title <- if (association == "CNV") "Correlations of CNV with mRNA expression" else {
      "Correlation between methylation and\nmRNA expression"
    }
    output_name <- if (length(associations) == 1L) name else paste0(name, "_", tolower(association))
    width <- max(10, min(20, length(levels(data$Cancer)) * 0.25 + 2))
    height <- max(5, length(levels(data$Gene)) * 0.24 + 1.8)
    theme <- ggplot2::theme_minimal(base_size = 10) + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
      axis.text.y = ggplot2::element_text(size = 12),
      panel.grid.major = ggplot2::element_line(colour = "grey80", linetype = "dashed"),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 1.2),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 8)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 8)),
      legend.key.size = grid::unit(0.6, "lines"),
      legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_text(size = 8),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.5),
      axis.ticks.length = grid::unit(0.1, "cm"),
      plot.caption = ggplot2::element_text(size = 8, hjust = 0)
    )
    make_plot <- function(frame) {
      ggplot2::ggplot(frame, ggplot2::aes(x = Cancer, y = Gene)) +
        ggplot2::geom_point(ggplot2::aes(size = point_size, fill = correlation, colour = sig),
                            shape = 21, stroke = 0.5) +
        ggplot2::scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                                      midpoint = 0, limits = c(-1, 1), breaks = c(-1, 0, 1)) +
        ggplot2::scale_colour_manual(values = c("<=0.05" = "black", ">0.05" = "grey"),
                                     drop = FALSE) + theme +
        ggplot2::labs(x = "Cancer type", y = "Symbol", fill = "Spearman cor.",
                       colour = "Source FDR", title = paste(label, title))
    }
    # First source variant: CNV_Methylation_mRNA.
    data$point_size <- pmax(data$correlation * 5, 2)
    plot <- make_plot(data) + ggplot2::guides(size = "none") +
      ggplot2::labs(caption = "Original size rule: pmax(rho x 5, 2). Outline uses source FDR; unavailable entries are omitted.")
    paths[[association]] <- .panfam_original_assoc_save(function() print(plot), directory,
                                                        output_name, width, height)

    # Second source variant: CNV_Methylation_mRNA1. Retain its FDR threshold
    # and 1.5-to-5 point sizes, while making legend labels match the true FDR
    # values rather than applying a second unrelated size rescaling.
    data$neg_log10_fdr <- -log10(data$fdr + 1e-300)
    maximum <- max(data$neg_log10_fdr)
    size_for <- function(values) {
      if (maximum <= 25) return(ifelse(values < 25, 1.5, 5))
      ifelse(values < 25, 1.5, 1.5 + 3.5 * pmin(1, (values - 25) / (maximum - 25)))
    }
    data$point_size <- ifelse(data$fdr <= 0.05, size_for(data$neg_log10_fdr), 1.5)
    cancer_order <- names(sort(tapply(data$neg_log10_fdr, data$Cancer, sum), decreasing = TRUE))
    gene_order <- names(sort(tapply(data$neg_log10_fdr, data$Gene, sum), decreasing = TRUE))
    data$Cancer <- factor(data$Cancer, levels = cancer_order)
    data$Gene <- factor(data$Gene, levels = gene_order)
    plot <- make_plot(data)
    breaks <- unique(c(25, 50, 75, 100, round(maximum, 1)))
    breaks <- sort(breaks[breaks >= 25 & breaks <= maximum])
    if (maximum > 25 && length(breaks)) {
      plot <- plot + ggplot2::scale_size_identity(name = "-log10(FDR)", guide = "legend",
                                                 breaks = size_for(breaks), labels = breaks)
    } else {
      plot <- plot + ggplot2::scale_size_identity(guide = "none")
    }
    plot <- plot + ggplot2::labs(caption = paste0(
      "Original FDR-size variant: minimum size below -log10(FDR) = 25; larger values scale to size 5.\n",
      "Rows and columns are ordered by summed -log10(source FDR). Unavailable entries are omitted."))
    paths[[paste0(association, "_fdr")]] <- .panfam_original_assoc_save(function() print(plot),
                                                                      directory, paste0(output_name, "_fdr"), width, height)
  }
  invisible(paths)
}

.panfam_original_response_bars <- function(data, title, facet = FALSE, caption = NULL) {
  # Keep the original Yes/No/SD categories when the raw source annotation is
  # available. The binary Response column remains unchanged for statistics.
  raw_response <- if ("Response_original" %in% names(data)) data$Response_original else data$Response
  raw_response <- toupper(trimws(as.character(raw_response)))
  data$ResponsePlot <- ifelse(raw_response %in% c("CR", "PR", "R", "CRPR", "CR/PR", "YES", "RESPONDER"), "Yes",
                             ifelse(raw_response %in% c("PD", "NR", "N", "NO", "NON-RESPONDER"), "No",
                                    ifelse(raw_response == "SD", "SD", NA_character_)))
  valid <- is.finite(data$Score) & !is.na(data$ScoreGroup) & !is.na(data$ResponsePlot)
  data <- data[valid, , drop = FALSE]
  if (!nrow(data)) return(NULL)
  data$ScoreGroup <- factor(data$ScoreGroup, levels = c("Low", "High"))
  count <- as.data.frame(table(Cohort = data$Cohort, ScoreGroup = data$ScoreGroup,
                               Response = factor(data$ResponsePlot, levels = c("Yes", "No", "SD"))),
                         stringsAsFactors = FALSE)
  count <- count[count$Freq > 0L, , drop = FALSE]
  denominators <- stats::ave(count$Freq, count$Cohort, count$ScoreGroup, FUN = sum)
  count$percent <- 100 * count$Freq / denominators
  palette <- c(Yes = "tomato", No = "deepskyblue", SD = "gray")
  if (!facet && !is.null(caption)) caption <- paste(strwrap(caption, width = 55), collapse = "\n")
  plot <- ggplot2::ggplot(count, ggplot2::aes(x = ScoreGroup, y = percent, fill = Response)) +
    ggplot2::geom_bar(stat = "identity", colour = NA) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(round(percent), "%")),
                        position = ggplot2::position_stack(vjust = 0.5), size = 5) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_y_continuous(limits = c(0, 100), expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(title = title, y = "Percent weight", x = NULL, fill = "Response", caption = caption) +
    ggplot2::theme_minimal() + ggplot2::theme(
      text = ggplot2::element_text(size = 14),
      plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5),
      panel.border = ggplot2::element_blank(),
      axis.ticks.length = grid::unit(5, "pt"),
      axis.ticks.y.left = ggplot2::element_line(colour = "black"),
      axis.ticks.x = ggplot2::element_blank(), axis.line = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(size = 8, hjust = 0)
    )
  if (facet) {
    plot <- plot + ggplot2::facet_wrap(~Cohort, nrow = 1)
  } else {
    # Match the original hand-drawn asymmetric panel borders from plot_km.
    plot <- plot + ggplot2::annotation_custom(grid::linesGrob(
      x = grid::unit(c(0, 0), "npc"), y = grid::unit(c(0, 1), "npc"), gp = grid::gpar(lwd = 2))) +
      ggplot2::annotation_custom(grid::linesGrob(
        x = grid::unit(c(0, 1), "npc"), y = grid::unit(c(1, 1), "npc"), gp = grid::gpar(lwd = 2))) +
      ggplot2::annotation_custom(grid::linesGrob(
        x = grid::unit(c(1, 1), "npc"), y = grid::unit(c(0, 1), "npc"), gp = grid::gpar(lwd = 1))) +
      ggplot2::annotation_custom(grid::linesGrob(
        x = grid::unit(c(0, 1), "npc"), y = grid::unit(c(0, 0), "npc"), gp = grid::gpar(lwd = 1)))
  }
  plot
}

.panfam_original_immunotherapy_plot <- function(sample_table, directory, label) {
  # Source: plot_km / plot_km1 / plot_km2, whose response figure is a stacked
  # response proportion by median-split score group, next to a KM curve.
  # The earlier modular response-by-score boxplot was not in that source.
  if (!nrow(sample_table) || !any(is.finite(sample_table$Score))) return(invisible(NULL))
  cohorts <- unique(as.character(sample_table$Cohort[!is.na(sample_table$Cohort)]))
  grouped <- list()
  paths <- list()
  for (cohort in cohorts) {
    data <- sample_table[!is.na(sample_table$Cohort) & sample_table$Cohort == cohort &
                           is.finite(sample_table$Score), , drop = FALSE]
    if (!nrow(data)) next
    if (!"ScoreGroup" %in% names(data)) {
      cutoff <- stats::median(data$Score)
      data$ScoreGroup <- factor(ifelse(data$Score > cutoff, "High", "Low"), levels = c("Low", "High"))
    } else {
      data$ScoreGroup <- factor(data$ScoreGroup, levels = c("Low", "High"))
    }
    grouped[[cohort]] <- data
    bar <- .panfam_original_response_bars(data, paste(label, cohort), caption =
      "Within-cohort median split: Low <= median, High > median. Only known groups/responses are counted.")
    if (!is.null(bar)) {
      safe <- gsub("[^[:alnum:]_-]", "_", cohort)
      paths[[paste0(cohort, "_response")]] <- .panfam_original_assoc_save(function() print(bar),
        file.path(directory, "immunotherapy_response_cohorts"), safe, width = 4.5, height = 4)
    }
    if (!all(c("Time", "Event") %in% names(data))) next
    units <- if ("Time_unit" %in% names(data)) unique(data$Time_unit[!is.na(data$Time_unit)]) else character()
    known_days <- length(units) == 1L && units == "days"
    if (known_days) data$Time <- data$Time / 365
    fitted <- .panfam_survival_fit(data)
    if (is.null(fitted$fit)) next
    endpoint <- if ("Endpoint" %in% names(data)) as.character(data$Endpoint[[1L]]) else "Survival"
    x_label <- if (known_days) "Time (Years)" else paste(endpoint, "time (source units)")
    theme <- survminer::theme_survminer() + ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 16),
      legend.title = ggplot2::element_text(face = "bold", size = 12),
      legend.text = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(face = "bold", size = 13),
      axis.text = ggplot2::element_text(size = 11),
      panel.grid.major = ggplot2::element_line(colour = "grey90"),
      panel.background = ggplot2::element_rect(fill = "white")
    )
    km <- survminer::ggsurvplot(
      fitted$fit, data = fitted$data, conf.int = FALSE,
      pval = paste0("Log-rank P = ", format.pval(fitted$km$p_value, digits = 3)),
      pval.method = FALSE, pval.size = 5,
      surv.median.line = "hv", legend.title = "Strata", legend.labs = c("Low", "High"),
      legend = c(0.8, 0.9), xlab = x_label, ylab = "Survival Probability",
      palette = c("#D95F02", "#1B9E77"), ggtheme = theme,
      font.tickslab = 11, risk.table = FALSE, title = paste(label, cohort, endpoint)
    )$plot
    complete_bar <- .panfam_original_response_bars(fitted$data, cohort)
    if (is.null(complete_bar)) next
    combined <- patchwork::wrap_plots(km, complete_bar, widths = c(2, 1)) +
      patchwork::plot_annotation(caption = paste0(
        "Complete survival cases; cohort median split. KM P is unadjusted log-rank P.\n",
        "Bars show known-response proportions in the same survival groups; missing responses excluded from denominators."),
        theme = ggplot2::theme(plot.caption = ggplot2::element_text(size = 8, hjust = 0)))
    safe <- gsub("[^[:alnum:]_-]", "_", cohort)
    paths[[paste0(cohort, "_KM_response")]] <- .panfam_original_assoc_save(function() print(combined),
      file.path(directory, "immunotherapy_KM_response"), safe, width = 9, height = 4.5)
  }
  if (length(grouped)) {
    all_data <- do.call(rbind, unname(grouped))
    plot <- .panfam_original_response_bars(all_data, paste(label, "immunotherapy response"), facet = TRUE,
      caption = paste0("Original plot_km response design. Within-cohort median split; known-response denominators.\n",
                       "Scores are computed within each cohort; absolute scores are not compared across cohorts."))
    if (!is.null(plot)) {
      paths$overview <- .panfam_original_assoc_save(function() print(plot), directory,
        "immunotherapy_response", width = max(8, length(grouped) * 3), height = 4.5)
    }
  }
  invisible(paths)
}
