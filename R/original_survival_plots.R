# Restore the source script's survival graphics without changing the fitted
# models, median groups, observations, or statistical result tables.

.panfam_original_survival_p <- function(value) {
  out <- rep("NA", length(value))
  valid <- is.finite(value)
  out[valid & value < 0.001] <- "<0.001"
  out[valid & value >= 0.001] <- sprintf("%.3f", value[valid & value >= 0.001])
  out
}

.panfam_original_survival_devices <- function(path, width, height, draw) {
  stem <- tools::file_path_sans_ext(path)
  paths <- c(pdf = paste0(stem, ".pdf"), png = paste0(stem, ".png"))
  dir.create(dirname(paths[["pdf"]]), recursive = TRUE, showWarnings = FALSE)
  render <- function(open_device) {
    open_device()
    device <- grDevices::dev.cur()
    tryCatch(draw(), finally = {
      if (device %in% grDevices::dev.list()) grDevices::dev.off(device)
    })
  }
  render(function() {
    if (capabilities("cairo")) {
      grDevices::cairo_pdf(paths[["pdf"]], width = width, height = height)
    } else {
      grDevices::pdf(paths[["pdf"]], width = width, height = height)
    }
  })
  render(function() {
    grDevices::png(paths[["png"]], width = width, height = height,
                   units = "in", res = 300, bg = "white")
  })
  invisible(paths)
}

.panfam_original_km_plot <- function(result, path, title, time_label = "Time (days)") {
  if (is.null(result$fit)) return(invisible(NULL))
  data <- result$data
  if (!all(c("Time", "Event", "ScoreGroup") %in% names(data))) {
    stop("The fitted survival result lacks Time, Event, or ScoreGroup.", call. = FALSE)
  }
  # Labels and colors follow the fitted strata. Reversing the displayed legend
  # does not change the Low/High factor or any observation's assigned group.
  groups <- sub("^.*=", "", names(result$fit$strata))
  if (!length(groups) || !all(groups %in% c("Low", "High"))) {
    stop("The fitted survival strata must be Low and High.", call. = FALSE)
  }
  colors <- c(Low = "#0066CC", High = "#FFCC33")
  unknown_units <- grepl("source|not specified|unspecified|unknown", time_label,
                         ignore.case = TRUE)
  if ("Time_unit" %in% names(data)) {
    unknown_units <- unknown_units || any(grepl("source|not specified|unspecified|unknown",
                                               data$Time_unit, ignore.case = TRUE), na.rm = TRUE)
  }
  days <- !unknown_units && grepl("\\bdays?\\b", time_label, ignore.case = TRUE)
  divisor <- if (days) 365 else 1
  axis_label <- if (days) "Time(years)" else time_label
  if (unknown_units && grepl("\\bdays?\\b", axis_label, ignore.case = TRUE)) {
    axis_label <- "Time (source units; unspecified)"
  }
  # xscale changes axis labels only; the original fit and the source time
  # values remain untouched. A year is displayed as 365 known source days.
  breaks <- pretty(c(0, max(data$Time, na.rm = TRUE) / divisor), n = 5)
  break_by <- if (length(breaks) > 1L) diff(breaks)[[1L]] * divisor else NULL
  cox <- result$cox
  group_model <- all(c("hr_high_vs_low", "group_ci_lower", "group_ci_upper", "group_p_value") %in% names(cox))
  plot_title <- paste0(title, "\nMedian split: High > median; Low <= median")
  km <- survminer::ggsurvplot(
    result$fit, data = data, title = plot_title,
    legend = c(0.9, 0.83), legend.title = " ", legend.labs = groups,
    xlab = axis_label, ylab = "Survival probability",
    xscale = divisor, break.time.by = break_by,
    palette = unname(colors[groups]), conf.int = FALSE,
    pval = FALSE, risk.table = TRUE, risk.table.height = 0.25,
    ggtheme = ggplot2::theme_bw(base_size = 15) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.80, size = 13),
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(colour = "black", linewidth = 1.8, fill = NA)
      )
  )
  plot <- km$plot
  if (identical(groups, c("Low", "High"))) {
    plot <- plot + ggplot2::guides(colour = ggplot2::guide_legend(reverse = TRUE),
                                  fill = ggplot2::guide_legend(reverse = TRUE))
  }
  number <- function(value) {
    if (length(value) == 1L && is.finite(value)) sprintf("%.2f", value) else "NA"
  }
  estimate <- if (group_model) cox$hr_high_vs_low else cox$hr_per_sd
  lower <- if (group_model) cox$group_ci_lower else cox$ci_lower
  upper <- if (group_model) cox$group_ci_upper else cox$ci_upper
  p_value <- if (group_model) cox$group_p_value else cox$p_value
  lines <- c(
    paste0(if (group_model) "HR (High vs Low): " else "HR per score SD: ", number(estimate)),
    paste0("95% CI: ", number(lower), " - ", number(upper)),
    paste0("Cox P: ", .panfam_original_survival_p(p_value))
  )
  fdr_field <- if (group_model) "group_fdr" else "fdr"
  if (fdr_field %in% names(cox) && length(cox[[fdr_field]]) == 1L && is.finite(cox[[fdr_field]])) {
    lines <- c(lines, paste0("BH FDR: ", .panfam_original_survival_p(cox[[fdr_field]])))
  }
  lines <- c(lines, paste0("Log-rank P: ", .panfam_original_survival_p(result$km$p_value)))
  annotation <- do.call(grid::grobTree, lapply(seq_along(lines), function(i) {
    grid::textGrob(lines[[i]], x = 0.035, y = 0.08 + (length(lines) - i) * 0.05,
                   just = "left", gp = grid::gpar(fontsize = 10, col = "black"))
  }))
  plot <- plot + ggplot2::annotation_custom(annotation)
  .panfam_original_survival_devices(path, width = 6, height = 5, draw = function() print(plot))
  km$plot <- plot
  risk_path <- paste0(tools::file_path_sans_ext(path), "_with_risk_table.pdf")
  .panfam_original_survival_devices(risk_path, width = 6, height = 7,
                                    draw = function() print(km, newpage = TRUE))
}

.panfam_original_survival_forest <- function(cox, directory, label) {
  required <- c("Endpoint", "Cancer", "n", "events")
  if (!nrow(cox)) return(invisible(NULL))
  if (!all(required %in% names(cox))) stop("The Cox table lacks required forest-plot columns.", call. = FALSE)
  group_model <- all(c("hr_high_vs_low", "group_ci_lower", "group_ci_upper", "group_p_value") %in% names(cox))
  fields <- if (group_model) {
    c(estimate = "hr_high_vs_low", lower = "group_ci_lower", upper = "group_ci_upper", p_value = "group_p_value", fdr = "group_fdr")
  } else {
    c(estimate = "hr_per_sd", lower = "ci_lower", upper = "ci_upper", p_value = "p_value", fdr = "fdr")
  }
  if (!all(fields[names(fields) != "fdr"] %in% names(cox))) {
    stop("The Cox table lacks fitted hazard-ratio, confidence-interval, or P-value columns.", call. = FALSE)
  }
  if (!fields[["fdr"]] %in% names(cox)) cox[[fields[["fdr"]]]] <- NA_real_
  endpoints <- unique(as.character(cox$Endpoint))
  paths <- list()
  for (endpoint in endpoints) {
    data <- cox[cox$Endpoint == endpoint, , drop = FALSE]
    estimate <- data[[fields[["estimate"]]]]
    lower <- data[[fields[["lower"]]]]
    upper <- data[[fields[["upper"]]]]
    p_value <- data[[fields[["p_value"]]]]
    fdr <- data[[fields[["fdr"]]]]
    valid <- is.finite(estimate) & estimate > 0 &
      is.finite(lower) & lower > 0 & is.finite(upper) & upper >= lower
    if (!any(valid)) next
    interval <- rep("Not estimable", nrow(data))
    interval[valid] <- sprintf("%.3f (%.3f - %.3f)", estimate[valid], lower[valid], upper[valid])
    tabletext <- cbind(
      c("Cancer", as.character(data$Cancer)),
      c("Total(N)", as.character(data$n)),
      c("Events", as.character(data$events)),
      c("HR (95% CI)", interval),
      c("P (raw)", .panfam_original_survival_p(p_value)),
      c("BH FDR", .panfam_original_survival_p(fdr))
    )
    row_colors <- rep("black", nrow(data))
    significant <- valid & is.finite(p_value) & p_value < 0.05
    row_colors[significant & estimate > 1] <- "#E64B35"
    row_colors[significant & estimate < 1] <- "#4DBBD5"
    mean <- estimate
    mean[!valid] <- lower[!valid] <- upper[!valid] <- NA_real_
    horizontal <- list(`1` = grid::gpar(lwd = 2, col = "black"),
                       `2` = grid::gpar(lwd = 2, col = "black"))
    horizontal[[as.character(nrow(data) + 2L)]] <- grid::gpar(lwd = 2, col = "black")
    plot <- forestplot::forestplot(
      labeltext = tabletext,
      title = paste0(label, " | ", endpoint, "\n",
                     if (group_model) "Median-split Cox: High vs Low" else "Univariate Cox; HR per score SD",
                     "; colors: raw P < 0.05"),
      mean = c(NA_real_, mean), lower = c(NA_real_, lower), upper = c(NA_real_, upper),
      hrzl_lines = horizontal, is.summary = rep(FALSE, nrow(data) + 1L),
      zero = 1, boxsize = 0.4, lineheight = grid::unit(0.2, "cm"), xlog = FALSE,
      col = forestplot::fpColors(box = "#4DBBD5", lines = "black", zero = "gray50",
                                 text = c("black", row_colors)),
      lwd.ci = 2.5, ci.vertices = TRUE, ci.vertices.height = 0.02,
      clip = c(0, 6), xticks = seq(0, 6, by = 2),
      graph.pos = 5, graphwidth = grid::unit(5, "cm"),
      xlab = if (group_model) "Hazard ratio: High vs Low (95% CI)" else "Hazard ratio per score SD (95% CI)",
      txt_gp = forestplot::fpTxtGp(
        title = grid::gpar(fontsize = 12, fontface = "bold"),
        label = grid::gpar(cex = 0.8), ticks = grid::gpar(cex = 0.9)
      )
    )
    name <- if (length(endpoints) == 1L) "survival_forest" else paste0("survival_forest_", endpoint)
    paths[[endpoint]] <- .panfam_original_survival_devices(
      file.path(directory, paste0(name, ".pdf")), width = 9,
      height = max(5, nrow(data) * 0.23 + 1.8),
      draw = function() { grid::grid.newpage(); print(plot, newpage = FALSE) }
    )
  }
  invisible(paths)
}

.panfam_original_gene_survival_heatmap <- function(cox, directory, label) {
  if (!nrow(cox)) return(invisible(NULL))
  required <- c("Endpoint", "Gene", "Cancer", "p_value")
  if (!all(required %in% names(cox))) stop("The gene Cox table lacks required heatmap columns.", call. = FALSE)
  estimate_field <- if ("hr_per_unit" %in% names(cox)) "hr_per_unit" else "hr_per_sd"
  p_field <- if ("unit_p_value" %in% names(cox)) "unit_p_value" else "p_value"
  if (!estimate_field %in% names(cox)) stop("The gene Cox table lacks a fitted hazard ratio.", call. = FALSE)
  cancer_order <- c("ACC", "BLCA", "BRCA", "CESC", "CHOL", "COAD", "DLBC", "ESCA",
                    "GBM", "HNSC", "KICH", "KIRC", "KIRP", "LAML", "LGG", "LIHC",
                    "LUAD", "LUSC", "MESO", "OV", "PAAD", "PCPG", "PRAD", "READ",
                    "SARC", "SKCM", "STAD", "TGCT", "THCA", "THYM", "UCEC", "UCS", "UVM")
  endpoints <- unique(as.character(cox$Endpoint))
  paths <- list()
  for (endpoint in endpoints) {
    data <- cox[cox$Endpoint == endpoint, , drop = FALSE]
    genes <- unique(as.character(data$Gene))
    cancers <- unique(as.character(data$Cancer))
    cancers <- c(intersect(cancer_order, cancers), setdiff(cancers, cancer_order))
    status <- matrix(NA_real_, nrow = length(genes), ncol = length(cancers),
                     dimnames = list(genes, cancers))
    estimate <- data[[estimate_field]]
    p_value <- data[[p_field]]
    valid <- is.finite(estimate) & estimate > 0 & is.finite(p_value)
    value <- rep(NA_real_, nrow(data))
    value[valid] <- 0
    value[valid & p_value < 0.05 & estimate > 1] <- 1
    value[valid & p_value < 0.05 & estimate < 1] <- -1
    status[cbind(match(data$Gene, genes), match(data$Cancer, cancers))] <- value
    # Preserve zero-score genes and distinguish unavailable estimates from
    # non-significant fits; the source's significant-only ordering dropped them.
    risky_score <- rowSums(status == 1, na.rm = TRUE) - rowSums(status == -1, na.rm = TRUE)
    row_order <- order(-risky_score, seq_along(risky_score))
    status <- status[row_order, , drop = FALSE]
    risky_score <- risky_score[row_order]
    annotation <- ComplexHeatmap::rowAnnotation(
      RiskyScore = ComplexHeatmap::anno_barplot(
        risky_score, baseline = 0, border = TRUE,
        gp = grid::gpar(fill = ifelse(risky_score > 0, "#E64B35", "#4DBBD5"),
                        col = "grey80", lwd = 0.5),
        width = grid::unit(2, "cm"), bar_width = 1,
        ylim = range(c(-1, 1, risky_score))
      ), annotation_name_gp = grid::gpar(fontsize = 10)
    )
    plot <- ComplexHeatmap::Heatmap(
      status, col = circlize::colorRamp2(c(-1, 0, 1), c("#4DBBD5", "white", "#E64B35")),
      na_col = "grey85", cluster_rows = FALSE, cluster_columns = FALSE,
      show_row_names = TRUE, show_column_names = TRUE, right_annotation = annotation,
      column_title = paste0(label, " | ", endpoint, "\nUnivariate Cox (gene expression); raw P < 0.05"),
      column_title_gp = grid::gpar(fontsize = 11),
      heatmap_legend_param = list(title = "", at = c(-1, 0, 1),
                                  labels = c("Protective", "P >= 0.05", "Risky")),
      cell_fun = function(j, i, x, y, width, height, fill) {
        grid::grid.rect(x, y, width, height, gp = grid::gpar(col = "black", fill = NA, lwd = 0.5))
      }
    )
    name <- if (length(endpoints) == 1L) "gene_survival_cox_heatmap" else paste0("gene_survival_cox_heatmap_", endpoint)
    paths[[endpoint]] <- .panfam_original_survival_devices(
      file.path(directory, paste0(name, ".pdf")), width = 10,
      height = max(4, min(16, length(genes) * 0.23 + 2)),
      draw = function() ComplexHeatmap::draw(plot)
    )
  }
  invisible(paths)
}
