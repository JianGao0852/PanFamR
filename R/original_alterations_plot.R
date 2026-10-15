# Restore the original Alterations() figure: event-count tiles accompanied by
# alteration-type proportions. Counts are aggregated from the modular tables;
# the Total row is supplied separately because one sample may have several
# altered genes/types and cannot be recovered by summing event counts.
.panfam_original_alterations_plot <- function(frequencies, sample_totals, directory, label) {
  if (!nrow(frequencies)) return(invisible(NULL))
  required <- c("Gene", "Alteration", "Cancer", "n_altered")
  if (!all(required %in% names(frequencies))) {
    stop("Alteration plot requires Gene, Alteration, Cancer and n_altered.", call. = FALSE)
  }
  values <- frequencies[is.finite(frequencies$n_altered) & frequencies$n_altered > 0 &
      !is.na(frequencies$Gene) & !is.na(frequencies$Cancer) & !is.na(frequencies$Alteration),
    required, drop = FALSE]
  if (!nrow(values)) return(invisible(NULL))
  if (!all(c("Cancer", "n_altered_samples") %in% names(sample_totals))) {
    stop("Total row requires Cancer and n_altered_samples from unique altered samples.", call. = FALSE)
  }
  if (anyDuplicated(as.character(sample_totals$Cancer))) {
    stop("Alteration sample totals must have one row per cancer.", call. = FALSE)
  }
  genes <- sort(unique(as.character(values$Gene)))
  cancers <- sort(unique(as.character(values$Cancer)))
  types <- sort(unique(trimws(as.character(values$Alteration))))
  values$Gene <- factor(as.character(values$Gene), levels = genes)
  values$Cancer <- factor(as.character(values$Cancer), levels = cancers)
  values$Alteration <- factor(trimws(as.character(values$Alteration)), levels = types)

  event_counts <- as.matrix(stats::xtabs(n_altered ~ Gene + Cancer, data = values))
  if ("n_profiled" %in% names(frequencies)) {
    coverage <- frequencies[is.finite(frequencies$n_profiled) & frequencies$n_profiled > 0L, , drop = FALSE]
    covered_pairs <- unique(paste(coverage$Gene, coverage$Cancer, sep = "\r"))
    for (cancer in cancers) {
      event_counts[!paste(genes, cancer, sep = "\r") %in% covered_pairs, cancer] <- NA_real_
    }
  }
  total <- sample_totals$n_altered_samples[match(cancers, as.character(sample_totals$Cancer))]
  total[!is.finite(total) | total < 0] <- NA_real_
  total_name <- if ("Total" %in% genes) "Total altered samples" else "Total"
  count_matrix <- rbind(event_counts, total)
  rownames(count_matrix)[nrow(count_matrix)] <- total_name
  colnames(count_matrix) <- cancers

  row_type_counts <- as.matrix(stats::xtabs(n_altered ~ Gene + Alteration, data = values))
  col_type_counts <- as.matrix(stats::xtabs(n_altered ~ Cancer + Alteration, data = values))
  row_type_proportions <- row_type_counts / rowSums(row_type_counts)
  col_type_proportions <- col_type_counts / rowSums(col_type_counts)
  row_type_proportions[!is.finite(row_type_proportions)] <- 0
  col_type_proportions[!is.finite(col_type_proportions)] <- 0
  row_type_proportions <- rbind(row_type_proportions, rep(0, length(types)))
  rownames(row_type_proportions)[nrow(row_type_proportions)] <- total_name
  colnames(row_type_proportions) <- types

  count_colors <- function(x) {
    colors <- rep("white", length(x))
    colors[is.finite(x) & x >= 10 & x < 20] <- "#ded3cf"
    colors[is.finite(x) & x >= 20 & x < 40] <- "#669dbb"
    colors[is.finite(x) & x >= 40 & x < 60] <- "#fdb1a4"
    colors[is.finite(x) & x >= 60] <- "#99253a"
    colors
  }
  stack_colors <- c(MUT = "#fc9a99", AMP = "#b2df8a", FUSION = "#34a12e", HOMDEL = "#1f78b4")
  missing_types <- setdiff(types, names(stack_colors))
  if (length(missing_types)) {
    stack_colors[missing_types] <- grDevices::rainbow(length(missing_types))
  }

  # Keep the original cell density, top/right stacked annotations, and legend.
  # Expand the device for large custom families or cancer panels rather than
  # shrinking labels or clipping them at the fixed original 10 x 6 inch size.
  label_width <- max(nchar(rownames(count_matrix), type = "width"), na.rm = TRUE)
  cancer_width <- max(nchar(cancers, type = "width"), na.rm = TRUE)
  width <- max(10, length(cancers) * 0.24 + label_width * 0.085 + 2.2)
  height <- max(6, nrow(count_matrix) * 0.24 + cancer_width * 0.085 + 1.9)
  draw <- function() {
    row_annotation <- ComplexHeatmap::rowAnnotation(
      " " = ComplexHeatmap::anno_barplot(row_type_proportions,
        gp = grid::gpar(fill = stack_colors[colnames(row_type_proportions)]),
        border = FALSE, bar_width = 0.8, axis = FALSE),
      show_annotation_name = FALSE)
    column_annotation <- ComplexHeatmap::columnAnnotation(
      " " = ComplexHeatmap::anno_barplot(col_type_proportions,
        gp = grid::gpar(fill = stack_colors[colnames(col_type_proportions)]),
        border = FALSE, bar_width = 0.8, axis = FALSE),
      show_annotation_name = FALSE)
    alteration_legend <- ComplexHeatmap::Legend(labels = names(stack_colors), title = "Alter",
      legend_gp = grid::gpar(fill = stack_colors), nrow = 1,
      title_gp = grid::gpar(fontsize = 11), labels_gp = grid::gpar(fontsize = 10),
      title_position = "leftcenter")
    heatmap <- ComplexHeatmap::Heatmap(count_matrix,
      col = count_colors, na_col = "grey90", cluster_rows = FALSE, cluster_columns = FALSE,
      top_annotation = column_annotation, right_annotation = row_annotation,
      cell_fun = function(j, i, x, y, width, height, fill) {
        text <- if (is.finite(count_matrix[i, j])) format(count_matrix[i, j], trim = TRUE, scientific = FALSE) else "NA"
        grid::grid.text(text, x, y, gp = grid::gpar(fontsize = 8))
        grid::grid.rect(x, y, width, height, gp = grid::gpar(col = "grey80", fill = NA, lwd = 0.5))
      },
      show_column_names = TRUE, show_row_names = TRUE, row_names_side = "left",
      row_names_gp = grid::gpar(fontsize = 11), column_names_gp = grid::gpar(fontsize = 11),
      show_heatmap_legend = FALSE,
      row_names_max_width = ComplexHeatmap::max_text_width(rownames(count_matrix), gp = grid::gpar(fontsize = 11)),
      column_names_max_height = ComplexHeatmap::max_text_width(colnames(count_matrix), gp = grid::gpar(fontsize = 11)))
    ComplexHeatmap::draw(heatmap, annotation_legend_list = list(alteration_legend),
      annotation_legend_side = "top", column_title = paste(label, "alteration counts"))
  }
  paths <- .panfam_original_assoc_save(draw, directory, "alteration_summary", width, height)
  invisible(list(paths = paths, event_counts = count_matrix,
    row_type_proportions = row_type_proportions, column_type_proportions = col_type_proportions))
}
