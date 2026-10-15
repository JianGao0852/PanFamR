.panfam_immune_analysis <- function(scores, expression, directory, method, panels = NULL) {
  datasets <- c(Immune_cells = "ssgsea_immucells", Immune_scores = "immuscore",
                TME_pathways = "TME_Pathaway_GSVA_Score", Hallmark_pathways = "HALLMARK_GSVA_Score")
  markers <- .panfam_marker_lists()
  available_panels <- c(names(datasets), names(markers))
  if (is.null(panels)) panels <- available_panels
  panels <- unique(match.arg(panels, available_panels, several.ok = TRUE))
  tables <- list()
  for (panel in intersect(panels, names(datasets))) {
    table <- .panfam_correlations(scores, .panfam_aux_object(datasets[[panel]]), panel, method)
    tables[[paste0(tolower(panel), "_correlations")]] <- table
    .panfam_correlation_plot(table, directory, paste0(tolower(panel), "_correlations"), panel)
  }
  for (panel in intersect(panels, names(markers))) {
    requested_genes <- markers[[panel]]
    genes <- intersect(requested_genes, names(expression))
    safe <- gsub(" ", "_", tolower(panel))
    tables[[paste0(safe, "_gene_coverage")]] <- data.frame(
      Panel = panel, Gene = requested_genes, Available = requested_genes %in% genes,
      reason = ifelse(requested_genes %in% genes, "", "Marker gene absent from expression input"),
      stringsAsFactors = FALSE
    )
    if (!length(genes)) next
    features <- expression[, genes, drop = FALSE]
    table <- .panfam_correlations(scores, features, panel, method)
    tables[[paste0(safe, "_correlations")]] <- table
    .panfam_correlation_plot(table, directory, paste0(safe, "_correlations"), panel)
  }
  tables
}

.panfam_tcga_gtex_analysis <- function(genes, cancers, directory, label) {
  data <- .panfam_aux_object("tcga_gtex")
  groups <- ifelse(data$sample_type == "TCGA_tumor", "Tumor",
                   ifelse(data$sample_type %in% c("TCGA_normal", "GTEx_normal"), "Normal", NA_character_))
  rows <- .panfam_bind(lapply(intersect(genes, names(data)), function(gene) {
    .panfam_bind(lapply(cancers, function(cancer) {
      keep <- data$project == cancer & !is.na(groups)
      value <- as.numeric(data[[gene]][keep])
      group <- groups[keep]
      cbind(Gene = gene, Cancer = cancer, .panfam_difference(value[group == "Tumor"], value[group == "Normal"]))
    }))
  }))
  table <- .panfam_fdr(rows)
  if (nrow(table)) {
    plot <- ggplot2::ggplot(table, ggplot2::aes(x = Cancer, y = Gene, fill = rank_biserial)) +
      ggplot2::geom_tile() + ggplot2::scale_fill_gradient2(low = "#397FB6", mid = "white", high = "#D65A45", limits = c(-1, 1), na.value = "grey90") +
      ggplot2::theme_minimal() + ggplot2::labs(title = paste(label, "TCGA/GTEx expression"), x = NULL, y = NULL, fill = "Rank effect") +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
    .panfam_save_plot(plot, directory, "tcga_gtex_differences", height = max(4, min(14, length(genes) * 0.22 + 2)))
  }
  list(tcga_gtex_gene_differences = table,
       tcga_gtex_gene_coverage = data.frame(Gene = genes, Available = genes %in% names(data)))
}

.panfam_event <- function(x) {
  value <- toupper(trimws(as.character(x)))
  ifelse(value %in% c("DEAD", "DECEASED", "1", "TRUE"), 1,
         ifelse(value %in% c("ALIVE", "LIVING", "0", "FALSE"), 0, NA_real_))
}

.panfam_response <- function(x) {
  value <- toupper(trimws(as.character(x)))
  ifelse(value %in% c("CR", "PR", "R", "CRPR", "CR/PR"), "Responder",
         ifelse(value %in% c("PD", "SD", "NR", "N"), "Non-responder", NA_character_))
}

.panfam_gse135222_input <- function() {
  raw <- utils::read.delim(panfam_data_path("auxiliary/GSE135222_GEO_RNA-seq_omicslab_exp.tsv"), check.names = FALSE)
  annotation <- utils::read.delim(panfam_data_path("auxiliary/gencode.gene.info.v22.tsv"), check.names = FALSE)
  metadata <- utils::read.delim(panfam_data_path("auxiliary/GSE135222_sample_metadata.tsv"), check.names = FALSE)
  keys <- sub("[.].*", "", annotation$gene_id)
  symbols <- annotation$gene_name[match(sub("[.].*", "", raw$gene_id), keys)]
  ids <- metadata$sample_id[match(names(raw)[-1], metadata$raw_sample)]
  if (anyNA(ids) || anyDuplicated(ids)) stop("GSE135222 expression sample mapping is incomplete.", call. = FALSE)
  input <- data.frame(GENE_SYMBOL = symbols, raw[, -1, drop = FALSE], check.names = FALSE)
  names(input)[-1] <- ids
  list(expression = input, metadata = metadata)
}

.panfam_immunotherapy_analysis <- function(genes, min_size, directory, label) {
  cohort_files <- list(
    GSE78220 = c("Melanoma-GSE78220.Response.Rds", "Melanoma-GSE78220.Response_clin.Rds"),
    GSE91061 = c("Melanoma-GSE91061.Response.Rds", "Melanoma-GSE91061.Response_clin.Rds"),
    GSE135222 = c("NSCLC_GSE135222.Response.Rds", "NSCLC_GSE135222.Response_clin.Rds"),
    Braun2020 = c("RCC-Braun_2020.Response.Rds", "RCC-Braun_2020.Response_clin.Rds")
  )
  samples <- response <- survival <- km <- coverage <- filtering <- list()
  for (cohort in names(cohort_files)) {
    files <- cohort_files[[cohort]]
    input <- as.data.frame(readRDS(panfam_data_path(paste0("auxiliary/", files[[1]]))))
    clinical <- as.data.frame(readRDS(panfam_data_path(paste0("auxiliary/", files[[2]]))))
    endpoint <- if (cohort == "GSE135222") "PFS" else "OS"
    if (cohort == "GSE135222") {
      native <- .panfam_gse135222_input()
      input <- native$expression
      idx <- match(as.character(clinical$sample_id), native$metadata$sample_id)
      clinical$PFS_time <- native$metadata$PFS_time[idx]
      clinical$PFS_event <- native$metadata$PFS_event[idx]
    }
    n_original <- nrow(clinical)
    treatment <- toupper(trimws(as.character(clinical$Treatment)))
    # Baseline samples are used; when stage is unannotated, preserve the sample
    # and explicitly count it in the filtering table.
    baseline <- is.na(treatment) | treatment == "" | treatment == "PRE"
    therapy <- toupper(as.character(clinical$Therapy))
    immune_therapy <- is.na(therapy) | therapy == "" | grepl("PD-?1|PD-?L1|CTLA", therapy)
    clinical <- clinical[baseline & immune_therapy, , drop = FALSE]
    clinical <- clinical[order(clinical$sample_id), , drop = FALSE]
    patient <- as.character(clinical$patient_name)
    missing_patient <- is.na(patient) | !nzchar(patient)
    patient[missing_patient] <- as.character(clinical$sample_id[missing_patient])
    clinical <- clinical[!duplicated(patient), , drop = FALSE]
    filtering[[cohort]] <- data.frame(Cohort = cohort, original_samples = n_original,
                                      baseline_immunotherapy_unique_patients = nrow(clinical),
                                      unknown_treatment_stage = sum(is.na(clinical$Treatment) | clinical$Treatment == ""),
                                      survival_endpoint = endpoint,
                                      expression_source = if (cohort == "GSE135222") "Native expression TSV + GEO sample mapping" else "Provided expression RDS",
                                      stringsAsFactors = FALSE)
    symbols <- trimws(as.character(input$GENE_SYMBOL))
    valid <- !is.na(symbols) & nzchar(symbols)
    matrix <- as.matrix(input[valid, setdiff(names(input), "GENE_SYMBOL"), drop = FALSE])
    suppressWarnings(storage.mode(matrix) <- "double")
    symbols <- symbols[valid]
    sums <- rowsum(matrix, group = symbols, reorder = FALSE)
    counts <- tabulate(match(symbols, rownames(sums)), nbins = nrow(sums))
    matrix <- sums / counts
    matched_columns <- match(as.character(clinical$sample_id), colnames(matrix))
    clinical <- clinical[!is.na(matched_columns), , drop = FALSE]
    matrix <- matrix[, matched_columns[!is.na(matched_columns)], drop = FALSE]
    complete_genes <- rowSums(!is.finite(matrix)) == 0L
    filtering[[cohort]]$excluded_nonfinite_genes <- sum(!complete_genes)
    matrix <- matrix[complete_genes, , drop = FALSE]
    available <- genes %in% rownames(matrix)
    coverage[[cohort]] <- data.frame(Cohort = cohort, Gene = genes, Available = available)
    if (sum(available) < min_size) {
      filtering[[cohort]]$reason <- "Too few family genes in cohort expression data"
      next
    }
    if (ncol(matrix) < 2L) {
      filtering[[cohort]]$reason <- "Fewer than two matched baseline patients"
      next
    }
    # Column-wise TPM rescaling/log transformation is monotone, so raw ranks
    # yield the same ssGSEA ordering without creating rounding-induced ties.
    score <- .panfam_score_matrix(matrix, stats::setNames(list(genes), label), min_size)
    filtering[[cohort]]$background_notes <- paste(attr(score, "background_notes"), collapse = "; ")
    time <- if (endpoint == "PFS") clinical$PFS_time else suppressWarnings(as.numeric(clinical[["overall survival (days)"]]))
    event <- if (endpoint == "PFS") clinical$PFS_event else .panfam_event(clinical[["vital status"]])
    unit <- if (endpoint == "PFS") "GEO source units (not specified)" else "days"
    data <- data.frame(Sample = as.character(clinical$sample_id), Cohort = cohort, Endpoint = endpoint,
                       Score = as.numeric(score[1, ]), Response = .panfam_response(clinical$response),
                       Time = time, Event = event, Time_unit = unit, stringsAsFactors = FALSE)
    samples[[cohort]] <- data
    test <- .panfam_difference(data$Score[data$Response %in% "Responder"], data$Score[data$Response %in% "Non-responder"])
    names(test)[names(test) == "n_tumor"] <- "n_responder"
    names(test)[names(test) == "n_normal"] <- "n_non_responder"
    response[[cohort]] <- cbind(Cohort = cohort, test)
    fit <- .panfam_survival_fit(data)
    survival[[cohort]] <- cbind(Cohort = cohort, Endpoint = endpoint, fit$cox)
    km[[cohort]] <- cbind(Cohort = cohort, Endpoint = endpoint, fit$km)
    .panfam_km_plot(fit, file.path(directory, "immunotherapy_KM", paste0(cohort, ".pdf")), paste(label, cohort, endpoint),
                    time_label = if (endpoint == "PFS") "PFS time (GEO source units)" else "Time (days)")
    filtering[[cohort]]$reason <- if (!any(is.finite(data$Time) & data$Event %in% c(0, 1))) "Survival outcome unavailable; response analysis retained" else ""
  }
  filtering <- lapply(filtering, function(x) {
    if (!"reason" %in% names(x)) x$reason <- ""
    if (!"background_notes" %in% names(x)) x$background_notes <- ""
    x
  })
  sample_table <- .panfam_bind(samples)
  if (nrow(sample_table) && any(!is.na(sample_table$Response))) {
    plot <- ggplot2::ggplot(sample_table[!is.na(sample_table$Response), , drop = FALSE], ggplot2::aes(x = Response, y = Score, fill = Response)) +
      ggplot2::geom_boxplot(outlier.shape = NA) + ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.15, seed = 1), size = 0.7) +
      ggplot2::facet_wrap(~Cohort, scales = "free_y", nrow = 1) + ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)) +
      ggplot2::labs(title = paste(label, "immunotherapy response"), x = NULL, y = "Within-cohort ssGSEA score")
    .panfam_save_plot(plot, directory, "immunotherapy_response", width = 13)
  }
  list(immunotherapy_samples = sample_table, immunotherapy_response = .panfam_fdr(.panfam_bind(response)),
       immunotherapy_survival = .panfam_fdr(.panfam_bind(survival)), immunotherapy_km = .panfam_fdr(.panfam_bind(km)),
       immunotherapy_gene_coverage = .panfam_bind(coverage), immunotherapy_filtering = .panfam_bind(filtering))
}

.panfam_alterations_analysis <- function(genes, cancers) {
  root <- system.file("extdata", "auxiliary", package = "PanFamR")
  files <- list.files(root, pattern = "^alterations_across_samples_.*[.]tsv$", full.names = TRUE)
  rows <- list()
  seen <- character()
  for (file in files) {
    table <- utils::read.delim(file, check.names = FALSE, quote = "", stringsAsFactors = FALSE)
    columns <- grep(": ", names(table), value = TRUE)
    columns <- columns[sub(": .*", "", columns) %in% genes & !columns %in% seen]
    if (!length(columns)) next
    seen <- c(seen, columns)
    cancer <- toupper(sub("_.*", "", table[["Study ID"]]))
    for (column in columns) {
      status <- as.character(table[[column]])
      keep <- !duplicated(paste(cancer, table[["Sample ID"]], sep = ":"))
      for (type in intersect(cancers, unique(cancer))) {
        values <- status[keep & cancer == type]
        profiled <- !is.na(values) & nzchar(values) & !values %in% c("not profiled", "NA")
        altered <- profiled & values != "no alteration"
        key <- paste(column, type)
        rows[[key]] <- data.frame(Gene = sub(": .*", "", column), Alteration = sub(".*: ", "", column),
                                  Cancer = type, n_profiled = sum(profiled), n_altered = sum(altered),
                                  fraction_altered = if (sum(profiled)) sum(altered) / sum(profiled) else NA_real_)
      }
    }
  }
  available <- unique(sub(": .*", "", seen))
  list(alteration_frequencies = .panfam_bind(rows), alteration_gene_coverage = data.frame(Gene = genes, Available = genes %in% available))
}

.panfam_cnv_methylation_analysis <- function(genes, cancers, directory, associations = c("CNV", "Methylation")) {
  associations <- unique(match.arg(associations, c("CNV", "Methylation"), several.ok = TRUE))
  root <- system.file("extdata", "family_tables", package = "PanFamR")
  files <- list.files(root, pattern = "^(CnvAndExpressionTable|ExpressionAndMethylationTable)[.]xlsx$", recursive = TRUE, full.names = TRUE)
  source_names <- c(CNV = "CnvAndExpressionTable.xlsx", Methylation = "ExpressionAndMethylationTable.xlsx")
  files <- files[basename(files) %in% source_names[associations]]
  rows <- list()
  for (file in files) {
    data <- as.data.frame(readxl::read_xlsx(file))
    data <- data[data$symbol %in% genes & toupper(data$cancertype) %in% cancers, , drop = FALSE]
    if (!nrow(data)) next
    rows[[file]] <- data.frame(Association = if (basename(file) == "CnvAndExpressionTable.xlsx") "CNV" else "Methylation",
                              Gene = data$symbol, Cancer = toupper(data$cancertype), correlation = data$spm, fdr = data$fdr,
                              source_table = paste(basename(dirname(file)), basename(file), sep = "/"))
  }
  table <- .panfam_bind(rows)
  if (nrow(table)) {
    table <- table[!duplicated(table[, c("Association", "Gene", "Cancer")]), , drop = FALSE]
    plot <- ggplot2::ggplot(table, ggplot2::aes(x = Cancer, y = Gene, fill = correlation)) +
      ggplot2::geom_tile() + ggplot2::facet_grid(Association ~ ., scales = "free_y", space = "free_y") +
      ggplot2::scale_fill_gradient2(low = "#397FB6", mid = "white", high = "#D65A45", limits = c(-1, 1), na.value = "grey90") +
      ggplot2::theme_minimal() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::labs(title = paste("Available source", paste(associations, collapse = "/"), "associations"), x = NULL, y = NULL)
    plot_name <- if (length(associations) == 2L) "cnv_methylation" else paste0(tolower(associations), "_associations")
    .panfam_save_plot(plot, directory, plot_name, height = max(5, min(16, length(genes) * 0.3 + 3)))
  }
  coverage <- .panfam_bind(lapply(associations, function(association) {
    available <- if (nrow(table)) table$Gene[table$Association == association] else character()
    data.frame(Association = association, Gene = genes, Available = genes %in% available)
  }))
  list(cnv_methylation_associations = table, cnv_methylation_gene_coverage = coverage)
}
