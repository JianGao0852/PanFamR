suppressPackageStartupMessages(library(PanFamR))

local({
  set.seed(20261014)
  expression <- data.frame(
    Cancer = rep(c("A", "B"), each = 12L),
    Group = rep(rep(c("Normal", "Tumor"), each = 6L), 2L),
    matrix(rexp(24L * 80L), nrow = 24L), check.names = FALSE
  )
  names(expression)[-(1:2)] <- paste0("G", seq_len(80L))
  rownames(expression) <- unlist(lapply(c("AA", "BB"), function(code) {
    paste0("TCGA-", code, "-", rep(sprintf("%04d", seq_len(6L)), 2L),
           rep(c("-11A", "-01A"), each = 6L))
  }))
  test_root <- tempfile("PanFamR-modular-")
  dir.create(test_root)
  assert_error <- function(expr, pattern) {
    err <- tryCatch({ force(expr); NULL }, error = function(e) e)
    stopifnot(inherits(err, "error"), grepl(pattern, conditionMessage(err), ignore.case = TRUE))
    invisible(err)
  }

  wd_before <- getwd()
  global_before <- ls(envir = .GlobalEnv, all.names = TRUE)
  lazy <- panfam_setup(work_dir = test_root, output_dir = "lazy",
                       genes = c("G1", "G2"), family_name = "Lazy")
  stopifnot(identical(getwd(), wd_before),
            identical(ls(envir = .GlobalEnv, all.names = TRUE), global_before),
            is.null(lazy$cache$expression), !lazy$cache$expression_ready,
            is.null(lazy$cache$scores), identical(lazy$cache$score_runs, 0L),
            nrow(panfam_status(lazy)) == 0L,
            identical(lazy$output_dir, normalizePath(file.path(test_root, "lazy"), winslash = "/")))
  assert_error(lazy$gene_sets <- list(Changed = "G3"), "locked binding")
  assert_error(lazy$cancers <- "A", "locked binding")
  stopifnot(environmentIsLocked(lazy), bindingIsLocked("gene_sets", lazy))
  lazy_report <- panfam_report(lazy)
  stopifnot(file.exists(lazy_report), lazy$cache$score_runs == 0L,
            nrow(panfam_status(lazy)) == 0L, is.null(lazy$cache$expression))
  cat("PASS: setup is lazy, preserves working/global state and locks configuration\n")

  project <- panfam_setup(work_dir = test_root, output_dir = "selected_steps",
                          genes = paste0("G", seq_len(5L)), family_name = "MyFamily",
                          expression = expression)
  stopifnot(project$cache$score_runs == 0L, nrow(panfam_status(project)) == 0L)
  difference <- panfam_expression(project)
  history <- panfam_status(project)
  stopifnot(inherits(difference, "panfam_step"),
            identical(history$Module, c("ssgsea", "expression")),
            !any(history$Status == "error"), project$cache$score_runs == 1L,
            nrow(difference$tables$MyFamily$score_differences) == 2L,
            nrow(difference$tables$MyFamily$gene_differences) == 10L,
            file.exists(file.path(difference$output_dirs[["MyFamily"]], "score_differences.csv")))
  scores_before <- project$cache$scores
  files_before <- sort(list.files(project$output_dir, recursive = TRUE, all.files = TRUE))
  scores_again <- panfam_ssgsea(project)
  stopifnot(identical(scores_again, scores_before), project$cache$score_runs == 1L,
            identical(panfam_status(project), history),
            identical(sort(list.files(project$output_dir, recursive = TRUE, all.files = TRUE)), files_before))
  cat("PASS: one expression step computes ssGSEA once; repeated scoring changes no outputs\n")

  paired <- panfam_paired(project)
  stopifnot(identical(panfam_status(project)$Module, c("ssgsea", "expression", "paired")),
            !any(paired$run_log$Status == "error"), project$cache$score_runs == 1L,
            nrow(paired$tables$MyFamily$paired_samples) == 12L,
            all(paired$tables$MyFamily$paired_score_differences$n_pairs == 6L))
  difference_again <- panfam_expression(project)
  stopifnot(project$cache$score_runs == 1L,
            difference_again$run_log$Run == 2L,
            difference_again$output_dirs[["MyFamily"]] != difference$output_dirs[["MyFamily"]],
            dir.exists(difference$output_dirs[["MyFamily"]]),
            file.exists(file.path(difference$output_dirs[["MyFamily"]], "result.rds")),
            identical(panfam_status(project)$Module, c("ssgsea", "expression", "paired", "expression")))
  history_before_report <- panfam_status(project)
  report <- panfam_report(project)
  report_text <- readLines(report, warn = FALSE, encoding = "UTF-8")
  stopifnot(file.exists(report), identical(panfam_status(project), history_before_report),
            project$cache$score_runs == 1L,
            any(grepl("expression_002", report_text, fixed = TRUE)))
  cat("PASS: independent steps preserve previous runs; reporting adds no analysis\n")

  unicode_names <- c("\u7b2c\u4e00\u7ec4", "\u7b2c\u4e8c\u7ec4")
  multiple <- panfam_setup(work_dir = test_root, output_dir = "unicode",
                           gene_sets = stats::setNames(list(c("G1", "G2"), c("G3", "G4")), unicode_names),
                           expression = expression)
  multiple_scores <- panfam_ssgsea(multiple)
  score_dirs <- panfam_status(multiple)$Output_dir
  family_csv <- utils::read.csv(file.path(multiple$output_dir, "gene_families.csv"), encoding = "UTF-8")
  log_csv <- utils::read.csv(file.path(multiple$output_dir, "run_log.csv"), encoding = "UTF-8")
  stopifnot(identical(names(multiple_scores), unicode_names),
            identical(unique(family_csv$Family), unicode_names),
            identical(log_csv$Family, unicode_names),
            !anyDuplicated(tolower(unname(multiple$family_directories))),
            length(unique(score_dirs)) == 2L,
            all(file.exists(file.path(score_dirs, "scores.csv"))), multiple$cache$score_runs == 1L)
  existing <- panfam_setup(work_dir = test_root, output_dir = "unicode",
                           genes = c("G1", "G2"), expression = expression)
  stopifnot(existing$output_dir != multiple$output_dir, dir.exists(multiple$output_dir))
  cat("PASS: multiple and Unicode families have distinct paths; existing projects are preserved\n")

  csv_path <- file.path(test_root, "csv-encoding.csv")
  csv_data <- data.frame(Label = c(unicode_names[[1L]], 'Quote " and comma, with\nnewline'),
                         Value = c(pi, NA_real_), Logical = c(TRUE, FALSE))
  PanFamR:::.panfam_write_csv(csv_data, csv_path)
  csv_again <- utils::read.csv(csv_path, encoding = "UTF-8", check.names = FALSE)
  stopifnot(identical(csv_again$Label, csv_data$Label),
            isTRUE(all.equal(csv_again$Value, csv_data$Value, tolerance = 1e-15)),
            identical(csv_again$Logical, csv_data$Logical))
  cat("PASS: CSV retains Unicode, quoted/newline fields, numeric precision and missing values\n")

  assert_error(panfam_expression(NULL), "panfam_setup")
  assert_error(panfam_status(list()), "panfam_setup")
  assert_error(panfam_setup(work_dir = test_root, genes = character()), "non-empty")
  assert_error(panfam_setup(work_dir = test_root, genes = c("G1", NA_character_)), "missing")
  assert_error(panfam_setup(work_dir = test_root, genes = "G1", min_size = 1.5), "positive integer")
  assert_error(panfam_setup(work_dir = test_root, genes = "G1", min_size = 0), "positive integer")
  assert_error(panfam_setup(work_dir = test_root, genes = "G1", cancers = NA_character_), "cancers")
  assert_error(panfam_survival(project, endpoints = "PFS"), "endpoints")
  assert_error(panfam_survival(project, endpoints = character()), "endpoints")
  assert_error(panfam_survival(project, endpoints = NA_character_), "endpoints")
  bad_cancer <- panfam_setup(work_dir = test_root, output_dir = "bad_cancer",
                             genes = "G1", expression = expression, cancers = "Unknown")
  assert_error(panfam_ssgsea(bad_cancer), "Cancer labels absent")
  stopifnot(bad_cancer$cache$score_runs == 0L, is.null(bad_cancer$cache$scores),
            nrow(panfam_status(bad_cancer)) == 0L)
  cat("PASS: invalid settings and unknown cancer labels produce informative errors\n")

  # Check public selector forwarding without accessing bundled clinical/cohort data.
  mock <- new.env(parent = asNamespace("PanFamR"))
  clone_function <- function(name) {
    fun <- get(name, envir = asNamespace("PanFamR"), inherits = FALSE)
    environment(fun) <- mock
    assign(name, fun, envir = mock)
  }
  for (name in c("panfam_survival", ".panfam_immune_step", "panfam_immune_cells",
                 "panfam_immune_scores", "panfam_tme", "panfam_hallmark", "panfam_checkpoints",
                 "panfam_chemokines", "panfam_chemokine_receptors", "panfam_immune_stimulators",
                 "panfam_immune_inhibitors", "panfam_cnv", "panfam_methylation",
                 "panfam_tcga_gtex", "panfam_alterations", "panfam_immunotherapy")) clone_function(name)
  mock$.panfam_run_step <- function(project, module, analyze, needs_scores = TRUE) {
    mock$module <- module
    mock$needs_scores <- needs_scores
    analyze(names(project$gene_sets)[[1L]], test_root)
  }
  mock$.panfam_project_expression <- function(project) {
    mock$expression_requests <- mock$expression_requests + 1L
    expression
  }
  mock$.panfam_project_cancers <- function(project) c("A", "B")
  mock$.panfam_survival_analysis <- function(scores, directory, label, endpoints) {
    mock$selected <- endpoints
    list()
  }
  mock$.panfam_immune_analysis <- function(scores, expression, directory, method, panels) {
    mock$selected <- panels
    list()
  }
  mock$.panfam_cnv_methylation_analysis <- function(genes, cancers, directory, associations) {
    mock$selected <- associations
    list()
  }
  mock$.panfam_tcga_gtex_analysis <- function(...) list()
  mock$.panfam_alterations_analysis <- function(...) list()
  mock$.panfam_immunotherapy_analysis <- function(...) list()
  mock$expression_requests <- 0L
  mock$panfam_survival(project, endpoints = "OS")
  stopifnot(identical(mock$module, "survival_OS"), identical(mock$selected, "OS"), mock$needs_scores)
  mock$panfam_survival(project, endpoints = c("PFI", "DSS", "PFI"))
  stopifnot(identical(mock$module, "survival_PFI_DSS"), identical(mock$selected, c("PFI", "DSS")))
  panels <- c(panfam_immune_cells = "Immune_cells", panfam_immune_scores = "Immune_scores",
              panfam_tme = "TME_pathways", panfam_hallmark = "Hallmark_pathways",
              panfam_checkpoints = "Immune checkpoints", panfam_chemokines = "Chemokines",
              panfam_chemokine_receptors = "Chemokine receptors",
              panfam_immune_stimulators = "Immune stimulators", panfam_immune_inhibitors = "Immune inhibitors")
  for (name in names(panels)) {
    get(name, envir = mock)(project)
    stopifnot(identical(mock$selected, unname(panels[[name]])), mock$needs_scores)
  }
  stopifnot(mock$expression_requests == 5L)
  mock$panfam_cnv(project)
  stopifnot(identical(mock$selected, "CNV"), identical(mock$module, "cnv"), !mock$needs_scores)
  mock$panfam_methylation(project)
  stopifnot(identical(mock$selected, "Methylation"), identical(mock$module, "methylation"), !mock$needs_scores)
  for (name in c("panfam_tcga_gtex", "panfam_alterations", "panfam_immunotherapy")) {
    get(name, envir = mock)(project)
    stopifnot(!mock$needs_scores)
  }
  stopifnot(identical(panfam_status(project), history_before_report), project$cache$score_runs == 1L)
  cat("PASS: survival, immune and alteration wrappers select only their requested analyses\n")
})
