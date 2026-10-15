suppressPackageStartupMessages(library(PanFamR))
set.seed(20261014)
expression <- data.frame(Cancer = rep(c("A", "B"), each = 10),
                         Group = rep(rep(c("Normal", "Tumor"), each = 5), 2),
                         matrix(rexp(20 * 80), nrow = 20), check.names = FALSE)
names(expression)[-(1:2)] <- paste0("G", seq_len(80))
rownames(expression) <- paste0("Sample_", seq_len(20))
sets <- list(Family_A = paste0("G", 1:5), Family_B = paste0("G", 6:12))
batch <- panfam_score(sets, expression = expression, min_size = 5)
for (name in names(sets)) {
  single <- ssGsea_Score(expression, sets[[name]], family_name = name)
  stopifnot(isTRUE(all.equal(single$Score, batch[[name]]$Score, tolerance = 1e-12)))
}
cat("PASS: shared ranking reproduces independent single-set normalization\n")
small <- panfam_score(genes = c("G1", "G2"), family_name = "SmallFamily", expression = expression)
stopifnot(all(is.finite(small$SmallFamily$Score)))
cat("PASS: custom two-gene families are supported\n")
fails <- function(expr) inherits(tryCatch({force(expr); NULL}, error = function(e) e), "error")
stopifnot(fails(panfam_score(genes = "NOT_A_REAL_GENE", expression = expression)),
          fails(panfam_score(list("A" = "G1", "A" = "G2"), expression = expression)),
          fails(panfam_score(genes = "G1", expression = expression, min_size = 2)))
bad <- expression
bad[1, "G1"] <- NA_real_
stopifnot(fails(panfam_score(genes = "G1", expression = bad)))
cat("PASS: unavailable genes, duplicate labels, minimum size, and missing expression are rejected\n")
match_ids <- PanFamR:::.panfam_match(c("TCGA-AA-0001-01A", "TCGA-BB-0002-01A", "nonTCGA"),
                                    c("TCGA-AA-0001-01", "TCGA-BB-0002", "nonTCGA"))
stopifnot(identical(match_ids, 1:3))
ambiguous <- PanFamR:::.panfam_match("TCGA-AA-0001-01A", c("TCGA-AA-0001-02A", "TCGA-AA-0001-03A"))
stopifnot(is.na(ambiguous))
stopifnot(identical(PanFamR:::.panfam_event(c("Alive", "Dead", "#N/A", "")), c(0, 1, NA_real_, NA_real_)),
          identical(PanFamR:::.panfam_response(c("CR", "PD", "UNK", "NE")), c("Responder", "Non-responder", NA_character_, NA_character_)))
cat("PASS: unambiguous TCGA matching and missing response/event handling\n")
data <- data.frame(Score = 1:20, Time = sample(1:100, 20), Event = rep(c(0, 1), 10))
fit <- PanFamR:::.panfam_survival_fit(data)
stopifnot(fit$km$n_low == 10L, fit$km$n_high == 10L,
          all(as.character(fit$data$ScoreGroup) == ifelse(fit$data$Score > 10.5, "High", "Low")))
cat("PASS: median groups exactly follow score values without reassignment\n")
result <- panfam_analyze(sets, expression = expression,
                         output_dir = file.path(tempdir(), "PanFamR-regression-results"),
                         analyses = c("expression", "paired"), min_size = 5)
stopifnot(inherits(result, "panfam_result"), file.exists(file.path(result$output_dir, "report.html")),
          !any(result$run_log$Status == "error"),
          all(result$run_log$Status[result$run_log$Module == "paired"] == "skipped"))
cat("PASS: two-family automated reports and empty-pair status\n")
unicode <- panfam_analyze(gene_sets = stats::setNames(list(c("G1", "G2"), c("G3", "G4")), c("\u7b2c\u4e00\u7ec4", "\u7b2c\u4e8c\u7ec4")),
                          expression = expression, analyses = "expression",
                          output_dir = file.path(tempdir(), "PanFamR-unicode-results"))
stopifnot(identical(names(unicode$scores), c("\u7b2c\u4e00\u7ec4", "\u7b2c\u4e8c\u7ec4")),
          length(list.files(unicode$output_dir, pattern="^scores[.]csv$",recursive=TRUE)) == 2L,
          !anyDuplicated(tolower(unname(unicode$settings$family_directories))))
cat("PASS: Chinese labels retain distinct result directories\n")
