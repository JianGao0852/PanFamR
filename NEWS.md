# PanFamR 0.3.1

* Restored the original plotting styles throughout the modular workflow: six-bin expression heatmaps with numeric labels, red/cyan TMB/MSI radar plots, blue/white/red ComplexHeatmap immune panels, alteration-count heatmaps with stacked type annotations, source CNV/methylation bubble plots, and immunotherapy response-composition bars with KM curves.
* Retained alteration frequency and coverage tables. Alteration heatmap cells count gene/cancer events, allowing multiple alteration types per sample; the Total row counts unique samples with at least one alteration in the user-defined input family.
* Restored Pearson as the default correlation method; Spearman remains available through project setup.
* Restored family survival plots to genuine fixed-median High-versus-Low Cox estimates, with group HR, confidence intervals, P values, and separately adjusted group FDR. Retained the additional continuous score-SD estimates in the result tables.
* Added individual-gene Cox estimates per original expression unit and restored the source Risky/Protective heatmap and count-difference annotation using raw P values.
* Restored yellow/blue ggsurvplot curves, table-style forest plots, and HR/CI/P annotations. Known day-based follow-up is displayed in years; unspecified source units remain unchanged.
* Save supported figures as PDF and 300 dpi PNG with a white background.
* Documented plotting quantities and raw-P significance conventions alongside BH-adjusted result tables. Source-script `logFC` means the tumor-minus-normal expression mean difference and is not automatically log2 fold change.
* Retained genuine median groups without manual sample reassignment or label reversal. Fixed-direction empirical ROC curves remain unsmoothed and can differ from inherited smoothed ROC output.
* Preserved custom gene-family definitions, independent analysis modules, bundled data, coverage reporting, and project output management.

# PanFamR 0.3.0

* Added a modular workflow: configure work and result paths and custom gene families once with `panfam_setup()`, then call each analysis separately.
* Added separate public interfaces for expression, paired comparisons, survival endpoints, TMB, MSI, immune panels, immunotherapy, alterations, CNV, and methylation.
* Added individual-gene survival summaries, tumor-score distributions, and fixed-direction exploratory ROC curves.
* Reuse ssGSEA scores within each project and preserve earlier results when repeating a step.
* Added status history and an HTML index containing only steps already executed.
* Write CSV output directly as UTF-8, retaining Unicode labels, quoted fields, numeric precision, and missing values.
* Updated the example script and English/Chinese documentation for step-by-step use in RStudio.
* Store large bundled inputs in 20 MiB pieces and automatically reconstruct expression and response data in the local cache.

# PanFamR 0.2.0

* Added a one-call workflow for arbitrary gene families with CSV, PDF/PNG, HTML reporting, and explicit data coverage.
* Added batch ssGSEA scoring compatible with GSVA 1.50+ and R 4.3+.
* Removed manual survival-group reassignment and high/low label swapping from inherited plots.
* Added fixed median survival groups, patient-aware matching, missing-outcome handling, and multiple-testing adjustment.
* Restored GSE135222 native expression and PFS annotation using bundled GEO sample mapping.
* Added automatic safe directory naming for Chinese family labels.

# PanFamR 0.1.0

* Added a portable package interface for ssGSEA scoring of named gene sets.
* Packaged the original downstream plotting functions.
