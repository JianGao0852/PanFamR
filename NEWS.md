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
