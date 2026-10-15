# Source-script plotting compatibility (PanFamR 0.3.1)

The modular functions use plotting adapters based on the supplied
`ssGSEA-Pan-Cancer.R`. Paths, family names, input selection, device management,
and unsupported-data handling are adapted for an installable package.

| Original function(s) | Modular entry point | Restored design |
| --- | --- | --- |
| `plot_diff_heatmap`, `plot_diff_heatmap1` | `panfam_expression`, `panfam_tcga_gtex` | Six expression-difference bins, original six-color palette, numeric labels, raw-P white masking, fixed cell aspect ratio |
| `ssGsea_Score_diff`, `ssGsea_Score_diff1` | `panfam_expression` | JCO grouped boxplots and Set1 boxplot/jitter variant |
| `ssGsea_Score_diff_Tumor` | `panfam_tumor_scores` | Median-ordered cancer-colored tumor boxplots |
| `ssGsea_Score_paired_boxplot`, `ssGsea_Score_one_paired_boxplot` | `panfam_paired` | NPG cyan/red paired boxes and patient lines; pooled and individual-cancer figures |
| `ssGsea_Score_pan_cancer_TMB`, `ssGsea_Score_pan_cancer_MSI` | `panfam_tmb`, `panfam_msi` | Red/cyan radar plots and raw-P labels |
| Immune, ICG, chemokine, receptor, stimulator/inhibitor, TME and Hallmark functions | Corresponding `panfam_*` modules | Unclustered ComplexHeatmap panels with observed-range blue/white/red colors and raw-P stars |
| `ssGsea_Score_pan_cancer_fores` | `panfam_survival` | Table-style forests using median-group High-versus-Low HR, confidence interval and raw P; separate group FDR retained |
| `plot_OS_KM_expression`, `ssGsea_Score_KM` | `panfam_survival` | Yellow/blue ggsurvplot curves, black border, group HR/CI/P labels; known days displayed in years |
| `plot_cox_heatmap1`, `plot_cox_heatmap2` | `panfam_gene_survival` | Risky/Protective discrete colors, cell borders, and signed count annotation |
| `Alterations` | `panfam_alterations` | Event-count heatmap with numeric cells and stacked alteration-type annotations |
| `CNV_Methylation_mRNA`, `CNV_Methylation_mRNA1` | `panfam_cnv`, `panfam_methylation` | Source bubble designs, source FDR outlines, plus an FDR-size variant |
| `plot_km`, `plot_km1`, `plot_km2` | `panfam_immunotherapy` | Within-cohort median-group response composition and KM/stacked-bar combinations |
| `ssGsea_Score_ROC` | `panfam_roc` | Original blue line, diagonal, AUC/CI annotation and individual-cancer layout |

## Statistical definitions

* Pearson is the default correlation method, matching the supplied script.
  Spearman remains an explicit alternative in project setup.
* The source expression heatmaps name `mean(tumor) - mean(normal)` as `logFC`.
  The package computes that same difference and labels it `Mean difference`;
  its units follow the expression input and are not assumed to be logarithmic.
* Source-style coloring/stars use raw P values; figure captions and module
  documentation distinguish these from adjusted values.
  Newly computed result tables additionally retain BH-adjusted values. Family
  grouped Cox uses `group_p_value`/`group_fdr`; continuous score-SD Cox uses
  `p_value`/`fdr`. These are separate models and adjustments.
* Family survival figures use the actual fixed median groups: Low <= median,
  High > median. `hr_high_vs_low` and `group_ci_lower`/`group_ci_upper` are the
  displayed grouped estimates. Continuous and per-unit estimates remain in
  explicitly named additional columns.
* Alteration heatmap gene/cancer cells count events, including different
  alteration types in the same sample. `Total` counts unique samples carrying
  at least one alteration among the custom input genes.
  Unprofiled cells and cancer totals remain unavailable, rather than zero.
  The source's combined COADREAD study is not assigned to COAD or READ without
  a separate cancer annotation.

## Explicit compatibility differences

The adapters do not move observations between score groups or reverse group
labels. Source routines containing such steps are not used to create the
modular survival figures. Missing or unestimable values remain marked rather
than being converted into evidence of no association.

ROC estimates retain a prespecified higher-score-predicts-tumor direction and
the empirical curve. They are not smoothed, so they can differ from the
inherited smoothed ROC output; an AUC below 0.5 is retained.

Unknown GEO follow-up units are not labelled as years. CNV/methylation values
and FDR are taken from the available source tables. Custom gene coverage is
checked before use, and undefined data are not synthesized to complete a plot.

Device handling saves a vector PDF and a white-background 300-dpi PNG. Figure
dimensions adapt when necessary for custom family sizes and labels.
