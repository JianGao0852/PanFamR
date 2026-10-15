# Bundled data notes

The files in this directory were supplied in the author's analysis workspace
and copied into PanFamR so users can reproduce the bundled pan-cancer
workflows. Core files include TCGA/GTEx expression matrices, TCGA survival
tables, TMB/MSI tables, immune-score tables, alteration summaries, and public
immunotherapy cohort inputs. `panfam_data_files()` lists every installed file.

GSE135222 uses the supplied native expression TSV with GENCODE v22 symbol
annotation. Its sample mapping and PFS indicator/time fields were prepared
from the official GEO series matrix retrieved on 2026-10-14:
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz
The original series matrix is bundled alongside the derived sample metadata.
GEO does not specify the unit of `pfs.time` in this annotation; it is retained
without rescaling and plotted as source units. The provided NSCLC response
RDS is retained for provenance but is not used as the expression input in the
default workflow because it contains extensive missing values.

The default immunotherapy workflow selects baseline samples where annotated,
excludes known non-immunotherapy controls, and uses one sample per patient.
Unknown response/event labels remain missing. Duplicate gene-symbol rows are
averaged; genes with non-finite values in selected samples are excluded and
recorded, without imputing expression. Constant rows retained by GSVA for
ranking are recorded in the cohort filtering table.

Mutation and CNV/methylation inputs are gene-restricted source tables; they do
not provide genome-wide coverage. The workflow reports gene coverage and
skips unavailable associations rather than treating them as negative results.

The bundled inputs occupy approximately 1.29 GiB. Large inputs
`all_tpm.Rdata`, `TCGA_GTEx_pancancer_mrna_pheno.rdata`, `paired_tpm.Rdata`,
and `RCC-Braun_2020.Response.Rds` are split into numbered `.partNN` pieces
of at most 20 MiB. PanFamR automatically reconstructs all split inputs in
the user's R cache directory when the corresponding dataset is first needed.

The source workspace did not include a complete provenance and license record
for every table. Before making a public GitHub release, the maintainer must
record the original source, citation, and redistribution terms for each
third-party dataset. The MIT license in the package applies to package code;
it does not grant rights to redistribute the bundled research data. MSigDB
GMT collections are intentionally not bundled; users may pass their own gene
sets to `panfam_analyze()`.
