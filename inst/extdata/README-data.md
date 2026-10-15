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

## Sources and references

See [the detailed provenance and reference record](../DATA_SOURCES.zh-CN.md)
for the file-by-file audit, original cohort accessions, 16 references, and
remaining gaps (reviewed 2026-10-15). The identifiable resources include TCGA,
GTEx, GEO, TIGER, MSigDB Hallmark, and the supplied GENCODE v22-derived
annotation. The maintainer confirmed TIGER as the source of the supplied
immunotherapy RDS collection; the GSE135222 default input separately uses
the native expression TSV and official GEO sample mapping described above.

TMB/MSI download sources and algorithms, the actual CNV/methylation export
portal, and some precomputed immune-score generation records remain
unconfirmed. These patient-level tables must not be attributed to MSigDB
merely because the original workspace also contains Hallmark signatures.
Database and method references do not establish a missing file-processing
chain or confirm all redistribution conditions. Original download dates and
release versions were not available for many supplied files.

The MIT license applies to package code; it does not grant rights to
redistribute third-party research data. Consult the applicable source terms
and attribution requirements. MSigDB GMT collections are intentionally not
bundled; users may pass their own gene sets to `panfam_analyze()`.
