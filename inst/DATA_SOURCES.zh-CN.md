# PanFamR 数据来源、处理边界与参考文献

核对日期：2026-10-15。适用版本：PanFamR 0.3.1。

本说明将**数据生成项目或原始研究**、**整理获取平台**和**分析方法文献**分别记录。数据库论文可以用于致谢和方法引用，但不能单独证明某个本地文件的下载版本、处理流程或再分发条件已经核实。

## 1. 当前已识别的资源

- **TCGA / TCGA PanCancer Atlas**：泛癌肿瘤及正常表达、临床结局和 TCGA 样本辅助表。依据包括 TCGA 样本条码、分组标签及变异文件中的 PanCancer Atlas study ID。原始整理文件大多未附具体下载门户或批次。[TCGA 官方项目](https://www.cancer.gov/ccg/research/genome-sequencing/tcga)、[PanCancer Atlas 数据入口](https://gdc.cancer.gov/about-data/publications/pancanatlas)。研究背景见文献[1]。
- **GTEx**：TCGA/GTEx 联合表达输入中的正常组织参考，文件明确保留 `GTEx_normal` 标签。具体 GTEx 发布版及联合矩阵处理流程尚未确认。文献[2]为 GTEx 项目参考，不能据此认定本地矩阵采用 v8。[GTEx Portal](https://www.gtexportal.org/home/)。
- **GEO**：GSE78220、GSE91061、GSE135222 三个可核对的公开系列。数据库参考见[3]，还需引用具体系列和原始研究[5–8]。
- **TIGER**：维护者于 2026-10-15 确认，免疫治疗表达/临床 RDS 整理文件经 TIGER 获取。平台参考见[4]。这是维护者补充的获取记录，原下载日期和本地转换脚本仍待完善。[TIGER 论文与入口说明](https://pmc.ncbi.nlm.nih.gov/articles/PMC10626175/)。**GSE135222 的默认计算另使用原生表达 TSV 与官方 GEO 样本映射，未使用保留的 TIGER 表达 RDS。**
- **MSigDB Hallmark**：源码和原工作目录中的 `h.all.v2025.1.Hs.symbols.gmt` 对应 MSigDB 的 Hallmark 集合；GMT 内含官方基因集链接。包内提供的是 50 项预计算通路评分，不分发 GMT。参考[10,11]。RData 未携带生成记录，因此未独立复核该评分文件的完整生成版本链。[MSigDB 官方引用说明](https://www.gsea-msigdb.org/gsea/msigdb)。
- **GENCODE v22 派生注释**：`gencode.gene.info.v22.tsv` 用于 Ensembl ID 与基因符号对应，参考[12]及 [GENCODE release 22 官方页面](https://www.gencodegenes.org/human/release_22)。尚未与官方原始文件校验和逐一匹配。

MSigDB 在本项目中是**基因集资源**。目前没有证据将患者级 TMB/MSI 或 CNV/甲基化相关性表归为 MSigDB 数据；这些表的实际下载渠道继续列为待核实。

## 2. 文件与出处对应表

下表中的文件名相对于安装包 `extdata` 目录。`.partNN` 是同一文件的存储分片，不是新的数据来源。

| 包内输入 / 模块 | 可确认的来源或证据 | 尚待补充的记录 |
| --- | --- | --- |
| `expression/all_tpm.Rdata.partNN`；ssGSEA、表达、相关分析 | 作者整理的 TCGA 表达输入；TCGA 条码及癌种/组别标签，11,236 个样本[1] | 下载门户、原始文件名和日期、定量流程、表达尺度及转换、注释版本。文件名含 tpm 不替代这些记录。 |
| `expression/TCGA_GTEx_pancancer_mrna_pheno.rdata.partNN`；正常参考 | 18,064 个样本，标签明确含 TCGA 肿瘤、TCGA 正常和 GTEx 正常[1,2] | GTEx 发布版、整合门户、归一化/组织匹配。尚不能确认为 UCSC Xena/Toil 文件。 |
| `expression/paired_tpm.Rdata.partNN` | 原项目保留的 TCGA 配对表达输入；模块化配对分析按患者条码匹配已计算评分 | 原文件下载及配对构建过程。 |
| `clinical/ALL_OS.Rdata`、`clinical/DSS_PFI_DFI.Rdata` | TCGA 样本的 OS、DSS、PFI、DFI 整理表；终点研究背景可参阅[16] | 原下载表、版本及整理脚本；仅凭终点名不能认定直接来自 TCGA-CDR。 |
| `auxiliary/All_TMB.Rdata` | 9,649 个 TCGA 样本的 TMB 辅助表 | 原始来源、算法、变异范围、分母、单位及版本。不得自行标注为某种指定 TMB 算法。 |
| `auxiliary/All_MSI.Rdata` | 9,202 个 TCGA 标识对应的连续 MSI 分值 | 原始论文/门户、算法、版本及阈值；不推定为 MANTIS 或 MSIsensor。 |
| `auxiliary/alterations_across_samples_*.tsv` | 7 份基因限制表含 `Study ID`，例如 `laml_tcga_pan_can_atlas_2018`，与 TCGA PanCancer Atlas 研究命名对应 | 实际导出平台、日期、病例与分子谱。cBioPortal 官方 DataHub 可核对研究 ID，但不等于已确认本地文件的导出渠道。 |
| `family_tables/*/CnvAndExpressionTable.xlsx` | 作者提供的 CNV–表达预计算相关系数/FDR 表 | 原下载平台及分析版本；GSCA 是待核对线索，当前不作为已确认来源。 |
| `family_tables/*/ExpressionAndMethylationTable.xlsx` | 作者提供的甲基化–表达预计算相关系数/FDR 表；包含探针标识 | 来源平台、甲基化平台、探针处理及 FDR 范围；不能归为 MSigDB。 |
| `family_tables` 下 SNV/CNV 汇总、差异甲基化 XLSX | 原项目保留的家族辅助表；当前 CNV/甲基化模块使用上述两个相关性表类型 | 原始来源和处理说明。当前包内共有 19 份家族 XLSX，其核心属性未给出来源链接。 |
| `auxiliary/ssgsea_immucells.Rdata` | 10,496×24 的预计算免疫细胞相关评分 | 基因签名来源、生成软件和参数。列名与常见 24 细胞面板相似，不足以确认 Bindea 原始签名，也不应将评分直接称为细胞比例。 |
| `auxiliary/immuscore.Rdata` | `StromalScore`、`ImmuneScore`、`ESTIMATEScore` 三列对应 ESTIMATE 体系，方法参考[15] | 实际下载来源、软件版本、输入与计算记录。方法引用不替代文件溯源。 |
| `auxiliary/TME_Pathaway_GSVA_Score.Rdata` | 17 项预计算评分；原工作目录的 `TME相关通路.txt` 有对应逐基因定义，源码保留 GSVA 计算框架[14] | 各签名最初整理来源及文献，不能将 17 项统一归因某篇论文。原定义文件尚未随包分发。 |
| `auxiliary/HALLMARK_GSVA_Score.Rdata` | 11,236×50；原脚本及原目录 GMT 对应 MSigDB Hallmark，原文件名指定 2025.1.Hs[10,11,14] | 当前 RData 的生成会话、软件版本和 GMT 校验和。 |
| `auxiliary/gencode.gene.info.v22.tsv` | 作者提供的 GENCODE v22 派生基因注释[12] | 派生脚本与官方原始文件校验和。 |
| 包代码中的免疫检查点、趋化因子、受体、刺激/抑制因子列表 | 沿用作者原脚本的基因面板 | 外部整理依据和引用；不据列表相似度将其直接归因 TISIDB 或其他数据库。 |

变异研究 ID 核对例：[cBioPortal DataHub 官方 LAML 元数据](https://raw.githubusercontent.com/cBioPortal/datahub/master/public/laml_tcga_pan_can_atlas_2018/meta_study.txt)。CNV/甲基化下载渠道核对线索：[GSCA 官方研究论文](https://doi.org/10.1093/bib/bbac558)。二者在此用于说明核对线索，**不表示本地文件已完成下载链匹配**。

## 3. 免疫治疗队列：平台与原始研究均应引用

| 本地队列 | 获取与原始研究 | 默认分析中的处理 |
| --- | --- | --- |
| GSE78220 | 维护者确认 RDS 经 TIGER 获取[4]；[GEO GSE78220](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78220)；Hugo 等原始研究[5] | 本地 28 条临床记录中 PRE 27、ON 1；优先保留治疗前样本，不能写为 28 个治疗前样本。 |
| GSE91061 | 维护者确认 RDS 经 TIGER 获取[4]；[GEO GSE91061](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE91061)；Riaz 等[6] | 原系列为 65 名患者的 109 个 RNA-seq 样本，其中治疗前 51、治疗中 58；实际纳入需看本包筛选表，不能把 109 写作独立患者数。 |
| GSE135222 | 原 RDS 经 TIGER 获取[4]并保留；默认表达采用 `GSE135222_GEO_RNA-seq_omicslab_exp.tsv`，样本映射及 PFS 注释依据[官方 GEO 系列](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE135222)；原始关联研究[7,8] | GEO 公开 27 名患者、27 个 RNA-seq 样本，本包按可核对样本映射使用。Kim2020 部分段落写 28，但 Methods 与 GEO 为 27，不强行统一。GEO 未注明 `pfs.time` 单位，保留来源单位。 |
| Braun2020 | 维护者确认 RDS 经 TIGER 获取[4]；Braun 等透明细胞肾癌研究[9]；不标为 GEO 队列 | 原研究整合 CheckMate-009/010/025。本地 311 条临床记录含抗 PD-1 181、everolimus 130；流程排除已知非免疫治疗对照，后续有效样本仍取决于表达和随访覆盖。 |

GSE135222 官方 series matrix 在 2026-10-14 获取，原文件与派生 `GSE135222_sample_metadata.tsv` 一并保留：
[NCBI 官方 series matrix](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135222/matrix/GSE135222_series_matrix.txt.gz)。其他既有 RDS 的历史下载日期未记录，不能把本次核对日期当作下载日期。

## 4. 方法文献与数据权利说明

- 家族评分使用 GSVA 包实现的 ssGSEA，方法参考[13,14]；Hallmark/MSigDB 另参考[10,11]。
- ESTIMATE[15]及 TCGA-CDR[16]在记录不完整处仅提供方法或终点背景，不能据此宣称当前预计算表的完整生成流程已经重现。
- 代码采用 MIT 许可证。第三方数据、签名及派生文件仍应遵循各自适用的条款和署名要求；新增引文不是对全部再分发条件的确认。
- GEO 建议同时引用数据库、系列编号及原研究：[官方引用说明](https://www.ncbi.nlm.nih.gov/geo/info/linking.html)。其[数据声明](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html)也保留了提交者可能主张的相关权利。
- MSigDB 的[官方许可说明](https://www.gsea-msigdb.org/gsea/license_terms_list.jsp)区分版本，并可能涉及附加条款。包内未分发 GMT；预计算评分的输入、引用和相关使用条件仍需记录。
- TCGA、GTEx、TIGER 等应按实际使用的数据类别、访问方式与原研究说明处理。本次未据论文开放获取状态推定所有底层个体数据均可任意再分发。

## 5. 待维护者补齐的最小记录

1. TMB/MSI 的原下载网页或论文附件、算法、数值单位及版本。
2. CNV/甲基化 XLSX 与变异 TSV 的实际导出网址、日期、所选病例/分子谱和处理参数。
3. 表达和临床文件的原始下载文件名、门户/版本、尺度及整理脚本。
4. 24 类免疫细胞签名、ESTIMATE 评分的生成记录，以及 17 项 TME 定义各自的外部引用。
5. TIGER 原下载文件清单和日期；所有再分发文件适用的原始条款或授权记录。

## 6. 参考文献

1. The Cancer Genome Atlas Research Network, Weinstein JN, Collisson EA, et al. The Cancer Genome Atlas Pan-Cancer analysis project. *Nature Genetics*. 2013;45:1113–1120. [DOI: 10.1038/ng.2764](https://doi.org/10.1038/ng.2764).

2. GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. *Science*. 2020;369:1318–1330. [DOI: 10.1126/science.aaz1776](https://doi.org/10.1126/science.aaz1776).

3. Barrett T, Wilhite SE, Ledoux P, et al. NCBI GEO: archive for functional genomics data sets—update. *Nucleic Acids Research*. 2013;41:D991–D995. [DOI: 10.1093/nar/gks1193](https://doi.org/10.1093/nar/gks1193).

4. Chen Z, Luo Z, Zhang D, et al. TIGER: A Web Portal of Tumor Immunotherapy Gene Expression Resource. *Genomics, Proteomics & Bioinformatics*. 2023;21:337–348. [DOI: 10.1016/j.gpb.2022.08.004](https://doi.org/10.1016/j.gpb.2022.08.004).

5. Hugo W, Zaretsky JM, Sun L, et al. Genomic and Transcriptomic Features of Response to Anti-PD-1 Therapy in Metastatic Melanoma. *Cell*. 2016;165:35–44. [DOI: 10.1016/j.cell.2016.02.065](https://doi.org/10.1016/j.cell.2016.02.065).

6. Riaz N, Havel JJ, Makarov V, et al. Tumor and Microenvironment Evolution during Immunotherapy with Nivolumab. *Cell*. 2017;171:934–949.e16. [DOI: 10.1016/j.cell.2017.09.028](https://doi.org/10.1016/j.cell.2017.09.028).

7. Jung H, Kim HS, Kim JY, et al. DNA methylation loss promotes immune evasion of tumours with high mutation and copy number load. *Nature Communications*. 2019;10:4278. [DOI: 10.1038/s41467-019-12159-9](https://doi.org/10.1038/s41467-019-12159-9).

8. Kim JY, Choi JK, Jung H. Genome-wide methylation patterns predict clinical benefit of immunotherapy in lung cancer. *Clinical Epigenetics*. 2020;12:119. [DOI: 10.1186/s13148-020-00907-4](https://doi.org/10.1186/s13148-020-00907-4).

9. Braun DA, Hou Y, Bakouny Z, et al. Interplay of somatic alterations and immune infiltration modulates response to PD-1 blockade in advanced clear cell renal cell carcinoma. *Nature Medicine*. 2020;26:909–918. [DOI: 10.1038/s41591-020-0839-y](https://doi.org/10.1038/s41591-020-0839-y).

10. Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles. *Proceedings of the National Academy of Sciences USA*. 2005;102:15545–15550. [DOI: 10.1073/pnas.0506580102](https://doi.org/10.1073/pnas.0506580102).

11. Liberzon A, Birger C, Thorvaldsdóttir H, et al. The Molecular Signatures Database (MSigDB) hallmark gene set collection. *Cell Systems*. 2015;1:417–425. [DOI: 10.1016/j.cels.2015.12.004](https://doi.org/10.1016/j.cels.2015.12.004).

12. Harrow J, Frankish A, Gonzalez JM, et al. GENCODE: The reference human genome annotation for The ENCODE Project. *Genome Research*. 2012;22:1760–1774. [DOI: 10.1101/gr.135350.111](https://doi.org/10.1101/gr.135350.111).

13. Barbie DA, Tamayo P, Boehm JS, et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. *Nature*. 2009;462:108–112. [DOI: 10.1038/nature08460](https://doi.org/10.1038/nature08460).

14. Hänzelmann S, Castelo R, Guinney J. GSVA: gene set variation analysis for microarray and RNA-Seq data. *BMC Bioinformatics*. 2013;14:7. [DOI: 10.1186/1471-2105-14-7](https://doi.org/10.1186/1471-2105-14-7).

15. Yoshihara K, Shahmoradgoli M, Martínez E, et al. Inferring tumour purity and stromal and immune cell admixture from expression data. *Nature Communications*. 2013;4:2612. [DOI: 10.1038/ncomms3612](https://doi.org/10.1038/ncomms3612).

16. Liu J, Lichtenberg T, Hoadley KA, et al. An Integrated TCGA Pan-Cancer Clinical Data Resource to Drive High-Quality Survival Outcome Analytics. *Cell*. 2018;173:400–416.e11. [DOI: 10.1016/j.cell.2018.02.052](https://doi.org/10.1016/j.cell.2018.02.052).
