# UF-HF Comorbidity Single-Cell Analysis — Reproducible Code

## Overview
This repository contains the analysis code for the study "CPQ-THBS axis links
ultrafiltration failure and heart failure via fibroblast/mesothelial remodeling".
The pipeline spans single-cell QC, Scissor phenotype projection, machine-learning
feature selection, cell-cell communication, pySCENIC GRN inference, and in-silico
perturbation.

## Software Versions
| Tool | Version |
|---|---|
| R | 4.4.3 |
| Python | 3.10 |
| Seurat | 5.4.0 |
| Scissor | 1.0.0 |
| Harmony | 1.2.4 |
| SingleR | 2.14.2 |
| scDblFinder | 1.25.0 |
| ClusterGVis | 0.99.9 |
| pySCENIC | 0.12.1 |
| CellChat | (current) |
| Geneformer | V2 (104M) |
| scTenifoldKnk | (current) |
| CellOracle | (current) |

## Script Execution Order
| Step | Script | Description |
|---|---|---|
| 01 | `01_qc_harmony_annotation.R` | QC (MAD 3×, percent.mt<15%, scDblFinder), LogNormalize, 2000 HVG, PCA 20, Harmony, UMAP, Louvain (res=0.25), SingleR (HPCA) + marker annotation |
| 02 | `02_clustergvis_modules.R` | ClusterGVis k-means modules across Monocle states S1–S6 + GO enrichment |
| 03 | `03_scissor_UF.R` / `03_scissor_HF.R` | Scissor phenotype projection (alpha=0.01, binomial): UF (long/short → 1/0) / HF (ICM+DCM → 1, Normal → 0) |
| 04 | `04_scissor_DEG_intersection.R` | Scissor+ vs Scissor− DEG (p.adj<0.05, log2FC>2); UF∩HF → 418 comorbid candidates |
| 05 | `05_ml_lasso_rf_svmrfe_shap.R` | LASSO (10-fold), RF (ntree=500, top30 Gini), SVM-RFE (10-fold, top10); intersection; SHAP (shapr + ranger probability forest) |
| 06 | `06_nomogram_calibration_dca.R` | 7-gene logistic nomogram (regplot), 1000× Bootstrap calibration, DCA, external validation (GSE116250/GSE141910) |
| 07 | `07_cellchat_CPQ.R` | CPQ High/Low CellChat in UF (Fibroblast, M-MCs, T_cell) and HF (Fibroblast); NMF patterns; THBS axis |
| 08 | `08_fibrosis_correlation_pseudobulk.R` | CPQ–fibrosis Spearman (5 HF cohorts), Fibrosis signature Pearson, pseudobulk E-MCs vs M-MCs Wilcoxon |
| 09 | `09_pyscenic_grn.py` | GRNBoost2 (multiprocessing) → cisTarget hg38 → AUCell → differential TF (Wilcoxon + BH, ΔMean>0) |
| 10a | `10a_geneformer.py` | Geneformer V2 embeddings + L2 logistic classifier (7:3, F1/confusion matrix) |
| 10b | `10b_sctenifoldknk.R` | Virtual KO (CPQ/EBF2/GLI3/THBS4 in HF; CPQ/HIF1A/SATB1/THBS1 in UF), q<0.01 & |Z|>1.8 |
| 10c | `10c_celloracle.py` | CellOracle GRN + virtual perturbation (HIF1A UF; RUNX1/RUNX2 HF; Random KO control) |
| 11 | `11_mediation_CPQ_THBS4_HF.R` | Mediation: X=CPQ, M=THBS4, Y=HF; Z-score; 1000× Bootstrap; GSE57338 discovery, GSE116250 validation |
| 12 | `12_statistics_utils.R` | Shared Wilcoxon/Pearson/Spearman/BH helpers + `theme_publication()` |

## Input Files Required
| File | Used by | Notes |
|---|---|---|
| `GSE248762_processed_v2.rds` | 01–08, 10b | Seurat v5 object, annotated |
| `GSE298023_processed_v2.rds` | 01–08, 10b | Seurat v5 object, annotated |
| `GSE125498_filtered_gene_level_expression_matrix.txt` | 03_UF, 04 | Bulk, genes × samples |
| `clinical_data.csv` | 03_UF | Columns: `Sample`, `Group` (`long`/`short`) |
| `GSE57338_expression_matrix.txt` + clinical | 03_HF, 05, 06, 11 | Discovery HF cohort |
| `GSE116250_expression_matrix.txt` + clinical | 06, 11 | Validation |
| `GSE141910_expression_matrix.txt` + clinical | 06, 08 | Validation |
| `fibrosis_signature.gmt` | 08 | **User-supplied fibrosis gene set** |
| `allTFs_hg38.txt` | 09 | Download: `https://resources.aertslab.org/cistarget/tf_lists/allTFs_hg38.txt` |
| cisTarget `.feather` databases (10kb, 500bp) | 09 | Download from `https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/` |
| `motifs-v9-nr.hgnc-m0.001-o0.0.tbl` | 09 | Download from `https://resources.aertslab.org/cistarget/motif2tf/` |
| Geneformer-V2-104M | 10a | Clone `https://huggingface.co/ctheodoris/Geneformer` |
| CellOracle human base GRN | 10c | Built-in via `co.data.load_human_promoter_base_GRN()` |

## Random Seeds
All scripts use `set.seed(1234)` (R) or `seed=1234` (Python) unless otherwise noted.
Bootstrap procedures (calibration, mediation) use 1000 resamples.

## Key Parameter Notes
- **Scissor alpha = 0.01**: balance between penalty strength and recall (no established default).
- **SingleR reference**: `celldex::HumanPrimaryCellAtlasData()` — chosen for broader mesenchymal/mesothelial coverage (713 primary samples, 157 subtypes) vs BlueprintEncode (259 samples).
- **pySCENIC GRN step**: uses `arboreto_with_multiprocessing.py` to avoid Dask version conflicts (official recommendation).
- **CellChat**: run independently per dataset (UF and HF); NMF `k=3` for outgoing/incoming patterns.
- **Fibrosis signature score**: gene set supplied by user (`fibrosis_signature.gmt`); method description must state the exact gene list.

## Methodological Notes for Manuscript
1. **Peritoneal injury score**: computed via ssGSEA/GSVA on GSE125498 using a user-defined gene set. **Add a sentence to Methods** specifying the gene set used.
2. **SHAP**: implemented in R via `shapr` + `ranger(probability=TRUE)`; "probability forest" is the correct term (not standard RF).
3. **S1–S6**: Monocle pseudotime states; S1 = early, S6 = Terminal. Confirm state numbering matches biological direction before running `02`.

## Contact
For questions, please open an issue or contact the corresponding author.