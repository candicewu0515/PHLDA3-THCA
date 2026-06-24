# PHLDA3 in Thyroid Carcinoma (THCA)

Reproducible analysis code and figures for the study:

> **PHLDA3 is induced in thyroid carcinoma cells and is an independent predictor of lymph-node metastasis and advanced stage: an integrated multi-cohort and single-cell analysis.**

PHLDA3 — a p53 target and AKT repressor — is over-expressed in thyroid carcinoma, discriminates tumour from normal tissue in two independent cohorts, independently predicts lymph-node metastasis and advanced stage, tracks with the BRAF^V600E^ subtype and an immunosuppressive microenvironment, and is induced cell-intrinsically in malignant thyroid cells at single-cell resolution.

The full manuscript draft is in [`PHLDA3_THCA_manuscript.md`](PHLDA3_THCA_manuscript.md).

## Environment

Python 3.13.11. Install dependencies:

```bash
pip install -r requirements.txt
```

Key versions: pandas 2.3.3, numpy 2.4.4, scipy 1.17.1, matplotlib 3.10.9, scikit-learn 1.8.0, statsmodels 0.14.6, lifelines 0.30.3, gseapy 1.3.0, cellxgene-census 1.18.0, xenaPython 1.0.14.

## Data sources

Large raw files are **not** tracked in this repository (size / public availability). Download them into the repository root before running the from-scratch scripts:

| File / source | Used by | Download |
|---|---|---|
| `TCGA-THCA.htseq_counts.tsv` (HTSeq counts) | most scripts | GDC / UCSC Xena (TCGA-THCA) |
| `THCA_clinical.tsv` (tracked) | clinical scripts | GDC |
| `GSE33630.txt.gz` (series matrix) | `geo_validation_thca.py` | [GEO GSE33630](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33630) |
| `gse193581_raw/` (per-sample UMI matrices from `GSE193581_RAW.tar`) | `sc_tumor_phlda3.py` | [GEO GSE193581](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193581) |
| Pan-cancer TCGA+GTEx | `pancancer_phlda3.py` (auto-fetch) | UCSC Xena Toil hub (via `xenaPython`) |
| Somatic mutations | `braf_ras_thca.py` (auto-fetch) | cBioPortal API (`thca_tcga_pan_can_atlas_2018`) |
| Normal-thyroid scRNA | `sc_normal_census.py` (auto-fetch) | CELLxGENE Census 2025-11-08 |
| Protein IHC images (tracked) | `hpa_ihc_thca.py` | Human Protein Atlas (ENSG00000174307) |

Intermediate caches that **are** tracked (so cache-based figures reproduce without large downloads): `thca_DE_full.csv`, `thca_DE_full_symbol.csv`, `pancancer_phlda3_raw.csv`, `thyroid_sc_phlda3.csv`, `PHLDA3_tumor_sc.csv`, and all result tables.

## Figures → scripts

| Figure | Script | Theme/style |
|---|---|---|
| 1 — diagnostic (volcano, tumour-vs-normal ROC, pan-cancer) | `volcano_thca.py`, `build_phlda3.py`, `pancancer_phlda3.py` | |
| 2 — external + protein validation | `geo_validation_thca.py`, `hpa_ihc_thca.py` | |
| 3 — clinicopathological associations | `clinical_phlda3.py` | |
| 4 — logistic regression + nomogram/DCA | `logistic_phlda3.py`, `nomogram_dca_thca.py` | |
| 5 — BRAF/RAS driver mutations | `braf_ras_thca.py` | |
| 6 — co-expression GO/KEGG/GSEA | `coexpr_enrich_thca.py` | |
| 7 — immune infiltration + checkpoints | `immune_thca.py` | |
| 8 — single-cell (normal vs tumour) | `combined_singlecell_nature.py` | Nature-style |

All figures share one palette/theme defined in [`nature_style.py`](nature_style.py) and are exported as SVG (editable), PDF, 600-dpi TIFF and PNG.

Upstream helpers: `screen_thca.py` (tumour-vs-normal DE screen → `thca_DE_full.csv`), `map_de_symbols.py` (Ensembl→symbol), `sc_normal_census.py` (normal-thyroid scRNA fetch), `sc_tumor_phlda3.py` (tumour scRNA → `PHLDA3_tumor_sc.csv`).

## Reproduce (cache-based figures, no large downloads)

```bash
pip install -r requirements.txt
python3 volcano_thca.py             # Fig 1a
python3 pancancer_phlda3.py         # Fig 1c   (auto-fetches via xenaPython)
python3 combined_singlecell_nature.py   # Fig 8
```

For the remaining figures, download the raw files listed above, then run the corresponding scripts.

## Data availability / attribution

TCGA (GDC), GEO (GSE33630, GSE193581), UCSC Xena, cBioPortal, Human Protein Atlas, and CZI CELLxGENE Census are publicly available; please cite the original data producers and methods papers listed in the manuscript references.
