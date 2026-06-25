# PHLDA3 in Thyroid Carcinoma (THCA)

Reproducible code, figures and manuscript for:

> **An integrated multi-omics and single-cell bioinformatics analysis PHLDA3 marks aggressive driver-defined thyroid carcinoma: an integrated multi-omics and single-cell analysis.**

PHLDA3 (a p53 target / AKT repressor) is induced cell-intrinsically in malignant thyroid cells — associated with promoter hypomethylation within a BRAF^V600E^/p53-driven, checkpoint-high immune context — and is comparably elevated in BRAF^V600E^ and RET/NTRK-fusion-positive tumours; its association with lymph-node metastasis is largely shared with this aggressive driver landscape rather than driver-independent.

## Repository layout

```
data/         # downloaded raw data + intermediate caches (large raw files git-ignored)
code/         # analysis scripts (Python + R) and the shared figure theme (nature_style.py)
figures/      # result figures (PDF)
manuscript/   # PHLDA3_THCA_manuscript.md
README.md · requirements.txt · .gitignore
```

**Run scripts from the repository root**, e.g. `python3 code/01_volcano.py` or `Rscript code/04_run_estimate.R`. Scripts read inputs from `data/` and write figures to `figures/`.

## Environment

Python 3.13.11 (`pip install -r requirements.txt`) and R 4.5 (packages: `estimate`; optional `infercnv`/`CellChat`). Key Python versions: pandas 2.3.3, numpy 2.4.4, scipy 1.17.1, matplotlib 3.10.9, scikit-learn 1.8.0, statsmodels 0.14.6, lifelines 0.30.3, gseapy 1.3.0, pingouin 0.6.1, xenaPython 1.0.14, cellxgene-census 1.18.0.

## Data (large raw files not tracked — download into `data/`)

| File / source | Used by | Download |
|---|---|---|
| `TCGA-THCA.htseq_counts.tsv`, `THCA_clinical.tsv` | most scripts | GDC / UCSC Xena |
| `TCGA-THCA.methylation450.tsv.gz`, `TCGA-THCA.mirna.tsv.gz` | methylation / miRNA | UCSC Xena (GDC hub) |
| `GSE33630.txt.gz` | `geo_validation_thca.py` | [GEO GSE33630](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33630) |
| `gse193581_raw/` + `gse193581_anno.txt.gz` | tumour single-cell | [GEO GSE193581](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193581) |
| `gdsc2/*.rds` | `run_gdsc2.R` | oncoPredict data (OSF c6tfx) |
| pan-cancer / mutations / census / ENCORI / L1000CDS2 | auto-fetch | UCSC Xena · cBioPortal · CELLxGENE · ENCORI · L1000CDS2 |

Tracked intermediate caches (so cache-based figures reproduce without large downloads): `data/thca_DE_full_symbol.csv`, `data/pancancer_phlda3_raw.csv`, `data/thyroid_sc_phlda3.csv`, `data/PHLDA3_tumor_sc.csv`, `data/PHLDA3_cmap_moa_clean.csv`, `data/gdsc2_phlda3_drugcorr.csv`, `data/cg04055835_betas.csv`, `data/thca_fusion_status.csv`, etc.

## Analyses (figure → script)

Scripts are numbered by manuscript figure (`0N_` = main Figure N; `SN_` = Supplementary Figure N). Run from the repository root.

| Manuscript figure | Scripts |
|---|---|
| Fig 1 — expression / validation / pan-cancer | `01_volcano.py`, `01_expression_roc.py`, `01_pancancer.py`, `01_geo_validation.py` |
| Fig 2 — driver-aware N1 (nested BRAF+fusion models, driver classes, incremental AUC) | `02_figure2_composite.py` (composite) ← `02_logistic_n1.py`, `02_braf_sensitivity.py`, `02_braf_stratified.py`, `02_expanded_n1_fusion.py`, `02_driver_robustness.py`, `02_nomogram_dca.py` |
| Fig 3 — driver landscape (4-class: expression, N1 frequency, promoter methylation) | `03_driver_landscape.py` (composite) ← `03_braf_ras.py` (driver calls) |
| Fig 4 — immune + tumour-purity correction | `04_immune.py`, `04_purity_partialcorr.py` (+ `04_run_estimate.R`) |
| Fig 5 — single cell (tumour, p53 module, cell-cell communication) | `05_figure5_composite.py` (composite) ← `05_sc_tumor.py`, `05_sc_p53_module.py`, `05_sc_sample_level.py`, `05_sc_combined.py`, `05_sc_liana.py` |
| S1 — clinicopathology | `S1_clinical.py` |
| S2 — mechanism (CNV + methylation) | `S2_cnv.py`, `S2_methylation.py` |
| S3 — co-expression GO/KEGG/GSEA | `S3_enrichment.py` |
| S4 — miRNA / ceRNA (exploratory) | `S4_mirna.py`, `S4_cerna.py` |
| S5 — drug repurposing | `S5_l1000cds2.py`, `S5_run_gdsc2.R` + `S5_gdsc2.py`, `S5_cmap.py` |
| S6 — protein (HPA) | `S6_hpa.py` |

Upstream helpers: `00_de_screen.py` (DE screen), `00_map_symbols.py` (Ensembl→symbol), `00_sc_normal_fetch.py` (normal-thyroid scRNA fetch). All figures share `code/nature_style.py` (Arial, unified palette, PDF export). Independent parallel cross-validation outputs (purity, single-cell incl. LIANA, immune) are kept under `results/`.

## Data availability / attribution

TCGA (GDC), GEO, UCSC Xena, cBioPortal, Human Protein Atlas, ENCORI, L1000CDS2, GDSC2 (oncoPredict) and CELLxGENE Census are publicly available; cite the original data producers and methods papers listed in the manuscript.
