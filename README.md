# PHLDA3 in Thyroid Carcinoma (THCA)

Reproducible code, figures and manuscript for:

> **An integrated multi-omics and single-cell bioinformatics analysis identifies PHLDA3 as a BRAF-associated transcriptomic marker linked to lymph-node metastasis in thyroid carcinoma.**

PHLDA3 (a p53 target / AKT repressor) is induced cell-intrinsically in malignant thyroid cells — associated with promoter hypomethylation within a BRAF^V600E^/p53-driven, immunosuppressive context — and independently correlates with lymph-node metastasis (robust to BRAF adjustment).

## Repository layout

```
data/         # downloaded raw data + intermediate caches (large raw files git-ignored)
code/         # analysis scripts (Python + R) and the shared figure theme (nature_style.py)
figures/      # result figures (PDF)
manuscript/   # PHLDA3_THCA_manuscript.md
README.md · requirements.txt · .gitignore
```

**Run scripts from the repository root**, e.g. `python3 code/volcano_thca.py` or `Rscript code/run_estimate.R`. Scripts read inputs from `data/` and write figures to `figures/`.

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

Tracked intermediate caches (so cache-based figures reproduce without large downloads): `data/thca_DE_full_symbol.csv`, `data/pancancer_phlda3_raw.csv`, `data/thyroid_sc_phlda3.csv`, `data/PHLDA3_tumor_sc.csv`, `data/PHLDA3_cmap_moa_clean.csv`, `data/gdsc2_phlda3_drugcorr.csv`, etc.

## Analyses (figure → script)

| Manuscript figure | Script |
|---|---|
| Expression / validation / pan-cancer | `volcano_thca.py`, `build_phlda3.py`, `pancancer_phlda3.py`, `geo_validation_thca.py` |
| N1 logistic + BRAF sensitivity + nomogram | `logistic_phlda3.py`, `braf_sensitivity_thca.py`, `nomogram_dca_thca.py` |
| BRAF/RAS drivers | `braf_ras_thca.py` |
| Mechanism (CNV + methylation) | `cnv_phlda3.py`, `meth_phlda3.py` |
| Co-expression GO/KEGG/GSEA | `coexpr_enrich_thca.py` |
| Immune + purity correction | `immune_thca.py`, `purity_partialcorr_thca.py` (+ `run_estimate.R`) |
| miRNA / ceRNA | `mirna_phlda3.py`, `cerna_phlda3.py` |
| Drug repurposing | `l1000cds2_phlda3.py`, `run_gdsc2.R` + `gdsc2_phlda3.py`, `cmap_clean_phlda3.py` |
| Single cell (tumour + p53 module) | `sc_tumor_phlda3.py`, `sc_p53_module_thca.py`, `combined_singlecell_nature.py` |
| Protein (HPA) | `hpa_ihc_thca.py` |

Upstream helpers: `screen_thca.py` (DE screen), `map_de_symbols.py` (Ensembl→symbol), `sc_normal_census.py` (normal-thyroid scRNA fetch). All figures share `code/nature_style.py` (Arial, unified palette, PDF export).

## Data availability / attribution

TCGA (GDC), GEO, UCSC Xena, cBioPortal, Human Protein Atlas, ENCORI, L1000CDS2, GDSC2 (oncoPredict) and CELLxGENE Census are publicly available; cite the original data producers and methods papers listed in the manuscript.
