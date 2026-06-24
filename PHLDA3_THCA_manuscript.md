# PHLDA3 is induced in thyroid carcinoma cells and is an independent predictor of lymph-node metastasis and advanced stage: an integrated multi-cohort and single-cell analysis

*(Working draft for submission to the Egyptian Journal of Medical Human Genetics / Journal of Genetic Engineering and Biotechnology / Bioinformation type venue. Author list, affiliations, and corresponding author to be completed by the submitting team.)*

---

## Abstract

**Background.** Pleckstrin-homology-like domain family A member 3 (PHLDA3) is a p53 target gene that represses AKT and behaves as a tumour suppressor in neuroendocrine tumours. Its role in thyroid carcinoma (THCA) is unknown; to our knowledge no dedicated study has examined PHLDA3 in thyroid cancer.

**Methods.** We integrated bulk RNA-seq and clinical data from The Cancer Genome Atlas (TCGA-THCA, 502 tumours), pan-cancer TCGA+GTEx expression, two independent Gene Expression Omnibus cohorts (GSE33630, GSE193581), somatic mutation calls from cBioPortal, the Human Protein Atlas, and two single-cell RNA-seq resources (CELLxGENE Census; GSE193581). We assessed differential expression, diagnostic performance (ROC), clinicopathological correlations, multivariable logistic regression, a nomogram with calibration and decision-curve analysis, driver-mutation association, co-expression-based functional enrichment (GO/KEGG/GSEA), immune-cell infiltration (ssGSEA) and immune-checkpoint correlation, and cell-type-resolved expression.

**Results.** PHLDA3 was significantly over-expressed in THCA versus normal thyroid (log2 fold-change +1.77, P = 4.9 × 10⁻³¹) with high diagnostic accuracy (TCGA AUC 0.965; external GSE33630 AUC 0.932, P = 4.5 × 10⁻¹⁴), ranking among the most up-regulated cancers pan-cancer. High PHLDA3 was associated with higher T stage, lymph-node metastasis (N1), distant metastasis and advanced AJCC stage. In multivariable logistic regression adjusting for age, sex and T stage, PHLDA3 independently predicted N1 metastasis (OR 1.42 per SD, 95% CI 1.15–1.75, P = 1.0 × 10⁻³) and advanced stage (OR 1.63, 95% CI 1.21–2.18, P = 1.2 × 10⁻³); a PHLDA3-based nomogram showed acceptable discrimination (10-fold cross-validated AUC 0.68) and net clinical benefit. PHLDA3 was highest in BRAF^V600E^-mutant tumours (P < 1 × 10⁻⁴). Co-expressed genes converged on the p53 signalling pathway, apoptosis, interferon responses and epithelial–mesenchymal transition. PHLDA3 correlated with an immunosuppressive microenvironment (regulatory T cells, dendritic cells) and with immune-checkpoint genes (CD276, LGALS9, HAVCR2, CD274, CTLA4; all P < 0.001). At single-cell resolution PHLDA3 was nearly absent from normal thyroid follicular cells (0.4% of cells) but became the predominant PHLDA3-expressing population in malignant cells of thyroid tumours (47.5%), indicating tumour-cell-intrinsic induction.

**Conclusions.** PHLDA3 is transcriptionally induced in thyroid cancer cells and marks nodal spread, advanced stage and the BRAF^V600E^ aggressive subtype, supporting its evaluation as a diagnostic and risk-stratification biomarker in THCA.

**Keywords:** PHLDA3; thyroid carcinoma; TCGA; biomarker; lymph-node metastasis; tumour microenvironment; single-cell RNA-seq

---

## 1. Background

Thyroid cancer is the most common endocrine malignancy and one of the fastest-rising cancers in incidence worldwide, a pattern largely attributed to over-diagnosis of indolent papillary tumours [4]. Although most patients have an excellent prognosis, a clinically important subset presents with lymph-node metastasis, extrathyroidal extension or progression to poorly differentiated and anaplastic disease. Papillary thyroid carcinoma (PTC) is driven predominantly by mutually exclusive mutations in *BRAF* (most often BRAF^V600E^) and *RAS* genes, and BRAF^V600E^ is associated with more aggressive clinicopathological features including lymph-node metastasis and advanced stage [3,5]. Robust molecular markers that identify aggressive behaviour at diagnosis remain a clinical need.

Pleckstrin-homology-like domain family A member 3 (PHLDA3) is a small phosphatidylinositol-binding protein and a direct transcriptional target of p53. By competing with AKT for membrane phosphoinositides, PHLDA3 represses AKT activation and, through this mechanism, functions as a tumour suppressor: in pancreatic neuroendocrine tumours the *PHLDA3* locus undergoes loss of heterozygosity and methylation, and its inactivation cooperates with *MEN1* loss to drive tumourigenesis [1]. Other PHLDA-family members similarly repress AKT and show tumour-suppressive activity [2]. PHLDA3 is therefore generally regarded as a context-dependent tumour suppressor, yet its expression and clinical relevance in thyroid carcinoma have not, to our knowledge, been reported—a search of PubMed returns no dedicated PHLDA3–thyroid study.

Here we systematically characterise PHLDA3 in THCA by integrating bulk transcriptomes, an external validation cohort, pan-cancer reference data, somatic mutation profiles, protein-level immunohistochemistry, functional enrichment, immune deconvolution and single-cell RNA-seq. We show that, contrary to its canonical tumour-suppressor behaviour in neuroendocrine tissue, PHLDA3 is transcriptionally induced in thyroid cancer cells and independently marks lymph-node metastasis, advanced stage and the BRAF^V600E^ subtype.

## 2. Materials and methods

### 2.1 Data sources
TCGA-THCA HTSeq counts and clinical annotations (502 primary tumours with paired expression and clinicopathological data; adjacent normal thyroid samples) were obtained from the Genomic Data Commons. Pan-cancer tumour-versus-normal comparison used the harmonised TCGA + GTEx expression matrix (RSEM, log2 norm-count) from the UCSC Xena Toil hub, with primary site used to pair each TCGA tumour to TCGA and GTEx normals. Survival/clinical endpoint definitions followed the TCGA Pan-Cancer Clinical Data Resource [11]. Two independent thyroid cohorts were retrieved from the Gene Expression Omnibus: GSE33630 (Affymetrix HG-U133 Plus 2.0; 45 normal, 49 PTC, 11 anaplastic thyroid carcinoma [ATC]) for external diagnostic validation, and GSE193581 (author-annotated single-cell RNA-seq of PTC, ATC and normal thyroid) for tumour-cell analysis. Somatic mutation calls (BRAF, HRAS, KRAS, NRAS) for the TCGA-THCA PanCancer Atlas cohort were queried from cBioPortal. Protein expression was examined in the Human Protein Atlas [12]. Normal-thyroid single-cell expression was queried from the CZI CELLxGENE Census.

### 2.2 Expression, diagnostic and clinicopathological analyses
HTSeq counts were converted to counts-per-million and log2-transformed; genes were collapsed to one value per patient. Tumour-versus-normal differences were tested by the Mann–Whitney *U* test and Benjamini–Hochberg false-discovery-rate (FDR) correction; diagnostic discrimination was summarised by the area under the receiver-operating-characteristic curve (AUC). PHLDA3-high and -low groups were defined by the median. Associations with age, sex, T stage, N stage, M stage and AJCC stage were tested by Fisher's exact / χ² tests and Mann–Whitney *U* tests.

### 2.3 Logistic regression, nomogram and decision-curve analysis
Multivariable logistic regression modelled lymph-node metastasis (N1 vs N0) and advanced stage (III/IV vs I/II) as a function of standardised PHLDA3 expression, age, sex and T stage. Odds ratios are reported per one standard deviation of PHLDA3. A nomogram was constructed from the N1 model (PHLDA3, age, T stage); discrimination was assessed by 10-fold cross-validated AUC, calibration by quantile calibration curves, and clinical utility by decision-curve analysis (net benefit across threshold probabilities) [8].

### 2.4 Mutation, enrichment, immune and single-cell analyses
PHLDA3 expression was compared across BRAF^V600E^, RAS-mutant and wild-type/other tumours (Mann–Whitney *U*; Kruskal–Wallis). Genes co-expressed with PHLDA3 across tumours (Pearson correlation) were submitted to GO Biological Process and KEGG over-representation analysis (Enrichr [15]); the full ranked correlation list was analysed by Gene Set Enrichment Analysis against the MSigDB Hallmark collection [7]. Immune-cell infiltration was estimated by single-sample GSEA [10] using the Danaher immune-cell signatures [14] and correlated with PHLDA3; immune-checkpoint genes were correlated directly. CIBERSORT-type deconvolution frameworks [6] and GEPIA [13] informed methodological choices. For single-cell analysis, normal-thyroid cells (CELLxGENE Census) and tumour cells (GSE193581, author cell-type labels, CP10k-normalised) were summarised as the percentage of cells expressing PHLDA3 and mean expression per cell type. All analyses were performed in Python 3.13.11 using pandas 2.3.3, numpy 2.4.4, scipy 1.17.1, matplotlib 3.10.9, scikit-learn 1.8.0, statsmodels 0.14.6, lifelines 0.30.3 and gseapy 1.3.0; pan-cancer and single-cell data were accessed with xenaPython 1.0.14 and cellxgene-census 1.18.0, and gene identifiers were mapped with mygene. Enrichment used the GO_Biological_Process_2021, KEGG_2021_Human and MSigDB_Hallmark_2020 libraries. P < 0.05 was considered significant; multiple testing was controlled by the Benjamini–Hochberg FDR where indicated.

## 3. Results

### 3.1 PHLDA3 is over-expressed in thyroid carcinoma and across cancers (Fig. 1)
In TCGA-THCA, PHLDA3 was significantly up-regulated in tumours relative to normal thyroid (log2 fold-change +1.77; FDR = 3.2 × 10⁻²⁸) and discriminated tumour from normal with high accuracy (AUC 0.965). In the harmonised TCGA + GTEx pan-cancer comparison, THCA showed the fourth-largest tumour-versus-normal up-regulation of PHLDA3 among 24 cancer types (P = 4.9 × 10⁻³¹), confirming that PHLDA3 induction is a prominent feature of thyroid cancer rather than a pan-cancer artefact.

### 3.2 External and protein-level validation (Fig. 2)
In the independent GSE33630 cohort, PHLDA3 was higher in PTC and ATC than in normal thyroid (Kruskal–Wallis P = 1.7 × 10⁻¹³) and separated tumour from normal with an AUC of 0.932 (Mann–Whitney P = 4.5 × 10⁻¹⁴), reproducing the TCGA finding in a different platform and patient population. In the Human Protein Atlas, PHLDA3 protein was detectable in thyroid tumour cells with cytoplasmic and nuclear localisation; immunohistochemical staining intensity was limited by small antibody-cohort size and is reported qualitatively.

### 3.3 PHLDA3 correlates with aggressive clinicopathology (Fig. 3)
PHLDA3-high tumours were enriched for higher T category (P = 7.5 × 10⁻³), lymph-node metastasis (N1; P = 6.5 × 10⁻³), distant metastasis (P = 0.043) and advanced AJCC stage (P = 4.4 × 10⁻³). Continuous PHLDA3 expression was higher in T3/T4 than T1/T2, in N1 than N0, and in stage III/IV than I/II tumours.

### 3.4 PHLDA3 independently predicts nodal metastasis and advanced stage (Fig. 4)
In multivariable logistic regression adjusting for age, sex and T stage, PHLDA3 remained an independent predictor of lymph-node metastasis (OR 1.42 per SD, 95% CI 1.15–1.75, P = 1.0 × 10⁻³; n = 450) and of advanced stage (OR 1.63, 95% CI 1.21–2.18, P = 1.2 × 10⁻³; n = 496). A PHLDA3-based nomogram for N1 metastasis (PHLDA3 + age + T stage) achieved a 10-fold cross-validated AUC of 0.68 with acceptable calibration, and decision-curve analysis demonstrated positive net benefit over "treat-all" and "treat-none" strategies across a clinically relevant range of threshold probabilities.

### 3.5 PHLDA3 is highest in BRAF^V600E^ tumours (Fig. 5)
Among 489 sequenced TCGA-THCA tumours (BRAF mutant n = 286, of which 283 BRAF^V600E^; RAS mutant n = 59), PHLDA3 expression was significantly higher in BRAF-mutant than wild-type tumours (P < 1 × 10⁻⁴) and lower in RAS-mutant tumours (P = 3.7 × 10⁻⁴). Across driver subtypes, PHLDA3 was highest in BRAF^V600E^ (median 5.82), intermediate in wild-type/other (5.36) and lowest in RAS-mutant tumours (5.24; Kruskal–Wallis P < 1 × 10⁻⁴), aligning PHLDA3 with the more aggressive BRAF^V600E^-driven phenotype.

### 3.6 PHLDA3 co-expression implicates p53, apoptosis, interferon and EMT programmes (Fig. 6)
Genes positively co-expressed with PHLDA3 were strongly enriched for the KEGG p53 signalling pathway and apoptosis, consistent with PHLDA3's identity as a p53 target; the top positively correlated genes included canonical p53 effectors (BBC3/PUMA, BAX, DDB2). GSEA against Hallmark gene sets showed positive enrichment for interferon-γ and interferon-α responses, the p53 pathway, inflammatory response, epithelial–mesenchymal transition and apoptosis.

### 3.7 PHLDA3 marks an immunosuppressive, checkpoint-high microenvironment (Fig. 7)
ssGSEA-estimated infiltration showed positive correlation of PHLDA3 with regulatory T cells (r = 0.34), dendritic cells, mast cells and exhausted CD8 T cells (all P < 0.001). PHLDA3 correlated positively with multiple immune-checkpoint genes, most strongly CD276/B7-H3 (r = 0.53), LGALS9/galectin-9, HAVCR2/TIM-3, CD274/PD-L1 and CTLA4 (all P < 0.001), indicating association with an immunosuppressive, checkpoint-enriched tumour microenvironment concordant with the interferon signature.

### 3.8 PHLDA3 is induced in malignant thyroid cells at single-cell resolution (Fig. 8)
In normal thyroid (CELLxGENE Census, 15,646 cells), PHLDA3 was expressed in epithelial and basal compartments but was nearly absent from thyroid follicular cells—the cell of origin of THCA—(0.4% of cells positive). By contrast, in thyroid tumours (GSE193581, 67,414 cells), malignant cells were the predominant PHLDA3-expressing population (47.5% positive, the highest of any cell type), and PHLDA3 was markedly higher in PTC and ATC malignant cells than in normal epithelial cells (P < 1 × 10⁻⁴). Together these data indicate that PHLDA3 over-expression in bulk tumour is tumour-cell-intrinsic and induced upon malignant transformation rather than reflecting a shift in cellular composition.

## 4. Discussion

Across seven complementary data modalities we find a consistent picture: PHLDA3 is transcriptionally induced in thyroid cancer cells, discriminates tumour from normal tissue with high accuracy in two independent cohorts, and independently marks lymph-node metastasis, advanced stage and the BRAF^V600E^ aggressive subtype. The single-cell data resolve a key ambiguity of bulk analyses by localising the induction to malignant cells: PHLDA3 is essentially silent in normal follicular cells but is the leading PHLDA3-expressing population in tumours, implying cell-intrinsic up-regulation during transformation.

This behaviour is, at first sight, paradoxical. PHLDA3 is a p53 target and AKT repressor that acts as a tumour suppressor in pancreatic neuroendocrine tumours, where its locus is lost and methylated [1,2]. Two non-exclusive explanations reconcile our findings with that biology. First, PHLDA3 induction may be a consequence of p53-pathway activation in thyroid tumours: its top co-expressed genes are canonical p53 effectors (PUMA, BAX, DDB2) and its co-expression signature is enriched for the p53 pathway and apoptosis, so high PHLDA3 may mark tumours with an engaged but ultimately insufficient p53 stress response rather than acting as a driver. Second, the tumour-suppressive function of PHLDA3 is tissue-context-dependent; a p53-inducible gene can be a passenger marker of aggressiveness in an epithelial cancer even where it restrains AKT in neuroendocrine tissue. The association of PHLDA3 with BRAF^V600E^—which activates MAPK signalling and is itself linked to nodal metastasis and advanced stage [3,5]—supports its interpretation as a marker of the aggressive, MAPK-driven thyroid phenotype.

The microenvironmental correlations add a further layer. PHLDA3-high tumours show an interferon-response signature and positive correlation with regulatory T cells and with several immune-checkpoint molecules, most strikingly CD276/B7-H3 and the TIM-3/galectin-9 axis. This immunosuppressive, checkpoint-enriched profile is consistent with the more aggressive thyroid phenotype—including the stromal and microenvironmental remodelling described in single-cell studies of anaplastic thyroid carcinoma [9]—and raises the hypothesis, requiring experimental testing, that PHLDA3-high tumours might be candidates for checkpoint-directed strategies.

Clinically, the value of PHLDA3 here is as a diagnostic and risk-stratification marker rather than a prognostic one. Its diagnostic AUC (0.93–0.97 across cohorts) is high, and it independently predicts nodal metastasis and advanced stage; the nomogram and decision-curve analysis frame this as potential clinical utility for pre-operative risk assessment, where prediction of occult nodal disease informs the extent of surgery. We deliberately do not claim a survival association: thyroid cancer has very few events in TCGA, and PHLDA3 was not associated with progression-free interval, consistent with the Human Protein Atlas classifying it as "unprognostic" for survival in THCA. Reporting this null result transparently avoids over-interpretation.

## 5. Limitations

This study is computational and association-based; it does not establish causality, and the inferred relationships (co-expression, immune correlations, single-cell detection rates) require experimental validation. Immune deconvolution and ssGSEA are estimates rather than direct measurements, and the correlation magnitudes, although highly significant at this sample size, are modest. The Human Protein Atlas thyroid-tumour cohort is small (n = 4), so protein-level evidence is qualitative. The normal-thyroid single-cell data contain no tumour cells and the follicular-cell detection rate may be partly affected by single-cell drop-out; we therefore report "percentage of cells in which PHLDA3 is detected" and corroborate the contrast with an independent tumour single-cell dataset. Functional and mechanistic experiments (knock-down/over-expression in thyroid cell lines, AKT-pathway read-outs, and analysis of patient outcomes in event-rich cohorts) are needed to determine whether PHLDA3 is a driver or a passenger marker of aggressiveness.

## 6. Conclusions

PHLDA3 is induced in thyroid carcinoma cells and is an accurate diagnostic marker and an independent predictor of lymph-node metastasis and advanced stage that tracks with the BRAF^V600E^ subtype and an immunosuppressive microenvironment. These integrated, multi-cohort and single-cell findings nominate PHLDA3 as a candidate biomarker for diagnosis and pre-operative risk stratification in thyroid cancer and motivate functional studies of its role downstream of p53 in the thyroid context.

## Declarations

**Ethics approval and consent to participate.** Not applicable; this study used de-identified, publicly available data from TCGA, GEO, the Human Protein Atlas, cBioPortal and CELLxGENE.
**Consent for publication.** Not applicable.
**Availability of data and materials.** All data are publicly available: TCGA-THCA (Genomic Data Commons, https://portal.gdc.cancer.gov), TCGA+GTEx (UCSC Xena, https://xenabrowser.net), GSE33630 and GSE193581 (GEO, https://www.ncbi.nlm.nih.gov/geo), cBioPortal (https://www.cbioportal.org), Human Protein Atlas (https://www.proteinatlas.org), CELLxGENE Census (https://cellxgene.cziscience.com). All analysis code, software versions and figure-generating scripts are openly available at https://github.com/candicewu0515/PHLDA3-THCA.
**Competing interests.** The authors declare no competing interests.
**Funding.** To be completed.
**Authors' contributions.** To be completed.

## References

1. Ohki R, Saito K, Chen Y, et al. PHLDA3 is a novel tumor suppressor of pancreatic neuroendocrine tumors. *Proc Natl Acad Sci USA*. 2014;111(23):E2404–E2413. doi:10.1073/pnas.1319962111
2. Chen Y, Takikawa M, Tsutsumi S, et al. PHLDA1, another PHLDA family protein that inhibits Akt. *Cancer Sci*. 2018;109(11):3532–3542. doi:10.1111/cas.13796
3. Cancer Genome Atlas Research Network. Integrated genomic characterization of papillary thyroid carcinoma. *Cell*. 2014;159(3):676–690. doi:10.1016/j.cell.2014.09.050
4. Pizzato M, Li M, Vignat J, et al. The epidemiological landscape of thyroid cancer worldwide: GLOBOCAN estimates for incidence and mortality rates in 2020. *Lancet Diabetes Endocrinol*. 2022;10(4):264–272. doi:10.1016/S2213-8587(22)00035-3
5. Wei X, Wang X, Xiong J, et al. Risk and prognostic factors for BRAF mutations in papillary thyroid carcinoma. *Biomed Res Int*. 2022;2022:9959649. doi:10.1155/2022/9959649
6. Newman AM, Liu CL, Green MR, et al. Robust enumeration of cell subsets from tissue expression profiles. *Nat Methods*. 2015;12(5):453–457. doi:10.1038/nmeth.3337
7. Subramanian A, Tamayo P, Mootha VK, et al. Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles. *Proc Natl Acad Sci USA*. 2005;102(43):15545–15550. doi:10.1073/pnas.0506580102
8. Vickers AJ, Elkin EB. Decision curve analysis: a novel method for evaluating prediction models. *Med Decis Making*. 2006;26(6):565–574. doi:10.1177/0272989X06295361
9. Pan Z, Xu T, Bao L, et al. CREB3L1 promotes tumor growth and metastasis of anaplastic thyroid carcinoma by remodeling the tumor microenvironment. *Mol Cancer*. 2022;21(1):190. doi:10.1186/s12943-022-01658-x
10. Barbie DA, Tamayo P, Boehm JS, et al. Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. *Nature*. 2009;462(7269):108–112. doi:10.1038/nature08460
11. Liu J, Lichtenberg T, Hoadley KA, et al. An integrated TCGA pan-cancer clinical data resource to drive high-quality survival outcome analytics. *Cell*. 2018;173(2):400–416.e11. doi:10.1016/j.cell.2018.02.052
12. Uhlén M, Fagerberg L, Hallström BM, et al. Tissue-based map of the human proteome. *Science*. 2015;347(6220):1260419. doi:10.1126/science.1260419
13. Tang Z, Li C, Kang B, et al. GEPIA: a web server for cancer and normal gene expression profiling and interactive analyses. *Nucleic Acids Res*. 2017;45(W1):W98–W102. doi:10.1093/nar/gkx247
14. Danaher P, Warren S, Dennis L, et al. Gene expression markers of tumor infiltrating leukocytes. *J Immunother Cancer*. 2017;5:18. doi:10.1186/s40425-017-0215-8
15. Kuleshov MV, Jones MR, Rouillard AD, et al. Enrichr: a comprehensive gene set enrichment analysis web server 2016 update. *Nucleic Acids Res*. 2016;44(W1):W90–W97. doi:10.1093/nar/gkw377

## Figure legends

**Figure 1. PHLDA3 is over-expressed in thyroid carcinoma and across cancers.** (a) Whole-transcriptome volcano plot of TCGA-THCA tumour-versus-normal differential expression; PHLDA3 is highlighted. (b) PHLDA3 expression in tumour versus normal thyroid and diagnostic ROC (AUC 0.965). (c) Pan-cancer tumour-versus-normal PHLDA3 expression (TCGA tumour vs TCGA + GTEx normal); THCA highlighted (***P < 0.001).

**Figure 2. External and protein-level validation.** (a) PHLDA3 across normal thyroid, PTC and ATC in GSE33630 and (b) diagnostic ROC (AUC 0.932). (c) Representative Human Protein Atlas immunohistochemistry of PHLDA3 in normal thyroid and thyroid cancer.

**Figure 3. Clinicopathological associations.** PHLDA3 expression by T category, lymph-node status (N) and AJCC stage; P values by Mann–Whitney U test.

**Figure 4. PHLDA3 independently predicts nodal metastasis and advanced stage.** (a) Multivariable logistic-regression forest plots for N1 metastasis and advanced stage. (b) Nomogram for predicting N1 metastasis with ROC, calibration curve and decision-curve analysis.

**Figure 5. PHLDA3 and driver mutations.** PHLDA3 expression by BRAF status, RAS status and driver subtype (BRAF^V600E^ / RAS / wild-type-other).

**Figure 6. Co-expression functional enrichment.** GO Biological Process and KEGG over-representation of PHLDA3-co-expressed genes, and GSEA of Hallmark gene sets ranked by correlation with PHLDA3.

**Figure 7. PHLDA3 and the immune microenvironment.** Correlation of PHLDA3 with ssGSEA immune-cell infiltration (Danaher signatures) and with immune-checkpoint genes.

**Figure 8. PHLDA3 is induced in malignant thyroid cells.** Percentage of cells expressing PHLDA3 across cell types in (a) normal thyroid (CELLxGENE Census) and (b) thyroid tumours (GSE193581); the cell of origin (follicular cells) and malignant cells are highlighted.
