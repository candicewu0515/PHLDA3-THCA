# Response to Reviewers

We thank the reviewer for a careful and constructive evaluation. The comments have substantially improved the rigour and tone of the manuscript. We have performed two additional analyses, added one expanded limitation, and revised the language throughout to avoid over-statement. Below we respond to each point; reviewer comments are in italics and our changes reference the revised manuscript and Supplementary figures. All analysis code is openly available (see Data Availability).

---

## Major comments

**Major 1 — The N1 model could be more clinically independent (histology, extrathyroidal extension, multifocality, tumour size, TERT promoter, RET/NTRK fusion).**

We agree. We examined the harmonised TCGA-THCA clinical table available to us and it does not contain tumour size, extrathyroidal extension, histological subtype or multifocality; AJCC T stage—already in the model—partly captures tumour size and extrathyroidal extension, which we now state explicitly. We attempted to incorporate the key molecular events but they were not available: TERT promoter mutations (C228T/C250T) are non-coding and absent from the standard TCGA-THCA whole-exome calls (cBioPortal returns only three coding TERT variants, no promoter mutations), and harmonised per-sample RET/NTRK fusion calls were not available through the queried TCGA/cBioPortal resources in a model-ready form. Rather than over-state the model, we have added an explicit limitation listing these unmeasured molecular and clinicopathological variables and recommending validation in cohorts with targeted/fusion sequencing and complete pathology (Limitations). We also strengthened the BRAF analysis as detailed in Major 2.

**Major 2 — BRAF should be stratified, not only adjusted (PHLDA3 vs N1 within BRAF-mutant and BRAF-wildtype; interaction).**

Done. We added BRAF-stratified N1 models and a PHLDA3×BRAF interaction term (new Methods 2.3; Results 3.4; shown as panel c of main Figure 2). The effect direction was consistent in both subgroups—BRAF-mutant OR 1.27 (95% CI 0.91–1.79) and BRAF-wildtype OR 1.27 (95% CI 0.92–1.75)—with no significant PHLDA3×BRAF interaction (P = 0.96). As expected, the subgroup analyses were individually underpowered, so we report them as directionally consistent sensitivity analyses rather than definitive subgroup evidence.

**Major 3 — Single-cell p-values risk pseudoreplication; the independent unit is the patient/sample.**

We fully agree and have addressed this directly. We re-tested the malignant-versus-normal-epithelial contrast after first aggregating expression per sample (new Methods 2.7; Results 3.11; shown as panels e,f of main Figure 5). At the sample level (14 malignant-cell samples vs 5 normal-epithelial samples), PHLDA3 remained higher in malignant cells (mean expression Mann–Whitney P = 0.002; percentage of cells detected P = 0.014). The text now states that the cell-level comparison is descriptive and that the conclusion rests on the effect size and detection proportion rather than the cell-level p-value. We applied the same sample-level safeguard to the within-malignant PHLDA3–p53 module-score correlation, which persisted as a pseudobulk correlation across malignant samples (Spearman r = 0.58, P = 0.030, n = 14); we also removed the cell-level p-value from the Abstract and now report the malignant-versus-normal difference there at the sample level, so no single-cell claim depends on pseudoreplication.

**Major 4 — The lead CpG is data-driven; avoid cherry-picking.**

We have clarified the selection procedure to make it transparent and rule-based (Methods 2.5). All 13 PHLDA3 HM450 probes were screened; promoter probes were defined by genomic position (within 1.5 kb of the Ensembl TSS) before identifying the lead site; cg04055835 was then reported as the promoter CpG with the strongest inverse expression correlation. The full per-probe results (genomic distance to TSS, promoter status, tumour-versus-normal methylation, and methylation–expression correlation) are shown in Supplementary Fig. S2 so the selection is fully visible.

**Major 5 — "via/through promoter hypomethylation" is too causal.**

Corrected throughout. We now write that PHLDA3 over-expression is "associated with" (not "via/through") promoter hypomethylation (Abstract, Background, Results 3.5, Discussion, Conclusions). We acknowledge that no methylation editing, 5-aza treatment or chromatin (ATAC/H3K27ac/H3K4me3) data were generated, so causality is not claimed.

**Major 6 — The ceRNA network is scientifically weak and should be further downgraded.**

We agree and have downgraded it. ceRNA has been removed from the Abstract and keywords; the Results section is shortened and explicitly states the network is "hypothesis-generating rather than mechanistic proof," listing the experiments required (miRNA mimic/inhibitor, lncRNA knockdown with miRNA-inhibitor rescue, luciferase-reporter 3′-UTR assays, RIP/RNA pull-down); the figure is in the Supplementary material (S4). We no longer use "sponging contributes to up-regulation"; for SNHG12 we state only that it was positively correlated with PHLDA3 and shared CLIP-supported miRNA regulators, suggesting a possible ceRNA relationship.

**Major 7 — L1000 drug repurposing reads as an add-on and may distract from the N1-biomarker focus.**

We have de-emphasised it. Drug repurposing has been removed from the Abstract; in the Results it is condensed and framed as exploratory in-silico analysis, with the softer wording "compounds predicted to oppose the PHLDA3-high co-expression signature" rather than "candidate to revert the PHLDA3-high state." The figure remains in the Supplementary material (S5), and we state that these predictions are hypothesis-generating and were not experimentally or clinically validated in thyroid cancer.

**Major 8 — "pre-operative" is too strong given surgical-tumour RNA-seq.**

Corrected. Because our expression data derive from resected tumour (not fine-needle-aspiration or blood), we now write "candidate risk-stratification marker" and "may inform future pre-operative risk-stratification studies," and the Discussion explicitly notes that pre-operative application would require validation in FNA/pre-operative samples.

**Major 9 — The external cohort validates expression, not the primary endpoint; avoid implying N1 was externally validated.**

Agreed. The Conclusion has been reworded to "Multi-cohort expression validation together with TCGA-based clinical modelling nominate PHLDA3 as a candidate marker linked to lymph-node metastasis," and external validation of the clinical (N1) association in annotated cohorts is named as a key next step. The Abstract and Limitations are consistent with this.

---

## Minor comments

- **Title.** Changed to "…identifies PHLDA3 as a BRAF-associated transcriptomic marker linked to lymph-node metastasis in thyroid carcinoma."
- **"independently predicts."** Retained in the Results where the regression supports it; in the Conclusions we now use "independently associated with."
- **Figure 1 ROC.** The legend and Methods 2.2 continue to state that the tumour-versus-normal ROC is confirmation of over-expression, not a clinical diagnostic test.
- **Nomogram AUC 0.68.** Now described as "modest discrimination"; we no longer use "acceptable" repeatedly.
- **"immunosuppressive microenvironment."** Changed to "checkpoint-high immune microenvironment" (Section 3.8 heading and text, Abstract, Conclusions), with a sentence noting that ssGSEA plus checkpoint-gene correlation cannot fully define immunosuppression.
- **HPA n = 4.** Presented only as qualitative protein-level support in the Supplementary material (S6); it does not carry the main argument.

---

We believe these revisions make the central claim appropriately precise and well supported, while clearly framing the regulatory, drug and ceRNA analyses as exploratory. We thank the reviewer for comments that materially improved the manuscript.
