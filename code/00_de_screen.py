#!/usr/bin/env python3
"""THCA tumor-vs-normal DE + diagnostic-AUC screen to find under-published candidate biomarker genes."""
import pandas as pd, numpy as np
from scipy.stats import mannwhitneyu

F = "data/TCGA-THCA.htseq_counts.tsv"
print("loading", F, "...", flush=True)
df = pd.read_csv(F, sep="\t", index_col=0)
print("matrix:", df.shape, flush=True)

# Xena htseq_counts are log2(count+1). Undo to raw-ish counts for CPM.
mat = (2**df - 1).clip(lower=0)

barcodes = df.columns
stype = pd.Series([b.split("-")[3][:2] for b in barcodes], index=barcodes)
tumor = barcodes[stype == "01"]
normal = barcodes[stype == "11"]
print(f"tumor(01)={len(tumor)}  normal(11)={len(normal)}", flush=True)

# CPM normalize then log2
cpm = mat.div(mat.sum(axis=0), axis=1) * 1e6
logcpm = np.log2(cpm + 1)

# expression filter: keep genes with mean CPM >= 1 in either group
keep = (cpm[tumor].mean(axis=1) >= 1) | (cpm[normal].mean(axis=1) >= 1)
logcpm = logcpm[keep]
print("genes after expression filter:", logcpm.shape[0], flush=True)

T = logcpm[tumor].values
N = logcpm[normal].values
log2fc = T.mean(axis=1) - N.mean(axis=1)

# Mann-Whitney U gives both p and AUC (U/(n1*n2))
n1, n2 = T.shape[1], N.shape[1]
aucs, pvals = np.empty(len(log2fc)), np.empty(len(log2fc))
for i in range(len(log2fc)):
    u, p = mannwhitneyu(T[i], N[i], alternative="two-sided")
    aucs[i] = u / (n1 * n2)   # AUC for tumor>normal direction
    pvals[i] = p

res = pd.DataFrame({
    "ensembl": [g.split(".")[0] for g in logcpm.index],
    "log2FC_tumor_vs_normal": log2fc,
    "AUC": aucs,                       # >0.5 up in tumor, <0.5 down in tumor
    "AUC_disc": np.abs(aucs - 0.5) + 0.5,  # discrimination regardless of direction
    "pval": pvals,
    "mean_logCPM_tumor": T.mean(axis=1),
})
# BH FDR
from scipy.stats import false_discovery_control
res["FDR"] = false_discovery_control(res["pval"].values)
res = res.sort_values(["AUC_disc", "log2FC_tumor_vs_normal"], ascending=[False, False])

# map top 300 to symbols
import mygene
top = res.head(300).copy()
mg = mygene.MyGeneInfo()
q = mg.querymany(top["ensembl"].tolist(), scopes="ensembl.gene",
                 fields="symbol,type_of_gene", species="human", verbose=False)
m = {d["query"]: (d.get("symbol"), d.get("type_of_gene")) for d in q if not d.get("notfound")}
top["symbol"] = top["ensembl"].map(lambda e: m.get(e, (None, None))[0])
top["gene_type"] = top["ensembl"].map(lambda e: m.get(e, (None, None))[1])

top.to_csv("data/thca_DE_top300.csv", index=False)
res.to_csv("data/thca_DE_full.csv", index=False)
print("\n=== TOP 30 by diagnostic discrimination (|AUC-0.5|), protein-coding ===", flush=True)
pc = top[top["gene_type"] == "protein-coding"].copy()
show = pc.head(30)[["symbol", "ensembl", "log2FC_tumor_vs_normal", "AUC", "FDR", "mean_logCPM_tumor"]]
pd.set_option("display.width", 140)
print(show.to_string(index=False), flush=True)
print("\nsaved: thca_DE_top300.csv, thca_DE_full.csv", flush=True)
