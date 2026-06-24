#!/usr/bin/env python3
"""PHLDA3 co-expression in THCA tumors -> GO/KEGG (Enrichr on top co-expressed genes)
+ GSEA preranked on the full PHLDA3-correlation ranking."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gseapy as gp
import nature_style as ns
ns.set_style(font_size=8)

GENE, ENS = "PHLDA3", "ENSG00000174307"

# ---------- tumor logCPM matrix, expression filter ----------
df = pd.read_csv("TCGA-THCA.htseq_counts.tsv", sep="\t", index_col=0)
mat = (2**df - 1).clip(lower=0)
stype = pd.Series([b.split("-")[3][:2] for b in df.columns], index=df.columns)
tum = df.columns[stype == "01"]
cpm = mat[tum].div(mat[tum].sum(axis=0), axis=1) * 1e6
keep = cpm.mean(axis=1) >= 1
logcpm = np.log2(cpm[keep] + 1)
logcpm.index = [g.split(".")[0] for g in logcpm.index]

# map symbols
m = pd.read_csv("thca_DE_full_symbol.csv")[["ensembl","symbol"]].dropna()
e2s = dict(zip(m.ensembl, m.symbol))
logcpm = logcpm[logcpm.index.isin(e2s)]
logcpm.index = [e2s[e] for e in logcpm.index]
logcpm = logcpm[~logcpm.index.duplicated()]
print("genes x tumors:", logcpm.shape)

# ---------- correlate every gene with PHLDA3 ----------
phl = logcpm.loc[GENE]
X = logcpm.sub(logcpm.mean(axis=1), axis=0)
r = X.values @ (phl - phl.mean()).values / (
        np.sqrt((X.values**2).sum(axis=1)) * np.sqrt(((phl-phl.mean())**2).sum()))
corr = pd.Series(r, index=logcpm.index).drop(GENE).sort_values(ascending=False)
corr.to_csv("PHLDA3_coexpression.csv", header=["pearson_r"])
print(f"top + : {list(corr.head(8).index)}")
print(f"top - : {list(corr.tail(8).index)}")

# ---------- Enrichr: GO BP + KEGG on top positively co-expressed genes ----------
top_pos = corr[corr > 0.4].head(300).index.tolist()
print(f"co-expressed genes (r>0.4, top300): {len(top_pos)}")
enr = gp.enrichr(gene_list=top_pos, organism="human",
                 gene_sets=["GO_Biological_Process_2021", "KEGG_2021_Human"],
                 outdir=None)
e = enr.results.copy()
e["nlogq"] = -np.log10(e["Adjusted P-value"].clip(lower=1e-300))
e["Term"] = e["Term"].str.replace(r"\(GO:\d+\)", "", regex=True).str.strip()
e.to_csv("PHLDA3_enrichr_GO_KEGG.csv", index=False)

# ---------- GSEA preranked on full correlation ranking ----------
rnk = corr.reset_index(); rnk.columns = ["gene", "r"]
pre = gp.prerank(rnk=rnk, gene_sets="MSigDB_Hallmark_2020",
                 min_size=10, max_size=500, permutation_num=1000,
                 seed=42, outdir=None, no_plot=True)
g = pre.res2d.copy()
for col in ["NES","FDR q-val","NOM p-val"]:
    g[col] = pd.to_numeric(g[col], errors="coerce")
g = g.sort_values("NES", ascending=False)
g.to_csv("PHLDA3_GSEA_hallmark.csv", index=False)
print(f"GSEA hallmark sets: {len(g)}, sig(FDR<0.25): {(g['FDR q-val']<0.25).sum()}")

# ================= figure: GO | KEGG | GSEA =================
fig, ax = plt.subplots(1, 3, figsize=(16, 5))

def barh(a, sub, title, color):
    sub = sub.sort_values("nlogq").tail(10)
    a.barh(range(len(sub)), sub["nlogq"], color=color, alpha=0.85, edgecolor="black", lw=0.4)
    a.set_yticks(range(len(sub)))
    a.set_yticklabels([t[:42] for t in sub["Term"]], fontsize=8)
    a.set_xlabel("-log10 adj.P"); a.set_title(title, fontweight="bold", fontsize=11)
    a.spines[["top","right"]].set_visible(False)

barh(ax[0], e[e["Gene_set"]=="GO_Biological_Process_2021"], "GO Biological Process", ns.C["accent"])
barh(ax[1], e[e["Gene_set"]=="KEGG_2021_Human"], "KEGG Pathway", ns.C["normal"])

# GSEA: top |NES| hallmark bars (signed)
gg = g.reindex(g["NES"].abs().sort_values().index).tail(12).sort_values("NES")
cols = [ns.C["tumor"] if v>0 else ns.C["normal"] for v in gg["NES"]]
ax[2].barh(range(len(gg)), gg["NES"], color=cols, alpha=0.85, edgecolor="black", lw=0.4)
ax[2].set_yticks(range(len(gg)))
ax[2].set_yticklabels([t.replace("_"," ").title()[:34] for t in gg["Term"]], fontsize=8)
ax[2].axvline(0, color="black", lw=0.8)
ax[2].set_xlabel("NES (PHLDA3-correlated)")
ax[2].set_title("GSEA Hallmark", fontweight="bold", fontsize=11)
ax[2].spines[["top","right"]].set_visible(False)

for i, lab in enumerate(["a","b","c"]):
    ns.panel_label(ax[i], lab)
fig.suptitle(f"{GENE} co-expression functional enrichment (TCGA-THCA tumors, n={logcpm.shape[1]})",
             fontweight="bold", fontsize=12, y=1.02)
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_enrichment")
print("+ 3 CSVs")
