#!/usr/bin/env python3
"""Composite main Figure 1 (over-expression + external validation) — Nature-style:
(a) TCGA-THCA tumour-vs-normal volcano (PHLDA3 highlighted), (b) pan-cancer PHLDA3
tumour-vs-normal delta (THCA highlighted), (c) external GSE33630 validation,
(d) tumour-vs-normal ROC (confirmation of over-expression, not a diagnostic test)."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, roc_auc_score
from scipy.stats import kruskal
import nature_style as ns
ns.set_style(font_size=7)

GENE, ENS = "PHLDA3", "ENSG00000174307"
fig = plt.figure(figsize=(7.2, 6.0)); gs = fig.add_gridspec(2, 2, hspace=0.5, wspace=0.42)
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1]); axC=fig.add_subplot(gs[1,0]); axD=fig.add_subplot(gs[1,1])

# (a) volcano
de = pd.read_csv("data/thca_DE_full_symbol.csv")
de["neglogfdr"] = -np.log10(de.FDR.clip(lower=1e-300))
sig = (de.FDR < 0.05) & (de.log2FC_tumor_vs_normal.abs() > 1)
axA.scatter(de.log2FC_tumor_vs_normal[~sig], de.neglogfdr[~sig], s=2, color=ns.C["ns"], alpha=.3, ec="none")
axA.scatter(de.log2FC_tumor_vs_normal[sig], de.neglogfdr[sig], s=2, color=ns.C["normal"], alpha=.4, ec="none")
ph = de[de.ensembl==ENS]
axA.scatter(ph.log2FC_tumor_vs_normal, ph.neglogfdr, s=40, color=ns.C["tumor"], ec="black", lw=.6, zorder=5)
axA.annotate("PHLDA3", (float(ph.log2FC_tumor_vs_normal), float(ph.neglogfdr)),
             textcoords="offset points", xytext=(6,2), fontsize=7, fontweight="bold", color=ns.C["accent2"])
axA.axvline(0, color=ns.C["grey_mid"], lw=.6, ls="--")
axA.set_xlabel("log2 fold-change (tumour vs normal)", fontsize=6.5); axA.set_ylabel("−log10 FDR", fontsize=6.5)
axA.set_title("TCGA-THCA differential expression", fontweight="bold", fontsize=7.5); axA.spines[["top","right"]].set_visible(False)

# (b) pan-cancer tumour-vs-normal delta, THCA highlighted (disease-paired stats cache)
pc = pd.read_csv("data/PHLDA3_pancancer_stats.csv").sort_values("fc").reset_index(drop=True)
cols = [ns.C["tumor"] if a=="THCA" else ns.C["ns"] for a in pc.abbr]
axB.barh(range(len(pc)), pc.fc.values, color=cols, ec="black", lw=.3)
axB.set_yticks(range(len(pc))); axB.set_yticklabels(pc.abbr, fontsize=5.2)
for lbl in axB.get_yticklabels():
    if lbl.get_text()=="THCA": lbl.set_color(ns.C["tumor"]); lbl.set_fontweight("bold")
rank = pc.sort_values("fc",ascending=False).reset_index(drop=True).index[pc.sort_values("fc",ascending=False).reset_index(drop=True).abbr=="THCA"][0]+1
axB.axvline(0, color="black", lw=.6); axB.set_xlabel("Δ median PHLDA3 (tumour − normal)", fontsize=6.5)
axB.set_title(f"Pan-cancer over-expression (THCA rank {rank}/{len(pc)})", fontweight="bold", fontsize=7.5); axB.spines[["top","right"]].set_visible(False)

# (c) external GSE33630
geo = pd.read_csv("data/PHLDA3_GSE33630.csv")
order=[g for g in ["Normal","PTC","ATC"] if g in geo.group.unique()]
gg=[geo[geo.group==g].expr.values for g in order]; pk=kruskal(*gg)[1]
cmap={"Normal":ns.C["normal"],"PTC":ns.C["tumor"],"ATC":ns.C["accent2"]}
bp=axC.boxplot(gg, widths=.6, patch_artist=True, showfliers=False, medianprops=dict(color="black",lw=1))
for patch,g in zip(bp["boxes"],order): patch.set_facecolor(cmap[g]); patch.set_alpha(.65)
for j,(v,g) in enumerate(zip(gg,order),1): axC.scatter(np.random.normal(j,.07,len(v)),v,s=8,color=cmap[g],ec="black",lw=.2,zorder=3)
axC.set_xticks(range(1,len(order)+1)); axC.set_xticklabels([f"{g}\n(n={len(v)})" for g,v in zip(order,gg)],fontsize=6.3)
axC.set_ylabel("PHLDA3 expression", fontsize=6.5)
axC.set_title(f"External validation GSE33630 (KW P={pk:.0e})", fontweight="bold", fontsize=7); axC.spines[["top","right"]].set_visible(False)

# (d) tumour-vs-normal ROC (confirmation only) — from htseq
df=pd.read_csv("data/TCGA-THCA.htseq_counts.tsv",sep="\t",index_col=0)
mat=(2**df-1).clip(lower=0); row=mat.index[mat.index.str.startswith(ENS)][0]
lc=np.log2(mat.loc[row]/mat.sum(axis=0)*1e6+1)
st=pd.Series([b.split("-")[3][:2] for b in df.columns],index=df.columns)
y=(st=="01").astype(int)[st.isin(["01","11"])]; x=lc[st.isin(["01","11"])]
auc=roc_auc_score(y,x); fpr,tpr,_=roc_curve(y,x)
axD.plot(fpr,tpr,color=ns.C["tumor"],lw=1.6); axD.plot([0,1],[0,1],"--",color=ns.C["grey_mid"],lw=.8)
axD.text(.55,.12,f"AUC = {auc:.3f}",fontsize=7,color=ns.C["accent2"])
axD.set_xlabel("1 − Specificity",fontsize=6.5); axD.set_ylabel("Sensitivity",fontsize=6.5)
axD.set_title("Tumour-vs-normal ROC\n(confirmation, not a diagnostic test)",fontweight="bold",fontsize=7); axD.spines[["top","right"]].set_visible(False)

for ax,l in zip([axA,axB,axC,axD],"abcd"): ns.panel_label(ax,l,x=-0.16,y=1.06,fs=10)
fig.suptitle("Figure 1. PHLDA3 is over-expressed in thyroid carcinoma and validated externally",
             fontweight="bold",fontsize=8,y=1.0)
ns.save_fig(fig,"PHLDA3_THCA_Figure1")
print(f"done Figure 1 composite; THCA tumour-vs-normal AUC={auc:.3f}")
