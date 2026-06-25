#!/usr/bin/env python3
"""Task 6 (single-cell deep) — does PHLDA3 track p53-program activity within
malignant thyroid cells? Per-cell PHLDA3 vs a p53-target module score
(GSE193581). A lightweight, database-free read-out of the SCENIC TP53-regulon
question, testing the bulk p53 link at single-cell resolution."""
import pandas as pd, numpy as np, glob, os, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import nature_style as ns
ns.set_style(font_size=8)

GENE = "PHLDA3"
# curated canonical p53 transcriptional targets (PHLDA3 itself excluded)
P53 = ["CDKN1A","MDM2","BBC3","BAX","GADD45A","SESN1","SESN2","FAS","TNFRSF10B",
       "RRM2B","BTG2","DDB2","POLH","TP53I3","ZMAT3","PLK3","TIGAR","AEN"]

anno = pd.read_csv("gse193581_anno.txt.gz", sep="\t", header=None, skiprows=1,
                   names=["barcode","sample","celltype"])
bc2ct = dict(zip(anno.barcode, anno.celltype))

recs = []
for f in sorted(glob.glob("gse193581_raw/*_UMI.txt.gz")):
    sample = re.search(r"_([A-Z]+\d+)_UMI", os.path.basename(f)).group(1)
    grp = "PTC" if sample.startswith("PTC") else "ATC" if sample.startswith("ATC") else "NORM"
    df = pd.read_csv(f, sep="\t", index_col=0)
    total = df.sum(axis=0).replace(0, np.nan)
    norm = np.log1p(df.div(total, axis=1) * 1e4)          # log CP10k
    if GENE not in norm.index: continue
    p53g = [g for g in P53 if g in norm.index]
    p53score = norm.loc[p53g].mean(axis=0)                # p53-target module score
    phl = norm.loc[GENE]
    for bc in df.columns:
        ct = bc2ct.get(bc)
        if ct is None: continue
        recs.append((sample, grp, ct, float(phl[bc]), float(p53score[bc])))
d = pd.DataFrame(recs, columns=["sample","grp","celltype","PHLDA3","p53"])
d.to_csv("PHLDA3_sc_p53.csv", index=False)

mal = d[d.celltype=="Malignant cell"]
rho_m, p_m = spearmanr(mal.PHLDA3, mal.p53)
rho_all, p_all = spearmanr(d.PHLDA3, d.p53)
print(f"malignant cells n={len(mal)}  PHLDA3 vs p53-score rho={rho_m:.2f} p={p_m:.1e}")
# per cell type
ct_rho = (d.groupby("celltype").apply(lambda x: spearmanr(x.PHLDA3,x.p53)[0] if len(x)>50 else np.nan)
          .dropna().sort_values())
print(ct_rho.round(2).to_string())

# ================= figure =================
fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
# A: malignant-cell density scatter PHLDA3 vs p53 score
m = mal.sample(min(8000,len(mal)), random_state=0)
ax[0].scatter(m.p53, m.PHLDA3, s=3, color=ns.C["tumor"], alpha=0.18, edgecolors="none")
b,a = np.polyfit(mal.p53, mal.PHLDA3, 1); xs=np.linspace(mal.p53.min(),mal.p53.max(),50)
ax[0].plot(xs, b*xs+a, color="black", lw=1.3)
ax[0].set_xlabel("p53-target module score (log CP10k)"); ax[0].set_ylabel(f"{GENE} (log CP10k)")
ax[0].set_title(f"Malignant cells: PHLDA3 vs p53 program\nSpearman r={rho_m:.2f}, P<1e-300" if p_m==0
                else f"Malignant cells: PHLDA3 vs p53 program\nSpearman r={rho_m:.2f}, P={p_m:.1e}",
                fontweight="bold", fontsize=9)
ax[0].spines[["top","right"]].set_visible(False)
# B: per-cell-type correlation bar (malignant highlighted)
cols=[ns.C["highlight"] if c=="Malignant cell" else ns.C["ns"] for c in ct_rho.index]
ax[1].barh(range(len(ct_rho)), ct_rho.values, color=cols, edgecolor="black", lw=0.4)
ax[1].set_yticks(range(len(ct_rho))); ax[1].set_yticklabels(ct_rho.index, fontsize=8)
ax[1].axvline(0,color="black",lw=0.8)
ax[1].set_xlabel("Spearman r (PHLDA3 vs p53 score)")
ax[1].set_title("PHLDA3–p53 coupling by cell type", fontweight="bold", fontsize=9)
ax[1].spines[["top","right"]].set_visible(False)
for a_,l in zip(ax,["a","b"]): ns.panel_label(a_,l)
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_sc_p53")
print("done task 6 (single-cell p53 module)")
