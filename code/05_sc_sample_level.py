#!/usr/bin/env python3
"""Reviewer Major-3 — sample-level (not cell-level) test of PHLDA3 in malignant
vs non-malignant epithelial cells (GSE193581). Aggregates per sample first to
avoid single-cell pseudoreplication; the cell-level comparison is descriptive."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import nature_style as ns
ns.set_style(font_size=8)

d = pd.read_csv("data/PHLDA3_tumor_sc.csv")   # barcode, sample, grp, celltype, phlda3, raw
# per-sample mean expression & % detected, malignant (tumour samples) vs epithelial (normal samples)
mal = (d[d.celltype=="Malignant cell"].groupby("sample")
         .agg(mean=("phlda3","mean"), pct=("raw",lambda x:(x>0).mean()*100), n=("phlda3","size")))
epi = (d[(d.celltype=="Epithelial cell")&(d.grp=="NORM")].groupby("sample")
         .agg(mean=("phlda3","mean"), pct=("raw",lambda x:(x>0).mean()*100), n=("phlda3","size")))
mal=mal[mal.n>=20]; epi=epi[epi.n>=20]
p_mean = mannwhitneyu(mal["mean"], epi["mean"], alternative="two-sided")[1]
p_pct  = mannwhitneyu(mal["pct"],  epi["pct"],  alternative="two-sided")[1]
print(f"sample-level: malignant samples n={len(mal)} (PHLDA3+ {mal.pct.mean():.0f}%), "
      f"normal-epithelial samples n={len(epi)} ({epi.pct.mean():.1f}%)")
print(f"  mean expr  MWU P={p_mean:.3f};  %detected MWU P={p_pct:.3f}")
pd.concat([mal.assign(group="Malignant"), epi.assign(group="Normal epithelial")]).to_csv("data/PHLDA3_sc_sample_level.csv")

# figure: sample-level dots (each point = one patient/sample)
fig, ax = plt.subplots(1, 2, figsize=(7.5, 4))
for a,(col,lab,pv) in zip(ax, [("pct","% cells PHLDA3+",p_pct),("mean","mean PHLDA3 (log CP10k)",p_mean)]):
    g=[epi[col].values, mal[col].values]; cols=[ns.C["normal"],ns.C["tumor"]]
    bp=a.boxplot(g,widths=.55,patch_artist=True,showfliers=False,medianprops=dict(color="black",lw=1.2))
    for patch,cc in zip(bp["boxes"],cols): patch.set_facecolor(cc); patch.set_alpha(.5)
    for j,(gg,cc) in enumerate(zip(g,cols),1):
        a.scatter(np.random.normal(j,.07,len(gg)),gg,s=22,color=cc,ec="black",lw=.4,zorder=3)
    a.set_xticks([1,2]); a.set_xticklabels([f"Normal epi\n({len(epi)} samples)",f"Malignant\n({len(mal)} samples)"])
    a.set_ylabel(lab); a.set_title(f"Sample-level  (MWU P={pv:.3f})",fontweight="bold",fontsize=8.5)
    a.spines[["top","right"]].set_visible(False)
fig.suptitle("PHLDA3 in malignant vs normal epithelial cells — per-sample (each point = one sample)",
             fontweight="bold",fontsize=9,y=1.0)
plt.tight_layout(); ns.save_fig(fig,"PHLDA3_THCA_sc_sample_level")
print("done sample-level single-cell")
