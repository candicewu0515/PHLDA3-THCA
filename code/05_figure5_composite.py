#!/usr/bin/env python3
"""Composite main Figure 5 (single cell) — Nature-style multi-panel:
(a) PHLDA3 by cell type, (b) malignant vs normal epithelial, (c) PHLDA3 vs p53
module in malignant cells, (d) PHLDA3-p53 coupling by cell type,
(e) sample-level % detected, (f) sample-level mean. All from cached per-cell CSVs."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import nature_style as ns
ns.set_style(font_size=7)

d   = pd.read_csv("data/PHLDA3_tumor_sc.csv")          # per cell: sample,grp,celltype,phlda3,raw
p53 = pd.read_csv("data/PHLDA3_sc_p53.csv")            # per cell: celltype,PHLDA3,p53
samp= pd.read_csv("data/PHLDA3_sc_sample_level.csv")   # per sample: mean,pct,n,group

fig = plt.figure(figsize=(7.2, 6.4))
gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 1, 1], hspace=0.6, wspace=0.42)
axa = fig.add_subplot(gs[0, 0]); axb = fig.add_subplot(gs[0, 1])
axc = fig.add_subplot(gs[1, 0]); axd = fig.add_subplot(gs[1, 1])
axe = fig.add_subplot(gs[2, 0]); axf = fig.add_subplot(gs[2, 1])

# (a) % expressing by cell type (dotplot, malignant highlighted)
g = (d.groupby("celltype").agg(n=("raw","size"), pct=("raw",lambda x:(x>0).mean()*100),
        mean=("phlda3","mean")).reset_index())
g = g[g.n>=30].sort_values("pct").reset_index(drop=True)
norm=Normalize(0,g["mean"].max()); cmap=ns.SEQ_RED
sizes=40+220*(g["n"]/g["n"].max())**0.5
for i,(_,r) in enumerate(g.iterrows()):
    mal=r.celltype=="Malignant cell"
    axa.scatter(r.pct,i,s=sizes.iloc[i],c=[cmap(norm(r["mean"]))],
                ec=ns.C["highlight"] if mal else "black",lw=1.8 if mal else .5,zorder=3)
axa.set_yticks(range(len(g))); axa.set_yticklabels([f"{c} ({n:,})" for c,n in zip(g.celltype,g.n)],fontsize=6)
for lbl in axa.get_yticklabels():
    if lbl.get_text().startswith("Malignant"): lbl.set_color(ns.C["highlight"]); lbl.set_fontweight("bold")
axa.set_xlabel("% cells PHLDA3+"); axa.set_title("Cell-type expression (GSE193581)",fontweight="bold",fontsize=7.5)
axa.grid(axis="x",ls=":",alpha=.5); axa.spines[["top","right"]].set_visible(False)

# (b) malignant PTC/ATC vs normal epithelial (cell-level box, descriptive)
mal=d[d.celltype=="Malignant cell"]
gb=[(d[(d.celltype=="Epithelial cell")&(d.grp=="NORM")].phlda3.values,"Normal\nepi",ns.C["normal"]),
    (mal[mal.grp=="PTC"].phlda3.values,"PTC\nmal",ns.C["tumor"]),
    (mal[mal.grp=="ATC"].phlda3.values,"ATC\nmal",ns.C["accent2"])]
bp=axb.boxplot([x[0] for x in gb],widths=.6,patch_artist=True,showfliers=False,medianprops=dict(color="black",lw=1))
for patch,x in zip(bp["boxes"],gb): patch.set_facecolor(x[2]); patch.set_alpha(.65)
axb.set_xticks([1,2,3]); axb.set_xticklabels([x[1] for x in gb],fontsize=6.5)
axb.set_ylabel("PHLDA3 (log CP10k)"); axb.set_title("Malignant vs normal (cell level)",fontweight="bold",fontsize=7.5)
axb.spines[["top","right"]].set_visible(False)

# (c) PHLDA3 vs p53 module in malignant cells
m=p53[p53.celltype=="Malignant cell"]; rho,_=spearmanr(m.PHLDA3,m.p53)
ms=m.sample(min(6000,len(m)),random_state=0)
axc.scatter(ms.p53,ms.PHLDA3,s=2,color=ns.C["tumor"],alpha=.15,ec="none")
b,a=np.polyfit(m.p53,m.PHLDA3,1); xs=np.linspace(m.p53.min(),m.p53.max(),40)
axc.plot(xs,b*xs+a,color="black",lw=1.1)
axc.set_xlabel("p53-target module score"); axc.set_ylabel("PHLDA3 (log CP10k)")
axc.set_title(f"Malignant cells: PHLDA3 vs p53 (r={rho:.2f})",fontweight="bold",fontsize=7.5)
axc.spines[["top","right"]].set_visible(False)

# (d) PHLDA3-p53 coupling by cell type
ct=(p53.groupby("celltype").apply(lambda x: spearmanr(x.PHLDA3,x.p53)[0] if len(x)>50 else np.nan)
    .dropna().sort_values())
cols=[ns.C["highlight"] if c=="Malignant cell" else ns.C["ns"] for c in ct.index]
axd.barh(range(len(ct)),ct.values,color=cols,ec="black",lw=.4)
axd.set_yticks(range(len(ct))); axd.set_yticklabels(ct.index,fontsize=6); axd.axvline(0,color="black",lw=.7)
axd.set_xlabel("Spearman r (PHLDA3 vs p53)"); axd.set_title("p53 coupling by cell type",fontweight="bold",fontsize=7.5)
axd.spines[["top","right"]].set_visible(False)

# (e,f) sample-level (each point = one sample)
mal_s=samp[samp.group=="Malignant"]; epi_s=samp[samp.group=="Normal epithelial"]
for ax,col,lab in [(axe,"pct","% cells PHLDA3+"),(axf,"mean","mean PHLDA3 (log CP10k)")]:
    gg=[epi_s[col].values,mal_s[col].values]; cc=[ns.C["normal"],ns.C["tumor"]]
    bp=ax.boxplot(gg,widths=.55,patch_artist=True,showfliers=False,medianprops=dict(color="black",lw=1))
    for patch,c in zip(bp["boxes"],cc): patch.set_facecolor(c); patch.set_alpha(.5)
    for j,(g2,c) in enumerate(zip(gg,cc),1): ax.scatter(np.random.normal(j,.06,len(g2)),g2,s=14,color=c,ec="black",lw=.3,zorder=3)
    p=mannwhitneyu(gg[0],gg[1])[1]
    ax.set_xticks([1,2]); ax.set_xticklabels([f"Normal\n({len(epi_s)})",f"Malig\n({len(mal_s)})"],fontsize=6.5)
    ax.set_ylabel(lab,fontsize=6.5); ax.set_title(f"Sample-level (P={p:.3f})",fontweight="bold",fontsize=7.5)
    ax.spines[["top","right"]].set_visible(False)

for ax,l in zip([axa,axb,axc,axd,axe,axf],"abcdef"): ns.panel_label(ax,l,x=-0.22,fs=10)
sm=ScalarMappable(norm=norm,cmap=cmap); sm.set_array([])
fig.colorbar(sm,ax=axa,fraction=.045,pad=.02).set_label("mean expr",fontsize=6)
fig.suptitle("Figure 5. PHLDA3 is induced in malignant thyroid cells and tracks p53-programme activity",
             fontweight="bold",fontsize=8.5,y=1.005)
ns.save_fig(fig,"PHLDA3_THCA_Figure5")
print("done Figure 5 composite")
