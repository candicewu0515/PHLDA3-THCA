#!/usr/bin/env python3
"""Composite main Figure 4 (immune microenvironment, purity-adjusted) — Nature-style:
(a) PHLDA3 vs immune-cell infiltration, (b) PHLDA3 vs immune-checkpoint genes,
each as raw vs tumour-purity-adjusted (ESTIMATE partial Spearman) coefficients,
associations surviving purity adjustment highlighted. All from cached CSVs."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nature_style as ns
ns.set_style(font_size=7)

pp = pd.read_csv("data/PHLDA3_purity_partialcorr.csv")   # item,r_raw,p_raw,r_partial,p_partial,survives,type
fig, ax = plt.subplots(1, 2, figsize=(7.2, 4.2))

def panel(a, sub, title, topn=12):
    d = sub.reindex(sub.r_partial.abs().sort_values().index).tail(topn)
    y = range(len(d))
    for i,(_,r) in enumerate(d.iterrows()):
        a.plot([r.r_raw, r.r_partial], [i, i], color=ns.C["ns"], lw=1.4, zorder=1)
        a.scatter(r.r_raw, i, s=22, facecolor="white", ec=ns.C["normal"], lw=1.0, zorder=2)
        col = ns.C["tumor"] if r.survives else ns.C["ns"]
        a.scatter(r.r_partial, i, s=34, color=col, ec="black", lw=.4, zorder=3)
    a.axvline(0, ls="--", color=ns.C["grey_mid"], lw=.8)
    a.set_yticks(list(y)); a.set_yticklabels(d.item, fontsize=6.3)
    for lbl, surv in zip(a.get_yticklabels(), d.survives):
        if surv: lbl.set_fontweight("bold")
    a.set_xlabel("Spearman r with PHLDA3", fontsize=6.5)
    a.set_title(title, fontweight="bold", fontsize=7.5); a.spines[["top","right"]].set_visible(False)

panel(ax[0], pp[pp.type=="immune_cell"], "Immune-cell infiltration (ssGSEA)")
panel(ax[1], pp[pp.type=="checkpoint"], "Immune-checkpoint genes")
# shared legend
from matplotlib.lines import Line2D
leg=[Line2D([0],[0],marker='o',color='w',markerfacecolor='white',markeredgecolor=ns.C["normal"],label='raw',ms=6),
     Line2D([0],[0],marker='o',color='w',markerfacecolor=ns.C["tumor"],markeredgecolor='black',label='purity-adjusted (survives)',ms=6),
     Line2D([0],[0],marker='o',color='w',markerfacecolor=ns.C["ns"],markeredgecolor='black',label='purity-adjusted (n.s.)',ms=6)]
ax[1].legend(handles=leg, frameon=False, fontsize=5.6, loc="lower right")
for a,l in zip(ax,"ab"): ns.panel_label(a,l,x=-0.2)
fig.suptitle("Figure 4. PHLDA3 and the immune microenvironment, adjusted for tumour purity",
             fontweight="bold", fontsize=8, y=1.0)
plt.tight_layout(); ns.save_fig(fig, "PHLDA3_THCA_Figure4")
print("done Figure 4 composite")
