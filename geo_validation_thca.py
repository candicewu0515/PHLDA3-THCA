#!/usr/bin/env python3
"""External validation of PHLDA3 up-regulation in thyroid cancer:
GEO GSE33630 (GPL570, normal vs PTC vs ATC). Boxplot + diagnostic ROC."""
import pandas as pd, numpy as np, gzip, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, kruskal
from sklearn.metrics import roc_curve, roc_auc_score
import nature_style as ns
ns.set_style(font_size=8)

GENE, PROBE = "PHLDA3", "218634_at"
F = "GSE33630.txt.gz"

# ---------- parse series matrix ----------
gsm, diag, tissue, expr_line = None, None, None, None
with gzip.open(F, "rt") as fh:
    in_table = False
    for ln in fh:
        if ln.startswith("!Sample_geo_accession"):
            gsm = [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
        elif ln.startswith("!Sample_characteristics_ch1"):
            vals = [x.strip('"') for x in ln.rstrip("\n").split("\t")[1:]]
            if any("diagnostic" in v for v in vals): diag = vals
            elif any("tissue type" in v for v in vals): tissue = vals
        elif ln.startswith("!series_matrix_table_begin"):
            in_table = True
        elif in_table and ln.startswith('"'+PROBE):
            expr_line = ln.rstrip("\n").split("\t"); break

vals = np.array([float(x.strip('"')) for x in expr_line[1:]])
tis = [t.split(":")[-1].strip().lower() for t in tissue]
dia = [d.split(":")[-1].strip() for d in diag]
def grp(t, d):
    if "non-tumor" in t or "non-tumor" in d.lower() or "normal" in t: return "Normal"
    if "anaplastic" in d.lower() or "ATC" in d: return "ATC"
    if "papillary" in d.lower() or "PTC" in d: return "PTC"
    return "Tumor" if "tumour" in t or "tumor" in t else "other"
g = pd.DataFrame({"gsm":gsm, "expr":vals, "group":[grp(t,d) for t,d in zip(tis,dia)]})
print(g["group"].value_counts().to_dict())
g.to_csv("PHLDA3_GSE33630.csv", index=False)

normal = g[g.group=="Normal"]["expr"].values
ptc    = g[g.group=="PTC"]["expr"].values
atc    = g[g.group=="ATC"]["expr"].values
tumor  = g[g.group.isin(["PTC","ATC"])]["expr"].values
p_tn = mannwhitneyu(tumor, normal, alternative="two-sided")[1]
p_kw = kruskal(normal, ptc, atc)[1]
# diagnostic ROC tumor vs normal
y = np.r_[np.ones(len(tumor)), np.zeros(len(normal))]
x = np.r_[tumor, normal]
auc = roc_auc_score(y, x); fpr,tpr,_ = roc_curve(y,x)
print(f"Normal={len(normal)} PTC={len(ptc)} ATC={len(atc)} | tumor-vs-normal MWU p={p_tn:.2e} AUC={auc:.3f}")

# ---------- figure ----------
fig, ax = plt.subplots(1, 2, figsize=(9.5, 4.4))
order=["Normal","PTC","ATC"]; cols=[ns.C["normal"], ns.C["tumor"], ns.C["accent2"]]
groups=[normal,ptc,atc]
bp=ax[0].boxplot(groups, widths=0.6, patch_artist=True, showfliers=False,
                 medianprops=dict(color="black",lw=1.3))
for patch,c in zip(bp["boxes"],cols): patch.set_facecolor(c); patch.set_alpha(0.6)
for j,(gg,c) in enumerate(zip(groups,cols),1):
    ax[0].scatter(np.random.normal(j,0.07,len(gg)), gg, s=10, color=c, alpha=0.5, edgecolors="none")
ax[0].set_xticks([1,2,3]); ax[0].set_xticklabels([f"{o}\n(n={len(gg)})" for o,gg in zip(order,groups)])
ax[0].set_ylabel(f"{GENE} expression (GPL570, log2)")
ax[0].set_title(f"GSE33630 external cohort\nKruskal p = {p_kw:.1e}", fontweight="bold", fontsize=11)
ax[0].spines[["top","right"]].set_visible(False)
ns.panel_label(ax[0], "a")

ax[1].plot(fpr,tpr,color=ns.C["tumor"],lw=2.2,label=f"AUC = {auc:.3f}")
ax[1].plot([0,1],[0,1],"--",color=ns.C["grey_mid"],lw=1)
ax[1].set_xlabel("1 - Specificity"); ax[1].set_ylabel("Sensitivity")
ax[1].set_title(f"Diagnostic ROC (tumor vs normal)\nMWU p = {p_tn:.1e}", fontweight="bold", fontsize=11)
ax[1].legend(loc="lower right", frameon=False); ax[1].spines[["top","right"]].set_visible(False)
ax[1].set_xlim(-0.02,1.02); ax[1].set_ylim(-0.02,1.02)
ns.panel_label(ax[1], "b")
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_geo_validation")
print("saved: PHLDA3_GSE33630.csv")
