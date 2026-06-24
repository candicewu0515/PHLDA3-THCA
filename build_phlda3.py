#!/usr/bin/env python3
"""PHLDA3 in TCGA-THCA: tumor-vs-normal boxplot + diagnostic ROC. Needs only expression."""
import pandas as pd, numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from sklearn.metrics import roc_curve, roc_auc_score
import nature_style as ns
ns.set_style(font_size=8)

GENE, ENS = "PHLDA3", "ENSG00000174307"
df = pd.read_csv("TCGA-THCA.htseq_counts.tsv", sep="\t", index_col=0)
mat = (2**df - 1).clip(lower=0)
libsize = mat.sum(axis=0)                      # per-sample library size (all genes)

row = mat.index[mat.index.str.startswith(ENS)][0]
cpm = mat.loc[row] / libsize * 1e6
logcpm = np.log2(cpm + 1)

stype = pd.Series([b.split("-")[3][:2] for b in df.columns], index=df.columns)
tumor = logcpm[stype == "01"]
normal = logcpm[stype == "11"]
u, p = mannwhitneyu(tumor, normal, alternative="two-sided")
print(f"{GENE}: tumor n={len(tumor)} median={tumor.median():.2f} | normal n={len(normal)} median={normal.median():.2f}")
print(f"log2FC(T-N)={tumor.mean()-normal.mean():+.2f}  MWU p={p:.2e}")

# labels: 1=tumor, 0=normal ; AUC
y = np.r_[np.ones(len(tumor)), np.zeros(len(normal))]
x = np.r_[tumor.values, normal.values]
auc = roc_auc_score(y, x)
fpr, tpr, _ = roc_curve(y, x)
print(f"diagnostic AUC = {auc:.3f}")

# ---- figure ----
fig, ax = plt.subplots(1, 2, figsize=(9, 4.2))

# (a) boxplot
bp = ax[0].boxplot([normal.values, tumor.values], widths=0.55, patch_artist=True,
                   showfliers=False, medianprops=dict(color="black", linewidth=1.4))
for patch, c in zip(bp["boxes"], [ns.C["normal"], ns.C["tumor"]]):
    patch.set_facecolor(c); patch.set_alpha(0.65)
for i, (d, c) in enumerate(zip([normal.values, tumor.values], [ns.C["normal"], ns.C["tumor"]]), start=1):
    ax[0].scatter(np.random.normal(i, 0.06, len(d)), d, s=7, color=c, alpha=0.45, zorder=3, edgecolors="none")
ax[0].set_xticks([1, 2]); ax[0].set_xticklabels([f"Normal\n(n={len(normal)})", f"Tumor\n(n={len(tumor)})"])
ax[0].set_ylabel(f"{GENE} expression  log2(CPM+1)")
ax[0].set_title(f"{GENE} in TCGA-THCA", fontweight="bold")
ymax = max(tumor.max(), normal.max())
ax[0].plot([1, 2], [ymax+0.4]*2, color="black", lw=1)
pstr = "****" if p < 1e-4 else ("***" if p < 1e-3 else ("**" if p < 1e-2 else "*"))
ax[0].text(1.5, ymax+0.5, pstr, ha="center", fontsize=14)
ax[0].spines[["top", "right"]].set_visible(False)
ns.panel_label(ax[0], "a")

# (b) ROC
ax[1].plot(fpr, tpr, color=ns.C["tumor"], lw=2.2, label=f"AUC = {auc:.3f}")
ax[1].plot([0, 1], [0, 1], "--", color=ns.C["grey_mid"], lw=1)
ax[1].set_xlabel("1 - Specificity"); ax[1].set_ylabel("Sensitivity")
ax[1].set_title(f"Diagnostic ROC: {GENE}", fontweight="bold")
ax[1].legend(loc="lower right", frameon=False)
ax[1].spines[["top", "right"]].set_visible(False)
ax[1].set_xlim(-0.02, 1.02); ax[1].set_ylim(-0.02, 1.02)
ns.panel_label(ax[1], "b")

plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_fig1")
