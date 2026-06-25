#!/usr/bin/env python3
"""THCA tumor-vs-normal whole-transcriptome volcano, PHLDA3 highlighted."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from adjustText import adjust_text
import nature_style as ns
ns.set_style(font_size=8)

res = pd.read_csv("data/thca_DE_full_symbol.csv")
res = res.dropna(subset=["log2FC_tumor_vs_normal", "FDR"]).copy()
res["fc"] = res["log2FC_tumor_vs_normal"]
res["nlogq"] = -np.log10(res["FDR"].clip(lower=1e-300))

FC_T, Q_T = 1.0, 0.05
res["cat"] = "ns"
res.loc[(res.fc >= FC_T) & (res.FDR < Q_T), "cat"] = "up"
res.loc[(res.fc <= -FC_T) & (res.FDR < Q_T), "cat"] = "down"
print(res["cat"].value_counts().to_dict())

COL = {"up": ns.C["tumor"], "down": ns.C["normal"], "ns": ns.C["ns"]}
fig, ax = plt.subplots(figsize=(6.4, 5.6))
for c in ["ns", "down", "up"]:
    s = res[res.cat == c]
    ax.scatter(s.fc, s.nlogq, s=7, c=COL[c], alpha=0.45 if c=="ns" else 0.6,
               edgecolors="none", rasterized=True, label=None)
ax.axvline(FC_T, ls="--", lw=0.8, color="grey"); ax.axvline(-FC_T, ls="--", lw=0.8, color="grey")
ax.axhline(-np.log10(Q_T), ls="--", lw=0.8, color="grey")

# label PHLDA3 + a handful of top hits each side
ph = res[res.symbol == "PHLDA3"].iloc[0]
ax.scatter(ph.fc, ph.nlogq, s=90, c="#000000", marker="*", zorder=6)
texts = [ax.annotate("PHLDA3", (ph.fc, ph.nlogq), fontsize=11, fontweight="bold", color="black")]
res["score"] = res["nlogq"] * res["fc"].abs()
for c in ["up", "down"]:
    for _, r in res[res.cat == c].nlargest(6, "score").iterrows():
        if r.symbol == "PHLDA3" or pd.isna(r.symbol): continue
        texts.append(ax.text(r.fc, r.nlogq, r.symbol, fontsize=8))
adjust_text(texts, arrowprops=dict(arrowstyle="-", color="grey", lw=0.5), ax=ax)

n_up = (res.cat=="up").sum(); n_dn = (res.cat=="down").sum()
ax.set_xlabel("log2 fold change (Tumor vs Normal)")
ax.set_ylabel("-log10 FDR")
ax.set_title(f"TCGA-THCA differential expression\n↑{n_up} up   ↓{n_dn} down  (|log2FC|≥1, FDR<0.05)",
             fontweight="bold", fontsize=11)
ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_volcano")
