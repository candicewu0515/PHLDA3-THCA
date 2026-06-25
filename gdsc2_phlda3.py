#!/usr/bin/env python3
"""Task 7 (drug sensitivity) — PHLDA3 expression vs GDSC2 drug response across 805
cancer cell lines. Spearman of PHLDA3 (RMA) vs log(IC50); negative = PHLDA3-high
more sensitive. Complements the L1000CDS2 signature-reversal analysis."""
import pandas as pd, numpy as np, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nature_style as ns
ns.set_style(font_size=8)

d = pd.read_csv("gdsc2_phlda3_drugcorr.csv")
d["name"] = d["drug"].str.replace(r"_\d+$", "", regex=True)
sig = d[d["fdr"] < 0.05].copy()
MAPK = {"Trametinib","Selumetinib","SCH772984","PD0325901","Refametinib","ERK_2440","ERK_6604","RO-4987655","Sapitinib","Dasatinib"}
def col(name, rho):
    if name in MAPK: return ns.C["accent"]
    return ns.C["normal"] if rho < 0 else ns.C["tumor"]

sens = sig.sort_values("rho").head(12)            # most sensitizing (rho<0)
resist = sig.sort_values("rho", ascending=False).head(8)  # most resistance (rho>0)
plot = pd.concat([resist[::-1], sens[::-1]])      # resistance top, sensitizing bottom

fig, ax = plt.subplots(figsize=(6.8, 5.6))
cols = [col(n, r) for n, r in zip(plot["name"], plot["rho"])]
ax.barh(range(len(plot)), plot["rho"], color=cols, edgecolor="black", lw=0.4)
ax.set_yticks(range(len(plot)))
ax.set_yticklabels(plot["name"], fontsize=7.5)
for lbl in ax.get_yticklabels():
    if lbl.get_text() in MAPK: lbl.set_color(ns.C["accent"]); lbl.set_fontweight("bold")
ax.axvline(0, color="black", lw=0.8)
ax.set_xlabel("Spearman r: PHLDA3 expression vs log(IC50)  (GDSC2, 805 cell lines)")
ax.set_title("PHLDA3-high cell lines are sensitive to MEK/ERK inhibitors\n"
             "(left = more sensitive · teal = MAPK-pathway inhibitor; FDR<0.05)",
             fontweight="bold", fontsize=8.5)
ax.text(0.02, 0.02, "← more sensitive", transform=ax.transAxes, fontsize=7, color=ns.C["normal"])
ax.text(0.98, 0.02, "more resistant →", transform=ax.transAxes, fontsize=7, color=ns.C["tumor"], ha="right")
ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_gdsc2")
print(f"sensitizing (FDR<0.05): {(sig.rho<0).sum()} | top MEK/ERK: Trametinib, SCH772984, Selumetinib")
print("done GDSC2")
