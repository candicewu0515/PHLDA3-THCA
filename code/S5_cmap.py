#!/usr/bin/env python3
"""Task 7 (CMap, perturbagen-class) — clue.io CMap connectivity of the PHLDA3-high
signature, MOA-class level (gsea_result.gct from the clue.io result archive).
norm_cs < 0 = class REVERSES the PHLDA3-high program. Confirms L1000CDS2 (HDAC)."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nature_style as ns
ns.set_style(font_size=8)

agg = pd.read_csv("data/PHLDA3_cmap_moa_clean.csv")
agg.columns = ["cls", "norm_cs"]
agg = agg.sort_values("norm_cs").reset_index(drop=True)

REV = {"HDAC_INHIBITOR","HSP_INHIBITOR"}     # therapeutic reversers of interest
MIM = {"MEK_INHIBITOR","ERK_INHIBITOR","PI3K_INHIBITOR","AKT_INHIBITOR","RAF_INHIBITOR"}
# panel: top reversing classes + the highlighted mechanistic classes
top_rev = agg.head(12)
focus = agg[agg.cls.isin(REV | MIM)]
panel = pd.concat([top_rev, focus]).drop_duplicates("cls").sort_values("norm_cs")
def color(c, v):
    if c in REV: return ns.C["tumor"]
    if c in MIM: return ns.C["accent"]
    return ns.C["normal"] if v < 0 else ns.C["ns"]
def nice(c): return c.replace("_"," ").title()

fig, ax = plt.subplots(figsize=(7, 5.6))
cols = [color(c, v) for c, v in zip(panel.cls, panel.norm_cs)]
ax.barh(range(len(panel)), panel.norm_cs, color=cols, edgecolor="black", lw=0.4)
ax.set_yticks(range(len(panel)))
ax.set_yticklabels([nice(c) for c in panel.cls], fontsize=7)
for lbl in ax.get_yticklabels():
    t = lbl.get_text().upper().replace(" ","_")
    if t in REV: lbl.set_color(ns.C["tumor"]); lbl.set_fontweight("bold")
    elif t in MIM: lbl.set_color(ns.C["accent"]); lbl.set_fontweight("bold")
ax.axvline(0, color="black", lw=0.8)
ax.set_xlabel("CMap connectivity (norm_cs):  ← reverses · mimics →")
ax.set_title("clue.io CMap (MOA-class): HDAC/HSP90 inhibitors reverse the\n"
             "PHLDA3-high program (red), consistent with L1000CDS2",
             fontweight="bold", fontsize=8.5)
hd = agg[agg.cls=="HDAC_INHIBITOR"]
ax.text(0.02, 0.02, f"HDAC inhibitor: norm_cs={hd.norm_cs.iloc[0]:.2f} (reverser, rank {agg.index[agg.cls=='HDAC_INHIBITOR'][0]+1}/{len(agg)})",
        transform=ax.transAxes, fontsize=6.5, color=ns.C["tumor"])
ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_cmap")
print("HDAC_INHIBITOR norm_cs:", round(hd.norm_cs.iloc[0],2), "| top reverser:", agg.cls.iloc[0], round(agg.norm_cs.iloc[0],2))
print("done CMap clean")
