#!/usr/bin/env python3
"""PHLDA3 vs immune landscape in THCA: ssGSEA immune-cell infiltration (Danaher
signatures) + immune-checkpoint gene correlation."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import gseapy as gp
import nature_style as ns
ns.set_style(font_size=8)

GENE = "PHLDA3"

# Danaher et al. 2017 immune cell signatures
IMMUNE = {
 "B cells": ["BLK","CD19","MS4A1","TNFRSF17","FCRL2","PNOC","SPIB","TCL1A"],
 "CD8 T cells": ["CD8A","CD8B"],
 "Cytotoxic cells": ["CTSW","GNLY","GZMA","GZMB","GZMH","KLRB1","KLRD1","KLRK1","NKG7","PRF1"],
 "DC": ["CCL13","CD209","HSD11B1"],
 "Exhausted CD8": ["CD244","EOMES","LAG3","PTGER4"],
 "Macrophages": ["CD163","CD68","CD84","MS4A4A"],
 "Mast cells": ["CPA3","HDC","MS4A2","TPSAB1","TPSB2"],
 "Neutrophils": ["CSF3R","S100A12","CEACAM3","FCGR3B","FFAR2"],
 "NK CD56dim": ["IL21R","KIR2DL3","KIR3DL1","KIR3DL2"],
 "NK cells": ["NCR1","XCL1","XCL2"],
 "T cells": ["CD3D","CD3E","CD3G","CD6","SH2D1A","TRAT1"],
 "Th1 cells": ["TBX21"],
 "Treg": ["FOXP3"],
}
CHECKPOINTS = ["PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT",
               "BTLA","IDO1","CD276","TNFRSF9","ICOS","CD27","LGALS9"]

# ---------- tumor logCPM matrix (symbols) ----------
df = pd.read_csv("data/TCGA-THCA.htseq_counts.tsv", sep="\t", index_col=0)
mat = (2**df - 1).clip(lower=0)
stype = pd.Series([b.split("-")[3][:2] for b in df.columns], index=df.columns)
tum = df.columns[stype == "01"]
cpm = mat[tum].div(mat[tum].sum(axis=0), axis=1) * 1e6
keep = cpm.mean(axis=1) >= 1
logcpm = np.log2(cpm[keep] + 1); logcpm.index = [g.split(".")[0] for g in logcpm.index]
e2s = pd.read_csv("data/thca_DE_full_symbol.csv")[["ensembl","symbol"]].dropna()
e2s = dict(zip(e2s.ensembl, e2s.symbol))
logcpm = logcpm[logcpm.index.isin(e2s)]; logcpm.index = [e2s[e] for e in logcpm.index]
logcpm = logcpm[~logcpm.index.duplicated()]
phl = logcpm.loc[GENE]
print("genes x tumors:", logcpm.shape)

# ---------- ssGSEA infiltration ----------
ss = gp.ssgsea(data=logcpm, gene_sets=IMMUNE, sample_norm_method="rank",
               min_size=1, outdir=None, no_plot=True)
nes = ss.res2d.pivot(index="Term", columns="Name", values="NES").astype(float)
nes = nes[logcpm.columns]
rows = []
for cell in nes.index:
    r, p = pearsonr(nes.loc[cell].values, phl.values)
    rows.append({"cell": cell, "r": r, "p": p})
inf = pd.DataFrame(rows).sort_values("r")
inf["sig"] = inf["p"].map(lambda p:"***" if p<1e-3 else "**" if p<1e-2 else "*" if p<.05 else "")
inf.to_csv("data/PHLDA3_immune_infiltration.csv", index=False)

# ---------- checkpoint correlation ----------
rows = []
for g in CHECKPOINTS:
    if g in logcpm.index:
        r, p = pearsonr(logcpm.loc[g].values, phl.values)
        rows.append({"gene": g, "r": r, "p": p})
cp = pd.DataFrame(rows).sort_values("r")
cp["sig"] = cp["p"].map(lambda p:"***" if p<1e-3 else "**" if p<1e-2 else "*" if p<.05 else "")
cp.to_csv("data/PHLDA3_immune_checkpoints.csv", index=False)
print("infiltration r range:", round(inf.r.min(),2), round(inf.r.max(),2))
print("checkpoint  r range:", round(cp.r.min(),2), round(cp.r.max(),2))

# ================= figure =================
fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.2))

def corr_bars(a, d, key, title):
    cols = [ns.C["tumor"] if v>0 else ns.C["normal"] for v in d["r"]]
    a.barh(range(len(d)), d["r"], color=cols, alpha=0.85, edgecolor="black", lw=0.4)
    a.set_yticks(range(len(d))); a.set_yticklabels(d[key], fontsize=9)
    for i,(v,s) in enumerate(zip(d["r"], d["sig"])):
        a.text(v + (0.01 if v>=0 else -0.01), i, s, va="center",
               ha="left" if v>=0 else "right", fontsize=9)
    a.axvline(0, color="black", lw=0.8)
    a.set_xlabel(f"Pearson r with {GENE}")
    a.set_title(title, fontweight="bold", fontsize=11)
    a.spines[["top","right"]].set_visible(False)
    m = max(0.05, d["r"].abs().max()*1.25); a.set_xlim(-m, m)

corr_bars(ax[0], inf, "cell", "Immune-cell infiltration (ssGSEA)")
corr_bars(ax[1], cp, "gene", "Immune-checkpoint genes")
for i, lab in enumerate(["a","b"]):
    ns.panel_label(ax[i], lab)
fig.suptitle(f"{GENE} and the immune microenvironment (TCGA-THCA, n={logcpm.shape[1]})",
             fontweight="bold", fontsize=12, y=1.0)
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_immune")
print("+ 2 CSVs")
