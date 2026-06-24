#!/usr/bin/env python3
"""PHLDA3 in THCA TUMOR single cells (GEO GSE193581: PTC/ATC/normal, author-annotated,
incl. 16,049 malignant cells). Confirms PHLDA3 expression in malignant thyroid cells."""
import pandas as pd, numpy as np, glob, os, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import nature_style as ns
ns.set_style(font_size=8)

GENE = "PHLDA3"

# ---------- annotation: barcode -> celltype ----------
anno = pd.read_csv("gse193581_anno.txt.gz", sep="\t", header=None, skiprows=1,
                   names=["barcode","sample","celltype"])
bc2ct = dict(zip(anno["barcode"], anno["celltype"]))
print("annotated cells:", len(anno), "| celltypes:", anno.celltype.nunique())

# ---------- per-sample: CP10k-normalized PHLDA3 per cell ----------
recs = []
for f in sorted(glob.glob("gse193581_raw/*_UMI.txt.gz")):
    sample = re.search(r"_([A-Z]+\d+)_UMI", os.path.basename(f)).group(1)
    grp = "PTC" if sample.startswith("PTC") else "ATC" if sample.startswith("ATC") else "NORM"
    df = pd.read_csv(f, sep="\t", index_col=0)
    total = df.sum(axis=0).replace(0, np.nan)              # UMI per cell
    if GENE not in df.index:
        print("  ! no PHLDA3 in", sample); continue
    norm = np.log1p(df.loc[GENE] / total * 1e4)
    for bc, v, raw in zip(df.columns, norm.values, df.loc[GENE].values):
        ct = bc2ct.get(bc)
        if ct is None: continue                            # keep only annotated (QC-passed) cells
        recs.append((bc, sample, grp, ct, v, int(raw)))
    print(f"  {sample:7s} grp={grp} cells={df.shape[1]} annotated={sum(1 for bc in df.columns if bc in bc2ct)}")

d = pd.DataFrame(recs, columns=["barcode","sample","grp","celltype","phlda3","raw"])
d.to_csv("PHLDA3_tumor_sc.csv", index=False)
print("merged cells:", len(d), "| by group:", d.grp.value_counts().to_dict())

# ================= figure =================
fig, ax = plt.subplots(1, 2, figsize=(14, 6))

# ---- A: dotplot PHLDA3 by cell type ----
g = d.groupby("celltype").agg(n=("phlda3","size"),
        pct=("raw", lambda x:(x>0).mean()*100), mean=("phlda3","mean")).reset_index()
g = g[g.n>=30].sort_values("pct").reset_index(drop=True)
norm = Normalize(vmin=0, vmax=g["mean"].max()); cmap=ns.SEQ_RED
sizes = 50 + 320*(g["n"]/g["n"].max())**0.5
for i,(_,r) in enumerate(g.iterrows()):
    mal = r.celltype=="Malignant cell"
    ax[0].scatter(r.pct, i, s=sizes.iloc[i], c=[cmap(norm(r["mean"]))],
                  edgecolors=ns.C["highlight"] if mal else "black", linewidths=2.2 if mal else 0.6, zorder=3)
ax[0].set_yticks(range(len(g)))
ax[0].set_yticklabels([f"{c}  (n={n})" for c,n in zip(g.celltype,g.n)], fontsize=9)
for lbl in ax[0].get_yticklabels():
    if lbl.get_text().startswith("Malignant"): lbl.set_color(ns.C["highlight"]); lbl.set_fontweight("bold")
ax[0].set_xlabel("% cells expressing PHLDA3"); ax[0].grid(axis="x", ls=":", alpha=0.5)
ax[0].set_title("PHLDA3 across cell types\n(GSE193581 thyroid tumors)", fontweight="bold", fontsize=11)
ax[0].spines[["top","right"]].set_visible(False)
sm=ScalarMappable(norm=norm,cmap=cmap); sm.set_array([])
fig.colorbar(sm, ax=ax[0], fraction=0.04, pad=0.02).set_label("Mean expr (log CP10k)")

# ---- B: PHLDA3 in malignant cells (PTC/ATC) vs normal epithelial ----
mal = d[d.celltype=="Malignant cell"]
groups, labels, cols = [], [], []
for sub,lab,c in [(mal[mal.grp=="PTC"],"PTC\nmalignant",ns.C["tumor"]),
                  (mal[mal.grp=="ATC"],"ATC\nmalignant",ns.C["accent2"]),
                  (d[(d.celltype=="Epithelial cell")&(d.grp=="NORM")],"Normal\nepithelial",ns.C["normal"])]:
    if len(sub)>=10: groups.append(sub.phlda3.values); labels.append(lab); cols.append(c)
bp=ax[1].boxplot(groups, widths=0.6, patch_artist=True, showfliers=False,
                 medianprops=dict(color="black",lw=1.3))
for patch,c in zip(bp["boxes"],cols): patch.set_facecolor(c); patch.set_alpha(0.65)
for j,(gg,c) in enumerate(zip(groups,cols),1):
    ax[1].scatter(np.random.normal(j,0.07,min(len(gg),800)),
                  np.random.choice(gg,min(len(gg),800),replace=False), s=4, color=c, alpha=0.25, edgecolors="none")
ax[1].set_xticks(range(1,len(groups)+1))
ax[1].set_xticklabels([f"{l}\n(n={len(g)})" for l,g in zip(labels,groups)])
ax[1].set_ylabel(f"{GENE} expression (log CP10k)")
# stats: PTC malignant vs normal epithelial
if len(groups)>=3:
    p=mannwhitneyu(groups[0],groups[2],alternative="two-sided")[1]
    pstr=f"PTC mal vs normal epi: p={p:.1e}" if p>=1e-4 else "PTC mal vs normal epi: p<1e-4"
    ax[1].set_title(f"PHLDA3 in malignant thyroid cells\n{pstr}", fontweight="bold", fontsize=11)
ax[1].spines[["top","right"]].set_visible(False)
for i, lab in enumerate(["a","b"]):
    ns.panel_label(ax[i], lab)
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_tumor_singlecell")
# print key stats
print("\n=== PHLDA3 by cell type (% expressing) ===")
print(g.sort_values("pct",ascending=False)[["celltype","n","pct","mean"]].round(2).to_string(index=False))
print("\nMalignant cells: n=%d, %%expr=%.1f%%, mean=%.2f" % (len(mal),(mal.raw>0).mean()*100, mal.phlda3.mean()))
print("+ PHLDA3_tumor_sc.csv")
