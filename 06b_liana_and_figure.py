#!/usr/bin/env python3
"""06b — memory-light: run liana L-R from PHLDA3-high malignant cells, then build the
combined Fig.06 from cached results (avoids holding the full matrix + scanpy together).
Run after 06_singlecell_deep.py has produced results/06_sc_panel_<group>.parquet."""
import os, gc, argparse
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr
import nature_style as ns
ns.set_style(font_size=8)
from importlib import import_module
m6 = import_module("06_singlecell_deep")
GENE, P53_TARGETS, EPI_MARKERS, OTHER_MARKERS = m6.GENE, m6.P53_TARGETS, m6.EPI_MARKERS, m6.OTHER_MARKERS

ap = argparse.ArgumentParser(); ap.add_argument("--group", default="TUMOR"); args = ap.parse_args()
group = args.group
cache = f"results/06_sc_panel_{group}.parquet"
tag = "" if group == "TUMOR" else f"_{group}"

# -------- liana (subsample first, free matrix, then import scanpy/liana) --------
def liana_step():
    import glob
    out = f"results/06_liana_phlda3hi_{group}.csv"
    rng = np.random.default_rng(0)
    # build subsample directly from small per-sample caches (never load full matrix)
    parts = []
    for f in sorted(glob.glob("results/06_cache/*.parquet")):
        p = pd.read_parquet(f)
        if not len(p): continue
        g = p.iloc[0]["grp"] if "grp" in p.columns else ""
        if group == "TUMOR" and g not in ("PTC","ATC"): continue
        if group in ("PTC","ATC") and g != group: continue
        k = min(len(p), 800)
        parts.append(p.iloc[rng.choice(len(p), k, replace=False)])
    sm = pd.concat(parts, ignore_index=True); del parts; gc.collect()
    gcols = [c for c in sm.columns if c not in ("barcode","celltype","sample","grp")]
    # cap to 600 per celltype
    keep = []
    for ct, sub in sm.groupby("celltype"):
        kk = min(len(sub), 600); keep += list(rng.choice(sub.index, kk, replace=False))
    sm = sm.loc[keep].copy()
    malv = sm.loc[sm.celltype == "Malignant cell", GENE]
    thr = malv[malv > 0].median() if (malv > 0).any() else 0
    print("liana subsample cells:", len(sm), "genes:", len(gcols))
    grp = sm.celltype.astype(str).copy()
    mm = sm.celltype == "Malignant cell"
    grp[mm & (sm[GENE] > thr)] = "Malignant_PHLDA3hi"
    grp[mm & (sm[GENE] <= thr)] = "Malignant_PHLDA3lo"
    import anndata as ad, liana as li
    A = ad.AnnData(X=sm[gcols].values.astype("float32"),
                   obs=pd.DataFrame({"label": grp.values}, index=sm.barcode.values),
                   var=pd.DataFrame(index=gcols))
    A.layers["normalized"] = A.X.copy()
    li.mt.rank_aggregate(A, groupby="label", expr_prop=0.1, use_raw=False,
                         layer="normalized", verbose=False, n_perms=50, seed=0)
    r = A.uns["liana_res"]
    r = r[r.source == "Malignant_PHLDA3hi"].copy()
    mcol = "magnitude_rank" if "magnitude_rank" in r.columns else "lr_means"
    r = r.sort_values(mcol).head(40); r.to_csv(out, index=False)
    print("liana saved:", out, "rows:", len(r)); return r

liana_csv = f"results/06_liana_phlda3hi_{group}.csv"
liana_top = pd.read_csv(liana_csv) if os.path.exists(liana_csv) else liana_step()

# -------- read only needed columns for p53 program (memory-light) --------
import pyarrow.parquet as pq
allcols = pq.ParquetFile(cache).schema.names
p53 = [g for g in P53_TARGETS if g in allcols]
dl = pd.read_parquet(cache, columns=["barcode","celltype",GENE] + p53)
dl["p53_program"] = dl[p53].mean(axis=1)
mal = dl[dl.celltype == "Malignant cell"].copy(); mal["pos"] = mal[GENE] > 0
rho, prho = spearmanr(mal[GENE], mal.p53_program)
pos = mal[mal.pos].p53_program; neg = mal[~mal.pos].p53_program
pmw = mannwhitneyu(pos, neg, alternative="greater")[1]
mk = pd.read_csv(f"results/06_malignant_markers_{group}.csv")
del dl; gc.collect()

# ================= figure =================
ncol = 4 if (liana_top is not None and len(liana_top)) else 3
fig, ax = plt.subplots(1, ncol, figsize=(3.5*ncol, 4.2))
y = np.arange(len(mk))[::-1]; h = 0.36
ax[0].barh(y+h/2, mk.PHLDA3hi_malignant_mean, height=h, color=ns.C["tumor"], alpha=0.9, label="PHLDA3-hi malignant")
ax[0].barh(y-h/2, mk.immune_cells_mean, height=h, color=ns.C["normal"], alpha=0.9, label="Immune cells")
ax[0].set_yticks(y); ax[0].set_yticklabels(mk.marker_set, fontsize=8)
ax[0].set_xlabel("Mean expr (log CP10k)"); ax[0].legend(fontsize=6.5, loc="lower right")
ax[0].set_title("Lineage identity of\nPHLDA3-high cells", fontweight="bold", fontsize=9.5)
ax[0].spines[["top","right"]].set_visible(False); ns.panel_label(ax[0], "a")

bp = ax[1].boxplot([neg.values, pos.values], widths=0.55, patch_artist=True, showfliers=False,
                   medianprops=dict(color="black", lw=1.3))
for patch, c in zip(bp["boxes"], [ns.C["normal"], ns.C["tumor"]]):
    patch.set_facecolor(c); patch.set_alpha(0.65)
ax[1].set_xticks([1,2]); ax[1].set_xticklabels([f"PHLDA3-\n(n={len(neg)})", f"PHLDA3+\n(n={len(pos)})"])
ax[1].set_ylabel("p53 target-gene program")
ax[1].set_title(f"p53 program in malignant cells\nMWU p={pmw:.1e}", fontweight="bold", fontsize=9.5)
ax[1].spines[["top","right"]].set_visible(False); ns.panel_label(ax[1], "b")

s = mal.sample(min(len(mal), 3000), random_state=0)
ax[2].scatter(s[GENE], s.p53_program, s=4, color=ns.C["tumor"], alpha=0.2, edgecolors="none")
mm, bb = np.polyfit(mal[GENE], mal.p53_program, 1)
xs = np.linspace(mal[GENE].min(), mal[GENE].max(), 50)
ax[2].plot(xs, mm*xs+bb, color=ns.C["grey_dark"], lw=1.4)
ax[2].set_xlabel(f"{GENE} (log CP10k)"); ax[2].set_ylabel("p53 target-gene program")
ax[2].set_title(f"p53 program vs {GENE}\n(malignant) rho={rho:+.2f} p={prho:.1e}", fontweight="bold", fontsize=9.5)
ax[2].spines[["top","right"]].set_visible(False); ns.panel_label(ax[2], "c")

if ncol == 4:
    r = liana_top.copy()
    r["pair"] = r["ligand_complex"] + "->" + r["receptor_complex"] + "\n(to " + r["target"].astype(str) + ")"
    mcol = "magnitude_rank" if "magnitude_rank" in r.columns else "lr_means"
    r = r.sort_values(mcol).head(10).iloc[::-1]
    score = -np.log10(r[mcol].clip(lower=1e-6)) if mcol == "magnitude_rank" else r[mcol]
    ax[3].barh(range(len(r)), score, color=ns.C["accent2"], alpha=0.9)
    ax[3].set_yticks(range(len(r))); ax[3].set_yticklabels(r["pair"], fontsize=6)
    ax[3].set_xlabel("-log10 magnitude_rank" if mcol == "magnitude_rank" else "LR mean")
    ax[3].set_title("Top signalling from\nPHLDA3-high malignant cells", fontweight="bold", fontsize=9.5)
    ax[3].spines[["top","right"]].set_visible(False); ns.panel_label(ax[3], "d")

fig.suptitle(f"PHLDA3 in malignant thyroid cells - p53 program & cell communication (GSE193581, {group})",
             fontweight="bold", fontsize=11, y=1.03)
plt.tight_layout()
ns.save_fig(fig, f"results/06_singlecell_deep{tag}")
print("done 06b")
