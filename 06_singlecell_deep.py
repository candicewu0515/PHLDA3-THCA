#!/usr/bin/env python3
"""06 — Deeper single-cell analysis of PHLDA3 in thyroid tumours (GSE193581).

Three questions:
  (A) Are PHLDA3-high cells malignant?  -> lineage-marker confirmation
      (inferCNV proper was attempted but is NOT runnable in this offline sandbox:
       it requires per-gene genomic coordinates from a GTF/BioMart download and an
       R/Bioconductor stack, neither reachable here. We instead confirm malignant
       identity from author cell-type labels + epithelial/thyroid lineage markers
       vs immune/stromal markers in PHLDA3+ cells.)
  (B) p53 target-gene program vs PHLDA3 (PHLDA3 is a direct p53 target).
      Full pySCENIC needs cisTarget feather databases (GBs, unreachable offline),
      so we score a curated canonical p53 target-gene program per cell
      (AUCell/score_genes-style mean of normalised targets, PHLDA3 excluded to
      avoid circularity) and relate it to PHLDA3.
  (C) Ligand-receptor communication of PHLDA3-high malignant cells (CellChat-style)
      using liana's bundled `consensus` resource (offline) — no R / CellChat needed.

Memory-light: only a gene PANEL (L-R genes + p53 targets + markers + PHLDA3) is kept.
Normalisation (CP10k+log1p) uses the FULL per-cell UMI total before subsetting.

Run:  python3 06_singlecell_deep.py                 # all tumours (PTC+ATC)
      python3 06_singlecell_deep.py --group PTC      # single tumour type
Outputs (results/): 06_sc_panel.parquet (cache), 06_p53_program.csv,
      06_malignant_markers.csv, 06_liana_phlda3hi.csv, 06_singlecell_deep.{svg,pdf,tiff,png}
"""
import argparse, glob, os, re
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr
import nature_style as ns
ns.set_style(font_size=8)

GENE = "PHLDA3"
P53_TARGETS = ["CDKN1A","MDM2","BAX","BBC3","PMAIP1","GADD45A","SFN","SESN1","SESN2",
    "RRM2B","DDB2","XPC","POLH","TP53I3","TP53INP1","FAS","TNFRSF10B","ZMAT3","BTG2",
    "FDXR","APAF1","PERP","RPS27L","CCNG1","GDF15","AEN","TRIAP1","EDA2R","TNFRSF10A",
    "EI24","PLK3","CD82"]                       # PHLDA3 intentionally EXCLUDED
EPI_MARKERS = ["EPCAM","KRT18","KRT19","KRT8","TG","TPO","TSHR","PAX8","CDH1"]
OTHER_MARKERS = {"Immune":["PTPRC","CD3D","LYZ","CD68"], "Fibroblast":["COL1A1","DCN","LUM"],
                 "Endothelial":["PECAM1","VWF","CDH5"]}


CACHEDIR = "results/06_cache"

def get_panel():
    """Build (once) the gene panel = L-R genes (liana consensus) + p53 targets + markers."""
    pf = "results/06_panel.txt"
    if os.path.exists(pf):
        return open(pf).read().split()
    import liana as li
    res = li.resource.select_resource("consensus")
    lr_genes = set()
    for col in ("ligand", "receptor"):
        for v in res[col].astype(str):
            lr_genes.update(v.split("_"))                 # split complexes
    markers = EPI_MARKERS + sum(OTHER_MARKERS.values(), [])
    panel = sorted(lr_genes | set(P53_TARGETS) | set(markers) | {GENE})
    open(pf, "w").write("\n".join(panel))
    return panel


def build_panel(group):
    """Resumable per-sample caching (each UMI file processed once, written to CACHEDIR)."""
    os.makedirs(CACHEDIR, exist_ok=True)
    panel = get_panel(); pset = set(panel)
    anno = pd.read_csv("gse193581_anno.txt.gz", sep="\t", header=None, skiprows=1,
                       names=["barcode","sample","celltype"])
    bc2ct = dict(zip(anno.barcode, anno.celltype))
    files = []
    for f in sorted(glob.glob("gse193581_raw/*_UMI.txt.gz")):
        sample = re.search(r"_([A-Z]+\d+)_UMI", os.path.basename(f)).group(1)
        grp = "PTC" if sample.startswith("PTC") else "ATC" if sample.startswith("ATC") else "NORM"
        if group != "ALL" and grp != group and not (group == "TUMOR" and grp in ("PTC","ATC")):
            continue
        files.append((f, sample, grp))
    for f, sample, grp in files:
        out = f"{CACHEDIR}/{sample}.parquet"
        if os.path.exists(out):
            print(f"  {sample:8s} cached -> skip", flush=True); continue
        total = None; rows = []
        for ch in pd.read_csv(f, sep="\t", index_col=0, chunksize=3000):
            s = ch.sum(axis=0)
            total = s if total is None else total.add(s, fill_value=0)
            hit = ch.index.intersection(pset)
            if len(hit): rows.append(ch.loc[hit])
        sub = pd.concat(rows) if rows else pd.DataFrame(columns=total.index)
        keep_cells = [c for c in total.index if c in bc2ct]
        if not keep_cells:
            pd.DataFrame().to_parquet(out); continue
        tot = total[keep_cells].replace(0, np.nan)
        sub = sub.reindex(columns=keep_cells)
        norm = np.log1p(sub.div(tot, axis=1) * 1e4).T
        norm["celltype"] = [bc2ct[c] for c in keep_cells]
        norm["sample"] = sample; norm["grp"] = grp
        norm.reset_index().rename(columns={"index":"barcode"}).to_parquet(out)
        print(f"  {sample:8s} {grp} cells={len(keep_cells)} panel_genes={sub.shape[0]} -> cached", flush=True)
    # combine
    parts = [pd.read_parquet(f"{CACHEDIR}/{s}.parquet") for _, s, _ in files]
    parts = [p for p in parts if len(p)]
    d = pd.concat(parts, ignore_index=True).fillna(0.0)
    return d, panel


def run(group):
    cache = f"results/06_sc_panel_{group}.parquet"
    if os.path.exists(cache):
        d = pd.read_parquet(cache); print("loaded cache", cache, d.shape)
    else:
        d, _ = build_panel(group)
        d.to_parquet(cache); print("cells x cols:", d.shape, "-> cached")
    gcols = [c for c in d.columns if c not in ("barcode","celltype","sample","grp")]

    # ---------- (B) p53 target-gene program ----------
    p53 = [g for g in P53_TARGETS if g in gcols]
    d["p53_program"] = d[p53].mean(axis=1)
    mal = d[d.celltype == "Malignant cell"].copy()
    mal["PHLDA3pos"] = mal[GENE] > 0
    rho, prho = spearmanr(mal[GENE], mal["p53_program"])
    pos = mal[mal.PHLDA3pos].p53_program; neg = mal[~mal.PHLDA3pos].p53_program
    u, pmw = mannwhitneyu(pos, neg, alternative="greater")
    pd.DataFrame([{"n_malignant":len(mal),"spearman_rho":rho,"spearman_p":prho,
                   "p53prog_PHLDA3pos_mean":pos.mean(),"p53prog_PHLDA3neg_mean":neg.mean(),
                   "MWU_p_pos_gt_neg":pmw,"n_p53_genes":len(p53)}]
                 ).to_csv(f"results/06_p53_program_{group}.csv", index=False)
    print(f"  p53 program vs PHLDA3 (malignant): rho={rho:+.3f} p={prho:.1e} | "
          f"PHLDA3+ {pos.mean():.3f} vs PHLDA3- {neg.mean():.3f} MWU p={pmw:.1e}")

    # ---------- (A) malignant lineage confirmation ----------
    # PHLDA3+ cells: composition across cell types + marker profile
    d["PHLDA3pos"] = d[GENE] > 0
    comp = (d[d.PHLDA3pos].celltype.value_counts(normalize=True) * 100).round(1)
    epi = [g for g in EPI_MARKERS if g in gcols]
    mk_rows = []
    malpos = d[(d.celltype == "Malignant cell") & d.PHLDA3pos]
    immune = d[d.celltype.isin(["T cell","B cell","Myeloid cell","NK cell"])]
    for label, genes in [("Epithelial/thyroid", epi)] + list(OTHER_MARKERS.items()):
        gs = [g for g in genes if g in gcols]
        mk_rows.append({"marker_set":label,
                        "PHLDA3hi_malignant_mean":malpos[gs].mean(axis=1).mean(),
                        "immune_cells_mean":immune[gs].mean(axis=1).mean()})
    mk = pd.DataFrame(mk_rows)
    mk.to_csv(f"results/06_malignant_markers_{group}.csv", index=False)
    comp.to_csv(f"results/06_PHLDA3pos_composition_{group}.csv")
    print("  PHLDA3+ cell composition (%):", comp.head(4).to_dict())

    # ---------- (C) liana L-R from PHLDA3-high malignant ----------
    liana_top = run_liana(d, gcols, group)

    make_figure(d, mal, mk, comp, liana_top, rho, prho, pos, neg, pmw, group)
    return dict(rho=rho, pmw=pmw, comp=comp, mk=mk, liana=liana_top)


def run_liana(d, gcols, group):
    try:
        import anndata as ad, scanpy as sc, liana as li
        rng = np.random.default_rng(0)
        # subsample per celltype for speed/memory
        idx = []
        for ct, sub in d.groupby("celltype"):
            n = min(len(sub), 600); idx += list(rng.choice(sub.index, n, replace=False))
        sm = d.loc[idx].copy()
        # split malignant by PHLDA3 status
        grp = sm.celltype.astype(str).copy()
        malmask = sm.celltype == "Malignant cell"
        med = d.loc[d.celltype == "Malignant cell", GENE]
        thr = med[med > 0].median() if (med > 0).any() else 0
        grp[malmask & (sm[GENE] > thr)] = "Malignant_PHLDA3hi"
        grp[malmask & (sm[GENE] <= thr)] = "Malignant_PHLDA3lo"
        A = ad.AnnData(X=sm[gcols].values.astype("float32"),
                       obs=pd.DataFrame({"label": grp.values}, index=sm.barcode.values),
                       var=pd.DataFrame(index=gcols))
        A.layers["normalized"] = A.X.copy()                 # already log1p CP10k
        li.mt.rank_aggregate(A, groupby="label", expr_prop=0.1, use_raw=False,
                             layer="normalized", verbose=False, n_perms=50, seed=0)
        r = A.uns["liana_res"]
        r = r[r.source == "Malignant_PHLDA3hi"].copy()
        # rank by interaction magnitude (lower magnitude_rank = stronger)
        mcol = "magnitude_rank" if "magnitude_rank" in r.columns else "lr_means"
        r = r.sort_values(mcol).head(40)
        r.to_csv(f"results/06_liana_phlda3hi_{group}.csv", index=False)
        print(f"  liana: {len(r)} top interactions from Malignant_PHLDA3hi saved")
        return r
    except Exception as e:
        print("  ! liana step failed:", type(e).__name__, str(e)[:160])
        return None


def make_figure(d, mal, mk, comp, liana_top, rho, prho, pos, neg, pmw, group):
    ncol = 4 if liana_top is not None and len(liana_top) else 3
    fig, ax = plt.subplots(1, ncol, figsize=(3.5*ncol, 4.2))

    # (a) malignant confirmation: marker means PHLDA3hi-malignant vs immune
    y = np.arange(len(mk))[::-1]; h = 0.36
    ax[0].barh(y+h/2, mk.PHLDA3hi_malignant_mean, height=h, color=ns.C["tumor"], alpha=0.9, label="PHLDA3-hi malignant")
    ax[0].barh(y-h/2, mk.immune_cells_mean, height=h, color=ns.C["normal"], alpha=0.9, label="Immune cells")
    ax[0].set_yticks(y); ax[0].set_yticklabels(mk.marker_set, fontsize=8)
    ax[0].set_xlabel("Mean expr (log CP10k)"); ax[0].legend(fontsize=6.5, loc="lower right")
    ax[0].set_title("Lineage identity of\nPHLDA3-high cells", fontweight="bold", fontsize=9.5)
    ax[0].spines[["top","right"]].set_visible(False); ns.panel_label(ax[0], "a")

    # (b) p53 program: PHLDA3+ vs PHLDA3- malignant
    bp = ax[1].boxplot([neg.values, pos.values], widths=0.55, patch_artist=True, showfliers=False,
                       medianprops=dict(color="black", lw=1.3))
    for patch, c in zip(bp["boxes"], [ns.C["normal"], ns.C["tumor"]]):
        patch.set_facecolor(c); patch.set_alpha(0.65)
    ax[1].set_xticks([1,2]); ax[1].set_xticklabels([f"PHLDA3−\n(n={len(neg)})", f"PHLDA3+\n(n={len(pos)})"])
    ax[1].set_ylabel("p53 target-gene program")
    ax[1].set_title(f"p53 program in malignant cells\nMWU p={pmw:.1e}", fontweight="bold", fontsize=9.5)
    ax[1].spines[["top","right"]].set_visible(False); ns.panel_label(ax[1], "b")

    # (c) p53 program vs PHLDA3 scatter (subsample for clarity)
    s = mal.sample(min(len(mal), 3000), random_state=0)
    ax[2].scatter(s[GENE], s.p53_program, s=4, color=ns.C["tumor"], alpha=0.2, edgecolors="none")
    m, b = np.polyfit(mal[GENE], mal.p53_program, 1)
    xs = np.linspace(mal[GENE].min(), mal[GENE].max(), 50)
    ax[2].plot(xs, m*xs+b, color=ns.C["grey_dark"], lw=1.4)
    ax[2].set_xlabel(f"{GENE} (log CP10k)"); ax[2].set_ylabel("p53 target-gene program")
    ax[2].set_title(f"p53 program vs {GENE}\nρ={rho:+.2f} p={prho:.1e}", fontweight="bold", fontsize=9.5)
    ax[2].spines[["top","right"]].set_visible(False); ns.panel_label(ax[2], "c")

    # (d) liana top interactions
    if ncol == 4:
        r = liana_top.copy()
        r["pair"] = r["ligand_complex"] + "→" + r["receptor_complex"] + "\n(to " + r["target"].astype(str) + ")"
        mcol = "magnitude_rank" if "magnitude_rank" in r.columns else "lr_means"
        r = r.sort_values(mcol).head(10).iloc[::-1]
        score = -np.log10(r[mcol].clip(lower=1e-6)) if mcol == "magnitude_rank" else r[mcol]
        ax[3].barh(range(len(r)), score, color=ns.C["accent2"], alpha=0.9)
        ax[3].set_yticks(range(len(r))); ax[3].set_yticklabels(r["pair"], fontsize=6)
        ax[3].set_xlabel("−log10 magnitude_rank" if mcol=="magnitude_rank" else "LR mean")
        ax[3].set_title("Top signalling from\nPHLDA3-high malignant cells", fontweight="bold", fontsize=9.5)
        ax[3].spines[["top","right"]].set_visible(False); ns.panel_label(ax[3], "d")

    fig.suptitle(f"PHLDA3 in malignant thyroid cells — p53 program & cell communication "
                 f"(GSE193581, {group})", fontweight="bold", fontsize=11, y=1.03)
    plt.tight_layout()
    tag = "" if group == "TUMOR" else f"_{group}"
    ns.save_fig(fig, f"results/06_singlecell_deep{tag}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--group", default="TUMOR", help="TUMOR (PTC+ATC, default) | PTC | ATC | ALL")
    args = ap.parse_args()
    run(args.group)
    print("\n06 done.")
