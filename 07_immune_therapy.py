#!/usr/bin/env python3
"""07 — PHLDA3, immune checkpoints and therapy relevance in TCGA-THCA.

Runnable offline (local expression only):
  * PHLDA3 vs an expanded immune-checkpoint / immunomodulator panel:
    Pearson + Spearman, PLUS purity-adjusted partial Spearman (covariate =
    ESTIMATEScore from 01_purity_estimate.py) -> shows checkpoint links are not
    purely a consequence of immune infiltrate.
  * PHLDA3-high vs -low checkpoint expression (Mann-Whitney) for the top hits.

NOT runnable in this offline sandbox (reported, not fabricated):
  * TMB: requires per-sample somatic mutation counts (GDC/cBioPortal) - network blocked.
  * Drug-sensitivity prediction (oncoPredict GDSC/CTRP) / CMap: require R + the GDSC/CTRP
    training matrices / CMap API - not installable/reachable here.

Run:  python3 07_immune_therapy.py               # THCA
      python3 07_immune_therapy.py --tissue ALL
Outputs (results/): 07_checkpoint_correlation.csv, 07_immune_therapy.{svg,pdf,tiff,png}
"""
import argparse, glob, os
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr, mannwhitneyu
import pingouin as pg
import nature_style as ns
ns.set_style(font_size=8)

GENE = "PHLDA3"
CHECKPOINTS = ["CD274","PDCD1","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT","BTLA",
    "IDO1","CD276","VTCN1","TNFRSF9","ICOS","CD27","LGALS9","VSIR","CD80","CD86",
    "CD40","CD70","TNFRSF18","TNFRSF4","CD200","NT5E"]


def load_tumor_logcpm_symbol(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    mat = (2**df - 1).clip(lower=0)
    stype = pd.Series([b.split("-")[3][:2] for b in df.columns], index=df.columns)
    tum = df.columns[stype == "01"]
    cpm = mat[tum].div(mat[tum].sum(axis=0), axis=1) * 1e6
    keep = cpm.mean(axis=1) >= 1
    logcpm = np.log2(cpm[keep] + 1)
    logcpm.index = [g.split(".")[0] for g in logcpm.index]
    e2s = pd.read_csv("thca_DE_full_symbol.csv")[["ensembl","symbol"]].dropna()
    e2s = dict(zip(e2s.ensembl, e2s.symbol))
    logcpm = logcpm[logcpm.index.isin(e2s)]; logcpm.index = [e2s[e] for e in logcpm.index]
    logcpm = logcpm[~logcpm.index.duplicated()]
    return logcpm


def run_tissue(path, tag):
    print(f"\n===== {tag} =====")
    logcpm = load_tumor_logcpm_symbol(path)
    if GENE not in logcpm.index:
        print("  ! no PHLDA3"); return None
    phl = logcpm.loc[GENE]
    # purity covariate from analysis 01 (reuse), matched by patient
    purity = None
    pf = f"results/01_purity_scores_{tag}.csv"
    if os.path.exists(pf):
        ps = pd.read_csv(pf)
        ps = ps[ps.tumor == 1].groupby("patient")["ESTIMATEScore"].mean()
        pat = pd.Series([c[:12] for c in phl.index], index=phl.index)
        purity = pat.map(ps)
        print(f"  purity covariate matched for {purity.notna().sum()}/{len(purity)} samples")

    rows = []
    for g in CHECKPOINTS:
        if g not in logcpm.index: continue
        x = logcpm.loc[g]
        r, p = pearsonr(x, phl); rs, ps_ = spearmanr(x, phl)
        rp, pp = np.nan, np.nan
        if purity is not None:
            dd = pd.DataFrame({"phl": phl.values, "cp": x.values, "pur": purity.values}).dropna()
            pc = pg.partial_corr(dd, x="phl", y="cp", covar="pur", method="spearman")
            pcol = "p_val" if "p_val" in pc.columns else "p-val"
            rp, pp = float(pc["r"].iloc[0]), float(pc[pcol].iloc[0])
        rows.append({"gene": g, "pearson_r": r, "pearson_p": p, "spearman_r": rs,
                     "spearman_p": ps_, "partial_r_purityadj": rp, "partial_p_purityadj": pp})
    cp = pd.DataFrame(rows).sort_values("spearman_r", ascending=False)
    cp.to_csv(f"results/07_checkpoint_correlation_{tag}.csv", index=False)
    print(cp[["gene","spearman_r","spearman_p","partial_r_purityadj","partial_p_purityadj"]]
          .head(12).round(3).to_string(index=False))

    # high/low PHLDA3 split for top-4 checkpoints
    med = phl.median(); hi = phl > med
    top4 = cp.head(4).gene.tolist()
    return dict(logcpm=logcpm, phl=phl, cp=cp, hi=hi, top4=top4, purity=purity, tag=tag)


def make_figure(res):
    logcpm, cp, hi, top4, tag = res["logcpm"], res["cp"], res["hi"], res["top4"], res["tag"]
    has_partial = res["purity"] is not None and cp["partial_r_purityadj"].notna().any()
    fig, ax = plt.subplots(1, 2, figsize=(12.5, 5.4),
                           gridspec_kw={"width_ratios": [1.25, 1]})

    # (a) correlation bars: spearman vs purity-adjusted
    d = cp.sort_values("spearman_r")
    y = np.arange(len(d))
    if has_partial:
        h = 0.4
        ax[0].barh(y+h/2, d.spearman_r, height=h, color=ns.C["grey_mid"], alpha=0.9, label="Spearman")
        ax[0].barh(y-h/2, d.partial_r_purityadj, height=h, color=ns.C["tumor"], alpha=0.9, label="Purity-adjusted")
        for yi, rr, pp in zip(y-h/2, d.partial_r_purityadj, d.partial_p_purityadj):
            s = "***" if pp < 1e-3 else "**" if pp < 1e-2 else "*" if pp < .05 else ""
            ax[0].text(rr + (0.006 if rr >= 0 else -0.006), yi, s, va="center",
                       ha="left" if rr >= 0 else "right", fontsize=7)
        ax[0].legend(fontsize=7, loc="lower right")
    else:
        ax[0].barh(y, d.spearman_r, color=[ns.C["tumor"] if v > 0 else ns.C["normal"] for v in d.spearman_r], alpha=0.9)
    ax[0].set_yticks(y); ax[0].set_yticklabels(d.gene, fontsize=8)
    ax[0].axvline(0, color="black", lw=0.8); ax[0].set_xlabel(f"Spearman rho with {GENE}")
    ax[0].set_title("PHLDA3 vs immune checkpoints\n(purity-adjusted)", fontweight="bold", fontsize=10.5)
    ax[0].spines[["top","right"]].set_visible(False); ns.panel_label(ax[0], "a")

    # (b) high/low boxplots for top-4 checkpoints
    pos = np.arange(len(top4))
    w = 0.36
    for i, g in enumerate(top4):
        lo = logcpm.loc[g][~hi].values; hh = logcpm.loc[g][hi].values
        for j, (vals, off, c) in enumerate([(lo, -w/2, ns.C["normal"]), (hh, w/2, ns.C["tumor"])]):
            bp = ax[1].boxplot([vals], positions=[i+off], widths=w*0.9, patch_artist=True,
                               showfliers=False, medianprops=dict(color="black", lw=1.1))
            bp["boxes"][0].set_facecolor(c); bp["boxes"][0].set_alpha(0.65)
        p = mannwhitneyu(lo, hh)[1]
        s = "***" if p < 1e-3 else "**" if p < 1e-2 else "*" if p < .05 else "ns"
        ax[1].text(i, max(hh.max(), lo.max())*1.0, s, ha="center", fontsize=9)
    ax[1].set_xticks(pos); ax[1].set_xticklabels(top4)
    ax[1].set_ylabel("Checkpoint expr log2(CPM+1)")
    from matplotlib.patches import Patch
    ax[1].legend(handles=[Patch(facecolor=ns.C["normal"], alpha=.65, label="PHLDA3 low"),
                          Patch(facecolor=ns.C["tumor"], alpha=.65, label="PHLDA3 high")], fontsize=7)
    ax[1].set_title("Checkpoints by PHLDA3 group", fontweight="bold", fontsize=10.5)
    ax[1].spines[["top","right"]].set_visible(False); ns.panel_label(ax[1], "b")

    fig.suptitle(f"PHLDA3 and immune-checkpoint landscape ({tag}, n={logcpm.shape[1]} tumours)",
                 fontweight="bold", fontsize=11, y=1.0)
    plt.tight_layout()
    tag2 = "" if tag == "THCA" else f"_{tag}"
    ns.save_fig(fig, f"results/07_immune_therapy{tag2}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--tissue", default="THCA")
    args = ap.parse_args()
    files = sorted(glob.glob("TCGA-*.htseq_counts.tsv")) if args.tissue.upper() == "ALL" \
            else [f"TCGA-{args.tissue}.htseq_counts.tsv"]
    for f in files:
        tag = f.split(".")[0].replace("TCGA-", "")
        res = run_tissue(f, tag)
        if res: make_figure(res)
    print("\nNOTE: TMB and oncoPredict/CMap drug prediction require external data/R and "
          "could not be run in this offline environment (see header).")
    print("07 done.")
