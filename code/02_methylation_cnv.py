#!/usr/bin/env python3
"""02 — Regulatory mechanism of PHLDA3 over-expression in TCGA-THCA:
   promoter DNA-methylation (HM450) vs expression, and copy-number (GISTIC) vs expression.

Data are pulled live from the cBioPortal REST API for the single gene PHLDA3
(entrezGeneId 23612) in study `thca_tcga_pan_can_atlas_2018`:
   * mRNA expression z-scores (RNA-seq v2, _rna_seq_v2_mrna_median_Zscores)
   * GISTIC copy-number calls       (_gistic, -2..+2)
   * DNA methylation (HM450)        (_methylation_hm450, beta values)
Then: Spearman correlation of expression vs methylation (expect negative if promoter
hypomethylation drives induction) and vs copy-number (expect positive if amplification
contributes), with boxplots by GISTIC category.

IMPORTANT — execution status in THIS repo:
   This script was NOT executed in the offline analysis sandbox used to add analyses
   01/06/07 (that sandbox allows only PyPI traffic, so cBioPortal/GDC/Xena are blocked).
   It is written to run unchanged on a networked machine. If no network is available
   it prints a clear message and exits WITHOUT producing fabricated results.

Run:  python3 02_methylation_cnv.py                 # THCA
      python3 02_methylation_cnv.py --study <cbioportal_study_id>
Outputs (results/): 02_phlda3_meth_cnv.csv, 02_methylation_cnv.{svg,pdf,tiff,png}
"""
import argparse, sys
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, kruskal
import nature_style as ns
ns.set_style(font_size=8)

GENE, ENTREZ = "PHLDA3", 23612
API = "https://www.cbioportal.org/api"


def fetch_profile(study, profile_suffix, sample_list_suffix):
    import requests
    prof = f"{study}_{profile_suffix}"
    slist = f"{study}_{sample_list_suffix}"
    url = f"{API}/molecular-profiles/{prof}/molecular-data"
    r = requests.get(url, params={"sampleListId": slist, "entrezGeneId": ENTREZ,
                                  "projection": "SUMMARY"}, timeout=60)
    r.raise_for_status()
    d = pd.DataFrame(r.json())
    if d.empty: return pd.Series(dtype=float)
    return pd.Series(d["value"].values, index=d["sampleId"].values, dtype=float)


def main(study):
    try:
        expr = fetch_profile(study, "rna_seq_v2_mrna_median_Zscores", "rna_seq_v2_mrna")
        cnv = fetch_profile(study, "gistic", "cna")
        meth = fetch_profile(study, "methylation_hm450", "methylation_hm450")
    except Exception as e:
        print("NETWORK/API UNAVAILABLE — cannot fetch cBioPortal data in this environment.")
        print("Reason:", type(e).__name__, str(e)[:140])
        print("Run this script on a networked machine to produce results (no fabrication).")
        sys.exit(0)

    if expr.empty:
        print("No expression returned; check study id."); sys.exit(0)
    df = pd.DataFrame({"expr": expr, "cnv": cnv, "meth": meth})
    df.to_csv("results/02_phlda3_meth_cnv.csv")
    print("samples with expr/cnv/meth:",
          df.expr.notna().sum(), df.cnv.notna().sum(), df.meth.notna().sum())

    fig, ax = plt.subplots(1, 2, figsize=(9.5, 4.3))
    # methylation vs expression
    m = df[["expr","meth"]].dropna()
    rho, p = spearmanr(m.meth, m.expr) if len(m) > 10 else (np.nan, np.nan)
    ax[0].scatter(m.meth, m.expr, s=10, color=ns.C["tumor"], alpha=0.5, edgecolors="none")
    if len(m) > 10:
        a, b = np.polyfit(m.meth, m.expr, 1); xs = np.linspace(m.meth.min(), m.meth.max(), 50)
        ax[0].plot(xs, a*xs+b, color=ns.C["grey_dark"], lw=1.4)
    ax[0].set_xlabel("PHLDA3 promoter methylation (HM450 beta)")
    ax[0].set_ylabel("PHLDA3 mRNA (z-score)")
    ax[0].set_title(f"Methylation vs expression\nrho={rho:+.2f} p={p:.1e}", fontweight="bold", fontsize=10)
    ax[0].spines[["top","right"]].set_visible(False); ns.panel_label(ax[0], "a")

    # CNV vs expression
    c = df[["expr","cnv"]].dropna()
    cats = sorted(c.cnv.unique())
    labels = {-2:"DeepDel",-1:"ShallowDel",0:"Diploid",1:"Gain",2:"Amp"}
    groups = [c[c.cnv == k].expr.values for k in cats]
    bp = ax[1].boxplot(groups, widths=0.6, patch_artist=True, showfliers=False,
                       medianprops=dict(color="black", lw=1.2))
    for patch in bp["boxes"]: patch.set_facecolor(ns.C["tumor"]); patch.set_alpha(0.6)
    ax[1].set_xticks(range(1, len(cats)+1)); ax[1].set_xticklabels([labels.get(k, k) for k in cats])
    rho2, p2 = spearmanr(c.cnv, c.expr) if len(c) > 10 else (np.nan, np.nan)
    ax[1].set_ylabel("PHLDA3 mRNA (z-score)")
    ax[1].set_title(f"Copy-number vs expression\nrho={rho2:+.2f} p={p2:.1e}", fontweight="bold", fontsize=10)
    ax[1].spines[["top","right"]].set_visible(False); ns.panel_label(ax[1], "b")

    pd.DataFrame([{"meth_expr_rho": rho, "meth_expr_p": p, "cnv_expr_rho": rho2, "cnv_expr_p": p2,
                   "n_meth": len(m), "n_cnv": len(c)}]).to_csv("results/02_summary.csv", index=False)
    fig.suptitle(f"PHLDA3 regulatory mechanism in {study}", fontweight="bold", fontsize=11, y=1.02)
    plt.tight_layout(); ns.save_fig(fig, "results/02_methylation_cnv")
    print("02 done.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--study", default="thca_tcga_pan_can_atlas_2018")
    main(ap.parse_args().study)
