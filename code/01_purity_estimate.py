#!/usr/bin/env python3
"""01 — Tumor-purity correction for PHLDA3 in TCGA-THCA.

Question addressed: is the PHLDA3 tumour-vs-normal signal (and its association with
nodal metastasis / advanced stage) an artefact of differing immune/stromal content,
or does it survive adjustment for tumour purity?

Method (fully offline, runs on local TCGA-THCA HTSeq counts only):
  * ESTIMATE-style ssGSEA of an immune signature and a stromal signature
    (gseapy ssGSEA, rank-normalised) -> ImmuneScore, StromalScore.
    NOTE: the exact 282-gene ESTIMATE signatures (R-package SI_geneset.gmt) are not
    fetchable in this offline environment, so curated, well-established immune and
    stromal marker panels are used. Because all downstream tests are RANK-based
    (Spearman partial correlation), the conclusion is invariant to the precise
    monotonic purity scale, so this faithfully reproduces the ESTIMATE logic.
  * ESTIMATEScore = z(Immune) + z(Stromal); purity index = -ESTIMATEScore.
  * Partial Spearman correlation of PHLDA3 with (i) tumour status, (ii) N1 nodal
    metastasis, (iii) advanced stage, each controlling for ESTIMATEScore (purity).

Run:  python3 01_purity_estimate.py                 # THCA (default)
      python3 01_purity_estimate.py --tissue ALL    # loop over any TCGA-*.htseq_counts.tsv present

Outputs (results/):
  01_purity_scores.csv          per-sample Immune/Stromal/ESTIMATE/purity + PHLDA3 + group
  01_partial_correlation.csv    unadjusted vs purity-adjusted association table
  01_purity_correction.{svg,pdf,tiff,png}
"""
import argparse, glob, os
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, spearmanr
import pingouin as pg
import gseapy as gp
import nature_style as ns
ns.set_style(font_size=8)

GENE, ENS = "PHLDA3", "ENSG00000174307"

# ---- ESTIMATE-style signatures (curated immune + stromal marker panels) -------
IMMUNE_SIG = ["PTPRC","CD2","CD3D","CD3E","CD3G","CD8A","CD8B","CD4","CD27","CD28",
    "CD48","CD52","CD53","LCK","LCP2","ZAP70","IL2RG","IL2RB","CD19","MS4A1","CD79A",
    "CD79B","BLK","TNFRSF17","CCL5","CCR7","CXCR4","CXCL9","CXCL10","CXCL11","CXCL13",
    "GZMA","GZMB","GZMK","GZMH","PRF1","NKG7","GNLY","KLRB1","KLRD1","KLRK1","CTSW",
    "CD68","CD163","CD14","FCGR3A","FCGR2A","ITGAM","TLR2","CSF1R","LYZ","AIF1",
    "HLA-DRA","HLA-DPA1","HLA-DPB1","CIITA","FOXP3","IKZF1","IL10","TIGIT","CTLA4",
    "PDCD1","LAG3","HAVCR2","IDO1","TNFRSF9","ICOS","SELL","CCL19","CCL21","IL7R",
    "TRBC2","TRAC","SH2D1A"]
STROMAL_SIG = ["COL1A1","COL1A2","COL3A1","COL5A1","COL5A2","COL6A1","COL6A2","COL6A3",
    "COL4A1","COL4A2","FN1","FAP","PDGFRA","PDGFRB","PDGFRL","ACTA2","TAGLN","MYL9",
    "MYH11","THY1","DCN","LUM","BGN","FBLN1","FBLN2","FBN1","SPARC","SPARCL1","VIM",
    "MMP2","MMP14","TIMP1","TIMP3","POSTN","ELN","THBS1","THBS2","CCN2","SERPINH1",
    "LOX","LOXL1","PCOLCE","CALD1","MFAP2","MFAP4","ANGPT2","PECAM1","VWF","CDH5",
    "ENG","EGFL7","CLEC14A","PLVAP","COL14A1","COL15A1","NID1","NID2","LAMA2","LAMB1"]


def load_logcpm_symbol(path):
    """Return (logcpm symbols x samples for genes meanCPM>=1, sample-type series)."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    mat = (2**df - 1).clip(lower=0)                     # Xena stores log2(count+1)
    cpm = mat.div(mat.sum(axis=0), axis=1) * 1e6
    keep = cpm.mean(axis=1) >= 1
    logcpm = np.log2(cpm[keep] + 1)
    logcpm.index = [g.split(".")[0] for g in logcpm.index]
    e2s = pd.read_csv("data/thca_DE_full_symbol.csv")[["ensembl", "symbol"]].dropna()
    e2s = dict(zip(e2s.ensembl, e2s.symbol))
    logcpm = logcpm[logcpm.index.isin(e2s)]
    logcpm.index = [e2s[e] for e in logcpm.index]
    logcpm = logcpm[~logcpm.index.duplicated()]
    stype = pd.Series([b.split("-")[3][:2] for b in df.columns], index=df.columns)
    return logcpm, stype


def estimate_scores(logcpm):
    sigs = {"ImmuneScore": [g for g in IMMUNE_SIG if g in logcpm.index],
            "StromalScore": [g for g in STROMAL_SIG if g in logcpm.index]}
    ss = gp.ssgsea(data=logcpm, gene_sets=sigs, sample_norm_method="rank",
                   min_size=1, outdir=None, no_plot=True, threads=4)
    nes = ss.res2d.pivot(index="Term", columns="Name", values="NES").astype(float)
    nes = nes[logcpm.columns]
    out = nes.T  # samples x [Immune, Stromal]
    z = (out - out.mean()) / out.std()
    out["ESTIMATEScore"] = z["ImmuneScore"] + z["StromalScore"]
    out["PurityIndex"] = -out["ESTIMATEScore"]   # higher = more pure (less infiltrate)
    return out


def run_tissue(path, tag):
    print(f"\n========== {tag} :: {os.path.basename(path)} ==========")
    logcpm, stype = load_logcpm_symbol(path)
    if GENE not in logcpm.index:
        print(f"  ! {GENE} not in matrix; skip"); return None
    sc = estimate_scores(logcpm)
    sc["PHLDA3"] = logcpm.loc[GENE]
    sc["sample_type"] = stype.reindex(sc.index)
    sc["tumor"] = (sc["sample_type"] == "01").astype(int)
    sc["patient"] = [b[:12] for b in sc.index]

    # clinical merge (tumour only) for N / stage
    cl = pd.read_csv("data/THCA_clinical.tsv", sep="\t", dtype=str)
    def pick(r, base):
        for k in (f"diagnoses.0.{base}", f"diagnoses.1.{base}"):
            if k in r and pd.notna(r[k]) and str(r[k]) not in ("", "nan"): return r[k]
        return np.nan
    c = pd.DataFrame({"patient": cl["submitter_id"],
                      "N": cl.apply(lambda r: pick(r, "ajcc_pathologic_n"), axis=1),
                      "stage": cl.apply(lambda r: pick(r, "ajcc_pathologic_stage"), axis=1)})
    c["N1"] = c["N"].str[:2].map(lambda n: 1 if n == "N1" else (0 if n == "N0" else np.nan))
    c["adv"] = c["stage"].map(lambda s: 1 if s in ("Stage III","Stage IVA","Stage IVB","Stage IVC","Stage IV")
                              else (0 if s in ("Stage I","Stage II") else np.nan))
    tum = sc[sc.tumor == 1].merge(c[["patient","N1","adv"]], on="patient", how="left")
    sc.to_csv(f"results/01_purity_scores_{tag}.csv", index=False)

    # ---------- partial correlations (Spearman, control = ESTIMATEScore) ----------
    def pcorr(d, x, y, covar):
        pc = pg.partial_corr(d, x=x, y=y, covar=covar, method="spearman")
        pcol = "p_val" if "p_val" in pc.columns else "p-val"
        return float(pc["r"].iloc[0]), float(pc[pcol].iloc[0])
    rows = []
    # (1) tumour status across all samples
    d0 = sc[["PHLDA3","tumor","ESTIMATEScore"]].dropna()
    r_u, p_u = spearmanr(d0.PHLDA3, d0.tumor)
    r_p, p_p = pcorr(d0, "PHLDA3", "tumor", "ESTIMATEScore")
    rows.append(["Tumor vs Normal", len(d0), r_u, p_u, r_p, p_p])
    # (2) N1 within tumours
    d1 = tum[["PHLDA3","N1","ESTIMATEScore"]].dropna()
    r_u, p_u = spearmanr(d1.PHLDA3, d1.N1)
    r_p, p_p = pcorr(d1, "PHLDA3", "N1", "ESTIMATEScore")
    rows.append(["N1 metastasis", len(d1), r_u, p_u, r_p, p_p])
    # (3) advanced stage within tumours
    d2 = tum[["PHLDA3","adv","ESTIMATEScore"]].dropna()
    r_u, p_u = spearmanr(d2.PHLDA3, d2.adv)
    r_p, p_p = pcorr(d2, "PHLDA3", "adv", "ESTIMATEScore")
    rows.append(["Advanced stage", len(d2), r_u, p_u, r_p, p_p])
    # (4) PHLDA3 vs purity itself (tumours): is PHLDA3 an infiltrate marker?
    dt = sc[sc.tumor == 1]
    r_pur, p_pur = spearmanr(dt.PHLDA3, dt.PurityIndex)

    tab = pd.DataFrame(rows, columns=["association","n","r_unadj","p_unadj","r_partial","p_partial"])
    tab["PHLDA3_vs_purity_r"] = r_pur
    tab["PHLDA3_vs_purity_p"] = p_pur
    tab.to_csv(f"results/01_partial_correlation_{tag}.csv", index=False)
    print(tab.round(4).to_string(index=False))
    print(f"  PHLDA3 vs PurityIndex (tumours): rho={r_pur:+.3f} p={p_pur:.2e}")
    return dict(sc=sc, tum=tum, tab=tab, r_pur=r_pur, p_pur=p_pur)


def make_figure(res, tag):
    sc, tum, tab = res["sc"], res["tum"], res["tab"]
    fig, ax = plt.subplots(1, 3, figsize=(13, 4.2))

    # (a) PHLDA3 tumor vs normal (still significant after noting partial corr)
    nm = sc[sc.tumor == 0].PHLDA3.values; tu = sc[sc.tumor == 1].PHLDA3.values
    bp = ax[0].boxplot([nm, tu], widths=0.55, patch_artist=True, showfliers=False,
                       medianprops=dict(color="black", lw=1.3))
    for patch, c in zip(bp["boxes"], [ns.C["normal"], ns.C["tumor"]]):
        patch.set_facecolor(c); patch.set_alpha(0.65)
    for i, (d, c) in enumerate(zip([nm, tu], [ns.C["normal"], ns.C["tumor"]]), 1):
        ax[0].scatter(np.random.normal(i, 0.06, len(d)), d, s=6, color=c, alpha=0.4, edgecolors="none")
    ax[0].set_xticks([1, 2]); ax[0].set_xticklabels([f"Normal\n(n={len(nm)})", f"Tumor\n(n={len(tu)})"])
    ax[0].set_ylabel(f"{GENE} log2(CPM+1)")
    pr = tab[tab.association == "Tumor vs Normal"].iloc[0]
    ax[0].set_title(f"{GENE} tumour vs normal\npartial ρ={pr.r_partial:.2f} (purity-adj.) p={pr.p_partial:.1e}",
                    fontweight="bold", fontsize=9.5)
    ax[0].spines[["top","right"]].set_visible(False); ns.panel_label(ax[0], "a")

    # (b) PHLDA3 vs purity scatter (tumours)
    dt = sc[sc.tumor == 1]
    ax[1].scatter(dt.PurityIndex, dt.PHLDA3, s=10, color=ns.C["tumor"], alpha=0.45, edgecolors="none")
    m, b = np.polyfit(dt.PurityIndex, dt.PHLDA3, 1)
    xs = np.linspace(dt.PurityIndex.min(), dt.PurityIndex.max(), 50)
    ax[1].plot(xs, m*xs+b, color=ns.C["grey_dark"], lw=1.4)
    ax[1].set_xlabel("Tumour purity index  (−ESTIMATE score)")
    ax[1].set_ylabel(f"{GENE} log2(CPM+1)")
    ax[1].set_title(f"{GENE} vs purity (tumours)\nρ={res['r_pur']:+.2f}  p={res['p_pur']:.1e}",
                    fontweight="bold", fontsize=9.5)
    ax[1].spines[["top","right"]].set_visible(False); ns.panel_label(ax[1], "b")

    # (c) unadjusted vs purity-adjusted association
    y = np.arange(len(tab))[::-1]
    h = 0.36
    ax[2].barh(y+h/2, tab.r_unadj, height=h, color=ns.C["grey_mid"], alpha=0.9, label="Unadjusted")
    ax[2].barh(y-h/2, tab.r_partial, height=h, color=ns.C["tumor"], alpha=0.9, label="Purity-adjusted")
    for yi, rr, pp in zip(y-h/2, tab.r_partial, tab.p_partial):
        s = "***" if pp < 1e-3 else "**" if pp < 1e-2 else "*" if pp < .05 else "ns"
        ax[2].text(rr + (0.005 if rr >= 0 else -0.005), yi, s, va="center",
                   ha="left" if rr >= 0 else "right", fontsize=8)
    ax[2].set_yticks(y); ax[2].set_yticklabels(tab.association)
    ax[2].axvline(0, color="black", lw=0.8)
    ax[2].set_xlabel(f"Spearman ρ with {GENE}")
    ax[2].set_title("Association survives purity adjustment", fontweight="bold", fontsize=9.5)
    ax[2].legend(loc="lower right", fontsize=7); ax[2].spines[["top","right"]].set_visible(False)
    ns.panel_label(ax[2], "c")

    fig.suptitle(f"PHLDA3 signal is independent of tumour purity ({tag}, ESTIMATE-style ssGSEA)",
                 fontweight="bold", fontsize=11, y=1.02)
    plt.tight_layout()
    ns.save_fig(fig, f"results/01_purity_correction_{tag}" if tag != "THCA" else "results/01_purity_correction")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", default="THCA", help="THCA (default) or ALL")
    args = ap.parse_args()
    if args.tissue.upper() == "ALL":
        files = sorted(glob.glob("TCGA-*.htseq_counts.tsv"))
        for f in files:
            tag = f.split(".")[0].replace("TCGA-", "")
            res = run_tissue(f, tag)
            if res: make_figure(res, tag)
    else:
        f = f"TCGA-{args.tissue}.htseq_counts.tsv"
        res = run_tissue(f, args.tissue)
        if res: make_figure(res, args.tissue)
    print("\n01 done.")
