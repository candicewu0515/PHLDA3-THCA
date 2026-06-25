#!/usr/bin/env python3
"""Task 1 — tumour-purity correction of PHLDA3 immune associations.
ESTIMATE (R) -> tumour purity; partial Spearman correlation of PHLDA3 with immune
checkpoints and Danaher immune-cell ssGSEA scores, controlling for purity.
Shows which immune associations survive purity adjustment."""
import pandas as pd, numpy as np, subprocess, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import pingouin as pg
import gseapy as gp
import nature_style as ns
ns.set_style(font_size=8)

GENE = "PHLDA3"
IMMUNE = {  # Danaher signatures (same as immune_thca.py)
 "B cells":["BLK","CD19","MS4A1","TNFRSF17","FCRL2","PNOC","SPIB","TCL1A"],
 "CD8 T cells":["CD8A","CD8B"],
 "Cytotoxic cells":["CTSW","GNLY","GZMA","GZMB","GZMH","KLRB1","KLRD1","KLRK1","NKG7","PRF1"],
 "DC":["CCL13","CD209","HSD11B1"],"Exhausted CD8":["CD244","EOMES","LAG3","PTGER4"],
 "Macrophages":["CD163","CD68","CD84","MS4A4A"],"Mast cells":["CPA3","HDC","MS4A2","TPSAB1","TPSB2"],
 "Neutrophils":["CSF3R","S100A12","CEACAM3","FCGR3B","FFAR2"],"NK CD56dim":["IL21R","KIR2DL3","KIR3DL1","KIR3DL2"],
 "NK cells":["NCR1","XCL1","XCL2"],"T cells":["CD3D","CD3E","CD3G","CD6","SH2D1A","TRAT1"],
 "Th1 cells":["TBX21"],"Treg":["FOXP3"]}
CHECKPOINTS = ["PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT","BTLA","IDO1","CD276","TNFRSF9","ICOS","CD27","LGALS9"]

# ---------- tumour logCPM (symbols) ----------
df = pd.read_csv("data/TCGA-THCA.htseq_counts.tsv", sep="\t", index_col=0)
mat = (2**df-1).clip(lower=0)
stype = pd.Series([b.split("-")[3][:2] for b in df.columns], index=df.columns)
tum = df.columns[stype=="01"]
cpm = mat[tum].div(mat[tum].sum(axis=0),axis=1)*1e6
logcpm = np.log2(cpm[cpm.mean(axis=1)>=1]+1); logcpm.index=[g.split(".")[0] for g in logcpm.index]
e2s = pd.read_csv("data/thca_DE_full_symbol.csv")[["ensembl","symbol"]].dropna()
e2s = dict(zip(e2s.ensembl,e2s.symbol))
logcpm = logcpm[logcpm.index.isin(e2s)]; logcpm.index=[e2s[e] for e in logcpm.index]
logcpm = logcpm[~logcpm.index.duplicated()]
phl = logcpm.loc[GENE]

# ---------- ESTIMATE purity (R) ----------
logcpm.rename_axis("GeneSymbol").to_csv("data/estimate_input.tsv", sep="\t")
subprocess.run(["Rscript","code/run_estimate.R","data/estimate_input.tsv","data/estimate_scores.gct"], check=True)
sc = pd.read_csv("data/estimate_scores.gct", sep="\t", skiprows=2, index_col=0).drop(columns=["Description"])
sc.columns = [c.replace(".","-") for c in sc.columns]
est = sc.loc["ESTIMATEScore"]
# affymetrix-calibrated purity formula (widely applied to RNA-seq ESTIMATE scores)
purity = np.cos(0.6049872018 + 0.0001467884*est.astype(float))
purity.index = est.index

# align samples (tumour barcodes)
common = [c for c in logcpm.columns if c in purity.index]
phl = phl[common]; purity = purity[common]; logcpm = logcpm[common]
print(f"samples with purity: {len(common)}  purity range {purity.min():.2f}-{purity.max():.2f}")

# ---------- ssGSEA immune-cell scores ----------
ss = gp.ssgsea(data=logcpm, gene_sets=IMMUNE, sample_norm_method="rank", min_size=1, outdir=None, no_plot=True)
nes = ss.res2d.pivot(index="Term", columns="Name", values="NES").astype(float)[logcpm.columns]

# ---------- raw vs partial (purity) correlation ----------
def raw_partial(values, label):
    d = pd.DataFrame({"PHLDA3":phl.values, "x":values, "purity":purity.values}).dropna()
    r_raw,p_raw = spearmanr(d.PHLDA3, d.x)
    pc = pg.partial_corr(d, x="PHLDA3", y="x", covar="purity", method="spearman")
    pcol = [c for c in pc.columns if c.lower().replace("-","").replace("_","") in ("pval","p")][0]
    return {"item":label, "r_raw":r_raw, "p_raw":p_raw,
            "r_partial":float(pc["r"].iloc[0]), "p_partial":float(pc[pcol].iloc[0])}

rows_cp = [raw_partial(logcpm.loc[g].values, g) for g in CHECKPOINTS if g in logcpm.index]
rows_ic = [raw_partial(nes.loc[c].values, c) for c in nes.index]
cp = pd.DataFrame(rows_cp).sort_values("r_raw"); ic = pd.DataFrame(rows_ic).sort_values("r_raw")
for d in (cp,ic):
    d["survives"] = (d["p_partial"]<0.05) & (np.sign(d["r_raw"])==np.sign(d["r_partial"]))
pd.concat([cp.assign(type="checkpoint"), ic.assign(type="immune_cell")]).to_csv("data/PHLDA3_purity_partialcorr.csv", index=False)
print(f"checkpoints surviving purity adj: {cp.survives.sum()}/{len(cp)}")
print(f"immune cells surviving purity adj: {ic.survives.sum()}/{len(ic)}")

# ---------- figure: raw vs partial r (lollipop pairs) ----------
fig, ax = plt.subplots(1, 2, figsize=(11, 5))
def paired(a, d, key, title):
    y = np.arange(len(d))
    for i,(_,r) in enumerate(d.iterrows()):
        a.plot([r.r_raw, r.r_partial],[i,i], color=ns.C["grey_light"], lw=1.2, zorder=1)
        a.scatter(r.r_raw, i, s=30, color=ns.C["ns"], edgecolors="black", lw=0.4, zorder=2)
        sig = r.p_partial<0.05
        a.scatter(r.r_partial, i, s=42,
                  color=ns.C["tumor"] if (sig and r.r_partial>0) else ns.C["normal"] if (sig and r.r_partial<0) else "white",
                  edgecolors="black", lw=0.5, zorder=3)
    a.axvline(0, color="black", lw=0.7)
    a.set_yticks(y); a.set_yticklabels(d[key], fontsize=7.5)
    a.set_xlabel("Spearman r with PHLDA3"); a.set_title(title, fontweight="bold", fontsize=9)
    a.spines[["top","right"]].set_visible(False)
paired(ax[0], cp, "item", "Immune checkpoints")
paired(ax[1], ic, "item", "Immune-cell infiltration (ssGSEA)")
from matplotlib.lines import Line2D
ax[0].legend(handles=[Line2D([0],[0],marker="o",color="w",mfc=ns.C["ns"],mec="black",label="raw r",ms=6),
                      Line2D([0],[0],marker="o",color="w",mfc=ns.C["tumor"],mec="black",label="purity-adjusted (P<0.05)",ms=6)],
             loc="lower right", fontsize=6.5)
for a,l in zip(ax,["a","b"]): ns.panel_label(a,l)
fig.suptitle(f"PHLDA3 immune associations survive tumour-purity correction (ESTIMATE, n={len(common)})",
             fontweight="bold", fontsize=10, y=1.0)
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_purity")
print("done task 1")
