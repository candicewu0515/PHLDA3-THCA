#!/usr/bin/env python3
"""Task 2b — PHLDA3 promoter methylation (HM450) vs expression in TCGA-THCA.
Tests whether PHLDA3 over-expression is associated with promoter hypomethylation
(epigenetic activation), complementing the CNV analysis (mostly diploid)."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, mannwhitneyu
import xenaPython as xena
import nature_style as ns
ns.set_style(font_size=8)

GENE, ENS = "PHLDA3", "ENSG00000174307"
TSS = 201469430   # PHLDA3 GRCh38 TSS (Ensembl: minus strand, gene end = TSS)

# ---------- methylation betas (13 PHLDA3 probes x samples) ----------
host, ds = "https://gdc.xenahubs.net", "TCGA-THCA.methylation450.tsv"
msamp = xena.dataset_samples(host, ds, None)
r = xena.dataset_gene_probes_values(host, ds, msamp, [GENE])
probes = r[0]["name"]; pos = [p["chromstart"] for p in r[0]["position"]]
beta = pd.DataFrame(np.array(r[1]), index=probes, columns=msamp).astype(float)
dist = {pr: abs(TSS - p) for pr, p in zip(probes, pos)}     # distance to TSS
# promoter probes = within 1500 bp of TSS
prom = [pr for pr in probes if dist[pr] <= 1500]
print(f"PHLDA3 probes: {len(probes)}, promoter(<=1500bp TSS): {len(prom)}")

# ---------- expression per patient ----------
df = pd.read_csv("TCGA-THCA.htseq_counts.tsv", sep="\t", index_col=0)
mat = (2**df-1).clip(lower=0); row = mat.index[mat.index.str.startswith(ENS)][0]
logcpm = np.log2(mat.loc[row]/mat.sum(axis=0)*1e6+1)
stype = pd.Series([b.split("-")[3][:2] for b in df.columns], index=df.columns)
expr = pd.DataFrame({"patient":[b[:12] for b in logcpm.index], "type":stype.values, GENE:logcpm.values})
expr_t = expr[expr.type=="01"].groupby("patient")[GENE].mean()

# methylation per patient/type
mcols = pd.DataFrame({"s":msamp, "patient":[s[:12] for s in msamp], "type":[s[12:15].strip("-")[:2] for s in msamp]})
def beta_patient(prlist, typ):
    cols = mcols[mcols.type==typ]["s"].tolist()
    sub = beta.loc[prlist, [c for c in cols if c in beta.columns]].mean(axis=0)
    out = pd.DataFrame({"patient":[c[:12] for c in sub.index], "beta":sub.values}).groupby("patient")["beta"].mean()
    return out
prom_t = beta_patient(prom, "01"); prom_n = beta_patient(prom, "11")

# ---------- per-probe: tumor-vs-normal + corr with expression ----------
rows=[]
for pr in probes:
    bt = beta_patient([pr],"01"); bn = beta_patient([pr],"11")
    d = pd.DataFrame({"beta":bt}).join(expr_t.rename("expr"), how="inner").dropna()
    rho,p = spearmanr(d.beta, d.expr) if len(d)>10 else (np.nan,np.nan)
    pmw = mannwhitneyu(bt.dropna(), bn.dropna())[1] if len(bn.dropna())>3 else np.nan
    rows.append({"probe":pr, "dist_TSS":dist[pr], "promoter":pr in prom,
                 "rho_expr":rho, "p_expr":p, "tumor_beta":bt.mean(), "normal_beta":bn.mean(), "p_TvsN":pmw})
pt = pd.DataFrame(rows).sort_values("dist_TSS")
pt.to_csv("PHLDA3_methylation_probes.csv", index=False)

# lead regulatory probe = promoter CpG (<=1500bp TSS) most inversely correlated with expression
cand = pt[pt.promoter & pt.rho_expr.notna()]
LEAD = cand.loc[cand.rho_expr.idxmin(), "probe"]
lead_t = beta_patient([LEAD], "01"); lead_n = beta_patient([LEAD], "11")
dprom = pd.DataFrame({"beta":lead_t}).join(expr_t.rename("expr"), how="inner").dropna()
rho_p, p_p = spearmanr(dprom.beta, dprom.expr)
p_tn = mannwhitneyu(lead_t.dropna(), lead_n.dropna())[1]
print(f"lead promoter CpG = {LEAD}")
print(f"  meth vs expr: rho={rho_p:.2f} p={p_p:.1e}")
print(f"  tumor={lead_t.mean():.2f} normal={lead_n.mean():.2f} p={p_tn:.1e}")
prom_t, prom_n = lead_t, lead_n

# ================= figure =================
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.4))
# A: per-probe Spearman r with expression (ordered by TSS distance)
cols = [ns.C["tumor"] if pr_ in prom else ns.C["ns"] for pr_ in pt["probe"]]
ax[0].barh(range(len(pt)), pt["rho_expr"], color=cols, edgecolor="black", lw=0.4)
ax[0].set_yticks(range(len(pt))); ax[0].set_yticklabels(pt["probe"], fontsize=6)
ax[0].axvline(0, color="black", lw=0.8); ax[0].set_xlabel("Spearman r (methylation vs expr)")
ax[0].set_title("Per-probe methylation–expression\n(red = promoter)", fontweight="bold", fontsize=9)
ax[0].spines[["top","right"]].set_visible(False)
# B: promoter methylation vs expression scatter
ax[1].scatter(dprom.beta, dprom.expr, s=9, color=ns.C["tumor"], alpha=0.5, edgecolors="none")
m,b = np.polyfit(dprom.beta, dprom.expr, 1); xs=np.linspace(dprom.beta.min(),dprom.beta.max(),50)
ax[1].plot(xs, m*xs+b, color="black", lw=1.2)
ax[1].set_xlabel(f"Promoter CpG {LEAD} methylation (β)"); ax[1].set_ylabel(f"{GENE} log2(CPM+1)")
ax[1].set_title(f"Promoter CpG vs expression\nSpearman r={rho_p:.2f}, P={p_p:.1e}", fontweight="bold", fontsize=9)
ax[1].spines[["top","right"]].set_visible(False)
# C: promoter methylation tumor vs normal
g=[prom_n.dropna().values, prom_t.dropna().values]
bp=ax[2].boxplot(g, widths=0.6, patch_artist=True, showfliers=False, medianprops=dict(color="black",lw=1.2))
for patch,c in zip(bp["boxes"],[ns.C["normal"],ns.C["tumor"]]): patch.set_facecolor(c); patch.set_alpha(0.7)
for j,(gg,c) in enumerate(zip(g,[ns.C["normal"],ns.C["tumor"]]),1):
    ax[2].scatter(np.random.normal(j,0.07,len(gg)), gg, s=7, color=c, alpha=0.5, edgecolors="none")
ax[2].set_xticks([1,2]); ax[2].set_xticklabels([f"Normal\n(n={len(g[0])})",f"Tumor\n(n={len(g[1])})"])
ax[2].set_ylabel(f"Promoter CpG {LEAD} methylation (β)")
ax[2].set_title(f"Promoter CpG: tumor vs normal\nMWU P={p_tn:.1e}", fontweight="bold", fontsize=9)
ax[2].spines[["top","right"]].set_visible(False)
for a,l in zip(ax,["a","b","c"]): ns.panel_label(a,l)
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_methylation")
print("done task 2b")
