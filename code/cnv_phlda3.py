#!/usr/bin/env python3
"""Task 2a — PHLDA3 copy-number (GISTIC, cBioPortal) vs expression in TCGA-THCA."""
import pandas as pd, numpy as np, requests
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import kruskal, spearmanr
import nature_style as ns
ns.set_style(font_size=8)

GENE, ENS, ENTREZ = "PHLDA3", "ENSG00000174307", 23612
BASE, STUDY = "https://www.cbioportal.org/api", "thca_tcga_pan_can_atlas_2018"

# expression per patient
df = pd.read_csv("data/TCGA-THCA.htseq_counts.tsv", sep="\t", index_col=0)
mat = (2**df-1).clip(lower=0); row = mat.index[mat.index.str.startswith(ENS)][0]
logcpm = np.log2(mat.loc[row]/mat.sum(axis=0)*1e6+1)
stype = pd.Series([b.split("-")[3][:2] for b in df.columns], index=df.columns)
tum = logcpm[stype=="01"]
expr = pd.DataFrame({"patient":[b[:12] for b in tum.index], GENE:tum.values}).groupby("patient")[GENE].mean()

# GISTIC discrete CNA from cBioPortal
r = requests.post(f"{BASE}/molecular-profiles/{STUDY}_gistic/molecular-data/fetch",
    json={"sampleListId":f"{STUDY}_all","entrezGeneIds":[ENTREZ]},
    headers={"Content-Type":"application/json"}, timeout=60).json()
cna = pd.DataFrame([{"patient":x["sampleId"][:12],"cna":x["value"]} for x in r])
d = cna.merge(expr.rename(GENE), on="patient", how="inner").dropna()
LAB = {-2:"Deep del",-1:"Shallow del",0:"Diploid",1:"Gain",2:"Amp"}
d["cat"] = d["cna"].map(LAB)
order = [LAB[k] for k in [-2,-1,0,1,2] if LAB[k] in set(d["cat"])]
rho,p_sp = spearmanr(d["cna"], d[GENE])
groups = [d[d.cat==o][GENE].values for o in order]
p_kw = kruskal(*[g for g in groups if len(g)>0])[1]
d.to_csv("data/PHLDA3_cnv.csv", index=False)
print(f"n={len(d)}  CNA dist={d['cat'].value_counts().to_dict()}")
print(f"Spearman(CNA,expr) rho={rho:.2f} p={p_sp:.2e}  Kruskal p={p_kw:.2e}")

# figure
fig, ax = plt.subplots(figsize=(5,4.4))
cols = {"Deep del":ns.C["normal"],"Shallow del":"#9FC0E0","Diploid":ns.C["ns"],
        "Gain":"#E2A6A1","Amp":ns.C["tumor"]}
bp = ax.boxplot(groups, widths=0.6, patch_artist=True, showfliers=False,
                medianprops=dict(color="black",lw=1.2))
for patch,o in zip(bp["boxes"],order): patch.set_facecolor(cols[o]); patch.set_alpha(0.75)
for j,(g,o) in enumerate(zip(groups,order),1):
    ax.scatter(np.random.normal(j,0.07,len(g)), g, s=7, color=cols[o], alpha=0.5, edgecolors="none")
ax.set_xticks(range(1,len(order)+1)); ax.set_xticklabels([f"{o}\n(n={len(g)})" for o,g in zip(order,groups)])
ax.set_ylabel(f"{GENE} log2(CPM+1)")
ax.set_title(f"PHLDA3 copy number vs expression\nSpearman rho={rho:.2f}, Kruskal P={p_kw:.1e}",
             fontweight="bold", fontsize=9)
ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_cnv")
print("done task 2a (CNV)")
