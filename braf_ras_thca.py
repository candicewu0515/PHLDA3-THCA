#!/usr/bin/env python3
"""PHLDA3 expression vs THCA driver mutations (BRAF V600E / RAS), TCGA via cBioPortal."""
import pandas as pd, numpy as np, requests
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, kruskal
import nature_style as ns

ns.set_style(font_size=8)

GENE, ENS = "PHLDA3", "ENSG00000174307"
BASE = "https://www.cbioportal.org/api"
STUDY = "thca_tcga_pan_can_atlas_2018"

# ---------- PHLDA3 expression per patient ----------
df = pd.read_csv("TCGA-THCA.htseq_counts.tsv", sep="\t", index_col=0)
mat = (2**df - 1).clip(lower=0)
row = mat.index[mat.index.str.startswith(ENS)][0]
logcpm = np.log2(mat.loc[row] / mat.sum(axis=0) * 1e6 + 1)
stype = pd.Series([b.split("-")[3][:2] for b in df.columns], index=df.columns)
tum = logcpm[stype == "01"]
expr = pd.DataFrame({"patient":[b[:12] for b in tum.index], GENE:tum.values}).groupby("patient")[GENE].mean()

# ---------- mutations from cBioPortal ----------
seq = requests.get(f"{BASE}/sample-lists/{STUDY}_sequenced/sample-ids", timeout=30).json()
ENTREZ = {673:"BRAF", 4893:"NRAS", 3265:"HRAS", 3845:"KRAS"}
muts = requests.post(f"{BASE}/molecular-profiles/{STUDY}_mutations/mutations/fetch",
    json={"sampleListId":f"{STUDY}_sequenced","entrezGeneIds":list(ENTREZ)},
    headers={"Content-Type":"application/json"}, timeout=60).json()

braf_v600e, braf_any, ras_any = set(), set(), set()
for m in muts:
    pid = m["sampleId"][:12]; g = ENTREZ.get(m["entrezGeneId"])
    if g == "BRAF":
        braf_any.add(pid)
        if str(m.get("proteinChange","")).startswith("V600"): braf_v600e.add(pid)
    else:
        ras_any.add(pid)
seq_pat = {s[:12] for s in seq}
print(f"sequenced patients={len(seq_pat)}  BRAF={len(braf_any)} (V600E={len(braf_v600e)})  RAS={len(ras_any)}")

d = expr.reset_index()
d = d[d.patient.isin(seq_pat)].copy()
d["BRAF"] = np.where(d.patient.isin(braf_any), "BRAF-mut", "BRAF-WT")
d["RAS"]  = np.where(d.patient.isin(ras_any),  "RAS-mut",  "RAS-WT")
def driver(p):
    if p in braf_v600e: return "BRAF V600E"
    if p in ras_any:    return "RAS"
    return "WT/other"
d["driver"] = d.patient.map(driver)
d.to_csv("PHLDA3_braf_ras.csv", index=False)

# ---------- figure ----------
fig, ax = plt.subplots(1, 3, figsize=(13, 4.4))
def box(a, groups, labels, colors, title, p, ptype="p"):
    bp=a.boxplot(groups, widths=0.55, patch_artist=True, showfliers=False,
                 medianprops=dict(color="black",lw=1.3))
    for patch,c in zip(bp["boxes"],colors): patch.set_facecolor(c); patch.set_alpha(0.6)
    for j,(g,c) in enumerate(zip(groups,colors),1):
        a.scatter(np.random.normal(j,0.06,len(g)), g, s=6, color=c, alpha=0.4, edgecolors="none")
    a.set_xticks(range(1,len(groups)+1))
    a.set_xticklabels([f"{l}\n(n={len(g)})" for l,g in zip(labels,groups)])
    a.set_ylabel(f"{GENE} log2(CPM+1)")
    pstr = f"{ptype} = {p:.1e}" if p>=1e-4 else f"{ptype} < 1e-4"
    a.set_title(f"{title}\n{pstr}", fontweight="bold", fontsize=11)
    a.spines[["top","right"]].set_visible(False)

# A: BRAF
g1=[d[d.BRAF=="BRAF-WT"][GENE].values, d[d.BRAF=="BRAF-mut"][GENE].values]
box(ax[0], g1, ["BRAF-WT","BRAF-mut"], [ns.C["normal"],ns.C["tumor"]], "BRAF mutation",
    mannwhitneyu(*g1,alternative="two-sided")[1])
# B: RAS
g2=[d[d.RAS=="RAS-WT"][GENE].values, d[d.RAS=="RAS-mut"][GENE].values]
box(ax[1], g2, ["RAS-WT","RAS-mut"], [ns.C["normal"],ns.C["accent"]], "RAS mutation",
    mannwhitneyu(*g2,alternative="two-sided")[1])
# C: driver groups
order=["WT/other","RAS","BRAF V600E"]
g3=[d[d.driver==o][GENE].values for o in order]
box(ax[2], g3, order, [ns.C["ns"],ns.C["accent"],ns.C["tumor"]], "Driver subtype",
    kruskal(*g3)[1], ptype="Kruskal p")
for i, lab in enumerate(["a","b","c"]):
    ns.panel_label(ax[i], lab)
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_braf_ras")
# print medians
for o in order: print(f"  {o:12s} n={len(d[d.driver==o])} median={np.median(d[d.driver==o][GENE]):.2f}")
print("saved table: PHLDA3_braf_ras.csv")
