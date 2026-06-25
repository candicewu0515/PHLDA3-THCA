#!/usr/bin/env python3
"""Reviewer Major-1 — expanded N1 model now INCLUDING RET/NTRK/ALK fusion status
(per-sample cBioPortal structural variants, data/thca_fusion_status.csv).
Does PHLDA3 stay associated with N1 after adjusting for BRAF AND driver fusions?
Also: PHLDA3 expression across the full driver landscape (BRAF/RAS/Fusion/WT)."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
from scipy.stats import kruskal, mannwhitneyu
import nature_style as ns
ns.set_style(font_size=8)

GENE, ENS = "PHLDA3", "ENSG00000174307"
# expression + clinical
df = pd.read_csv("data/TCGA-THCA.htseq_counts.tsv", sep="\t", index_col=0)
mat=(2**df-1).clip(lower=0); row=mat.index[mat.index.str.startswith(ENS)][0]
lc=np.log2(mat.loc[row]/mat.sum(axis=0)*1e6+1); stt=pd.Series([b.split("-")[3][:2] for b in df.columns],index=df.columns)
expr=pd.DataFrame({"patient":[b[:12] for b in lc[stt=="01"].index],GENE:lc[stt=="01"].values}).groupby("patient")[GENE].mean()
cl=pd.read_csv("data/THCA_clinical.tsv",sep="\t",dtype=str)
def pick(r,bb):
    for k in [f"diagnoses.0.{bb}",f"diagnoses.1.{bb}"]:
        if k in r and pd.notna(r[k]) and str(r[k]) not in("","nan"): return r[k]
    return np.nan
c=pd.DataFrame({"patient":cl["submitter_id"],"age":pd.to_numeric(cl["demographic.age_at_index"],errors="coerce"),
  "male":(cl["demographic.gender"]=="male").astype(float),
  "T":cl.apply(lambda r:pick(r,"ajcc_pathologic_t"),axis=1),"N":cl.apply(lambda r:pick(r,"ajcc_pathologic_n"),axis=1)})
c["T34"]=c["T"].str[:2].map({"T1":0,"T2":0,"T3":1,"T4":1}); c["N1"]=c["N"].str[:2].map({"N0":0,"N1":1})
# drivers: BRAF / RAS (from cache) + fusion (new)
br=pd.read_csv("data/PHLDA3_braf_ras.csv")[["patient","BRAF","RAS"]]
fus=pd.read_csv("data/thca_fusion_status.csv"); fset=set(fus.patient)
d=c.merge(expr.rename(GENE),on="patient").merge(br,on="patient")
d["BRAFmut"]=(d.BRAF=="BRAF-mut").astype(float); d["RASmut"]=(d.RAS=="RAS-mut").astype(float)
d["Fusion"]=d.patient.isin(fset).astype(float)
d["z"]=(d[GENE]-d[GENE].mean())/d[GENE].std()
# driver class (mutually exclusive priority: BRAF > RAS > Fusion > WT)
def cls(r):
    if r.BRAFmut: return "BRAF V600E"
    if r.RASmut:  return "RAS-mut"
    if r.Fusion:  return "Fusion+"
    return "WT/other"
d["class"]=d.apply(cls,axis=1)

# ---- expanded N1 model: +BRAF +Fusion ----
dm=d.dropna(subset=["N1","z","age","male","T34","BRAFmut","Fusion"])
def fit(cols):
    m=sm.Logit(dm["N1"].astype(int), sm.add_constant(dm[cols].astype(float))).fit(disp=0)
    return m
m_base=fit(["z","age","male","T34"])
m_braf=fit(["z","age","male","T34","BRAFmut"])
m_full=fit(["z","age","male","T34","BRAFmut","Fusion"])
def orr(m): return np.exp(m.params["z"]), *np.exp(m.conf_int().loc["z"]), m.pvalues["z"]
print("=== PHLDA3 OR for N1 (per 1 SD) ===")
for lab,m in [("base (PHLDA3+age+sex+T)",m_base),("+BRAF",m_braf),("+BRAF+Fusion",m_full)]:
    o=orr(m); print(f"  {lab:26s} OR={o[0]:.2f} [{o[1]:.2f}-{o[2]:.2f}] p={o[3]:.2e}  (n={int(dm.shape[0])})")
print(f"  Fusion+ N1 covariate OR={np.exp(m_full.params['Fusion']):.2f} p={m_full.pvalues['Fusion']:.3f}")
pd.DataFrame([dict(model=l,OR=orr(m)[0],lo=orr(m)[1],hi=orr(m)[2],p=orr(m)[3],n=dm.shape[0])
   for l,m in [("PHLDA3+age+sex+T",m_base),("+BRAF",m_braf),("+BRAF+Fusion",m_full)]]).to_csv("data/PHLDA3_expanded_n1.csv",index=False)

# ---- PHLDA3 by driver class ----
order=["BRAF V600E","WT/other","Fusion+","RAS-mut"]
groups=[d[d["class"]==k][GENE].values for k in order]
H,pk=kruskal(*groups)
print("\n=== PHLDA3 by driver class (median) ===")
for k,g in zip(order,groups): print(f"  {k:12s} n={len(g):3d}  median={np.median(g):.2f}")
print(f"  Kruskal-Wallis p={pk:.2e};  Fusion+ vs WT MWU p={mannwhitneyu(d[d['class']=='Fusion+'][GENE],d[d['class']=='WT/other'][GENE])[1]:.3f}")

# figure
fig,ax=plt.subplots(1,2,figsize=(7.4,3.2))
cols={"BRAF V600E":ns.C["tumor"],"WT/other":ns.C["ns"],"Fusion+":ns.C["accent"],"RAS-mut":ns.C["normal"]}
bp=ax[0].boxplot(groups,widths=.62,patch_artist=True,showfliers=False,medianprops=dict(color="black",lw=1.1))
for patch,k in zip(bp["boxes"],order): patch.set_facecolor(cols[k]); patch.set_alpha(.65)
ax[0].set_xticks(range(1,5)); ax[0].set_xticklabels([f"{k}\n(n={len(g)})" for k,g in zip(order,groups)],fontsize=7)
ax[0].set_ylabel(f"{GENE} log2(CPM+1)"); ax[0].set_title(f"PHLDA3 across driver classes\nKruskal–Wallis P={pk:.1e}",fontweight="bold",fontsize=8.5)
ax[0].spines[["top","right"]].set_visible(False)
# forest of PHLDA3 OR across nested models
rows=[("+BRAF+Fusion",orr(m_full)),("+BRAF",orr(m_braf)),("base",orr(m_base))]
for i,(lab,o) in enumerate(rows):
    ax[1].plot([o[1],o[2]],[i,i],color=ns.C["tumor"],lw=2); ax[1].scatter(o[0],i,s=55,color=ns.C["tumor"],ec="black",lw=.6,zorder=3)
    ax[1].text(o[2],i,f"  {o[0]:.2f} ({o[1]:.2f}-{o[2]:.2f})",va="center",fontsize=7.5)
ax[1].axvline(1,ls="--",color=ns.C["grey_mid"],lw=1); ax[1].set_yticks(range(3)); ax[1].set_yticklabels([r[0] for r in rows])
ax[1].set_xscale("log"); ax[1].set_xlabel("PHLDA3 OR for N1 (per 1 SD)"); ax[1].margins(x=.5,y=.25)
ax[1].set_title("PHLDA3–N1 OR attenuates to n.s.\nafter adding BRAF + RET/NTRK fusion",fontweight="bold",fontsize=8.5)
ax[1].spines[["top","right"]].set_visible(False)
for a,l in zip(ax,"ab"): ns.panel_label(a,l)
plt.tight_layout(); ns.save_fig(fig,"PHLDA3_THCA_expanded_n1_fusion")
print("done expanded N1 + fusion")
