#!/usr/bin/env python3
"""Composite main Figure 2 (driver-aware N1 analysis) — Nature-style:
(a) base multivariable N1 forest, (b) PHLDA3-OR attenuation across nested models
(base -> +BRAF -> +BRAF+fusion), (c) PHLDA3 across mutually exclusive driver
classes, (d) incremental 10-fold CV-AUC (PHLDA3 adds ~0 over driver status).
The nomogram/ROC/calibration/DCA are exploratory and shown in Supplementary."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
from sklearn.model_selection import cross_val_predict
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from scipy.stats import kruskal
import nature_style as ns
ns.set_style(font_size=7)

GENE, ENS = "PHLDA3", "ENSG00000174307"
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
br=pd.read_csv("data/PHLDA3_braf_ras.csv")[["patient","BRAF","RAS"]]
fset=set(pd.read_csv("data/thca_fusion_status.csv").patient)
d=c.merge(expr.rename(GENE),on="patient").merge(br,on="patient")
d["BRAFmut"]=(d.BRAF=="BRAF-mut").astype(float); d["RASmut"]=(d.RAS=="RAS-mut").astype(float)
d["Fusion"]=d.patient.isin(fset).astype(float); d["z"]=(d[GENE]-d[GENE].mean())/d[GENE].std()
# mutually exclusive driver class: fusion first, then BRAF, then RAS, else WT (overlaps verified ~0)
def cls(r):
    if r.Fusion: return "Fusion+"
    if r.BRAFmut: return "BRAF V600E"
    if r.RASmut:  return "RAS-mut"
    return "WT/other"
d["class"]=d.apply(cls,axis=1)
dm=d.dropna(subset=["N1","z","age","male","T34","BRAFmut","Fusion"])

def OR(cols,var="z"):
    m=sm.Logit(dm.N1.astype(int),sm.add_constant(dm[cols].astype(float))).fit(disp=0)
    return np.exp(m.params[var]),*np.exp(m.conf_int().loc[var]),m.pvalues[var]
def cvauc(cols):
    return roc_auc_score(dm.N1, cross_val_predict(LogisticRegression(max_iter=2000),
            dm[cols].values,dm.N1.values,cv=10,method="predict_proba")[:,1])

fig=plt.figure(figsize=(7.2,6.0)); gs=fig.add_gridspec(2,2,hspace=0.55,wspace=0.5)
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1]); axC=fig.add_subplot(gs[1,0]); axD=fig.add_subplot(gs[1,1])

# (a) base multivariable forest
mb=sm.Logit(dm.N1.astype(int),sm.add_constant(dm[["z","age","male","T34"]].astype(float))).fit(disp=0)
labs={"z":"PHLDA3","age":"Age","male":"Male","T34":"T3/T4"}
fr=[(labs[v],np.exp(mb.params[v]),*np.exp(mb.conf_int().loc[v])) for v in ["z","age","male","T34"]][::-1]
for i,(l,o,lo,hi) in enumerate(fr):
    h=l=="PHLDA3"; col=ns.C["tumor"] if h else ns.C["normal"]
    axA.plot([lo,hi],[i,i],color=col,lw=1.6); axA.scatter(o,i,s=34 if h else 22,color=col,ec="black",lw=.4,zorder=3)
axA.axvline(1,ls="--",color=ns.C["grey_mid"],lw=.9); axA.set_yticks(range(4)); axA.set_yticklabels([f[0] for f in fr],fontsize=6.5)
axA.set_xscale("log"); axA.set_xlabel("OR for N1 (95% CI)",fontsize=6.5); axA.margins(x=.3)
axA.set_title("Base N1 model (PHLDA3+age+sex+T)",fontweight="bold",fontsize=7); axA.spines[["top","right"]].set_visible(False)

# (b) PHLDA3-OR attenuation across nested models
nest=[("base",["z","age","male","T34"]),("+BRAF",["z","age","male","T34","BRAFmut"]),
      ("+BRAF+fusion",["z","age","male","T34","BRAFmut","Fusion"])]
rows=[(l,*OR(cols)) for l,cols in nest][::-1]
for i,(l,o,lo,hi,p) in enumerate(rows):
    sig=p<0.05; col=ns.C["tumor"] if sig else ns.C["ns"]
    axB.plot([lo,hi],[i,i],color=col,lw=1.8); axB.scatter(o,i,s=42,color=col,ec="black",lw=.5,zorder=3)
    ptxt="P<0.001" if p<0.001 else f"P={p:.2f}"
    axB.text(hi,i,f"  {o:.2f} ({ptxt})",va="center",fontsize=6.3)
axB.axvline(1,ls="--",color=ns.C["grey_mid"],lw=.9); axB.set_yticks(range(3)); axB.set_yticklabels([r[0] for r in rows],fontsize=6.5)
axB.set_xscale("log"); axB.set_xlabel("PHLDA3 OR for N1 (per 1 SD)",fontsize=6.5); axB.margins(x=.55,y=.22)
axB.set_title("PHLDA3–N1 OR attenuates to n.s.\nadding BRAF + RET/NTRK/ALK fusion",fontweight="bold",fontsize=7); axB.spines[["top","right"]].set_visible(False)

# (c) PHLDA3 by mutually exclusive driver class
order=["BRAF V600E","Fusion+","RAS-mut","WT/other"]; groups=[d[d["class"]==k][GENE].values for k in order]
pk=kruskal(*groups)[1]; cols={"BRAF V600E":ns.C["tumor"],"Fusion+":ns.C["accent"],"RAS-mut":ns.C["normal"],"WT/other":ns.C["ns"]}
bp=axC.boxplot(groups,widths=.62,patch_artist=True,showfliers=False,medianprops=dict(color="black",lw=1))
for patch,k in zip(bp["boxes"],order): patch.set_facecolor(cols[k]); patch.set_alpha(.65)
axC.set_xticks(range(1,5)); axC.set_xticklabels([f"{k}\n(n={len(g)})" for k,g in zip(order,groups)],fontsize=6)
axC.set_ylabel(f"{GENE} log2(CPM+1)",fontsize=6.5); axC.set_title(f"PHLDA3 by driver class (KW P={pk:.0e})",fontweight="bold",fontsize=7)
axC.spines[["top","right"]].set_visible(False)

# (d) incremental CV-AUC
base=["age","male","T34"]; drv=base+["BRAFmut","Fusion"]
a_drv=cvauc(drv); a_drvg=cvauc(drv+["z"])
aucs=[("Clinical\n(age+sex+T)",cvauc(base)),("+PHLDA3",cvauc(base+["z"])),
      ("Driver-aware\n(+BRAF+fusion)",a_drv),("Driver-aware\n+PHLDA3",a_drvg)]
bcol=[ns.C["ns"],ns.C["tumor"],ns.C["normal"],ns.C["accent"]]
axD.bar(range(4),[a[1] for a in aucs],color=bcol,ec="black",lw=.5,width=.7)
for i,a in enumerate(aucs): axD.text(i,a[1]+.005,f"{a[1]:.3f}",ha="center",fontsize=6)
axD.set_xticks(range(4)); axD.set_xticklabels([a[0] for a in aucs],fontsize=5.6); axD.set_ylim(0.5,0.78)
axD.set_ylabel("10-fold CV AUC",fontsize=6.5); axD.axhline(0.5,color="black",lw=.6,ls=":")
axD.set_title(f"PHLDA3 adds ΔAUC={a_drvg-a_drv:+.3f} over driver status",fontweight="bold",fontsize=6.8)
axD.spines[["top","right"]].set_visible(False)

for ax,l in zip([axA,axB,axC,axD],"abcd"): ns.panel_label(ax,l,x=-0.16,y=1.06,fs=10)
fig.suptitle("Figure 2. PHLDA3's lymph-node-metastasis association is shared with the kinase-driver landscape",
             fontweight="bold",fontsize=8,y=1.0)
ns.save_fig(fig,"PHLDA3_THCA_Figure2")
print(f"done Figure 2 (driver-aware); driver-aware AUC {a_drv:.3f} -> +PHLDA3 {a_drvg:.3f}")
