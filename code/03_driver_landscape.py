#!/usr/bin/env python3
"""Main Figure 3 — PHLDA3 across the full THCA driver landscape (4 mutually
exclusive classes: BRAF V600E, RET/NTRK/ALK fusion-positive, RAS-mutant,
WT/other). Goes beyond Fig 2c (model interpretation) to driver biology:
(a) PHLDA3 expression, (b) N1 frequency, (c) promoter methylation cg04055835,
(d) integrated PHLDA3 vs promoter methylation, coloured by driver class."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import kruskal, spearmanr
import nature_style as ns
ns.set_style(font_size=7)

GENE, ENS = "PHLDA3", "ENSG00000174307"
df=pd.read_csv("data/TCGA-THCA.htseq_counts.tsv",sep="\t",index_col=0)
mat=(2**df-1).clip(lower=0); row=mat.index[mat.index.str.startswith(ENS)][0]
lc=np.log2(mat.loc[row]/mat.sum(axis=0)*1e6+1); stt=pd.Series([b.split("-")[3][:2] for b in df.columns],index=df.columns)
expr=pd.DataFrame({"patient":[b[:12] for b in lc[stt=="01"].index],GENE:lc[stt=="01"].values}).groupby("patient")[GENE].mean()
cl=pd.read_csv("data/THCA_clinical.tsv",sep="\t",dtype=str)
def pick(r,bb):
    for k in [f"diagnoses.0.{bb}",f"diagnoses.1.{bb}"]:
        if k in r and pd.notna(r[k]) and str(r[k]) not in("","nan"): return r[k]
    return np.nan
N=pd.DataFrame({"patient":cl["submitter_id"],"N":cl.apply(lambda r:pick(r,"ajcc_pathologic_n"),axis=1)})
N["N1"]=N["N"].str[:2].map({"N0":0,"N1":1})
br=pd.read_csv("data/PHLDA3_braf_ras.csv")[["patient","BRAF","RAS"]]
fset=set(pd.read_csv("data/thca_fusion_status.csv").patient)
me=pd.read_csv("data/cg04055835_betas.csv"); me["patient"]=me["sample"].str[:12]
meth=me.groupby("patient").beta.mean()

d=expr.rename(GENE).reset_index().merge(br,on="patient").merge(N[["patient","N1"]],on="patient",how="left")
d["BRAFmut"]=(d.BRAF=="BRAF-mut"); d["RASmut"]=(d.RAS=="RAS-mut"); d["Fusion"]=d.patient.isin(fset)
def cls(r):
    if r.Fusion: return "Fusion+"
    if r.BRAFmut: return "BRAF V600E"
    if r.RASmut:  return "RAS-mut"
    return "WT/other"
d["class"]=d.apply(cls,axis=1)
d=d.merge(meth.rename("meth"),on="patient",how="left")
order=["BRAF V600E","Fusion+","RAS-mut","WT/other"]
COL={"BRAF V600E":ns.C["tumor"],"Fusion+":ns.C["accent"],"RAS-mut":ns.C["normal"],"WT/other":ns.C["ns"]}

fig=plt.figure(figsize=(7.2,5.8)); gs=fig.add_gridspec(2,2,hspace=0.5,wspace=0.42)
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1]); axC=fig.add_subplot(gs[1,0]); axD=fig.add_subplot(gs[1,1])

def box(ax,col,ylab,title):
    g=[d[d["class"]==k][col].dropna().values for k in order]; p=kruskal(*g)[1]
    bp=ax.boxplot(g,widths=.62,patch_artist=True,showfliers=False,medianprops=dict(color="black",lw=1))
    for patch,k in zip(bp["boxes"],order): patch.set_facecolor(COL[k]); patch.set_alpha(.65)
    for j,(gg,k) in enumerate(zip(g,order),1): ax.scatter(np.random.normal(j,.07,len(gg)),gg,s=3,color=COL[k],alpha=.35,ec="none",zorder=3)
    ax.set_xticks(range(1,5)); ax.set_xticklabels([f"{k}\n(n={len(g[i])})" for i,k in enumerate(order)],fontsize=5.8)
    ax.set_ylabel(ylab,fontsize=6.5); ax.set_title(f"{title} (KW P={p:.0e})",fontweight="bold",fontsize=7)
    ax.spines[["top","right"]].set_visible(False); return p

# (a) PHLDA3 expression
box(axA,GENE,f"{GENE} log2(CPM+1)","PHLDA3 expression by driver class")
# (b) N1 frequency
nf=d.dropna(subset=["N1"]).groupby("class").N1.agg(["mean","size"]).reindex(order)
axB.bar(range(4),nf["mean"]*100,color=[COL[k] for k in order],ec="black",lw=.5,width=.7)
for i,(m,n) in enumerate(zip(nf["mean"],nf["size"])): axB.text(i,m*100+1.5,f"{m*100:.0f}%",ha="center",fontsize=6)
axB.set_xticks(range(4)); axB.set_xticklabels([f"{k}\n(n={int(n)})" for k,n in zip(order,nf["size"])],fontsize=5.8)
axB.set_ylabel("N1 (lymph-node met.) frequency (%)",fontsize=6.3); axB.set_ylim(0,100)
axB.set_title("Nodal-metastasis frequency by class",fontweight="bold",fontsize=7); axB.spines[["top","right"]].set_visible(False)
# (c) promoter methylation cg04055835
box(axC,"meth","cg04055835 β (promoter)","Promoter methylation by driver class")
# (d) integrated PHLDA3 vs promoter methylation
dd=d.dropna(subset=["meth"]); rho,_=spearmanr(dd.meth,dd[GENE])
for k in order: s=dd[dd["class"]==k]; axD.scatter(s.meth,s[GENE],s=12,color=COL[k],alpha=.7,ec="none",label=k)
b,a=np.polyfit(dd.meth,dd[GENE],1); xs=np.linspace(dd.meth.min(),dd.meth.max(),40)
axD.plot(xs,b*xs+a,color="black",lw=1,ls=":")
axD.set_xlabel("cg04055835 β (promoter methylation)",fontsize=6.5); axD.set_ylabel(f"{GENE} log2(CPM+1)",fontsize=6.5)
axD.set_title(f"Aggressive classes: high PHLDA3, low methylation (r={rho:.2f})",fontweight="bold",fontsize=6.6)
axD.legend(frameon=False,fontsize=5.5,loc="upper right",handletextpad=.2); axD.spines[["top","right"]].set_visible(False)

for ax,l in zip([axA,axB,axC,axD],"abcd"): ns.panel_label(ax,l,x=-0.16,y=1.06,fs=10)
fig.suptitle("Figure 3. PHLDA3 marks the aggressive driver landscape of thyroid carcinoma",fontweight="bold",fontsize=8,y=1.0)
ns.save_fig(fig,"PHLDA3_THCA_Figure3")
print("done Figure 3 (driver landscape)")
print(d.groupby("class")[GENE].median().reindex(order).to_string())
print("N1 freq:", (nf["mean"]*100).round(0).to_dict())
