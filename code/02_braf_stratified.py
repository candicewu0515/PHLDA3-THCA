#!/usr/bin/env python3
"""Reviewer Major-2 — BRAF-stratified N1 analysis + PHLDA3xBRAF interaction.
Does PHLDA3 predict N1 within BRAF-mutant and within BRAF-wildtype tumours?
Outputs a subgroup forest (Supplementary) + interaction test."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
import nature_style as ns
ns.set_style(font_size=8)

GENE, ENS = "PHLDA3", "ENSG00000174307"
df = pd.read_csv("data/TCGA-THCA.htseq_counts.tsv", sep="\t", index_col=0)
mat=(2**df-1).clip(lower=0); row=mat.index[mat.index.str.startswith(ENS)][0]
logcpm=np.log2(mat.loc[row]/mat.sum(axis=0)*1e6+1)
st=pd.Series([b.split("-")[3][:2] for b in df.columns],index=df.columns); tum=logcpm[st=="01"]
expr=pd.DataFrame({"patient":[b[:12] for b in tum.index],GENE:tum.values}).groupby("patient")[GENE].mean()
cl=pd.read_csv("data/THCA_clinical.tsv",sep="\t",dtype=str)
def pick(r,b):
    for k in [f"diagnoses.0.{b}",f"diagnoses.1.{b}"]:
        if k in r and pd.notna(r[k]) and str(r[k]) not in("","nan"): return r[k]
    return np.nan
c=pd.DataFrame({"patient":cl["submitter_id"],"age":pd.to_numeric(cl["demographic.age_at_index"],errors="coerce"),
  "male":(cl["demographic.gender"]=="male").astype(float),
  "T":cl.apply(lambda r:pick(r,"ajcc_pathologic_t"),axis=1),"N":cl.apply(lambda r:pick(r,"ajcc_pathologic_n"),axis=1)})
c["T34"]=c["T"].str[:2].map({"T1":0,"T2":0,"T3":1,"T4":1}); c["N1"]=c["N"].str[:2].map({"N0":0,"N1":1})
br=pd.read_csv("data/PHLDA3_braf_ras.csv")[["patient","BRAF"]]; br["BRAFmut"]=(br["BRAF"]=="BRAF-mut").astype(float)
d=c.merge(expr.rename(GENE),on="patient").merge(br[["patient","BRAFmut"]],on="patient")
d["z"]=(d[GENE]-d[GENE].mean())/d[GENE].std()
d=d.dropna(subset=["N1","z","age","male","T34","BRAFmut"])

def fit(sub):
    m=sm.Logit(sub["N1"].astype(int), sm.add_constant(sub[["z","age","male","T34"]].astype(float))).fit(disp=0)
    return np.exp(m.params["z"]), *np.exp(m.conf_int().loc["z"]), m.pvalues["z"], len(sub)
mut=fit(d[d.BRAFmut==1]); wt=fit(d[d.BRAFmut==0]); ov=fit(d)
# interaction
mi=smf.logit("N1 ~ z*BRAFmut + age + male + T34", data=d).fit(disp=0)
pint=mi.pvalues["z:BRAFmut"]
print("Overall   : OR=%.2f [%.2f-%.2f] p=%.2e n=%d"%ov)
print("BRAF-mut  : OR=%.2f [%.2f-%.2f] p=%.2e n=%d"%mut)
print("BRAF-wt   : OR=%.2f [%.2f-%.2f] p=%.2e n=%d"%wt)
print("PHLDA3xBRAF interaction p = %.3f"%pint)
pd.DataFrame([dict(zip(["OR","lo","hi","p","n"],x))|{"group":g} for g,x in
   [("Overall",ov),("BRAF-mutant",mut),("BRAF-wildtype",wt)]]).to_csv("data/PHLDA3_braf_stratified.csv",index=False)

# forest
fig,ax=plt.subplots(figsize=(6.2,2.8))
rows=[("BRAF-wildtype",wt),("BRAF-mutant",mut),("Overall",ov)]
for i,(lab,r) in enumerate(rows):
    OR,lo,hi,p,n=r
    ax.plot([lo,hi],[i,i],color=ns.C["tumor"],lw=2); ax.scatter(OR,i,s=60,color=ns.C["tumor"],ec="black",lw=.6,zorder=3)
    ax.text(ax.get_xlim()[1] if i==0 else hi,i,f"  OR {OR:.2f} ({lo:.2f}-{hi:.2f}), P={p:.0e}, n={n}",va="center",fontsize=7.5)
ax.axvline(1,ls="--",color=ns.C["grey_mid"],lw=1)
ax.set_yticks(range(3)); ax.set_yticklabels([r[0] for r in rows]); ax.set_xscale("log")
ax.set_xlabel(f"{GENE} odds ratio for N1 (per 1 SD)"); ax.set_ylim(-.5,2.5); ax.margins(x=.55)
ax.set_title(f"PHLDA3 predicts N1 within BRAF subgroups\nPHLDA3xBRAF interaction P = {pint:.2f} (no significant interaction)",
             fontweight="bold",fontsize=8.5); ax.spines[["top","right"]].set_visible(False)
plt.tight_layout(); ns.save_fig(fig,"PHLDA3_THCA_braf_stratified")
print("done BRAF-stratified")
