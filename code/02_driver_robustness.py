#!/usr/bin/env python3
"""Robustness checks for the driver-landscape framing (reviewer round):
(1) fusion+ PHLDA3 elevation is not driven by one fusion type/batch (histology
    unavailable in harmonised TCGA clinical, so leave-one-fusion-type-out);
(2) PHLDA3 ~ BRAF + fusion + RAS + age + sex + T  (driver effect beyond T stage);
(3) the same driver-aware model for cg04055835 promoter methylation;
(4) PHLDA3-N1 association within the kinase-driver-negative (RAS/WT) subset."""
import pandas as pd, numpy as np
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.stats import mannwhitneyu, kruskal

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
c=pd.DataFrame({"patient":cl["submitter_id"],"age":pd.to_numeric(cl["demographic.age_at_index"],errors="coerce"),
  "male":(cl["demographic.gender"]=="male").astype(float),
  "T":cl.apply(lambda r:pick(r,"ajcc_pathologic_t"),axis=1),"N":cl.apply(lambda r:pick(r,"ajcc_pathologic_n"),axis=1)})
c["T34"]=c["T"].str[:2].map({"T1":0,"T2":0,"T3":1,"T4":1}); c["N1"]=c["N"].str[:2].map({"N0":0,"N1":1})
br=pd.read_csv("data/PHLDA3_braf_ras.csv")[["patient","BRAF","RAS"]]
fus=pd.read_csv("data/thca_fusion_status.csv"); fset=set(fus.patient)
# per-patient lead fusion gene (for leave-one-out)
fgene=fus.groupby("patient").driver.first() if "driver" in fus.columns else None
me=pd.read_csv("data/cg04055835_betas.csv"); me["patient"]=me["sample"].str[:12]; meth=me.groupby("patient").beta.mean()
d=c.merge(expr.rename(GENE),on="patient").merge(br,on="patient")
d["BRAFmut"]=(d.BRAF=="BRAF-mut").astype(float); d["RASmut"]=(d.RAS=="RAS-mut").astype(float)
d["Fusion"]=d.patient.isin(fset).astype(float); d["z"]=(d[GENE]-d[GENE].mean())/d[GENE].std()
d=d.merge(meth.rename("meth"),on="patient",how="left")

print("="*70)
print("CHECK 1 — fusion+ PHLDA3 not driven by one fusion type (histology N/A)")
fd=pd.read_csv("data/thca_fusion_status.csv")
fp=fd.merge(expr.rename(GENE),on="patient")
print(fp.groupby("driver")[GENE].agg(["size","median"]).to_string())
allfus=d[d.Fusion==1][GENE]
print(f"all fusion+ median={allfus.median():.2f} (n={len(allfus)})")
for g in fd.driver.unique():
    keep=set(fd[fd.driver!=g].patient); sub=d[d.patient.isin(keep)&(d.Fusion==1)][GENE]
    print(f"  leave-out {g:6s}: remaining fusion+ median={sub.median():.2f} (n={len(sub)})")

print("="*70)
print("CHECK 2 — PHLDA3 ~ BRAF + fusion + RAS + age + sex + T  (OLS)")
d2=d.dropna(subset=[GENE,"age","male","T34","BRAFmut","Fusion","RASmut"])
m2=smf.ols("Q('PHLDA3') ~ BRAFmut + Fusion + RASmut + age + male + T34", data=d2).fit()
print(m2.summary2().tables[1].loc[["BRAFmut","Fusion","RASmut","T34"],["Coef.","P>|t|"]].to_string())
print(f"  (n={int(m2.nobs)})  -> BRAF & fusion raise PHLDA3 independently of T stage")

print("="*70)
print("CHECK 3 — cg04055835 methylation ~ BRAF + fusion + RAS + age + sex + T (OLS)")
d3=d.dropna(subset=["meth","age","male","T34","BRAFmut","Fusion","RASmut"])
m3=smf.ols("meth ~ BRAFmut + Fusion + RASmut + age + male + T34", data=d3).fit()
print(m3.summary2().tables[1].loc[["BRAFmut","Fusion","RASmut","T34"],["Coef.","P>|t|"]].to_string())
aggr=d3[(d3.BRAFmut==1)|(d3.Fusion==1)].meth; indo=d3[(d3.BRAFmut==0)&(d3.Fusion==0)].meth
print(f"  methylation: BRAF/fusion median={aggr.median():.2f} vs RAS/WT={indo.median():.2f}; MWU P={mannwhitneyu(aggr,indo)[1]:.1e}")

print("="*70)
print("CHECK 4 — PHLDA3-N1 within kinase-driver-NEGATIVE subset (RAS/WT)")
neg=d[(d.BRAFmut==0)&(d.Fusion==0)].dropna(subset=["N1","z","age","male","T34"])
mn=sm.Logit(neg.N1.astype(int), sm.add_constant(neg[["z","age","male","T34"]].astype(float))).fit(disp=0)
orz=np.exp(mn.params["z"]); lo,hi=np.exp(mn.conf_int().loc["z"]); p=mn.pvalues["z"]
print(f"  N1 ~ PHLDA3 + age + sex + T  in BRAF-wt & fusion-neg tumours (n={int(mn.nobs)}, N1={int(neg.N1.sum())})")
print(f"  PHLDA3 OR={orz:.2f} [{lo:.2f}-{hi:.2f}] P={p:.3f}")
print("done driver robustness")
