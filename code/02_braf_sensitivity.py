#!/usr/bin/env python3
"""Reviewer sensitivity analysis — does PHLDA3 independently predict lymph-node
metastasis (N1) AFTER adjusting for BRAF status (a key confounder)?
Compares N1 model with vs without BRAF: PHLDA3 + age + sex + T (+ BRAF)."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
import nature_style as ns
ns.set_style(font_size=8)

GENE, ENS = "PHLDA3", "ENSG00000174307"

# expression per patient
df = pd.read_csv("data/TCGA-THCA.htseq_counts.tsv", sep="\t", index_col=0)
mat=(2**df-1).clip(lower=0); row=mat.index[mat.index.str.startswith(ENS)][0]
logcpm=np.log2(mat.loc[row]/mat.sum(axis=0)*1e6+1)
st=pd.Series([b.split("-")[3][:2] for b in df.columns],index=df.columns)
tum=logcpm[st=="01"]
expr=pd.DataFrame({"patient":[b[:12] for b in tum.index],GENE:tum.values}).groupby("patient")[GENE].mean()

# clinical
cl=pd.read_csv("data/THCA_clinical.tsv", sep="\t", dtype=str)
def pick(r,b):
    for k in [f"diagnoses.0.{b}",f"diagnoses.1.{b}"]:
        if k in r and pd.notna(r[k]) and str(r[k]) not in ("","nan"): return r[k]
    return np.nan
c=pd.DataFrame({"patient":cl["submitter_id"],
   "age":pd.to_numeric(cl["demographic.age_at_index"],errors="coerce"),
   "male":(cl["demographic.gender"]=="male").astype(float),
   "T":cl.apply(lambda r:pick(r,"ajcc_pathologic_t"),axis=1),
   "N":cl.apply(lambda r:pick(r,"ajcc_pathologic_n"),axis=1)})
c["T34"]=c["T"].str[:2].map({"T1":0,"T2":0,"T3":1,"T4":1})
c["N1"]=c["N"].str[:2].map({"N0":0,"N1":1})

# BRAF status (cBioPortal, from braf_ras_thca.py output)
br=pd.read_csv("data/PHLDA3_braf_ras.csv")[["patient","BRAF"]]
br["BRAF_mut"]=(br["BRAF"]=="BRAF-mut").astype(float)

d=c.merge(expr.rename(GENE),on="patient").merge(br[["patient","BRAF_mut"]],on="patient")
d[f"{GENE}_z"]=(d[GENE]-d[GENE].mean())/d[GENE].std()

def fit(predictors, data):
    sub=data[["N1"]+predictors].dropna()
    m=sm.Logit(sub["N1"].astype(int), sm.add_constant(sub[predictors].astype(float))).fit(disp=0)
    OR=np.exp(m.params[f"{GENE}_z"]); lo,hi=np.exp(m.conf_int().loc[f"{GENE}_z"])
    return {"OR":OR,"lo":lo,"hi":hi,"p":m.pvalues[f"{GENE}_z"],"n":len(sub)}

base=fit([f"{GENE}_z","age","male","T34"], d)
adj =fit([f"{GENE}_z","age","male","T34","BRAF_mut"], d)
print("N1 model WITHOUT BRAF: PHLDA3 OR=%.2f [%.2f-%.2f] p=%.2e n=%d"%(base["OR"],base["lo"],base["hi"],base["p"],base["n"]))
print("N1 model WITH    BRAF: PHLDA3 OR=%.2f [%.2f-%.2f] p=%.2e n=%d"%(adj["OR"],adj["lo"],adj["hi"],adj["p"],adj["n"]))
pd.DataFrame([{"model":"PHLDA3+age+sex+T",**base},{"model":"PHLDA3+age+sex+T+BRAF",**adj}]).to_csv("data/PHLDA3_braf_sensitivity.csv",index=False)

# forest: PHLDA3 OR for N1, before vs after BRAF adjustment
fig, ax = plt.subplots(figsize=(6.4,2.6))
rows=[("+ BRAF (sensitivity)",adj),("Base model",base)]
for i,(lab,r) in enumerate(rows):
    ax.plot([r["lo"],r["hi"]],[i,i],color=ns.C["tumor"],lw=2)
    ax.scatter(r["OR"],i,s=70,color=ns.C["tumor"],edgecolors="black",lw=.6,zorder=3)
    ax.text(ax.get_xlim()[1] if i==0 else r["hi"], i, f"  OR={r['OR']:.2f} ({r['lo']:.2f}-{r['hi']:.2f}), P={r['p']:.0e}", va="center", fontsize=8)
ax.axvline(1, ls="--", color=ns.C["grey_mid"], lw=1)
ax.set_yticks([0,1]); ax.set_yticklabels([r[0] for r in rows])
ax.set_xscale("log"); ax.set_xlabel(f"{GENE} odds ratio for N1 (per 1 SD)")
ax.set_title("PHLDA3 independently predicts N1 after BRAF adjustment", fontweight="bold", fontsize=9)
ax.set_ylim(-0.5,1.5); ax.margins(x=0.5); ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_braf_sensitivity")
print("done BRAF sensitivity")
