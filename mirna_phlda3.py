#!/usr/bin/env python3
"""Task 5 — PHLDA3-targeting miRNAs (ENCORI/AGO2-CLIP) vs PHLDA3 expression in
TCGA-THCA. CLIP-supported targeting + inverse correlation = candidate repressor
miRNAs. Output also feeds a lncRNA-miRNA-PHLDA3 ceRNA hypothesis."""
import pandas as pd, numpy as np, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import nature_style as ns
ns.set_style(font_size=8)

GENE, ENS = "PHLDA3", "ENSG00000174307"

# ---------- ENCORI: CLIP-supported miRNAs targeting PHLDA3 ----------
enc = pd.read_csv("encori_phlda3_mirna.txt", sep="\t", comment="#")
enc = enc.groupby("miRNAname").agg(clip=("clipExpNum","max")).reset_index()
def to_precursor(m):  # hsa-miR-15a-5p -> hsa-mir-15a
    b = re.sub(r"-(3p|5p)$","", m); return b.replace("miR","mir")
enc["precursor"] = enc["miRNAname"].map(to_precursor)
print(f"ENCORI PHLDA3-targeting miRNAs: {len(enc)}")

# ---------- TCGA miRNA expression ----------
mi = pd.read_csv("TCGA-THCA.mirna.tsv.gz", sep="\t", index_col=0)
mi = mi[[c for c in mi.columns if c.split("-")[3][:2]=="01"]]      # tumour
mi.columns = [c[:12] for c in mi.columns]
mi = mi.groupby(level=0, axis=1).mean()

# ---------- PHLDA3 expression per patient ----------
df = pd.read_csv("TCGA-THCA.htseq_counts.tsv", sep="\t", index_col=0)
mat=(2**df-1).clip(lower=0); row=mat.index[mat.index.str.startswith(ENS)][0]
logcpm=np.log2(mat.loc[row]/mat.sum(axis=0)*1e6+1)
st=pd.Series([b.split("-")[3][:2] for b in df.columns],index=df.columns)
expr=pd.DataFrame({"patient":[b[:12] for b in logcpm[st=="01"].index],GENE:logcpm[st=="01"].values}).groupby("patient")[GENE].mean()

# ---------- match precursor -> TCGA miRNA rows, correlate ----------
rows=[]
for _,r in enc.iterrows():
    hits=[ix for ix in mi.index if ix==r.precursor or ix.startswith(r.precursor+"-")]
    if not hits: continue
    mir_expr = mi.loc[hits].mean(axis=0)
    d = pd.DataFrame({"mir":mir_expr}).join(expr.rename("phl"), how="inner").dropna()
    if len(d)<30: continue
    rho,p = spearmanr(d.mir, d.phl)
    rows.append({"miRNA":r.miRNAname,"clip":r["clip"],"rho":rho,"p":p,"n":len(d)})
res=pd.DataFrame(rows).sort_values("rho")
res["sig"]=res["p"]<0.05
res.to_csv("PHLDA3_mirna_corr.csv", index=False)
neg = res[(res.rho<0)&(res.sig)]
print(f"miRNAs matched & tested: {len(res)} | CLIP-supported negative (candidate repressors): {len(neg)}")
print(neg.sort_values('rho').head(10)[['miRNA','clip','rho','p']].round(3).to_string(index=False))

# ---------- figure: top candidate repressor miRNAs ----------
top = res[res.sig].sort_values("rho").head(18)[::-1]
cols=[ns.C["normal"] if v<0 else ns.C["tumor"] for v in top["rho"]]
fig, ax = plt.subplots(figsize=(6.6,5.2))
ax.barh(range(len(top)), top["rho"], color=cols, edgecolor="black", lw=0.4)
ax.set_yticks(range(len(top)))
ax.set_yticklabels([f"{m}  (CLIP {int(c)})" for m,c in zip(top["miRNA"],top["clip"])], fontsize=7)
ax.axvline(0,color="black",lw=0.8)
ax.set_xlabel(f"Spearman r with {GENE} (TCGA-THCA tumours)")
ax.set_title("AGO2-CLIP-supported miRNAs targeting PHLDA3\n(blue = inverse → candidate repressors)",
             fontweight="bold", fontsize=8.5)
ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_mirna")
print("done task 5 (miRNA core)")
