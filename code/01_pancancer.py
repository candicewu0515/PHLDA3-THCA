#!/usr/bin/env python3
"""Pan-cancer PHLDA3 expression (TCGA tumor vs TCGA+GTEx normal, UCSC Xena toil).
GEPIA-style tumor-vs-normal across cancers, THCA highlighted."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu
import nature_style as ns
ns.set_style(font_size=8)

ABBR = {  # TCGA primary disease -> abbreviation
 "thyroid carcinoma":"THCA","breast invasive carcinoma":"BRCA","lung adenocarcinoma":"LUAD",
 "lung squamous cell carcinoma":"LUSC","kidney clear cell carcinoma":"KIRC",
 "kidney papillary cell carcinoma":"KIRP","kidney chromophobe":"KICH","liver hepatocellular carcinoma":"LIHC",
 "stomach adenocarcinoma":"STAD","colon adenocarcinoma":"COAD","rectum adenocarcinoma":"READ",
 "bladder urothelial carcinoma":"BLCA","prostate adenocarcinoma":"PRAD","uterine corpus endometrioid carcinoma":"UCEC",
 "head & neck squamous cell carcinoma":"HNSC","esophageal carcinoma":"ESCA","pancreatic adenocarcinoma":"PAAD",
 "cervical & endocervical cancer":"CESC","glioblastoma multiforme":"GBM","brain lower grade glioma":"LGG",
 "skin cutaneous melanoma":"SKCM","ovarian serous cystadenocarcinoma":"OV","testicular germ cell tumor":"TGCT",
 "cholangiocarcinoma":"CHOL","adrenocortical cancer":"ACC","pheochromocytoma & paraganglioma":"PCPG",
 "uterine carcinosarcoma":"UCS","sarcoma":"SARC","mesothelioma":"MESO","thymoma":"THYM",
 "diffuse large b-cell lymphoma":"DLBC","acute myeloid leukemia":"LAML","uveal melanoma":"UVM"}

df = pd.read_csv("data/pancancer_phlda3_raw.csv")
df["PHLDA3"] = pd.to_numeric(df["PHLDA3"], errors="coerce")
df = df.dropna(subset=["PHLDA3"])
df["d"] = df["disease"].str.lower().str.strip()
tum_all = df[(df.study=="TCGA") & (df.stype=="Primary Tumor")]
norm_all = df[df.stype.isin(["Normal Tissue","Solid Tissue Normal"])]   # GTEx + TCGA normal

rows=[]
for d, sub in tum_all.groupby("d"):
    site = sub["site"].mode().iloc[0] if not sub["site"].mode().empty else None
    nrm = norm_all[norm_all.site==site]
    if len(nrm) < 5 or d not in ABBR: continue
    t, n = sub.PHLDA3.values, nrm.PHLDA3.values
    p = mannwhitneyu(t, n, alternative="two-sided")[1]
    rows.append({"abbr":ABBR[d], "tumor":t, "normal":n, "nt":len(t), "nn":len(n),
                 "fc":np.median(t)-np.median(n), "p":p})
res = pd.DataFrame(rows).sort_values("fc", ascending=False).reset_index(drop=True)
pd.DataFrame([{k:r[k] for k in ["abbr","nt","nn","fc","p"]} for _,r in res.iterrows()]
            ).to_csv("data/PHLDA3_pancancer_stats.csv", index=False)
print(f"cancers plotted: {len(res)}")
print(res[["abbr","nt","nn","fc","p"]].to_string(index=False))

# ---------- figure: grouped tumor/normal boxplot ----------
fig, ax = plt.subplots(figsize=(17, 5.6))
pos=0; xticks=[]; labels=[]
for _, r in res.iterrows():
    is_thca = r["abbr"]=="THCA"
    ct = ns.C["tumor"] if not is_thca else ns.C["highlight"]
    cn = ns.C["normal"]
    bt=ax.boxplot(r["tumor"], positions=[pos], widths=0.62, patch_artist=True, showfliers=False,
                  medianprops=dict(color="black",lw=1))
    bn=ax.boxplot(r["normal"], positions=[pos+0.7], widths=0.62, patch_artist=True, showfliers=False,
                  medianprops=dict(color="black",lw=1))
    bt["boxes"][0].set_facecolor(ct); bt["boxes"][0].set_alpha(0.85 if is_thca else 0.6)
    bn["boxes"][0].set_facecolor(cn); bn["boxes"][0].set_alpha(0.6)
    sig = "***" if r.p<1e-3 else "**" if r.p<1e-2 else "*" if r.p<0.05 else ""
    ymax=max(np.percentile(r["tumor"],98), np.percentile(r["normal"],98))
    if sig: ax.text(pos+0.35, ymax+0.2, sig, ha="center", fontsize=8)
    xticks.append(pos+0.35)
    labels.append(r["abbr"]); pos+=2
# highlight THCA tick
ax.set_xticks(xticks)
ax.set_xticklabels(labels, rotation=90, fontsize=9)
for lbl in ax.get_xticklabels():
    if lbl.get_text()=="THCA": lbl.set_color(ns.C["highlight"]); lbl.set_fontweight("bold")
ax.set_ylabel("PHLDA3 expression  log2(norm_count+1)")
ax.set_title("Pan-cancer PHLDA3: tumor vs normal (TCGA tumor vs TCGA+GTEx normal)  —  THCA highlighted",
             fontweight="bold", fontsize=12)
from matplotlib.patches import Patch
ax.legend(handles=[Patch(fc=ns.C["tumor"],alpha=0.6,label="Tumor"),Patch(fc=ns.C["normal"],alpha=0.6,label="Normal")],
          loc="upper right", frameon=False)
ax.spines[["top","right"]].set_visible(False)
ax.margins(x=0.01)
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_pancancer")
print("+ PHLDA3_pancancer_stats.csv")
