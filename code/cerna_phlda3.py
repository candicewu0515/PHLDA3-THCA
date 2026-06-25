#!/usr/bin/env python3
"""Task 5 (ceRNA) — lncRNA-miRNA-PHLDA3 competing-endogenous-RNA network.
lncRNAs that (i) sponge PHLDA3-repressing miRNAs (ENCORI/AGO2-CLIP) and
(ii) positively correlate with PHLDA3 -> candidate ceRNA regulators."""
import pandas as pd, numpy as np, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nature_style as ns
ns.set_style(font_size=8)

# lead repressor miRNAs (CLIP-supported, inverse with PHLDA3)
mcorr = pd.read_csv("data/PHLDA3_mirna_corr.csv")
LEAD = ["hsa-miR-7-5p","hsa-miR-195-5p","hsa-miR-424-5p","hsa-miR-497-5p","hsa-miR-103a-3p"]
mrho = dict(zip(mcorr.miRNA, mcorr.rho))

# ENCORI lncRNA-miRNA interactions
rawcols=None; pairs=[]
for ln in open("data/encori_lncrna_raw.txt"):
    if ln.startswith("#"): continue
    f=ln.rstrip("\n").split("\t")
    if f[0]=="miRNAid": rawcols=f; continue
    if rawcols and len(f)>=4:
        pairs.append({"miRNA":f[1], "lncRNA":re.sub(r"\.\d+$","",f[3])})
pp=pd.DataFrame(pairs).drop_duplicates()
pp=pp[pp.miRNA.isin(LEAD)]

# PHLDA3 co-expression (positively correlated lncRNAs = ceRNA candidates)
co=pd.read_csv("data/PHLDA3_coexpression.csv"); co.columns=["gene","r"]
co["gene"]=co["gene"].astype(str).str.replace(r"\.\d+$","",regex=True)
cor=dict(zip(co.gene, co.r))
pp["lnc_r"]=pp["lncRNA"].map(cor)
cand=pp[pp.lnc_r>0.2].dropna()       # positive correlation with PHLDA3
# rank lncRNAs by correlation, keep top
top_lnc=(cand.groupby("lncRNA")["lnc_r"].max().sort_values(ascending=False).head(8).index.tolist())
net=cand[cand.lncRNA.isin(top_lnc)]
net.to_csv("data/PHLDA3_cerna_network.csv", index=False)
print(f"ceRNA lncRNAs (sponge lead miRNA + r>0.2 with PHLDA3): {len(top_lnc)}")
print(net[['lncRNA','miRNA','lnc_r']].drop_duplicates().to_string(index=False))

# ================= 3-layer network figure =================
mids=sorted(net.miRNA.unique(), key=lambda m: mrho.get(m,0))
lncs=sorted(top_lnc, key=lambda l: cand[cand.lncRNA==l].lnc_r.max())
fig, ax = plt.subplots(figsize=(8.5, 6))
xL, xM, xR = 0.0, 1.0, 2.0
yL={l:i*(len(mids)/max(len(lncs),1)) for i,l in enumerate(lncs)}
yM={m:i for i,m in enumerate(mids)}
yPH=(len(mids)-1)/2
# edges lncRNA-miRNA
for _,e in net.iterrows():
    ax.plot([xL,xM],[yL[e.lncRNA],yM[e.miRNA]], color=ns.C["grey_light"], lw=0.8, zorder=1)
# edges miRNA-PHLDA3 (repression)
for m in mids:
    ax.plot([xM,xR],[yM[m],yPH], color=ns.C["normal"], lw=1.4, ls=(0,(4,2)), zorder=1)
# nodes
for l,y in yL.items():
    ax.scatter(xL,y,s=420,color=ns.C["tumor"],alpha=.85,edgecolors="black",lw=.5,zorder=3)
    ax.text(xL-0.06,y,l,ha="right",va="center",fontsize=7)
for m,y in yM.items():
    ax.scatter(xM,y,s=300,color=ns.C["accent"],alpha=.9,edgecolors="black",lw=.5,zorder=3)
    ax.text(xM,y+0.16,m.replace("hsa-",""),ha="center",va="bottom",fontsize=7,
            color=ns.C["grey_dark"],fontweight="bold")
ax.scatter(xR,yPH,s=1500,color=ns.C["highlight"],edgecolors="black",lw=.8,zorder=3)
ax.text(xR,yPH,"PHLDA3",ha="center",va="center",fontsize=9,fontweight="bold")
ax.text(xL,len(mids)+0.2,"lncRNA (ceRNA)\n+corr with PHLDA3",ha="center",fontsize=7.5,fontweight="bold",color=ns.C["tumor"])
ax.text(xM,len(mids)+0.2,"miRNA\n(AGO2-CLIP, repressor)",ha="center",fontsize=7.5,fontweight="bold",color=ns.C["accent"])
ax.text(xR,len(mids)+0.2,"target",ha="center",fontsize=7.5,fontweight="bold",color=ns.C["highlight"])
ax.set_xlim(-0.7,2.5); ax.set_ylim(-0.7,len(mids)+0.7); ax.axis("off")
ax.set_title("PHLDA3 ceRNA network: lncRNA — miRNA — PHLDA3\n(TCGA-THCA + ENCORI AGO2-CLIP)",
             fontweight="bold", fontsize=10)
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_cerna")
print("done task 5 (ceRNA)")
