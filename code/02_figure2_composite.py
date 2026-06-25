#!/usr/bin/env python3
"""Composite main Figure 2 (N1 prediction) — Nature-style multi-panel:
(a) N1 multivariable forest, (b) BRAF sensitivity, (c) BRAF-stratified subgroups,
(d) nomogram, (e) ROC, (f) calibration, (g) decision curve. Forests from cached
tables; nomogram/ROC/calibration/DCA recomputed from the N1 model."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
from sklearn.model_selection import cross_val_predict
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.calibration import calibration_curve
import nature_style as ns
ns.set_style(font_size=7)

GENE, ENS = "PHLDA3", "ENSG00000174307"
# ---- data + N1 model (PHLDA3 + age + T34) for nomogram/ROC/cal/DCA ----
df = pd.read_csv("data/TCGA-THCA.htseq_counts.tsv", sep="\t", index_col=0)
mat=(2**df-1).clip(lower=0); row=mat.index[mat.index.str.startswith(ENS)][0]
lc=np.log2(mat.loc[row]/mat.sum(axis=0)*1e6+1)
st=pd.Series([b.split("-")[3][:2] for b in df.columns],index=df.columns)
expr=pd.DataFrame({"patient":[b[:12] for b in lc[st=="01"].index],GENE:lc[st=="01"].values}).groupby("patient")[GENE].mean()
cl=pd.read_csv("data/THCA_clinical.tsv",sep="\t",dtype=str)
def pick(r,bb):
    for k in [f"diagnoses.0.{bb}",f"diagnoses.1.{bb}"]:
        if k in r and pd.notna(r[k]) and str(r[k]) not in("","nan"): return r[k]
    return np.nan
c=pd.DataFrame({"patient":cl["submitter_id"],"age":pd.to_numeric(cl["demographic.age_at_index"],errors="coerce"),
  "T":cl.apply(lambda r:pick(r,"ajcc_pathologic_t"),axis=1),"N":cl.apply(lambda r:pick(r,"ajcc_pathologic_n"),axis=1)})
d=c.merge(expr.rename(GENE),on="patient"); d["T34"]=d["T"].str[:2].map({"T1":0,"T2":0,"T3":1,"T4":1})
d["N1"]=d["N"].str[:2].map({"N0":0,"N1":1}); d=d.dropna(subset=[GENE,"age","T34","N1"]).reset_index(drop=True)
PRED=[GENE,"age","T34"]; X,y=d[PRED].astype(float),d["N1"].astype(int)
m=sm.Logit(y,sm.add_constant(X)).fit(disp=0); b=m.params
cv=cross_val_predict(LogisticRegression(max_iter=1000),X.values,y.values,cv=10,method="predict_proba")[:,1]
auc_app=roc_auc_score(y,m.predict(sm.add_constant(X))); auc_cv=roc_auc_score(y,cv)

fig=plt.figure(figsize=(7.2,8.6))
gs=fig.add_gridspec(3,3,height_ratios=[0.85,1.05,0.9],hspace=0.65,wspace=0.5)
axA=fig.add_subplot(gs[0,0]); axB=fig.add_subplot(gs[0,1]); axC=fig.add_subplot(gs[0,2])
axN=fig.add_subplot(gs[1,:]); axR=fig.add_subplot(gs[2,0]); axCa=fig.add_subplot(gs[2,1]); axD=fig.add_subplot(gs[2,2])

def forest(ax,rows,title,hero=GENE):
    rows=rows[::-1]
    for i,r in enumerate(rows):
        h=r["lab"]==hero or "PHLDA3" in str(r["lab"])
        col=ns.C["tumor"] if h else ns.C["normal"]
        ax.plot([r["lo"],r["hi"]],[i,i],color=col,lw=1.6); ax.scatter(r["OR"],i,s=34 if h else 22,color=col,ec="black",lw=.4,zorder=3)
    ax.axvline(1,ls="--",color=ns.C["grey_mid"],lw=.9); ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([r["lab"] for r in rows],fontsize=6); ax.set_xscale("log")
    ax.set_xlabel("OR (95% CI)",fontsize=6.5); ax.set_title(title,fontweight="bold",fontsize=7); ax.margins(x=.3)
    ax.spines[["top","right"]].set_visible(False); ax.tick_params(labelsize=6)

# (a) N1 multivariable forest (cache)
lt=pd.read_csv("data/PHLDA3_logistic_table.csv")
PR={"PHLDA3_z":"PHLDA3","age":"Age","age_c":"Age","male":"Male","T34":"T3/T4"}
nv=lt[(lt.outcome=="N1")&(lt.model=="Multivariable")]
forest(axA,[{"lab":PR.get(r["var"],r["var"]),"OR":r.OR,"lo":r.lo,"hi":r.hi} for _,r in nv.iterrows()],"N1 multivariable")
# (b) BRAF sensitivity (cache)
bs=pd.read_csv("data/PHLDA3_braf_sensitivity.csv")
forest(axB,[{"lab":"PHLDA3 +BRAF" if "BRAF" in r["model"] else "PHLDA3 base","OR":r.OR,"lo":r.lo,"hi":r.hi} for _,r in bs.iterrows()],"BRAF sensitivity")
# (c) BRAF-stratified (cache)
bst=pd.read_csv("data/PHLDA3_braf_stratified.csv")
forest(axC,[{"lab":r["group"],"OR":r.OR,"lo":r.lo,"hi":r.hi} for _,r in bst.iterrows()],"BRAF-stratified")

# (d) nomogram
LB={GENE:f"{GENE} log2(CPM+1)","age":"Age (years)","T34":"T stage"}
rng={v:(X[v].min(),X[v].max()) for v in PRED}; gg=lambda v,x:b[v]*x
ref={v:min(gg(v,rng[v][0]),gg(v,rng[v][1])) for v in PRED}; span={v:abs(b[v])*(rng[v][1]-rng[v][0]) for v in PRED}
maxc=max(span.values()); pts=lambda v,x:100*(gg(v,x)-ref[v])/maxc
order=sorted(PRED,key=lambda v:-span[v]); yrow={v:i for i,v in enumerate(["points"]+order+["total","lp"])}
ytop=len(yrow)-1; Y=lambda nm:ytop-yrow[nm]; LX=-4
axN.plot([0,100],[Y("points")]*2,"k-",lw=.9)
for t in range(0,101,20): axN.plot([t,t],[Y("points")-.08,Y("points")+.08],"k-",lw=.7); axN.text(t,Y("points")+.22,str(t),ha="center",fontsize=5.5)
axN.text(LX,Y("points"),"Points",ha="right",va="center",fontweight="bold",fontsize=6.5)
for v in order:
    lo,hi=rng[v]
    if v=="T34": tk=[0,1]; lb=["T1/2","T3/4"]
    elif v=="age": tk=np.arange(20,hi+1,15); lb=[f"{t:.0f}" for t in tk]
    else: tk=np.linspace(lo,hi,5); lb=[f"{t:.1f}" for t in tk]
    pp=[pts(v,t) for t in tk]; axN.plot([min(pp),max(pp)],[Y(v)]*2,"k-",lw=.9)
    for t,q,l in zip(tk,pp,lb): axN.plot([q,q],[Y(v)-.07,Y(v)+.07],"k-",lw=.7); axN.text(q,Y(v)-.26,l,ha="center",fontsize=5.5)
    axN.text(LX,Y(v),LB[v],ha="right",va="center",fontsize=6.5)
tot=sum(max(pts(v,rng[v][0]),pts(v,rng[v][1])) for v in PRED); axN.plot([0,tot],[Y("total")]*2,"k-",lw=.9)
for t in np.linspace(0,tot,6): axN.plot([t,t],[Y("total")-.08,Y("total")+.08],"k-",lw=.7); axN.text(t,Y("total")+.22,f"{t:.0f}",ha="center",fontsize=5.5)
axN.text(LX,Y("total"),"Total Points",ha="right",va="center",fontweight="bold",fontsize=6.5)
base=b["const"]+sum(ref[v] for v in PRED); p2t=lambda p:(100/maxc)*(np.log(p/(1-p))-base)
vis=[(p,p2t(p)) for p in [.1,.2,.3,.4,.5,.6,.7,.8] if -2<=p2t(p)<=tot+5]
axN.plot([vis[0][1],vis[-1][1]],[Y("lp")]*2,"k-",lw=.9)
for p,xp in vis: axN.plot([xp,xp],[Y("lp")-.08,Y("lp")+.08],"k-",lw=.7); axN.text(xp,Y("lp")-.26,f"{p:.1f}",ha="center",fontsize=5.5)
axN.text(LX,Y("lp"),"Pr(N1)",ha="right",va="center",fontweight="bold",fontsize=6.5)
axN.set_xlim(-26,106); axN.set_ylim(-.7,ytop+.8); axN.axis("off")
axN.set_title("Nomogram for predicting lymph-node metastasis (N1)",fontweight="bold",fontsize=7.5)

# (e) ROC
fpr,tpr,_=roc_curve(y,m.predict(sm.add_constant(X))); fc,tc,_=roc_curve(y,cv)
axR.plot(fpr,tpr,color=ns.C["tumor"],lw=1.6,label=f"Apparent {auc_app:.2f}")
axR.plot(fc,tc,color=ns.C["normal"],lw=1.3,label=f"10-fold CV {auc_cv:.2f}")
axR.plot([0,1],[0,1],"--",color=ns.C["grey_mid"],lw=.8); axR.set_xlabel("1 - Specificity",fontsize=6.5); axR.set_ylabel("Sensitivity",fontsize=6.5)
axR.set_title("ROC",fontweight="bold",fontsize=7.5); axR.legend(frameon=False,fontsize=5.5,loc="lower right"); axR.spines[["top","right"]].set_visible(False)
# (f) calibration
pt,pp=calibration_curve(y,cv,n_bins=8,strategy="quantile")
axCa.plot([0,1],[0,1],"--",color=ns.C["grey_mid"],lw=.8); axCa.plot(pp,pt,"o-",color=ns.C["tumor"],lw=1.4,ms=3)
axCa.set_xlabel("Predicted",fontsize=6.5); axCa.set_ylabel("Observed",fontsize=6.5); axCa.set_title("Calibration",fontweight="bold",fontsize=7.5)
axCa.set_xlim(0,1); axCa.set_ylim(0,1); axCa.spines[["top","right"]].set_visible(False)
# (g) DCA
def nb(yt,p,th): pr=p>=th; n=len(yt); return np.sum(pr&(yt==1))/n - np.sum(pr&(yt==0))/n*(th/(1-th))
ths=np.linspace(.01,.8,80); nbm=[nb(y.values,cv,t) for t in ths]; prev=y.mean(); nba=[prev-(1-prev)*(t/(1-t)) for t in ths]
axD.plot(ths,nbm,color=ns.C["tumor"],lw=1.6,label="Nomogram"); axD.plot(ths,nba,color=ns.C["grey_mid"],lw=1,label="Treat all")
axD.axhline(0,color="black",lw=.9,label="Treat none"); axD.set_ylim(-.05,max(nbm)*1.15)
axD.set_xlabel("Threshold",fontsize=6.5); axD.set_ylabel("Net benefit",fontsize=6.5); axD.set_title("Decision curve",fontweight="bold",fontsize=7.5)
axD.legend(frameon=False,fontsize=5.5,loc="upper right"); axD.spines[["top","right"]].set_visible(False)

for ax,l in zip([axA,axB,axC,axN,axR,axCa,axD],"abcdefg"): ns.panel_label(ax,l,x=-0.18,y=1.05,fs=10)
fig.suptitle("Figure 2. PHLDA3 independently predicts lymph-node metastasis, including after BRAF adjustment",
             fontweight="bold",fontsize=8.5,y=1.005)
ns.save_fig(fig,"PHLDA3_THCA_Figure2")
print(f"AUC apparent={auc_app:.2f} CV={auc_cv:.2f}; done Figure 2 composite")
