#!/usr/bin/env python3
"""PHLDA3-based clinical prediction model for lymph-node metastasis (N1):
nomogram + ROC + calibration curve + decision-curve analysis (DCA).
Predictors: PHLDA3 expression, age, T stage. CV used for honest calibration/DCA."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
from sklearn.model_selection import cross_val_predict
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.calibration import calibration_curve
import nature_style as ns

ns.set_style(font_size=8)

GENE, ENS = "PHLDA3", "ENSG00000174307"

# ---------- merge expression + clinical (same pipeline as logistic) ----------
df = pd.read_csv("data/TCGA-THCA.htseq_counts.tsv", sep="\t", index_col=0)
mat = (2**df - 1).clip(lower=0)
row = mat.index[mat.index.str.startswith(ENS)][0]
logcpm = np.log2(mat.loc[row] / mat.sum(axis=0) * 1e6 + 1)
stype = pd.Series([b.split("-")[3][:2] for b in df.columns], index=df.columns)
tum = logcpm[stype == "01"]
expr = pd.DataFrame({"patient":[b[:12] for b in tum.index], GENE:tum.values}).groupby("patient")[GENE].mean()
cl = pd.read_csv("data/THCA_clinical.tsv", sep="\t", dtype=str)
def pick(r, base):
    for k in [f"diagnoses.0.{base}", f"diagnoses.1.{base}"]:
        if k in r and pd.notna(r[k]) and str(r[k]) not in ("","nan"): return r[k]
    return np.nan
c = pd.DataFrame({"patient":cl["submitter_id"],
    "age":pd.to_numeric(cl["demographic.age_at_index"],errors="coerce"),
    "T":cl.apply(lambda r:pick(r,"ajcc_pathologic_t"),axis=1),
    "N":cl.apply(lambda r:pick(r,"ajcc_pathologic_n"),axis=1)})
d = c.merge(expr.rename(GENE), on="patient", how="inner")
d["T34"] = d["T"].str[:2].map({"T1":0,"T2":0,"T3":1,"T4":1})
d["N1"]  = d["N"].str[:2].map({"N0":0,"N1":1})
d = d.dropna(subset=[GENE,"age","T34","N1"]).reset_index(drop=True)
PRED = [GENE, "age", "T34"]
X, y = d[PRED].astype(float), d["N1"].astype(int)
print(f"model n={len(d)}, N1 events={y.sum()}")

# ---------- fit (statsmodels for clean coefficients) ----------
Xc = sm.add_constant(X)
m = sm.Logit(y, Xc).fit(disp=0)
b = m.params
print(m.summary2().tables[1][["Coef.","P>|z|"]].round(4).to_string())

# ---------- honest CV predictions for AUC / calibration / DCA ----------
clf = LogisticRegression(max_iter=1000)
cv_prob = cross_val_predict(clf, X.values, y.values, cv=10, method="predict_proba")[:,1]
auc_cv = roc_auc_score(y, cv_prob)
auc_app = roc_auc_score(y, m.predict(Xc))
print(f"AUC apparent={auc_app:.3f}  10-fold CV={auc_cv:.3f}")

# ================= figure: nomogram | ROC | calibration | DCA =================
fig = plt.figure(figsize=(15, 9))
gs = fig.add_gridspec(2, 3, height_ratios=[1.05, 1])
axN = fig.add_subplot(gs[0, :])          # nomogram spans top
axR = fig.add_subplot(gs[1, 0])
axC = fig.add_subplot(gs[1, 1])
axD = fig.add_subplot(gs[1, 2])

# ---- nomogram (points always 0..100, lowest-risk level = 0 points) ----
LABELS = {GENE:f"{GENE} log2(CPM+1)", "age":"Age (years)", "T34":"T stage"}
ranges = {v:(X[v].min(), X[v].max()) for v in PRED}
g = lambda v,x: b[v]*x                                  # log-odds contribution
ref = {v: min(g(v,ranges[v][0]), g(v,ranges[v][1])) for v in PRED}   # least-risk
span = {v: abs(b[v])*(ranges[v][1]-ranges[v][0]) for v in PRED}
maxc = max(span.values())
def pts(v, x): return 100*(g(v,x)-ref[v])/maxc          # always >= 0
order = sorted(PRED, key=lambda v:-span[v])
yrow = {v:i for i,v in enumerate(["points"]+order+["total","lp"])}
ytop = len(yrow)-1
def Y(name): return ytop - yrow[name]
LX = -6                                                 # left label x

# Points ruler 0-100
axN.plot([0,100],[Y("points")]*2, "k-", lw=1)
for t in range(0,101,10):
    axN.plot([t,t],[Y("points")-0.08,Y("points")+0.08],"k-",lw=0.8)
    axN.text(t, Y("points")+0.20, str(t), ha="center", fontsize=8)
axN.text(LX, Y("points"), "Points", ha="right", va="center", fontweight="bold")

# predictor axes
for v in order:
    lo,hi = ranges[v]
    if v=="T34":
        ticks=[0,1]; labs=["T1/T2","T3/T4"]
    elif v=="age":
        ticks=np.arange(20, hi+1, 15); labs=[f"{t:.0f}" for t in ticks]
    else:
        ticks=np.linspace(lo,hi,6); labs=[f"{t:.1f}" for t in ticks]
    p = [pts(v,t) for t in ticks]
    axN.plot([min(p),max(p)],[Y(v)]*2,"k-",lw=1)
    for t,pp,lb in zip(ticks,p,labs):
        axN.plot([pp,pp],[Y(v)-0.07,Y(v)+0.07],"k-",lw=0.8)
        axN.text(pp, Y(v)-0.24, lb, ha="center", fontsize=7.5)
    axN.text(LX, Y(v), LABELS[v], ha="right", va="center")

# total points
tot_hi = sum(max(pts(v,ranges[v][0]),pts(v,ranges[v][1])) for v in PRED)
axN.plot([0,tot_hi],[Y("total")]*2,"k-",lw=1)
for t in np.linspace(0,tot_hi,7):
    axN.plot([t,t],[Y("total")-0.08,Y("total")+0.08],"k-",lw=0.8)
    axN.text(t, Y("total")+0.20, f"{t:.0f}", ha="center", fontsize=8)
axN.text(LX, Y("total"), "Total Points", ha="right", va="center", fontweight="bold")

# probability axis: tot = (100/maxc)*(lp - const - sum ref)
base = b["const"] + sum(ref[v] for v in PRED)
def prob_to_tot(p):
    lp = np.log(p/(1-p)); return (100/maxc)*(lp - base)
probs=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8]
xs=[prob_to_tot(p) for p in probs]
vis=[(p,x) for p,x in zip(probs,xs) if -2<=x<=tot_hi+5]
axN.plot([vis[0][1],vis[-1][1]],[Y("lp")]*2,"k-",lw=1)
for p,xp in vis:
    axN.plot([xp,xp],[Y("lp")-0.08,Y("lp")+0.08],"k-",lw=0.8)
    axN.text(xp, Y("lp")-0.24, f"{p:.1f}", ha="center", fontsize=7.5)
axN.text(LX, Y("lp"), "Pr(N1 metastasis)", ha="right", va="center", fontweight="bold")
axN.set_xlim(-34,108); axN.set_ylim(-0.7, ytop+0.8); axN.axis("off")
axN.set_title("Nomogram for predicting lymph-node metastasis (N1)", fontweight="bold", fontsize=12)

# ---- ROC ----
fpr,tpr,_ = roc_curve(y, m.predict(Xc))
axR.plot(fpr,tpr,color=ns.C["tumor"],lw=2,label=f"Apparent AUC={auc_app:.3f}")
fprc,tprc,_ = roc_curve(y, cv_prob)
axR.plot(fprc,tprc,color=ns.C["normal"],lw=1.6,ls="-",label=f"10-fold CV AUC={auc_cv:.3f}")
axR.plot([0,1],[0,1],"--",color=ns.C["grey_mid"],lw=1)
axR.set_xlabel("1 - Specificity"); axR.set_ylabel("Sensitivity")
axR.set_title("Model ROC", fontweight="bold"); axR.legend(loc="lower right", frameon=False, fontsize=8)
axR.spines[["top","right"]].set_visible(False)

# ---- calibration ----
pt,pp = calibration_curve(y, cv_prob, n_bins=8, strategy="quantile")
axC.plot([0,1],[0,1],"--",color=ns.C["grey_mid"],lw=1,label="Ideal")
axC.plot(pp,pt,"o-",color=ns.C["tumor"],lw=1.8,ms=5,label="Observed (CV)")
axC.set_xlabel("Predicted probability"); axC.set_ylabel("Observed frequency")
axC.set_title("Calibration", fontweight="bold"); axC.legend(loc="upper left", frameon=False, fontsize=8)
axC.set_xlim(0,1); axC.set_ylim(0,1); axC.spines[["top","right"]].set_visible(False)

# ---- DCA ----
def net_benefit(yt, p, th):
    pred = p >= th; n=len(yt)
    tp=np.sum(pred & (yt==1)); fp=np.sum(pred & (yt==0))
    return tp/n - fp/n*(th/(1-th))
ths=np.linspace(0.01,0.8,80)
nb_model=[net_benefit(y.values,cv_prob,t) for t in ths]
prev=y.mean()
nb_all=[prev - (1-prev)*(t/(1-t)) for t in ths]
axD.plot(ths,nb_model,color=ns.C["tumor"],lw=2,label="Nomogram")
axD.plot(ths,nb_all,color=ns.C["grey_mid"],lw=1.2,label="Treat all")
axD.axhline(0,color="black",lw=1,label="Treat none")
axD.set_ylim(-0.05, max(nb_model)*1.15)
axD.set_xlabel("Threshold probability"); axD.set_ylabel("Net benefit")
axD.set_title("Decision curve (DCA)", fontweight="bold"); axD.legend(loc="upper right", frameon=False, fontsize=8)
axD.spines[["top","right"]].set_visible(False)

ns.panel_label(axN, "a")
ns.panel_label(axR, "b")
ns.panel_label(axC, "c")
ns.panel_label(axD, "d")
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_nomogram")

# save coefficient table
ct = m.summary2().tables[1].reset_index().rename(columns={"index":"term"})
ct.to_csv("data/PHLDA3_nomogram_model.csv", index=False)
print("saved table: PHLDA3_nomogram_model.csv")
