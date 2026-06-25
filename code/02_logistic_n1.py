#!/usr/bin/env python3
"""PHLDA3 x THCA: logistic regression — does PHLDA3 INDEPENDENTLY predict
lymph-node metastasis (N1) and late stage (III/IV)? Univariate + multivariable,
with an OR forest plot for the manuscript."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import statsmodels.api as sm
import nature_style as ns

ns.set_style(font_size=8)

GENE, ENS = "PHLDA3", "ENSG00000174307"

# ---------- expression (tumor only, per-patient) ----------
df = pd.read_csv("data/TCGA-THCA.htseq_counts.tsv", sep="\t", index_col=0)
mat = (2**df - 1).clip(lower=0)
libsize = mat.sum(axis=0)
row = mat.index[mat.index.str.startswith(ENS)][0]
logcpm = np.log2(mat.loc[row] / libsize * 1e6 + 1)
stype = pd.Series([b.split("-")[3][:2] for b in df.columns], index=df.columns)
tum = logcpm[stype == "01"]
expr = pd.DataFrame({"patient": [b[:12] for b in tum.index], GENE: tum.values}) \
        .groupby("patient")[GENE].mean()

# ---------- clinical ----------
cl = pd.read_csv("data/THCA_clinical.tsv", sep="\t", dtype=str)
def pick(r, base):
    for k in [f"diagnoses.0.{base}", f"diagnoses.1.{base}"]:
        if k in r and pd.notna(r[k]) and str(r[k]) not in ("", "nan"): return r[k]
    return np.nan
c = pd.DataFrame({
    "patient": cl["submitter_id"],
    "age": pd.to_numeric(cl["demographic.age_at_index"], errors="coerce"),
    "gender": cl["demographic.gender"],
    "T": cl.apply(lambda r: pick(r, "ajcc_pathologic_t"), axis=1),
    "N": cl.apply(lambda r: pick(r, "ajcc_pathologic_n"), axis=1),
    "stage": cl.apply(lambda r: pick(r, "ajcc_pathologic_stage"), axis=1),
})
d = c.merge(expr.rename(GENE), on="patient", how="inner")

# ---------- variable engineering ----------
# Predictor of interest: PHLDA3 standardized -> OR is "per 1 SD increase"
d["PHLDA3_z"] = (d[GENE] - d[GENE].mean()) / d[GENE].std()
med = d[GENE].median()
d["PHLDA3_high"] = (d[GENE] > med).astype(int)        # High(1) vs Low(0)
# Covariates
d["age_c"] = d["age"]                                  # continuous (years)
d["male"] = (d["gender"] == "male").astype(int)
d["T34"]  = d["T"].str[:2].map({"T1":0,"T2":0,"T3":1,"T4":1})
# Outcomes
d["N1"]   = d["N"].str[:2].map({"N0":0,"N1":1})
d["late"] = d["stage"].map(lambda s: 0 if s in ("Stage I","Stage II")
            else (1 if s in ("Stage III","Stage IV","Stage IVA","Stage IVB","Stage IVC") else np.nan))

print(f"patients merged: {len(d)}")
print(f"N1 outcome: {d['N1'].value_counts().to_dict()}  | late stage: {d['late'].value_counts().to_dict()}")

# ---------- logistic regression helper ----------
def fit(data, outcome, predictors):
    sub = data[[outcome] + predictors].dropna()
    X = sm.add_constant(sub[predictors].astype(float))
    y = sub[outcome].astype(int)
    m = sm.Logit(y, X).fit(disp=0)
    out = []
    for v in predictors:
        OR = np.exp(m.params[v]); lo, hi = np.exp(m.conf_int().loc[v])
        out.append({"outcome": outcome, "var": v, "n": len(sub),
                    "OR": OR, "lo": lo, "hi": hi, "p": m.pvalues[v]})
    return pd.DataFrame(out)

PRETTY = {"PHLDA3_z": f"{GENE} (per 1 SD)", "PHLDA3_high": f"{GENE} High vs Low",
          "age_c": "Age (per year)", "male": "Male vs Female", "T34": "T3/T4 vs T1/T2"}

results = []
# ---- Outcome 1: lymph-node metastasis (N1) ----
results.append(("Univariate",   fit(d, "N1", ["PHLDA3_z"])))
results.append(("Multivariable",fit(d, "N1", ["PHLDA3_z", "age_c", "male", "T34"])))
# ---- Outcome 2: late stage (III/IV) ----
results.append(("Univariate",   fit(d, "late", ["PHLDA3_z"])))
results.append(("Multivariable",fit(d, "late", ["PHLDA3_z", "age_c", "male", "T34"])))

tbl = pd.concat([r.assign(model=lbl) for lbl, r in results], ignore_index=True)
tbl["sig"] = tbl["p"].map(lambda p: "***" if p<1e-3 else "**" if p<1e-2 else "*" if p<0.05 else "ns")
tbl_out = tbl[["outcome","model","var","n","OR","lo","hi","p","sig"]].copy()
tbl_out.to_csv("data/PHLDA3_logistic_table.csv", index=False)
print("\n=== logistic regression (OR [95% CI], p) ===")
for _, r in tbl_out.iterrows():
    print(f"{r['outcome']:5s} {r['model']:13s} {PRETTY.get(r['var'],r['var']):22s} "
          f"OR={r['OR']:.2f} [{r['lo']:.2f}, {r['hi']:.2f}]  p={r['p']:.2e} {r['sig']}")

# ---------- forest plot (multivariable models) ----------
OUT = {"N1": "Lymph-node metastasis (N1 vs N0)", "late": "Late stage (III/IV vs I/II)"}
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
for pidx, (ax, oc) in enumerate(zip(axes, ["N1", "late"])):
    mv = tbl[(tbl.outcome==oc) & (tbl.model=="Multivariable")].iloc[::-1].reset_index(drop=True)
    ys = np.arange(len(mv))
    for i, r in mv.iterrows():
        is_gene = r["var"] == "PHLDA3_z"
        col = ns.C["tumor"] if is_gene else ns.C["normal"]
        ax.plot([r.lo, r.hi], [i, i], color=col, lw=1.8, zorder=2)
        ax.scatter(r.OR, i, s=70 if is_gene else 45, color=col, zorder=3,
                   edgecolors="black", linewidths=0.6)
        ax.text(ax.get_xlim()[1], i, f"  {r.OR:.2f} ({r.lo:.2f}–{r.hi:.2f}) {r.sig}",
                va="center", fontsize=8.5)
    ax.axvline(1, ls="--", color="grey", lw=1, zorder=1)
    ax.set_yticks(ys); ax.set_yticklabels([PRETTY.get(v, v) for v in mv["var"]])
    ax.set_xscale("log")
    ax.set_xlabel("Odds ratio (95% CI), log scale")
    ax.set_title(f"{OUT[oc]}\n(multivariable, n={int(mv['n'].iloc[0])})",
                 fontweight="bold", fontsize=10)
    ax.spines[["top","right"]].set_visible(False)
    ax.margins(x=0.45)
    ns.panel_label(ax, "a" if pidx == 0 else "b")
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_logistic")
print("saved table: PHLDA3_logistic_table.csv")
