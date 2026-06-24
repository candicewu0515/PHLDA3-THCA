#!/usr/bin/env python3
"""PHLDA3 x THCA clinical: baseline table (high vs low) + clinicopathologic correlation + survival feasibility."""
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, kruskal, fisher_exact, chi2_contingency
import nature_style as ns
ns.set_style(font_size=8)

GENE, ENS = "PHLDA3", "ENSG00000174307"

# ---- expression (tumor only) ----
df = pd.read_csv("TCGA-THCA.htseq_counts.tsv", sep="\t", index_col=0)
mat = (2**df - 1).clip(lower=0)
libsize = mat.sum(axis=0)
row = mat.index[mat.index.str.startswith(ENS)][0]
logcpm = np.log2(mat.loc[row] / libsize * 1e6 + 1)
stype = pd.Series([b.split("-")[3][:2] for b in df.columns], index=df.columns)
tum = logcpm[stype == "01"]
expr = pd.DataFrame({"patient": [b[:12] for b in tum.index], GENE: tum.values})
expr = expr.groupby("patient")[GENE].mean()   # collapse dup samples per patient

# ---- clinical ----
cl = pd.read_csv("THCA_clinical.tsv", sep="\t", dtype=str)
def pick(r, base):  # prefer diagnoses.0, fallback .1
    for k in [f"diagnoses.0.{base}", f"diagnoses.1.{base}"]:
        if k in r and pd.notna(r[k]) and str(r[k]) not in ("", "nan"): return r[k]
    return np.nan
c = pd.DataFrame({
    "patient": cl["submitter_id"],
    "age": pd.to_numeric(cl["demographic.age_at_index"], errors="coerce"),
    "gender": cl["demographic.gender"],
    "vital": cl["demographic.vital_status"],
    "T": cl.apply(lambda r: pick(r, "ajcc_pathologic_t"), axis=1),
    "N": cl.apply(lambda r: pick(r, "ajcc_pathologic_n"), axis=1),
    "M": cl.apply(lambda r: pick(r, "ajcc_pathologic_m"), axis=1),
    "stage": cl.apply(lambda r: pick(r, "ajcc_pathologic_stage"), axis=1),
})
d = c.merge(expr.rename(GENE), on="patient", how="inner")
print(f"merged tumor patients with expression+clinical: {len(d)}")
print("vital_status:", d["vital"].value_counts().to_dict())

# ---- group definitions ----
d["T_grp"] = d["T"].str[:2].map(lambda t: "T1/T2" if t in ("T1","T2") else ("T3/T4" if t in ("T3","T4") else np.nan))
d["N_grp"] = d["N"].str[:2].map(lambda n: "N0" if n=="N0" else ("N1" if n=="N1" else np.nan))
d["stage_grp"] = d["stage"].map(lambda s: "I/II" if s in ("Stage I","Stage II") else ("III/IV" if s in ("Stage III","Stage IVA","Stage IVB","Stage IVC","Stage IV") else np.nan))
med = d[GENE].median()
d["PHLDA3_grp"] = np.where(d[GENE] > med, "High", "Low")
d["age_grp"] = np.where(d["age"] > d["age"].median(), ">median", "<=median")

# ---- baseline table: PHLDA3 High vs Low across features ----
rows = []
for feat in ["age_grp","gender","T_grp","N_grp","M","stage_grp"]:
    sub = d[[feat,"PHLDA3_grp"]].dropna()
    ct = pd.crosstab(sub[feat], sub["PHLDA3_grp"])
    if ct.shape == (2,2):
        _, p = fisher_exact(ct)
    else:
        _, p, _, _ = chi2_contingency(ct)
    rows.append({"feature": feat, "p_value": p, "levels": "/".join(ct.index.astype(str))})
bt = pd.DataFrame(rows)
bt.to_csv("PHLDA3_baseline_table.csv", index=False)
print("\n=== baseline table (PHLDA3 High vs Low) ===")
print(bt.to_string(index=False))

# ---- clinicopathologic correlation: expression vs groups ----
def mwu(a, b): return mannwhitneyu(a, b, alternative="two-sided")[1]
panels = [("T_grp", ["T1/T2","T3/T4"]), ("N_grp", ["N0","N1"]), ("stage_grp", ["I/II","III/IV"])]
fig, ax = plt.subplots(1, 3, figsize=(11, 4))
for i,(feat, order) in enumerate(panels):
    sub = d[[feat, GENE]].dropna()
    groups = [sub[sub[feat]==o][GENE].values for o in order]
    p = mwu(*groups)
    bp = ax[i].boxplot(groups, widths=0.55, patch_artist=True, showfliers=False,
                       medianprops=dict(color="black", lw=1.3))
    for patch,cc in zip(bp["boxes"], [ns.C["normal"], ns.C["tumor"]]):
        patch.set_facecolor(cc); patch.set_alpha(0.6)
    for j,(g,cc) in enumerate(zip(groups,[ns.C["normal"], ns.C["tumor"]]),1):
        ax[i].scatter(np.random.normal(j,0.06,len(g)), g, s=6, color=cc, alpha=0.4, edgecolors="none")
    ax[i].set_xticks([1,2]); ax[i].set_xticklabels([f"{o}\n(n={len(g)})" for o,g in zip(order,groups)])
    ax[i].set_ylabel(f"{GENE} log2(CPM+1)" if i==0 else "")
    pstr = f"p = {p:.1e}" if p>=1e-4 else "p < 1e-4"
    ax[i].set_title(f"{feat.replace('_grp','')}   {pstr}", fontweight="bold", fontsize=11)
    ax[i].spines[["top","right"]].set_visible(False)
    ns.panel_label(ax[i], "abc"[i])
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_clinical")
print("\nsaved: PHLDA3_baseline_table.csv")
