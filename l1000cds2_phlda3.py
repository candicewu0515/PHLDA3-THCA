#!/usr/bin/env python3
"""Task 7 (drug repurposing) — L1000CDS2 query with the PHLDA3-high signature.
Reverse mode: find L1000 perturbagens whose signature OPPOSES the PHLDA3-high
program (candidate compounds to revert the aggressive/de-differentiated state).
Independent, API-based complement to the clue.io CMap query."""
import requests, json, collections
import pandas as pd, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nature_style as ns
ns.set_style(font_size=8)

up = [g.strip() for g in open("PHLDA3_cmap_up.txt") if g.strip()]
dn = [g.strip() for g in open("PHLDA3_cmap_down.txt") if g.strip()]
print(f"signature: {len(up)} up / {len(dn)} down")

URL = "https://maayanlab.cloud/L1000CDS2/query"
payload = {"data": {"upGenes": up, "dnGenes": dn},
           "config": {"aggravate": False, "searchMethod": "geneSet",
                      "share": False, "combination": False, "db-version": "latest"},
           "metadata": []}
res = requests.post(URL, json=payload, timeout=120).json()
top = res.get("topMeta", [])
print(f"top reversing signatures returned: {len(top)}")

rows = [{"drug": t.get("pert_desc") or t.get("pert_id"),
         "pert_id": t.get("pert_id"), "score": t.get("score"),
         "cell": t.get("cell_id"), "dose": t.get("pert_dose"), "time": t.get("pert_time")}
        for t in top]
df = pd.DataFrame(rows)
df = df[df["drug"].notna() & (df["drug"] != "-666")]
df.to_csv("PHLDA3_L1000CDS2_hits.csv", index=False)

# aggregate by drug: best (lowest) score + frequency across cell lines/doses
agg = (df.groupby("drug")
         .agg(best_score=("score","min"), n_sig=("score","size"))
         .reset_index().sort_values(["best_score"]).head(20))
agg.to_csv("PHLDA3_L1000CDS2_drugs.csv", index=False)
print("\n=== top reversing compounds (lower score = stronger reversal) ===")
print(agg.head(15).to_string(index=False))

# ---------- figure: recurrence of reversing compounds (frequency = robustness) ----------
# scores are near-tied at the top, so recurrence across cell lines/doses is the signal
HDAC = {"trichostatin A","vorinostat","panobinostat","entinostat","mocetinostat","romidepsin","scriptaid"}
HSP90 = {"geldanamycin","tanespimycin","alvespimycin"}
freq = agg[agg.n_sig>=1].sort_values(["n_sig","best_score"], ascending=[True,False]).tail(15)
def cls_color(d):
    dl=d.lower()
    if dl in {x.lower() for x in HDAC}: return ns.C["tumor"]
    if dl in {x.lower() for x in HSP90}: return ns.C["accent"]
    return ns.C["ns"]
cols=[cls_color(d) for d in freq["drug"]]
fig, ax = plt.subplots(figsize=(6.8, 5))
ax.barh(range(len(freq)), freq["n_sig"], color=cols, edgecolor="black", lw=0.4)
ax.set_yticks(range(len(freq))); ax.set_yticklabels(freq["drug"], fontsize=7.5)
ax.set_xlabel("Recurrence among top-50 reversing L1000 signatures")
ax.set_title("Compounds reversing the PHLDA3-high signature (L1000CDS2)\nred = HDAC inhibitor · teal = HSP90 inhibitor",
             fontweight="bold", fontsize=8.5)
ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
ns.save_fig(fig, "PHLDA3_THCA_L1000CDS2")
print("done L1000CDS2")
