#!/usr/bin/env python3
"""Nature-style combined single-cell figure.
Core conclusion: PHLDA3 is near-absent in normal thyroid follicular cells but is the
dominant PHLDA3+ population in malignant tumor cells -> induced upon transformation.
Backend: Python/matplotlib (exclusive). Source data: thyroid_sc_phlda3.csv (normal,
CELLxGENE Census) + PHLDA3_tumor_sc.csv (tumor, GSE193581)."""
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, numpy as np

# ── MANDATORY editable-text + Nature rcParams ──────────────────────────────
mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "svg.fonttype": "none",     # editable SVG text
    "pdf.fonttype": 42,         # editable PDF text
    "font.size": 7,
    "axes.linewidth": 0.6,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,  "ytick.major.size": 2.5,
    "legend.frameon": False,
})
NEUTRAL, STEM, HERO = "#9A9A9A", "#CFCECE", "#B64342"   # neutral dot, stem, signal red

# ── data ───────────────────────────────────────────────────────────────────
def summarize(df, ctcol, valcol, nmin):
    g = (df.groupby(ctcol)
           .agg(n=(valcol, "size"), pct=(valcol, lambda x: (x > 0).mean()*100))
           .reset_index().rename(columns={ctcol: "celltype"}))
    return g[g.n >= nmin].sort_values("pct").reset_index(drop=True)

gn = summarize(pd.read_csv("thyroid_sc_phlda3.csv"), "cell_type", "PHLDA3", 50)
gt = summarize(pd.read_csv("PHLDA3_tumor_sc.csv"),  "celltype",  "raw",    30)
xmax = max(gn.pct.max(), gt.pct.max()) * 1.18

# ── figure: Nature double-column 183 mm ─────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3))
plt.subplots_adjust(left=0.20, right=0.985, top=0.84, bottom=0.16, wspace=1.05)

def lollipop(ax, g, hero, title, subtitle):
    for i, r in g.iterrows():
        is_h = r.celltype == hero
        ax.hlines(i, 0, r.pct, color=STEM, lw=1.1, zorder=1)
        ax.plot(r.pct, i, "o", ms=5.5 if is_h else 3.6,
                color=HERO if is_h else NEUTRAL,
                mec="black", mew=0.5, zorder=3)
        ax.text(r.pct + xmax*0.02, i, f"{r.pct:.1f}%", va="center", ha="left",
                fontsize=6, color=HERO if is_h else "0.45",
                fontweight="bold" if is_h else "normal")
    ax.set_yticks(range(len(g)))
    ax.set_yticklabels([f"{c}  ({n:,})" for c, n in zip(g.celltype, g.n)], fontsize=6.3)
    for lbl in ax.get_yticklabels():
        if lbl.get_text().startswith(hero): lbl.set_color(HERO); lbl.set_fontweight("bold")
    ax.set_xlim(0, xmax); ax.set_ylim(-0.6, len(g)-0.4)
    ax.set_xlabel("Cells expressing PHLDA3 (%)", fontsize=7)
    ax.set_title(title, fontsize=7.5, fontweight="bold", pad=12)
    ax.annotate(subtitle, xy=(0.5, 1.0), xycoords="axes fraction",
                xytext=(0, 3), textcoords="offset points",
                ha="center", va="bottom", fontsize=5.8, color="0.4")

# hero row faint highlight
def hero_band(ax, g, hero):
    i = g.index[g.celltype == hero]
    if len(i): ax.axhspan(i[0]-0.45, i[0]+0.45, color=HERO, alpha=0.07, zorder=0)

lollipop(axes[0], gn, "thyroid follicular cell", "Normal thyroid",
         "CELLxGENE Census · 15,646 cells")
lollipop(axes[1], gt, "Malignant cell", "Thyroid tumour",
         "GSE193581 · 67,414 cells · PTC/ATC")
hero_band(axes[0], gn, "thyroid follicular cell")
hero_band(axes[1], gt, "Malignant cell")

# panel labels a / b
for ax, lab in zip(axes, ["a", "b"]):
    ax.text(-0.55, 1.06, lab, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top", ha="left")

# ── export: editable SVG/PDF + 600-dpi TIFF ────────────────────────────────
base = "PHLDA3_THCA_singlecell_combined_nature"
fig.savefig(f"{base}.svg", bbox_inches="tight")
fig.savefig(f"{base}.pdf", bbox_inches="tight")
fig.savefig(f"{base}.tiff", dpi=600, bbox_inches="tight", pil_kwargs={"compression":"tiff_lzw"})
fig.savefig(f"{base}.png", dpi=600, bbox_inches="tight")
plt.close(fig)
print("follicular (normal) %:", round(gn.loc[gn.celltype=='thyroid follicular cell','pct'].iloc[0],2))
print("malignant (tumour)  %:", round(gt.loc[gt.celltype=='Malignant cell','pct'].iloc[0],2))
print(f"saved: {base}.svg/.pdf/.tiff/.png")
