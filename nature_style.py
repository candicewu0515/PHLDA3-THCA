#!/usr/bin/env python3
"""Shared Nature-style theme for ALL PHLDA3-THCA figures — one palette, one font,
one export path. Import this in every figure script:

    import nature_style as ns
    ns.set_style()                 # call once, before plotting
    ... use ns.C["tumor"] etc ...
    ns.save_fig(fig, "PHLDA3_THCA_xxx")   # writes .svg/.pdf/.tiff/.png
"""
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# ── unified SEMANTIC palette (use these names everywhere) ───────────────────
C = {
    # signal axis: tumour / high-expression / up-regulated / hero  = red
    "tumor":  "#B64342", "high": "#B64342", "up": "#B64342", "hero": "#B64342",
    # reference axis: normal / low-expression / down-regulated     = blue
    "normal": "#3775BA", "low": "#3775BA", "down": "#3775BA", "ref": "#3775BA",
    # third / intermediate category (e.g. RAS, DC)                 = teal
    "accent": "#42949E",
    # extreme category (e.g. ATC, BRAF-V600E emphasis)             = dark red
    "accent2": "#7E2A2A",
    # neutral / not-significant
    "ns": "#CFCECE", "neutral": "#9A9A9A", "grey_light": "#CFCECE",
    "grey_mid": "#9A9A9A", "grey_dark": "#4D4D4D", "black": "#272727",
    # special highlight (single called-out item, e.g. THCA in pan-cancer)
    "highlight": "#E08214",
}
# legacy ad-hoc hex -> semantic (for mechanical replacement in old scripts)
LEGACY = {"#C44E52": C["tumor"], "#4C72B0": C["normal"], "#55A868": C["accent"],
          "#8C2D2D": C["accent2"], "#BFBFBF": C["ns"], "#E41A1C": C["highlight"],
          "#888": C["grey_mid"], "#8C2D2D2": C["accent2"]}

# sequential colormap for continuous fills (white -> signal red)
SEQ_RED = LinearSegmentedColormap.from_list("seq_red", ["#FBF0EF", C["tumor"], "#5E1F1F"])
# diverging (down blue - white - up red)
DIV = LinearSegmentedColormap.from_list("div_br", [C["normal"], "#F2F2F2", C["tumor"]])

def set_style(font_size=8):
    """Nature rcParams: Arial, editable vector text, thin despined axes."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "svg.fonttype": "none",     # editable text in SVG
        "pdf.fonttype": 42,         # editable TrueType in PDF
        "font.size": font_size,
        "axes.titlesize": font_size + 1,
        "axes.labelsize": font_size,
        "axes.linewidth": 0.7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.major.width": 0.7, "ytick.major.width": 0.7,
        "xtick.major.size": 2.5,  "ytick.major.size": 2.5,
        "xtick.labelsize": font_size - 1, "ytick.labelsize": font_size - 1,
        "legend.frameon": False, "legend.fontsize": font_size - 1,
        "figure.dpi": 120,
    })

def panel_label(ax, lab, x=-0.12, y=1.04, fs=11):
    ax.text(x, y, lab, transform=ax.transAxes, fontsize=fs, fontweight="bold",
            va="bottom", ha="left")

def save_fig(fig, base, dpi=600, formats=("pdf",)):
    """Submission figure export. Default: PDF only (editable vector, submission-
    ready). Pass formats=("pdf","png") for a preview raster, or add "tiff"/"svg"
    only if a journal specifically requires them at final submission."""
    for fmt in formats:
        kw = {"bbox_inches": "tight"}
        if fmt in ("png", "tiff"): kw["dpi"] = dpi
        if fmt == "tiff": kw["pil_kwargs"] = {"compression": "tiff_lzw"}
        fig.savefig(f"{base}.{fmt}", **kw)
    print(f"saved: {base}." + "/.".join(formats))
