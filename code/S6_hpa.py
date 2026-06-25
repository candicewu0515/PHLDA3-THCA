#!/usr/bin/env python3
"""HPA protein-level: representative PHLDA3 IHC in normal thyroid vs thyroid cancer.
Qualitative confirmation of protein expression / subcellular localization (NOT a
quantitative tumor-vs-normal claim: HPA cancer n=4, all 'Medium')."""
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

panels = [
    ("data/hpa_normal_thyroid.jpg", "Normal thyroid gland",
     "Glandular cells: Medium\nIntensity: moderate"),
    ("data/hpa_thyroid_cancer.jpg", "Thyroid cancer (papillary)",
     "Tumor cells: Medium (4/4 patients)\nLocation: cytoplasmic + nuclear"),
]
plt.rcParams.update({"font.size": 11})
fig, ax = plt.subplots(1, 2, figsize=(10, 5.4))
for a, (f, title, note) in zip(ax, panels):
    im = Image.open(f)
    s = min(im.size); im = im.crop(((im.width-s)//2,(im.height-s)//2,(im.width+s)//2,(im.height+s)//2))
    a.imshow(im); a.set_xticks([]); a.set_yticks([])
    a.set_title(title, fontweight="bold", fontsize=12)
    a.text(0.02, 0.02, note, transform=a.transAxes, fontsize=9, va="bottom", ha="left",
           bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="grey", alpha=0.85))
fig.suptitle("PHLDA3 protein (HPA IHC, antibody HPA017385)", fontweight="bold", fontsize=13, y=0.99)
plt.tight_layout()
plt.savefig("PHLDA3_THCA_hpa_ihc.png", dpi=200, bbox_inches="tight")
plt.savefig("PHLDA3_THCA_hpa_ihc.pdf", bbox_inches="tight")
print("saved: PHLDA3_THCA_hpa_ihc.png/.pdf")
