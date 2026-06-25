#!/usr/bin/env python3
"""Fetch normal-thyroid single-cell PHLDA3 expression from the CZI CELLxGENE Census
-> thyroid_sc_phlda3.csv  (consumed by combined_singlecell_nature.py).
Requires network access; census_version pinned for reproducibility."""
import cellxgene_census, numpy as np, pandas as pd

CENSUS_VERSION = "2025-11-08"   # pinned stable release used in the manuscript
with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
    ad = cellxgene_census.get_anndata(
        census, organism="Homo sapiens",
        var_value_filter="feature_name == 'PHLDA3'",
        obs_value_filter="tissue == 'thyroid gland' and is_primary_data == True",
        obs_column_names=["cell_type", "disease", "assay"])
x = np.asarray(ad[:, ad.var.feature_name == "PHLDA3"].X.todense()).ravel()
df = ad.obs.copy(); df["PHLDA3"] = x
df.to_csv("data/thyroid_sc_phlda3.csv", index=False)
print(f"saved thyroid_sc_phlda3.csv: {len(df)} cells, {df['cell_type'].nunique()} cell types")
