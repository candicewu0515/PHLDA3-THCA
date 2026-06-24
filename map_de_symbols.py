#!/usr/bin/env python3
"""Map Ensembl IDs in the differential-expression table to gene symbols
(thca_DE_full.csv -> thca_DE_full_symbol.csv) via mygene.
Consumed by volcano_thca.py, coexpr_enrich_thca.py, immune_thca.py.
Requires network access (mygene.info)."""
import pandas as pd, mygene

res = pd.read_csv("thca_DE_full.csv")
mg = mygene.MyGeneInfo()
q = mg.querymany(res["ensembl"].tolist(), scopes="ensembl.gene",
                 fields="symbol,type_of_gene", species="human", verbose=False)
m = {d["query"]: (d.get("symbol"), d.get("type_of_gene")) for d in q if not d.get("notfound")}
res["symbol"]    = res["ensembl"].map(lambda e: m.get(e, (None, None))[0])
res["gene_type"] = res["ensembl"].map(lambda e: m.get(e, (None, None))[1])
res.to_csv("thca_DE_full_symbol.csv", index=False)
print(f"mapped {res['symbol'].notna().sum()}/{len(res)} genes -> thca_DE_full_symbol.csv")
