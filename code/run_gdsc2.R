#!/usr/bin/env Rscript
# GDSC2 drug-sensitivity: Spearman of PHLDA3 expression vs log(IC50) across cell
# lines -> data/gdsc2_phlda3_drugcorr.csv (consumed by gdsc2_phlda3.py).
# Run from repo root: Rscript code/run_gdsc2.R
expr <- readRDS("data/gdsc2/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
res  <- readRDS("data/gdsc2/GDSC2_Res.rds")
ph <- as.numeric(expr["PHLDA3", ]); names(ph) <- colnames(expr)
cells <- intersect(names(ph), rownames(res)); ph <- ph[cells]; R <- res[cells, ]
out <- data.frame(drug = colnames(R), rho = NA_real_, p = NA_real_, n = NA_integer_)
for (i in seq_len(ncol(R))) {
  x <- R[, i]; ok <- !is.na(x)
  ct <- suppressWarnings(cor.test(ph[ok], x[ok], method = "spearman"))
  out$rho[i] <- ct$estimate; out$p[i] <- ct$p.value; out$n[i] <- sum(ok)
}
out$fdr <- p.adjust(out$p, "BH")
write.csv(out, "data/gdsc2_phlda3_drugcorr.csv", row.names = FALSE)
cat("saved data/gdsc2_phlda3_drugcorr.csv:", nrow(out), "drugs\n")
