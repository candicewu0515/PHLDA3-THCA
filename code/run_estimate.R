#!/usr/bin/env Rscript
# ESTIMATE stromal/immune/ESTIMATE scores from a gene-symbol x sample expression matrix.
suppressMessages(library(estimate))
args <- commandArgs(trailingOnly = TRUE)             # input.tsv  output_scores.gct
filterCommonGenes(input.f = args[1], output.f = "data/estimate_filtered.gct", id = "GeneSymbol")
estimateScore("data/estimate_filtered.gct", args[2], platform = "illumina")
cat("ESTIMATE done\n")
