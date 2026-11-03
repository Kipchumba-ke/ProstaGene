#!/usr/bin/env Rscript
# Checks that the pipeline recovers the effects built into the synthetic dataset.
suppressPackageStartupMessages(library(data.table))
fail <- character()
check <- function(ok, msg) { cat(if (ok) "PASS " else "FAIL ", msg, "\n", sep = ""); if (!ok) fail <<- c(fail, msg) }

de <- fread("results/de/deseq2_all_genes.tsv")
h  <- fread("results/gsea/fgsea_hallmark.tsv")
nes <- function(p) h[pathway == p, NES]; padj <- function(p) h[pathway == p, padj]
dir <- function(g) de[gene_name == g, direction][1]
called <- function(g, want) {
  d <- dir(g)
  if (is.na(d)) check(FALSE, paste(g, "missing from DE results (filtered out at QC?)"))
  else check(d == want, sprintf("%s called %s in tumor (got: %s)", g, want, d))
}

for (g in c("AMACR", "PCA3", "HPN", "DLX1", "TDRD1")) called(g, "up")
for (g in c("KRT5", "KRT15", "TP63")) called(g, "down")
check(nes("HALLMARK_MYC_TARGETS_V1") > 0 && padj("HALLMARK_MYC_TARGETS_V1") < 0.05, "MYC targets V1 enriched up")
check(nes("HALLMARK_E2F_TARGETS") > 0 && padj("HALLMARK_E2F_TARGETS") < 0.05, "E2F targets enriched up")
check(nes("HALLMARK_MYOGENESIS") < 0 && padj("HALLMARK_MYOGENESIS") < 0.05, "Myogenesis enriched down")
check(file.exists("results/qc/qc_sample_flags.tsv"), "QC flag table written")
check(file.exists("results/de/volcano.png"), "volcano plot written")

if (length(fail)) stop(length(fail), " check(s) failed: ", paste(fail, collapse = "; "))
cat("All checks passed.\n")
