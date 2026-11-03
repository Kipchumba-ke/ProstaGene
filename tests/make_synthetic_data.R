#!/usr/bin/env Rscript
# ---------------------------------------------------------------
# make_synthetic_data.R — small synthetic paired dataset for CI
#
# 8 patients x (normal, tumor) over the MSigDB Hallmark genes plus the
# QC marker genes, negative-binomial counts with KNOWN effects:
#   MYC targets V1 and E2F targets  up   in tumor (log2FC ~ +1)
#   Myogenesis                      down in tumor (log2FC ~ -1.2)
#   AMACR, PCA3, HPN, DLX1, TDRD1   up   (log2FC +3)
#   KRT5, KRT15, TP63, MYLK, CNN1   down (log2FC -2)
# Written where the pipeline expects the paired SummarizedExperiment,
# so CI can run qc -> deseq2 -> gsea without downloading TCGA.
# tests/check_results.R then checks those effects are recovered.
# ---------------------------------------------------------------
suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(msigdbr)
  library(yaml)
})
set.seed(1)
cfg <- yaml::read_yaml("config/config.yaml")

new_api <- "collection" %in% names(formals(msigdbr::msigdbr))
h <- if (new_api) msigdbr(species = "Homo sapiens", collection = "H") else msigdbr(species = "Homo sapiens", category = "H")
up_tumor   <- c("AMACR", "PCA3", "HPN", "DLX1", "TDRD1")
down_tumor <- c("KRT5", "KRT15", "TP63", "MYLK", "CNN1")
genes <- sort(unique(c(h$gene_symbol, up_tumor, down_tumor)))
G <- length(genes)

lfc <- setNames(rep(0, G), genes)
set_lfc <- function(set, value, sd) { g <- intersect(unique(h$gene_symbol[h$gs_name == set]), genes)
  lfc[g] <<- value + rnorm(length(g), 0, sd) }
set_lfc("HALLMARK_MYC_TARGETS_V1", 1.0, 0.3)
set_lfc("HALLMARK_E2F_TARGETS", 1.0, 0.3)
set_lfc("HALLMARK_MYOGENESIS", -1.2, 0.3)
lfc[up_tumor] <- 3; lfc[down_tumor] <- -2

n_pairs <- 8
patients <- sprintf("TCGA-SY-%04d", seq_len(n_pairs))
cd <- data.frame(
  patient   = rep(patients, each = 2),
  condition = factor(rep(c("normal", "tumor"), n_pairs), levels = c("normal", "tumor")),
  stringsAsFactors = FALSE)
cd$tss   <- ifelse(seq_len(nrow(cd)) %% 4 < 2, "SY", "SZ")
cd$plate <- ifelse(rep(seq_len(n_pairs), each = 2) %% 2 == 0, "A000", "A001")
rownames(cd) <- paste0(cd$patient, ifelse(cd$condition == "tumor", "-01A", "-11A"), "-11R-", cd$plate, "-07")

base_mu <- exp(rnorm(G, 5, 1.5))
size_factor <- exp(rnorm(nrow(cd), 0, 0.15))
patient_effect <- matrix(rnorm(G * n_pairs, 0, 0.25), G, n_pairs)[, match(cd$patient, patients)]
mu <- base_mu * 2^(outer(lfc, as.numeric(cd$condition == "tumor"))) * exp(patient_effect) *
      matrix(size_factor, G, nrow(cd), byrow = TRUE)
counts <- matrix(rnbinom(length(mu), size = 20, mu = mu), G, nrow(cd), dimnames = list(sprintf("SYN%05d", seq_len(G)), rownames(cd)))
storage.mode(counts) <- "integer"

se <- SummarizedExperiment(assays = setNames(list(counts), cfg$count_assay),
                           rowData = data.frame(gene_id = rownames(counts), gene_name = genes, gene_type = "protein_coding"),
                           colData = cd)
dir.create(dirname(cfg$outputs$se_paired), recursive = TRUE, showWarnings = FALSE)
saveRDS(se, cfg$outputs$se_paired)
write.table(data.frame(barcode = rownames(cd), cd), cfg$outputs$sample_sheet, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Synthetic dataset: %d genes x %d samples (%d pairs) -> %s\n", G, ncol(se), n_pairs, cfg$outputs$se_paired))
