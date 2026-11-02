#!/usr/bin/env Rscript
# ---------------------------------------------------------------
# 06_geo_validation.R — replicate TCGA-PRAD findings in GSE104131
#
# GSE104131: high-grade (Gleason >= 4+3) prostate tumors + matched
# adjacent normal, 8 African-American and 8 European-American men.
# Counts: NCBI-generated raw counts (GRCh38.p13, Entrez Gene IDs).
#
#  1. Build sample sheet from the series matrix; keep complete pairs
#  2. DESeq2 ~ pair + condition, apeglm shrinkage
#  3. Gene-level concordance with TCGA (sign agreement, correlation)
#  4. TCGA signature score: does it separate tumor from normal here?
#  5. Hallmark GSEA in GEO; NES concordance with TCGA
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(ggrepel)
})

if (exists("snakemake")) {
  inp <- snakemake@input; out <- snakemake@output; p <- snakemake@params
  threads <- snakemake@threads
  con <- file(snakemake@log[[1]], open = "wt"); sink(con); sink(con, type = "message")
} else {
  inp <- list(counts = "data/raw/geo/GSE104131_raw_counts_GRCh38.p13_NCBI.tsv.gz",
              matrix = "data/raw/geo/GSE104131_series_matrix.txt.gz",
              tcga_de = "results/de/deseq2_all_genes.tsv",
              tcga_hallmark = "results/gsea/fgsea_hallmark.tsv")
  v <- "results/validation"
  out <- list(sample_sheet = "data/processed/geo_sample_sheet.tsv",
              dds_fit = "data/processed/geo_dds_fitted.rds",
              de = file.path(v, "geo_deseq2_all_genes.tsv"),
              hallmark = file.path(v, "geo_fgsea_hallmark.tsv"),
              summary = file.path(v, "validation_summary.tsv"),
              lfc_plot = file.path(v, "lfc_concordance.png"),
              score_plot = file.path(v, "tcga_signature_score.png"),
              nes_plot = file.path(v, "hallmark_nes_concordance.png"))
  p <- list(fdr = 0.05, lfc = 1, min_count = 10, sig_size = 100, seed = 42)
  threads <- 4
}
for (f in unlist(out)) dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 11))
fdr <- p[["fdr"]]; lfc_cut <- p[["lfc"]]
set.seed(p[["seed"]])

# ---- 1. sample sheet ---------------------------------------------
pd <- Biobase::pData(getGEO(filename = inp[["matrix"]], getGPL = FALSE))
meta <- data.table(
  gsm       = pd$geo_accession,
  title     = pd$title,
  pair      = paste0("P", sub("^Sample_([0-9]+)_.*$", "\\1", pd$title)),
  condition = ifelse(grepl("tumor", pd$characteristics_ch1, ignore.case = TRUE), "tumor", "normal"),
  ethnicity = trimws(sub("^ethnicity:\\s*", "", pd$characteristics_ch1.1))
)
meta[, ancestry := fifelse(grepl("African", ethnicity), "AA",
                   fifelse(grepl("European", ethnicity), "EA", NA_character_))]

cts <- fread(inp[["counts"]])
gene_ids <- as.character(cts$GeneID)
cts <- as.matrix(cts[, -"GeneID"]); rownames(cts) <- gene_ids

meta[, in_counts := gsm %in% colnames(cts)]
missing <- meta[in_counts == FALSE]
message("Samples in series matrix: ", nrow(meta), " | in NCBI counts: ", sum(meta$in_counts))
if (nrow(missing)) { message("Missing from counts:"); print(missing[, .(gsm, pair, condition, ancestry)]) }

complete <- meta[in_counts == TRUE, .N, by = pair][N == 2, pair]
meta[, complete_pair := pair %in% complete]
fwrite(meta, out[["sample_sheet"]], sep = "\t")

m <- meta[complete_pair == TRUE][order(pair, condition)]
message("Complete pairs: ", length(complete), "  [AA ",
        uniqueN(m[ancestry == "AA", pair]), ", EA ", uniqueN(m[ancestry == "EA", pair]), "]")
stopifnot(length(complete) >= 5)

cd <- data.frame(pair = factor(m$pair),
                 condition = relevel(factor(m$condition), ref = "normal"),
                 ancestry = factor(m$ancestry),
                 row.names = m$gsm)
cts <- cts[, m$gsm]

# ---- 2. DESeq2 ---------------------------------------------------
sym <- mapIds(org.Hs.eg.db, keys = rownames(cts), keytype = "ENTREZID",
              column = "SYMBOL", multiVals = "first")
dds <- DESeqDataSetFromMatrix(cts, cd, design = ~ pair + condition)
mcols(dds)$symbol <- unname(sym[rownames(dds)])

min_n <- length(complete)                       # = samples per condition
keep <- rowSums(counts(dds) >= p[["min_count"]]) >= min_n
message("Genes before filtering: ", nrow(dds), " | after: ", sum(keep))
dds <- dds[keep, ]

bp <- if (threads > 1) BiocParallel::MulticoreParam(threads) else BiocParallel::SerialParam()
dds <- DESeq(dds, parallel = threads > 1, BPPARAM = bp, quiet = TRUE)
coef_name <- "condition_tumor_vs_normal"
r  <- results(dds, name = coef_name, alpha = fdr)
s  <- lfcShrink(dds, coef = coef_name, type = "apeglm", quiet = TRUE)
geo <- data.table(entrez = rownames(r), symbol = mcols(dds)$symbol,
                  baseMean = r$baseMean, log2FC = r$log2FoldChange,
                  log2FC_shrunk = s$log2FoldChange, stat = r$stat,
                  pvalue = r$pvalue, padj = r$padj)[order(padj)]
geo[, direction := fifelse(!is.na(padj) & padj < fdr & abs(log2FC_shrunk) >= lfc_cut,
                           fifelse(log2FC_shrunk > 0, "up", "down"), "ns")]
fwrite(geo, out[["de"]], sep = "\t")
saveRDS(dds, out[["dds_fit"]])
message("GEO DEGs (padj < ", fdr, ", |shrunk log2FC| >= ", lfc_cut, "): ",
        sum(geo$direction != "ns"), " [up ", sum(geo$direction == "up"),
        ", down ", sum(geo$direction == "down"), "]")

# ---- 3. gene-level concordance with TCGA -------------------------
tcga <- fread(inp[["tcga_de"]])[!is.na(gene_name) & gene_name != ""]
tcga <- tcga[order(padj)][!duplicated(gene_name)]
geo_u <- geo[!is.na(symbol)][order(padj)][!duplicated(symbol)]
cmp <- merge(tcga[, .(symbol = gene_name, tcga_lfc = log2FC_shrunk, tcga_padj = padj, tcga_dir = direction)],
             geo_u[, .(symbol, geo_lfc = log2FC_shrunk, geo_padj = padj, geo_stat = stat)],
             by = "symbol")
message("Genes shared with TCGA: ", nrow(cmp))

deg <- cmp[tcga_dir != "ns" & !is.na(geo_lfc)]
same_sign <- sign(deg$geo_lfc) == ifelse(deg$tcga_dir == "up", 1, -1)
rep_sig   <- same_sign & !is.na(deg$geo_padj) & deg$geo_padj < fdr
bt <- binom.test(sum(same_sign), length(same_sign), p = 0.5)
rho <- cor(cmp$tcga_lfc, cmp$geo_lfc, method = "spearman", use = "complete.obs")
rho_deg <- cor(deg$tcga_lfc, deg$geo_lfc, method = "spearman", use = "complete.obs")

cmp[, class := fifelse(tcga_dir == "ns", "TCGA ns",
               fifelse(!is.na(geo_padj) & geo_padj < fdr &
                         sign(geo_lfc) == sign(tcga_lfc), "replicated (GEO padj<0.05)", "TCGA DEG only"))]
lab <- head(cmp[class == "replicated (GEO padj<0.05)"][order(-abs(tcga_lfc) * abs(geo_lfc))], 15)
g <- ggplot(cmp, aes(tcga_lfc, geo_lfc, colour = class)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_text_repel(data = lab, aes(label = symbol), size = 2.8, colour = "black", max.overlaps = 30) +
  scale_colour_manual(values = c(`TCGA ns` = "grey80", `TCGA DEG only` = "#E3A587",
                                 `replicated (GEO padj<0.05)` = "#6A3D9A"), name = NULL) +
  labs(x = "TCGA-PRAD log2FC (shrunk)", y = "GSE104131 log2FC (shrunk)",
       title = "Gene-level replication: TCGA-PRAD vs GSE104131",
       subtitle = sprintf("Spearman rho = %.2f (all shared genes), %.2f (TCGA DEGs); sign agreement %.0f%%",
                          rho, rho_deg, 100 * mean(same_sign))) +
  theme(legend.position = "bottom")
ggsave(out[["lfc_plot"]], g, width = 7, height = 7, dpi = 150)

# ---- 4. TCGA signature score in GEO ------------------------------
vsd <- vst(dds, blind = TRUE)
vm  <- assay(vsd); rownames(vm) <- mcols(dds)$symbol
vm  <- vm[!is.na(rownames(vm)) & !duplicated(rownames(vm)), ]
n_sig <- p[["sig_size"]]
up_sig   <- intersect(head(tcga[direction == "up"][order(padj)]$gene_name,   n_sig), rownames(vm))
down_sig <- intersect(head(tcga[direction == "down"][order(padj)]$gene_name, n_sig), rownames(vm))
z <- t(scale(t(vm)))
score <- colMeans(z[up_sig, , drop = FALSE]) - colMeans(z[down_sig, , drop = FALSE])
sc <- data.table(gsm = names(score), score = score)[m, on = "gsm"]
auc <- with(sc, {
  t <- score[condition == "tumor"]; n <- score[condition == "normal"]
  mean(outer(t, n, ">") + 0.5 * outer(t, n, "=="))
})
wide <- dcast(sc, pair + ancestry ~ condition, value.var = "score")
pair_frac <- mean(wide$tumor > wide$normal)
message(sprintf("TCGA signature (%d up, %d down genes): AUC = %.3f; tumor > normal in %d/%d pairs",
                length(up_sig), length(down_sig), auc, sum(wide$tumor > wide$normal), nrow(wide)))

g <- ggplot(sc, aes(condition, score)) +
  geom_line(aes(group = pair), colour = "grey70") +
  geom_point(aes(colour = ancestry), size = 2.6) +
  scale_colour_manual(values = c(AA = "#B5651D", EA = "#3A6EA5")) +
  labs(x = NULL, y = "TCGA signature score (z, up minus down)",
       title = "TCGA-derived signature in GSE104131",
       subtitle = sprintf("Top %d up / %d down TCGA DEGs; AUC = %.2f; lines join pairs",
                          length(up_sig), length(down_sig), auc))
ggsave(out[["score_plot"]], g, width = 5.5, height = 5.5, dpi = 150)

# ---- 5. Hallmark GSEA in GEO -------------------------------------
new_api <- "collection" %in% names(formals(msigdbr::msigdbr))
hdf <- if (new_api) msigdbr(species = "Homo sapiens", collection = "H") else
                    msigdbr(species = "Homo sapiens", category = "H")
hall <- split(hdf$gene_symbol, hdf$gs_name)
rk <- geo_u[!is.na(stat)]
ranks <- sort(setNames(rk$stat, rk$symbol), decreasing = TRUE)
fh <- fgsea(hall, ranks, minSize = 15, maxSize = 500, eps = 0, nproc = threads)[order(padj)]
fwrite(copy(fh)[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")],
       out[["hallmark"]], sep = "\t")

th <- fread(inp[["tcga_hallmark"]])
nes <- merge(th[, .(pathway, tcga_NES = NES, tcga_padj = padj)],
             fh[, .(pathway, geo_NES = NES, geo_padj = padj)], by = "pathway")
r_nes <- cor(nes$tcga_NES, nes$geo_NES, method = "pearson")
tsig  <- nes[tcga_padj < fdr]
nes_dir <- mean(sign(tsig$tcga_NES) == sign(tsig$geo_NES))
nes_rep <- sum(sign(tsig$tcga_NES) == sign(tsig$geo_NES) & tsig$geo_padj < fdr)
message(sprintf("Hallmark NES: Pearson r = %.2f; of %d TCGA-significant sets, %.0f%% same direction, %d significant in GEO",
                r_nes, nrow(tsig), 100 * nes_dir, nes_rep))

nes[, label := gsub("_", " ", sub("^HALLMARK_", "", pathway))]
nes[, status := fifelse(tcga_padj < fdr & geo_padj < fdr & sign(tcga_NES) == sign(geo_NES),
                        "significant in both", fifelse(tcga_padj < fdr, "TCGA only", "neither / GEO only"))]
g <- ggplot(nes, aes(tcga_NES, geo_NES, colour = status)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 2.4) +
  geom_text_repel(data = nes[status == "significant in both"], aes(label = label),
                  size = 2.6, colour = "black", max.overlaps = 40) +
  scale_colour_manual(values = c(`significant in both` = "#6A3D9A", `TCGA only` = "#E3A587",
                                 `neither / GEO only` = "grey70"), name = NULL) +
  labs(x = "TCGA-PRAD NES", y = "GSE104131 NES", title = "Hallmark pathway replication",
       subtitle = sprintf("Pearson r = %.2f; %d/%d TCGA-significant sets replicated", r_nes, nes_rep, nrow(tsig))) +
  theme(legend.position = "bottom")
ggsave(out[["nes_plot"]], g, width = 7.5, height = 7, dpi = 150)

# ---- summary -----------------------------------------------------
summ <- data.table(
  metric = c("complete_pairs", "pairs_AA", "pairs_EA", "geo_genes_tested", "geo_DEGs",
             "shared_genes", "tcga_DEGs_in_geo", "sign_agreement_tcga_DEGs", "sign_agreement_binom_p",
             "replicated_tcga_DEGs_geo_padj05", "spearman_lfc_all", "spearman_lfc_tcga_DEGs",
             "signature_auc", "signature_pairs_tumor_gt_normal",
             "hallmark_nes_pearson", "hallmark_tcga_sig_same_direction", "hallmark_replicated"),
  value = c(length(complete), uniqueN(m[ancestry == "AA", pair]), uniqueN(m[ancestry == "EA", pair]),
            nrow(geo), sum(geo$direction != "ns"), nrow(cmp), nrow(deg),
            round(mean(same_sign), 3), signif(bt$p.value, 3), sum(rep_sig),
            round(rho, 3), round(rho_deg, 3), round(auc, 3), round(pair_frac, 3),
            round(r_nes, 3), round(nes_dir, 3), nes_rep))
fwrite(summ, out[["summary"]], sep = "\t")
message("Validation summary:"); print(summ)
message("Done.")
