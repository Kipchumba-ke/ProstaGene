#!/usr/bin/env Rscript
# ---------------------------------------------------------------
# 04_deseq2.R — paired differential expression, tumor vs normal
#
#  Main:        all 52 pairs, design ~ patient + condition
#  Shrinkage:   apeglm LFC shrinkage for ranking / plotting
#  Sensitivity: refit without patients whose tumor was flagged in QC,
#               compare fold changes and DEG overlap
#  Purity:      published purity estimates for the flagged tumors
#               (TCGA values use decimal commas, e.g. "0,83")
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(BiocParallel)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

if (exists("snakemake")) {
  inp <- snakemake@input; out <- snakemake@output; p <- snakemake@params
  threads <- snakemake@threads
  con <- file(snakemake@log[[1]], open = "wt"); sink(con); sink(con, type = "message")
} else {
  inp <- list(dds = "data/processed/prad_dds_filtered.rds",
              flags = "results/qc/qc_sample_flags.tsv",
              se_paired = "data/processed/prad_se_paired.rds")
  de <- "results/de"
  out <- list(dds_fit = "data/processed/prad_dds_fitted.rds",
              table_all = file.path(de, "deseq2_all_genes.tsv"),
              table_sig = file.path(de, "deseq2_significant.tsv"),
              volcano = file.path(de, "volcano.png"),
              ma = file.path(de, "ma_plot.png"),
              heatmap = file.path(de, "top50_heatmap.png"),
              sens_plot = file.path(de, "sensitivity_concordance.png"),
              sens_summary = file.path(de, "sensitivity_summary.tsv"),
              purity = file.path(de, "purity_check.tsv"))
  p <- list(fdr = 0.05, lfc = 1)
  threads <- 4
}
for (f in unlist(out)) dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
bp <- if (threads > 1) MulticoreParam(threads) else SerialParam()
theme_set(theme_bw(base_size = 11))
fdr <- p[["fdr"]]; lfc_cut <- p[["lfc"]]
coef_name <- "condition_tumor_vs_normal"

# ---- prepare -----------------------------------------------------
dds <- readRDS(inp[["dds"]])
# Patient IDs: replace '-' so factor levels are syntactic (silences the DESeq2 note)
dds$patient   <- factor(gsub("-", "_", as.character(dds$patient)))
dds$condition <- relevel(factor(dds$condition), ref = "normal")
design(dds)   <- ~ patient + condition

gene_info <- as.data.frame(mcols(dds)[, intersect(c("gene_name", "gene_type"),
                                                  colnames(mcols(dds))), drop = FALSE])
gene_info$gene_id <- rownames(dds)

fit_and_extract <- function(d, shrink = TRUE) {
  d <- DESeq(d, parallel = threads > 1, BPPARAM = bp, quiet = TRUE)
  r <- results(d, name = coef_name, alpha = fdr)
  tab <- data.frame(gene_id = rownames(r), baseMean = r$baseMean,
                    log2FC = r$log2FoldChange, lfcSE = r$lfcSE,
                    stat = r$stat, pvalue = r$pvalue, padj = r$padj)
  if (shrink) {
    s <- lfcShrink(d, coef = coef_name, type = "apeglm",
                   parallel = threads > 1, BPPARAM = bp, quiet = TRUE)
    tab$log2FC_shrunk <- s$log2FoldChange
  }
  list(dds = d, tab = tab)
}

# ---- main fit: all pairs -----------------------------------------
message("Main fit: ", nlevels(dds$patient), " pairs, ", nrow(dds), " genes")
main <- fit_and_extract(dds, shrink = TRUE)
tab  <- merge(gene_info, main$tab, by = "gene_id")
tab  <- tab[order(tab$padj, -abs(tab$log2FC_shrunk), na.last = TRUE), ]

tab$direction <- "ns"
is_sig <- !is.na(tab$padj) & tab$padj < fdr & abs(tab$log2FC_shrunk) >= lfc_cut
tab$direction[is_sig & tab$log2FC_shrunk > 0] <- "up"
tab$direction[is_sig & tab$log2FC_shrunk < 0] <- "down"

message("Significant (padj < ", fdr, ", |shrunk log2FC| >= ", lfc_cut, "): ",
        sum(is_sig), "  [up ", sum(tab$direction == "up"),
        ", down ", sum(tab$direction == "down"), "]")
message("padj < ", fdr, " at any fold change: ", sum(tab$padj < fdr, na.rm = TRUE))

write.table(tab, out[["table_all"]], sep = "\t", quote = FALSE, row.names = FALSE)
write.table(tab[is_sig, ], out[["table_sig"]], sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(main$dds, out[["dds_fit"]])

lab <- if ("gene_name" %in% names(tab)) tab$gene_name else tab$gene_id
message("Top 15 by padj:")
print(head(data.frame(gene = lab, log2FC_shrunk = round(tab$log2FC_shrunk, 2),
                      padj = signif(tab$padj, 3)), 15), row.names = FALSE)

# ---- volcano -----------------------------------------------------
dir_cols <- c(up = "#C8553D", down = "#2E6F9E", ns = "grey75")
v <- tab[!is.na(tab$padj), ]
v$label <- if ("gene_name" %in% names(v)) v$gene_name else v$gene_id
min_p <- min(v$padj[v$padj > 0])
v$neglog10 <- -log10(pmax(v$padj, min_p))
sv <- v[v$direction != "ns", ]
top_lab <- head(sv[order(-abs(sv$log2FC_shrunk) * sv$neglog10), ], 20)
g <- ggplot(v, aes(log2FC_shrunk, neglog10, colour = direction)) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = -log10(fdr), linetype = "dashed", colour = "grey50") +
  geom_text_repel(data = top_lab, aes(label = label), size = 3, colour = "black",
                  max.overlaps = 30, min.segment.length = 0) +
  scale_colour_manual(values = dir_cols) +
  labs(x = "log2 fold change, tumor vs normal (apeglm-shrunk)",
       y = expression(-log[10]~adjusted~p),
       title = "TCGA-PRAD: tumor vs matched normal",
       subtitle = paste0(nlevels(dds$patient), " pairs; ", sum(v$direction == "up"),
                         " up, ", sum(v$direction == "down"), " down"))
ggsave(out[["volcano"]], g, width = 7.5, height = 6.5, dpi = 150)

# ---- MA plot -----------------------------------------------------
g <- ggplot(v, aes(log10(baseMean), log2FC_shrunk, colour = direction)) +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_hline(yintercept = 0, colour = "grey30") +
  scale_colour_manual(values = dir_cols) +
  labs(x = "log10 mean normalised count", y = "log2 fold change (shrunk)", title = "MA plot")
ggsave(out[["ma"]], g, width = 7, height = 5, dpi = 150)

# ---- top-50 heatmap (VST, design-aware) --------------------------
vsd <- vst(main$dds, blind = FALSE)
top <- head(tab$gene_id[is_sig], 50)
if (length(top) >= 2) {
  m <- assay(vsd)[top, ]
  m <- m - rowMeans(m)
  nm <- lab[match(top, tab$gene_id)]
  nm[is.na(nm) | nm == ""] <- top[is.na(nm) | nm == ""]
  rownames(m) <- make.unique(nm)
  ann <- data.frame(condition = vsd$condition, row.names = colnames(m))
  pheatmap(m, annotation_col = ann,
           annotation_colors = list(condition = c(normal = "#4C9F70", tumor = "#C8553D")),
           show_colnames = FALSE, fontsize_row = 7,
           color = colorRampPalette(c("#2E6F9E", "white", "#C8553D"))(101),
           breaks = seq(-max(abs(m)), max(abs(m)), length.out = 102),
           main = "Top 50 DEGs (VST, row-centred)",
           filename = out[["heatmap"]], width = 9, height = 10)
}

# ---- purity of tumors (decimal-comma values) ---------------------
se <- readRDS(inp[["se_paired"]])
cd <- as.data.frame(colData(se))
to_num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
pur_cols <- intersect(c("paper_ABSOLUTE_purity", "paper_avgRNA_purity",
                        "paper_avgDNA_purity", "paper_Tumor_cellularity_pathology"),
                      colnames(cd))
flags <- read.delim(inp[["flags"]])
flagged_bc <- flags$barcode[flags$flagged %in% c(TRUE, "TRUE")]
pur <- data.frame(barcode = rownames(cd)[cd$condition == "tumor"])
for (k in pur_cols) pur[[k]] <- to_num(cd[pur$barcode, k])
pur$qc_flagged <- pur$barcode %in% flagged_bc
write.table(pur, out[["purity"]], sep = "\t", quote = FALSE, row.names = FALSE)
message("Tumor purity — flagged tumors:")
print(pur[pur$qc_flagged, ], row.names = FALSE)
message("Tumor purity — all paired tumors (summary):")
print(summary(pur[, pur_cols, drop = FALSE]))

# ---- sensitivity: drop pairs of flagged tumors -------------------
drop_pat <- unique(gsub("-", "_", substr(flagged_bc, 1, 12)))
if (length(drop_pat)) {
  d2 <- dds[, !(dds$patient %in% drop_pat)]
  d2$patient <- droplevels(d2$patient)
  message("Sensitivity fit: dropping ", length(drop_pat), " pairs -> ",
          nlevels(d2$patient), " pairs")
  sens <- fit_and_extract(d2, shrink = FALSE)$tab

  cmp <- merge(main$tab[, c("gene_id", "log2FC", "padj")],
               sens[, c("gene_id", "log2FC", "padj")],
               by = "gene_id", suffixes = c("_main", "_sens"))
  sig_m <- with(cmp, !is.na(padj_main) & padj_main < fdr & abs(log2FC_main) >= lfc_cut)
  sig_s <- with(cmp, !is.na(padj_sens) & padj_sens < fdr & abs(log2FC_sens) >= lfc_cut)
  both  <- sig_m & sig_s
  summ <- data.frame(
    pairs_main        = nlevels(dds$patient),
    pairs_sensitivity = nlevels(d2$patient),
    sig_main          = sum(sig_m),
    sig_sensitivity   = sum(sig_s),
    sig_both          = sum(both),
    jaccard           = round(sum(both) / sum(sig_m | sig_s), 3),
    same_direction_in_overlap = round(mean(sign(cmp$log2FC_main[both]) ==
                                           sign(cmp$log2FC_sens[both])), 3),
    spearman_log2FC_all = round(cor(cmp$log2FC_main, cmp$log2FC_sens,
                                    method = "spearman", use = "complete.obs"), 3)
  )
  message("Sensitivity summary (MLE log2FC, same thresholds):"); print(summ, row.names = FALSE)
  write.table(summ, out[["sens_summary"]], sep = "\t", quote = FALSE, row.names = FALSE)

  cmp$class <- ifelse(both, "both", ifelse(sig_m, "main only", ifelse(sig_s, "sensitivity only", "neither")))
  g <- ggplot(cmp, aes(log2FC_main, log2FC_sens, colour = class)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey40") +
    geom_point(size = 0.6, alpha = 0.6) +
    scale_colour_manual(values = c(both = "#6A3D9A", `main only` = "#C8553D",
                                   `sensitivity only` = "#2E6F9E", neither = "grey80")) +
    coord_equal() +
    labs(x = paste0("log2FC, all ", summ$pairs_main, " pairs"),
         y = paste0("log2FC, ", summ$pairs_sensitivity, " pairs"),
         title = "Sensitivity: excluding QC-borderline tumors",
         subtitle = paste0("Spearman rho = ", summ$spearman_log2FC_all,
                           "; DEG Jaccard = ", summ$jaccard))
  ggsave(out[["sens_plot"]], g, width = 6.5, height = 6.5, dpi = 150)
} else {
  message("No flagged tumors; sensitivity fit skipped.")
  write.table(data.frame(note = "no flagged tumors"), out[["sens_summary"]],
              sep = "\t", quote = FALSE, row.names = FALSE)
  file.create(out[["sens_plot"]])
}

writeLines(capture.output(sessionInfo()), sub("\\.tsv$", "_sessionInfo.txt", out[["table_all"]]))
message("Done.")
