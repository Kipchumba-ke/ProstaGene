#!/usr/bin/env Rscript
# ---------------------------------------------------------------
# 03_qc.R — QC and exploratory analysis of the paired TCGA-PRAD set
#
#  1. Library size / detected-gene summary
#  2. Low-count gene filtering
#  3. VST (blind) -> PCA on top variable genes, coloured by
#     condition, tissue source site (TSS) and plate
#  4. Sample-to-sample distance heatmap
#  5. Tumor purity sanity check: tumors that sit with the normals,
#     by PCA position and by a tumor-vs-basal marker score
#
# Outputs a filtered DESeqDataSet + VST object for M3, figures,
# and a per-sample flag table. Nothing is removed automatically —
# flagged samples are reviewed before DE.
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})

# ---- I/O: Snakemake object if present, else defaults -------------
if (exists("snakemake")) {
  in_se  <- snakemake@input[["se_paired"]]
  out    <- snakemake@output
  p      <- snakemake@params
  assay_name <- snakemake@config[["count_assay"]]
  con <- file(snakemake@log[[1]], open = "wt"); sink(con); sink(con, type = "message")
} else {
  in_se <- "data/processed/prad_se_paired.rds"
  qc <- "results/qc"
  out <- list(
    dds = "data/processed/prad_dds_filtered.rds",
    vsd = "data/processed/prad_vst.rds",
    flags = file.path(qc, "qc_sample_flags.tsv"),
    libsize = file.path(qc, "library_sizes.png"),
    pca_condition = file.path(qc, "pca_condition.png"),
    pca_tss = file.path(qc, "pca_tss.png"),
    pca_plate = file.path(qc, "pca_plate.png"),
    heatmap = file.path(qc, "sample_distance_heatmap.png"),
    markers = file.path(qc, "purity_marker_score.png")
  )
  p <- list(min_count = 10, n_top_genes = 500)
  assay_name <- "unstranded"
}
for (f in unlist(out)) dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)

theme_set(theme_bw(base_size = 11))
cond_cols <- c(normal = "#4C9F70", tumor = "#C8553D")

# ---- build DESeqDataSet ------------------------------------------
se <- readRDS(in_se)
cd <- as.data.frame(colData(se))[, c("patient", "condition", "tss", "plate")]
cd$condition <- factor(cd$condition, levels = c("normal", "tumor"))
cd$patient   <- droplevels(factor(cd$patient))
cd$tss       <- factor(cd$tss)
cd$plate     <- factor(cd$plate)

dds <- DESeqDataSetFromMatrix(
  countData = assay(se, assay_name),
  colData   = cd,
  design    = ~ patient + condition
)
rd_cols <- intersect(c("gene_id", "gene_name", "gene_type"), colnames(rowData(se)))
mcols(dds) <- cbind(mcols(dds), rowData(se)[, rd_cols, drop = FALSE])

# ---- 1. library sizes --------------------------------------------
lib <- data.frame(
  barcode   = colnames(dds),
  condition = cd$condition,
  lib_size  = colSums(counts(dds)),
  detected  = colSums(counts(dds) > 0)
)
message("Library size (millions): median ",
        round(median(lib$lib_size) / 1e6, 1), ", range ",
        paste(round(range(lib$lib_size) / 1e6, 1), collapse = "–"))

g <- ggplot(lib, aes(reorder(barcode, lib_size), lib_size / 1e6, fill = condition)) +
  geom_col() + scale_fill_manual(values = cond_cols) +
  labs(x = NULL, y = "Assigned reads (millions)", title = "Library sizes") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(out[["libsize"]], g, width = 9, height = 4, dpi = 150)

# ---- 2. filter low-count genes -----------------------------------
# Keep genes with >= min_count reads in at least as many samples as
# the smaller condition group (52 here), so genes expressed only in
# tumors or only in normals are retained.
min_n <- min(table(cd$condition))
keep  <- rowSums(counts(dds) >= p[["min_count"]]) >= min_n
message("Genes before filtering: ", nrow(dds), " | after: ", sum(keep))
if ("gene_type" %in% colnames(mcols(dds))) {
  message("Gene types retained:")
  print(head(sort(table(mcols(dds)$gene_type[keep]), decreasing = TRUE), 6))
}
dds <- dds[keep, ]
dds <- estimateSizeFactors(dds)

# ---- 3. VST + PCA ------------------------------------------------
vsd <- vst(dds, blind = TRUE)
mat <- assay(vsd)

ntop <- min(p[["n_top_genes"]], nrow(mat))
top  <- order(MatrixGenerics::rowVars(mat), decreasing = TRUE)[seq_len(ntop)]
pca  <- prcomp(t(mat[top, ]))
pve  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
message("PCA variance explained (PC1–PC5): ", paste0(pve[1:5], "%", collapse = ", "))

pc <- cbind(cd, as.data.frame(pca$x[, 1:3]))
pc$barcode <- rownames(pc)

pca_plot <- function(colour_by, title, legend = "right") {
  ggplot(pc, aes(PC1, PC2)) +
    geom_line(aes(group = patient), colour = "grey75", linewidth = 0.3) +  # links each pair
    geom_point(aes(colour = .data[[colour_by]], shape = condition), size = 2.4) +
    labs(x = paste0("PC1 (", pve[1], "%)"), y = paste0("PC2 (", pve[2], "%)"),
         title = title, subtitle = paste0("VST, top ", ntop, " variable genes; grey lines join patient pairs")) +
    theme(legend.position = legend)
}
ggsave(out[["pca_condition"]],
       pca_plot("condition", "PCA by condition") + scale_colour_manual(values = cond_cols),
       width = 7, height = 5.5, dpi = 150)
ggsave(out[["pca_tss"]],   pca_plot("tss",   "PCA by tissue source site", "bottom"),
       width = 8, height = 7, dpi = 150)
ggsave(out[["pca_plate"]], pca_plot("plate", "PCA by sequencing plate", "bottom"),
       width = 8, height = 7, dpi = 150)

# Rough check: does PC1 separate condition, and is TSS/plate associated
# with PC1/PC2? (Kruskal–Wallis p-values; descriptive only.)
assoc <- sapply(c("condition", "tss", "plate"), function(v)
  sapply(c("PC1", "PC2"), function(k) signif(kruskal.test(pc[[k]] ~ pc[[v]])$p.value, 3)))
message("Kruskal–Wallis p-values, PC vs covariate:"); print(assoc)

# ---- 4. sample distance heatmap ----------------------------------
d  <- dist(t(mat))
dm <- as.matrix(d)
ann <- data.frame(condition = cd$condition, row.names = colnames(mat))
pheatmap(dm,
         clustering_distance_rows = d, clustering_distance_cols = d,
         annotation_col = ann, annotation_row = ann,
         annotation_colors = list(condition = cond_cols),
         show_rownames = FALSE, show_colnames = FALSE,
         color = colorRampPalette(c("#08306B", "#6BAED6", "white"))(100),
         main = "Sample-to-sample Euclidean distance (VST)",
         filename = out[["heatmap"]], width = 9, height = 8)

# ---- 5. purity sanity check --------------------------------------
# (a) PCA: is each tumor nearer the normal or the tumor centroid (PC1–PC2)?
cent <- aggregate(cbind(PC1, PC2) ~ condition, data = pc, FUN = median)
dist_to <- function(lbl) sqrt((pc$PC1 - cent$PC1[cent$condition == lbl])^2 +
                              (pc$PC2 - cent$PC2[cent$condition == lbl])^2)
pc$nearest_centroid <- ifelse(dist_to("tumor") < dist_to("normal"), "tumor", "normal")

# (b) Marker score: tumor-up genes minus normal basal-cell genes
#     (basal cells are lost in prostate adenocarcinoma).
tumor_markers  <- c("AMACR", "PCA3", "HPN", "DLX1", "TDRD1")
normal_markers <- c("KRT5", "KRT15", "TP63", "MYLK", "CNN1")
gn <- if ("gene_name" %in% colnames(mcols(vsd))) mcols(vsd)$gene_name else rownames(vsd)
score_set <- function(genes) {
  idx <- which(gn %in% genes)
  if (!length(idx)) return(rep(NA_real_, ncol(mat)))
  message("  markers found: ", paste(unique(gn[idx]), collapse = ", "))
  colMeans(t(scale(t(mat[idx, , drop = FALSE]))))
}
message("Tumor markers:");  s_t <- score_set(tumor_markers)
message("Normal markers:"); s_n <- score_set(normal_markers)
pc$marker_score <- s_t - s_n

normal_q90 <- quantile(pc$marker_score[pc$condition == "normal"], 0.9, na.rm = TRUE)
pc$flag_pca    <- pc$condition == "tumor" & pc$nearest_centroid == "normal"
pc$flag_marker <- pc$condition == "tumor" & pc$marker_score < normal_q90
pc$flagged     <- pc$flag_pca | pc$flag_marker

g <- ggplot(pc, aes(condition, marker_score)) +
  geom_line(aes(group = patient), colour = "grey80", linewidth = 0.3) +
  geom_point(aes(colour = flagged), size = 2) +
  geom_hline(yintercept = normal_q90, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(values = c(`FALSE` = "grey30", `TRUE` = "#D62728")) +
  labs(x = NULL, y = "Tumor-marker minus basal-marker score (z)",
       title = "Purity sanity check",
       subtitle = "Dashed line: 90th percentile of normals; red = flagged tumor")
ggsave(out[["markers"]], g, width = 6, height = 5.5, dpi = 150)

flags <- merge(lib, pc[, c("barcode", "patient", "tss", "plate", "PC1", "PC2", "PC3",
                           "nearest_centroid", "marker_score",
                           "flag_pca", "flag_marker", "flagged")],
               by = "barcode")
flags <- flags[order(!flags$flagged, flags$patient), ]
write.table(flags, out[["flags"]], sep = "\t", quote = FALSE, row.names = FALSE)

message("Flagged tumors: ", sum(flags$flagged))
if (any(flags$flagged))
  print(flags[flags$flagged, c("barcode", "nearest_centroid", "marker_score")], row.names = FALSE)

# ---- save for M3 -------------------------------------------------
saveRDS(dds, out[["dds"]])
saveRDS(vsd, out[["vsd"]])
message("Done.")
