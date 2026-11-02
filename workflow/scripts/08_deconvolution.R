#!/usr/bin/env Rscript
# ---------------------------------------------------------------
# 08_deconvolution.R — tissue composition (M8, part 1)
#
# Bulk RNA-seq mixes tumor cells with stroma and immune cells.
# This step estimates composition and asks which findings survive
# adjusting for it.
#
#  1. Scores for all TCGA-PRAD samples (log2 TPM+1, gene symbols):
#       ESTIMATE stromal / immune (tidyestimate, Yoshihara 2013)
#       MCP-counter immune & stromal populations (Becht 2016)
#       Marker scores: smooth muscle, basal epithelium
#  2. Paired tumor vs normal composition (Wilcoxon signed-rank)
#  3. Composition-adjusted paired DE:
#       ~ patient + stromal_z + immune_z + condition  -> Hallmark GSEA
#     Which tumor-vs-normal pathways are composition-driven?
#  4. Black vs White tumors: composition differences, and whether the
#     interferon signal survives adjustment for stromal/immune content
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(BiocParallel)
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})
for (pkg in c("tidyestimate", "MCPcounter"))
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is missing. Run once: Rscript workflow/scripts/install_extras.R", call. = FALSE)

if (exists("snakemake")) {
  inp <- snakemake@input; out <- snakemake@output; p <- snakemake@params
  threads <- snakemake@threads
  con <- file(snakemake@log[[1]], open = "wt"); sink(con); sink(con, type = "message")
} else {
  inp <- list(se_full = "data/processed/prad_se_full.rds",
              se_paired = "data/processed/prad_se_paired.rds",
              dds = "data/processed/prad_dds_filtered.rds",
              hallmark = "results/gsea/fgsea_hallmark.tsv",
              race_sheet = "results/ancestry/tcga_race_sample_sheet.tsv",
              race_hallmark = "results/ancestry/tcga_race_hallmark.tsv")
  d <- "results/deconv"
  out <- list(scores = file.path(d, "cell_scores.tsv"),
              paired_tests = file.path(d, "paired_composition_tests.tsv"),
              paired_plot = file.path(d, "paired_composition.png"),
              adj_hallmark = file.path(d, "adjusted_paired_hallmark.tsv"),
              adj_plot = file.path(d, "adjusted_paired_hallmark.png"),
              race_comp = file.path(d, "race_composition.tsv"),
              race_plot = file.path(d, "race_composition.png"),
              race_adj_hallmark = file.path(d, "race_adjusted_hallmark.tsv"),
              summary = file.path(d, "deconv_summary.tsv"))
  p <- list(fdr = 0.05, min_count = 10, seed = 42)
  threads <- 4
}
for (f in unlist(out)) dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 11))
fdr <- p[["fdr"]]
set.seed(p[["seed"]])
bp <- if (threads > 1) MulticoreParam(threads) else SerialParam()
short <- function(x) gsub("_", " ", sub("^HALLMARK_", "", x))
zs <- function(x) as.numeric(scale(x))

new_api <- "collection" %in% names(formals(msigdbr::msigdbr))
hdf  <- if (new_api) msigdbr(species = "Homo sapiens", collection = "H") else
                     msigdbr(species = "Homo sapiens", category = "H")
hall <- split(hdf$gene_symbol, hdf$gs_name)
gsea_h <- function(stat, symbol) {
  d <- data.table(stat, symbol)[!is.na(stat) & !is.na(symbol) & symbol != ""]
  d <- d[order(-abs(stat))][!duplicated(symbol)]
  r <- sort(setNames(d$stat, d$symbol), decreasing = TRUE)
  fgsea(hall, r, minSize = 15, maxSize = 500, eps = 0, nproc = threads)[, .(pathway, NES, padj)]
}

# =================================================================
# 1. Composition scores for every sample
# =================================================================
se  <- readRDS(inp[["se_full"]])
rd  <- as.data.frame(rowData(se))
tpm <- assay(se, "tpm_unstrand")
sym <- rd$gene_name
keep <- !is.na(sym) & sym != ""
tpm <- tpm[keep, ]; sym <- sym[keep]
o   <- order(-rowMeans(tpm))                       # one row per symbol: highest-expressed ID
tpm <- tpm[o, ]; sym <- sym[o]
tpm <- tpm[!duplicated(sym), ]; rownames(tpm) <- sym[!duplicated(sym)]
expr <- log2(tpm + 1)
message("Expression matrix for scoring: ", nrow(expr), " symbols x ", ncol(expr), " samples")

# ESTIMATE (relative scores for RNA-seq; purity formula is Affymetrix-only, so not used)
est_in <- tidyestimate::filter_common_genes(as.data.frame(expr), id = "hgnc_symbol", tidy = FALSE,
                                            tell_missing = FALSE, find_alias = FALSE)
est <- as.data.table(tidyestimate::estimate_score(est_in, is_affymetrix = FALSE))
setnames(est, c("sample", "stromal", "immune", "estimate"), c("barcode", "estimate_stromal", "estimate_immune", "estimate_total"))

# MCP-counter
mcp <- MCPcounter::MCPcounter.estimate(expr, featuresType = "HUGO_symbols")
mcp <- as.data.table(t(mcp), keep.rownames = "barcode")
setnames(mcp, names(mcp)[-1], paste0("mcp_", gsub("[^A-Za-z0-9]+", "_", names(mcp)[-1])))

# Marker scores (mean z across all samples)
marker_score <- function(genes) {
  g <- intersect(genes, rownames(expr))
  colMeans(t(scale(t(expr[g, , drop = FALSE]))))
}
mk <- data.table(barcode = colnames(expr),
                 smooth_muscle = marker_score(c("ACTG2", "DES", "CNN1", "MYH11", "ACTA2", "TAGLN")),
                 basal_epithelium = marker_score(c("KRT5", "KRT15", "TP63")))

sc <- Reduce(function(a, b) merge(a, b, by = "barcode"), list(est, mcp, mk))
sc[, `:=`(patient = substr(barcode, 1, 12), code = substr(barcode, 14, 15))]
sc[, condition := fifelse(code == "11", "normal", fifelse(code == "01", "tumor", "other"))]
fwrite(sc, out[["scores"]], sep = "\t")
score_cols <- setdiff(names(sc), c("barcode", "patient", "code", "condition", "estimate_total"))
message("Scores: ", paste(score_cols, collapse = ", "))

# =================================================================
# 2. Paired composition: tumor vs normal
# =================================================================
sp  <- readRDS(inp[["se_paired"]])
pm  <- data.table(barcode = colnames(sp), patient = as.character(sp$patient),
                  condition = as.character(sp$condition))
ps  <- merge(pm, sc[, c("barcode", score_cols), with = FALSE], by = "barcode")
pl  <- melt(ps, id.vars = c("barcode", "patient", "condition"), variable.name = "score")
wide <- dcast(pl, patient + score ~ condition, value.var = "value")
pt <- wide[, .(median_normal = median(normal), median_tumor = median(tumor),
               tumor_higher_pairs = sum(tumor > normal), pairs = .N,
               p = wilcox.test(tumor, normal, paired = TRUE)$p.value), by = score]
pt[, padj := p.adjust(p, "BH")][order(p)]
fwrite(pt, out[["paired_tests"]], sep = "\t")
message("Paired composition, tumor vs normal (Wilcoxon signed-rank):")
print(pt[order(p), .(score, median_normal = round(median_normal, 2), median_tumor = round(median_tumor, 2),
                     tumor_higher = paste0(tumor_higher_pairs, "/", pairs), padj = signif(padj, 3))])

show <- c("estimate_stromal", "estimate_immune", "mcp_Fibroblasts", "smooth_muscle",
          "basal_epithelium", "mcp_Endothelial_cells", "mcp_T_cells", "mcp_Monocytic_lineage")
show <- intersect(show, score_cols)
pls <- pl[score %in% show]
pls[, score := factor(score, levels = intersect(show, unique(score)))]
g <- ggplot(pls, aes(condition, value)) +
  geom_line(aes(group = patient), colour = "grey75", linewidth = 0.3) +
  geom_point(aes(colour = condition), size = 1.2) +
  scale_colour_manual(values = c(normal = "#4C9F70", tumor = "#C8553D"), guide = "none") +
  facet_wrap(~ score, scales = "free_y", nrow = 2) +
  labs(x = NULL, y = "Score", title = "Tissue composition, tumor vs matched normal (52 pairs)",
       subtitle = "ESTIMATE, MCP-counter and marker scores; lines join each patient's samples")
ggsave(out[["paired_plot"]], g, width = 11, height = 6, dpi = 150)

# =================================================================
# 3. Composition-adjusted paired DE
# =================================================================
dds <- readRDS(inp[["dds"]])
cmp <- sc[match(colnames(dds), barcode)]
stopifnot(all(cmp$barcode == colnames(dds)))
dds$patient   <- factor(gsub("-", "_", as.character(dds$patient)))
dds$condition <- relevel(factor(dds$condition), ref = "normal")
dds$stromal_z <- zs(cmp$estimate_stromal)
dds$immune_z  <- zs(cmp$estimate_immune)
r_cs <- cor(as.numeric(dds$condition == "tumor"), dds$stromal_z)
r_ci <- cor(as.numeric(dds$condition == "tumor"), dds$immune_z)
message(sprintf("Correlation of tumor status with stromal score r = %.2f, immune score r = %.2f", r_cs, r_ci))

design(dds) <- ~ patient + stromal_z + immune_z + condition
dds <- DESeq(dds, parallel = threads > 1, BPPARAM = bp, quiet = TRUE)
ra  <- results(dds, name = "condition_tumor_vs_normal", alpha = fdr)
gn  <- mcols(dds)$gene_name
n_adj_any <- sum(ra$padj < fdr, na.rm = TRUE)
n_adj_deg <- sum(ra$padj < fdr & abs(ra$log2FoldChange) >= 1, na.rm = TRUE)
message(sprintf("Composition-adjusted paired DE: %d genes padj < %.2f; %d also |log2FC| >= 1 (unshrunk)",
                n_adj_any, fdr, n_adj_deg))

ha <- gsea_h(ra$stat, gn)
hm <- fread(inp[["hallmark"]])[, .(pathway, NES_main = NES, padj_main = padj)]
ah <- merge(hm, ha[, .(pathway, NES_adj = NES, padj_adj = padj)], by = "pathway")
ah[, class := fifelse(padj_main < fdr & padj_adj < fdr & sign(NES_main) == sign(NES_adj), "robust to composition",
              fifelse(padj_main < fdr, "composition-sensitive", fifelse(padj_adj < fdr, "emerges after adjustment", "neither")))]
ah <- ah[order(padj_main)]
fwrite(ah, out[["adj_hallmark"]], sep = "\t")
message("Hallmark tumor-vs-normal sets after composition adjustment:")
print(ah[padj_main < fdr | padj_adj < fdr][order(class, -NES_main),
         .(pathway = short(pathway), NES_main = round(NES_main, 2), NES_adj = round(NES_adj, 2),
           padj_adj = signif(padj_adj, 3), class)])

g <- ggplot(ah, aes(NES_main, NES_adj, colour = class)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 2.6) +
  ggrepel::geom_text_repel(data = ah[class != "neither"], aes(label = short(pathway)),
                           size = 2.6, colour = "black", max.overlaps = 40) +
  scale_colour_manual(values = c(`robust to composition` = "#6A3D9A", `composition-sensitive` = "#E3A587",
                                 `emerges after adjustment` = "#2E6F9E", neither = "grey75"), name = NULL) +
  labs(x = "NES, main paired model", y = "NES, adjusted for ESTIMATE stromal + immune scores",
       title = "Which tumor-vs-normal pathways reflect tissue composition?",
       subtitle = "~ patient + stromal + immune + condition") +
  theme(legend.position = "bottom")
ggsave(out[["adj_plot"]], g, width = 8, height = 7.5, dpi = 150)

# =================================================================
# 4. Black vs White tumors
# =================================================================
tb <- fread(inp[["race_sheet"]])
tb[, race := factor(race, levels = c("White", "Black"))]
tb[, gleason_group := factor(gleason_group, levels = c("GS6", "GS7", "GS8_10"))]
tr <- merge(tb, sc[, c("barcode", score_cols), with = FALSE], by = "barcode")

rc <- rbindlist(lapply(score_cols, function(s) {
  fit <- lm(zs(tr[[s]]) ~ race + gleason_group + age_s, data = tr)
  co  <- summary(fit)$coefficients["raceBlack", ]
  data.table(score = s, diff_sd_black_vs_white = co[["Estimate"]], se = co[["Std. Error"]], p = co[["Pr(>|t|)"]])
}))
rc[, padj := p.adjust(p, "BH")]
rc <- rc[order(p)]
fwrite(rc, out[["race_comp"]], sep = "\t")
message("Composition, Black vs White tumors (SD units; adjusted for Gleason group and age):")
print(rc[, .(score, diff_sd = round(diff_sd_black_vs_white, 2), padj = signif(padj, 3))])

rl <- melt(tr[, c("barcode", "race", show), with = FALSE], id.vars = c("barcode", "race"), variable.name = "score")
rl[, score := factor(score, levels = show)]
rl[, label := sprintf("%s\nadj p = %s", score, signif(rc$padj[match(score, rc$score)], 2))]
rl[, label := factor(label, levels = unique(label[order(score)]))]
g <- ggplot(rl, aes(race, value, fill = race)) +
  geom_boxplot(outlier.size = 0.6, width = 0.6) +
  scale_fill_manual(values = c(White = "#3A6EA5", Black = "#B5651D"), guide = "none") +
  facet_wrap(~ label, scales = "free_y", nrow = 2) +
  labs(x = NULL, y = "Score", title = "Tumor composition by self-reported race (TCGA-PRAD)",
       subtitle = sprintf("%d Black, %d White tumors; p from linear model adjusted for Gleason group and age",
                          sum(tr$race == "Black"), sum(tr$race == "White"))) +
  theme(strip.text = element_text(size = 8))
ggsave(out[["race_plot"]], g, width = 11, height = 6.5, dpi = 150)

# Does the interferon signal survive composition adjustment?
ct  <- assay(se, "unstranded")[, tr$barcode]
cdt <- data.frame(tr[, .(race, gleason_group = droplevels(gleason_group), age_s,
                         stromal_z = zs(estimate_stromal), immune_z = zs(estimate_immune))],
                  row.names = tr$barcode)
dr <- DESeqDataSetFromMatrix(ct, cdt, design = ~ gleason_group + age_s + stromal_z + immune_z + race)
dr <- dr[rowSums(counts(dr) >= p[["min_count"]]) >= sum(tr$race == "Black"), ]
dr <- DESeq(dr, parallel = threads > 1, BPPARAM = bp, quiet = TRUE)
rr <- results(dr, name = "race_Black_vs_White", alpha = fdr)
gt <- rd[rownames(dr), "gene_type"]
rep_ok <- gt %in% c("protein_coding", "lncRNA")
hr <- gsea_h(rr$stat[rep_ok], rd[rownames(dr), "gene_name"][rep_ok])
mr <- fread(inp[["race_hallmark"]])[, .(pathway, NES_main = NES, padj_main = padj)]
ra2 <- merge(mr, hr[, .(pathway, NES_adj = NES, padj_adj = padj)], by = "pathway")[order(padj_main)]
fwrite(ra2, out[["race_adj_hallmark"]], sep = "\t")
message("Race Hallmark sets before vs after stromal + immune adjustment:")
print(ra2[padj_main < fdr | padj_adj < fdr, .(pathway = short(pathway), NES_main = round(NES_main, 2),
          padj_main = signif(padj_main, 3), NES_adj = round(NES_adj, 2), padj_adj = signif(padj_adj, 3))])
ifn <- ra2[grepl("INTERFERON", pathway)]

# ---- summary -----------------------------------------------------
summ <- data.table(
  metric = c("paired_stromal_tumor_lower_pairs", "paired_immune_padj", "cor_tumor_stromal", "cor_tumor_immune",
             "adj_paired_genes_padj05", "adj_paired_genes_padj05_lfc1",
             "hallmark_robust", "hallmark_composition_sensitive", "hallmark_emerging",
             "race_immune_diff_sd", "race_immune_padj",
             "race_ifna_nes_adj", "race_ifna_padj_adj", "race_ifng_nes_adj", "race_ifng_padj_adj"),
  value = c(pt[score == "estimate_stromal", pairs - tumor_higher_pairs], signif(pt[score == "estimate_immune", padj], 3),
            round(r_cs, 3), round(r_ci, 3), n_adj_any, n_adj_deg,
            sum(ah$class == "robust to composition"), sum(ah$class == "composition-sensitive"),
            sum(ah$class == "emerges after adjustment"),
            round(rc[score == "estimate_immune", diff_sd_black_vs_white], 3), signif(rc[score == "estimate_immune", padj], 3),
            round(ifn[pathway == "HALLMARK_INTERFERON_ALPHA_RESPONSE", NES_adj], 3),
            signif(ifn[pathway == "HALLMARK_INTERFERON_ALPHA_RESPONSE", padj_adj], 3),
            round(ifn[pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE", NES_adj], 3),
            signif(ifn[pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE", padj_adj], 3)))
fwrite(summ, out[["summary"]], sep = "\t")
message("Deconvolution summary:"); print(summ)
message("Done.")
