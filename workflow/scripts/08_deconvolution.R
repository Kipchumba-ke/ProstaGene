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
#  3. Composition-adjusted paired DE -> Hallmark GSEA, two models:
#       ~ patient + stromal_z + immune_z + condition            (ESTIMATE)
#       ~ patient + stromal_z + immune_z + smooth_z + basal_z + condition
#     Which tumor-vs-normal pathways are composition-driven?
#  4. Black vs White tumors: composition differences, and whether the
#     interferon signal survives adjustment, using ESTIMATE scores and,
#     separately, MCP-counter scores built from lineage marker genes
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
dds$smooth_z  <- zs(cmp$smooth_muscle)
dds$basal_z   <- zs(cmp$basal_epithelium)
is_tum <- as.numeric(dds$condition == "tumor")
r_cs <- cor(is_tum, dds$stromal_z); r_ci <- cor(is_tum, dds$immune_z)
r_sm <- cor(is_tum, dds$smooth_z);  r_ba <- cor(is_tum, dds$basal_z)
message(sprintf("Correlation of tumor status with scores: ESTIMATE stromal r = %.2f, immune r = %.2f, smooth muscle r = %.2f, basal r = %.2f",
                r_cs, r_ci, r_sm, r_ba))
# ESTIMATE's stromal signature is fibroblast/ECM-based; prostate stroma is mostly smooth muscle,
# so a second model also adjusts for smooth-muscle and basal-cell scores.

gn <- mcols(dds)$gene_name
fit_paired <- function(f) {
  d <- dds; design(d) <- f
  d <- DESeq(d, parallel = threads > 1, BPPARAM = bp, quiet = TRUE)
  results(d, name = "condition_tumor_vs_normal", alpha = fdr)
}
ra  <- fit_paired(~ patient + stromal_z + immune_z + condition)
ra2 <- fit_paired(~ patient + stromal_z + immune_z + smooth_z + basal_z + condition)
n_adj_any <- sum(ra$padj < fdr, na.rm = TRUE)
n_adj_deg <- sum(ra$padj < fdr & abs(ra$log2FoldChange) >= 1, na.rm = TRUE)
n_full_any <- sum(ra2$padj < fdr, na.rm = TRUE)
message(sprintf("Paired DE adjusted for ESTIMATE: %d genes padj < %.2f (%d with |log2FC| >= 1); + smooth muscle + basal: %d genes padj < %.2f",
                n_adj_any, fdr, n_adj_deg, n_full_any, fdr))

classify <- function(nm, pm, na, pa) fifelse(pm < fdr & pa < fdr & sign(nm) == sign(na), "robust to composition",
                                     fifelse(pm < fdr, "composition-sensitive",
                                     fifelse(pa < fdr, "emerges after adjustment", "neither")))
hm <- fread(inp[["hallmark"]])[, .(pathway, NES_main = NES, padj_main = padj)]
ah <- Reduce(function(x, y) merge(x, y, by = "pathway"),
             list(hm, gsea_h(ra$stat, gn)[, .(pathway, NES_adj = NES, padj_adj = padj)],
                  gsea_h(ra2$stat, gn)[, .(pathway, NES_full = NES, padj_full = padj)]))
ah[, class := classify(NES_main, padj_main, NES_adj, padj_adj)]
ah[, class_full := classify(NES_main, padj_main, NES_full, padj_full)]
ah <- ah[order(padj_main)]
fwrite(ah, out[["adj_hallmark"]], sep = "\t")
message("Hallmark tumor-vs-normal sets: main vs ESTIMATE-adjusted vs fully adjusted (+ smooth muscle + basal):")
print(ah[padj_main < fdr | padj_adj < fdr | padj_full < fdr][order(class_full, -NES_main),
         .(pathway = short(pathway), NES_main = round(NES_main, 2), NES_adj = round(NES_adj, 2),
           NES_full = round(NES_full, 2), padj_full = signif(padj_full, 3), class, class_full)])
message("Note: the smooth-muscle score shares genes with HALLMARK_MYOGENESIS, so its loss under full adjustment is partly circular.")

al <- rbind(ah[, .(pathway, NES_main, NES = NES_adj, cls = class, model = "Adjusted: ESTIMATE stromal + immune")],
            ah[, .(pathway, NES_main, NES = NES_full, cls = class_full, model = "Adjusted: + smooth muscle + basal")])
g <- ggplot(al, aes(NES_main, NES, colour = cls)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 2.2) +
  ggrepel::geom_text_repel(data = al[cls %in% c("composition-sensitive", "emerges after adjustment")],
                           aes(label = short(pathway)), size = 2.4, colour = "black", max.overlaps = 40) +
  scale_colour_manual(values = c(`robust to composition` = "#6A3D9A", `composition-sensitive` = "#E3A587",
                                 `emerges after adjustment` = "#2E6F9E", neither = "grey75"), name = NULL) +
  facet_wrap(~ model) +
  labs(x = "NES, main paired model", y = "NES, adjusted model",
       title = "Which tumor-vs-normal pathways reflect tissue composition?",
       subtitle = "Labelled: sets whose significance changes with adjustment") +
  theme(legend.position = "bottom")
ggsave(out[["adj_plot"]], g, width = 12, height = 6.5, dpi = 150)

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
# Two adjustments, because they can remove different things:
#   ESTIMATE: stromal + immune scores. The immune signature shares genes with the
#             Hallmark interferon sets, so adjusting for it can remove interferon signal itself.
#   MCP:      T-cell, B-lineage and fibroblast scores, built from lineage markers.
imm_sig <- tidyestimate::gene_sets$immune_signature
ifn_genes <- unique(c(hall$HALLMARK_INTERFERON_ALPHA_RESPONSE, hall$HALLMARK_INTERFERON_GAMMA_RESPONSE))
ov_ifn <- intersect(imm_sig, ifn_genes)
message(sprintf("ESTIMATE immune signature genes also in Hallmark interferon sets: %d of %d (%s)",
                length(ov_ifn), length(imm_sig), paste(head(ov_ifn, 15), collapse = ", ")))

ct  <- assay(se, "unstranded")[, tr$barcode]
gt  <- rd[rownames(ct), "gene_type"]
fit_race_h <- function(covs) {
  cdt <- data.frame(tr[, .(race, gleason_group = droplevels(gleason_group), age_s)], row.names = tr$barcode)
  for (cv in covs) cdt[[paste0(cv, "_z")]] <- zs(tr[[cv]])
  f <- as.formula(paste("~ gleason_group + age_s +", paste(paste0(covs, "_z"), collapse = " + "), "+ race"))
  d <- DESeqDataSetFromMatrix(ct, cdt, design = f)
  keep <- rowSums(counts(d) >= p[["min_count"]]) >= sum(tr$race == "Black")
  d <- DESeq(d[keep, ], parallel = threads > 1, BPPARAM = bp, quiet = TRUE)
  r <- results(d, name = "race_Black_vs_White", alpha = fdr)
  ok <- rd[rownames(d), "gene_type"] %in% c("protein_coding", "lncRNA")
  gsea_h(r$stat[ok], rd[rownames(d), "gene_name"][ok])
}
h_est <- fit_race_h(c("estimate_stromal", "estimate_immune"))
h_mcp <- fit_race_h(c("mcp_T_cells", "mcp_B_lineage", "mcp_Fibroblasts"))
mr <- fread(inp[["race_hallmark"]])[, .(pathway, NES_main = NES, padj_main = padj)]
rh <- Reduce(function(x, y) merge(x, y, by = "pathway"),
             list(mr, h_est[, .(pathway, NES_est = NES, padj_est = padj)],
                  h_mcp[, .(pathway, NES_mcp = NES, padj_mcp = padj)]))[order(padj_main)]
fwrite(rh, out[["race_adj_hallmark"]], sep = "\t")
shift_est <- median(rh$NES_est - rh$NES_main); shift_mcp <- median(rh$NES_mcp - rh$NES_main)
message(sprintf("Median NES shift across all 50 sets after adjustment: ESTIMATE %+.2f, MCP %+.2f (large uniform shifts suggest over-adjustment)",
                shift_est, shift_mcp))
message("Race Hallmark sets: main vs ESTIMATE-adjusted vs MCP-adjusted:")
print(rh[padj_main < fdr | padj_est < fdr | padj_mcp < fdr,
         .(pathway = short(pathway), NES_main = round(NES_main, 2), padj_main = signif(padj_main, 2),
           NES_est = round(NES_est, 2), padj_est = signif(padj_est, 2),
           NES_mcp = round(NES_mcp, 2), padj_mcp = signif(padj_mcp, 2))])
ifn <- rh[grepl("INTERFERON", pathway)]
gv  <- function(pw, col) { v <- ifn[pathway == pw][[col]]; if (length(v)) v else NA_real_ }
IFA <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"; IFG <- "HALLMARK_INTERFERON_GAMMA_RESPONSE"

# ---- summary -----------------------------------------------------
summ <- data.table(
  metric = c("paired_stromal_tumor_lower_pairs", "paired_immune_padj",
             "cor_tumor_stromal", "cor_tumor_immune", "cor_tumor_smooth_muscle", "cor_tumor_basal",
             "adj_paired_genes_padj05", "adj_paired_genes_padj05_lfc1", "full_adj_paired_genes_padj05",
             "hallmark_robust", "hallmark_composition_sensitive", "hallmark_emerging",
             "hallmark_robust_full", "hallmark_composition_sensitive_full", "hallmark_emerging_full",
             "race_immune_diff_sd", "race_immune_padj", "estimate_immune_ifn_overlap_genes",
             "race_nes_median_shift_estimate", "race_nes_median_shift_mcp",
             "race_ifna_nes_est", "race_ifna_padj_est", "race_ifng_nes_est", "race_ifng_padj_est",
             "race_ifna_nes_mcp", "race_ifna_padj_mcp", "race_ifng_nes_mcp", "race_ifng_padj_mcp"),
  value = c(pt[score == "estimate_stromal", pairs - tumor_higher_pairs], signif(pt[score == "estimate_immune", padj], 3),
            round(r_cs, 3), round(r_ci, 3), round(r_sm, 3), round(r_ba, 3),
            n_adj_any, n_adj_deg, n_full_any,
            sum(ah$class == "robust to composition"), sum(ah$class == "composition-sensitive"),
            sum(ah$class == "emerges after adjustment"),
            sum(ah$class_full == "robust to composition"), sum(ah$class_full == "composition-sensitive"),
            sum(ah$class_full == "emerges after adjustment"),
            round(rc[score == "estimate_immune", diff_sd_black_vs_white], 3), signif(rc[score == "estimate_immune", padj], 3),
            length(ov_ifn), round(shift_est, 3), round(shift_mcp, 3),
            round(gv(IFA, "NES_est"), 3), signif(gv(IFA, "padj_est"), 3),
            round(gv(IFG, "NES_est"), 3), signif(gv(IFG, "padj_est"), 3),
            round(gv(IFA, "NES_mcp"), 3), signif(gv(IFA, "padj_mcp"), 3),
            round(gv(IFG, "NES_mcp"), 3), signif(gv(IFG, "padj_mcp"), 3)))
fwrite(summ, out[["summary"]], sep = "\t")
message("Deconvolution summary:"); print(summ)
message("Done.")
