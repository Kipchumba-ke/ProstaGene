#!/usr/bin/env Rscript
# ---------------------------------------------------------------
# 07_ancestry.R — exploratory ancestry analyses (M6), revised
#
# Part A (GSE104131, 7 AA + 6 EA pairs)
#   - Gene level: nested DESeq2 design (pairs within ancestry),
#     ~ ancestry + ancestry:pair_n + ancestry:condition
#   - Pathway level: PATIENT is the unit of replication.
#     Per pair, tumor-minus-normal VST change; per Hallmark set,
#     mean change of its genes; AA (7) vs EA (6) by exact Wilcoxon.
#     (Pre-ranked GSEA on the interaction statistic was dropped:
#     it treats genes as independent and is badly anti-conservative
#     with 13 patients.)
#   - Pair QC: each pair's tumor effect vs TCGA (Spearman rho)
#
# Part B (TCGA-PRAD tumors, Black vs White, self-reported race)
#   - ~ gleason_group + age + race
#   - DEG reporting restricted to protein-coding + lncRNA genes, with
#     BH recomputed in that set. Mitochondrial pseudogene hits are
#     written separately: likely a mapping artifact of mtDNA
#     haplogroup divergence from the (European) reference.
#   - Site-matched sensitivity: sites with >= 3 Black and >= 3 White
#     patients, adding tissue source site to the model.
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(BiocParallel)
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
  inp <- list(geo_dds = "data/processed/geo_dds_fitted.rds",
              se_full = "data/processed/prad_se_full.rds",
              tcga_de = "results/de/deseq2_all_genes.tsv")
  a <- "results/ancestry"
  out <- list(geo_table = file.path(a, "geo_ancestry_tumor_effects.tsv"),
              geo_scatter = file.path(a, "geo_ancestry_lfc_scatter.png"),
              geo_hallmark = file.path(a, "geo_ancestry_hallmark_pairlevel.tsv"),
              geo_hallmark_plot = file.path(a, "geo_ancestry_hallmark_pairlevel.png"),
              geo_pair_qc = file.path(a, "geo_pair_qc.tsv"),
              tcga_sheet = file.path(a, "tcga_race_sample_sheet.tsv"),
              tcga_table = file.path(a, "tcga_race_deseq2.tsv"),
              tcga_mt = file.path(a, "tcga_race_mt_pseudogene_hits.tsv"),
              tcga_volcano = file.path(a, "tcga_race_volcano.png"),
              tcga_hallmark = file.path(a, "tcga_race_hallmark.tsv"),
              tcga_hallmark_plot = file.path(a, "tcga_race_hallmark.png"),
              tcga_site_hallmark = file.path(a, "tcga_race_site_matched_hallmark.tsv"),
              tcga_site_plot = file.path(a, "tcga_race_site_matched_nes.png"),
              summary = file.path(a, "ancestry_summary.tsv"))
  p <- list(fdr = 0.05, lfc = 1, min_count = 10, min_per_site = 3, seed = 42)
  threads <- 4
}
for (f in unlist(out)) dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 11))
fdr <- p[["fdr"]]; lfc_cut <- p[["lfc"]]
set.seed(p[["seed"]])
bp <- if (threads > 1) MulticoreParam(threads) else SerialParam()
anc_cols  <- c(AA = "#B5651D", EA = "#3A6EA5")
race_cols <- c(`higher in Black` = "#B5651D", `higher in White` = "#3A6EA5", ns = "grey75")

# ---- helpers -----------------------------------------------------
new_api <- "collection" %in% names(formals(msigdbr::msigdbr))
hdf  <- if (new_api) msigdbr(species = "Homo sapiens", collection = "H") else
                     msigdbr(species = "Homo sapiens", category = "H")
hall <- split(hdf$gene_symbol, hdf$gs_name)
short <- function(x) gsub("_", " ", sub("^HALLMARK_", "", x))

gsea_h <- function(stat, symbol) {
  d <- data.table(stat, symbol)[!is.na(stat) & !is.na(symbol) & symbol != ""]
  d <- d[order(-abs(stat))][!duplicated(symbol)]
  r <- sort(setNames(d$stat, d$symbol), decreasing = TRUE)
  fgsea(hall, r, minSize = 15, maxSize = 500, eps = 0, nproc = threads)[, .(pathway, NES, pval, padj, size)]
}
tidy_res <- function(r, prefix) {
  d <- data.table(log2FC = r$log2FoldChange, stat = r$stat, pvalue = r$pvalue, padj = r$padj)
  setnames(d, paste0(prefix, "_", names(d))); d
}

# =================================================================
# Part A — GSE104131
# =================================================================
dg  <- readRDS(inp[["geo_dds"]])
cdg <- as.data.frame(colData(dg))[, c("pair", "condition", "ancestry")]
cdg$pair      <- as.character(cdg$pair)
cdg$condition <- relevel(factor(cdg$condition), ref = "normal")
cdg$ancestry  <- relevel(factor(cdg$ancestry), ref = "EA")
cdg$pair_n    <- factor(ave(cdg$pair, cdg$ancestry, FUN = function(x) as.integer(factor(x))))

# ---- A1. gene level: nested design -------------------------------
mm <- model.matrix(~ ancestry + ancestry:pair_n + ancestry:condition, cdg)
mm <- mm[, colSums(mm != 0) > 0]                 # drop empty columns (7 vs 6 pairs)
colnames(mm) <- make.names(colnames(mm))
da <- DESeqDataSetFromMatrix(counts(dg), cdg, design = ~ 1)
da <- DESeq(da, full = mm, parallel = threads > 1, BPPARAM = bp, quiet = TRUE)
rn <- resultsNames(da)
aa <- grep("^ancestryAA\\.conditiontumor$", rn, value = TRUE)
ea <- grep("^ancestryEA\\.conditiontumor$", rn, value = TRUE)
stopifnot(length(aa) == 1, length(ea) == 1)
r_aa  <- results(da, name = aa, alpha = fdr)
r_ea  <- results(da, name = ea, alpha = fdr)
r_int <- results(da, contrast = list(aa, ea), alpha = fdr)

gt <- cbind(data.table(entrez = rownames(da), symbol = mcols(dg)$symbol, baseMean = r_aa$baseMean),
            tidy_res(r_aa, "AA"), tidy_res(r_ea, "EA"), tidy_res(r_int, "diff"))[order(diff_padj)]
fwrite(gt, out[["geo_table"]], sep = "\t")

n_deg <- function(lfc, padj) sum(!is.na(padj) & padj < fdr & abs(lfc) >= lfc_cut)
n_aa  <- n_deg(gt$AA_log2FC, gt$AA_padj)
n_ea  <- n_deg(gt$EA_log2FC, gt$EA_padj)
n_int <- sum(!is.na(gt$diff_padj) & gt$diff_padj < fdr)
rho_ae <- cor(gt$AA_log2FC, gt$EA_log2FC, method = "spearman", use = "complete.obs")
message(sprintf("GEO tumor-vs-normal DEGs: AA %d, EA %d | genes whose effect differs by ancestry (padj < %.2f): %d | AA vs EA log2FC rho = %.2f",
                n_aa, n_ea, fdr, n_int, rho_ae))
if (n_int) print(gt[!is.na(diff_padj) & diff_padj < fdr,
                    .(symbol, AA_log2FC = round(AA_log2FC, 2), EA_log2FC = round(EA_log2FC, 2),
                      diff_padj = signif(diff_padj, 3))])

gt[, class := fifelse(!is.na(diff_padj) & diff_padj < fdr, "differs by ancestry",
              fifelse((!is.na(AA_padj) & AA_padj < fdr) | (!is.na(EA_padj) & EA_padj < fdr),
                      "tumor effect (either group)", "ns"))]
g <- ggplot(gt, aes(EA_log2FC, AA_log2FC, colour = class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey50") +
  geom_point(size = 0.6, alpha = 0.6) +
  geom_text_repel(data = gt[class == "differs by ancestry"], aes(label = symbol),
                  size = 3, colour = "black") +
  scale_colour_manual(values = c(ns = "grey85", `tumor effect (either group)` = "#9DB4C0",
                                 `differs by ancestry` = "#C8553D"), name = NULL) +
  coord_cartesian(xlim = c(-8, 8), ylim = c(-8, 8)) +
  labs(x = "Tumor vs normal log2FC, European-American (6 pairs)",
       y = "Tumor vs normal log2FC, African-American (7 pairs)",
       title = "GSE104131: is the tumor effect ancestry-specific?",
       subtitle = sprintf("Spearman rho = %.2f; %d genes differ at padj < %.2f", rho_ae, n_int, fdr)) +
  theme(legend.position = "bottom")
ggsave(out[["geo_scatter"]], g, width = 7, height = 7.3, dpi = 150)

# ---- A2. pathway level, patient as unit --------------------------
vs  <- vst(dg, blind = TRUE)
vm  <- assay(vs)
sym <- mcols(dg)$symbol
ok  <- !is.na(sym) & sym != "" & !duplicated(sym)
vm  <- vm[ok, ]; rownames(vm) <- sym[ok]

pairs <- sort(unique(cdg$pair))
pair_anc <- as.character(cdg$ancestry[match(pairs, cdg$pair)])
delta <- sapply(pairs, function(pp) {
  tu <- rownames(cdg)[cdg$pair == pp & cdg$condition == "tumor"]
  no <- rownames(cdg)[cdg$pair == pp & cdg$condition == "normal"]
  vm[, tu] - vm[, no]
})
colnames(delta) <- pairs

# Pair QC: how TCGA-like is each pair's tumor effect?
tcga <- fread(inp[["tcga_de"]])[!is.na(gene_name) & gene_name != ""][order(padj)][!duplicated(gene_name)]
shared <- intersect(rownames(delta), tcga$gene_name)
tl <- setNames(tcga$log2FC_shrunk, tcga$gene_name)[shared]
pq <- data.table(pair = pairs, ancestry = pair_anc,
                 rho_vs_tcga = sapply(pairs, function(pp) cor(delta[shared, pp], tl, method = "spearman")))
pq[, low_concordance := rho_vs_tcga < median(rho_vs_tcga) - 2 * mad(rho_vs_tcga)]
pq <- pq[order(rho_vs_tcga)]
fwrite(pq, out[["geo_pair_qc"]], sep = "\t")
message("Per-pair concordance with TCGA tumor effect (Spearman rho):"); print(pq)

# Hallmark score per pair = mean tumor-minus-normal VST change of set genes
sets_ok <- Filter(function(s) length(intersect(s, rownames(delta))) >= 15, hall)
sc <- t(sapply(sets_ok, function(s) colMeans(delta[intersect(s, rownames(delta)), , drop = FALSE])))
is_aa <- pair_anc == "AA"
ph <- data.table(
  pathway  = rownames(sc),
  n_genes  = sapply(sets_ok, function(s) length(intersect(s, rownames(delta)))),
  median_AA = apply(sc[, is_aa, drop = FALSE], 1, median),
  median_EA = apply(sc[, !is_aa, drop = FALSE], 1, median),
  p_wilcox = apply(sc, 1, function(x) wilcox.test(x[is_aa], x[!is_aa], exact = TRUE)$p.value))
ph[, diff_AA_minus_EA := median_AA - median_EA]
ph[, padj := p.adjust(p_wilcox, "BH")]

# Descriptive within-group GSEA (standard tumor-vs-normal ranking, one per group)
h_aa <- gsea_h(gt$AA_stat, gt$symbol)[, .(pathway, AA_NES = NES, AA_gsea_padj = padj)]
h_ea <- gsea_h(gt$EA_stat, gt$symbol)[, .(pathway, EA_NES = NES, EA_gsea_padj = padj)]
ph <- merge(merge(ph, h_aa, by = "pathway", all.x = TRUE), h_ea, by = "pathway", all.x = TRUE)[order(p_wilcox)]
fwrite(ph, out[["geo_hallmark"]], sep = "\t")
min_p <- 2 / choose(length(pairs), sum(is_aa))       # smallest attainable two-sided exact p
message(sprintf("Pair-level Hallmark test (AA %d vs EA %d pairs; smallest possible p = %.4f): %d sets padj < %.2f, %d raw p < 0.05",
                sum(is_aa), sum(!is_aa), min_p, sum(ph$padj < fdr), fdr, sum(ph$p_wilcox < 0.05)))
print(head(ph[, .(pathway = short(pathway), median_AA = round(median_AA, 3), median_EA = round(median_EA, 3),
                  p_wilcox = signif(p_wilcox, 3), padj = signif(padj, 3))], 12))

top <- head(ph$pathway, 8)
sl <- data.table(pathway = rep(rownames(sc), ncol(sc)), pair = rep(colnames(sc), each = nrow(sc)),
                 score = as.vector(sc))[pathway %in% top]
sl[, ancestry := pair_anc[match(pair, pairs)]]
sl[, low_qc := pair %in% pq[low_concordance == TRUE]$pair]
sl[, facet := factor(sprintf("%s\np = %.3f", short(pathway), ph$p_wilcox[match(pathway, ph$pathway)]),
                     levels = sprintf("%s\np = %.3f", short(top), ph$p_wilcox[match(top, ph$pathway)]))]
g <- ggplot(sl, aes(ancestry, score, colour = ancestry)) +
  geom_hline(yintercept = 0, colour = "grey70") +
  geom_boxplot(outlier.shape = NA, width = 0.5, colour = "grey50") +
  geom_point(aes(shape = low_qc), size = 2.6, position = position_jitter(width = 0.08, seed = 1)) +
  scale_colour_manual(values = anc_cols, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4), labels = c("ok", "low TCGA concordance"), name = "Pair") +
  facet_wrap(~ facet, nrow = 2, scales = "free_y") +
  labs(x = NULL, y = "Mean tumor-minus-normal change (VST)",
       title = "Hallmark tumor effects per patient pair, by ancestry (GSE104131)",
       subtitle = "Each point is one patient; exact Wilcoxon AA vs EA (raw p shown)") +
  theme(legend.position = "bottom", strip.text = element_text(size = 8))
ggsave(out[["geo_hallmark_plot"]], g, width = 11, height = 6.5, dpi = 150)

# =================================================================
# Part B — TCGA-PRAD tumors, Black vs White
# =================================================================
se <- readRDS(inp[["se_full"]])
bc <- colnames(se)
tm <- data.table(barcode = bc, patient = substr(bc, 1, 12), code = substr(bc, 14, 15),
                 vial = substr(bc, 16, 16), analyte = substr(bc, 20, 20), plate = substr(bc, 22, 25),
                 tss = substr(bc, 6, 7))[code == "01"]
tm[, arank := match(analyte, c("H", "R", "T"), nomatch = 4L)]
tm <- tm[order(patient, vial, arank, -rank(plate))][!duplicated(patient)]

cd <- as.data.frame(colData(se))[tm$barcode, ]
race_raw <- tolower(as.character(cd$race))
tm[, race := fifelse(race_raw == "black or african american", "Black",
             fifelse(race_raw == "white", "White", NA_character_))]
num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(x))))
gl <- if (all(c("primary_gleason_grade", "secondary_gleason_grade") %in% names(cd)))
  num(cd$primary_gleason_grade) + num(cd$secondary_gleason_grade) else rep(NA_real_, nrow(cd))
tm[, gleason := gl]
tm[, gleason_group := factor(fifelse(gleason <= 6, "GS6", fifelse(gleason == 7, "GS7", "GS8_10")),
                             levels = c("GS6", "GS7", "GS8_10"))]
tm[, age := if ("age_at_index" %in% names(cd)) num(cd$age_at_index) else num(cd$age_at_diagnosis) / 365.25]

tb <- tm[!is.na(race) & !is.na(gleason_group) & !is.na(age)]
tb[, race := factor(race, levels = c("White", "Black"))]
tb[, age_s := as.numeric(scale(age))]
message("TCGA tumors analysed: ", nrow(tb), " [Black ", sum(tb$race == "Black"), ", White ", sum(tb$race == "White"), "]")
print(table(race = tb$race, gleason = tb$gleason_group))
fwrite(tb, out[["tcga_sheet"]], sep = "\t")

ri <- as.data.frame(rowData(se))
fit_race <- function(sub, design) {
  ct <- assay(se, "unstranded")[, sub$barcode]
  cdt <- data.frame(sub[, .(race, gleason_group = droplevels(gleason_group), age_s, tss = factor(tss))],
                    row.names = sub$barcode)
  d <- DESeqDataSetFromMatrix(ct, cdt, design = design)
  d <- d[rowSums(counts(d) >= p[["min_count"]]) >= sum(sub$race == "Black"), ]
  d <- DESeq(d, parallel = threads > 1, BPPARAM = bp, quiet = TRUE)
  r <- results(d, name = "race_Black_vs_White", alpha = fdr)
  list(dds = d, res = r)
}
annotate <- function(d, r, shrunk = NULL) {
  x <- data.table(gene_id = rownames(d), gene_name = ri[rownames(d), "gene_name"],
                  gene_type = ri[rownames(d), "gene_type"], baseMean = r$baseMean,
                  log2FC = r$log2FoldChange, stat = r$stat, pvalue = r$pvalue, padj_all = r$padj)
  if (!is.null(shrunk)) x[, log2FC_shrunk := shrunk$log2FoldChange]
  x[, reportable := gene_type %in% c("protein_coding", "lncRNA")]
  x[, mt_pseudogene := grepl("pseudogene", gene_type) & grepl("^MT", gene_name)]
  x[reportable == TRUE & !is.na(pvalue), padj := p.adjust(pvalue, "BH")]
  x[order(padj_all)]
}

# ---- B1. main model ----------------------------------------------
main <- fit_race(tb, ~ gleason_group + age_s + race)
sh <- lfcShrink(main$dds, coef = "race_Black_vs_White", type = "apeglm",
                parallel = threads > 1, BPPARAM = bp, quiet = TRUE)
tt <- annotate(main$dds, main$res, sh)
tt[, direction := fifelse(reportable & !is.na(padj) & padj < fdr & abs(log2FC_shrunk) >= lfc_cut,
                          fifelse(log2FC_shrunk > 0, "higher in Black", "higher in White"), "ns")]
fwrite(tt, out[["tcga_table"]], sep = "\t")

mt <- tt[mt_pseudogene == TRUE & !is.na(padj_all) & padj_all < fdr]
fwrite(mt, out[["tcga_mt"]], sep = "\t")
message(sprintf("Mitochondrial pseudogene hits (padj_all < %.2f): %d, of which %d higher in Black; %d of the top 25 genes overall",
                fdr, nrow(mt), sum(mt$log2FC > 0), sum(head(tt, 25)$mt_pseudogene)))

rep_tt <- tt[reportable == TRUE]
n_any <- sum(rep_tt$padj < fdr, na.rm = TRUE); n_deg_r <- sum(rep_tt$direction != "ns")
message(sprintf("Protein-coding + lncRNA genes tested: %d | padj < %.2f: %d | also |shrunk log2FC| >= %g: %d [Black %d, White %d]",
                nrow(rep_tt), fdr, n_any, lfc_cut, n_deg_r,
                sum(rep_tt$direction == "higher in Black"), sum(rep_tt$direction == "higher in White")))
print(head(rep_tt[order(padj), .(gene_name, gene_type, log2FC_shrunk = round(log2FC_shrunk, 2),
                                 padj = signif(padj, 3))], 15))

v <- rep_tt[!is.na(padj)]
v[, nl := -log10(pmax(padj, min(padj[padj > 0])))]
g <- ggplot(v, aes(log2FC_shrunk, nl, colour = direction)) +
  geom_point(size = 0.7, alpha = 0.7) +
  geom_hline(yintercept = -log10(fdr), linetype = "dashed", colour = "grey50") +
  geom_text_repel(data = head(v[padj < fdr][order(padj)], 20), aes(label = gene_name),
                  size = 2.8, colour = "black", max.overlaps = 30) +
  scale_colour_manual(values = race_cols, name = NULL) +
  labs(x = "log2FC, Black vs White tumors (apeglm-shrunk)", y = expression(-log[10]~adjusted~p),
       title = "TCGA-PRAD tumors by self-reported race",
       subtitle = sprintf("%d Black vs %d White; Gleason group + age adjusted; protein-coding + lncRNA only",
                          sum(tb$race == "Black"), sum(tb$race == "White"))) +
  theme(legend.position = "bottom")
ggsave(out[["tcga_volcano"]], g, width = 7.5, height = 6.5, dpi = 150)

ht <- gsea_h(rep_tt$stat, rep_tt$gene_name)[order(padj)]
fwrite(ht, out[["tcga_hallmark"]], sep = "\t")
message("TCGA Hallmark, Black vs White (padj < 0.05):")
print(ht[padj < fdr][order(-NES), .(pathway = short(pathway), NES = round(NES, 2), padj = signif(padj, 3))])

hp <- ht[padj < 0.25]
if (nrow(hp)) {
  hp[, label := factor(short(pathway), levels = short(pathway[order(NES)]))]
  g <- ggplot(hp, aes(NES, label, fill = NES > 0, alpha = padj < fdr)) +
    geom_col() +
    scale_fill_manual(values = c(`TRUE` = "#B5651D", `FALSE` = "#3A6EA5"),
                      labels = c(`TRUE` = "higher in Black", `FALSE` = "higher in White"), name = NULL) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4),
                       labels = c(`TRUE` = "padj < 0.05", `FALSE` = "padj < 0.25"), name = NULL) +
    labs(x = "NES", y = NULL, title = "Hallmark GSEA: TCGA-PRAD tumors, Black vs White",
         subtitle = "Wald-statistic ranking; Gleason group + age adjusted")
  ggsave(out[["tcga_hallmark_plot"]], g, width = 8.5, height = max(4, 0.3 * nrow(hp) + 1.5), dpi = 150)
} else file.create(out[["tcga_hallmark_plot"]])

# ---- B2. site-matched sensitivity --------------------------------
st <- tb[, .(black = sum(race == "Black"), white = sum(race == "White")), by = tss]
sites <- st[black >= p[["min_per_site"]] & white >= p[["min_per_site"]], tss]
ts <- tb[tss %in% sites]
message(sprintf("Site-matched subset: sites %s | %d Black, %d White",
                paste(sites, collapse = ", "), sum(ts$race == "Black"), sum(ts$race == "White")))
print(st[tss %in% sites])
site <- fit_race(ts, ~ tss + gleason_group + age_s + race)
ts_tab <- annotate(site$dds, site$res)[reportable == TRUE]
hs <- gsea_h(ts_tab$stat, ts_tab$gene_name)
cmp <- merge(ht[, .(pathway, NES_main = NES, padj_main = padj)],
             hs[, .(pathway, NES_site = NES, padj_site = padj)], by = "pathway")[order(padj_main)]
fwrite(cmp, out[["tcga_site_hallmark"]], sep = "\t")
r_site <- cor(cmp$NES_main, cmp$NES_site)
sig_main <- cmp[padj_main < fdr]
same_dir <- mean(sign(sig_main$NES_main) == sign(sig_main$NES_site))
message(sprintf("Site-matched vs main Hallmark NES: r = %.2f; %.0f%% of main-significant sets keep direction; %d significant in site-matched model",
                r_site, 100 * same_dir, sum(cmp$padj_site < fdr)))
print(cmp[grepl("INTERFERON|ALLOGRAFT|INFLAMMATORY|TNFA", pathway),
          .(pathway = short(pathway), NES_main = round(NES_main, 2), padj_main = signif(padj_main, 3),
            NES_site = round(NES_site, 2), padj_site = signif(padj_site, 3))])

cmp[, status := fifelse(padj_main < fdr & padj_site < fdr & sign(NES_main) == sign(NES_site), "significant in both",
                fifelse(padj_main < fdr, "main only", fifelse(padj_site < fdr, "site-matched only", "neither")))]
g <- ggplot(cmp, aes(NES_main, NES_site, colour = status)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 2.4) +
  geom_text_repel(data = cmp[status == "significant in both"], aes(label = short(pathway)),
                  size = 2.6, colour = "black", max.overlaps = 40) +
  scale_colour_manual(values = c(`significant in both` = "#6A3D9A", `main only` = "#E3A587",
                                 `site-matched only` = "#9DB4C0", neither = "grey75"), name = NULL) +
  labs(x = "NES, all sites (Gleason + age)", y = "NES, site-matched (site + Gleason + age)",
       title = "Black vs White Hallmark NES: robustness to tissue source site",
       subtitle = sprintf("Sites %s; r = %.2f", paste(sites, collapse = ", "), r_site)) +
  theme(legend.position = "bottom")
ggsave(out[["tcga_site_plot"]], g, width = 7.5, height = 7, dpi = 150)

# ---- summary -----------------------------------------------------
summ <- data.table(
  metric = c("geo_pairs_AA", "geo_pairs_EA", "geo_DEGs_AA", "geo_DEGs_EA", "geo_interaction_genes",
             "geo_AA_vs_EA_lfc_spearman", "geo_pairlevel_hallmark_padj05", "geo_pairlevel_hallmark_rawp05",
             "geo_low_concordance_pairs", "tcga_tumors_black", "tcga_tumors_white",
             "tcga_mt_pseudogene_hits", "tcga_race_genes_padj05", "tcga_race_genes_padj05_lfc1",
             "tcga_race_hallmark_padj05", "site_matched_black", "site_matched_white",
             "site_matched_nes_r", "site_matched_same_direction", "site_matched_hallmark_padj05"),
  value = c(sum(is_aa), sum(!is_aa), n_aa, n_ea, n_int, round(rho_ae, 3),
            sum(ph$padj < fdr), sum(ph$p_wilcox < 0.05), sum(pq$low_concordance),
            sum(tb$race == "Black"), sum(tb$race == "White"), nrow(mt), n_any, n_deg_r,
            sum(ht$padj < fdr), sum(ts$race == "Black"), sum(ts$race == "White"),
            round(r_site, 3), round(same_dir, 3), sum(cmp$padj_site < fdr)))
fwrite(summ, out[["summary"]], sep = "\t")
message("Ancestry summary:"); print(summ)
message("Done.")