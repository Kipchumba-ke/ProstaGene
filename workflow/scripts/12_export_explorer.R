#!/usr/bin/env Rscript
# ---------------------------------------------------------------
# 12_export_explorer.R — export pipeline results for the web explorer
#
# Writes compact JSON files to explorer/public/data/ that the React
# app reads at run time. Only summary-level results are exported
# (gene statistics, pathway scores, model outputs); per-sample data
# are limited to de-identified TCGA expression for the top genes and
# the survival/score table, all from open-access TCGA data.
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
  library(jsonlite)
})

if (exists("snakemake")) {
  inp <- snakemake@input; out_dir <- dirname(snakemake@output[["summary"]]); p <- snakemake@params
  con <- file(snakemake@log[[1]], open = "wt"); sink(con); sink(con, type = "message")
} else {
  inp <- list(
    sheet = "data/processed/sample_sheet_paired.tsv", de = "results/de/deseq2_all_genes.tsv",
    sens = "results/de/sensitivity_summary.tsv", vsd = "data/processed/prad_vst.rds",
    hallmark = "results/gsea/fgsea_hallmark.tsv", validation = "results/validation/validation_summary.tsv",
    geo_de = "results/validation/geo_deseq2_all_genes.tsv", geo_hallmark = "results/validation/geo_fgsea_hallmark.tsv",
    ancestry = "results/ancestry/ancestry_summary.tsv", race_de = "results/ancestry/tcga_race_deseq2.tsv",
    race_hallmark = "results/ancestry/tcga_race_hallmark.tsv",
    deconv = "results/deconv/deconv_summary.tsv", adj_hallmark = "results/deconv/adjusted_paired_hallmark.tsv",
    paired_comp = "results/deconv/paired_composition_tests.tsv",
    surv_summary = "results/survival/survival_summary.tsv", surv_data = "results/survival/survival_dataset.tsv",
    cox = "results/survival/cox_results.tsv",
    ga = "results/genetic_ancestry/genetic_ancestry_summary.tsv", afr_de = "results/genetic_ancestry/afr_deseq2.tsv",
    afr_hallmark = "results/genetic_ancestry/afr_hallmark.tsv", germline = "results/genetic_ancestry/germline_gene_tests.tsv",
    sb = "results/subtypes/subtype_summary.tsv", sb_sheet = "results/subtypes/subtype_sample_sheet.tsv",
    sb_hallmark = "results/subtypes/subtype_hallmark.tsv", sb_erg = "results/subtypes/erg_vs_ancestry.tsv",
    sb_cox = "results/subtypes/subtype_cox.tsv", sb_paired = "results/subtypes/paired_by_erg_hallmark.tsv")
  out_dir <- "explorer/public/data"
  p <- list(n_expr_genes = 300)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
w <- function(x, name) {
  f <- file.path(out_dir, name)
  write_json(x, f, auto_unbox = TRUE, digits = NA, na = "null", null = "null")
  message(sprintf("  %-18s %7.1f KB", name, file.size(f) / 1024))
}
r3 <- function(x) round(x, 3)
s3 <- function(x) signif(x, 3)
metrics <- function(path) {
  d <- fread(path, colClasses = "character")
  v <- lapply(d$value, function(z) { n <- suppressWarnings(as.numeric(z)); if (!is.na(n)) n else z })
  setNames(v, d$metric)
}

# ---- summary ------------------------------------------------------
sheet <- fread(inp[["sheet"]])
de    <- fread(inp[["de"]])
sens  <- fread(inp[["sens"]])
hall  <- fread(inp[["hallmark"]])
summary <- list(
  generated = format(Sys.time(), "%Y-%m-%d"),
  design = list(pairs = uniqueN(sheet$patient), samples = nrow(sheet), genes_tested = nrow(de)),
  de = list(degs = sum(de$direction != "ns"), up = sum(de$direction == "up"), down = sum(de$direction == "down"),
            any_fc = sum(de$padj < 0.05, na.rm = TRUE),
            top = I(head(de[direction != "ns" & !is.na(gene_name)][order(padj)]$gene_name, 10))),
  sensitivity = as.list(sens[1]),
  pathways = list(hallmark_sig = sum(hall$padj < 0.05), up = sum(hall$padj < 0.05 & hall$NES > 0),
                  down = sum(hall$padj < 0.05 & hall$NES < 0)),
  validation = metrics(inp[["validation"]]),
  ancestry = metrics(inp[["ancestry"]]),
  composition = metrics(inp[["deconv"]]),
  survival = metrics(inp[["surv_summary"]]),
  genetic_ancestry = metrics(inp[["ga"]]),
  subtypes = metrics(inp[["sb"]]))
message("Writing explorer data to ", out_dir)
w(summary, "summary.json")

# ---- genes (columnar, one row per symbol) -------------------------
g <- de[!is.na(gene_name) & gene_name != ""][order(padj)][!duplicated(gene_name)]
geo  <- fread(inp[["geo_de"]])[!is.na(symbol)][order(padj)][!duplicated(symbol)]
race <- fread(inp[["race_de"]])[!is.na(gene_name)][order(padj_all)][!duplicated(gene_name)]
afr  <- fread(inp[["afr_de"]])[!is.na(gene_name)][order(pvalue)][!duplicated(gene_name)]
germ <- fread(inp[["germline"]])
m <- function(tab, key, col, keys) tab[[col]][match(keys, tab[[key]])]
genes <- list(
  symbol = g$gene_name, type = g$gene_type, baseMean = r3(log10(g$baseMean + 1)),
  lfc = r3(g$log2FC_shrunk), padj = s3(g$padj), dir = g$direction,
  geo_lfc = r3(m(geo, "symbol", "log2FC_shrunk", g$gene_name)), geo_padj = s3(m(geo, "symbol", "padj", g$gene_name)),
  race_lfc = r3(m(race, "gene_name", "log2FC_shrunk", g$gene_name)), race_padj = s3(m(race, "gene_name", "padj", g$gene_name)),
  afr_lfc = r3(m(afr, "gene_name", "log2FC_per_full_AFR", g$gene_name)), afr_padj = s3(m(afr, "gene_name", "padj", g$gene_name)),
  germline = I(germ$gene[germ$likely_germline %in% TRUE]))
w(genes, "genes.json")

# ---- pathways ----------------------------------------------------
pw <- hall[, .(pathway, NES = r3(NES), padj = s3(padj), size)]
add <- function(base, path, cols, new) {
  d <- fread(path)[, c("pathway", cols), with = FALSE]; setnames(d, cols, new); merge(base, d, by = "pathway", all.x = TRUE)
}
pw <- add(pw, inp[["adj_hallmark"]], c("NES_full", "padj_full", "class_full"), c("NES_adj", "padj_adj", "composition"))
pw <- add(pw, inp[["geo_hallmark"]], c("NES", "padj"), c("NES_geo", "padj_geo"))
pw <- add(pw, inp[["race_hallmark"]], c("NES", "padj"), c("NES_race", "padj_race"))
pw <- add(pw, inp[["afr_hallmark"]], c("NES_afr", "padj_afr"), c("NES_afr", "padj_afr"))
pw <- add(pw, inp[["sb_paired"]], c("NES_ERGpos", "NES_ERGneg", "p_wilcox"), c("NES_ergpos", "NES_ergneg", "p_erg_diff"))
num_cols <- setdiff(names(pw), c("pathway", "composition", "size"))
for (cc in num_cols) pw[, (cc) := if (grepl("^NES", cc)) r3(get(cc)) else s3(get(cc))]
pw[, label := gsub("_", " ", sub("^HALLMARK_", "", pathway))]
w(pw[order(padj)], "hallmark.json")

# ---- expression for top genes (paired samples) --------------------
vsd <- readRDS(inp[["vsd"]])
gsym <- if ("gene_name" %in% names(mcols(vsd))) mcols(vsd)$gene_name else rownames(vsd)
want <- unique(c(head(g[direction != "ns"]$gene_name, p[["n_expr_genes"]]),
                 c("AMACR", "PCA3", "HPN", "DLX1", "TDRD1", "ERG", "KLK3", "AR", "KRT5", "TP63", "MYC")))
idx <- match(want, gsym); want <- want[!is.na(idx)]; idx <- idx[!is.na(idx)]
mat <- round(assay(vsd)[idx, , drop = FALSE], 2)
expr <- list(samples = data.table(patient = as.character(vsd$patient), condition = as.character(vsd$condition)),
             genes = setNames(lapply(seq_along(want), function(i) unname(mat[i, ])), want))
w(expr, "expression.json")

# ---- survival ----------------------------------------------------
sd <- fread(inp[["surv_data"]])
score_cols <- intersect(c("tcga_tumor_signature", "hm_E2F_TARGETS", "hm_G2M_CHECKPOINT", "hm_ANDROGEN_RESPONSE",
                          "hm_INTERFERON_GAMMA_RESPONSE", "hm_MYOGENESIS", "estimate_stromal", "estimate_immune"), names(sd))
surv <- list(
  rows = sd[, c("pfi_months", "pfi", "gleason_group", score_cols), with = FALSE][, pfi_months := round(pfi_months, 1)],
  scores = I(score_cols),
  cox = fread(inp[["cox"]])[, .(score, model, HR = r3(HR_per_SD), lower = r3(lower95), upper = r3(upper95),
                                padj = s3(padj), c_index = r3(c_index))])
for (cc in score_cols) surv$rows[, (cc) := r3(get(cc))]
w(surv, "survival.json")

# ---- subtypes -----------------------------------------------------
ss <- fread(inp[["sb_sheet"]])
sub <- list(
  counts = ss[!is.na(subtype_group), .N, by = subtype_group][order(-N)],
  by_ancestry = ss[!is.na(subtype_group) & !is.na(afr)][, afr_group := ifelse(afr >= 0.5, "AFR >= 0.5", "AFR < 0.5")][
    , .N, by = .(afr_group, subtype_group)],
  erg_ancestry = fread(inp[["sb_erg"]]),
  hallmark = fread(inp[["sb_hallmark"]])[padj < 0.05][, .(contrast, pathway, NES = r3(NES), padj = s3(padj))],
  cox = fread(inp[["sb_cox"]])[, .(model, term, HR = r3(HR), lower = r3(lower95), upper = r3(upper95), p = s3(p), n, events)])
w(sub, "subtypes.json")

comp <- fread(inp[["paired_comp"]])[, .(score, median_normal = r3(median_normal), median_tumor = r3(median_tumor),
                                         tumor_higher_pairs, pairs, padj = s3(padj))]
w(comp, "composition.json")
message("Done.")
