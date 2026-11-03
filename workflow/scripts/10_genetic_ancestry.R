#!/usr/bin/env Rscript
# ---------------------------------------------------------------
# 10_genetic_ancestry.R — genetic ancestry instead of self-reported race (M9)
#
# Source: TCGA Ancestry & Molecular Correlates working group
# (Carrot-Zhang et al. 2020, Cancer Cell), GDC publication page:
#   Admixture_by_sample.txt  (per-sample ancestry proportions)
#   UCSF_Ancestry_Calls.csv  (UCSF PCA-based calls; fallback)
#
#  1. Genetic African-ancestry fraction (AFR) per patient; agreement
#     with self-reported race
#  2. Tumors from two-way AFR/EUR patients: DESeq2
#       ~ gleason_group + age + AFR   (AFR on a 0-1 scale)
#     Hallmark GSEA; compare with the self-reported-race result (M6)
#  3. Germline test for the top race-associated genes: does expression
#     track AFR fraction (a) in all tumors, (b) WITHIN self-identified
#     Black patients with AFR >= 0.5 (admixture dose-response), and (c) in matched
#     NORMAL tissue? Tracking AFR in normals and within Black patients
#     points to inherited variation rather than tumor biology.
#  4. Interferon-gamma score vs AFR, with and without lymphocyte scores
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

if (exists("snakemake")) {
  inp <- snakemake@input; out <- snakemake@output; p <- snakemake@params
  threads <- snakemake@threads
  con <- file(snakemake@log[[1]], open = "wt"); sink(con); sink(con, type = "message")
} else {
  inp <- list(admix = "data/raw/ancestry/Admixture_by_sample.txt",
              ucsf = "data/raw/ancestry/UCSF_Ancestry_Calls.csv",
              se_full = "data/processed/prad_se_full.rds",
              se_paired = "data/processed/prad_se_paired.rds",
              race_table = "results/ancestry/tcga_race_deseq2.tsv",
              race_hallmark = "results/ancestry/tcga_race_hallmark.tsv",
              scores = "results/deconv/cell_scores.tsv")
  a <- "results/genetic_ancestry"
  out <- list(sheet = file.path(a, "ancestry_sample_sheet.tsv"),
              race_plot = file.path(a, "race_vs_genetic_ancestry.png"),
              afr_table = file.path(a, "afr_deseq2.tsv"),
              afr_hallmark = file.path(a, "afr_hallmark.tsv"),
              nes_plot = file.path(a, "afr_vs_race_nes.png"),
              germline = file.path(a, "germline_gene_tests.tsv"),
              germline_plot = file.path(a, "germline_gene_scatter.png"),
              ifn = file.path(a, "ifn_vs_afr.tsv"),
              summary = file.path(a, "genetic_ancestry_summary.tsv"))
  p <- list(fdr = 0.05, min_count = 10, min_two_way = 0.9, n_genes = 30, seed = 42)
  threads <- 4
}
for (f in unlist(out)) dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 11))
fdr <- p[["fdr"]]
set.seed(p[["seed"]])
bp <- if (threads > 1) MulticoreParam(threads) else SerialParam()
short <- function(x) gsub("_", " ", sub("^HALLMARK_", "", x))
num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.eE-]", "", as.character(x))))   # ancestry proportions (may use 1e-05)
num_clin <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(x)))) # clinical text, e.g. "Pattern 3"
race_cols <- c(White = "#3A6EA5", Black = "#B5651D", `Other/unknown` = "grey55")

# =================================================================
# 1. Read ancestry proportions (column names detected, then printed)
# =================================================================
read_ancestry <- function(path) {
  if (!file.exists(path)) return(NULL)
  x <- fread(path, sep = "auto", header = TRUE, fill = TRUE)
  message("\n", basename(path), ": ", nrow(x), " rows; columns: ", paste(names(x), collapse = ", "))
  print(head(x, 3))
  idcol <- names(x)[sapply(x, function(v) mean(grepl("^TCGA-", as.character(v))) > 0.5)][1]
  find  <- function(pat) {
    cand <- grep(pat, names(x), ignore.case = TRUE, value = TRUE)
    if (!length(cand)) return(NA_character_)
    cand <- cand[vapply(cand, function(cn) { v <- num(x[[cn]]); mean(!is.na(v)) > 0.5 && max(v, na.rm = TRUE) <= 100 }, logical(1))]
    if (length(cand)) cand[1] else NA_character_
  }
  afr <- find("afr"); eur <- find("eur"); eas <- find("eas")
  if (is.na(idcol) || is.na(afr) || is.na(eur)) {
    message("Could not identify ID / AFR / EUR columns in ", basename(path)); return(NULL)
  }
  d <- data.table(id = as.character(x[[idcol]]), afr = num(x[[afr]]), eur = num(x[[eur]]),
                  eas = if (!is.na(eas)) num(x[[eas]]) else NA_real_)
  if (max(d$afr, na.rm = TRUE) > 1.5) d[, `:=`(afr = afr / 100, eur = eur / 100, eas = eas / 100)]
  callcol <- grep("consensus|call|ancestry|ethnic", names(x), ignore.case = TRUE, value = TRUE)
  callcol <- callcol[!callcol %in% c(idcol, afr, eur, eas)][1]
  d[, call := if (!is.na(callcol)) as.character(x[[callcol]]) else NA_character_]
  d[, patient := substr(id, 1, 12)]
  message(sprintf("Using ID = '%s', AFR = '%s', EUR = '%s', EAS = '%s', call = '%s'",
                  idcol, afr, eur, eas, callcol))
  d[!is.na(afr)][order(patient)][!duplicated(patient)]
}
anc <- read_ancestry(inp[["admix"]])
src <- "Admixture_by_sample.txt"
if (is.null(anc) || !nrow(anc)) { anc <- read_ancestry(inp[["ucsf"]]); src <- "UCSF_Ancestry_Calls.csv" }
if (is.null(anc) || !nrow(anc)) stop("Neither ancestry file could be parsed; see the column listing above.")
message("Ancestry proportions taken from: ", src, " (", nrow(anc), " patients pan-cancer)")

# ---- one primary tumor per patient, clinical covariates ----------
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
             fifelse(race_raw == "white", "White", "Other/unknown"))]
gl <- num_clin(cd$primary_gleason_grade) + num_clin(cd$secondary_gleason_grade)
tm[, gleason_group := factor(fifelse(gl <= 6, "GS6", fifelse(gl == 7, "GS7", "GS8_10")),
                             levels = c("GS6", "GS7", "GS8_10"))]
tm[, age := if ("age_at_index" %in% names(cd)) num_clin(cd$age_at_index) else num_clin(cd$age_at_diagnosis) / 365.25]

dt <- merge(tm, anc[, .(patient, afr, eur, eas, call)], by = "patient")
dt[, two_way := (afr + eur) >= p[["min_two_way"]]]
fwrite(dt, out[["sheet"]], sep = "\t")
message(sprintf("PRAD primary tumors with genetic ancestry: %d of %d (Gleason group known: %d; age known: %d)",
                nrow(dt), nrow(tm), sum(!is.na(dt$gleason_group)), sum(!is.na(dt$age))))
print(dt[, .(n = .N, median_AFR = round(median(afr), 3), min_AFR = round(min(afr), 3),
             max_AFR = round(max(afr), 3), AFR_ge_0.5 = sum(afr >= 0.5)), by = race])
if (!all(is.na(dt$call))) { message("Published ancestry call by self-reported race:"); print(table(race = dt$race, call = dt$call)) }

g <- ggplot(dt, aes(race, afr, colour = race)) +
  geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.8) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = race_cols, guide = "none") +
  labs(x = "Self-reported race", y = "Genetic African-ancestry fraction (AFR)",
       title = "Self-reported race vs genetic ancestry (TCGA-PRAD)",
       subtitle = paste0("Ancestry proportions: Carrot-Zhang et al. 2020 (", src, ")"))
ggsave(out[["race_plot"]], g, width = 6.5, height = 5, dpi = 150)

# =================================================================
# 2. Genome-wide: expression vs AFR fraction (two-way admixed tumors)
# =================================================================
new_api <- "collection" %in% names(formals(msigdbr::msigdbr))
hdf  <- if (new_api) msigdbr(species = "Homo sapiens", collection = "H") else
                     msigdbr(species = "Homo sapiens", category = "H")
hall <- split(hdf$gene_symbol, hdf$gs_name)
rd <- as.data.frame(rowData(se))

ta <- dt[two_way == TRUE & !is.na(gleason_group) & !is.na(age)]
if (nrow(ta) < 50) stop("Only ", nrow(ta), " tumors with two-way ancestry and complete Gleason/age; check the clinical parsing above.")
ta[, age_s := as.numeric(scale(age))]
n_afr_hi <- sum(ta$afr >= 0.5)
message(sprintf("Genome-wide model: %d two-way (AFR+EUR >= %.2f) tumors, %d with AFR >= 0.5",
                nrow(ta), p[["min_two_way"]], n_afr_hi))
cdt <- data.frame(ta[, .(gleason_group = droplevels(gleason_group), age_s, afr)], row.names = ta$barcode)
d <- DESeqDataSetFromMatrix(assay(se, "unstranded")[, ta$barcode], cdt, design = ~ gleason_group + age_s + afr)
d <- d[rowSums(counts(d) >= p[["min_count"]]) >= max(10, n_afr_hi), ]
d <- DESeq(d, parallel = threads > 1, BPPARAM = bp, quiet = TRUE)
r <- results(d, name = "afr", alpha = fdr)
at <- data.table(gene_id = rownames(d), gene_name = rd[rownames(d), "gene_name"],
                 gene_type = rd[rownames(d), "gene_type"], baseMean = r$baseMean,
                 log2FC_per_full_AFR = r$log2FoldChange, stat = r$stat, pvalue = r$pvalue)
at[, reportable := gene_type %in% c("protein_coding", "lncRNA")]
at[reportable == TRUE & !is.na(pvalue), padj := p.adjust(pvalue, "BH")]
at <- at[order(padj)]
fwrite(at, out[["afr_table"]], sep = "\t")
message(sprintf("Genes associated with AFR fraction (protein-coding + lncRNA, FDR < %.2f): %d",
                fdr, sum(at$padj < fdr, na.rm = TRUE)))
print(head(at[reportable == TRUE, .(gene_name, log2FC_per_full_AFR = round(log2FC_per_full_AFR, 2),
                                    padj = signif(padj, 3))], 15))

ra <- at[reportable == TRUE & !is.na(stat)][order(-abs(stat))][!duplicated(gene_name)]
ha <- fgsea(hall, sort(setNames(ra$stat, ra$gene_name), decreasing = TRUE),
            minSize = 15, maxSize = 500, eps = 0, nproc = threads)[, .(pathway, NES_afr = NES, padj_afr = padj)]
rh <- fread(inp[["race_hallmark"]])[, .(pathway, NES_race = NES, padj_race = padj)]
hh <- merge(rh, ha, by = "pathway")[order(padj_afr)]
fwrite(hh, out[["afr_hallmark"]], sep = "\t")
r_nes <- cor(hh$NES_race, hh$NES_afr)
message(sprintf("Hallmark NES, AFR fraction vs self-reported race (M6): r = %.2f", r_nes))
print(hh[padj_afr < fdr | padj_race < fdr, .(pathway = short(pathway), NES_race = round(NES_race, 2),
          padj_race = signif(padj_race, 2), NES_afr = round(NES_afr, 2), padj_afr = signif(padj_afr, 2))])

hh[, status := fifelse(padj_race < fdr & padj_afr < fdr & sign(NES_race) == sign(NES_afr), "significant in both",
               fifelse(padj_race < fdr, "race only", fifelse(padj_afr < fdr, "AFR only", "neither")))]
g <- ggplot(hh, aes(NES_race, NES_afr, colour = status)) +
  geom_hline(yintercept = 0, colour = "grey70") + geom_vline(xintercept = 0, colour = "grey70") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 2.4) +
  ggrepel::geom_text_repel(data = hh[status == "significant in both"], aes(label = short(pathway)),
                           size = 2.6, colour = "black", max.overlaps = 40) +
  scale_colour_manual(values = c(`significant in both` = "#6A3D9A", `race only` = "#E3A587",
                                 `AFR only` = "#9DB4C0", neither = "grey75"), name = NULL) +
  labs(x = "NES, Black vs White (self-reported, M6)", y = "NES, per unit genetic AFR fraction",
       title = "Self-reported race vs genetic ancestry: Hallmark pathways",
       subtitle = sprintf("Both adjusted for Gleason group and age; r = %.2f", r_nes)) +
  theme(legend.position = "bottom")
ggsave(out[["nes_plot"]], g, width = 7.5, height = 7, dpi = 150)

# =================================================================
# 3. Germline test for the top race-associated genes
# =================================================================
expr_sym <- function(sx) {
  tpm <- assay(sx, "tpm_unstrand"); sym <- rd[rownames(sx), "gene_name"]
  ok <- !is.na(sym) & sym != ""; tpm <- tpm[ok, , drop = FALSE]; sym <- sym[ok]
  o <- order(-rowMeans(tpm)); tpm <- tpm[o, , drop = FALSE]; sym <- sym[o]
  tpm <- tpm[!duplicated(sym), , drop = FALSE]; rownames(tpm) <- sym[!duplicated(sym)]
  log2(tpm + 1)
}
race_tab <- fread(inp[["race_table"]])
top <- race_tab[reportable == TRUE & !is.na(padj)][order(padj)]
top <- head(top[!duplicated(gene_name)]$gene_name, p[["n_genes"]])

ex_t <- expr_sym(se[, dt$barcode])
sp   <- readRDS(inp[["se_paired"]])
nb   <- colnames(sp)[sp$condition == "normal"]
nrm  <- data.table(barcode = nb, patient = substr(nb, 1, 12))[anc[, .(patient, afr)], on = "patient", nomatch = 0]
ex_n <- expr_sym(se[, nrm$barcode])
# Within-group dose-response uses self-identified Black patients WITH substantial African
# ancestry (AFR >= 0.5): the few with near-zero AFR would otherwise act as high-leverage
# points and turn the test back into a two-group comparison.
blk  <- dt[race == "Black" & afr >= 0.5 & !is.na(gleason_group) & !is.na(age)]
message(sprintf("Self-identified Black patients excluded from within-group tests for AFR < 0.5: %d",
                nrow(dt[race == "Black" & afr < 0.5])))
message(sprintf("Germline tests: %d top race genes; tumors n = %d; Black-only tumors n = %d (AFR %.2f-%.2f); normals n = %d (%d with AFR >= 0.5)",
                length(top), nrow(dt[!is.na(gleason_group) & !is.na(age)]), nrow(blk),
                min(blk$afr), max(blk$afr), nrow(nrm), sum(nrm$afr >= 0.5)))

lm_afr <- function(y, d, adj = TRUE) {
  if (all(is.na(y)) || sd(y, na.rm = TRUE) == 0) return(c(NA, NA))
  f <- if (adj) y ~ afr + gleason_group + age else y ~ afr
  dd <- cbind(d, y = y)
  if (adj) dd[, gleason_group := droplevels(gleason_group)]
  co <- tryCatch(summary(lm(f, data = dd))$coefficients, error = function(e) NULL)
  if (is.null(co)) return(c(NA, NA))
  if (!"afr" %in% rownames(co)) return(c(NA, NA))
  co["afr", c("Estimate", "Pr(>|t|)")]
}
dtc <- dt[!is.na(gleason_group) & !is.na(age)]
gtest <- rbindlist(lapply(top, function(gname) {
  yt <- if (gname %in% rownames(ex_t)) ex_t[gname, dtc$barcode] else rep(NA_real_, nrow(dtc))
  yb <- if (gname %in% rownames(ex_t)) ex_t[gname, blk$barcode] else rep(NA_real_, nrow(blk))
  yn <- if (gname %in% rownames(ex_n)) ex_n[gname, nrm$barcode] else rep(NA_real_, nrow(nrm))
  a1 <- lm_afr(yt, dtc); a2 <- lm_afr(yb, blk); a3 <- lm_afr(yn, nrm, adj = FALSE)
  data.table(gene = gname,
             race_log2FC = race_tab[gene_name == gname]$log2FC_shrunk[1],
             beta_all_tumors = a1[1], p_all_tumors = a1[2],
             beta_within_black = a2[1], p_within_black = a2[2],
             beta_normals = a3[1], p_normals = a3[2])
}))
for (cc in c("p_all_tumors", "p_within_black", "p_normals"))
  gtest[, (sub("^p_", "padj_", cc)) := p.adjust(get(cc), "BH")]
gtest[, likely_germline := padj_all_tumors < fdr & sign(beta_all_tumors) == sign(race_log2FC) &
                            ((p_within_black < 0.05 & sign(beta_within_black) == sign(race_log2FC)) |
                             (p_normals < 0.05 & sign(beta_normals) == sign(race_log2FC)))]
fwrite(gtest, out[["germline"]], sep = "\t")
message("Top race-associated genes vs genetic AFR fraction (beta = change in log2 TPM from 0 to 1 AFR):")
print(gtest[, .(gene, race_log2FC = round(race_log2FC, 2), beta_tum = round(beta_all_tumors, 2),
                padj_tum = signif(padj_all_tumors, 2), beta_blk = round(beta_within_black, 2),
                p_blk = signif(p_within_black, 2), beta_norm = round(beta_normals, 2),
                p_norm = signif(p_normals, 2), likely_germline)])

show_g <- head(gtest[order(p_all_tumors)]$gene, 6)
sl <- rbindlist(lapply(show_g, function(gname) rbind(
  data.table(gene = gname, tissue = "Tumor", afr = dt$afr, race = dt$race,
             y = if (gname %in% rownames(ex_t)) ex_t[gname, dt$barcode] else NA_real_),
  data.table(gene = gname, tissue = "Normal", afr = nrm$afr,
             race = dt$race[match(nrm$patient, dt$patient)],
             y = if (gname %in% rownames(ex_n)) ex_n[gname, nrm$barcode] else NA_real_))))
sl[is.na(race), race := "Other/unknown"]
sl[, gene := factor(gene, levels = show_g)]
g <- ggplot(sl, aes(afr, y, colour = race)) +
  geom_point(size = 1, alpha = 0.7) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = "black", linewidth = 0.5) +
  scale_colour_manual(values = race_cols, name = "Self-reported") +
  facet_grid(tissue ~ gene, scales = "free_y") +
  labs(x = "Genetic African-ancestry fraction", y = "log2(TPM + 1)",
       title = "Do the top race-associated genes track genetic ancestry?",
       subtitle = "Top: tumors. Bottom: matched normal tissue. Black line: linear fit") +
  theme(legend.position = "bottom", strip.text = element_text(size = 8))
ggsave(out[["germline_plot"]], g, width = 13, height = 6, dpi = 150)

# =================================================================
# 4. Interferon-gamma score vs AFR fraction
# =================================================================
zt <- t(scale(t(ex_t[rowMeans(ex_t) > 1, ])))
ifng_score <- colMeans(zt[intersect(hall$HALLMARK_INTERFERON_GAMMA_RESPONSE, rownames(zt)), , drop = FALSE])
cs <- fread(inp[["scores"]])
ii <- merge(dtc[, .(barcode, afr, race, gleason_group, age)], cs[, .(barcode, mcp_T_cells, mcp_B_lineage)], by = "barcode")
ii[, ifng := as.numeric(scale(ifng_score[barcode]))]
fit_i <- function(f, d) { co <- summary(lm(f, data = d))$coefficients["afr", ]
  data.table(beta_sd_per_full_AFR = co[["Estimate"]], p = co[["Pr(>|t|)"]]) }
ifn <- rbind(
  cbind(model = "all tumors: + Gleason + age", fit_i(ifng ~ afr + gleason_group + age, ii)),
  cbind(model = "all tumors: + Gleason + age + T/B cells", fit_i(ifng ~ afr + gleason_group + age + mcp_T_cells + mcp_B_lineage, ii)),
  cbind(model = "Black patients with AFR >= 0.5 only: + Gleason + age",
        fit_i(ifng ~ afr + gleason_group + age, ii[race == "Black" & afr >= 0.5])))
fwrite(ifn, out[["ifn"]], sep = "\t")
message("Interferon-gamma score vs AFR fraction (SD change from 0 to 1 AFR):"); print(ifn)

# ---- summary -----------------------------------------------------
v1 <- function(x) if (length(x)) x[1] else NA
summ <- data.table(
  metric = c("ancestry_source", "tumors_with_ancestry", "black_median_afr", "black_min_afr",
             "white_max_afr", "other_unknown_afr_ge_0.5", "two_way_tumors", "two_way_afr_ge_0.5",
             "afr_genes_fdr05", "afr_hallmark_fdr05", "nes_r_afr_vs_race",
             "ifna_nes_afr", "ifna_padj_afr", "ifng_nes_afr", "ifng_padj_afr",
             "top_genes_tested", "likely_germline_genes", "likely_germline_list",
             "normals_with_afr", "normals_afr_ge_0.5", "black_within_group_n", "black_excluded_low_afr",
             "ifng_afr_beta", "ifng_afr_p", "ifng_afr_beta_adj_TB", "ifng_afr_p_adj_TB",
             "ifng_afr_beta_black_only", "ifng_afr_p_black_only"),
  value = c(src, nrow(dt), round(median(dt[race == "Black"]$afr), 3), round(min(dt[race == "Black"]$afr), 3),
            round(max(dt[race == "White"]$afr), 3), sum(dt[race == "Other/unknown"]$afr >= 0.5),
            nrow(ta), n_afr_hi, sum(at$padj < fdr, na.rm = TRUE), sum(hh$padj_afr < fdr),
            round(r_nes, 3),
            v1(round(hh[pathway == "HALLMARK_INTERFERON_ALPHA_RESPONSE"]$NES_afr, 3)),
            v1(signif(hh[pathway == "HALLMARK_INTERFERON_ALPHA_RESPONSE"]$padj_afr, 3)),
            v1(round(hh[pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE"]$NES_afr, 3)),
            v1(signif(hh[pathway == "HALLMARK_INTERFERON_GAMMA_RESPONSE"]$padj_afr, 3)),
            nrow(gtest), sum(gtest$likely_germline, na.rm = TRUE),
            paste(gtest[likely_germline == TRUE]$gene, collapse = ";"),
            nrow(nrm), sum(nrm$afr >= 0.5), nrow(blk), nrow(dt[race == "Black" & afr < 0.5]),
            round(ifn$beta_sd_per_full_AFR[1], 3), signif(ifn$p[1], 3),
            round(ifn$beta_sd_per_full_AFR[2], 3), signif(ifn$p[2], 3),
            round(ifn$beta_sd_per_full_AFR[3], 3), signif(ifn$p[3], 3)))
fwrite(summ, out[["summary"]], sep = "\t")
message("Genetic ancestry summary:"); print(summ)
message("Done.")
