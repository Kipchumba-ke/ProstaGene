#!/usr/bin/env Rscript
# ---------------------------------------------------------------
# 11_subtypes.R — split by molecular subtype (M10)
#
# Subtype labels: TCGA PRAD 2015 (Cell) mutation/fusion classes, added
# by TCGAbiolinks as paper_Subtype (ERG, ETV1, ETV4, FLI1, SPOP, FOXA1,
# IDH1, other) for ~333 tumors. ERG-fusion tumors overexpress ERG, so
# ERG status is also called from expression for every tumor and checked
# against the paper labels.
#
#  1. Subtype frequencies, incl. ERG by genetic African-ancestry fraction
#  2. Paired tumor vs normal, ERG+ vs ERG- pairs: nested DESeq2 design;
#     pathway differences tested with the patient as the unit (Wilcoxon)
#  3. Tumor-only subtype programmes: ~ gleason + age + subtype group
#  4. Progression-free interval by ERG status and subtype group
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(DESeq2)
  library(BiocParallel)
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(survival)
  library(ggplot2)
})

if (exists("snakemake")) {
  inp <- snakemake@input; out <- snakemake@output; p <- snakemake@params
  threads <- snakemake@threads
  con <- file(snakemake@log[[1]], open = "wt"); sink(con); sink(con, type = "message")
} else {
  inp <- list(se_full = "data/processed/prad_se_full.rds",
              se_paired = "data/processed/prad_se_paired.rds",
              ancestry = "results/genetic_ancestry/ancestry_sample_sheet.tsv",
              survival = "results/survival/survival_dataset.tsv")
  s <- "results/subtypes"
  out <- list(sheet = file.path(s, "subtype_sample_sheet.tsv"),
              erg_plot = file.path(s, "erg_expression_call.png"),
              counts_plot = file.path(s, "subtype_by_ancestry.png"),
              erg_ancestry = file.path(s, "erg_vs_ancestry.tsv"),
              paired = file.path(s, "paired_by_erg_hallmark.tsv"),
              paired_plot = file.path(s, "paired_by_erg_hallmark.png"),
              subtype_hallmark = file.path(s, "subtype_hallmark.tsv"),
              subtype_plot = file.path(s, "subtype_hallmark_heatmap.png"),
              cox = file.path(s, "subtype_cox.tsv"),
              summary = file.path(s, "subtype_summary.tsv"))
  p <- list(fdr = 0.05, min_count = 10, seed = 42)
  threads <- 4
}
for (f in unlist(out)) dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 11))
fdr <- p[["fdr"]]
set.seed(p[["seed"]])
bp <- if (threads > 1) MulticoreParam(threads) else SerialParam()
short <- function(x) gsub("_", " ", sub("^HALLMARK_", "", x))
num_clin <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(x))))
grp_cols <- c(ERG = "#5E8C7A", ETS_other = "#8FB8A8", SPOP = "#C99A6B", FOXA1 = "#D9B98C", Other = "#B7B7AE")

new_api <- "collection" %in% names(formals(msigdbr::msigdbr))
hdf  <- if (new_api) msigdbr(species = "Homo sapiens", collection = "H") else
                     msigdbr(species = "Homo sapiens", category = "H")
hall <- split(hdf$gene_symbol, hdf$gs_name)
gsea_h <- function(stat, symbol) {
  d <- data.table(stat, symbol)[!is.na(stat) & !is.na(symbol) & symbol != ""]
  d <- d[order(-abs(stat))][!duplicated(symbol)]
  fgsea(hall, sort(setNames(d$stat, d$symbol), decreasing = TRUE),
        minSize = 15, maxSize = 500, eps = 0, nproc = threads)[, .(pathway, NES, padj)]
}

se <- readRDS(inp[["se_full"]])
rd <- as.data.frame(rowData(se))
erg_id <- rownames(rd)[which(rd$gene_name == "ERG")[1]]
log_tpm <- function(ids, cols) log2(assay(se, "tpm_unstrand")[ids, cols, drop = FALSE] + 1)

# =================================================================
# 1. Subtype labels + expression-based ERG call
# =================================================================
bc <- colnames(se)
tm <- data.table(barcode = bc, patient = substr(bc, 1, 12), code = substr(bc, 14, 15),
                 vial = substr(bc, 16, 16), analyte = substr(bc, 20, 20), plate = substr(bc, 22, 25))[code == "01"]
tm[, arank := match(analyte, c("H", "R", "T"), nomatch = 4L)]
tm <- tm[order(patient, vial, arank, -rank(plate))][!duplicated(patient)]
cd <- as.data.frame(colData(se))[tm$barcode, ]

scol <- intersect(c("paper_Subtype", "paper_subtype", "paper_SUBTYPE"), names(cd))[1]
if (is.na(scol)) {
  cand <- grep("^paper_", names(cd), value = TRUE)
  hit  <- cand[sapply(cand, function(cn) any(grepl("ERG", as.character(cd[[cn]]))))]
  scol <- hit[1]
}
message("Subtype column: ", scol)
raw <- if (!is.na(scol)) as.character(cd[[scol]]) else rep(NA_character_, nrow(cd))
print(table(raw, useNA = "ifany"))
up <- toupper(raw)
tm[, paper_subtype := raw]
tm[, subtype_group := fifelse(grepl("ERG", up), "ERG",
                      fifelse(grepl("ETV1|ETV4|FLI1", up), "ETS_other",
                      fifelse(grepl("SPOP", up), "SPOP",
                      fifelse(grepl("FOXA1", up), "FOXA1",
                      fifelse(!is.na(raw) & raw != "" & !grepl("^NA$", raw), "Other", NA_character_)))))]
tm[, subtype_group := factor(subtype_group, levels = names(grp_cols))]

# ERG expression in tumors, normals for reference
nb <- bc[substr(bc, 14, 15) == "11"]
tm[, erg_log2tpm := as.numeric(log_tpm(erg_id, barcode))]
erg_norm <- as.numeric(log_tpm(erg_id, nb))
km <- kmeans(tm$erg_log2tpm, centers = quantile(tm$erg_log2tpm, c(0.2, 0.9)))
hi_cluster <- which.max(km$centers)
thr <- mean(c(max(tm$erg_log2tpm[km$cluster != hi_cluster]), min(tm$erg_log2tpm[km$cluster == hi_cluster])))
tm[, erg_expr_positive := erg_log2tpm > thr]
lab <- tm[!is.na(subtype_group)]
conc <- mean(lab$erg_expr_positive == (lab$subtype_group == "ERG"))
message(sprintf("ERG expression threshold (k-means, 2 groups): log2(TPM+1) = %.2f; ERG+ by expression: %d of %d tumors; agreement with paper ERG label: %.1f%% (n = %d)",
                thr, sum(tm$erg_expr_positive), nrow(tm), 100 * conc, nrow(lab)))
print(table(paper_ERG = lab$subtype_group == "ERG", expression_ERG = lab$erg_expr_positive))

g <- ggplot(rbind(data.table(tissue = "Tumor", v = tm$erg_log2tpm), data.table(tissue = "Normal", v = erg_norm)),
            aes(v, fill = tissue)) +
  geom_histogram(bins = 50, alpha = 0.75, position = "identity") +
  geom_vline(xintercept = thr, linetype = "dashed") +
  scale_fill_manual(values = c(Tumor = "#5E8C7A", Normal = "#C9C3B6")) +
  labs(x = "ERG log2(TPM + 1)", y = "Samples", fill = NULL,
       title = "Calling ERG-fusion status from ERG expression",
       subtitle = sprintf("Dashed line: threshold %.2f; %.0f%% agreement with TCGA 2015 labels (n = %d)", thr, 100 * conc, nrow(lab)))
ggsave(out[["erg_plot"]], g, width = 7, height = 4.5, dpi = 150)

# Clinical + ancestry
gl <- num_clin(cd$primary_gleason_grade) + num_clin(cd$secondary_gleason_grade)
tm[, gleason_group := factor(fifelse(gl <= 6, "GS6", fifelse(gl == 7, "GS7", "GS8_10")), levels = c("GS6", "GS7", "GS8_10"))]
tm[, age := if ("age_at_index" %in% names(cd)) num_clin(cd$age_at_index) else num_clin(cd$age_at_diagnosis) / 365.25]
anc <- fread(inp[["ancestry"]])[, .(patient, afr, race)]
tm <- merge(tm, anc, by = "patient", all.x = TRUE)
fwrite(tm, out[["sheet"]], sep = "\t")
message("Subtype groups (TCGA 2015 labels):"); print(table(tm$subtype_group, useNA = "ifany"))

# ERG vs African ancestry
ea <- tm[!is.na(afr) & !is.na(gleason_group) & !is.na(age)]
ea[, afr_group := factor(fifelse(afr >= 0.5, "AFR >= 0.5", "AFR < 0.5"), levels = c("AFR < 0.5", "AFR >= 0.5"))]
tab_e <- ea[, .(n = .N, erg_pos = sum(erg_expr_positive), pct_erg = round(100 * mean(erg_expr_positive), 1)), by = afr_group][order(afr_group)]
ft <- fisher.test(table(ea$afr_group, ea$erg_expr_positive))
lr <- summary(glm(erg_expr_positive ~ afr + gleason_group + age, family = binomial, data = ea))$coefficients["afr", ]
erg_anc <- cbind(tab_e, fisher_or = c(NA, unname(ft$estimate)), fisher_p = c(NA, ft$p.value),
                 logit_or_per_full_afr = c(NA, exp(lr[["Estimate"]])), logit_p = c(NA, lr[["Pr(>|z|)"]]))
fwrite(erg_anc, out[["erg_ancestry"]], sep = "\t")
message("ERG-positive (expression) by genetic African ancestry:"); print(erg_anc)

cnt <- tm[!is.na(subtype_group) & !is.na(afr)][, afr_group := fifelse(afr >= 0.5, "AFR >= 0.5", "AFR < 0.5")]
cnt <- cnt[, .N, by = .(afr_group, subtype_group)][, pct := 100 * N / sum(N), by = afr_group]
g <- ggplot(cnt, aes(afr_group, pct, fill = subtype_group)) +
  geom_col(width = 0.6) +
  scale_fill_manual(values = grp_cols, name = "Subtype") +
  labs(x = "Genetic African-ancestry fraction", y = "% of tumors",
       title = "Molecular subtype by genetic ancestry (TCGA 2015 labels)",
       subtitle = paste0("n = ", paste(cnt[, .(n = sum(N)), by = afr_group][order(afr_group)]$n, collapse = " vs ")))
ggsave(out[["counts_plot"]], g, width = 6, height = 5, dpi = 150)

# =================================================================
# 2. Paired tumor vs normal, ERG+ vs ERG- pairs
# =================================================================
sp <- readRDS(inp[["se_paired"]])
cdp <- data.frame(patient = as.character(sp$patient), condition = as.character(sp$condition), row.names = colnames(sp))
erg_by_patient <- setNames(tm$erg_expr_positive, tm$patient)
pair_bc <- data.table(barcode = colnames(sp), patient = cdp$patient)
tum_bc  <- colnames(sp)[cdp$condition == "tumor"]
erg_pair <- setNames(as.numeric(log_tpm(erg_id, tum_bc)) > thr, substr(tum_bc, 1, 12))
cdp$erg <- factor(ifelse(erg_pair[cdp$patient], "ERGpos", "ERGneg"), levels = c("ERGneg", "ERGpos"))
cdp$condition <- factor(cdp$condition, levels = c("normal", "tumor"))
cdp$pair_n <- factor(ave(cdp$patient, cdp$erg, FUN = function(x) as.integer(factor(x))))
n_pos <- length(unique(cdp$patient[cdp$erg == "ERGpos"])); n_neg <- length(unique(cdp$patient[cdp$erg == "ERGneg"]))
message(sprintf("Paired set: %d ERG+ pairs, %d ERG- pairs", n_pos, n_neg))

dds <- DESeqDataSetFromMatrix(assay(sp, "unstranded"), cdp, design = ~ 1)
dds <- dds[rowSums(counts(dds) >= p[["min_count"]]) >= min(n_pos, n_neg), ]
mm <- model.matrix(~ erg + erg:pair_n + erg:condition, cdp)
mm <- mm[, colSums(mm != 0) > 0]; colnames(mm) <- make.names(colnames(mm))
dds <- DESeq(dds, full = mm, parallel = threads > 1, BPPARAM = bp, quiet = TRUE)
rn <- resultsNames(dds)
cpos <- grep("^ergERGpos\\.conditiontumor$", rn, value = TRUE); cneg <- grep("^ergERGneg\\.conditiontumor$", rn, value = TRUE)
r_pos <- results(dds, name = cpos, alpha = fdr); r_neg <- results(dds, name = cneg, alpha = fdr)
r_dif <- results(dds, contrast = list(cpos, cneg), alpha = fdr)
gnm <- rd[rownames(dds), "gene_name"]
n_dif <- sum(r_dif$padj < fdr, na.rm = TRUE)
message(sprintf("Tumor-vs-normal DEGs (padj < %.2f, |LFC| >= 1): ERG+ %d, ERG- %d; genes whose tumor effect differs by ERG status: %d",
                fdr, sum(r_pos$padj < fdr & abs(r_pos$log2FoldChange) >= 1, na.rm = TRUE),
                sum(r_neg$padj < fdr & abs(r_neg$log2FoldChange) >= 1, na.rm = TRUE), n_dif))
top_dif <- data.table(gene = gnm, ERGpos_lfc = r_pos$log2FoldChange, ERGneg_lfc = r_neg$log2FoldChange, padj = r_dif$padj)[order(padj)]
print(head(top_dif[!is.na(gene)], 12))

# Patient-level pathway test: per-pair tumor-minus-normal change, ERG+ vs ERG-
vs <- assay(vst(dds, blind = TRUE)); sym <- gnm
ok <- !is.na(sym) & !duplicated(sym); vs <- vs[ok, ]; rownames(vs) <- sym[ok]
pats <- unique(cdp$patient)
delta <- sapply(pats, function(pt) vs[, rownames(cdp)[cdp$patient == pt & cdp$condition == "tumor"]] -
                                   vs[, rownames(cdp)[cdp$patient == pt & cdp$condition == "normal"]])
is_pos <- erg_pair[pats]
sets <- Filter(function(g) length(intersect(g, rownames(delta))) >= 15, hall)
sc <- t(sapply(sets, function(g) colMeans(delta[intersect(g, rownames(delta)), , drop = FALSE])))
pp <- data.table(pathway = rownames(sc),
                 median_ERGpos = apply(sc[, is_pos, drop = FALSE], 1, median),
                 median_ERGneg = apply(sc[, !is_pos, drop = FALSE], 1, median),
                 p_wilcox = apply(sc, 1, function(x) wilcox.test(x[is_pos], x[!is_pos])$p.value))
pp[, padj_wilcox := p.adjust(p_wilcox, "BH")]
hp <- gsea_h(r_pos$stat, gnm)[, .(pathway, NES_ERGpos = NES, padj_ERGpos = padj)]
hn <- gsea_h(r_neg$stat, gnm)[, .(pathway, NES_ERGneg = NES, padj_ERGneg = padj)]
pp <- merge(merge(pp, hp, by = "pathway"), hn, by = "pathway")[order(p_wilcox)]
fwrite(pp, out[["paired"]], sep = "\t")
message("Hallmark tumor effect by ERG status (patient-level Wilcoxon on per-pair change):")
print(head(pp[, .(pathway = short(pathway), NES_ERGpos = round(NES_ERGpos, 2), NES_ERGneg = round(NES_ERGneg, 2),
                  p_wilcox = signif(p_wilcox, 2), padj_wilcox = signif(padj_wilcox, 2))], 12))

pl <- melt(pp[, .(pathway, `ERG+` = NES_ERGpos, `ERG-` = NES_ERGneg, padj_wilcox)],
           id.vars = c("pathway", "padj_wilcox"), variable.name = "group", value.name = "NES")
pl[, label := factor(short(pathway), levels = short(pp[order(pp$NES_ERGpos + pp$NES_ERGneg)]$pathway))]
g <- ggplot(pl, aes(NES, label)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_line(aes(group = label), colour = "grey80") +
  geom_point(aes(colour = group, shape = padj_wilcox < fdr), size = 2.4) +
  scale_colour_manual(values = c(`ERG+` = "#5E8C7A", `ERG-` = "#C99A6B"), name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), labels = c("difference ns", "differs (patient-level FDR < 0.05)"), name = NULL) +
  labs(x = "NES, tumor vs matched normal", y = NULL,
       title = "Tumor-vs-normal pathway changes in ERG+ and ERG- tumors",
       subtitle = sprintf("%d ERG+ and %d ERG- patient pairs", n_pos, n_neg)) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 7))
ggsave(out[["paired_plot"]], g, width = 8.5, height = 10, dpi = 150)

# =================================================================
# 3. Tumor-only subtype programmes (paper labels)
# =================================================================
ts <- tm[!is.na(subtype_group) & !is.na(gleason_group) & !is.na(age)]
ts[, subtype_group := relevel(droplevels(subtype_group), ref = "Other")]
ts[, age_s := as.numeric(scale(age))]
message("Tumor-only subtype model (reference = Other): "); print(table(ts$subtype_group))
cdt <- data.frame(ts[, .(subtype_group, gleason_group = droplevels(gleason_group), age_s)], row.names = ts$barcode)
dt <- DESeqDataSetFromMatrix(assay(se, "unstranded")[, ts$barcode], cdt, design = ~ gleason_group + age_s + subtype_group)
dt <- dt[rowSums(counts(dt) >= p[["min_count"]]) >= 10, ]
dt <- DESeq(dt, parallel = threads > 1, BPPARAM = bp, quiet = TRUE)
coefs <- grep("^subtype_group_", resultsNames(dt), value = TRUE)
gn2 <- rd[rownames(dt), "gene_name"]
sh <- rbindlist(lapply(coefs, function(cf) {
  r <- results(dt, name = cf, alpha = fdr)
  top <- head(data.table(gene = gn2, lfc = r$log2FoldChange, padj = r$padj)[!is.na(padj)][order(padj)], 3)
  message(sprintf("  %s: %d genes padj < %.2f; top: %s", cf, sum(r$padj < fdr, na.rm = TRUE), fdr,
                  paste(sprintf("%s (%.1f)", top$gene, top$lfc), collapse = ", ")))
  cbind(contrast = sub("^subtype_group_", "", cf), gsea_h(r$stat, gn2))
}))
fwrite(sh, out[["subtype_hallmark"]], sep = "\t")
keep_pw <- sh[padj < fdr, unique(pathway)]
hm <- sh[pathway %in% keep_pw]
hm[, label := short(pathway)]
hm[, contrast := gsub("_vs_Other", "", contrast)]
g <- ggplot(hm, aes(contrast, label, fill = NES)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(padj < fdr, "•", "")), size = 4) +
  scale_fill_gradient2(low = "#C99A6B", mid = "#F7F5F0", high = "#4F7F6F", midpoint = 0) +
  labs(x = "Subtype vs 'Other'", y = NULL, title = "Hallmark programmes by molecular subtype (tumors only)",
       subtitle = "Adjusted for Gleason group and age; dot = FDR < 0.05") +
  theme(axis.text.y = element_text(size = 7))
ggsave(out[["subtype_plot"]], g, width = 7, height = max(4, 0.22 * length(keep_pw) + 2), dpi = 150)

# =================================================================
# 4. Progression-free interval by ERG status / subtype
# =================================================================
sv <- fread(inp[["survival"]])[, .(barcode, pfi, pfi_months, t_stage)]
cx <- merge(tm, sv, by = "barcode")[!is.na(gleason_group) & !is.na(t_stage) & !is.na(age)]
cx[, gleason_group := factor(gleason_group, levels = c("GS6", "GS7", "GS8_10"))]
# Subtype groups with fewer than 3 progression events give infinite hazard ratios
# (e.g. FOXA1: 9 tumors, 0 events); merge them into "Other" for the Cox model.
ev <- cx[!is.na(subtype_group), .(events = sum(pfi), n = .N), by = subtype_group]
small <- as.character(ev[events < 3 & subtype_group != "Other", subtype_group])
cx[, cox_group := as.character(subtype_group)]
cx[cox_group %in% small, cox_group := "Other"]
message("Progression events by subtype group:"); print(ev[order(subtype_group)])
if (length(small)) message("Merged into 'Other' for the Cox model (< 3 events): ", paste(small, collapse = ", "))
cox_row <- function(f, d, term) {
  fit <- coxph(f, data = d); s <- summary(fit)
  rows <- grep(term, rownames(s$coefficients), value = TRUE)
  data.table(term = rows, HR = s$conf.int[rows, "exp(coef)"], lower95 = s$conf.int[rows, "lower .95"],
             upper95 = s$conf.int[rows, "upper .95"], p = s$coefficients[rows, "Pr(>|z|)"], n = s$n, events = s$nevent)
}
cox <- rbind(
  cbind(model = "ERG (expression) unadjusted", cox_row(Surv(pfi_months, pfi) ~ erg_expr_positive, cx, "erg")),
  cbind(model = "ERG (expression) + Gleason + T + age",
        cox_row(Surv(pfi_months, pfi) ~ erg_expr_positive + gleason_group + t_stage + age, cx, "erg")),
  cbind(model = "Subtype group (paper) + Gleason + T + age",
        cox_row(Surv(pfi_months, pfi) ~ relevel(factor(cox_group), ref = "Other") + gleason_group + t_stage + age,
                cx[!is.na(cox_group)], "cox_group")))
cox[, term := gsub("relevel\\(factor\\(cox_group\\), ref = \"Other\"\\)", "subtype:", term)]
fwrite(cox, out[["cox"]], sep = "\t")
message("Progression-free interval:"); print(cox[, .(model, term, HR = round(HR, 2), CI = sprintf("%.2f-%.2f", lower95, upper95), p = signif(p, 2), n, events)])

# ---- summary -----------------------------------------------------
v1 <- function(x) if (length(x)) x[1] else NA
summ <- data.table(
  metric = c("subtype_column", "tumors", "tumors_with_paper_subtype", "erg_threshold_log2tpm",
             "erg_pos_expression", "erg_call_agreement_pct", "erg_pct_afr_low", "erg_pct_afr_high",
             "n_afr_high", "erg_afr_logit_or", "erg_afr_logit_p", "paired_erg_pos", "paired_erg_neg",
             "genes_tumor_effect_differs_by_erg", "pathways_differ_by_erg_fdr05",
             "erg_pfi_hr_adj", "erg_pfi_p_adj", "cox_groups_merged_into_other"),
  value = c(scol, nrow(tm), sum(!is.na(tm$subtype_group)), round(thr, 3), sum(tm$erg_expr_positive),
            round(100 * conc, 1), v1(tab_e[afr_group == "AFR < 0.5", pct_erg]), v1(tab_e[afr_group == "AFR >= 0.5", pct_erg]),
            v1(tab_e[afr_group == "AFR >= 0.5", n]), round(exp(lr[["Estimate"]]), 3), signif(lr[["Pr(>|z|)"]], 3),
            n_pos, n_neg, n_dif, sum(pp$padj_wilcox < fdr),
            round(cox[model == "ERG (expression) + Gleason + T + age", HR], 3),
            signif(cox[model == "ERG (expression) + Gleason + T + age", p], 3),
            if (length(small)) paste(small, collapse = ";") else "none"))
fwrite(summ, out[["summary"]], sep = "\t")
message("Subtype summary:"); print(summ)
message("Done.")
