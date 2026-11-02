#!/usr/bin/env Rscript
# ---------------------------------------------------------------
# 09_survival.R — association with progression (M8, part 2)
#
# Outcome: progression-free interval (PFI) from the TCGA Pan-Cancer
# Clinical Data Resource (TCGA-CDR; Liu et al. 2018), the endpoint
# recommended for TCGA-PRAD (few deaths, so OS is underpowered).
#
# One primary tumor per patient. Per-tumor scores (z across tumors):
#   - TCGA tumor signature (top 100 up minus top 100 down DEGs)
#   - Hallmark set scores (mean gene z): proliferation, MYC, androgen,
#     interferon, EMT, myogenesis
#   - Composition: ESTIMATE stromal/immune, MCP-counter populations
# Cox models per SD: univariable, and adjusted for Gleason group,
# pathological T stage and age. C-index gain over clinical-only.
# Kaplan-Meier curves (median split) for four key scores.
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(data.table)
  library(survival)
  library(readxl)
  library(msigdbr)
  library(ggplot2)
})

if (exists("snakemake")) {
  inp <- snakemake@input; out <- snakemake@output; p <- snakemake@params
  con <- file(snakemake@log[[1]], open = "wt"); sink(con); sink(con, type = "message")
} else {
  inp <- list(se_full = "data/processed/prad_se_full.rds",
              cdr = "data/raw/clinical/TCGA-CDR-SupplementalTableS1.xlsx",
              de = "results/de/deseq2_all_genes.tsv",
              scores = "results/deconv/cell_scores.tsv")
  s <- "results/survival"
  out <- list(dataset = file.path(s, "survival_dataset.tsv"),
              cox = file.path(s, "cox_results.tsv"),
              km = file.path(s, "km_plots.png"),
              forest = file.path(s, "forest_plot.png"),
              summary = file.path(s, "survival_summary.tsv"))
  p <- list(fdr = 0.05, sig_size = 100)
}
for (f in unlist(out)) dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 11))
fdr <- p[["fdr"]]
zs <- function(x) as.numeric(scale(x))
num <- function(x) suppressWarnings(as.numeric(gsub("[^0-9.]", "", as.character(x))))

# ---- 1. one primary tumor per patient + clinical covariates ------
se <- readRDS(inp[["se_full"]])
bc <- colnames(se)
tm <- data.table(barcode = bc, patient = substr(bc, 1, 12), code = substr(bc, 14, 15),
                 vial = substr(bc, 16, 16), analyte = substr(bc, 20, 20), plate = substr(bc, 22, 25))[code == "01"]
tm[, arank := match(analyte, c("H", "R", "T"), nomatch = 4L)]
tm <- tm[order(patient, vial, arank, -rank(plate))][!duplicated(patient)]
cd <- as.data.frame(colData(se))[tm$barcode, ]

gl <- if (all(c("primary_gleason_grade", "secondary_gleason_grade") %in% names(cd)))
  num(cd$primary_gleason_grade) + num(cd$secondary_gleason_grade) else rep(NA_real_, nrow(cd))
tm[, gleason_group := factor(fifelse(gl <= 6, "GS6", fifelse(gl == 7, "GS7", "GS8_10")),
                             levels = c("GS6", "GS7", "GS8_10"))]
tcol <- intersect(c("ajcc_pathologic_t", "paper_Pathologic_T", "pathologic_T"), names(cd))[1]
tm[, t_stage := if (!is.na(tcol)) {
  t <- toupper(substr(as.character(cd[[tcol]]), 1, 2))
  factor(fifelse(t == "T4", "T3_4", fifelse(t == "T3", "T3_4", fifelse(t == "T2", "T2", NA_character_))),
         levels = c("T2", "T3_4"))
} else factor(NA_character_)]
tm[, age := if ("age_at_index" %in% names(cd)) num(cd$age_at_index) else num(cd$age_at_diagnosis) / 365.25]
race_raw <- tolower(as.character(cd$race))
tm[, race := fifelse(race_raw == "black or african american", "Black",
             fifelse(race_raw == "white", "White", "Other/unknown"))]
message("T stage taken from column: ", tcol)

# ---- 2. TCGA-CDR progression-free interval -----------------------
sheets <- excel_sheets(inp[["cdr"]])
cdr <- NULL
for (sh in sheets) {
  x <- as.data.table(read_excel(inp[["cdr"]], sheet = sh, na = c("NA", "#N/A", "[Not Available]")))
  if (all(c("bcr_patient_barcode", "PFI", "PFI.time") %in% names(x))) { cdr <- x; break }
}
if (is.null(cdr)) stop("No sheet with bcr_patient_barcode / PFI / PFI.time in ", inp[["cdr"]])
cdr <- cdr[type == "PRAD", .(patient = bcr_patient_barcode, pfi = as.numeric(PFI), pfi_days = as.numeric(PFI.time))]
dt <- merge(tm, cdr, by = "patient")
dt <- dt[!is.na(pfi) & !is.na(pfi_days) & pfi_days > 0]
dt[, pfi_months := pfi_days / 30.44]
message(sprintf("Tumors with PFI: %d | progression events: %d | median follow-up %.0f months",
                nrow(dt), sum(dt$pfi), median(dt$pfi_months)))

# ---- 3. expression scores ----------------------------------------
rd  <- as.data.frame(rowData(se))
tpm <- assay(se, "tpm_unstrand")[, dt$barcode]
sym <- rd$gene_name
ok  <- !is.na(sym) & sym != ""
tpm <- tpm[ok, ]; sym <- sym[ok]
o   <- order(-rowMeans(tpm)); tpm <- tpm[o, ]; sym <- sym[o]
tpm <- tpm[!duplicated(sym), ]; rownames(tpm) <- sym[!duplicated(sym)]
ex  <- log2(tpm + 1)
ex  <- ex[rowMeans(ex) > 1, ]                         # drop barely expressed genes
z   <- t(scale(t(ex)))
set_score <- function(genes) { g <- intersect(genes, rownames(z)); colMeans(z[g, , drop = FALSE]) }

de <- fread(inp[["de"]])[!is.na(gene_name) & gene_name != ""][order(padj)]
up   <- head(de[direction == "up"]$gene_name,   p[["sig_size"]])
down <- head(de[direction == "down"]$gene_name, p[["sig_size"]])

new_api <- "collection" %in% names(formals(msigdbr::msigdbr))
hdf  <- if (new_api) msigdbr(species = "Homo sapiens", collection = "H") else
                     msigdbr(species = "Homo sapiens", category = "H")
hall <- split(hdf$gene_symbol, hdf$gs_name)
hsets <- c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_MYC_TARGETS_V1",
           "HALLMARK_ANDROGEN_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
           "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
           "HALLMARK_MYOGENESIS")

dt[, tcga_tumor_signature := zs(set_score(up) - set_score(down))]
for (h in hsets) dt[, (sub("^HALLMARK_", "hm_", h)) := zs(set_score(hall[[h]]))]

cs <- fread(inp[["scores"]])
comp_cols <- intersect(c("estimate_stromal", "estimate_immune", "mcp_Fibroblasts",
                         "mcp_CD8_T_cells", "mcp_Cytotoxic_lymphocytes", "smooth_muscle"), names(cs))
dt <- merge(dt, cs[, c("barcode", comp_cols), with = FALSE], by = "barcode")
for (cc in comp_cols) dt[, (cc) := zs(get(cc))]
score_cols <- c("tcga_tumor_signature", sub("^HALLMARK_", "hm_", hsets), comp_cols)
fwrite(dt, out[["dataset"]], sep = "\t")

# ---- 4. Cox models -----------------------------------------------
clin_terms <- c("gleason_group", "t_stage", "age")
clin_terms <- clin_terms[sapply(clin_terms, function(v) uniqueN(na.omit(dt[[v]])) > 1)]
cc_set <- dt[complete.cases(dt[, c("pfi", "pfi_months", clin_terms, score_cols), with = FALSE])]
cc_set[, age := zs(age)]
message(sprintf("Complete-case set for Cox models: %d tumors, %d events; clinical covariates: %s",
                nrow(cc_set), sum(cc_set$pfi), paste(clin_terms, collapse = ", ")))
print(cc_set[, .(n = .N, events = sum(pfi)), by = gleason_group][order(gleason_group)])

f_clin <- as.formula(paste("Surv(pfi_months, pfi) ~", paste(clin_terms, collapse = " + ")))
fit_clin <- coxph(f_clin, data = cc_set)
c_clin <- unname(summary(fit_clin)$concordance[1])
message(sprintf("Clinical-only model C-index: %.3f", c_clin))

cox_row <- function(s, adjusted) {
  rhs <- if (adjusted) paste(c(s, clin_terms), collapse = " + ") else s
  fit <- coxph(as.formula(paste("Surv(pfi_months, pfi) ~", rhs)), data = cc_set)
  co  <- summary(fit)$coefficients[s, ]
  ci  <- summary(fit)$conf.int[s, ]
  data.table(score = s, model = if (adjusted) "adjusted" else "univariable",
             HR_per_SD = ci[["exp(coef)"]], lower95 = ci[["lower .95"]], upper95 = ci[["upper .95"]],
             p = co[["Pr(>|z|)"]], c_index = unname(summary(fit)$concordance[1]),
             ph_global_p = tryCatch(cox.zph(fit)$table["GLOBAL", "p"], error = function(e) NA_real_))
}
cox <- rbindlist(c(lapply(score_cols, cox_row, adjusted = FALSE), lapply(score_cols, cox_row, adjusted = TRUE)))
cox[, padj := p.adjust(p, "BH"), by = model]
cox[model == "adjusted", c_gain_vs_clinical := c_index - c_clin]
cox <- cox[order(model, p)]
fwrite(cox, out[["cox"]], sep = "\t")
message("Cox models (HR per SD; adjusted = + ", paste(clin_terms, collapse = " + "), "):")
print(cox[, .(score, model, HR = round(HR_per_SD, 2), CI = sprintf("%.2f-%.2f", lower95, upper95),
              padj = signif(padj, 3), C = round(c_index, 3))])

# Exploratory: progression by race (few events in Black patients)
rt <- dt[race %in% c("Black", "White"), .(n = .N, events = sum(pfi)), by = race]
message("PFI events by self-reported race:"); print(rt)

# ---- 5. Kaplan-Meier, median split -------------------------------
km_scores <- intersect(c("tcga_tumor_signature", "hm_E2F_TARGETS", "hm_INTERFERON_GAMMA_RESPONSE",
                         "estimate_stromal"), score_cols)
km <- rbindlist(lapply(km_scores, function(s) {
  d <- copy(dt); d[, grp := factor(fifelse(get(s) > median(get(s)), "high", "low"), levels = c("low", "high"))]
  sf <- survfit(Surv(pfi_months, pfi) ~ grp, data = d)
  lr <- survdiff(Surv(pfi_months, pfi) ~ grp, data = d)
  lp <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
  st <- rep(sub("^grp=", "", names(sf$strata)), sf$strata)
  k  <- data.table(score = s, group = st, time = sf$time, surv = sf$surv)
  k  <- rbind(data.table(score = s, group = c("low", "high"), time = 0, surv = 1), k)
  k[, panel := sprintf("%s\nlog-rank p = %s", s, signif(lp, 2))]
  k
}))
km[, group := factor(group, levels = c("low", "high"))]
g <- ggplot(km, aes(time, surv, colour = group)) +
  geom_step(linewidth = 0.8) +
  scale_colour_manual(values = c(low = "#2E6F9E", high = "#C8553D"), name = "Score (median split)") +
  facet_wrap(~ panel, nrow = 1) +
  coord_cartesian(ylim = c(0.4, 1)) +
  labs(x = "Months", y = "Progression-free probability",
       title = "Progression-free interval by score (TCGA-PRAD primary tumors)",
       subtitle = sprintf("%d tumors, %d progression events; unadjusted", nrow(dt), sum(dt$pfi))) +
  theme(legend.position = "bottom", strip.text = element_text(size = 8))
ggsave(out[["km"]], g, width = 12, height = 4.6, dpi = 150)

# ---- 6. forest plot (adjusted HRs) -------------------------------
fp <- cox[model == "adjusted"]
fp[, label := factor(score, levels = score[order(HR_per_SD)])]
fp[, sig := padj < fdr]
g <- ggplot(fp, aes(HR_per_SD, label, colour = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbar(aes(xmin = lower95, xmax = upper95), width = 0.2, orientation = "y") +
  geom_point(size = 2.4) +
  scale_x_log10() +
  scale_colour_manual(values = c(`TRUE` = "#C8553D", `FALSE` = "grey40"),
                      labels = c(`TRUE` = paste0("FDR < ", fdr), `FALSE` = "ns"), name = NULL) +
  labs(x = "Hazard ratio per SD (log scale)", y = NULL,
       title = "Scores and progression-free interval",
       subtitle = paste0("Cox models adjusted for ", paste(clin_terms, collapse = ", "),
                         sprintf("; %d tumors, %d events", nrow(cc_set), sum(cc_set$pfi)))) +
  theme(legend.position = "bottom")
ggsave(out[["forest"]], g, width = 8, height = 6, dpi = 150)

# ---- summary -----------------------------------------------------
best <- cox[model == "adjusted"][order(p)][1]
v1 <- function(x) if (length(x)) x[1] else NA
summ <- data.table(
  metric = c("tumors_with_pfi", "pfi_events", "median_followup_months", "cox_complete_cases", "cox_events",
             "c_index_clinical", "adjusted_scores_fdr05", "univariable_scores_fdr05",
             "best_adjusted_score", "best_adjusted_hr", "best_adjusted_padj", "best_adjusted_c_gain",
             "signature_adjusted_hr", "signature_adjusted_padj",
             "events_black", "n_black", "events_white", "n_white"),
  value = c(nrow(dt), sum(dt$pfi), round(median(dt$pfi_months), 1), nrow(cc_set), sum(cc_set$pfi),
            round(c_clin, 3), sum(cox[model == "adjusted"]$padj < fdr), sum(cox[model == "univariable"]$padj < fdr),
            best$score, round(best$HR_per_SD, 3), signif(best$padj, 3), round(best$c_gain_vs_clinical, 3),
            v1(round(cox[model == "adjusted" & score == "tcga_tumor_signature", HR_per_SD], 3)),
            v1(signif(cox[model == "adjusted" & score == "tcga_tumor_signature", padj], 3)),
            v1(rt[race == "Black", events]), v1(rt[race == "Black", n]),
            v1(rt[race == "White", events]), v1(rt[race == "White", n])))
fwrite(summ, out[["summary"]], sep = "\t")
message("Survival summary:"); print(summ)
message("Done.")
