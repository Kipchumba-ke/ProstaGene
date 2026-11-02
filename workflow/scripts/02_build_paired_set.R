#!/usr/bin/env Rscript
# ---------------------------------------------------------------
# 02_build_paired_set.R
# Keep only patients with BOTH a primary tumor (01) and a solid
# tissue normal (11) sample, one aliquot each, and build the sample
# sheet used by DESeq2 (design ~ patient + condition).
#
# TCGA barcode layout, e.g. TCGA-2A-A8VL-01A-21R-A37L-07
#   1-12  patient     (TCGA-2A-A8VL)
#   14-15 sample code (01 tumor, 11 normal)
#   16    vial        (A, B, ...)
#   20    analyte     (R, H, T ...)
#   22-25 plate       (A37L)
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(yaml)
})

if (exists("snakemake")) {
  cfg       <- snakemake@config
  in_se     <- snakemake@input[["se_full"]]
  out_se    <- snakemake@output[["se_paired"]]
  out_sheet <- snakemake@output[["sample_sheet"]]
  con <- file(snakemake@log[[1]], open = "wt"); sink(con); sink(con, type = "message")
} else {
  args <- commandArgs(trailingOnly = TRUE)
  cfg  <- yaml::read_yaml(if (length(args)) args[1] else "config/config.yaml")
  in_se     <- cfg$outputs$se_full
  out_se    <- cfg$outputs$se_paired
  out_sheet <- cfg$outputs$sample_sheet
}

se <- readRDS(in_se)
bc <- colnames(se)

meta <- data.frame(
  barcode     = bc,
  patient     = substr(bc, 1, 12),
  sample_code = substr(bc, 14, 15),
  vial        = substr(bc, 16, 16),
  analyte     = substr(bc, 20, 20),
  plate       = substr(bc, 22, 25),
  tss         = substr(bc, 6, 7),       # tissue source site (batch proxy)
  stringsAsFactors = FALSE
)
meta <- meta[meta$sample_code %in% c("01", "11"), ]
meta$condition <- ifelse(meta$sample_code == "01", "tumor", "normal")

# ---- deduplicate: one aliquot per patient x condition -------------
prio <- unlist(cfg$analyte_priority)
meta$analyte_rank <- match(meta$analyte, prio)
meta$analyte_rank[is.na(meta$analyte_rank)] <- length(prio) + 1
meta <- meta[order(meta$patient, meta$condition, meta$vial,
                   meta$analyte_rank, -rank(meta$plate)), ]
n_before <- nrow(meta)
meta <- meta[!duplicated(meta[, c("patient", "condition")]), ]
message("Duplicate aliquots dropped: ", n_before - nrow(meta))

# ---- keep patients with both tumor and normal ---------------------
tab <- table(meta$patient, meta$condition)
paired_patients <- rownames(tab)[tab[, "tumor"] == 1 & tab[, "normal"] == 1]
meta <- meta[meta$patient %in% paired_patients, ]
message("Paired patients: ", length(paired_patients),
        "  (samples: ", nrow(meta), ")")
stopifnot(length(paired_patients) >= 10)

# ---- carry over useful clinical columns if present ----------------
cd <- as.data.frame(colData(se))
keep_clin <- intersect(
  c("race", "ethnicity", "age_at_diagnosis", "primary_gleason_grade",
    "secondary_gleason_grade", "ajcc_pathologic_t", "paper_Reviewed_Gleason",
    "vital_status"),
  colnames(cd)
)
meta <- cbind(meta, cd[meta$barcode, keep_clin, drop = FALSE])

meta$condition <- factor(meta$condition, levels = c("normal", "tumor"))
meta$patient   <- factor(meta$patient)
meta$analyte_rank <- NULL

se_paired <- se[, meta$barcode]
# GDCprepare's colData already has a 'patient' column; overwrite
# explicitly rather than cbind, to avoid duplicate column names.
se_paired$patient   <- meta$patient
se_paired$condition <- meta$condition
se_paired$tss       <- meta$tss
se_paired$plate     <- meta$plate

saveRDS(se_paired, out_se)
write.table(meta, out_sheet, sep = "\t", quote = FALSE, row.names = FALSE)

message("Race breakdown among paired patients:")
if ("race" %in% colnames(meta)) print(table(meta$race[meta$condition == "tumor"], useNA = "ifany"))
message("Done.")
