#!/usr/bin/env Rscript
# ---------------------------------------------------------------
# 01_download_tcga.R
# Query GDC for TCGA-PRAD STAR counts (Primary Tumor + Solid Tissue
# Normal), download, and assemble a RangedSummarizedExperiment.
#
# Run via Snakemake (preferred) or standalone:
#   Rscript workflow/scripts/01_download_tcga.R config/config.yaml
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(yaml)
})

# ---- config: Snakemake object if present, else CLI arg ----------
if (exists("snakemake")) {
  cfg      <- snakemake@config
  out_se   <- snakemake@output[["se_full"]]
  out_man  <- snakemake@output[["manifest"]]
  out_prov <- snakemake@output[["provenance"]]
  log_file <- snakemake@log[[1]]
  con <- file(log_file, open = "wt"); sink(con); sink(con, type = "message")
} else {
  args <- commandArgs(trailingOnly = TRUE)
  cfg  <- yaml::read_yaml(if (length(args)) args[1] else "config/config.yaml")
  out_se   <- cfg$outputs$se_full
  out_man  <- cfg$outputs$manifest
  out_prov <- cfg$outputs$provenance
}
dir.create(dirname(out_se), recursive = TRUE, showWarnings = FALSE)

# ---- record which GDC data release we are pulling --------------
gdc_info <- tryCatch(getGDCInfo(), error = function(e) NULL)
message("GDC data release: ", if (is.null(gdc_info)) "unknown" else gdc_info$data_release)

# ---- query -------------------------------------------------------
query <- GDCquery(
  project       = cfg$project,
  data.category = cfg$gdc$data_category,
  data.type     = cfg$gdc$data_type,
  workflow.type = cfg$gdc$workflow_type,
  sample.type   = unlist(cfg$gdc$sample_types)
)

res <- getResults(query)
message("Files returned by query: ", nrow(res))
print(table(res$sample_type))

# File-level manifest (file IDs, barcodes) for reproducibility
manifest_cols <- intersect(
  c("id", "file_name", "cases", "sample.submitter_id", "sample_type",
    "md5sum", "file_size", "updated_datetime"),
  colnames(res)
)
write.table(res[, manifest_cols], out_man, sep = "\t",
            quote = FALSE, row.names = FALSE)

# ---- download (resumable: already-present files are skipped) -----
GDCdownload(
  query,
  method          = "api",
  directory       = cfg$gdc$download_dir,
  files.per.chunk = cfg$gdc$files_per_chunk
)

# ---- assemble SummarizedExperiment -------------------------------
# GDCprepare attaches clinical + sample metadata (incl. race/ethnicity,
# Gleason, stage) to colData, and gene_name/gene_type to rowData.
se <- GDCprepare(query, directory = cfg$gdc$download_dir, save = FALSE)

message("Assays available: ", paste(assayNames(se), collapse = ", "))
stopifnot(cfg$count_assay %in% assayNames(se))
message("Dimensions (genes x samples): ", paste(dim(se), collapse = " x "))

saveRDS(se, out_se)

# ---- provenance --------------------------------------------------
prov <- list(
  date_downloaded   = as.character(Sys.time()),
  gdc_data_release  = if (is.null(gdc_info)) NA else gdc_info$data_release,
  project           = cfg$project,
  workflow_type     = cfg$gdc$workflow_type,
  n_files           = nrow(res),
  n_by_sample_type  = as.list(table(res$sample_type)),
  tcgabiolinks      = as.character(packageVersion("TCGAbiolinks")),
  r_version         = R.version.string
)
yaml::write_yaml(prov, out_prov)
message("Done.")
