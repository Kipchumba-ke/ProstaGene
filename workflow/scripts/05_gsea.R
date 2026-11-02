#!/usr/bin/env Rscript
# ---------------------------------------------------------------
# 05_gsea.R — pathway analysis of the paired DESeq2 results
#
#  1. Rank all tested genes by DESeq2 Wald statistic
#     (one entry per gene symbol: keep the max |stat| Ensembl ID)
#  2. fgsea on MSigDB Hallmark and Reactome (via msigdbr)
#  3. Collapse redundant Reactome pathways (fgsea::collapsePathways)
#  4. Cross-check: over-representation analysis (clusterProfiler)
#     of up- and down-regulated DEGs against Hallmark
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(ggplot2)
})

if (exists("snakemake")) {
  inp <- snakemake@input; out <- snakemake@output; p <- snakemake@params
  threads <- snakemake@threads
  con <- file(snakemake@log[[1]], open = "wt"); sink(con); sink(con, type = "message")
} else {
  inp <- list(de = "results/de/deseq2_all_genes.tsv")
  g <- "results/gsea"
  out <- list(ranks = file.path(g, "ranked_genes.rnk"),
              hallmark = file.path(g, "fgsea_hallmark.tsv"),
              reactome = file.path(g, "fgsea_reactome.tsv"),
              ora = file.path(g, "ora_hallmark.tsv"),
              concordance = file.path(g, "hallmark_gsea_vs_ora.tsv"),
              hallmark_plot = file.path(g, "hallmark_nes.png"),
              reactome_plot = file.path(g, "reactome_top_dotplot.png"),
              enrichment_pdf = file.path(g, "hallmark_enrichment_plots.pdf"))
  p <- list(fdr = 0.05, min_size = 15, max_size = 500, seed = 42)
  threads <- 4
}
for (f in unlist(out)) dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
theme_set(theme_bw(base_size = 11))
fdr <- p[["fdr"]]; min_size <- p[["min_size"]]; max_size <- p[["max_size"]]
set.seed(p[["seed"]])

# ---- 1. ranked gene list -----------------------------------------
de <- fread(inp[["de"]])
de <- de[!is.na(stat) & !is.na(gene_name) & gene_name != ""]
de <- de[order(-abs(stat))][!duplicated(gene_name)]          # one row per symbol
ranks <- setNames(de$stat, de$gene_name)
ranks <- sort(ranks, decreasing = TRUE)
message("Ranked genes (unique symbols): ", length(ranks))
fwrite(data.table(gene = names(ranks), stat = ranks), out[["ranks"]],
       sep = "\t", col.names = FALSE)                          # GSEA-desktop .rnk format

# ---- 2. gene sets ------------------------------------------------
# msigdbr >= 10 renamed category/subcategory -> collection/subcollection
# and moved non-Hallmark collections to the 'msigdbdf' package.
get_sets <- function(coll, sub = NULL) {
  new_api <- "collection" %in% names(formals(msigdbr::msigdbr))
  df <- tryCatch(
    if (new_api) msigdbr(species = "Homo sapiens", collection = coll, subcollection = sub)
    else         msigdbr(species = "Homo sapiens", category   = coll, subcategory   = sub),
    error = function(e) stop(
      "msigdbr could not load ", coll, " ", sub, ": ", conditionMessage(e),
      "\nIf this mentions 'msigdbdf', install it with:\n",
      "  install.packages('msigdbdf', repos = c('https://igordot.r-universe.dev', ",
      "'https://cloud.r-project.org'))", call. = FALSE))
  df <- unique(as.data.table(df)[, .(gs_name, gene_symbol)])
  list(list = split(df$gene_symbol, df$gs_name), t2g = df)
}
hall <- get_sets("H")
reac <- get_sets("C2", "CP:REACTOME")
message("Gene sets: Hallmark ", length(hall$list), ", Reactome ", length(reac$list))

run_fgsea <- function(sets) {
  r <- fgsea(pathways = sets, stats = ranks, minSize = min_size, maxSize = max_size,
             eps = 0, nproc = threads)
  r[order(padj, -abs(NES))]
}
flatten_le <- function(r) {
  r <- copy(r)
  r[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
  r
}

# ---- 3a. Hallmark ------------------------------------------------
fh <- run_fgsea(hall$list)
message("Hallmark sets tested: ", nrow(fh), " | padj < ", fdr, ": ", sum(fh$padj < fdr),
        " [up ", sum(fh$padj < fdr & fh$NES > 0), ", down ", sum(fh$padj < fdr & fh$NES < 0), "]")
fwrite(flatten_le(fh), out[["hallmark"]], sep = "\t")
message("Hallmark, significant (by NES):")
print(fh[padj < fdr][order(-NES), .(pathway = sub("^HALLMARK_", "", pathway),
                                    NES = round(NES, 2), padj = signif(padj, 3), size)])

# ---- 3b. Reactome + redundancy collapse --------------------------
fr <- run_fgsea(reac$list)
fr[, main_pathway := NA]
sig_r <- fr[padj < fdr]
if (nrow(sig_r) > 1) {
  coll <- tryCatch(
    collapsePathways(sig_r[order(pval)], reac$list, ranks),
    error = function(e) { message("collapsePathways failed: ", conditionMessage(e)); NULL })
  if (!is.null(coll)) fr[padj < fdr, main_pathway := pathway %in% coll$mainPathways]
}
message("Reactome sets tested: ", nrow(fr), " | padj < ", fdr, ": ", nrow(sig_r),
        " | non-redundant: ", sum(fr$main_pathway %in% TRUE))
fwrite(flatten_le(fr), out[["reactome"]], sep = "\t")

# ---- plots -------------------------------------------------------
pretty <- function(x, prefix) {
  x <- gsub("_", " ", sub(prefix, "", x))
  ifelse(nchar(x) > 60, paste0(substr(x, 1, 57), "..."), x)
}

hp <- fh[padj < fdr]
if (nrow(hp)) {
  hp[, label := pretty(pathway, "^HALLMARK_")]
  hp[, label := factor(label, levels = label[order(NES)])]
  g <- ggplot(hp, aes(NES, label, fill = NES > 0)) +
    geom_col() +
    geom_text(aes(label = formatC(padj, format = "e", digits = 1),
                  hjust = ifelse(NES > 0, -0.1, 1.1)), size = 2.8) +
    scale_fill_manual(values = c(`TRUE` = "#C8553D", `FALSE` = "#2E6F9E"),
                      labels = c(`TRUE` = "up in tumor", `FALSE` = "down in tumor"), name = NULL) +
    scale_x_continuous(expand = expansion(mult = 0.25)) +
    labs(x = "Normalised enrichment score", y = NULL,
         title = "Hallmark GSEA: tumor vs matched normal",
         subtitle = paste0("fgsea on DESeq2 Wald statistic; padj < ", fdr, " shown (labels = padj)"))
  ggsave(out[["hallmark_plot"]], g, width = 9, height = max(4, 0.28 * nrow(hp) + 1.5), dpi = 150)
} else file.create(out[["hallmark_plot"]])

rp <- fr[padj < fdr & (main_pathway %in% TRUE | all(is.na(main_pathway)))]
rp <- rbind(head(rp[NES > 0][order(-NES)], 15), head(rp[NES < 0][order(NES)], 15))
if (nrow(rp)) {
  rp[, label := pretty(pathway, "^REACTOME_")]
  rp[, label := factor(label, levels = unique(label[order(NES)]))]
  g <- ggplot(rp, aes(NES, label, size = size, colour = -log10(padj))) +
    geom_vline(xintercept = 0, colour = "grey60") +
    geom_point() +
    scale_colour_viridis_c(option = "plasma", end = 0.9) +
    labs(x = "Normalised enrichment score", y = NULL, size = "Set size",
         colour = expression(-log[10]~padj),
         title = "Reactome GSEA: top non-redundant pathways",
         subtitle = "Up to 15 up and 15 down; redundancy collapsed with fgsea::collapsePathways")
  ggsave(out[["reactome_plot"]], g, width = 10, height = max(4, 0.3 * nrow(rp) + 1.5), dpi = 150)
} else file.create(out[["reactome_plot"]])

# Running-score plots: top 3 up and top 3 down Hallmark sets
pick <- c(head(fh[padj < fdr & NES > 0][order(-NES)]$pathway, 3),
          head(fh[padj < fdr & NES < 0][order(NES)]$pathway, 3))
pdf(out[["enrichment_pdf"]], width = 6, height = 4)
for (pw in pick) {
  r <- fh[pathway == pw]
  print(plotEnrichment(hall$list[[pw]], ranks) +
          labs(title = pretty(pw, "^HALLMARK_"),
               subtitle = sprintf("NES = %.2f, padj = %.1e", r$NES, r$padj),
               x = "Gene rank (tumor-up -> tumor-down)", y = "Enrichment score"))
}
dev.off()

# ---- 4. ORA cross-check (Hallmark) -------------------------------
universe <- names(ranks)
up   <- intersect(de[direction == "up"]$gene_name,   universe)
down <- intersect(de[direction == "down"]$gene_name, universe)
t2g  <- as.data.frame(hall$t2g)
ora1 <- function(genes, label) {
  e <- enricher(genes, universe = universe, TERM2GENE = t2g,
                minGSSize = min_size, maxGSSize = max_size,
                pvalueCutoff = 1, qvalueCutoff = 1)
  if (is.null(e)) return(NULL)
  d <- as.data.table(as.data.frame(e)); d[, direction := label]; d
}
ora <- rbind(ora1(up, "up"), ora1(down, "down"), fill = TRUE)
if (is.null(ora)) ora <- data.table()
fwrite(ora, out[["ora"]], sep = "\t")

ora_up   <- if (nrow(ora)) ora[direction == "up",   .(pathway = ID, padj_ora_up = p.adjust)]   else data.table(pathway = character(), padj_ora_up = numeric())
ora_down <- if (nrow(ora)) ora[direction == "down", .(pathway = ID, padj_ora_down = p.adjust)] else data.table(pathway = character(), padj_ora_down = numeric())
conc <- merge(fh[, .(pathway, NES, padj_gsea = padj)], ora_up, by = "pathway", all.x = TRUE)
conc <- merge(conc, ora_down, by = "pathway", all.x = TRUE)
conc[, agrees := NA]
conc[padj_gsea < fdr & NES > 0, agrees := !is.na(padj_ora_up)   & padj_ora_up   < fdr]
conc[padj_gsea < fdr & NES < 0, agrees := !is.na(padj_ora_down) & padj_ora_down < fdr]
conc <- conc[order(padj_gsea)]
fwrite(conc, out[["concordance"]], sep = "\t")
n_sig <- sum(!is.na(conc$agrees))
message("ORA cross-check: of ", n_sig, " significant Hallmark GSEA sets, ",
        sum(conc$agrees, na.rm = TRUE), " are also significant by ORA in the same direction")

writeLines(capture.output(sessionInfo()), sub("\\.tsv$", "_sessionInfo.txt", out[["hallmark"]]))
message("Done.")
