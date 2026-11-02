# ProstaGene: Differential Expression & Pathway Analysis in Prostate Cancer

Paired tumor vs. matched normal RNA-seq analysis of TCGA prostate adenocarcinoma (TCGA-PRAD).
It includes pathway enrichment, validation in an independent GEO cohort, and an exploratory
ancestry analysis. Everything runs as one reproducible Snakemake pipeline and renders to a Quarto report.

**📄 [Read the full report](https://Kipchumba-ke.github.io/ProstaGene/docs/index.html)**

**Author:** Kipchumba, Nursing & Public Health · Bioinformatics

---

## Key findings

| | Result |
|---|---|
| **Design** | 52 TCGA-PRAD patients with tumor + matched normal (104 samples); DESeq2 `~ patient + condition` |
| **Differential expression** | 2,233 DEGs (FDR < 0.05, ≥ 2-fold): 763 up, 1,470 down in tumor. Top genes: DLX1, AMACR, PCA3, HPN, TDRD1 |
| **Pathways** | 27 / 50 Hallmark sets enriched. MYC targets, E2F/G2M proliferation and androgen response up; myogenesis, EMT and inflammatory sets down, largely reflecting stromal loss |
| **Robustness** | Excluding 4 borderline-purity pairs: fold-change ρ = 0.99, DEG Jaccard = 0.87 |
| **Validation (GSE104131, 13 pairs)** | 92% of TCGA DEGs same direction; TCGA signature AUC = 0.98; 20 / 27 Hallmark sets replicate (NES r = 0.73) |
| **Ancestry (exploratory)** | Tumors from Black patients (n = 57) show higher interferon-α/γ signalling than tumors from White patients (n = 384), robust to tissue-source-site adjustment. No ancestry-specific tumor effect detectable in the small GEO cohort |

<p align="center">
  <img src="results/gsea/hallmark_nes.png" width="48%">
  <img src="results/validation/hallmark_nes_concordance.png" width="40%">
</p>

## Methodological notes worth reading

- **Purity:** four tumors sat near the normal cluster on PCA. They were kept, with a sensitivity analysis, rather than dropped based on how the data looked.
- **A mapping artifact:** the strongest "race-associated" genes were mitochondrial pseudogenes, 14 of 15 higher in Black patients. The likely cause is African mtDNA haplogroups diverging from the European reference mitochondrial genome. They were excluded and reported separately.
- **An anti-conservative test caught and replaced:** pre-ranked GSEA on an ancestry interaction statistic flagged 41 / 50 pathways from just 7 vs 6 patients. A patient-level test (exact Wilcoxon on per-pair pathway scores) found none.
- **Confounding:** in TCGA, race is confounded with tissue source site and Gleason grade. Site-matched models are reported alongside the main ones.

## Pipeline

```
download_tcga → build_paired_set → qc → deseq2 → gsea → geo_validation → ancestry → report
```

| Step | Script | Tools |
|---|---|---|
| Download | `01_download_tcga.R` | TCGAbiolinks (GDC STAR counts) |
| Pairing | `02_build_paired_set.R` | one aliquot per patient × sample type |
| QC | `03_qc.R` | VST, PCA, sample distances, purity check |
| DGE | `04_deseq2.R` | DESeq2, apeglm, sensitivity analysis |
| Pathways | `05_gsea.R` | fgsea + msigdbr (Hallmark, Reactome), clusterProfiler ORA |
| Validation | `06_geo_validation.R` | GSE104131, NCBI counts, org.Hs.eg.db |
| Ancestry | `07_ancestry.R` | nested DESeq2 design, patient-level Wilcoxon, site-matched models |
| Report | `reports/report.qmd` | Quarto |

## Running it

Requires Linux or macOS (on Windows, use WSL2) and about 3 GB of disk.

```bash
micromamba create -n prad-dge -f workflow/envs/r-bioc.yaml
micromamba activate prad-dge
snakemake -n            # dry run
snakemake --cores 4     # full run
```

With conda instead of micromamba, `conda env create -f workflow/envs/r-bioc.yaml` works the same way.

**GEO files (manual download).** NCBI's FTP server can refuse automated downloads, so the two GSE104131 files are fetched by hand into `data/raw/geo/`:

```bash
mkdir -p data/raw/geo && cd data/raw/geo
curl -fL -A "Mozilla/5.0" -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104131/matrix/GSE104131_series_matrix.txt.gz
curl -fL -o GSE104131_raw_counts_GRCh38.p13_NCBI.tsv.gz \
  "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE104131&format=file&file=GSE104131_raw_counts_GRCh38.p13_NCBI.tsv.gz"
```

If the FTP link returns 403, download the series matrix from the GEO page for GSE104131 in a browser.

The GDC data release used is recorded in `data/processed/download_provenance.yaml`.

## Repository layout

```
config/config.yaml              # parameters
workflow/Snakefile              # pipeline
workflow/envs/r-bioc.yaml       # single environment for the whole project
workflow/scripts/01–07_*.R      # one script per step
data/processed/                 # sample sheets, GDC manifest, provenance (.rds git-ignored)
results/{qc,de,gsea,validation,ancestry}/   # tables and figures
reports/report.qmd              # report source
docs/index.html                 # rendered report (GitHub Pages)
```

## Limitations

Bulk tissue mixes cell-state and composition effects. Purity estimates are incomplete. The validation cohort is small and high-grade only. Self-reported race is not genetic ancestry, and African-ancestry samples are few. All ancestry results are exploratory. See the report for details.

## Data use

TCGA open-access data come from the NCI Genomic Data Commons and are subject to the
[GDC Data Use Agreement](https://gdc.cancer.gov/access-data/data-access-processes-and-tools).
GSE104131 is publicly available from NCBI GEO.

## License

MIT
