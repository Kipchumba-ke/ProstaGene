# Differential Expression & Pathway Analysis in Prostate Cancer (TCGA-PRAD)

Paired tumor vs. matched normal RNA-seq analysis of TCGA prostate adenocarcinoma, with
pathway enrichment and validation in an independent GEO cohort. Built as a fully
reproducible Snakemake pipeline with a Quarto report.

**Author:** Kipchumba · Nursing & Public Health · Bioinformatics

## Questions
1. Which genes are differentially expressed between prostate tumors and matched normal tissue, controlling for patient?
2. Which Hallmark and Reactome pathways are enriched among those changes?
3. Do the findings replicate in an independent cohort?
4. *(Exploratory)* What can and can't be said about ancestry differences, given TCGA's limited African-ancestry representation?

## Approach
| Step | Tooling |
|---|---|
| Data | GDC harmonized STAR counts via `TCGAbiolinks` (Primary Tumor + Solid Tissue Normal) |
| Pairing | Patients with both sample types; one aliquot per sample type |
| DGE | `DESeq2`, design `~ patient + condition`, `apeglm` LFC shrinkage |
| Pathways | `fgsea` / `clusterProfiler` on MSigDB Hallmark + Reactome (`msigdbr`) |
| Validation | Independent paired prostate RNA-seq cohort from GEO |
| Reproducibility | Snakemake + conda, Quarto report |

## Repository layout
```
config/config.yaml          # all parameters
workflow/Snakefile          # pipeline
workflow/envs/r-bioc.yaml   # R/Bioconductor environment
workflow/scripts/           # one R script per step
data/raw/                   # GDC downloads (git-ignored)
data/processed/             # SummarizedExperiments, sample sheet, provenance
results/                    # tables and figures
reports/                    # Quarto report
```

## Running it
Requires Linux or macOS (on Windows, use WSL2). About 2–3 GB of disk for the raw counts.

```bash
conda env create -f environment.yml
conda activate prad-pipeline
snakemake -n                          # dry run: shows planned steps
snakemake --use-conda --cores 4       # full run
```

**Using micromamba instead of conda:** Snakemake's `--use-conda` needs a real
`conda` (or `mamba`) binary, and shell aliases are not visible to it. Build one
combined environment and run without `--use-conda`:

```bash
micromamba create -n prad-dge -f workflow/envs/r-bioc.yaml   # pins python=3.12
micromamba install -n prad-dge -c conda-forge -c bioconda "python=3.12" "snakemake>=8" quarto
micromamba activate prad-dge
snakemake -n
snakemake --cores 4
```

Steps can also be run standalone: `Rscript workflow/scripts/01_download_tcga.R config/config.yaml`.

The GDC data release used is recorded in `data/processed/download_provenance.yaml`.

## Status
- [x] M1 Project structure + data download — 553 samples; 52 tumor/normal pairs (104 samples)
- [x] M2 QC & exploratory analysis
- [ ] M3 Paired DESeq2
- [ ] M4 GSEA (Hallmark, Reactome)
- [ ] M5 GEO validation
- [ ] M6 Ancestry-stratified comparison (optional)
- [ ] M7 Full pipeline + Quarto report

## Data use
TCGA open-access data are from the NCI Genomic Data Commons and are subject to
the [GDC Data Use Agreement](https://gdc.cancer.gov/access-data/data-access-processes-and-tools).
