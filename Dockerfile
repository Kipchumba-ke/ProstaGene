# ProstaGene pipeline image: R/Bioconductor + Snakemake + Quarto from workflow/envs/r-bioc.yaml
#
#   docker build -t prostagene .
#   docker run --rm -v "$PWD/data:/work/data" -v "$PWD/results:/work/results" prostagene snakemake --cores 4
#
ARG MAMBA_TAG=1.5.10
FROM mambaorg/micromamba:${MAMBA_TAG}

LABEL org.opencontainers.image.title="ProstaGene" \
      org.opencontainers.image.description="TCGA-PRAD paired differential expression and pathway pipeline" \
      org.opencontainers.image.licenses="MIT"

# 1. Environment (cached unless the env file changes)
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/r-bioc.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && micromamba clean --all --yes

# Activate the env for the following RUN steps and at run time
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# 2. Packages not reliably on conda channels (tidyestimate from CRAN, MCPcounter from GitHub)
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/scripts/install_extras.R /tmp/install_extras.R
RUN Rscript /tmp/install_extras.R

# 3. Project code (the working directory must be writable by the non-root image user)
USER root
RUN mkdir -p /work && chown "$MAMBA_USER:$MAMBA_USER" /work
USER $MAMBA_USER
WORKDIR /work
COPY --chown=$MAMBA_USER:$MAMBA_USER . /work

CMD ["snakemake", "--cores", "4"]
