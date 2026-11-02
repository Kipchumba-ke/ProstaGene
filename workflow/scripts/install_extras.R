#!/usr/bin/env Rscript
# One-time install of packages not reliably available from conda channels.
# Run inside the activated environment:  Rscript workflow/scripts/install_extras.R

repos <- "https://cloud.r-project.org"

if (!requireNamespace("tidyestimate", quietly = TRUE))
  install.packages("tidyestimate", repos = repos)

if (!requireNamespace("MCPcounter", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = repos)
  remotes::install_github("ebecht/MCPcounter", ref = "master", subdir = "Source", upgrade = "never")
}

for (pkg in c("tidyestimate", "MCPcounter", "survival", "readxl"))
  cat(sprintf("%-14s %s\n", pkg, if (requireNamespace(pkg, quietly = TRUE)) "OK" else "MISSING"))
