#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
script <- sub("^--file=", "", grep("^--file=", args, value = TRUE)[1])
root <- normalizePath(dirname(script), mustWork = TRUE)
Sys.setenv(RADNET_CODE_PACKAGE = root)
source(file.path(root, "R", "common.R"))

steps <- c("01_validate_source_data.R", "02_reproduce_figures_tables.R", "03_refit_models.R", "05_refit_figure6.R", "06_finalize_v3_displays.R", "04_session_info.R")
for (step in steps) {
  cat("\n=== ", step, " ===\n", sep = "")
  source(file.path(root, "R", step), chdir = FALSE)
}
cat("\nAll packaged reproduction steps completed successfully.\n")
