#!/usr/bin/env Rscript
# Run from the repository root: Rscript display_update_scripts/v4_2026-10-30/run_all.R
ROOT <- normalizePath(getwd())
stopifnot(file.exists('display_update_scripts/v4_2026-10-30/export_helpers.R'))
for(d in c('data/display_assets/v4_2026-10-30/figures','data/display_assets/v4_2026-10-30/qa/audits','data/display_assets/v4_2026-10-30/qa/previews','data/display_assets/v4_2026-10-30/qa/logs','data/display_assets/v4_2026-10-30/source_data','data/display_assets/v4_2026-10-30/.tmp','data/display_assets/v4_2026-10-30/qa/cache'))dir.create(d,recursive=TRUE,showWarnings=FALSE)
Sys.setenv(XDG_CACHE_HOME=file.path(ROOT,'data/display_assets/v4_2026-10-30/qa/cache'),R_USER_CACHE_DIR=file.path(ROOT,'data/display_assets/v4_2026-10-30/qa/cache'),
 NATURE_FIGURE_CONFIG=file.path(ROOT,'data/display_assets/v4_2026-10-30/qa/nature-figure-backend.json'))
source('display_update_scripts/v4_2026-10-30/export_helpers.R')
SOURCE_DIR <- file.path(ROOT,'source_data_manuscript_v5')
OUTPUT_DIR <- ROOT
combustion <- c('Coal','Natural gas','Oil','Other fossil','Biomass','Solid waste')
args <- commandArgs(trailingOnly=TRUE)
selected <- function(x)length(args)==0 || x %in% args
if(selected('Figure_1'))source('display_update_scripts/v4_2026-10-30/Figure1_display_v3.R')
if(selected('quantitative'))source('display_update_scripts/v4_2026-10-30/02_reproduce_figures_tables.R')
if(selected('Figure_6')){
 source('display_update_scripts/v4_2026-10-30/Figure6_display_v3.R')
 render_figure6_v3(fread('source_data_manuscript_v5/Figure_6_internal_validity_robustness.csv'),
 'data/display_assets/v4_2026-10-30/figures/Figure_6','data/display_assets/v4_2026-10-30/qa/audits','data/display_assets/v4_2026-10-30/qa/audits',skill_root='display_update_scripts/v4_2026-10-30')
}
for(id in c('S1','S2','S6','S7'))if(selected(id))source(paste0('display_update_scripts/v4_2026-10-30/Supplementary_Figure_',id,'.R'))
if(selected('Figure_2')) {
 py <- Sys.getenv('QA_PYTHON',unset='python3')
 status <- system2(py,shQuote('display_update_scripts/v4_2026-10-30/extract_Figure2.py'))
 if(status!=0)stop('Figure 2 extraction failed')
}
writeLines(capture.output(sessionInfo()),'data/display_assets/v4_2026-10-30/qa/R_session_info.txt')
cat('R exports complete.\n')
