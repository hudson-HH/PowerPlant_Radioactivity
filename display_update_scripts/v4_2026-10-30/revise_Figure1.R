#!/usr/bin/env Rscript
# Optional Figure 1-only wrapper; run from the repository root.
status <- system2(file.path(R.home('bin'),'Rscript'),
 c('display_update_scripts/v4_2026-10-30/run_all.R','Figure_1'))
if(status!=0)stop('Figure 1 redraw failed')
