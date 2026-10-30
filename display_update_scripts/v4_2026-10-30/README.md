# Figure revision 2026-10-30 redraw scripts

Run from the repository root:

```sh
Rscript display_update_scripts/v4_2026-10-30/run_all.R
```

This redraws Figures 1–6 and Supplementary Figures S1–S8 without fitting models.
It writes PNG, PDF, TIFF and SVG files to
`data/display_assets/v4_2026-10-30/figures/`, using the manuscript filenames in
`figure_names.json`. Existing v4 artwork is overwritten when this command runs.
The supplied revision artwork is committed unchanged; redraw file hashes can
vary with fonts, graphics libraries and PDF metadata.

Requirements: R packages `data.table`, `ggplot2` (with inside-legend support),
`patchwork`, `jsonlite`, `svglite`, `ragg`, `maps`, and `sf`; Cairo PDF support
and Helvetica or a compatible font. Python 3 requires PyMuPDF (`fitz`) and
Pillow. Set `QA_PYTHON` to that Python interpreter if `python3` is unsuitable.
No software is installed and no upstream data are downloaded by these scripts.

Optional selections, also from the repository root:

```sh
Rscript display_update_scripts/v4_2026-10-30/run_all.R Figure_1 S1
Rscript display_update_scripts/v4_2026-10-30/run_all.R quantitative
Rscript display_update_scripts/v4_2026-10-30/run_all.R Figure_2 Figure_6 S2 S6 S7
```

`quantitative` redraws Figures 3–5 and S3–S5/S8 together. S1 reads the supplied
Figure 1 combustion-role table, so it can also run alone. All paths are resolved
from the repository root; v5 sources and frozen_v3 artwork remain read-only.
Generated QA, previews, logs, caches and sizing files stay under
`data/display_assets/v4_2026-10-30/qa/` and its sibling `.tmp/`.

| Script | Produces / purpose |
| --- | --- |
| `run_all.R` | Initializes output directories and runs the selected drawing modules, including Figure 2 extraction; records the R session |
| `Figure1_display_v3.R` | Figure 1 combustion-only map and the three new Figure 1 location/count files |
| `revise_Figure1.R` | Figure 1-only wrapper around `run_all.R` |
| `extract_Figure2.py` | Figure 2 in four formats, cropping the title band from frozen_v3; writes crop validation |
| `02_reproduce_figures_tables.R` | Figures 4/5 and S8; orchestrates Figure 3 and S3–S5 (no table or model regeneration) |
| `02a_redraw_overall_inventory.R` | Figure 3, with its wind legend inside panel a |
| `02b_redraw_supplementary_diagnostics.R` | S3/S4 monitor-FE estimates and S5 leave-one-monitor-out results |
| `Figure6_display_v3.R` | Figure 6 from unchanged v5 estimates |
| `Supplementary_Figure_S1.R` | Six non-exclusive combustion-role panels and three S1 location/count files |
| `Supplementary_Figure_S2.R` | Three regional maps and a regional-count recheck |
| `Supplementary_Figure_S6.R` | S6 non-combustion-adjusted estimates |
| `Supplementary_Figure_S7.R` | S7 current versus spline estimates |
| `export_helpers.R` | Common four-format exports, exact plotted-layer comparisons (tolerance 0) and rendered alignment measurements |
| `preflight.py` | Static module checks with shared export context |
| `compare_csv_sources.R` | Exact v5 versus archived reproduction-source checks, plus Figure 1/S1 subset and count verification |
| `make_comparisons.R` | Before/after PNGs for all 14 figures (requires R `magick`) |
| `compare_Supplementary_Figure_S1.R` | S1-only before/after comparison (requires R `magick`) |
| `run_audits.py` | PDF/raster export checks and collision/font/alignment diagnostics |
| `qa_helpers/` | Bundled alignment, collision, effective-font, text and static figure auditors |
| `run_all.sh` | Optional zsh wrapper running preflight, all redraws, comparisons and audits |

For the optional full wrapper:

```sh
zsh display_update_scripts/v4_2026-10-30/run_all.sh
```

`run_audits.py` retains raw font diagnostics; use
`qa_helpers/effective_font_audit.py` for Cairo text-matrix-adjusted font sizes.
The wrapper does not regenerate the supplied historical QA report.

All direct numerical inputs come from `source_data_manuscript_v5/`; Figure 2
uses `data/display_assets/frozen_v3/figures/`. Figure/source mappings are listed
in the main README. The optional CSV audit additionally reads the existing
`data/source_data/` reproduction exports for independent comparisons, without
copying them. S5 ordered indices are reconstructed by sorting unchanged v5
estimates within each fuel/metric; its full-sample reference comes from the
20-km rows of `Figure_4_estimates.csv`. S2 uses unique plant IDs/coordinates from
`Figure_1A_powerplant_locations.csv`; these exactly match the original physical
plant coordinate set.

The input delivery's `source_data/original/`, `source_data/reproduction/` and
`scripts/original/` staging copies are omitted. Its one-time staging/code-generation
tools (`stage_inputs.py`, `prepare_scripts.py`) and report/follow-up mutation
tools (`write_report.py`, `qa_Supplementary_Figure_S1.py`) are also omitted:
they depend on the external delivery workspace and its unbundled historical
QA artifacts, and are not needed to redraw the figures. Existing historical
scripts remain under `legacy_scripts/` and the original repository code paths.
