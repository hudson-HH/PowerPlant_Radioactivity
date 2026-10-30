> Imported QA report supplied with the 2026-10-30 figure revision. Its visual
> review and numerical findings describe the supplied artwork. References below
> to audit JSON/RDS files and before/after previews identify historical delivery
> artifacts; those artifacts are not bundled in this repository. Current artwork:
> `data/display_assets/v4_2026-10-30/figures/`; new display tables:
> `data/display_assets/v4_2026-10-30/source_data/`. Portable redraw commands are in
> `display_update_scripts/v4_2026-10-30/README.md`. The original staging folders
> and one-time setup commands described at the end are historical, superseded
> for this repository by that README.

# Presentation-only figure revision QA

All 14 figures were changed. None needed no change. No model was refitted and no estimate, confidence interval, analytical sample, category colour, axis transformation, axis range or scientific reference line was changed.

## Delivered files

`figures/` contains 14 PNGs at 600 dpi, 14 vector PDFs with embedded fonts, 14 LZW TIFFs at 600 dpi, and 14 SVGs. `qa/` contains one side-by-side frozen/new PNG for every figure. Scientific plots contain no raster image objects in their PDFs. Figure 2 retains its original slide raster effects alongside the unchanged vectors and embedded text.

## Data and axis checks

- `audits/csv_source_comparison.json`: all shared fields in the reproduction display inputs match the manuscript v5 CSVs exactly; maximum absolute difference 0. This includes estimates, CI endpoints, keys, wind references and available sample-size columns.
- `audits/numerical_invariance.json` and `*_plot_signatures.rds`: exact before/after comparison of every built scientific layer column, facet layout, axis numeric range, breaks and tick labels. All 12 applicable R comparisons passed at tolerance 0. S2 excludes only the expressly removed label layer from this comparison.
- Figure 1 retains its permitted original role subset; S1 now exactly matches that Figure 1 role table, with unchanged coordinates. All retained Figure 1 points are inside the original geographic bounds.
- S2 independently reproduces the frozen regional monitor/plant counts, original geographic limits, projection, buffers and point data.
- Figure 2 retains exact cropped PNG pixels and unchanged box text; only the title band/object was removed.
- Decimal tick formatting was explicitly retained with the original `scipen=999` setting.

## Figure 1 combustion counts

Unique physical plants: **5,068**. Non-exclusive plant–role points: **8,326**. RadNet locations: **168**, unchanged.

| Role | Plants / role points |
|---|---:|
| Coal | 798 |
| Natural gas | 3,500 |
| Oil | 2,675 |
| Other fossil | 108 |
| Biomass | 1,057 |
| Solid waste | 188 |

Original role rows: 18,205; non-combustion role rows removed: 9,879. Plants with multiple combustion roles retain one plotted point per role.

## S1 non-exclusive panel counts

S1 now uses the exact Figure 1 combustion-role location table. All 8,326 distinct plant–role rows across 5,068 physical plants are included; no rows are excluded. A multi-role plant appears in every applicable panel. All six expected counts were independently verified against the input table and exported source data.

| Panel | Plants / role points |
|---|---:|
| Coal | 798 |
| Natural gas | 3,500 |
| Oil | 2,675 |
| Other fossil | 108 |
| Biomass | 1,057 |
| Solid waste | 188 |

Six populated panels, ordered Coal, Natural gas, Oil / Other fossil, Biomass, Solid waste. No RadNet points. The earlier exclusive single-role rule is superseded.

## Per-figure changes and verification

### Figure 1

Single combustion-only map; original non-exclusive roles, colours, shapes, alpha, marker sizes, draw order, geographic limits and all 168 RadNet triangles retained. Biomass remains magenta diamonds. Former panel b removed. 5,068 unique plants; 8,326 plant–role points.

Removed text:

- Power-plant inventory and RadNet locations
- a  Energy types (non-exclusive)
- b  Combustion-based plant classes

Export size: 7.20333 × 2.91333 inches. PNG/TIFF: 600 dpi; TIFF: LZW. All PDF fonts embedded; effective minimum glyph size 6.199 pt. Collision audit: PASS. Measured alignment audit: PASS (single-panel maps are not applicable).

Historical comparison artifact (not bundled): `Figure_1_before_after.png`

### Figure 2

Original title band cropped by 192 pixels at 600 dpi. Every retained PNG pixel and all box text are unchanged. Vector objects, embedded fonts and the original slide raster effects are preserved.

Removed text:

- Figure 2. Study design and analytical framework

Export size: 7.20500 × 3.73333 inches. PNG/TIFF: 600 dpi; TIFF: LZW. All PDF fonts embedded; effective minimum glyph size 8.639 pt. Collision audit: PASS. Alignment: not applicable to the original flowchart.

Historical comparison artifact (not bundled): `Figure_2_before_after.png`

### Figure 3

Original six-entry wind legend moved inside panel a at top right, clear of curves and grid lines. Panels b/c use the original plot objects. All estimates, CIs, reference lines, axis ranges and tick labels unchanged.

Removed text:

- Overall plant inventory and ambient gross-beta activity
- 140 monitors; 558,915 quasi-daily rows; directional comparisons evaluated at mean wind 2.986 m/s

Export size: 11.00000 × 5.80000 inches. PNG/TIFF: 600 dpi; TIFF: LZW. All PDF fonts embedded; effective minimum glyph size 8.000 pt. Collision audit: PASS. Measured alignment audit: PASS (single-panel maps are not applicable).

Historical comparison artifact (not bundled): `Figure_3_before_after.png`

### Figure 4

Only overall title/subtitle removed; all eight facets retained.

Removed text:

- Fuel-specific upwind actual-operation associations
- Monitor × year fixed effects; effects evaluated at each model's sample-median wind speed

Export size: 11.00000 × 6.70000 inches. PNG/TIFF: 600 dpi; TIFF: LZW. All PDF fonts embedded; effective minimum glyph size 8.000 pt. Collision audit: PASS. Measured alignment audit: PASS (single-panel maps are not applicable).

Historical comparison artifact (not bundled): `Figure_4_before_after.png`

### Figure 5

Only overall title/subtitle removed; all eight facets retained.

Removed text:

- Directional contrast: upwind minus downwind
- Joint directional models with monitor × year fixed effects; sample-median wind speed

Export size: 11.00000 × 6.70000 inches. PNG/TIFF: 600 dpi; TIFF: LZW. All PDF fonts embedded; effective minimum glyph size 8.000 pt. Collision audit: PASS. Measured alignment audit: PASS (single-panel maps are not applicable).

Historical comparison artifact (not bundled): `Figure_5_before_after.png`

### Figure 6

Only overall title/subtitle/caption removed; panel a/b/c titles, six data facets, model legends and all 56 estimates/CIs retained.

Removed text:

- Internal validity and robustness
- 20 km upwind exposure; shared baseline sample-median wind; 95% confidence intervals
- Count: 3.154 m/s (all fuels). Heat: coal 3.153, gas 3.139, oil 3.136, biomass 3.154 m/s.
- Reference medians are shared across a-c; model-specific sample sizes are provided in Source Data.

Export size: 7.20333 × 7.87333 inches. PNG/TIFF: 600 dpi; TIFF: LZW. All PDF fonts embedded; effective minimum glyph size 6.199 pt. Collision audit: PASS. Measured alignment audit: PASS (single-panel maps are not applicable).

Historical comparison artifact (not bundled): `Figure_6_before_after.png`

### Supplementary Figure S1

Six non-exclusive combustion-role panels in a 3-column × 2-row grid, using every row and unchanged coordinate from the Figure 1 role table. Multi-role plants appear in each applicable panel. Grey state outlines and black points (size 0.45, alpha 0.55) retained; no overall title or RadNet points.

Export size: **13 × 5.125 inches** (previously 13 × 7). The height is derived from the fixed-aspect map grid and title strips, with small uniform 5.76-pt outer margins; geographic bounds and projection are unchanged. PNG/TIFF: 7,800 × 3,075 pixels, 600 dpi; TIFF: LZW. PDF: vector layers and embedded fonts. Effective minimum glyph size: 7.199 pt. Six editable SVG text elements.

Rendered alignment: PASS, 1.5-pt tolerance, no exemptions. Collision audit: PASS. The existing raw Cairo `Tf=1` diagnostic is resolved by the effective transformed font-size check. All six panels and the complete R preview were visually inspected for title clearance, clipping, point visibility, and removal of the empty bands.

Data QA: S1 coordinates and row order exactly match Figure 1. The S1 plot-signature comparison verifies the current role map before/after title removal; it does not claim invariance to the superseded exclusive-role figure. See `audits/Supplementary_Figure_S1_followup_validation.json`.

Other figures and their source/preview/comparison assets: 229 protected files checked by SHA-256; none changed.

Historical comparison artifact (not bundled): `Supplementary_Figure_S1_before_after.png`

### Supplementary Figure S2

Only overall title, explanatory subtitle/caption and in-panel sample-count labels removed. Three original regional maps, geographic bounds, plant/monitor positions and 20/50-km buffers retained. Regional counts independently reproduce 30/1,091; 35/838; 68/3,150.

Removed text:

- Supplementary Figure S2. Regional distribution of RadNet monitors and nearby power plants
- Solid circles indicate 20-km monitor-centered buffers; dashed circles indicate 50-km buffers.
- Black triangles represent RadNet monitors; grey points represent unique nearby power plants. Regional divisions are descriptive and are not statistical analysis strata.
- Regional plant counts represent unique power plants located within 50 km of at least one RadNet monitor shown in that panel. A plant may contribute to more than one descriptive panel if it lies within 50 km of monitors on opposite sides of a regional boundary.
- 30 RadNet monitors
- 1,091 unique power plants within 50 km
- 35 RadNet monitors
- 838 unique power plants within 50 km
- 68 RadNet monitors
- 3,150 unique power plants within 50 km

Export size: 18.00000 × 7.20000 inches. PNG/TIFF: 600 dpi; TIFF: LZW. All PDF fonts embedded; effective minimum glyph size 11.000 pt. Collision audit: PASS. Measured alignment audit: PASS (single-panel maps are not applicable).

Historical comparison artifact (not bundled): `Supplementary_Figure_S2_before_after.png`

### Supplementary Figure S3

Only overall title/subtitle removed; all eight facets and 56 estimates/CIs retained.

Removed text:

- Supplementary Figure S3. Upwind associations: monitor fixed effects
- Monitor fixed effects; linear and quadratic calendar time; sample-median wind

Export size: 11.00000 × 6.70000 inches. PNG/TIFF: 600 dpi; TIFF: LZW. All PDF fonts embedded; effective minimum glyph size 8.000 pt. Collision audit: PASS. Measured alignment audit: PASS (single-panel maps are not applicable).

Historical comparison artifact (not bundled): `Supplementary_Figure_S3_before_after.png`

### Supplementary Figure S4

Only overall title/subtitle removed; all eight facets and 56 directional contrasts/CIs retained.

Removed text:

- Supplementary Figure S4. Directional contrasts: monitor fixed effects
- Monitor fixed effects; linear and quadratic calendar time; sample-median wind

Export size: 11.00000 × 6.70000 inches. PNG/TIFF: 600 dpi; TIFF: LZW. All PDF fonts embedded; effective minimum glyph size 8.000 pt. Collision audit: PASS. Measured alignment audit: PASS (single-panel maps are not applicable).

Historical comparison artifact (not bundled): `Supplementary_Figure_S4_before_after.png`

### Supplementary Figure S5

Only overall title/subtitle removed; nine labelled panels, 1,204 omission estimates/CIs and both reference-line types retained.

Removed text:

- Supplementary Figure S5. Leave-one-monitor-out influence
- 20 km, 90°; 140 inventory omissions and 133 per fuel–metric; full-sample wind references

Export size: 11.00000 × 9.00000 inches. PNG/TIFF: 600 dpi; TIFF: LZW. All PDF fonts embedded; effective minimum glyph size 8.000 pt. Collision audit: PASS. Measured alignment audit: PASS (single-panel maps are not applicable).

Historical comparison artifact (not bundled): `Supplementary_Figure_S5_before_after.png`

### Supplementary Figure S6

Only overall title/subtitle removed; original 2 × 2 grid, fuel colours/shapes, samples, axis scales, zero reference lines and 112 adjusted estimates/CIs retained. Vector PDF newly exported from original plotting code.

Removed text:

- Combustion associations after adjustment for non-combustion energy
- Wind, solar, hydro, nuclear and geothermal upwind Count simultaneously adjusted.

Export size: 13.00000 × 8.50000 inches. PNG/TIFF: 600 dpi; TIFF: LZW. All PDF fonts embedded; effective minimum glyph size 8.000 pt. Collision audit: PASS. Measured alignment audit: PASS (single-panel maps are not applicable).

Historical comparison artifact (not bundled): `Supplementary_Figure_S6_before_after.png`

### Supplementary Figure S7

Only overall title/subtitle removed; original 2 × 4 grid, current/spline encodings and 112 estimates/CIs retained. Vector PDF newly exported from original plotting code.

Removed text:

- Fuel-specific associations: current vs natural-spline covariate models
- 20–50 km; Monitor × year FE; natural cubic splines with df=3 for continuous background covariates

Export size: 13.00000 × 7.50000 inches. PNG/TIFF: 600 dpi; TIFF: LZW. All PDF fonts embedded; effective minimum glyph size 8.000 pt. Collision audit: PASS. Measured alignment audit: PASS (single-panel maps are not applicable).

Historical comparison artifact (not bundled): `Supplementary_Figure_S7_before_after.png`

### Supplementary Figure S8

Only overall title/subtitle removed; original two metric panels, model legend, 28 estimates/CIs and zero reference lines retained.

Removed text:

- Biomass associations before and after adjustment for GFED fire-sourced PM2.5
- 2001–2023 matched samples; monitor × year fixed effects; effects evaluated at sample-median wind speed

Export size: 10.00000 × 8.50000 inches. PNG/TIFF: 600 dpi; TIFF: LZW. All PDF fonts embedded; effective minimum glyph size 8.000 pt. Collision audit: PASS. Measured alignment audit: PASS (single-panel maps are not applicable).

Historical comparison artifact (not bundled): `Supplementary_Figure_S8_before_after.png`

## Alignment and collision audit interpretation

All final PDFs pass the rendered collision audit with 0 failures and 0 warnings. R measures actual forced PDF viewport rectangles at the export dimensions; all applicable alignment audits pass with the original 1.5-point tolerance. Figure 3 preserves the supplied 1.35:1:1 panel-width ratio. S2 retains its original region-specific map aspects; the existing Eastern-map centering and map widths/gutters are explicitly exempted rather than changing geographic bounds.

The copied collision auditor was corrected locally to split Cairo text traces into actual glyph runs and to honour PDF clipping rectangles. Without those corrections, a combined multi-line trace falsely appears as one huge label and invisible out-of-panel map paths falsely cross panel titles. No collision thresholds were relaxed. Raw `Tf=1` operands in Cairo PDFs are transformed by text matrices; `effective_fonts.json` measures the actual rendered point sizes. The raw text audit logs are diagnostic and do not replace this transformed-font check.

`source_preflight_with_context.json` checks each plotting module together with its shared export helper: 0 FAIL findings. Remaining warnings request a real R parse and journal-width verification; all delivered R scripts parse successfully, and existing figure widths are preserved as instructed. Standalone module scans cannot discover shared export functions and are retained only as raw diagnostics.

## Visual review

Every panel and each assembled figure was inspected from the final R preview and the before/after comparisons. Panel titles, axes, legends, uncertainty bands/error bars and zero/full-estimate reference lines remain legible. Figure 3’s inset legend occupies the upper-right empty part of panel a, with all six original entries and no block below the panels. Figure 1 has only the six requested role keys plus RadNet; biomass remains magenta diamonds.

## Read-only sources and reproducibility

All 90 initially recorded external input files were rehashed at completion and remain unchanged. See `input_manifest.json` and `audits/read_only_input_check.json`. All task outputs, intermediate files, logs and caches are within this directory.

Run from this directory with `zsh scripts/run_all.sh`, or use the commands in `README.md`. The drawing code reads local relative-path input copies; it never writes to the source repositories. `scripts/original/` preserves supplied code for provenance and is not executed. Setup scripts `stage_inputs.py` and `prepare_scripts.py` are historical staging tools; do not rerun them over edited display copies.

No unresolved problems.

## Repository integration verification (2026-10-30)

The portable revision entry point was run successfully from the repository root:
`Rscript display_update_scripts/v4_2026-10-30/run_all.R`, using a Python interpreter
with PyMuPDF and Pillow via `QA_PYTHON`. No model-fitting script was invoked.

- All 14 regenerated PNGs were pixel-identical to the supplied revision artwork.
- All six regenerated Figure 1/S1 display source files were byte-identical to the supplied files.
- All 28 original manuscript-v5 staging copies matched existing repository v5 files byte for byte.
- All 12 applicable scientific plot signatures passed exact comparison at tolerance 0.
- The CSV/subset audit passed all 15 entries with maximum numerical difference 0.
- All 14 PDF/raster exports passed embedded-font, 600-dpi, LZW-TIFF and vector-layer checks; Figure 2 retains the documented slide raster effects.
- All 14 collision audits and all 13 applicable rendered alignment audits returned exit 0; Figure 2 has no panel alignment manifest.
- All 14 effective-font/SVG audits passed. Shared-context preflight checks had 0 failures; their retained warnings concern R parsing and journal-width review. All copied R/Python scripts parsed successfully.
- All 14 before/after comparison images and the separate S1 comparison were generated successfully.

After verification, every committed figure was restored from the read-only input
and verified byte-identical in all four formats. Generated QA artifacts, previews,
caches and sizing files were removed from the commit; future runs write these
under the v4 asset directory. The supplied historical report above remains
distinct from these repository integration checks. `frozen_v3/` and all existing
scientific source files remain unchanged.
