# Particle radioactivity near US combustion power plants

Huai Hu, Fenghua Wang, Jinghan Wang, Yeqing Gu, Hongbing Xu*, Qiang Liu*  
Institute of Radiation Medicine, Chinese Academy of Medical Sciences & Peking Union Medical College (CAMS & PUMC)

This study examines associations between US combustion power plant activity
and ambient particle radioactivity using RadNet air-filter gross-beta
measurements, plant inventories, daily fuel-specific operating records, wind
direction and speed, and environmental covariates. Prepared filter-period
matrices support analyses of coal, natural gas, oil, and biomass plant counts
and heat input, directional contrasts, and sensitivity analyses. Gross-beta
activity is an integrated particle measurement rather than an isotope-specific
or PM2.5-specific measurement; the observational associations do not establish
causation. Heat input measures fuel energy, not electricity generation or
radioactive emissions.

Citation placeholder: Hu H, et al. Particle radioactivity near US combustion power plants. Manuscript in preparation.

| Location | Contents |
| --- | --- |
| `run_all.R` | Single-entry validation, redraw, refit, final display, and session-recording sequence |
| `R/` | Publication reproduction scripts and final Figure 1/6 rendering modules |
| `data/prepared_model_data/` | Both published compressed analysis matrices |
| `data/source_data/` | Frozen model estimates, numerical display sources, and supplementary results |
| `data/display_assets/frozen_v3/` | Existing manuscript Figures 1–6, Supplementary Figures S1–S8, Table 1, and Tables S1a–S7 |
| `source_data_manuscript_v5/` | All 28 direct manuscript-v5 display source files |
| `data_dictionary.csv` | Matrix variables, units, and definitions where supplied |
| `display_update_scripts/` | Additional Figure 1/2/6 revision, extraction, and audit scripts |
| `legacy_scripts/` | Historical upstream preparation and analysis scripts; external inputs required |
| `provenance/` | Input manifests, historical verification evidence, filename mapping, and path-edit report |
| `tools/nature-figure/scripts/` | Bundled alignment and collision audit helpers |
| `docs/package_status.md` | Translated archive status and current limitations |
| `outputs/`, `validation/` | Supplied historical outputs are retained; new intermediate/work directories and validation reports are ignored by Git |
| `SHA256SUMS.txt` | Final release-file hashes, excluding itself and `.git/` |

Data sources are linked below using the manuscript's Data Availability URLs,
except AQS, for which an official EPA access page is supplied because the
manuscript did not contain a URL.

| Source | Access |
| --- | --- |
| EPA RadNet air-filter measurements and metadata | [RadNet](https://www.epa.gov/radnet), [Envirofacts](https://enviro.epa.gov) |
| EIA-860 plant and generator inventories | [EIA-860](https://www.eia.gov/electricity/data/eia860/) |
| EIA-923 operation and fuel records | [EIA-923](https://www.eia.gov/electricity/data/eia923/) |
| EPA Clean Air Markets Program Data (CAMPD) | [CAMPD](https://campd.epa.gov), [power-sector data](https://www.epa.gov/power-sector/power-sector-data) |
| North American Regional Reanalysis (NARR) | [NOAA NARR](https://psl.noaa.gov/data/gridded/data.narr.html) |
| HYSPLIT atmospheric transport | [NOAA HYSPLIT](https://www.ready.noaa.gov/HYSPLIT.php) |
| EPA Air Quality System (AQS) | [Obtaining AQS data](https://www.epa.gov/aqs/obtaining-aqs-data) |
| USGS terrestrial radioactivity | [USGS Open-File Report 2005-1413](https://pubs.usgs.gov/of/2005/1413/) |
| SILSO International Sunspot Number | [WDC-SILSO, Royal Observatory of Belgium](https://www.sidc.be/SILSO/) |
| GFED fire-sourced PM2.5 sensitivity product | [Zenodo dataset](https://doi.org/10.5281/zenodo.15493914) |

The fire-product link discrepancy in the manuscript is recorded in
`docs/package_status.md`. Original upstream downloads, gridded meteorological
fields, complete trajectories, and the full plant–monitor construction are
not bundled. Original providers' terms and attribution requirements continue
to apply to their products, including SILSO.

Both files in `data/prepared_model_data/` are included in this release by the
authors' decision. Counts below are from `provenance/prepared_data_manifest.csv`.

| File | Rows | Columns | Description |
| --- | ---: | ---: | --- |
| `filter_period_primary_2001_2024.csv.gz` | 107,407 | 136 | Prepared eligible filter-period matrix, 133 monitors, outcomes, environmental covariates, and four-fuel upwind/downwind exposures |
| `filter_period_biomass_GFED_2001_2023.csv.gz` | 103,196 | 33 | Overlap matrix with biomass exposures and fire-sourced PM2.5 coverage; paired fits select matching complete observations |

These are prepared analytical matrices, not raw EPA/EIA downloads. Row counts
describe eligible inputs, not necessarily actual model sample sizes. The
`data_dictionary.csv` lists matrix variables and supplied definitions; entries
marked “See source variable / README” remain incomplete descriptions. Keys
include monitor identifiers and filter-period dates; collection intervals are
half-open `[start,end)`. CAMPD counts are distinct eligible reporting plants
without an additional positive-heat threshold; the operational biomass
category maps CAMPD `WOOD`.

The reference environment is **R 4.5.3**, **fixest 0.14.2**, and **lme4 2.0-1**.
Scripts additionally use `data.table`, `Matrix`, `digest`, `ggplot2`,
`patchwork`, `svglite`, `maps`, `jsonlite`, and `ragg`, plus R's standard
`splines`, `grid`, `grDevices`, `parallel`, and `tools` packages. Historical
session details are retained at
`provenance/revision_v3/legacy_validation/sessionInfo.txt`; they are not a new
validation of this repository. Install requirements in your own environment;
the entry point does not install software or download data.

Final display export requires Cairo PDF support and available Helvetica or a
compatible font. PDF alignment, collision, and font audits use **Python 3.10+**
and **PyMuPDF** (`fitz`); the optional Figure 2 extraction also uses **Pillow**
and **lxml**, with a LibreOffice PDF export of the original presentation.
Bundled QA helpers are used by default. `NATURE_FIGURE_QA_PYTHON` can select a
Python interpreter; `NATURE_FIGURE_SKILL` can select a compatible helper tree.

From the repository root:

```sh
Rscript run_all.R
```

The entry point resolves the root from its own script location. It writes
derived files under `outputs/` and `validation/`. The existing redraw/refit
stages use `outputs/v1_intermediates/`; final Figure 1/6 displays use
`outputs/figures/`. Source data and frozen display assets are retained.
Existing historical files under `outputs/figures/` and `outputs/tables/` are
also included; their presence does not document a run of this revised entry point.

| Script or asset | Figures/tables and operation |
| --- | --- |
| `R/01_validate_source_data.R` | Check source schemas, expected rows, Figure 4/6 winds, and matched GFED samples |
| `R/02_reproduce_figures_tables.R` | Redraw Figures 3–6, S3–S5 and S8 from estimates; export Table 1 and Tables S2–S4 plus environmental, spline, and GFED source tables corresponding to S5–S7 |
| `R/02a_redraw_overall_inventory.R` | Figure 3 redraw from numerical exports; sourced by stage 02 |
| `R/02b_redraw_supplementary_diagnostics.R` | S3/S4 monitor-FE and S5 leave-one-monitor-out redraws; sourced by stage 02 |
| `R/03_refit_models.R` | Refit Figure 4 and S7 spline models (56 each), Figure 5 joint-direction models (56), S8 paired GFED comparisons (28 fits), and the environmental model underlying Table S5; compare frozen estimates |
| `R/05_refit_figure6.R` | Refit Figure 6 models, recover actual panel sample sizes, and compare all 56 displayed estimates |
| `R/06_finalize_v3_displays.R` | Source and invoke `R/Figure1_display_v3.R` and `R/Figure6_display_v3.R` on frozen v5 sources for final Figures 1/6; no models are refitted in this step |
| `R/04_session_info.R` | Write the runtime session record after preceding steps complete |
| `data/display_assets/frozen_v3/figures/` | Frozen Figures 1–6 and S1–S8; Figure 2 and S1/S2/S6/S7 final artwork are distributed assets rather than regenerated by the entry point |
| `data/display_assets/frozen_v3/tables/` | Frozen Table 1 and Tables S1a, S1b, S2, S3, S4, S5, S6, S7 with manuscript-facing names |
| `display_update_scripts/revise_Figure1.R` | Optional descriptive Figure 1 revision from supplied plant/monitor sources |
| `display_update_scripts/extract_Figure2.py` | Optional Figure 2 vector extraction; original PPTX/PDF inputs are not bundled |
| `display_update_scripts/revise_Figure6.R` | Optional Figure 6 revision from the prepared matrix and archived 1 m/s reference; outside `run_all.R` |
| `R/verify_Figure6_covariance.R` | Optional independent verification from archived coefficient/covariance records; no model refit |
| `legacy_scripts/` | Historical producing scripts for upstream construction, maps, descriptive/exposure support (S1a/S1b), period-specific results (S3), environmental (S5), spline (S6), GFED (S7), and other supplementary analyses |

Figures 3, S3/S4, and S5 are redrawn without refitting their original models.
Figure 3 uses a separate 140-monitor quasi-daily inventory sample and must not
be refitted with the 133-monitor filter-period matrix. Monitor-only models,
leave-one-monitor-out fits, period-specific fits, and non-combustion adjustment
fits are not all rerun by the entry point. The optional curator script
`R/00_prepare_model_matrix.R` requires the original external research tree via
`RADNET_RESEARCH_ROOT`; its upstream data are not part of this repository.

Figure 3 uses sample-mean wind (2.985606 m/s); Figures 4–5 use model-specific
sample-median wind. Figure 6B uses the **shared baseline-sample median wind
for each fuel–metric**, before fixed-effect singleton removal: **3.154 m/s**
for all count metrics and biomass heat input, **3.153 m/s** for coal heat input,
**3.139 m/s** for gas heat input, and **3.136 m/s** for oil heat input. Full
precision and samples are in
`source_data_manuscript_v5/Figure_6_reference_winds_and_samples.csv`. Historical
1 m/s comparisons remain in provenance. Contrasts and confidence intervals
include the fitted covariance. The 20-km biomass exposure occurs at only two
monitors, limiting effective exposure replication and interval reliability.
Analyses were iterative and do not apply a multiplicity adjustment.

**Validation status:** Figure 6 covariance was independently verified; the
archive includes coefficients/covariances and its verification script. A fresh
read-only check of all 56 v5 rows passed with maximum absolute error
3.47e-16; see `provenance/release_audit/Figure6_covariance_verification.json`.
**Full single-entry reproduction has not yet been run end-to-end for this
repository.** Historical reports are retained in `provenance/revision_v3/`,
and do not establish current entry-point success. Model checks use 1e-7 for
OLS and 1e-5 for Figure 6 mixed models. The historical H16-versus-H15 FAIL is
retained because samples differ; historical input N must not be confused with
actual `nobs`. The environmental model's actual N is 107,367, whereas the
prepared primary matrix has 107,407 rows. Original finite-cluster-t p values
and the package's normal-Wald presentation are distinguished in provenance.

The code is released under the [MIT License](LICENSE). This repository includes
both prepared datasets and all supplied data files. MIT does not supersede
third-party data terms. There is currently no public remote, archival DOI, or
versioned public deposit. `SHA256SUMS.txt` identifies this local tree.

Corresponding authors: Hongbing Xu, xuhongbing@irm-cams.ac.cn;
Qiang Liu, liuqiang@irm-cams.ac.cn.
