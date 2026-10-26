# Package version and validation status

## English translation of the supplied status note (6 October 2026)

The latest existing code directory, originally named
`code_availability_v3_待验证/`, is a complete copy containing 300 files, excluding
`.DS_Store`. Its source was the existing Project_007 Figure 1/2/6 revision v3
delivery, under `05_Deliverables/2026-10-01_图1图2图6更新_v3/code_availability`
in the original project tree.

English and Chinese manuscript v5 changed language and structure only; their
analytical figures and tables match v3/v4. That packaging operation organized
and copied existing code and added updated display scripts. It did not refit
models or change the existing analysis programs.

The following distinctions applied to the supplied archive:

1. `run_all.R` called `R/06_finalize_v3_displays.R`, which had not been created.
   The complete revised entry point could not be considered validated, and
   `validation/` contained no complete run record for that version.
2. The original README retained the older Figure 6B reference wind of 1 m/s.
   The manuscript, `R/05_refit_figure6.R`, and updated source tables instead used
   the shared baseline-sample median wind for each fuel–metric combination.
3. Independent focused Figure 6 verification had passed. It did not validate
   the complete package entry point.
4. Original `R/` scripts, prepared model matrices, display source tables,
   provenance, and historical scripts were retained. Historical scripts depend
   on their original research directories and external inputs. Additional
   display-update scripts also document the revision history.
5. An older archive with the same code-availability name was not the current
   v3 directory. A separate archive was created containing this status note
   and the actual materials assembled at that time.
6. A public repository, code licence, and persistent identifier still awaited
   author action. This was a local draft and could not support a promise that
   every latest result was fully reproducible from a single entry point.

Additional scripts, now in `display_update_scripts/`, were:

| Script | Purpose |
| --- | --- |
| `revise_Figure1.R` | Revised Figure 1 mutually exclusive process classification and maps |
| `extract_Figure2.py` | Extract Figure 2 from the third slide of the original presentation |
| `revise_Figure6.R` | Figure 6 models and shared-median-wind revision |
| `Figure6_display_v3.R` | Revised Figure 6 rendering |
| `verify_Figure6_covariance.R` | Independent full-covariance Figure 6 verification |
| `audit_effective_pdf_fonts.py` | Figure 6 PDF font-geometry audit |

`source_data_manuscript_v5/` retains all 28 direct display source files for
manuscript v5, distinguished from historical package outputs. Supplementary
figures and tables were also stored separately in the submission directory;
the original source files and earlier manuscripts remained in Project_007.

## Changes for this local repository

The missing final-display entry-point module now sources and invokes the
existing Figure 1 and Figure 6 renderers with repository-relative paths. It
renders the frozen manuscript source tables; it does not refit models. The
README now describes the correct Figure 6B wind and the actual reproduction
scope. Active display-update paths now resolve within the repository. The
historical external-tree builder `R/00_prepare_model_matrix.R` still requires
`RADNET_RESEARCH_ROOT` and is outside the entry point.

English filename mappings are in `provenance/renamed_paths.csv`; changed path
expressions are in `provenance/path_changes.csv`. In that report, `[USER_HOME]`
redacts the original personal home prefix. Provenance metadata uses
`<ORIGINAL_RESEARCH_ROOT>`, `<ORIGINAL_PROJECT_ROOT>`, `<ORIGINAL_FIGURE_SKILL>`,
`<ORIGINAL_QA_PYTHON>`, and `<ORIGINAL_USER_HOME>` for historical locations. Historical script bodies
and serialized model records retain their original provenance. Numerical
data, estimates, models, and formulas have not been changed.

MIT licensing has been added for the software. Both prepared analysis matrices
are included for publication as instructed by the authors. Third-party source
products retain their provider terms and attribution requirements. The
repository remains local, with no remote, public deposit, or DOI. The complete
`Rscript run_all.R` pipeline has not been run end-to-end for this repository.
An independent read-only covariance-vector check was run on the archived model
records against all 56 v5 source rows. It passed (maximum error 3.47e-16),
without refitting models; the report is in `provenance/release_audit/`.

The optional Figure 2 extraction script expects `data/external_inputs/study_design.pptx`
and its LibreOffice PDF export `data/external_inputs/study_design.pdf`. Neither
original input is bundled; frozen Figure 2 exports and slide XML are retained.
The optional `revise_Figure6.R` uses the bundled 1 m/s historical reference for
its original revision comparisons and writes into `outputs/display_updates/`;
it is not called by `run_all.R`.

The manuscript's Data Availability paragraph links the fire-PM2.5 product to
`https://doi.org/10.5281/zenodo.15493914`, whereas its embedded reference includes
a different QFED record (`10.5281/zenodo.15496596`). The README follows the
Data Availability GFED link; this source-record discrepancy remains for author
review. No AQS URL was present in the supplied English manuscript; the README
uses the official EPA [AQS access page](https://www.epa.gov/aqs/obtaining-aqs-data).
