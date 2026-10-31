# Supplementary-table revision 2026-10-31

This workflow reformats the eight Supplementary Information (SI) tables for
print and repairs figure layout in external Word files. No estimate,
analytical sample or model changed. Manuscript and SI Word files are not part
of this repository and must be supplied as command-line arguments.

The exact printed tables are the eight CSVs in
`data/display_assets/v5_2026-10-31/tables/`. Full-precision values and all
exported columns remain in `data/display_assets/frozen_v3/tables/` and
`source_data_manuscript_v5/Supplementary_Figure_S8_GFED_estimates.csv`.
Table 1 is unchanged.

## Scripts and usage

Requires Python 3 and `lxml`. Examples below run from the repository root;
replace uppercase filenames with paths to your external Word documents.
The scripts also work from other working directories when invoked by their
script paths. Input documents are read only; each editing script requires a
separate output document. Output parent directories must already exist.

| Script | Purpose |
| --- | --- |
| `rebuild_tables.py` | Rebuild all eight Word tables from the original SI's full-precision cells; cross-check those source cells against repository CSVs, update notes and table layout, and write a JSON report. The repository root is resolved from the script location. |
| `verify_tables.py` | Independently compare every printed numeric cell with the original full-precision SI tables, including the paired S6 models. Read only; prints counts, mismatches and example rows. |
| `fix_figs.py` | Set figure paragraphs to single spacing, remove extra spacing and inappropriate heading styles, keep images with legends, split a shared image/caption paragraph, and restore the SI S1 bookmark. Also accepts the main manuscript. |
| `fix_extents.py` | Use each inline PNG's aspect ratio and existing frame width to correct frame height; set inner picture `a:ext` dimensions equal to `wp:extent` and clear effect extents. Also accepts the main manuscript. Image bytes are unchanged. |
| `table_red.py` | Prepare a reviewer copy with table headers red and table body text black; text outside tables is untouched. |

```sh
python3 display_update_scripts/v5_2026-10-31/rebuild_tables.py SOURCE_SI.docx TABLES_SI.docx --report TABLES_SI.report.json
python3 display_update_scripts/v5_2026-10-31/verify_tables.py SOURCE_SI.docx TABLES_SI.docx
python3 display_update_scripts/v5_2026-10-31/fix_figs.py TABLES_SI.docx SPACED_SI.docx
python3 display_update_scripts/v5_2026-10-31/fix_extents.py SPACED_SI.docx FINAL_SI.docx
python3 display_update_scripts/v5_2026-10-31/table_red.py FINAL_SI.docx REVIEWER_SI.docx
```

Omitting `--report` writes the report beside the rebuilt document with suffix
`.report.json`. Rebuilding and verification require the original SI with its
eight full-precision tables, original column names and expected notes, rather
than an already reformatted SI. Inspect the rebuild report's `csv_mismatches`
and `problems`, and require zero numeric mismatches from verification. The
original checks and formatting logic are retained; the verifier reports
numeric mismatches in its output, so its exit status alone is insufficient.
The rebuild report notes that S5 t statistics vary and were explicitly omitted.
These scripts do not refit models or regenerate the committed display CSVs.

Two supplied-script details are retained to preserve the original logic:
the S7 coefficient header built by `rebuild_tables.py` omits the unit
`(mBq/m³ per µg/m³)` present in the committed printed CSV; and the S7 CSV
uses `sample median wind` where the original DOCX schema expects
`sample_median`, which can appear as a metadata mismatch in its cross-check
report. Neither difference changes a numeric value. Match the S7 header to
the printed CSV when preparing the final Word document and review metadata
mismatches separately from numeric mismatches.

## Print formatting

- Use readable column names with units and value labels: `coal` → Coal,
  `count` → Plant count, `heat` → Daily heat input. Write S5 covariate names
  and units in full.
- Use thousands separators for counts; three significant figures for
  estimates, confidence limits, SEs, coefficients and percentages; three
  decimal places for wind speed; and three decimal places for P values,
  with values below 0.001 shown as `<0.001`. Use the true minus sign (−).
- Combine each estimate and its 95% CI in one cell. All eight tables are
  landscape; the original SI already has landscape sections except S1b,
  whose section is changed by the rebuild script.
- Omit the following constant or unnecessary printed columns while retaining
  them in the full-precision sources:

| Table | Omitted columns or layout change |
| --- | --- |
| S2 | Wind reference |
| S3 | Radius, status, wind reference |
| S4 | Radius, wind speed |
| S5 | t, model N, monitors, outcome unit, N_pre, note |
| S6 | N_input, monitors, SE, P; constant wind reference and model-specific spline metadata are described in notes. Place the main and natural-spline models side by side because model N and wind match in all 56 pairs. |
| S7 | Monitors, wind speed (3.152 m/s), absolute change, attenuation_pct; the constant wind reference is described in the note. |

## Word figure layout

Figure paragraphs are single-spaced and kept with their legends. Main
Figures 2–6 had been vertically squeezed by 1–11% in Word because the inner
picture dimensions differed from the frame; matching `a:ext` to `wp:extent`
fixes this without changing the figure files. Apply the spacing and extent
helpers to the main manuscript separately when needed. The v4 workflow has
no image-replacement script and no code that sizes pictures in a DOCX; no
v4 script required an extent change.

Repository integration checked Python compilation and CSV copying. It did
not run these transformations or numeric verification on the external Word
documents, which are not bundled here.
