# Local repository preparation record

Date: 26 October 2026. Associated research project: Project_007.
The task-specific instruction placed extraction and the local Git repository
in the current workspace; the original research tree and submission ZIP were
read only. No remote was created and no push was performed.

ZIP member names lacking the UTF-8 flag were decoded by encoding their
ZIP-decoded cp437 names back to bytes and decoding those bytes as UTF-8.
The v3 code directory became the repository root, the v5 display sources and
revision scripts became English-named subfolders, and the Chinese status
note was translated into `package_status.md`. All input files are accounted
for in `provenance/release_audit/input_file_integrity.csv`. Input hashes refer
to the supplied archive; final hashes are in `SHA256SUMS.txt`.

Edits to existing active code are restricted to root discovery, input/output
paths, Python interpreter/helper lookup, and creation of output directories.
`R/06_finalize_v3_displays.R` adds entry-point wiring to the existing Figure 1/6
rendering functions. Neither `run_all.R` nor model formulas were changed.
The final display module uses the existing frozen v5 sources after preceding
refit checks; it does not construct or alter analytical values.

All original data files, frozen artwork, historical numerical outputs,
serialized model records, and legacy script bodies were retained byte for
byte. Changed provenance metadata replaces private home prefixes with
documented placeholders and updates packaged filename references. Legacy
external-tree path references remain historical. Each edited path expression
is recorded in `provenance/path_changes.csv`; `[USER_HOME]` masks the original
personal prefix in that public report.

Verification included active-R parsing without evaluation, Python AST parsing,
prepared-matrix dimensions and source hashes, data/legacy hash preservation,
entry-point and helper-file presence, English path names, the 50 MB file limit,
and a text/CSV-gzip credential, email, and personal-path scan. The independent
Figure 6 covariance check used archived model records only. The R analysis
pipeline, model refits, and figure rendering were not run.

Preparation-tool failures were corrected before proceeding: Python rejected
an invalid encoding alias during an initial scan, and the patch tool rejected
a single patch that both deleted and added the README at the same path.
Neither failed operation changed input data or analytical code. The final scan
also caught eight JSON-escaped historical project paths that the first literal
prefix replacement missed; their personal home prefixes were redacted and the
scan was repeated. Numerical JSON fields were preserved. The corrected scan
and file writes completed successfully.

The pre-commit inventory check also caught a broad `outputs/` ignore rule that
would have omitted supplied historical artwork and tables. It was narrowed
to new intermediate/work directories before committing. Every supplied output
is retained and tracked; the final checksum manifest was regenerated afterward.

Future reuse: decode and check ZIP paths before extraction; inventory and hash
all inputs before editing; use separate path/rename manifests; preserve data
bytes; distinguish static checks and archived-model verification from a full
entry-point reproduction. No global research files were modified.
