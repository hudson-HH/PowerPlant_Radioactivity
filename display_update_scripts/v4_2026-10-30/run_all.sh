#!/bin/zsh
set -euo pipefail
cd "${0:A:h}/../.."
ASSETS="data/display_assets/v4_2026-10-30"
SCRIPTS="display_update_scripts/v4_2026-10-30"
mkdir -p "$ASSETS/.tmp" "$ASSETS/qa/cache" "$ASSETS/qa/logs" "$ASSETS/qa/audits"
export TMPDIR="$PWD/$ASSETS/.tmp"
export XDG_CACHE_HOME="$PWD/$ASSETS/qa/cache"
export R_USER_CACHE_DIR="$PWD/$ASSETS/qa/cache"
export PYTHONDONTWRITEBYTECODE=1
R_BIN="${R_BIN:-Rscript}"
QA_PYTHON="${QA_PYTHON:-python3}"
export QA_PYTHON
"$QA_PYTHON" "$SCRIPTS/preflight.py"
"$R_BIN" "$SCRIPTS/run_all.R"
"$R_BIN" "$SCRIPTS/compare_csv_sources.R"
"$R_BIN" "$SCRIPTS/make_comparisons.R"
"$QA_PYTHON" "$SCRIPTS/run_audits.py"
"$QA_PYTHON" "$SCRIPTS/qa_helpers/effective_font_audit.py"
