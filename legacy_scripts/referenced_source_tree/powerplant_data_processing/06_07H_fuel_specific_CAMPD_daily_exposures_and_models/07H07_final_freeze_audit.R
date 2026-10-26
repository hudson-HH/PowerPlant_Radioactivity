#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999
)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("缺少data.table。", call. = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
})

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

HROOT <- file.path(
  ROOT,
  "发电站数据处理",
  "06_07H_燃料特异CAMPD逐日暴露与模型"
)

H06 <- file.path(
  HROOT,
  "07H06_主分析与距离加权敏感性综合审计"
)

H07 <- file.path(
  HROOT,
  "07H07_论文结果表与图形准备"
)

SCRIPT07 <- file.path(
  HROOT,
  "07H07_论文结果表与图形准备.R"
)

FILE_H06_HARD <- file.path(
  H06,
  "10_07H06_hard_checks.csv"
)

FILE_H06_RW <- file.path(
  H06,
  "02_07H06_raw与distance_weighted稳健性比较.csv"
)

FILE_H07_HARD <- file.path(
  H07,
  "07H07_hard_checks.csv"
)

FILE_TABLE03 <- file.path(
  H07,
  "01_论文结果表",
  "03_raw与distance_weighted稳健性汇总.csv"
)

FILE_REPORT <- file.path(
  H07,
  "07H07_论文结果与图表说明.md"
)

OUT_AUDIT <- file.path(
  H07,
  "07H07_final_freeze_audit.csv"
)

OUT_MARKER <- file.path(
  H07,
  "07H07_FINAL_FREEZE.md"
)

needed <- c(
  FILE_H06_HARD,
  FILE_H06_RW,
  FILE_H07_HARD,
  FILE_TABLE03,
  FILE_REPORT,
  SCRIPT07
)

if (any(!file.exists(needed))) {
  cat("缺失文件：\n")
  print(needed[!file.exists(needed)])
  stop("冻结审计输入不完整。", call. = FALSE)
}

H06H <- fread(
  FILE_H06_HARD,
  showProgress = FALSE
)

H07H <- fread(
  FILE_H07_HARD,
  showProgress = FALSE
)

RW <- fread(
  FILE_H06_RW,
  showProgress = FALSE
)

T03 <- fread(
  FILE_TABLE03,
  showProgress = FALSE
)

n_h06_fail <- H06H[
  status != "PASS",
  .N
]

n_h07_fail <- H07H[
  status != "PASS",
  .N
]

n_total <- nrow(RW)

n_comparable <- RW[
  raw_status == "PASS" &
  weighted_status == "PASS",
  .N
]

n_incomplete <- RW[
  robustness_pattern ==
    "PAIR_INCOMPLETE",
  .N
]

n_same <- RW[
  raw_status == "PASS" &
  weighted_status == "PASS" &
  same_sign_fixed_increment %in% TRUE,
  .N
]

n_opposite <- RW[
  raw_status == "PASS" &
  weighted_status == "PASS" &
  same_sign_fixed_increment %in% FALSE,
  .N
]

n_both_ci <- RW[
  robustness_pattern %in%
    c(
      "BOTH_CI_EXCLUDE_ZERO_SAME_SIGN",
      "BOTH_CI_EXCLUDE_ZERO_OPPOSITE_SIGN"
    ),
  .N
]

n_raw_only <- RW[
  robustness_pattern ==
    "RAW_ONLY_CI_EXCLUDES_ZERO",
  .N
]

n_weighted_only <- RW[
  robustness_pattern ==
    "WEIGHTED_ONLY_CI_EXCLUDES_ZERO",
  .N
]

n_neither <- RW[
  robustness_pattern %in%
    c(
      "NEITHER_CI_EXCLUDES_ZERO_SAME_SIGN",
      "NEITHER_CI_EXCLUDES_ZERO_OPPOSITE_SIGN"
    ),
  .N
]

n_ci_partition <- (
  n_both_ci +
  n_raw_only +
  n_weighted_only +
  n_neither
)

OVERALL <- T03[
  comparison_family == "ALL" &
  branch == "ALL" &
  fuel == "ALL"
]

if (nrow(OVERALL) != 1L) {
  stop(
    paste0(
      "Table03总体行数量异常：",
      nrow(OVERALL)
    ),
    call. = FALSE
  )
}

TABLE_DIR <- file.path(
  H07,
  "01_论文结果表"
)

PLOT_DIR <- file.path(
  H07,
  "02_绘图数据"
)

FIG_DIR <- file.path(
  H07,
  "03_正式图形"
)

TABLE_FILES <- file.path(
  TABLE_DIR,
  c(
    "01_论文主结果_upwind单暴露_全部预设半径.csv",
    "02_论文补充结果_Joint_Count_Heat_upwind_全部半径.csv",
    "03_raw与distance_weighted稳健性汇总.csv",
    "04_upwind_downwind方向负对照汇总.csv",
    "05_quasi_filter时间尺度一致性汇总.csv",
    "06_模型覆盖与SKIP情况.csv",
    "07_全部预设累计半径结果模式.csv",
    "08_Joint_Count_Heat共线性解释表.csv",
    "09_Σ1d与ΣHeatd共线性审计.csv"
  )
)

PLOT_FILES <- file.path(
  PLOT_DIR,
  c(
    "01_upwind四类单暴露_1SD绘图数据.csv",
    "02_Count_vs_Σ1d_1SD绘图数据.csv",
    "03_Heat_vs_ΣHeatd_1SD绘图数据.csv",
    "04_upwind_downwind四类单暴露_1SD绘图数据.csv",
    "05_quasi_filter_upwind四类单暴露_1SD绘图数据.csv"
  )
)

FIG_BASE <- c(
  "Fig1_raw主分析_upwind_Count_Heat_1SD",
  "Fig2_Count_vs_Σ1d_距离加权敏感性_1SD",
  "Fig3_Heat_vs_ΣHeatd_距离加权敏感性_1SD",
  "FigS1_upwind_downwind方向负对照_四类单暴露_1SD",
  "FigS2_quasi_filter时间尺度一致性_upwind_1SD",
  "FigS3_Joint_Count_Heat_upwind_1SD"
)

PNG_FILES <- file.path(
  FIG_DIR,
  paste0(
    FIG_BASE,
    ".png"
  )
)

PDF_FILES <- file.path(
  FIG_DIR,
  paste0(
    FIG_BASE,
    ".pdf"
  )
)

SCRIPT_TEXT <- readLines(
  SCRIPT07,
  warn = FALSE,
  encoding = "UTF-8"
)

REPORT_TEXT <- readLines(
  FILE_REPORT,
  warn = FALSE,
  encoding = "UTF-8"
)

syntax_ok <- tryCatch(
  {
    invisible(
      parse(
        file = SCRIPT07
      )
    )
    TRUE
  },
  error = function(e) {
    FALSE
  }
)

get1 <- function(name) {

  if (!(name %in% names(OVERALL))) {
    return(NA_character_)
  }

  as.character(
    OVERALL[[name]][1]
  )
}

AUDIT <- data.table(
  check = c(
    "07H06_hard_checks_all_pass",
    "07H07_hard_checks_all_pass",
    "raw_weighted_total_pairs_480",
    "raw_weighted_comparable_pairs_426",
    "raw_weighted_pair_incomplete_54",
    "raw_weighted_same_sign_401",
    "raw_weighted_opposite_sign_25",
    "both_CI_exclude_zero_191",
    "raw_only_CI_excludes_zero_45",
    "weighted_only_CI_excludes_zero_12",
    "neither_CI_excludes_zero_178",
    "CI_partition_equals_comparable_426",
    "total_partition_equals_480",
    "Table03_overall_n_pairs_480",
    "Table03_overall_n_both_pass_426",
    "Table03_overall_n_pair_incomplete_54",
    "Table03_overall_n_same_sign_401",
    "Table03_overall_n_opposite_sign_25",
    "Table03_overall_both_CI_191",
    "Table03_overall_raw_only_45",
    "Table03_overall_weighted_only_12",
    "Table03_overall_neither_178",
    "paper_tables_9_exist",
    "plot_data_files_5_exist",
    "png_files_6_exist",
    "pdf_files_6_exist",
    "png_files_6_nonzero",
    "pdf_files_6_nonzero",
    "script_syntax_pass",
    "script_contains_PAIR_INCOMPLETE_fix",
    "script_contains_sigma_1d",
    "script_contains_sigma_heat",
    "script_uses_cairo_pdf",
    "R_X11_capability_true",
    "R_cairo_capability_true",
    "report_contains_PASS",
    "report_contains_correct_neither_178",
    "report_contains_pair_incomplete_54"
  ),

  observed = c(
    as.character(n_h06_fail),
    as.character(n_h07_fail),
    as.character(n_total),
    as.character(n_comparable),
    as.character(n_incomplete),
    as.character(n_same),
    as.character(n_opposite),
    as.character(n_both_ci),
    as.character(n_raw_only),
    as.character(n_weighted_only),
    as.character(n_neither),
    as.character(n_ci_partition),
    as.character(n_comparable + n_incomplete),
    get1("n_pairs"),
    get1("n_both_pass"),
    get1("n_pair_incomplete"),
    get1("n_same_sign"),
    get1("n_opposite_sign"),
    get1("n_both_ci_exclude_zero"),
    get1("n_raw_only_ci_excludes_zero"),
    get1("n_weighted_only_ci_excludes_zero"),
    get1("n_neither_ci_excludes_zero"),
    as.character(sum(file.exists(TABLE_FILES))),
    as.character(sum(file.exists(PLOT_FILES))),
    as.character(sum(file.exists(PNG_FILES))),
    as.character(sum(file.exists(PDF_FILES))),
    as.character(
      sum(
        file.exists(PNG_FILES) &
        file.info(PNG_FILES)$size > 0
      )
    ),
    as.character(
      sum(
        file.exists(PDF_FILES) &
        file.info(PDF_FILES)$size > 0
      )
    ),
    as.character(as.integer(syntax_ok)),
    as.character(
      as.integer(
        any(
          grepl(
            "N_PAIR_INCOMPLETE_CI <- RW[",
            SCRIPT_TEXT,
            fixed = TRUE
          )
        )
      )
    ),
    as.character(
      as.integer(
        any(
          grepl(
            "Σ1/d",
            SCRIPT_TEXT,
            fixed = TRUE
          )
        )
      )
    ),
    as.character(
      as.integer(
        any(
          grepl(
            "ΣHeat/d",
            SCRIPT_TEXT,
            fixed = TRUE
          )
        )
      )
    ),
    as.character(
      as.integer(
        any(
          grepl(
            "cairo_pdf",
            SCRIPT_TEXT,
            fixed = TRUE
          )
        )
      )
    ),
    as.character(
      as.integer(
        isTRUE(
          capabilities("X11")
        )
      )
    ),
    as.character(
      as.integer(
        isTRUE(
          capabilities("cairo")
        )
      )
    ),
    as.character(
      as.integer(
        any(
          grepl(
            "07H07 PAPER-RESULT PREPARATION: PASS",
            REPORT_TEXT,
            fixed = TRUE
          )
        )
      )
    ),
    as.character(
      as.integer(
        any(
          grepl(
            "Neither 95% CI excludes zero among comparable pairs: 178 pairs.",
            REPORT_TEXT,
            fixed = TRUE
          )
        )
      )
    ),
    as.character(
      as.integer(
        any(
          grepl(
            "Pair-incomplete (at least one model not estimable): 54 pairs.",
            REPORT_TEXT,
            fixed = TRUE
          )
        )
      )
    )
  ),

  expected = c(
    "0",
    "0",
    "480",
    "426",
    "54",
    "401",
    "25",
    "191",
    "45",
    "12",
    "178",
    "426",
    "480",
    "480",
    "426",
    "54",
    "401",
    "25",
    "191",
    "45",
    "12",
    "178",
    "9",
    "5",
    "6",
    "6",
    "6",
    "6",
    "1",
    "1",
    "1",
    "1",
    "1",
    "1",
    "1",
    "1",
    "1",
    "1"
  )
)

AUDIT[
  ,
  status :=
    fifelse(
      observed == expected,
      "PASS",
      "FAIL"
    )
]

fwrite(
  AUDIT,
  OUT_AUDIT
)

cat(
  "\n============================================================\n",
  "07H07 FINAL FREEZE AUDIT\n",
  "============================================================\n",
  sep = ""
)

print(
  AUDIT,
  nrows = Inf
)

N_FAIL <- AUDIT[
  status == "FAIL",
  .N
]

if (N_FAIL > 0L) {

  if (file.exists(OUT_MARKER)) {
    unlink(OUT_MARKER)
  }

  cat(
    "\n============================================================\n",
    "07H07 FINAL FREEZE: FAIL\n",
    "============================================================\n",
    "Hard FAIL = ",
    N_FAIL,
    "\n",
    sep = ""
  )

  quit(
    status = 1L
  )
}

MARKER <- c(
  "# 07H07 FINAL FREEZE",
  "",
  "**07H07 FINAL FREEZE: PASS**",
  "",
  paste0(
    "Frozen at: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  "",
  "## Frozen scope",
  "",
  "- 07H07 paper-result tables, plot data, and formal figures are frozen.",
  "- 07H01–07H06 upstream analyses were not modified.",
  "- No statistical model was refitted during the reporting-layer correction.",
  "",
  "## Final raw-versus-weighted classification",
  "",
  "- Total pairs: 480.",
  "- Both models estimable: 426.",
  "- Pair-incomplete: 54.",
  "- Same effect direction: 401/426 (94.1%).",
  "- Opposite effect direction: 25/426 (5.9%).",
  "- Both 95% CIs exclude zero: 191.",
  "- Raw only 95% CI excludes zero: 45.",
  "- Weighted only 95% CI excludes zero: 12.",
  "- Neither 95% CI excludes zero among comparable pairs: 178.",
  "",
  "## Reporting correction",
  "",
  "- SKIP specifications are no longer classified as Neither.",
  "- The correction affects reporting categories only.",
  "- Model estimates, standard errors, confidence intervals, P values, standardized effects, and fitted models are unchanged.",
  "",
  "## Figure output",
  "",
  "- Unicode labels Σ1/d and ΣHeat/d are retained.",
  "- PDF output uses cairo_pdf.",
  "- X11 and Cairo capabilities passed.",
  "- Six PNG and six PDF figures exist and are non-empty.",
  "",
  "**Any future extension should be implemented as a subsequent analysis stage rather than by modifying the frozen 07H07 outputs in place.**"
)

writeLines(
  MARKER,
  OUT_MARKER,
  useBytes = TRUE
)

cat(
  "\n============================================================\n",
  "07H07 FINAL FREEZE: PASS\n",
  "============================================================\n",
  "Freeze audit: ",
  OUT_AUDIT,
  "\n",
  "Freeze marker: ",
  OUT_MARKER,
  "\n",
  sep = ""
)
