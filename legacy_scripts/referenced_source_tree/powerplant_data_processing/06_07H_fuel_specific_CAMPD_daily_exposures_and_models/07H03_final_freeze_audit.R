#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999
)

suppressPackageStartupMessages({
  library(data.table)
})


# ============================================================
# 07H03 final freeze audit
#
# 只读07H02/07H03结果。
# 不重新拟合模型。
# 不修改07H01/07H02/07H03主脚本。
# ============================================================


stopf <- function(...) {
  stop(
    paste0(...),
    call. = FALSE
  )
}


ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

PP <- file.path(
  ROOT,
  "发电站数据处理"
)

H_ROOT <- file.path(
  PP,
  "06_07H_燃料特异CAMPD逐日暴露与模型"
)

H03_ROOT <- file.path(
  H_ROOT,
  "07H03_模型结果审计与汇总"
)

H03_SCRIPT <- file.path(
  H_ROOT,
  "07H03_燃料特异Count_Heat模型结果审计与汇总.R"
)

MASTER_FILE <- file.path(
  H03_ROOT,
  "01_07H03_全部模型结果主表.csv"
)

COL_FILE <- file.path(
  H03_ROOT,
  "05_07H03_Count_Heat共线性审计.csv"
)

HARD_FILE <- file.path(
  H03_ROOT,
  "11_07H03_hard_checks.csv"
)

FREEZE_AUDIT_FILE <- file.path(
  H03_ROOT,
  "07H03_final_freeze_audit.csv"
)

FREEZE_MARKER <- file.path(
  H03_ROOT,
  "07H03_FINAL_FREEZE.md"
)


# ============================================================
# Required files
# ============================================================

EXPECTED_OUTPUTS <- file.path(
  H03_ROOT,
  c(
    "01_07H03_全部模型结果主表.csv",
    "02_07H03_upwind主分析结果.csv",
    "03_07H03_downwind方向负对照结果.csv",
    "04_07H03_SKIP审计.csv",
    "05_07H03_Count_Heat共线性审计.csv",
    "06_07H03_模型覆盖情况.csv",
    "07_07H03_累计半径结果模式.csv",
    "08_07H03_upwind_downwind配对比较.csv",
    "09_07H03_quasi_filter配对比较.csv",
    "10_07H03_重点结果索引.csv",
    "11_07H03_hard_checks.csv",
    "12_07H03_模型结果审计与汇总报告.md",
    "13_sessionInfo.txt"
  )
)


if (!file.exists(H03_SCRIPT)) {
  stopf(
    "找不到冻结候选07H03主脚本：\n",
    H03_SCRIPT
  )
}


missing_outputs <- EXPECTED_OUTPUTS[
  !file.exists(EXPECTED_OUTPUTS)
]


if (length(missing_outputs) > 0L) {
  stopf(
    "07H03缺少正式输出：\n",
    paste(
      missing_outputs,
      collapse = "\n"
    )
  )
}


# ============================================================
# Read
# ============================================================

M <- fread(
  MASTER_FILE,
  showProgress = FALSE
)

C <- fread(
  COL_FILE,
  showProgress = FALSE
)

H <- fread(
  HARD_FILE,
  showProgress = FALSE
)


required_master <- c(
  "spec_id",
  "branch",
  "fuel",
  "direction",
  "radius_km",
  "model_type",
  "term",
  "status"
)

required_cor <- c(
  "branch",
  "fuel",
  "direction",
  "radius_km",
  "two_predictor_vif",
  "collinearity_flag_07H03"
)

required_hard <- c(
  "check",
  "observed",
  "expected",
  "status"
)


for (z in list(
  list(D = M, req = required_master, name = "07H03 master"),
  list(D = C, req = required_cor, name = "07H03 collinearity"),
  list(D = H, req = required_hard, name = "07H03 hard checks")
)) {

  miss <- setdiff(
    z$req,
    names(z$D)
  )

  if (length(miss) > 0L) {
    stopf(
      z$name,
      "缺少字段：",
      paste(
        miss,
        collapse = ", "
      )
    )
  }
}


# ============================================================
# Exact frozen dimensions
# ============================================================

N_EFFECT <- nrow(M)

N_SPECS <- uniqueN(
  M$spec_id
)

N_PASS <- M[
  status == "PASS",
  .N
]

N_SKIP <- M[
  status == "SKIP",
  .N
]

N_FAIL <- M[
  status == "FAIL",
  .N
]

N_COR <- nrow(C)

N_INF <- C[
  is.infinite(
    as.numeric(
      two_predictor_vif
    )
  ),
  .N
]

N_GE10 <- C[
  is.finite(
    as.numeric(
      two_predictor_vif
    )
  ) &
    as.numeric(
      two_predictor_vif
    ) >= 10,
  .N
]

N_5_10 <- C[
  is.finite(
    as.numeric(
      two_predictor_vif
    )
  ) &
    as.numeric(
      two_predictor_vif
    ) >= 5 &
    as.numeric(
      two_predictor_vif
    ) < 10,
  .N
]

N_H03_HARD_FAIL <- H[
  status != "PASS",
  .N
]

DUP_EFFECT_KEY <- M[
  ,
  .N,
  by = .(
    branch,
    fuel,
    direction,
    radius_km,
    model_type,
    term
  )
][
  N > 1L,
  .N
]

DUP_COR_KEY <- C[
  ,
  .N,
  by = .(
    branch,
    fuel,
    direction,
    radius_km
  )
][
  N > 1L,
  .N
]


AUDIT <- data.table(

  check = c(
    "07H03_main_script_exists",
    "07H03_formal_outputs_13",
    "07H03_hard_fail_zero",
    "master_effect_rows_960",
    "model_specs_720",
    "pass_effect_rows_846",
    "skip_effect_rows_114",
    "fail_effect_rows_zero",
    "duplicate_effect_keys_zero",
    "correlation_rows_240",
    "duplicate_correlation_keys_zero",
    "VIF_inf_6",
    "finite_VIF_ge10_14",
    "VIF_5_to_10_56"
  ),

  observed = c(
    as.character(file.exists(H03_SCRIPT)),
    as.character(sum(file.exists(EXPECTED_OUTPUTS))),
    as.character(N_H03_HARD_FAIL),
    as.character(N_EFFECT),
    as.character(N_SPECS),
    as.character(N_PASS),
    as.character(N_SKIP),
    as.character(N_FAIL),
    as.character(DUP_EFFECT_KEY),
    as.character(N_COR),
    as.character(DUP_COR_KEY),
    as.character(N_INF),
    as.character(N_GE10),
    as.character(N_5_10)
  ),

  expected = c(
    "TRUE",
    "13",
    "0",
    "960",
    "720",
    "846",
    "114",
    "0",
    "0",
    "240",
    "0",
    "6",
    "14",
    "56"
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
  FREEZE_AUDIT_FILE
)


print(
  AUDIT,
  nrows = Inf
)


N_FREEZE_FAIL <- AUDIT[
  status == "FAIL",
  .N
]


if (N_FREEZE_FAIL > 0L) {

  stopf(
    "07H03 final freeze audit存在",
    N_FREEZE_FAIL,
    "项FAIL。"
  )
}


FREEZE_TEXT <- c(

  "# 07H03 FINAL FREEZE",

  "",

  "- Status: PASS",

  paste0(
    "- Freeze time: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),

  paste0(
    "- Main script: ",
    H03_SCRIPT
  ),

  "- 07H03 main script is frozen and must not be modified for subsequent sensitivity analyses.",

  "- Effect rows: 960.",

  "- Model specs: 720.",

  "- PASS effect rows: 846.",

  "- SKIP effect rows: 114.",

  "- FAIL effect rows: 0.",

  "- Count–Heat correlation rows: 240.",

  "- Perfect collinearity VIF=Inf: 6.",

  "- Finite VIF>=10: 14.",

  "- 5<=VIF<10: 56.",

  "- All 07H03 hard checks passed.",

  "",

  "**07H03 FINAL FREEZE: PASS**"
)


writeLines(
  FREEZE_TEXT,
  FREEZE_MARKER,
  useBytes = TRUE
)


cat(
  "\n============================================================\n",
  "07H03 FINAL FREEZE: PASS\n",
  "============================================================\n",
  "冻结主脚本：\n",
  H03_SCRIPT,
  "\n\n冻结标记：\n",
  FREEZE_MARKER,
  "\n",
  sep = ""
)

