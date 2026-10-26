#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999
)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("缺少R包 data.table。", call. = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
})


# ============================================================
# 07H06
# 主分析与距离加权敏感性结果综合审计
#
# 只读：
#   07H03 frozen main results
#   07H05 distance-weighted sensitivity results
#
# 不重新拟合模型。
# 不修改07H01–07H05。
#
# 综合6类结果：
#   1 raw_count_only
#   2 raw_heat_only
#   3 joint_count
#   4 joint_heat
#   5 weighted_proximity
#   6 weighted_heat_proximity
#
# 不进行：
#   - significance-based radius selection
#   - best fuel ranking
#   - best exposure ranking
#   - upwind vs downwind formal coefficient-difference test
#
# upwind/downwind comparison is descriptive directional-control audit only.
# ============================================================


# ============================================================
# 1. Helpers
# ============================================================

stopf <- function(...) {
  stop(
    paste0(...),
    call. = FALSE
  )
}


safe_num <- function(x) {
  suppressWarnings(
    as.numeric(x)
  )
}


assert_file <- function(
  file,
  label
) {

  if (!file.exists(file)) {
    stopf(
      "找不到",
      label,
      "：\n",
      file
    )
  }
}


assert_cols <- function(
  D,
  required,
  label
) {

  miss <- setdiff(
    required,
    names(D)
  )

  if (length(miss) > 0L) {
    stopf(
      label,
      "缺少字段：",
      paste(
        miss,
        collapse = ", "
      )
    )
  }
}


assert_exact_set <- function(
  observed,
  expected,
  label
) {

  a <- sort(
    unique(
      as.character(
        observed
      )
    )
  )

  b <- sort(
    unique(
      as.character(
        expected
      )
    )
  )

  if (!identical(a, b)) {
    stopf(
      label,
      "取值集合错误。\nobserved=",
      paste(a, collapse = "|"),
      "\nexpected=",
      paste(b, collapse = "|")
    )
  }
}


ci_class <- function(
  status,
  low,
  high
) {

  out <- rep(
    NA_character_,
    length(status)
  )

  out[
    status == "SKIP"
  ] <- "SKIP"

  out[
    status == "FAIL"
  ] <- "FAIL"

  idx <- status == "PASS"

  out[
    idx &
      (
        !is.finite(low) |
        !is.finite(high)
      )
  ] <- "CI_NOT_AVAILABLE"

  out[
    idx &
      is.finite(low) &
      is.finite(high) &
      low > 0
  ] <- "CI_ABOVE_ZERO"

  out[
    idx &
      is.finite(low) &
      is.finite(high) &
      high < 0
  ] <- "CI_BELOW_ZERO"

  out[
    idx &
      is.finite(low) &
      is.finite(high) &
      low <= 0 &
      high >= 0
  ] <- "CI_INCLUDES_ZERO"

  out
}


effect_sign <- function(x) {

  fifelse(
    is.na(x),
    NA_character_,

    fifelse(
      x > 0,
      "positive",

      fifelse(
        x < 0,
        "negative",
        "zero"
      )
    )
  )
}


ci_excludes_zero <- function(
  status,
  low,
  high
) {

  status == "PASS" &
    is.finite(low) &
    is.finite(high) &
    (
      low > 0 |
      high < 0
    )
}


collapse_radii <- function(
  radius,
  keep
) {

  x <- sort(
    unique(
      radius[
        keep %in% TRUE &
        !is.na(radius)
      ]
    )
  )

  if (length(x) == 0L) {
    return("")
  }

  paste(
    x,
    collapse = ";"
  )
}


pair_ci_pattern <- function(
  status1,
  low1,
  high1,
  status2,
  low2,
  high2,
  label1,
  label2
) {

  n <- length(status1)

  out <- rep(
    NA_character_,
    n
  )

  p1 <- status1 == "PASS"

  p2 <- status2 == "PASS"

  c1 <- ci_excludes_zero(
    status1,
    low1,
    high1
  )

  c2 <- ci_excludes_zero(
    status2,
    low2,
    high2
  )

  out[
    !(p1 & p2)
  ] <- "PAIR_INCOMPLETE"

  out[
    p1 & p2 & c1 & c2
  ] <- "BOTH_CI_EXCLUDE_ZERO"

  out[
    p1 & p2 & c1 & !c2
  ] <- paste0(
    label1,
    "_ONLY_CI_EXCLUDES_ZERO"
  )

  out[
    p1 & p2 & !c1 & c2
  ] <- paste0(
    label2,
    "_ONLY_CI_EXCLUDES_ZERO"
  )

  out[
    p1 & p2 & !c1 & !c2
  ] <- "NEITHER_CI_EXCLUDES_ZERO"

  out
}


# ============================================================
# 2. Formal paths
# ============================================================

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

H05_ROOT <- file.path(
  H_ROOT,
  "07H05_model"
)


H03_MASTER_FILE <- file.path(
  H03_ROOT,
  "01_07H03_全部模型结果主表.csv"
)

H03_HARD_FILE <- file.path(
  H03_ROOT,
  "11_07H03_hard_checks.csv"
)


H05_EFFECT_FILE <- file.path(
  H05_ROOT,
  "01_07H05_weighted_model_effects_all.csv"
)

H05_COR_FILE <- file.path(
  H05_ROOT,
  "03_07H05_weighted_metric_collinearity.csv"
)

H05_COMPARE_FILE <- file.path(
  H05_ROOT,
  "05_07H05_vs_07H02_raw_sensitivity_comparison.csv"
)

H05_DISTANCE_FILE <- file.path(
  H05_ROOT,
  "07_07H05_distance_floor_audit.csv"
)

H05_HARD_FILE <- file.path(
  H05_ROOT,
  "08_07H05_hard_checks.csv"
)


OUT <- file.path(
  H_ROOT,
  "07H06_主分析与距离加权敏感性综合审计"
)


dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)


OUT_LONG <- file.path(
  OUT,
  "01_07H06_全部主分析与敏感性结果长表.csv"
)

OUT_RW <- file.path(
  OUT,
  "02_07H06_raw与distance_weighted稳健性比较.csv"
)

OUT_DIR <- file.path(
  OUT,
  "03_07H06_upwind_downwind全部结果配对.csv"
)

OUT_BRANCH <- file.path(
  OUT,
  "04_07H06_quasi_filter全部结果配对.csv"
)

OUT_RADIUS <- file.path(
  OUT,
  "05_07H06_累计半径结果模式.csv"
)

OUT_JOINT <- file.path(
  OUT,
  "06_07H06_joint_Count_Heat共线性解释表.csv"
)

OUT_WCOR <- file.path(
  OUT,
  "07_07H06_weighted指标共线性审计.csv"
)

OUT_COVERAGE <- file.path(
  OUT,
  "08_07H06_模型覆盖情况.csv"
)

OUT_UPWIND <- file.path(
  OUT,
  "09_07H06_upwind主结果综合摘要.csv"
)

OUT_HARD <- file.path(
  OUT,
  "10_07H06_hard_checks.csv"
)

OUT_REPORT <- file.path(
  OUT,
  "11_07H06_主分析与敏感性综合报告.md"
)

OUT_SESSION <- file.path(
  OUT,
  "12_sessionInfo.txt"
)


# ============================================================
# 3. Frozen design
# ============================================================

FUELS <- c(
  "coal",
  "gas",
  "oil",
  "other_fossil",
  "biomass",
  "solid_waste"
)

BRANCHES <- c(
  "filter_period",
  "quasi_daily"
)

DIRECTIONS <- c(
  "upwind",
  "downwind"
)

RADII <- c(
  20L,
  25L,
  30L,
  35L,
  40L,
  45L,
  50L,
  100L,
  150L,
  200L
)

COMPONENTS <- c(
  "raw_count_only",
  "raw_heat_only",
  "joint_count",
  "joint_heat",
  "weighted_proximity",
  "weighted_heat_proximity"
)


EXPECTED_H03_ROWS <- 960L

EXPECTED_H03_PASS <- 846L

EXPECTED_H03_SKIP <- 114L

EXPECTED_H03_FAIL <- 0L


EXPECTED_H05_ROWS <- 480L

EXPECTED_H05_PASS <- 426L

EXPECTED_H05_SKIP <- 54L

EXPECTED_H05_FAIL <- 0L


EXPECTED_LONG_ROWS <- 1440L

EXPECTED_RW_ROWS <- 480L

EXPECTED_DIR_PAIR_ROWS <- 720L

EXPECTED_BRANCH_PAIR_ROWS <- 720L

EXPECTED_RADIUS_GROUPS <- 144L

EXPECTED_JOINT_ROWS <- 480L

EXPECTED_WEIGHTED_COR_ROWS <- 240L


# ============================================================
# 4. Input existence
# ============================================================

assert_file(
  H03_MASTER_FILE,
  "07H03 master"
)

assert_file(
  H03_HARD_FILE,
  "07H03 hard checks"
)

assert_file(
  H05_EFFECT_FILE,
  "07H05 weighted effects"
)

assert_file(
  H05_COR_FILE,
  "07H05 weighted collinearity"
)

assert_file(
  H05_COMPARE_FILE,
  "07H05 raw-vs-weighted comparison"
)

assert_file(
  H05_DISTANCE_FILE,
  "07H05 distance-floor audit"
)

assert_file(
  H05_HARD_FILE,
  "07H05 hard checks"
)


# ============================================================
# 5. Read
# ============================================================

M03 <- fread(
  H03_MASTER_FILE,
  showProgress = FALSE
)

H03 <- fread(
  H03_HARD_FILE,
  showProgress = FALSE
)

M05 <- fread(
  H05_EFFECT_FILE,
  showProgress = FALSE
)

C05 <- fread(
  H05_COR_FILE,
  showProgress = FALSE
)

RW05 <- fread(
  H05_COMPARE_FILE,
  showProgress = FALSE
)

DIST05 <- fread(
  H05_DISTANCE_FILE,
  showProgress = FALSE
)

H05 <- fread(
  H05_HARD_FILE,
  showProgress = FALSE
)


# ============================================================
# 6. Schema audit
# ============================================================

REQ03 <- c(
  "branch",
  "fuel",
  "direction",
  "radius_km",
  "result_component",
  "status",
  "estimate_mBq_m3",
  "se_mBq_m3",
  "ci95_low_mBq_m3",
  "ci95_high_mBq_m3",
  "p_wald",
  "sd_raw_exposure",
  "estimate_per_1SD_mBq_m3",
  "raw_increment",
  "raw_unit",
  "pearson_r",
  "spearman_rho",
  "two_predictor_vif",
  "joint_interpretability",
  "note"
)


REQ05 <- c(
  "branch",
  "fuel",
  "direction",
  "radius_km",
  "metric",
  "status",
  "estimate_mBq_m3",
  "se_mBq_m3",
  "ci95_low_mBq_m3",
  "ci95_high_mBq_m3",
  "p_wald",
  "sd_raw_exposure",
  "estimate_per_1SD_mBq_m3",
  "ci95_low_per_1SD_mBq_m3",
  "ci95_high_per_1SD_mBq_m3",
  "raw_increment",
  "raw_unit",
  "singular",
  "convergence_message",
  "warnings",
  "note",
  "spec_id"
)


REQ_COR <- c(
  "branch",
  "fuel",
  "direction",
  "radius_km",
  "pearson_r",
  "spearman_rho",
  "two_predictor_vif",
  "collinearity_flag"
)


REQ_RW <- c(
  "branch",
  "fuel",
  "direction",
  "radius_km",
  "comparison_family",
  "raw_status",
  "weighted_status",
  "raw_estimate_mBq_m3",
  "weighted_estimate_mBq_m3",
  "raw_ci95_low_mBq_m3",
  "raw_ci95_high_mBq_m3",
  "weighted_ci95_low_mBq_m3",
  "weighted_ci95_high_mBq_m3",
  "raw_estimate_per_1SD_mBq_m3",
  "weighted_estimate_per_1SD_mBq_m3"
)


assert_cols(
  M03,
  REQ03,
  "07H03 master"
)

assert_cols(
  M05,
  REQ05,
  "07H05 weighted effects"
)

assert_cols(
  C05,
  REQ_COR,
  "07H05 weighted correlation"
)

assert_cols(
  RW05,
  REQ_RW,
  "07H05 raw-vs-weighted comparison"
)

assert_cols(
  H03,
  c(
    "check",
    "status"
  ),
  "07H03 hard checks"
)

assert_cols(
  H05,
  c(
    "check",
    "status"
  ),
  "07H05 hard checks"
)


# ============================================================
# 7. Frozen input hard checks
# ============================================================

N_H03_HARD_FAIL <- H03[
  status != "PASS",
  .N
]

N_H05_HARD_FAIL <- H05[
  status != "PASS",
  .N
]


if (
  N_H03_HARD_FAIL > 0L
) {

  print(
    H03[
      status != "PASS"
    ]
  )

  stopf(
    "07H03存在未PASS hard check。"
  )
}


if (
  N_H05_HARD_FAIL > 0L
) {

  print(
    H05[
      status != "PASS"
    ]
  )

  stopf(
    "07H05存在未PASS hard check。"
  )
}


# ============================================================
# 8. Exact input dimensions
# ============================================================

M03[
  ,
  radius_km :=
    as.integer(
      radius_km
    )
]

M05[
  ,
  radius_km :=
    as.integer(
      radius_km
    )
]

C05[
  ,
  radius_km :=
    as.integer(
      radius_km
    )
]

RW05[
  ,
  radius_km :=
    as.integer(
      radius_km
    )
]


N03_PASS <- M03[
  status == "PASS",
  .N
]

N03_SKIP <- M03[
  status == "SKIP",
  .N
]

N03_FAIL <- M03[
  status == "FAIL",
  .N
]


N05_PASS <- M05[
  status == "PASS",
  .N
]

N05_SKIP <- M05[
  status == "SKIP",
  .N
]

N05_FAIL <- M05[
  status == "FAIL",
  .N
]


if (
  nrow(M03) != EXPECTED_H03_ROWS ||
  N03_PASS != EXPECTED_H03_PASS ||
  N03_SKIP != EXPECTED_H03_SKIP ||
  N03_FAIL != EXPECTED_H03_FAIL
) {

  stopf(
    "07H03冻结结果维度不一致：rows=",
    nrow(M03),
    "; PASS=",
    N03_PASS,
    "; SKIP=",
    N03_SKIP,
    "; FAIL=",
    N03_FAIL
  )
}


if (
  nrow(M05) != EXPECTED_H05_ROWS ||
  N05_PASS != EXPECTED_H05_PASS ||
  N05_SKIP != EXPECTED_H05_SKIP ||
  N05_FAIL != EXPECTED_H05_FAIL
) {

  stopf(
    "07H05冻结结果维度不一致：rows=",
    nrow(M05),
    "; PASS=",
    N05_PASS,
    "; SKIP=",
    N05_SKIP,
    "; FAIL=",
    N05_FAIL
  )
}


if (
  nrow(C05) !=
    EXPECTED_WEIGHTED_COR_ROWS
) {

  stopf(
    "07H05 weighted correlation不是240行。"
  )
}


if (
  nrow(RW05) !=
    EXPECTED_RW_ROWS
) {

  stopf(
    "07H05 raw-vs-weighted comparison不是480行。"
  )
}


assert_exact_set(
  M03$branch,
  BRANCHES,
  "07H03 branch"
)

assert_exact_set(
  M03$fuel,
  FUELS,
  "07H03 fuel"
)

assert_exact_set(
  M03$direction,
  DIRECTIONS,
  "07H03 direction"
)

assert_exact_set(
  M03$radius_km,
  RADII,
  "07H03 radius"
)


assert_exact_set(
  M05$branch,
  BRANCHES,
  "07H05 branch"
)

assert_exact_set(
  M05$fuel,
  FUELS,
  "07H05 fuel"
)

assert_exact_set(
  M05$direction,
  DIRECTIONS,
  "07H05 direction"
)

assert_exact_set(
  M05$radius_km,
  RADII,
  "07H05 radius"
)


# ============================================================
# 9. Convert 07H03 to unified long table
# ============================================================

MAIN_COMPONENT <- fcase(

  M03$result_component ==
    "count_only",
  "raw_count_only",

  M03$result_component ==
    "heat_only",
  "raw_heat_only",

  M03$result_component ==
    "joint_count",
  "joint_count",

  M03$result_component ==
    "joint_heat",
  "joint_heat",

  default =
    NA_character_
)


if (
  any(
    is.na(
      MAIN_COMPONENT
    )
  )
) {

  stopf(
    "07H03存在无法识别的result_component。"
  )
}


L03 <- data.table(

  analysis_layer =
    "main_raw",

  branch =
    M03$branch,

  fuel =
    M03$fuel,

  direction =
    M03$direction,

  radius_km =
    M03$radius_km,

  component =
    MAIN_COMPONENT,

  source_component =
    M03$result_component,

  status =
    M03$status,

  raw_increment =
    safe_num(
      M03$raw_increment
    ),

  raw_unit =
    as.character(
      M03$raw_unit
    ),

  estimate_mBq_m3 =
    safe_num(
      M03$estimate_mBq_m3
    ),

  se_mBq_m3 =
    safe_num(
      M03$se_mBq_m3
    ),

  ci95_low_mBq_m3 =
    safe_num(
      M03$ci95_low_mBq_m3
    ),

  ci95_high_mBq_m3 =
    safe_num(
      M03$ci95_high_mBq_m3
    ),

  p_wald =
    safe_num(
      M03$p_wald
    ),

  sd_raw_exposure =
    safe_num(
      M03$sd_raw_exposure
    ),

  estimate_per_1SD_mBq_m3 =
    safe_num(
      M03$estimate_per_1SD_mBq_m3
    ),

  ci95_low_per_1SD_mBq_m3 =
    NA_real_,

  ci95_high_per_1SD_mBq_m3 =
    NA_real_,

  count_heat_pearson_r =
    safe_num(
      M03$pearson_r
    ),

  count_heat_spearman_rho =
    safe_num(
      M03$spearman_rho
    ),

  count_heat_vif =
    safe_num(
      M03$two_predictor_vif
    ),

  interpretability_flag =
    as.character(
      M03$joint_interpretability
    ),

  singular =
    NA,

  convergence_message =
    "",

  warnings =
    "",

  note =
    as.character(
      M03$note
    )
)


L03[
  ,
  ci_status :=
    ci_class(
      status,
      ci95_low_mBq_m3,
      ci95_high_mBq_m3
    )
]


L03[
  ,
  effect_sign :=
    effect_sign(
      estimate_mBq_m3
    )
]


# ============================================================
# 10. Convert 07H05 to unified long table
# ============================================================

WEIGHTED_COMPONENT <- fifelse(

  M05$metric ==
    "proximity_invkm",

  "weighted_proximity",

  fifelse(
    M05$metric ==
      "heat_proximity_mmbtu_per_km",

    "weighted_heat_proximity",

    NA_character_
  )
)


if (
  any(
    is.na(
      WEIGHTED_COMPONENT
    )
  )
) {

  stopf(
    "07H05存在无法识别的metric。"
  )
}


L05 <- data.table(

  analysis_layer =
    "distance_weighted_sensitivity",

  branch =
    M05$branch,

  fuel =
    M05$fuel,

  direction =
    M05$direction,

  radius_km =
    M05$radius_km,

  component =
    WEIGHTED_COMPONENT,

  source_component =
    M05$metric,

  status =
    M05$status,

  raw_increment =
    safe_num(
      M05$raw_increment
    ),

  raw_unit =
    as.character(
      M05$raw_unit
    ),

  estimate_mBq_m3 =
    safe_num(
      M05$estimate_mBq_m3
    ),

  se_mBq_m3 =
    safe_num(
      M05$se_mBq_m3
    ),

  ci95_low_mBq_m3 =
    safe_num(
      M05$ci95_low_mBq_m3
    ),

  ci95_high_mBq_m3 =
    safe_num(
      M05$ci95_high_mBq_m3
    ),

  p_wald =
    safe_num(
      M05$p_wald
    ),

  sd_raw_exposure =
    safe_num(
      M05$sd_raw_exposure
    ),

  estimate_per_1SD_mBq_m3 =
    safe_num(
      M05$estimate_per_1SD_mBq_m3
    ),

  ci95_low_per_1SD_mBq_m3 =
    safe_num(
      M05$ci95_low_per_1SD_mBq_m3
    ),

  ci95_high_per_1SD_mBq_m3 =
    safe_num(
      M05$ci95_high_per_1SD_mBq_m3
    ),

  count_heat_pearson_r =
    NA_real_,

  count_heat_spearman_rho =
    NA_real_,

  count_heat_vif =
    NA_real_,

  interpretability_flag =
    "SINGLE_WEIGHTED_MODEL",

  singular =
    M05$singular,

  convergence_message =
    as.character(
      M05$convergence_message
    ),

  warnings =
    as.character(
      M05$warnings
    ),

  note =
    as.character(
      M05$note
    )
)


L05[
  ,
  ci_status :=
    ci_class(
      status,
      ci95_low_mBq_m3,
      ci95_high_mBq_m3
    )
]


L05[
  ,
  effect_sign :=
    effect_sign(
      estimate_mBq_m3
    )
]


# ============================================================
# 11. Unified 1440-row result table
# ============================================================

LONG <- rbindlist(
  list(
    L03,
    L05
  ),
  use.names = TRUE,
  fill = TRUE
)


assert_exact_set(
  LONG$component,
  COMPONENTS,
  "07H06 component"
)


if (
  nrow(LONG) !=
    EXPECTED_LONG_ROWS
) {

  stopf(
    "07H06 unified long table不是1440行；实际=",
    nrow(LONG)
  )
}


N_LONG_DUP <- LONG[
  ,
  .N,
  by = .(
    branch,
    fuel,
    direction,
    radius_km,
    component
  )
][
  N > 1L,
  .N
]


if (
  N_LONG_DUP > 0L
) {

  stopf(
    "统一结果长表存在重复主键：",
    N_LONG_DUP
  )
}


setorder(
  LONG,
  branch,
  fuel,
  direction,
  radius_km,
  component
)


fwrite(
  LONG,
  OUT_LONG
)


# ============================================================
# 12. Raw vs weighted robustness
#
# Count-only  <-> Σ1/d
# Heat-only   <-> ΣHeat/d
# ============================================================

RW <- copy(
  RW05
)


NUM_RW <- c(
  "raw_estimate_mBq_m3",
  "weighted_estimate_mBq_m3",
  "raw_ci95_low_mBq_m3",
  "raw_ci95_high_mBq_m3",
  "weighted_ci95_low_mBq_m3",
  "weighted_ci95_high_mBq_m3",
  "raw_estimate_per_1SD_mBq_m3",
  "weighted_estimate_per_1SD_mBq_m3"
)


for (
  cc in NUM_RW
) {

  RW[
    ,
    (cc) :=
      safe_num(
        get(cc)
      )
  ]
}


RW[
  ,
  `:=`(

    raw_effect_sign =
      effect_sign(
        raw_estimate_mBq_m3
      ),

    weighted_effect_sign =
      effect_sign(
        weighted_estimate_mBq_m3
      ),

    raw_1SD_sign =
      effect_sign(
        raw_estimate_per_1SD_mBq_m3
      ),

    weighted_1SD_sign =
      effect_sign(
        weighted_estimate_per_1SD_mBq_m3
      ),

    same_sign_fixed_increment =
      fifelse(
        raw_status == "PASS" &
        weighted_status == "PASS" &
        is.finite(
          raw_estimate_mBq_m3
        ) &
        is.finite(
          weighted_estimate_mBq_m3
        ),

        sign(
          raw_estimate_mBq_m3
        ) ==
          sign(
            weighted_estimate_mBq_m3
          ),

        NA
      ),

    same_sign_1SD =
      fifelse(
        raw_status == "PASS" &
        weighted_status == "PASS" &
        is.finite(
          raw_estimate_per_1SD_mBq_m3
        ) &
        is.finite(
          weighted_estimate_per_1SD_mBq_m3
        ),

        sign(
          raw_estimate_per_1SD_mBq_m3
        ) ==
          sign(
            weighted_estimate_per_1SD_mBq_m3
          ),

        NA
      ),

    raw_ci_excludes_zero =
      ci_excludes_zero(
        raw_status,
        raw_ci95_low_mBq_m3,
        raw_ci95_high_mBq_m3
      ),

    weighted_ci_excludes_zero =
      ci_excludes_zero(
        weighted_status,
        weighted_ci95_low_mBq_m3,
        weighted_ci95_high_mBq_m3
      ),

    oneSD_effect_difference_mBq_m3 =
      weighted_estimate_per_1SD_mBq_m3 -
      raw_estimate_per_1SD_mBq_m3,

    oneSD_absolute_difference_mBq_m3 =
      abs(
        weighted_estimate_per_1SD_mBq_m3 -
        raw_estimate_per_1SD_mBq_m3
      )
  )
]


RW[
  ,
  robustness_pattern :=
    fcase(

      raw_status != "PASS" |
        weighted_status != "PASS",
      "PAIR_INCOMPLETE",

      raw_ci_excludes_zero &
        weighted_ci_excludes_zero &
        same_sign_fixed_increment %in% TRUE,
      "BOTH_CI_EXCLUDE_ZERO_SAME_SIGN",

      raw_ci_excludes_zero &
        weighted_ci_excludes_zero &
        same_sign_fixed_increment %in% FALSE,
      "BOTH_CI_EXCLUDE_ZERO_OPPOSITE_SIGN",

      raw_ci_excludes_zero &
        !weighted_ci_excludes_zero,
      "RAW_ONLY_CI_EXCLUDES_ZERO",

      !raw_ci_excludes_zero &
        weighted_ci_excludes_zero,
      "WEIGHTED_ONLY_CI_EXCLUDES_ZERO",

      !raw_ci_excludes_zero &
        !weighted_ci_excludes_zero &
        same_sign_fixed_increment %in% TRUE,
      "NEITHER_CI_EXCLUDES_ZERO_SAME_SIGN",

      !raw_ci_excludes_zero &
        !weighted_ci_excludes_zero &
        same_sign_fixed_increment %in% FALSE,
      "NEITHER_CI_EXCLUDES_ZERO_OPPOSITE_SIGN",

      default =
        "OTHER"
    )
]


N_RW_DUP <- RW[
  ,
  .N,
  by = .(
    branch,
    fuel,
    direction,
    radius_km,
    comparison_family
  )
][
  N > 1L,
  .N
]


if (
  N_RW_DUP > 0L
) {

  stopf(
    "raw-vs-weighted比较存在重复主键。"
  )
}


setorder(
  RW,
  branch,
  fuel,
  direction,
  radius_km,
  comparison_family
)


fwrite(
  RW,
  OUT_RW
)


# ============================================================
# 13. Upwind vs downwind comparison for ALL 6 components
#
# Descriptive negative-control comparison only.
# This is NOT a formal test of coefficient difference.
# ============================================================

UP <- LONG[
  direction == "upwind"
]


DOWN <- LONG[
  direction == "downwind"
]


UP2 <- UP[
  ,
  .(
    branch,
    fuel,
    radius_km,
    component,
    analysis_layer,

    upwind_status =
      status,

    upwind_estimate =
      estimate_mBq_m3,

    upwind_ci_low =
      ci95_low_mBq_m3,

    upwind_ci_high =
      ci95_high_mBq_m3,

    upwind_p =
      p_wald,

    upwind_1SD =
      estimate_per_1SD_mBq_m3,

    upwind_ci_status =
      ci_status,

    upwind_interpretability =
      interpretability_flag
  )
]


DOWN2 <- DOWN[
  ,
  .(
    branch,
    fuel,
    radius_km,
    component,
    analysis_layer,

    downwind_status =
      status,

    downwind_estimate =
      estimate_mBq_m3,

    downwind_ci_low =
      ci95_low_mBq_m3,

    downwind_ci_high =
      ci95_high_mBq_m3,

    downwind_p =
      p_wald,

    downwind_1SD =
      estimate_per_1SD_mBq_m3,

    downwind_ci_status =
      ci_status,

    downwind_interpretability =
      interpretability_flag
  )
]


DIRPAIR <- merge(
  UP2,
  DOWN2,

  by = c(
    "branch",
    "fuel",
    "radius_km",
    "component",
    "analysis_layer"
  ),

  all = TRUE,
  sort = FALSE
)


DIRPAIR[
  ,
  `:=`(

    both_estimable =
      upwind_status == "PASS" &
      downwind_status == "PASS",

    same_sign =
      fifelse(
        upwind_status == "PASS" &
        downwind_status == "PASS" &
        is.finite(
          upwind_estimate
        ) &
        is.finite(
          downwind_estimate
        ),

        sign(
          upwind_estimate
        ) ==
          sign(
            downwind_estimate
          ),

        NA
      ),

    same_sign_1SD =
      fifelse(
        upwind_status == "PASS" &
        downwind_status == "PASS" &
        is.finite(
          upwind_1SD
        ) &
        is.finite(
          downwind_1SD
        ),

        sign(
          upwind_1SD
        ) ==
          sign(
            downwind_1SD
          ),

        NA
      ),

    ci_pattern =
      pair_ci_pattern(
        upwind_status,
        upwind_ci_low,
        upwind_ci_high,

        downwind_status,
        downwind_ci_low,
        downwind_ci_high,

        "UPWIND",
        "DOWNWIND"
      )
  )
]


setorder(
  DIRPAIR,
  branch,
  fuel,
  radius_km,
  component
)


fwrite(
  DIRPAIR,
  OUT_DIR
)


# ============================================================
# 14. Quasi-daily vs filter-period comparison
# ============================================================

Q <- LONG[
  branch ==
    "quasi_daily"
]


F <- LONG[
  branch ==
    "filter_period"
]


Q2 <- Q[
  ,
  .(
    fuel,
    direction,
    radius_km,
    component,
    analysis_layer,

    quasi_status =
      status,

    quasi_estimate =
      estimate_mBq_m3,

    quasi_ci_low =
      ci95_low_mBq_m3,

    quasi_ci_high =
      ci95_high_mBq_m3,

    quasi_p =
      p_wald,

    quasi_1SD =
      estimate_per_1SD_mBq_m3,

    quasi_ci_status =
      ci_status,

    quasi_interpretability =
      interpretability_flag
  )
]


F2 <- F[
  ,
  .(
    fuel,
    direction,
    radius_km,
    component,
    analysis_layer,

    filter_status =
      status,

    filter_estimate =
      estimate_mBq_m3,

    filter_ci_low =
      ci95_low_mBq_m3,

    filter_ci_high =
      ci95_high_mBq_m3,

    filter_p =
      p_wald,

    filter_1SD =
      estimate_per_1SD_mBq_m3,

    filter_ci_status =
      ci_status,

    filter_interpretability =
      interpretability_flag
  )
]


BRPAIR <- merge(
  Q2,
  F2,

  by = c(
    "fuel",
    "direction",
    "radius_km",
    "component",
    "analysis_layer"
  ),

  all = TRUE,
  sort = FALSE
)


BRPAIR[
  ,
  `:=`(

    both_estimable =
      quasi_status == "PASS" &
      filter_status == "PASS",

    same_sign =
      fifelse(
        quasi_status == "PASS" &
        filter_status == "PASS" &
        is.finite(
          quasi_estimate
        ) &
        is.finite(
          filter_estimate
        ),

        sign(
          quasi_estimate
        ) ==
          sign(
            filter_estimate
          ),

        NA
      ),

    same_sign_1SD =
      fifelse(
        quasi_status == "PASS" &
        filter_status == "PASS" &
        is.finite(
          quasi_1SD
        ) &
        is.finite(
          filter_1SD
        ),

        sign(
          quasi_1SD
        ) ==
          sign(
            filter_1SD
          ),

        NA
      ),

    ci_pattern =
      pair_ci_pattern(
        quasi_status,
        quasi_ci_low,
        quasi_ci_high,

        filter_status,
        filter_ci_low,
        filter_ci_high,

        "QUASI",
        "FILTER"
      )
  )
]


setorder(
  BRPAIR,
  fuel,
  direction,
  radius_km,
  component
)


fwrite(
  BRPAIR,
  OUT_BRANCH
)


# ============================================================
# 15. Radius-pattern summary
#
# All 10 radii retained.
# No best-radius selection.
# ============================================================

RADIUS <- LONG[
  ,
  .(

    n_radii =
      uniqueN(
        radius_km
      ),

    n_pass =
      sum(
        status == "PASS"
      ),

    n_skip =
      sum(
        status == "SKIP"
      ),

    n_fail =
      sum(
        status == "FAIL"
      ),

    n_positive =
      sum(
        status == "PASS" &
        estimate_mBq_m3 > 0,
        na.rm = TRUE
      ),

    n_negative =
      sum(
        status == "PASS" &
        estimate_mBq_m3 < 0,
        na.rm = TRUE
      ),

    n_ci_above_zero =
      sum(
        ci_status ==
          "CI_ABOVE_ZERO",
        na.rm = TRUE
      ),

    n_ci_below_zero =
      sum(
        ci_status ==
          "CI_BELOW_ZERO",
        na.rm = TRUE
      ),

    n_ci_include_zero =
      sum(
        ci_status ==
          "CI_INCLUDES_ZERO",
        na.rm = TRUE
      ),

    positive_radii_km =
      collapse_radii(
        radius_km,
        status == "PASS" &
        estimate_mBq_m3 > 0
      ),

    negative_radii_km =
      collapse_radii(
        radius_km,
        status == "PASS" &
        estimate_mBq_m3 < 0
      ),

    ci_above_zero_radii_km =
      collapse_radii(
        radius_km,
        ci_status ==
          "CI_ABOVE_ZERO"
      ),

    ci_below_zero_radii_km =
      collapse_radii(
        radius_km,
        ci_status ==
          "CI_BELOW_ZERO"
      ),

    all_estimable_same_sign = {

      z <- sign(
        estimate_mBq_m3[
          status == "PASS" &
          is.finite(
            estimate_mBq_m3
          )
        ]
      )

      z <- z[
        z != 0
      ]

      if (
        length(z) == 0L
      ) {
        NA
      } else {
        uniqueN(z) == 1L
      }
    }
  ),

  by = .(
    analysis_layer,
    branch,
    fuel,
    direction,
    component
  )
]


setorder(
  RADIUS,
  branch,
  fuel,
  direction,
  component
)


fwrite(
  RADIUS,
  OUT_RADIUS
)


# ============================================================
# 16. Main joint Count/Heat interpretation table
# ============================================================

JOINT <- LONG[
  component %in%
    c(
      "joint_count",
      "joint_heat"
    )
]


if (
  nrow(JOINT) !=
    EXPECTED_JOINT_ROWS
) {

  stopf(
    "joint结果不是480行；实际=",
    nrow(JOINT)
  )
}


JOINT[
  ,
  joint_interpretation_note :=
    fifelse(

      is.infinite(
        count_heat_vif
      ),

      "完全共线：joint条件效应不应作独立效应解释",

      fifelse(

        is.finite(
          count_heat_vif
        ) &
        count_heat_vif >= 10,

        "高共线：joint条件效应需高度谨慎解释",

        fifelse(

          is.finite(
            count_heat_vif
          ) &
          count_heat_vif >= 5,

          "中度共线：joint条件效应需谨慎解释",

          fifelse(
            is.finite(
              count_heat_vif
            ),

            "未触发VIF>=5共线性警示",

            "VIF不可估计"
          )
        )
      )
    )
]


setorder(
  JOINT,
  branch,
  fuel,
  direction,
  radius_km,
  component
)


fwrite(
  JOINT,
  OUT_JOINT
)


# ============================================================
# 17. Weighted-metric collinearity
# ============================================================

WCOR <- copy(
  C05
)


setorder(
  WCOR,
  branch,
  fuel,
  direction,
  radius_km
)


fwrite(
  WCOR,
  OUT_WCOR
)


N_W_INF <- WCOR[
  is.infinite(
    two_predictor_vif
  ),
  .N
]


N_W_GE10 <- WCOR[
  is.finite(
    two_predictor_vif
  ) &
  two_predictor_vif >= 10,
  .N
]


N_W_5_10 <- WCOR[
  is.finite(
    two_predictor_vif
  ) &
  two_predictor_vif >= 5 &
  two_predictor_vif < 10,
  .N
]


# ============================================================
# 18. Model coverage
# ============================================================

COVERAGE <- LONG[
  ,
  .(

    n_rows =
      .N,

    n_pass =
      sum(
        status == "PASS"
      ),

    n_skip =
      sum(
        status == "SKIP"
      ),

    n_fail =
      sum(
        status == "FAIL"
      ),

    n_ci_above_zero =
      sum(
        ci_status ==
          "CI_ABOVE_ZERO",
        na.rm = TRUE
      ),

    n_ci_below_zero =
      sum(
        ci_status ==
          "CI_BELOW_ZERO",
        na.rm = TRUE
      ),

    n_ci_include_zero =
      sum(
        ci_status ==
          "CI_INCLUDES_ZERO",
        na.rm = TRUE
      )
  ),

  by = .(
    analysis_layer,
    branch,
    fuel,
    direction,
    component
  )
]


setorder(
  COVERAGE,
  branch,
  fuel,
  direction,
  component
)


fwrite(
  COVERAGE,
  OUT_COVERAGE
)


# ============================================================
# 19. Upwind integrated summary
#
# This is the paper-oriented index.
# Still preserves ALL fuels and ALL radii.
# ============================================================

UP_SUMMARY <- LONG[
  direction ==
    "upwind",
  .(

    analysis_layer,
    branch,
    fuel,
    radius_km,
    component,

    status,

    raw_increment,
    raw_unit,

    estimate_mBq_m3,
    ci95_low_mBq_m3,
    ci95_high_mBq_m3,
    p_wald,

    sd_raw_exposure,
    estimate_per_1SD_mBq_m3,

    ci_status,
    effect_sign,

    count_heat_vif,
    interpretability_flag,

    note
  )
]


setorder(
  UP_SUMMARY,
  fuel,
  branch,
  radius_km,
  component
)


fwrite(
  UP_SUMMARY,
  OUT_UPWIND
)


# ============================================================
# 20. Final dimension checks
# ============================================================

N_DIR_DUP <- DIRPAIR[
  ,
  .N,
  by = .(
    branch,
    fuel,
    radius_km,
    component
  )
][
  N > 1L,
  .N
]


N_BRANCH_DUP <- BRPAIR[
  ,
  .N,
  by = .(
    fuel,
    direction,
    radius_km,
    component
  )
][
  N > 1L,
  .N
]


N_RADIUS_DUP <- RADIUS[
  ,
  .N,
  by = .(
    branch,
    fuel,
    direction,
    component
  )
][
  N > 1L,
  .N
]


N_PAIR_INCOMPLETE_RW <- RW[
  raw_status != "PASS" |
  weighted_status != "PASS",
  .N
]


N_RW_SAME_SIGN <- RW[
  raw_status == "PASS" &
  weighted_status == "PASS" &
  same_sign_fixed_increment %in% TRUE,
  .N
]


N_RW_OPPOSITE_SIGN <- RW[
  raw_status == "PASS" &
  weighted_status == "PASS" &
  same_sign_fixed_increment %in% FALSE,
  .N
]


# ============================================================
# 21. Hard checks
# ============================================================

HARD <- data.table(

  check = c(

    "07H03_hard_checks_all_pass",

    "07H05_hard_checks_all_pass",

    "07H03_rows_960",

    "07H03_PASS_846",

    "07H03_SKIP_114",

    "07H03_FAIL_zero",

    "07H05_rows_480",

    "07H05_PASS_426",

    "07H05_SKIP_54",

    "07H05_FAIL_zero",

    "unified_long_rows_1440",

    "unified_long_duplicate_keys_zero",

    "six_components_present",

    "raw_weighted_rows_480",

    "raw_weighted_duplicate_keys_zero",

    "direction_pair_rows_720",

    "direction_pair_duplicate_keys_zero",

    "branch_pair_rows_720",

    "branch_pair_duplicate_keys_zero",

    "radius_pattern_groups_144",

    "radius_pattern_duplicate_keys_zero",

    "joint_rows_480",

    "weighted_correlation_rows_240",

    "weighted_VIF_inf_zero",

    "weighted_finite_VIF_ge10_26",

    "weighted_VIF_5_to_10_55"
  ),

  observed = c(

    as.character(
      N_H03_HARD_FAIL
    ),

    as.character(
      N_H05_HARD_FAIL
    ),

    as.character(
      nrow(M03)
    ),

    as.character(
      N03_PASS
    ),

    as.character(
      N03_SKIP
    ),

    as.character(
      N03_FAIL
    ),

    as.character(
      nrow(M05)
    ),

    as.character(
      N05_PASS
    ),

    as.character(
      N05_SKIP
    ),

    as.character(
      N05_FAIL
    ),

    as.character(
      nrow(LONG)
    ),

    as.character(
      N_LONG_DUP
    ),

    as.character(
      uniqueN(
        LONG$component
      )
    ),

    as.character(
      nrow(RW)
    ),

    as.character(
      N_RW_DUP
    ),

    as.character(
      nrow(DIRPAIR)
    ),

    as.character(
      N_DIR_DUP
    ),

    as.character(
      nrow(BRPAIR)
    ),

    as.character(
      N_BRANCH_DUP
    ),

    as.character(
      nrow(RADIUS)
    ),

    as.character(
      N_RADIUS_DUP
    ),

    as.character(
      nrow(JOINT)
    ),

    as.character(
      nrow(WCOR)
    ),

    as.character(
      N_W_INF
    ),

    as.character(
      N_W_GE10
    ),

    as.character(
      N_W_5_10
    )
  ),

  expected = c(

    "0",
    "0",

    "960",
    "846",
    "114",
    "0",

    "480",
    "426",
    "54",
    "0",

    "1440",
    "0",
    "6",

    "480",
    "0",

    "720",
    "0",

    "720",
    "0",

    "144",
    "0",

    "480",
    "240",

    "0",
    "26",
    "55"
  )
)


HARD[
  ,
  status :=
    fifelse(
      observed == expected,
      "PASS",
      "FAIL"
    )
]


fwrite(
  HARD,
  OUT_HARD
)


N_HARD_FAIL <- HARD[
  status == "FAIL",
  .N
]


# ============================================================
# 22. Console summaries
# ============================================================

cat(
  "\n============================================================\n",
  "A. ALL COMPONENT STATUS\n",
  "============================================================\n",
  sep = ""
)


print(
  LONG[
    ,
    .(
      PASS =
        sum(
          status == "PASS"
        ),

      SKIP =
        sum(
          status == "SKIP"
        ),

      FAIL =
        sum(
          status == "FAIL"
        )
    ),
    by = .(
      analysis_layer,
      component
    )
  ],
  nrows = Inf
)


cat(
  "\n============================================================\n",
  "B. RAW vs DISTANCE-WEIGHTED ROBUSTNESS\n",
  "============================================================\n",
  "Both PASS + same sign = ",
  N_RW_SAME_SIGN,
  "\n",
  "Both PASS + opposite sign = ",
  N_RW_OPPOSITE_SIGN,
  "\n",
  "Pair incomplete = ",
  N_PAIR_INCOMPLETE_RW,
  "\n",
  sep = ""
)


print(
  RW[
    ,
    .N,
    by = .(
      comparison_family,
      robustness_pattern
    )
  ][
    order(
      comparison_family,
      robustness_pattern
    )
  ],
  nrows = Inf
)


cat(
  "\n============================================================\n",
  "C. UPWIND vs DOWNWIND DESCRIPTIVE PATTERN\n",
  "============================================================\n",
  sep = ""
)


print(
  DIRPAIR[
    ,
    .N,
    by = .(
      component,
      ci_pattern
    )
  ][
    order(
      component,
      ci_pattern
    )
  ],
  nrows = Inf
)


cat(
  "\n============================================================\n",
  "D. QUASI vs FILTER DESCRIPTIVE PATTERN\n",
  "============================================================\n",
  sep = ""
)


print(
  BRPAIR[
    ,
    .N,
    by = .(
      component,
      ci_pattern
    )
  ][
    order(
      component,
      ci_pattern
    )
  ],
  nrows = Inf
)


cat(
  "\n============================================================\n",
  "E. WEIGHTED COLLINEARITY\n",
  "============================================================\n",
  "VIF=Inf = ",
  N_W_INF,
  "\n",
  "Finite VIF>=10 = ",
  N_W_GE10,
  "\n",
  "5<=VIF<10 = ",
  N_W_5_10,
  "\n",
  sep = ""
)


cat(
  "\n============================================================\n",
  "F. HARD CHECKS\n",
  "============================================================\n",
  sep = ""
)


print(
  HARD,
  nrows = Inf
)


# ============================================================
# 23. Report
# ============================================================

REPORT <- c(

  "# 07H06 主分析与距离加权敏感性综合审计",

  "",

  paste0(
    "- Generated: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),

  "",

  "## Analysis layers",

  "",

  "- Main raw exposure results: Count-only, Heat-only, Joint Count, Joint Heat.",

  "- Spatial-definition sensitivity results: inverse-distance plant proximity and inverse-distance Heat proximity.",

  "- Distance-weighted results remain sensitivity analyses and do not replace the frozen raw Count/Heat primary analyses.",

  "",

  "## Frozen model status",

  "",

  "- 07H03: 960 effect rows; 846 PASS; 114 SKIP; 0 FAIL.",

  "- 07H05: 480 weighted sensitivity models; 426 PASS; 54 SKIP; 0 FAIL.",

  "- 07H05 PASS fits had no singular fit, convergence message, or warning in the completed model audit.",

  "",

  "## Raw versus distance-weighted robustness",

  "",

  "- Count-only is paired with Σ1/d.",

  "- Heat-only is paired with ΣHeat/d.",

  "- All prespecified fuels, directions, branches, and cumulative radii are retained.",

  "- One-SD effects are included because raw Count/Heat and distance-weighted metrics have different physical units.",

  "- Agreement is described by direction, confidence-interval pattern, and standardized effect pattern.",

  "- No exposure definition is selected according to P value.",

  "",

  "## Directional-control interpretation",

  "",

  "- Upwind is the prespecified main direction.",

  "- Downwind is treated as a directional negative-control comparison.",

  "- Separate upwind/downwind confidence intervals do not constitute a formal statistical test of the difference between the two coefficients.",

  "",

  "## Time-scale interpretation",

  "",

  "- Quasi-daily and filter-period branches are compared descriptively for effect-direction and confidence-interval consistency.",

  "- Neither branch is selected according to statistical significance.",

  "",

  "## Radius interpretation",

  "",

  "- All 10 cumulative radii are retained.",

  "- Radius patterns are summarized across 20, 25, 30, 35, 40, 45, 50, 100, 150, and 200 km.",

  "- No best radius is selected.",

  "",

  "## Joint Count + Heat",

  "",

  "- Joint Count and Joint Heat remain part of the frozen main analysis.",

  "- Their interpretation is conditional on the other raw exposure component.",

  "- High or infinite Count–Heat VIF combinations are explicitly flagged and require cautious interpretation.",

  "",

  "## Weighted joint model decision support",

  "",

  paste0(
    "- Weighted Σ1/d versus ΣHeat/d VIF=Inf combinations: ",
    N_W_INF,
    "."
  ),

  paste0(
    "- Finite weighted VIF>=10 combinations: ",
    N_W_GE10,
    "."
  ),

  paste0(
    "- Weighted 5<=VIF<10 combinations: ",
    N_W_5_10,
    "."
  ),

  "- No weighted joint model is added in 07H06.",

  "",

  "## Distance-floor audit",

  "",

  paste0(
    "- 0–200 km plant-monitor links: ",
    DIST05$n_links_0_200km[1],
    "."
  ),

  paste0(
    "- Links with original distance <1 km: ",
    DIST05$n_links_distance_lt1km[1],
    " (",
    sprintf(
      "%.4f",
      DIST05$pct_links_distance_lt1km[1]
    ),
    "%)."
  ),

  "- The 1-km floor affects only the inverse-distance denominator; cumulative-radius membership uses the original frozen 07A distance.",

  "",

  if (
    N_HARD_FAIL == 0L
  ) {

    "**07H06 INTEGRATED RESULT AUDIT: PASS**"

  } else {

    "**07H06 INTEGRATED RESULT AUDIT: FAIL**"
  }
)


writeLines(
  REPORT,
  OUT_REPORT,
  useBytes = TRUE
)


capture.output(
  sessionInfo(),
  file = OUT_SESSION
)


# ============================================================
# 24. Finish
# ============================================================

cat(
  "\n============================================================\n"
)


if (
  N_HARD_FAIL == 0L
) {

  cat(
    "07H06 INTEGRATED RESULT AUDIT: PASS\n"
  )

} else {

  cat(
    "07H06 INTEGRATED RESULT AUDIT: FAIL\n"
  )
}


cat(
  "============================================================\n",
  "Unified rows = ",
  nrow(LONG),
  " / 1440\n",
  "Raw-vs-weighted rows = ",
  nrow(RW),
  " / 480\n",
  "Direction-pair rows = ",
  nrow(DIRPAIR),
  " / 720\n",
  "Branch-pair rows = ",
  nrow(BRPAIR),
  " / 720\n",
  "Radius-pattern groups = ",
  nrow(RADIUS),
  " / 144\n",
  "Joint rows = ",
  nrow(JOINT),
  " / 480\n",
  "Weighted correlation rows = ",
  nrow(WCOR),
  " / 240\n",
  "Hard FAIL = ",
  N_HARD_FAIL,
  "\n",
  sep = ""
)


if (
  N_HARD_FAIL > 0L
) {

  stopf(
    "07H06存在hard failure。"
  )
}

