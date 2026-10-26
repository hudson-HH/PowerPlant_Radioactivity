#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999
)


# ============================================================
# 07H07
# 论文结果表、绘图数据和正式图形准备
#
# 输入：
#   已PASS的07H06综合审计结果
#
# 本阶段：
#   - 不重新拟合模型
#   - 不修改07H01–07H06
#   - 不按P值选择fuel
#   - 不按P值选择radius
#   - 所有10个累计半径完整保留
#
# 论文主结果：
#   raw Count-only
#   raw Heat-only
#
# 空间定义敏感性：
#   Σ1/d
#   ΣHeat/d
#
# 补充条件效应：
#   Joint Count
#   Joint Heat
#
# 跨暴露定义比较效应大小时：
#   优先使用每1 SD暴露变化对应的gross beta变化
# ============================================================


# ============================================================
# 1. Packages
# ============================================================

if (
  !requireNamespace(
    "data.table",
    quietly = TRUE
  )
) {
  stop(
    "缺少R包 data.table。",
    call. = FALSE
  )
}


if (
  !requireNamespace(
    "ggplot2",
    quietly = TRUE
  )
) {
  stop(
    "缺少R包 ggplot2。",
    call. = FALSE
  )
}


suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})


# ============================================================
# 2. Helpers
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

  if (
    !file.exists(
      file
    )
  ) {

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


  if (
    length(miss) >
      0L
  ) {

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


effect_sign <- function(x) {

  fifelse(
    is.na(x),
    NA_character_,

    fifelse(
      x >
        0,
      "positive",

      fifelse(
        x <
          0,
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

  status ==
    "PASS" &
  is.finite(
    low
  ) &
  is.finite(
    high
  ) &
  (
    low >
      0 |
    high <
      0
  )
}


p_text <- function(p) {

  out <- rep(
    "",
    length(p)
  )


  idx <- is.finite(p)


  out[
    idx &
    p <
      0.001
  ] <- "<0.001"


  out[
    idx &
    p >=
      0.001
  ] <- sprintf(
    "%.3f",
    p[
      idx &
      p >=
        0.001
    ]
  )


  out
}


effect_ci_text <- function(
  status,
  est,
  low,
  high,
  digits = 4
) {

  out <- as.character(
    status
  )


  idx <-
    status ==
      "PASS" &
    is.finite(
      est
    ) &
    is.finite(
      low
    ) &
    is.finite(
      high
    )


  fmt <- paste0(
    "%.",
    digits,
    "f (%.",
    digits,
    "f, %.",
    digits,
    "f)"
  )


  out[idx] <- sprintf(
    fmt,
    est[idx],
    low[idx],
    high[idx]
  )


  out
}


component_label_fun <- function(x) {

  fcase(

    x ==
      "raw_count_only",
    "Count",

    x ==
      "raw_heat_only",
    "Heat",

    x ==
      "joint_count",
    "Joint Count",

    x ==
      "joint_heat",
    "Joint Heat",

    x ==
      "weighted_proximity",
    "Σ1/d",

    x ==
      "weighted_heat_proximity",
    "ΣHeat/d",

    default =
      x
  )
}


fuel_label_fun <- function(x) {

  fcase(

    x ==
      "coal",
    "Coal",

    x ==
      "gas",
    "Gas",

    x ==
      "oil",
    "Oil",

    x ==
      "other_fossil",
    "Other fossil",

    x ==
      "biomass",
    "Biomass",

    x ==
      "solid_waste",
    "Solid waste",

    default =
      x
  )
}


branch_label_fun <- function(x) {

  fcase(

    x ==
      "quasi_daily",
    "Quasi-daily",

    x ==
      "filter_period",
    "Filter-period",

    default =
      x
  )
}


direction_label_fun <- function(x) {

  fcase(

    x ==
      "upwind",
    "Upwind",

    x ==
      "downwind",
    "Downwind",

    default =
      x
  )
}


save_plot <- function(
  p,
  basefile,
  width,
  height
) {

  ggsave(
    filename = paste0(basefile, ".png"),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  ggsave(
    filename = paste0(basefile, ".pdf"),
    plot = p,
    width = width,
    height = height,
    units = "in",
    device = grDevices::cairo_pdf
  )
}


# ============================================================
# 3. Paths
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


H06 <- file.path(
  H_ROOT,
  "07H06_主分析与距离加权敏感性综合审计"
)


FILE_LONG <- file.path(
  H06,
  "01_07H06_全部主分析与敏感性结果长表.csv"
)


FILE_RW <- file.path(
  H06,
  "02_07H06_raw与distance_weighted稳健性比较.csv"
)


FILE_DIR <- file.path(
  H06,
  "03_07H06_upwind_downwind全部结果配对.csv"
)


FILE_BRANCH <- file.path(
  H06,
  "04_07H06_quasi_filter全部结果配对.csv"
)


FILE_RADIUS <- file.path(
  H06,
  "05_07H06_累计半径结果模式.csv"
)


FILE_JOINT <- file.path(
  H06,
  "06_07H06_joint_Count_Heat共线性解释表.csv"
)


FILE_WCOR <- file.path(
  H06,
  "07_07H06_weighted指标共线性审计.csv"
)


FILE_HARD <- file.path(
  H06,
  "10_07H06_hard_checks.csv"
)


OUT <- file.path(
  H_ROOT,
  "07H07_论文结果表与图形准备"
)


TABLE_DIR <- file.path(
  OUT,
  "01_论文结果表"
)


PLOT_DATA_DIR <- file.path(
  OUT,
  "02_绘图数据"
)


FIG_DIR <- file.path(
  OUT,
  "03_正式图形"
)


dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  PLOT_DATA_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


OUT_HARD <- file.path(
  OUT,
  "07H07_hard_checks.csv"
)


OUT_REPORT <- file.path(
  OUT,
  "07H07_论文结果与图表说明.md"
)


OUT_SESSION <- file.path(
  OUT,
  "sessionInfo.txt"
)


# ============================================================
# 4. Input files
# ============================================================

assert_file(
  FILE_LONG,
  "07H06 long table"
)

assert_file(
  FILE_RW,
  "07H06 raw-weighted comparison"
)

assert_file(
  FILE_DIR,
  "07H06 direction comparison"
)

assert_file(
  FILE_BRANCH,
  "07H06 branch comparison"
)

assert_file(
  FILE_RADIUS,
  "07H06 radius pattern"
)

assert_file(
  FILE_JOINT,
  "07H06 joint table"
)

assert_file(
  FILE_WCOR,
  "07H06 weighted correlation"
)

assert_file(
  FILE_HARD,
  "07H06 hard checks"
)


# ============================================================
# 5. Read
# ============================================================

LONG <- fread(
  FILE_LONG,
  showProgress = FALSE
)

RW <- fread(
  FILE_RW,
  showProgress = FALSE
)

DIRPAIR <- fread(
  FILE_DIR,
  showProgress = FALSE
)

BRPAIR <- fread(
  FILE_BRANCH,
  showProgress = FALSE
)

RADIUS <- fread(
  FILE_RADIUS,
  showProgress = FALSE
)

JOINT <- fread(
  FILE_JOINT,
  showProgress = FALSE
)

WCOR <- fread(
  FILE_WCOR,
  showProgress = FALSE
)

HARD06 <- fread(
  FILE_HARD,
  showProgress = FALSE
)


# ============================================================
# 6. 07H06 freeze-status audit
# ============================================================

assert_cols(
  HARD06,
  c(
    "check",
    "status"
  ),
  "07H06 hard checks"
)


N_H06_FAIL <- HARD06[
  status !=
    "PASS",
  .N
]


if (
  N_H06_FAIL >
    0L
) {

  print(
    HARD06[
      status !=
        "PASS"
    ]
  )


  stopf(
    "07H06仍存在未PASS hard check，禁止进入论文整理阶段。"
  )
}


# ============================================================
# 7. Required long-table schema
# ============================================================

REQ_LONG <- c(

  "analysis_layer",

  "branch",
  "fuel",
  "direction",
  "radius_km",
  "component",

  "status",

  "raw_increment",
  "raw_unit",

  "estimate_mBq_m3",
  "se_mBq_m3",

  "ci95_low_mBq_m3",
  "ci95_high_mBq_m3",

  "p_wald",

  "sd_raw_exposure",
  "estimate_per_1SD_mBq_m3",

  "ci95_low_per_1SD_mBq_m3",
  "ci95_high_per_1SD_mBq_m3",

  "count_heat_vif",
  "interpretability_flag",

  "note"
)


assert_cols(
  LONG,
  REQ_LONG,
  "07H06 long table"
)


# ============================================================
# 8. Exact dimensions
# ============================================================

if (
  nrow(LONG) !=
    1440L
) {

  stopf(
    "07H06 long table不是1440行；实际=",
    nrow(LONG)
  )
}


if (
  nrow(RW) !=
    480L
) {

  stopf(
    "07H06 raw-weighted table不是480行。"
  )
}


if (
  nrow(DIRPAIR) !=
    720L
) {

  stopf(
    "07H06 direction pairs不是720行。"
  )
}


if (
  nrow(BRPAIR) !=
    720L
) {

  stopf(
    "07H06 branch pairs不是720行。"
  )
}


if (
  nrow(RADIUS) !=
    144L
) {

  stopf(
    "07H06 radius groups不是144行。"
  )
}


if (
  nrow(JOINT) !=
    480L
) {

  stopf(
    "07H06 joint table不是480行。"
  )
}


if (
  nrow(WCOR) !=
    240L
) {

  stopf(
    "07H06 weighted correlation不是240行。"
  )
}


# ============================================================
# 9. Numeric normalization
# ============================================================

NUM_LONG <- c(

  "radius_km",

  "raw_increment",

  "estimate_mBq_m3",
  "se_mBq_m3",

  "ci95_low_mBq_m3",
  "ci95_high_mBq_m3",

  "p_wald",

  "sd_raw_exposure",

  "estimate_per_1SD_mBq_m3",

  "ci95_low_per_1SD_mBq_m3",
  "ci95_high_per_1SD_mBq_m3",

  "count_heat_vif"
)


for (
  cc in NUM_LONG
) {

  LONG[
    ,
    (cc) :=
      safe_num(
        get(cc)
      )
  ]
}


# ============================================================
# 10. Derive standardized 1-SD CI
#
# For a linear contrast:
#
# effect(1 SD)
# =
# effect(fixed increment)
# × SD / fixed increment
#
# 07H05 already contains 1-SD CI.
# Frozen raw models may not contain their 1-SD CI, so derive them
# from the same linear contrast.
# ============================================================

LONG[
  ,
  std_ratio :=
    fifelse(
      is.finite(
        sd_raw_exposure
      ) &
      is.finite(
        raw_increment
      ) &
      raw_increment >
        0,

      sd_raw_exposure /
        raw_increment,

      NA_real_
    )
]


LONG[
  ,
  derived_estimate_1SD :=
    estimate_mBq_m3 *
    std_ratio
]


LONG[
  ,
  derived_ci_low_1SD :=
    ci95_low_mBq_m3 *
    std_ratio
]


LONG[
  ,
  derived_ci_high_1SD :=
    ci95_high_mBq_m3 *
    std_ratio
]


# Use original saved 1-SD estimate when available,
# otherwise use exactly equivalent derived estimate.
LONG[
  ,
  estimate_1SD :=
    fifelse(
      is.finite(
        estimate_per_1SD_mBq_m3
      ),

      estimate_per_1SD_mBq_m3,

      derived_estimate_1SD
    )
]


LONG[
  ,
  ci_low_1SD :=
    fifelse(
      is.finite(
        ci95_low_per_1SD_mBq_m3
      ),

      ci95_low_per_1SD_mBq_m3,

      derived_ci_low_1SD
    )
]


LONG[
  ,
  ci_high_1SD :=
    fifelse(
      is.finite(
        ci95_high_per_1SD_mBq_m3
      ),

      ci95_high_per_1SD_mBq_m3,

      derived_ci_high_1SD
    )
]


# ============================================================
# 11. Standardization audit
# ============================================================

STD_EST_MISMATCH <- LONG[
  status ==
    "PASS" &
  is.finite(
    estimate_per_1SD_mBq_m3
  ) &
  is.finite(
    derived_estimate_1SD
  ) &
  abs(
    estimate_per_1SD_mBq_m3 -
    derived_estimate_1SD
  ) >
    1e-8 *
    pmax(
      1,
      abs(
        estimate_per_1SD_mBq_m3
      )
    ),
  .N
]


STD_CI_LOW_MISMATCH <- LONG[
  status ==
    "PASS" &
  is.finite(
    ci95_low_per_1SD_mBq_m3
  ) &
  is.finite(
    derived_ci_low_1SD
  ) &
  abs(
    ci95_low_per_1SD_mBq_m3 -
    derived_ci_low_1SD
  ) >
    1e-8 *
    pmax(
      1,
      abs(
        ci95_low_per_1SD_mBq_m3
      )
    ),
  .N
]


STD_CI_HIGH_MISMATCH <- LONG[
  status ==
    "PASS" &
  is.finite(
    ci95_high_per_1SD_mBq_m3
  ) &
  is.finite(
    derived_ci_high_1SD
  ) &
  abs(
    ci95_high_per_1SD_mBq_m3 -
    derived_ci_high_1SD
  ) >
    1e-8 *
    pmax(
      1,
      abs(
        ci95_high_per_1SD_mBq_m3
      )
    ),
  .N
]


if (
  STD_EST_MISMATCH >
    0L
) {

  stopf(
    "1-SD estimate与固定增量线性换算不一致：",
    STD_EST_MISMATCH
  )
}


if (
  STD_CI_LOW_MISMATCH >
    0L |
  STD_CI_HIGH_MISMATCH >
    0L
) {

  stopf(
    "1-SD CI与固定增量线性换算不一致。"
  )
}


# ============================================================
# 12. Paper-friendly labels
# ============================================================

LONG[
  ,
  `:=`(

    component_label =
      component_label_fun(
        component
      ),

    fuel_label =
      fuel_label_fun(
        fuel
      ),

    branch_label =
      branch_label_fun(
        branch
      ),

    direction_label =
      direction_label_fun(
        direction
      ),

    effect_sign =
      effect_sign(
        estimate_mBq_m3
      ),

    ci_excludes_zero =
      ci_excludes_zero(
        status,
        ci95_low_mBq_m3,
        ci95_high_mBq_m3
      ),

    ci_1SD_excludes_zero =
      ci_excludes_zero(
        status,
        ci_low_1SD,
        ci_high_1SD
      ),

    p_text =
      p_text(
        p_wald
      )
  )
]


LONG[
  ,
  effect_fixed_95CI :=
    effect_ci_text(
      status,
      estimate_mBq_m3,
      ci95_low_mBq_m3,
      ci95_high_mBq_m3,
      digits = 4
    )
]


LONG[
  ,
  effect_1SD_95CI :=
    effect_ci_text(
      status,
      estimate_1SD,
      ci_low_1SD,
      ci_high_1SD,
      digits = 4
    )
]


# ============================================================
# 13. Component definitions
# ============================================================

MAIN_SINGLE <- c(
  "raw_count_only",
  "raw_heat_only",
  "weighted_proximity",
  "weighted_heat_proximity"
)


RAW_MAIN <- c(
  "raw_count_only",
  "raw_heat_only"
)


JOINT_COMPONENTS <- c(
  "joint_count",
  "joint_heat"
)


# ============================================================
# 14. Paper Table 1
# Upwind single-exposure results, ALL radii
# ============================================================

TABLE1 <- LONG[
  direction ==
    "upwind" &
  component %in%
    MAIN_SINGLE,
  .(

    branch,
    branch_label,

    fuel,
    fuel_label,

    radius_km,

    component,
    component_label,

    analysis_layer,

    status,

    raw_increment,
    raw_unit,

    sd_raw_exposure,

    estimate_fixed =
      estimate_mBq_m3,

    ci_low_fixed =
      ci95_low_mBq_m3,

    ci_high_fixed =
      ci95_high_mBq_m3,

    effect_fixed_95CI,

    estimate_1SD,
    ci_low_1SD,
    ci_high_1SD,

    effect_1SD_95CI,

    p_wald,
    p_text,

    effect_sign,
    ci_excludes_zero,

    note
  )
]


setorder(
  TABLE1,
  fuel,
  branch,
  component,
  radius_km
)


if (
  nrow(TABLE1) !=
    480L
) {

  stopf(
    "论文upwind单暴露表不是480行；实际=",
    nrow(TABLE1)
  )
}


fwrite(
  TABLE1,
  file.path(
    TABLE_DIR,
    "01_论文主结果_upwind单暴露_全部预设半径.csv"
  )
)


# ============================================================
# 15. Paper Table 2
# Joint Count / Heat conditional effects
# ============================================================

TABLE2 <- LONG[
  direction ==
    "upwind" &
  component %in%
    JOINT_COMPONENTS,
  .(

    branch,
    branch_label,

    fuel,
    fuel_label,

    radius_km,

    component,
    component_label,

    status,

    raw_increment,
    raw_unit,

    estimate_fixed =
      estimate_mBq_m3,

    ci_low_fixed =
      ci95_low_mBq_m3,

    ci_high_fixed =
      ci95_high_mBq_m3,

    effect_fixed_95CI,

    estimate_1SD,
    ci_low_1SD,
    ci_high_1SD,

    effect_1SD_95CI,

    p_wald,
    p_text,

    count_heat_vif,
    interpretability_flag,

    note
  )
]


setorder(
  TABLE2,
  fuel,
  branch,
  component,
  radius_km
)


if (
  nrow(TABLE2) !=
    240L
) {

  stopf(
    "upwind joint表不是240行；实际=",
    nrow(TABLE2)
  )
}


fwrite(
  TABLE2,
  file.path(
    TABLE_DIR,
    "02_论文补充结果_Joint_Count_Heat_upwind_全部半径.csv"
  )
)


# ============================================================
# 16. Raw vs weighted robustness summary
# ============================================================

assert_cols(
  RW,
  c(
    "comparison_family",
    "branch",
    "fuel",
    "direction",
    "radius_km",
    "raw_status",
    "weighted_status",
    "same_sign_fixed_increment",
    "same_sign_1SD",
    "raw_ci_excludes_zero",
    "weighted_ci_excludes_zero",
    "robustness_pattern"
  ),
  "raw-vs-weighted"
)


ROBUST_SUMMARY <- RW[
  ,
  .(

    n_pairs =
      .N,

    n_both_pass =
      sum(
        raw_status == "PASS" &
        weighted_status == "PASS"
      ),

    n_pair_incomplete =
      sum(
        robustness_pattern == "PAIR_INCOMPLETE"
      ),

    n_same_sign =
      sum(
        raw_status == "PASS" &
        weighted_status == "PASS" &
        same_sign_fixed_increment %in% TRUE,
        na.rm = TRUE
      ),

    n_opposite_sign =
      sum(
        raw_status == "PASS" &
        weighted_status == "PASS" &
        same_sign_fixed_increment %in% FALSE,
        na.rm = TRUE
      ),

    same_sign_pct_among_comparable =
      100 *
      sum(
        raw_status == "PASS" &
        weighted_status == "PASS" &
        same_sign_fixed_increment %in% TRUE,
        na.rm = TRUE
      ) /
      sum(
        raw_status == "PASS" &
        weighted_status == "PASS"
      ),

    n_both_ci_exclude_zero =
      sum(
        robustness_pattern %in% c(
          "BOTH_CI_EXCLUDE_ZERO_SAME_SIGN",
          "BOTH_CI_EXCLUDE_ZERO_OPPOSITE_SIGN"
        )
      ),

    n_raw_only_ci_excludes_zero =
      sum(
        robustness_pattern ==
          "RAW_ONLY_CI_EXCLUDES_ZERO"
      ),

    n_weighted_only_ci_excludes_zero =
      sum(
        robustness_pattern ==
          "WEIGHTED_ONLY_CI_EXCLUDES_ZERO"
      ),

    n_neither_ci_excludes_zero =
      sum(
        robustness_pattern %in% c(
          "NEITHER_CI_EXCLUDES_ZERO_SAME_SIGN",
          "NEITHER_CI_EXCLUDES_ZERO_OPPOSITE_SIGN"
        )
      )
  ),

  by = .(
    comparison_family,
    branch,
    fuel
  )
]


ROBUST_OVERALL <- RW[
  ,
  .(

    comparison_family =
      "ALL",

    branch =
      "ALL",

    fuel =
      "ALL",

    n_pairs =
      .N,

    n_both_pass =
      sum(
        raw_status == "PASS" &
        weighted_status == "PASS"
      ),

    n_pair_incomplete =
      sum(
        robustness_pattern == "PAIR_INCOMPLETE"
      ),

    n_same_sign =
      sum(
        raw_status == "PASS" &
        weighted_status == "PASS" &
        same_sign_fixed_increment %in% TRUE,
        na.rm = TRUE
      ),

    n_opposite_sign =
      sum(
        raw_status == "PASS" &
        weighted_status == "PASS" &
        same_sign_fixed_increment %in% FALSE,
        na.rm = TRUE
      ),

    same_sign_pct_among_comparable =
      100 *
      sum(
        raw_status == "PASS" &
        weighted_status == "PASS" &
        same_sign_fixed_increment %in% TRUE,
        na.rm = TRUE
      ) /
      sum(
        raw_status == "PASS" &
        weighted_status == "PASS"
      ),

    n_both_ci_exclude_zero =
      sum(
        robustness_pattern %in% c(
          "BOTH_CI_EXCLUDE_ZERO_SAME_SIGN",
          "BOTH_CI_EXCLUDE_ZERO_OPPOSITE_SIGN"
        )
      ),

    n_raw_only_ci_excludes_zero =
      sum(
        robustness_pattern ==
          "RAW_ONLY_CI_EXCLUDES_ZERO"
      ),

    n_weighted_only_ci_excludes_zero =
      sum(
        robustness_pattern ==
          "WEIGHTED_ONLY_CI_EXCLUDES_ZERO"
      ),

    n_neither_ci_excludes_zero =
      sum(
        robustness_pattern %in% c(
          "NEITHER_CI_EXCLUDES_ZERO_SAME_SIGN",
          "NEITHER_CI_EXCLUDES_ZERO_OPPOSITE_SIGN"
        )
      )
  )
]


ROBUST_SUMMARY <- rbindlist(
  list(
    ROBUST_OVERALL,
    ROBUST_SUMMARY
  ),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  ROBUST_SUMMARY,
  file.path(
    TABLE_DIR,
    "03_raw与distance_weighted稳健性汇总.csv"
  )
)


# ============================================================
# 17. Directional negative-control summary
# ============================================================

assert_cols(
  DIRPAIR,
  c(
    "branch",
    "fuel",
    "radius_km",
    "component",
    "upwind_status",
    "downwind_status",
    "same_sign",
    "same_sign_1SD",
    "ci_pattern"
  ),
  "direction pair table"
)


DIR_SUMMARY <- DIRPAIR[
  ,
  .(

    n_radius_pairs =
      .N,

    n_both_estimable =
      sum(
        upwind_status ==
          "PASS" &
        downwind_status ==
          "PASS"
      ),

    n_pair_incomplete =
      sum(
        ci_pattern ==
          "PAIR_INCOMPLETE"
      ),

    n_both_ci_exclude_zero =
      sum(
        ci_pattern ==
          "BOTH_CI_EXCLUDE_ZERO"
      ),

    n_upwind_only_ci_excludes_zero =
      sum(
        ci_pattern ==
          "UPWIND_ONLY_CI_EXCLUDES_ZERO"
      ),

    n_downwind_only_ci_excludes_zero =
      sum(
        ci_pattern ==
          "DOWNWIND_ONLY_CI_EXCLUDES_ZERO"
      ),

    n_neither_ci_excludes_zero =
      sum(
        ci_pattern ==
          "NEITHER_CI_EXCLUDES_ZERO"
      ),

    n_same_sign =
      sum(
        same_sign %in%
          TRUE,
        na.rm = TRUE
      ),

    n_opposite_sign =
      sum(
        same_sign %in%
          FALSE,
        na.rm = TRUE
      )
  ),

  by = .(
    branch,
    fuel,
    component
  )
]


DIR_SUMMARY[
  ,
  `:=`(
    fuel_label =
      fuel_label_fun(
        fuel
      ),

    component_label =
      component_label_fun(
        component
      )
  )
]


setorder(
  DIR_SUMMARY,
  branch,
  fuel,
  component
)


fwrite(
  DIR_SUMMARY,
  file.path(
    TABLE_DIR,
    "04_upwind_downwind方向负对照汇总.csv"
  )
)


# ============================================================
# 18. Quasi vs filter consistency summary
# ============================================================

assert_cols(
  BRPAIR,
  c(
    "fuel",
    "direction",
    "radius_km",
    "component",
    "quasi_status",
    "filter_status",
    "same_sign",
    "same_sign_1SD",
    "ci_pattern"
  ),
  "branch pair table"
)


BR_SUMMARY <- BRPAIR[
  ,
  .(

    n_radius_pairs =
      .N,

    n_both_estimable =
      sum(
        quasi_status ==
          "PASS" &
        filter_status ==
          "PASS"
      ),

    n_pair_incomplete =
      sum(
        ci_pattern ==
          "PAIR_INCOMPLETE"
      ),

    n_both_ci_exclude_zero =
      sum(
        ci_pattern ==
          "BOTH_CI_EXCLUDE_ZERO"
      ),

    n_quasi_only_ci_excludes_zero =
      sum(
        ci_pattern ==
          "QUASI_ONLY_CI_EXCLUDES_ZERO"
      ),

    n_filter_only_ci_excludes_zero =
      sum(
        ci_pattern ==
          "FILTER_ONLY_CI_EXCLUDES_ZERO"
      ),

    n_neither_ci_excludes_zero =
      sum(
        ci_pattern ==
          "NEITHER_CI_EXCLUDES_ZERO"
      ),

    n_same_sign =
      sum(
        same_sign %in%
          TRUE,
        na.rm = TRUE
      ),

    n_opposite_sign =
      sum(
        same_sign %in%
          FALSE,
        na.rm = TRUE
      )
  ),

  by = .(
    fuel,
    direction,
    component
  )
]


BR_SUMMARY[
  ,
  `:=`(
    fuel_label =
      fuel_label_fun(
        fuel
      ),

    component_label =
      component_label_fun(
        component
      )
  )
]


setorder(
  BR_SUMMARY,
  fuel,
  direction,
  component
)


fwrite(
  BR_SUMMARY,
  file.path(
    TABLE_DIR,
    "05_quasi_filter时间尺度一致性汇总.csv"
  )
)


# ============================================================
# 19. Model coverage / SKIP table
# ============================================================

COVERAGE <- LONG[
  ,
  .(

    n_results =
      .N,

    n_pass =
      sum(
        status ==
          "PASS"
      ),

    n_skip =
      sum(
        status ==
          "SKIP"
      ),

    n_fail =
      sum(
        status ==
          "FAIL"
      )
  ),

  by = .(
    analysis_layer,
    branch,
    fuel,
    component
  )
]


COVERAGE[
  ,
  `:=`(
    fuel_label =
      fuel_label_fun(
        fuel
      ),

    component_label =
      component_label_fun(
        component
      )
  )
]


setorder(
  COVERAGE,
  branch,
  fuel,
  component
)


fwrite(
  COVERAGE,
  file.path(
    TABLE_DIR,
    "06_模型覆盖与SKIP情况.csv"
  )
)


# ============================================================
# 20. Radius-pattern table
# ============================================================

RADIUS_OUT <- copy(
  RADIUS
)


RADIUS_OUT[
  ,
  `:=`(
    fuel_label =
      fuel_label_fun(
        fuel
      ),

    component_label =
      component_label_fun(
        component
      ),

    branch_label =
      branch_label_fun(
        branch
      ),

    direction_label =
      direction_label_fun(
        direction
      )
  )
]


setorder(
  RADIUS_OUT,
  branch,
  fuel,
  direction,
  component
)


fwrite(
  RADIUS_OUT,
  file.path(
    TABLE_DIR,
    "07_全部预设累计半径结果模式.csv"
  )
)


# ============================================================
# 21. Joint VIF supplementary table
# ============================================================

JOINT_OUT <- copy(
  JOINT
)


JOINT_OUT[
  ,
  `:=`(
    fuel_label =
      fuel_label_fun(
        fuel
      ),

    component_label =
      component_label_fun(
        component
      )
  )
]


fwrite(
  JOINT_OUT,
  file.path(
    TABLE_DIR,
    "08_Joint_Count_Heat共线性解释表.csv"
  )
)


# ============================================================
# 22. Weighted collinearity supplementary table
# ============================================================

fwrite(
  WCOR,
  file.path(
    TABLE_DIR,
    "09_Σ1d与ΣHeatd共线性审计.csv"
  )
)


# ============================================================
# 23. Plot data
# ============================================================

PLOT_UPWIND <- LONG[
  direction ==
    "upwind" &
  component %in%
    MAIN_SINGLE,
  .(

    branch,
    branch_label,

    fuel,
    fuel_label,

    radius_km,

    component,
    component_label,

    status,

    estimate_1SD,
    ci_low_1SD,
    ci_high_1SD,

    p_wald,
    ci_1SD_excludes_zero
  )
]


fwrite(
  PLOT_UPWIND,
  file.path(
    PLOT_DATA_DIR,
    "01_upwind四类单暴露_1SD绘图数据.csv"
  )
)


PLOT_COUNT_RW <- PLOT_UPWIND[
  component %in%
    c(
      "raw_count_only",
      "weighted_proximity"
    )
]


PLOT_COUNT_RW[
  ,
  exposure_definition :=
    fifelse(
      component ==
        "raw_count_only",

      "Count",

      "Σ1/d"
    )
]


fwrite(
  PLOT_COUNT_RW,
  file.path(
    PLOT_DATA_DIR,
    "02_Count_vs_Σ1d_1SD绘图数据.csv"
  )
)


PLOT_HEAT_RW <- PLOT_UPWIND[
  component %in%
    c(
      "raw_heat_only",
      "weighted_heat_proximity"
    )
]


PLOT_HEAT_RW[
  ,
  exposure_definition :=
    fifelse(
      component ==
        "raw_heat_only",

      "Heat",

      "ΣHeat/d"
    )
]


fwrite(
  PLOT_HEAT_RW,
  file.path(
    PLOT_DATA_DIR,
    "03_Heat_vs_ΣHeatd_1SD绘图数据.csv"
  )
)


PLOT_DIRECTION <- LONG[
  component %in%
    MAIN_SINGLE,
  .(

    branch,
    branch_label,

    fuel,
    fuel_label,

    direction,
    direction_label,

    radius_km,

    component,
    component_label,

    status,

    estimate_1SD,
    ci_low_1SD,
    ci_high_1SD,

    ci_1SD_excludes_zero
  )
]


fwrite(
  PLOT_DIRECTION,
  file.path(
    PLOT_DATA_DIR,
    "04_upwind_downwind四类单暴露_1SD绘图数据.csv"
  )
)


PLOT_BRANCH <- LONG[
  direction ==
    "upwind" &
  component %in%
    MAIN_SINGLE,
  .(

    branch,
    branch_label,

    fuel,
    fuel_label,

    radius_km,

    component,
    component_label,

    status,

    estimate_1SD,
    ci_low_1SD,
    ci_high_1SD,

    ci_1SD_excludes_zero
  )
]


fwrite(
  PLOT_BRANCH,
  file.path(
    PLOT_DATA_DIR,
    "05_quasi_filter_upwind四类单暴露_1SD绘图数据.csv"
  )
)


# ============================================================
# 24. Factor order for figures
# ============================================================

FUEL_ORDER <- c(
  "Coal",
  "Gas",
  "Oil",
  "Other fossil",
  "Biomass",
  "Solid waste"
)


BRANCH_ORDER <- c(
  "Quasi-daily",
  "Filter-period"
)


PLOT_UPWIND[
  ,
  fuel_label :=
    factor(
      fuel_label,
      levels =
        FUEL_ORDER
    )
]


PLOT_UPWIND[
  ,
  branch_label :=
    factor(
      branch_label,
      levels =
        BRANCH_ORDER
    )
]


PLOT_COUNT_RW[
  ,
  fuel_label :=
    factor(
      fuel_label,
      levels =
        FUEL_ORDER
    )
]


PLOT_COUNT_RW[
  ,
  branch_label :=
    factor(
      branch_label,
      levels =
        BRANCH_ORDER
    )
]


PLOT_HEAT_RW[
  ,
  fuel_label :=
    factor(
      fuel_label,
      levels =
        FUEL_ORDER
    )
]


PLOT_HEAT_RW[
  ,
  branch_label :=
    factor(
      branch_label,
      levels =
        BRANCH_ORDER
    )
]


PLOT_DIRECTION[
  ,
  fuel_label :=
    factor(
      fuel_label,
      levels =
        FUEL_ORDER
    )
]


PLOT_DIRECTION[
  ,
  branch_label :=
    factor(
      branch_label,
      levels =
        BRANCH_ORDER
    )
]


PLOT_BRANCH[
  ,
  fuel_label :=
    factor(
      fuel_label,
      levels =
        FUEL_ORDER
    )
]


PLOT_BRANCH[
  ,
  branch_label :=
    factor(
      branch_label,
      levels =
        BRANCH_ORDER
    )
]


# ============================================================
# 25. Common figure theme
# ============================================================

paper_theme <- theme_bw(
  base_size = 11
) +
  theme(
    panel.grid.minor =
      element_blank(),

    legend.position =
      "bottom",

    strip.background =
      element_rect(
        fill = "white"
      ),

    axis.title =
      element_text(
        face = "bold"
      ),

    plot.title =
      element_text(
        face = "bold"
      )
  )


# ============================================================
# 26. Figure 1
# Raw Count / Heat upwind main analysis
# ============================================================

FIG1_DATA <- PLOT_UPWIND[
  status ==
    "PASS" &
  component %in%
    RAW_MAIN
]


FIG1_DATA[
  ,
  component_label :=
    factor(
      component_label,
      levels =
        c(
          "Count",
          "Heat"
        )
    )
]


P1 <- ggplot(
  FIG1_DATA,
  aes(
    x =
      radius_km,
    y =
      estimate_1SD,
    group =
      fuel_label,
    colour =
      fuel_label
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +

  geom_errorbar(
    aes(
      ymin =
        ci_low_1SD,
      ymax =
        ci_high_1SD
    ),
    width = 1.5,
    alpha = 0.35
  ) +

  geom_line(
    linewidth = 0.7
  ) +

  geom_point(
    size = 1.8
  ) +

  facet_wrap(
    vars(
      branch_label,
      component_label
    ),
    ncol = 2,
    scales = "free_y"
  ) +

  scale_x_continuous(
    breaks =
      c(
        20,
        50,
        100,
        150,
        200
      )
  ) +

  labs(
    title =
      "Fuel-specific upwind associations: raw Count and Heat",

    x =
      "Cumulative radius (km)",

    y =
      "Change in gross beta (mBq/m³) per 1 SD exposure",

    colour =
      "Fuel"
  ) +

  paper_theme


save_plot(
  P1,
  file.path(
    FIG_DIR,
    "Fig1_raw主分析_upwind_Count_Heat_1SD"
  ),
  width = 11,
  height = 8
)


# ============================================================
# 27. Figure 2
# Count vs Σ1/d
# ============================================================

FIG2_DATA <- PLOT_COUNT_RW[
  status ==
    "PASS"
]


P2 <- ggplot(
  FIG2_DATA,
  aes(
    x =
      radius_km,
    y =
      estimate_1SD,
    group =
      exposure_definition,
    colour =
      exposure_definition
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +

  geom_errorbar(
    aes(
      ymin =
        ci_low_1SD,
      ymax =
        ci_high_1SD
    ),
    width = 1.5,
    alpha = 0.35
  ) +

  geom_line(
    linewidth = 0.7
  ) +

  geom_point(
    size = 1.8
  ) +

  facet_wrap(
    vars(
      branch_label,
      fuel_label
    ),
    ncol = 6,
    scales = "free_y"
  ) +

  scale_x_continuous(
    breaks =
      c(
        20,
        50,
        100,
        150,
        200
      )
  ) +

  labs(
    title =
      "Sensitivity to inverse-distance weighting: Count versus Σ1/d",

    x =
      "Cumulative radius (km)",

    y =
      "Change in gross beta (mBq/m³) per 1 SD exposure",

    colour =
      "Exposure definition"
  ) +

  paper_theme


save_plot(
  P2,
  file.path(
    FIG_DIR,
    "Fig2_Count_vs_Σ1d_距离加权敏感性_1SD"
  ),
  width = 15,
  height = 7
)


# ============================================================
# 28. Figure 3
# Heat vs ΣHeat/d
# ============================================================

FIG3_DATA <- PLOT_HEAT_RW[
  status ==
    "PASS"
]


P3 <- ggplot(
  FIG3_DATA,
  aes(
    x =
      radius_km,
    y =
      estimate_1SD,
    group =
      exposure_definition,
    colour =
      exposure_definition
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +

  geom_errorbar(
    aes(
      ymin =
        ci_low_1SD,
      ymax =
        ci_high_1SD
    ),
    width = 1.5,
    alpha = 0.35
  ) +

  geom_line(
    linewidth = 0.7
  ) +

  geom_point(
    size = 1.8
  ) +

  facet_wrap(
    vars(
      branch_label,
      fuel_label
    ),
    ncol = 6,
    scales = "free_y"
  ) +

  scale_x_continuous(
    breaks =
      c(
        20,
        50,
        100,
        150,
        200
      )
  ) +

  labs(
    title =
      "Sensitivity to inverse-distance weighting: Heat versus ΣHeat/d",

    x =
      "Cumulative radius (km)",

    y =
      "Change in gross beta (mBq/m³) per 1 SD exposure",

    colour =
      "Exposure definition"
  ) +

  paper_theme


save_plot(
  P3,
  file.path(
    FIG_DIR,
    "Fig3_Heat_vs_ΣHeatd_距离加权敏感性_1SD"
  ),
  width = 15,
  height = 7
)


# ============================================================
# 29. Figure 4
# Upwind vs downwind directional-control comparison
# Single exposures only
# ============================================================

FIG4_DATA <- copy(
  PLOT_DIRECTION[
    status ==
      "PASS"
  ]
)


FIG4_DATA[
  ,
  direction_label :=
    factor(
      direction_label,
      levels =
        c(
          "Upwind",
          "Downwind"
        )
    )
]


P4 <- ggplot(
  FIG4_DATA,
  aes(
    x =
      radius_km,
    y =
      estimate_1SD,
    group =
      interaction(
        direction_label,
        branch_label
      ),
    colour =
      direction_label,
    linetype =
      branch_label
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.35
  ) +

  geom_line(
    linewidth = 0.65
  ) +

  geom_point(
    size = 1.3
  ) +

  facet_grid(
    component_label ~ fuel_label,
    scales = "free_y"
  ) +

  scale_x_continuous(
    breaks =
      c(
        20,
        100,
        200
      )
  ) +

  labs(
    title =
      "Directional-control comparison across exposure definitions",

    x =
      "Cumulative radius (km)",

    y =
      "Change in gross beta (mBq/m³) per 1 SD exposure",

    colour =
      "Direction",

    linetype =
      "Analysis branch"
  ) +

  paper_theme +

  theme(
    axis.text.x =
      element_text(
        size = 7
      ),

    strip.text =
      element_text(
        size = 8
      )
  )


save_plot(
  P4,
  file.path(
    FIG_DIR,
    "FigS1_upwind_downwind方向负对照_四类单暴露_1SD"
  ),
  width = 18,
  height = 12
)


# ============================================================
# 30. Figure 5
# Quasi vs filter upwind comparison
# ============================================================

FIG5_DATA <- copy(
  PLOT_BRANCH[
    status ==
      "PASS"
  ]
)


P5 <- ggplot(
  FIG5_DATA,
  aes(
    x =
      radius_km,
    y =
      estimate_1SD,
    group =
      branch_label,
    colour =
      branch_label
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.35
  ) +

  geom_line(
    linewidth = 0.65
  ) +

  geom_point(
    size = 1.3
  ) +

  facet_grid(
    component_label ~ fuel_label,
    scales = "free_y"
  ) +

  scale_x_continuous(
    breaks =
      c(
        20,
        100,
        200
      )
  ) +

  labs(
    title =
      "Quasi-daily versus filter-period consistency in upwind analyses",

    x =
      "Cumulative radius (km)",

    y =
      "Change in gross beta (mBq/m³) per 1 SD exposure",

    colour =
      "Analysis branch"
  ) +

  paper_theme +

  theme(
    axis.text.x =
      element_text(
        size = 7
      ),

    strip.text =
      element_text(
        size = 8
      )
  )


save_plot(
  P5,
  file.path(
    FIG_DIR,
    "FigS2_quasi_filter时间尺度一致性_upwind_1SD"
  ),
  width = 18,
  height = 12
)


# ============================================================
# 31. Figure 6
# Joint Count / Heat, upwind
# Supplementary only
# ============================================================

FIG6_DATA <- LONG[
  direction ==
    "upwind" &
  component %in%
    JOINT_COMPONENTS &
  status ==
    "PASS"
]


FIG6_DATA[
  ,
  `:=`(

    component_label =
      factor(
        component_label_fun(
          component
        ),
        levels =
          c(
            "Joint Count",
            "Joint Heat"
          )
      ),

    fuel_label =
      factor(
        fuel_label_fun(
          fuel
        ),
        levels =
          FUEL_ORDER
      ),

    branch_label =
      factor(
        branch_label_fun(
          branch
        ),
        levels =
          BRANCH_ORDER
      )
  )
]


P6 <- ggplot(
  FIG6_DATA,
  aes(
    x =
      radius_km,
    y =
      estimate_1SD,
    group =
      fuel_label,
    colour =
      fuel_label
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.4
  ) +

  geom_errorbar(
    aes(
      ymin =
        ci_low_1SD,
      ymax =
        ci_high_1SD
    ),
    width = 1.5,
    alpha = 0.35
  ) +

  geom_line(
    linewidth = 0.7
  ) +

  geom_point(
    size = 1.8
  ) +

  facet_wrap(
    vars(
      branch_label,
      component_label
    ),
    ncol = 2,
    scales = "free_y"
  ) +

  scale_x_continuous(
    breaks =
      c(
        20,
        50,
        100,
        150,
        200
      )
  ) +

  labs(
    title =
      "Conditional Joint Count and Joint Heat associations",

    subtitle =
      "Interpret jointly with Count–Heat VIF diagnostics",

    x =
      "Cumulative radius (km)",

    y =
      "Change in gross beta (mBq/m³) per 1 SD exposure",

    colour =
      "Fuel"
  ) +

  paper_theme


save_plot(
  P6,
  file.path(
    FIG_DIR,
    "FigS3_Joint_Count_Heat_upwind_1SD"
  ),
  width = 11,
  height = 8
)


# ============================================================
# 32. Key aggregate numbers for manuscript/report
# ============================================================

N_RW_COMPARABLE <- RW[
  raw_status == "PASS" &
  weighted_status == "PASS",
  .N
]


N_PAIR_INCOMPLETE_CI <- RW[
  robustness_pattern ==
    "PAIR_INCOMPLETE",
  .N
]


N_RW_SAME <- RW[
  raw_status == "PASS" &
  weighted_status == "PASS" &
  same_sign_fixed_increment %in% TRUE,
  .N
]


N_RW_OPPOSITE <- RW[
  raw_status == "PASS" &
  weighted_status == "PASS" &
  same_sign_fixed_increment %in% FALSE,
  .N
]


PCT_RW_SAME <- 100 *
  N_RW_SAME /
  N_RW_COMPARABLE


N_BOTH_CI <- RW[
  robustness_pattern %in% c(
    "BOTH_CI_EXCLUDE_ZERO_SAME_SIGN",
    "BOTH_CI_EXCLUDE_ZERO_OPPOSITE_SIGN"
  ),
  .N
]


N_RAW_ONLY_CI <- RW[
  robustness_pattern ==
    "RAW_ONLY_CI_EXCLUDES_ZERO",
  .N
]


N_WEIGHTED_ONLY_CI <- RW[
  robustness_pattern ==
    "WEIGHTED_ONLY_CI_EXCLUDES_ZERO",
  .N
]


N_NEITHER_CI <- RW[
  robustness_pattern %in% c(
    "NEITHER_CI_EXCLUDES_ZERO_SAME_SIGN",
    "NEITHER_CI_EXCLUDES_ZERO_OPPOSITE_SIGN"
  ),
  .N
]


if (
  N_BOTH_CI +
  N_RAW_ONLY_CI +
  N_WEIGHTED_ONLY_CI +
  N_NEITHER_CI !=
    N_RW_COMPARABLE
) {
  stopf(
    "可比较CI分类数量不守恒：",
    N_BOTH_CI +
    N_RAW_ONLY_CI +
    N_WEIGHTED_ONLY_CI +
    N_NEITHER_CI,
    " != ",
    N_RW_COMPARABLE
  )
}


if (
  N_RW_COMPARABLE +
  N_PAIR_INCOMPLETE_CI !=
    nrow(RW)
) {
  stopf(
    "总配对数量不守恒：",
    N_RW_COMPARABLE +
    N_PAIR_INCOMPLETE_CI,
    " != ",
    nrow(RW)
  )
}

weighted_proximity_dir <- DIRPAIR[
  component ==
    "weighted_proximity",
  .N,
  by = ci_pattern
]


get_pattern_n <- function(
  D,
  pattern
) {

  z <- D[
    ci_pattern ==
      pattern,
    N
  ]

  if (
    length(z) ==
      0L
  ) {
    return(0L)
  }

  as.integer(z[1])
}


WP_UP_ONLY <- get_pattern_n(
  weighted_proximity_dir,
  "UPWIND_ONLY_CI_EXCLUDES_ZERO"
)


WP_DOWN_ONLY <- get_pattern_n(
  weighted_proximity_dir,
  "DOWNWIND_ONLY_CI_EXCLUDES_ZERO"
)


WP_BOTH <- get_pattern_n(
  weighted_proximity_dir,
  "BOTH_CI_EXCLUDE_ZERO"
)


# ============================================================
# 33. Hard checks
# ============================================================

N_TABLE1_DUP <- TABLE1[
  ,
  .N,
  by = .(
    branch,
    fuel,
    radius_km,
    component
  )
][
  N >
    1L,
  .N
]


N_TABLE2_DUP <- TABLE2[
  ,
  .N,
  by = .(
    branch,
    fuel,
    radius_km,
    component
  )
][
  N >
    1L,
  .N
]


EXPECTED_FIG_BASE <- c(

  "Fig1_raw主分析_upwind_Count_Heat_1SD",

  "Fig2_Count_vs_Σ1d_距离加权敏感性_1SD",

  "Fig3_Heat_vs_ΣHeatd_距离加权敏感性_1SD",

  "FigS1_upwind_downwind方向负对照_四类单暴露_1SD",

  "FigS2_quasi_filter时间尺度一致性_upwind_1SD",

  "FigS3_Joint_Count_Heat_upwind_1SD"
)


EXPECTED_PNG <- file.path(
  FIG_DIR,
  paste0(
    EXPECTED_FIG_BASE,
    ".png"
  )
)


EXPECTED_PDF <- file.path(
  FIG_DIR,
  paste0(
    EXPECTED_FIG_BASE,
    ".pdf"
  )
)


HARD07 <- data.table(

  check = c(

    "07H06_hard_checks_all_pass",

    "input_long_rows_1440",

    "input_raw_weighted_rows_480",

    "input_direction_pairs_720",

    "input_branch_pairs_720",

    "input_radius_groups_144",

    "input_joint_rows_480",

    "input_weighted_correlation_rows_240",

    "standardized_effect_mismatch_zero",

    "standardized_CI_low_mismatch_zero",

    "standardized_CI_high_mismatch_zero",

    "upwind_single_table_rows_480",

    "upwind_single_duplicate_keys_zero",

    "upwind_joint_table_rows_240",

    "upwind_joint_duplicate_keys_zero",

    "raw_weighted_comparable_pairs_426",

    "raw_weighted_same_sign_401",

    "raw_weighted_opposite_sign_25",

    "all_expected_png_figures_exist",

    "all_expected_pdf_figures_exist"
  ),

  observed = c(

    as.character(
      N_H06_FAIL
    ),

    as.character(
      nrow(LONG)
    ),

    as.character(
      nrow(RW)
    ),

    as.character(
      nrow(DIRPAIR)
    ),

    as.character(
      nrow(BRPAIR)
    ),

    as.character(
      nrow(RADIUS)
    ),

    as.character(
      nrow(JOINT)
    ),

    as.character(
      nrow(WCOR)
    ),

    as.character(
      STD_EST_MISMATCH
    ),

    as.character(
      STD_CI_LOW_MISMATCH
    ),

    as.character(
      STD_CI_HIGH_MISMATCH
    ),

    as.character(
      nrow(TABLE1)
    ),

    as.character(
      N_TABLE1_DUP
    ),

    as.character(
      nrow(TABLE2)
    ),

    as.character(
      N_TABLE2_DUP
    ),

    as.character(
      N_RW_COMPARABLE
    ),

    as.character(
      N_RW_SAME
    ),

    as.character(
      N_RW_OPPOSITE
    ),

    as.character(
      sum(
        file.exists(
          EXPECTED_PNG
        )
      )
    ),

    as.character(
      sum(
        file.exists(
          EXPECTED_PDF
        )
      )
    )
  ),

  expected = c(

    "0",

    "1440",

    "480",

    "720",

    "720",

    "144",

    "480",

    "240",

    "0",

    "0",

    "0",

    "480",

    "0",

    "240",

    "0",

    "426",

    "401",

    "25",

    "6",

    "6"
  )
)


HARD07[
  ,
  status :=
    fifelse(
      observed ==
        expected,
      "PASS",
      "FAIL"
    )
]


fwrite(
  HARD07,
  OUT_HARD
)


N_HARD_FAIL <- HARD07[
  status ==
    "FAIL",
  .N
]


# ============================================================
# 34. Manuscript-oriented report
# ============================================================

REPORT <- c(

  "# 07H07 论文结果表与图形准备",

  "",

  paste0(
    "- Generated: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),

  "",

  "## 1. Analysis status",

  "",

  "- No new statistical model was fitted in 07H07.",

  "- 07H07 reads only the PASS results from 07H06.",

  "- All six fuels and all ten prespecified cumulative radii are retained.",

  "- No fuel or radius is selected according to statistical significance.",

  "",

  "## 2. Recommended main-paper results",

  "",

  "- Main exposure results: upwind Count-only and Heat-only models.",

  "- Primary effect-size display: change in gross beta per 1-SD exposure increase.",

  "- Spatial-definition sensitivity: Count versus Σ1/d and Heat versus ΣHeat/d.",

  "- Directional-control evidence: upwind versus downwind.",

  "- Temporal-scale robustness: quasi-daily versus filter-period.",

  "- Joint Count and Joint Heat are retained as conditional supplementary analyses and must be interpreted with VIF.",

  "",

  "## 3. Raw versus inverse-distance weighting",

  "",

  paste0(
    "- Comparable raw-weighted pairs: ",
    N_RW_COMPARABLE,
    "."
  ),

  paste0(
    "- Same effect direction: ",
    N_RW_SAME,
    "/",
    N_RW_COMPARABLE,
    " (",
    sprintf(
      "%.1f",
      PCT_RW_SAME
    ),
    "%)."
  ),

  paste0(
    "- Opposite effect direction: ",
    N_RW_OPPOSITE,
    "/",
    N_RW_COMPARABLE,
    "."
  ),

  paste0(
    "- Both raw and weighted 95% CI exclude zero: ",
    N_BOTH_CI,
    " pairs."
  ),

  paste0(
    "- Raw only 95% CI excludes zero: ",
    N_RAW_ONLY_CI,
    " pairs."
  ),

  paste0(
    "- Weighted only 95% CI excludes zero: ",
    N_WEIGHTED_ONLY_CI,
    " pairs."
  ),

  paste0(
    "- Neither 95% CI excludes zero among comparable pairs: ",
    N_NEITHER_CI,
    " pairs."
  ),

  paste0(
    "- Pair-incomplete (at least one model not estimable): ",
    N_PAIR_INCOMPLETE_CI,
    " pairs."
  ),

  "",

  "## 4. Directional-control example: Σ1/d",

  "",

  paste0(
    "- Upwind-only CI exclusion: ",
    WP_UP_ONLY,
    " radius/branch/fuel combinations."
  ),

  paste0(
    "- Downwind-only CI exclusion: ",
    WP_DOWN_ONLY,
    " combinations."
  ),

  paste0(
    "- Both upwind and downwind CI exclusion: ",
    WP_BOTH,
    " combinations."
  ),

  "- These are descriptive directional-control patterns, not formal tests of the difference between upwind and downwind coefficients.",

  "",

  "## 5. Main figures",

  "",

  "- Fig1: raw Count and Heat upwind results across all radii.",

  "- Fig2: Count versus Σ1/d spatial-definition sensitivity.",

  "- Fig3: Heat versus ΣHeat/d spatial-definition sensitivity.",

  "",

  "## 6. Supplementary figures",

  "",

  "- FigS1: upwind versus downwind directional-control patterns.",

  "- FigS2: quasi-daily versus filter-period consistency.",

  "- FigS3: Joint Count and Joint Heat conditional effects.",

  "",

  "## 7. Tables",

  "",

  "- Table 01: all upwind single-exposure results at all prespecified radii.",

  "- Table 02: all upwind Joint Count / Joint Heat results.",

  "- Table 03: raw-versus-distance-weighted robustness summary.",

  "- Table 04: directional negative-control summary.",

  "- Table 05: quasi-versus-filter-period consistency summary.",

  "- Table 06: model coverage and SKIP distribution.",

  "- Tables 07–09: radius patterns and collinearity diagnostics.",

  "",

  "## 8. Interpretation rules",

  "",

  "- Fixed-increment effects remain available for exposure-specific interpretation.",

  "- Cross-exposure comparisons should preferentially use 1-SD effects because Count, Heat, Σ1/d and ΣHeat/d have different units.",

  "- A CI excluding zero is reported descriptively; no best fuel or radius is selected from P values.",

  "- Joint Count/Heat estimates with elevated VIF should not be treated as clean independent source effects.",

  "- Association results do not by themselves establish atmospheric transport mechanism or causality.",

  "",

  if (
    N_HARD_FAIL ==
      0L
  ) {

    "**07H07 PAPER-RESULT PREPARATION: PASS**"

  } else {

    "**07H07 PAPER-RESULT PREPARATION: FAIL**"
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
# 35. Console
# ============================================================

cat(
  "\n============================================================\n",
  "A. PAPER RESULT PREPARATION\n",
  "============================================================\n",
  sep = ""
)


cat(
  "Upwind single-exposure rows = ",
  nrow(TABLE1),
  "\n",
  "Upwind joint rows = ",
  nrow(TABLE2),
  "\n",
  "Raw-weighted comparable pairs = ",
  N_RW_COMPARABLE,
  "\n",
  "Raw-weighted same sign = ",
  N_RW_SAME,
  " / ",
  N_RW_COMPARABLE,
  " (",
  sprintf(
    "%.1f",
    PCT_RW_SAME
  ),
  "%)\n",
  sep = ""
)


cat(
  "\n============================================================\n",
  "B. RAW vs WEIGHTED CI PATTERN\n",
  "============================================================\n",
  "Both CI exclude zero = ",
  N_BOTH_CI,
  "\n",
  "Raw only = ",
  N_RAW_ONLY_CI,
  "\n",
  "Weighted only = ",
  N_WEIGHTED_ONLY_CI,
  "\n",
  "Neither (both PASS) = ",
  N_NEITHER_CI,
  "\n",
  "Pair incomplete = ",
  N_PAIR_INCOMPLETE_CI,
  "\n",
  sep = ""
)


cat(
  "\n============================================================\n",
  "C. Σ1/d DIRECTIONAL PATTERN\n",
  "============================================================\n",
  "Upwind only = ",
  WP_UP_ONLY,
  "\n",
  "Downwind only = ",
  WP_DOWN_ONLY,
  "\n",
  "Both = ",
  WP_BOTH,
  "\n",
  sep = ""
)


cat(
  "\n============================================================\n",
  "D. HARD CHECKS\n",
  "============================================================\n",
  sep = ""
)


print(
  HARD07,
  nrows = Inf
)


cat(
  "\n============================================================\n"
)


if (
  N_HARD_FAIL ==
    0L
) {

  cat(
    "07H07 PAPER-RESULT PREPARATION: PASS\n"
  )

} else {

  cat(
    "07H07 PAPER-RESULT PREPARATION: FAIL\n"
  )
}


cat(
  "============================================================\n",
  "Paper tables directory:\n",
  TABLE_DIR,
  "\n\n",
  "Plot-data directory:\n",
  PLOT_DATA_DIR,
  "\n\n",
  "Figures directory:\n",
  FIG_DIR,
  "\n\n",
  "Hard FAIL = ",
  N_HARD_FAIL,
  "\n",
  sep = ""
)


if (
  N_HARD_FAIL >
    0L
) {

  stopf(
    "07H07存在hard failure。"
  )
}

