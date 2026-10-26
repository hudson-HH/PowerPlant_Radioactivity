#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999
)


# ============================================================
# 07H03
# 燃料特异 Count + Heat 模型结果审计与汇总
#
# 作用：
#   1. 读取已经完成的07H02模型结果
#   2. 不重新拟合任何模型
#   3. 不修改07H01/07H02
#   4. 不修改07A–07G冻结文件
#   5. 整理upwind主分析
#   6. 整理downwind方向负对照
#   7. 比较quasi-daily与filter-period
#   8. 汇总10个预设累计半径
#   9. 审计SKIP
#  10. 审计Count–Heat共线性
#
# 注意：
#   不根据P值选择“最佳半径”
#   不因为结果显著/不显著修改暴露定义
# ============================================================


# ============================================================
# 0. Package
# ============================================================

if (
  !requireNamespace(
    "data.table",
    quietly = TRUE
  )
) {

  stop(
    "缺少R包 data.table。当前项目此前已经使用该包；请先确认R环境。",
    call. = FALSE
  )
}


suppressPackageStartupMessages({

  library(
    data.table
  )
})


# ============================================================
# 1. Helpers
# ============================================================

stopf <- function(...) {

  stop(
    paste0(...),
    call. = FALSE
  )
}


fmt_int <- function(x) {

  format(
    x,
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
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
    length(
      miss
    ) >
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


assert_exact_set <- function(
  observed,
  expected,
  label
) {

  observed <- sort(
    unique(
      as.character(
        observed
      )
    )
  )

  expected <- sort(
    unique(
      as.character(
        expected
      )
    )
  )


  if (
    !identical(
      observed,
      expected
    )
  ) {

    stopf(
      label,
      "取值集合不一致。\nobserved = ",
      paste(
        observed,
        collapse = "|"
      ),
      "\nexpected = ",
      paste(
        expected,
        collapse = "|"
      )
    )
  }
}


collapse_radii <- function(
  radius,
  condition
) {

  keep <- condition %in%
    TRUE


  z <- radius[
    keep
  ]


  z <- sort(
    unique(
      z[
        !is.na(z)
      ]
    )
  )


  if (
    length(z) ==
      0L
  ) {

    return("")
  }


  paste(
    z,
    collapse = ";"
  )
}


effect_direction_fun <- function(x) {

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


ci_class <- function(
  status,
  low,
  high
) {

  fifelse(
    status !=
      "PASS",

    status,

    fifelse(
      is.na(low) |
        is.na(high),

      "CI_NOT_AVAILABLE",

      fifelse(
        low >
          0,

        "CI_ABOVE_ZERO",

        fifelse(
          high <
            0,

          "CI_BELOW_ZERO",

          "CI_INCLUDES_ZERO"
        )
      )
    )
  )
}


component_label <- function(
  model_type,
  term
) {

  fifelse(
    model_type ==
      "count_only" &
      term ==
        "count",

    "count_only",

    fifelse(
      model_type ==
        "heat_only" &
        term ==
          "heat",

      "heat_only",

      fifelse(
        model_type ==
          "joint_count_heat" &
          term ==
            "count",

        "joint_count",

        fifelse(
          model_type ==
            "joint_count_heat" &
            term ==
              "heat",

          "joint_heat",

          NA_character_
        )
      )
    )
  )
}


component_cn <- function(x) {

  fifelse(
    x ==
      "count_only",

    "Count-only：每增加1座实际运行CAMPD电厂",

    fifelse(
      x ==
        "heat_only",

      "Heat-only：每增加10,000 MMBtu Heat Input",

      fifelse(
        x ==
          "joint_count",

        "Joint Count：控制Heat后每增加1座实际运行CAMPD电厂",

        fifelse(
          x ==
            "joint_heat",

          "Joint Heat：控制Count后每增加10,000 MMBtu Heat Input",

          NA_character_
        )
      )
    )
  )
}


joint_qc_class <- function(
  model_type,
  vif
) {

  out <- rep(
    "NOT_JOINT_MODEL",
    length(vif)
  )


  idx <- model_type ==
    "joint_count_heat"


  out[
    idx &
      is.infinite(vif)
  ] <-
    "PERFECT_COLLINEARITY_VIF_INF"


  out[
    idx &
      is.finite(vif) &
      vif >=
        10
  ] <-
    "HIGH_VIF_GE10"


  out[
    idx &
      is.finite(vif) &
      vif >=
        5 &
      vif <
        10
  ] <-
    "MODERATE_VIF_5_TO_10"


  out[
    idx &
      is.finite(vif) &
      vif <
        5
  ] <-
    "NO_HIGH_VIF_FLAG"


  out[
    idx &
      is.na(vif)
  ] <-
    "VIF_NOT_ESTIMABLE"


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


H02 <- file.path(
  H_ROOT,
  "07H02_model"
)


OUT <- file.path(
  H_ROOT,
  "07H03_模型结果审计与汇总"
)


dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 3. Frozen 07H02 inputs
# ============================================================

FILE_EFFECT <- file.path(
  H02,
  "01_07H02_model_effects_all.csv"
)


FILE_STATUS <- file.path(
  H02,
  "02_07H02_model_status_all.csv"
)


FILE_COR <- file.path(
  H02,
  "03_07H02_count_heat_collinearity.csv"
)


FILE_HARD <- file.path(
  H02,
  "05_07H02_hard_checks.csv"
)


# ============================================================
# 4. Outputs
# ============================================================

OUT_MASTER <- file.path(
  OUT,
  "01_07H03_全部模型结果主表.csv"
)


OUT_UPWIND <- file.path(
  OUT,
  "02_07H03_upwind主分析结果.csv"
)


OUT_DOWNWIND <- file.path(
  OUT,
  "03_07H03_downwind方向负对照结果.csv"
)


OUT_SKIP <- file.path(
  OUT,
  "04_07H03_SKIP审计.csv"
)


OUT_COLLINEARITY <- file.path(
  OUT,
  "05_07H03_Count_Heat共线性审计.csv"
)


OUT_COVERAGE <- file.path(
  OUT,
  "06_07H03_模型覆盖情况.csv"
)


OUT_RADIUS <- file.path(
  OUT,
  "07_07H03_累计半径结果模式.csv"
)


OUT_DIRECTION <- file.path(
  OUT,
  "08_07H03_upwind_downwind配对比较.csv"
)


OUT_BRANCH <- file.path(
  OUT,
  "09_07H03_quasi_filter配对比较.csv"
)


OUT_KEY <- file.path(
  OUT,
  "10_07H03_重点结果索引.csv"
)


OUT_HARD <- file.path(
  OUT,
  "11_07H03_hard_checks.csv"
)


OUT_REPORT <- file.path(
  OUT,
  "12_07H03_模型结果审计与汇总报告.md"
)


OUT_SESSION <- file.path(
  OUT,
  "13_sessionInfo.txt"
)


# ============================================================
# 5. Frozen design constants
# ============================================================

FUELS <- c(
  "coal",
  "gas",
  "oil",
  "other_fossil",
  "biomass",
  "solid_waste"
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


BRANCHES <- c(
  "filter_period",
  "quasi_daily"
)


DIRECTIONS <- c(
  "upwind",
  "downwind"
)


MODEL_TYPES <- c(
  "count_only",
  "heat_only",
  "joint_count_heat"
)


EXPECTED_EFFECT_ROWS <- 960L

EXPECTED_MODEL_SPECS <- 720L

EXPECTED_PASS <- 846L

EXPECTED_SKIP <- 114L

EXPECTED_FAIL <- 0L

EXPECTED_COR_ROWS <- 240L

EXPECTED_VIF_INF <- 6L

EXPECTED_VIF_GE10 <- 14L

EXPECTED_VIF_5_10 <- 56L


# ============================================================
# 6. Input existence audit
# ============================================================

assert_file(
  FILE_EFFECT,
  "07H02 effect结果"
)


assert_file(
  FILE_STATUS,
  "07H02 model status"
)


assert_file(
  FILE_COR,
  "07H02 Count–Heat correlation"
)


assert_file(
  FILE_HARD,
  "07H02 hard checks"
)


cat(
  "\n============================================================\n",
  "07H03 INPUT PATH AUDIT\n",
  "============================================================\n",
  "ROOT = ",
  ROOT,
  "\n",
  "07H02 = ",
  H02,
  "\n",
  "07H03 output = ",
  OUT,
  "\n",
  sep = ""
)


# ============================================================
# 7. Read frozen 07H02 outputs
# ============================================================

E <- fread(
  FILE_EFFECT,
  showProgress = FALSE
)


S <- fread(
  FILE_STATUS,
  showProgress = FALSE
)


C <- fread(
  FILE_COR,
  showProgress = FALSE
)


H <- fread(
  FILE_HARD,
  showProgress = FALSE
)


# ============================================================
# 8. Schema audit
# ============================================================

E_REQUIRED <- c(
  "spec_id",
  "branch",
  "fuel",
  "direction",
  "radius_km",
  "model_type",
  "term",
  "status",
  "n_rows",
  "n_monitors",
  "estimate_mBq_m3",
  "se_mBq_m3",
  "ci95_low_mBq_m3",
  "ci95_high_mBq_m3",
  "p_wald",
  "note"
)


S_REQUIRED <- c(
  "spec_id",
  "branch",
  "model_type",
  "singular_any",
  "convergence_message_any",
  "warnings_any"
)


C_REQUIRED <- c(
  "branch",
  "fuel",
  "direction",
  "radius_km",
  "n_both_finite",
  "pearson_r",
  "spearman_rho",
  "two_predictor_vif",
  "collinearity_flag"
)


H_REQUIRED <- c(
  "check",
  "observed",
  "expected",
  "status"
)


assert_cols(
  E,
  E_REQUIRED,
  "01_07H02_model_effects_all.csv"
)


assert_cols(
  S,
  S_REQUIRED,
  "02_07H02_model_status_all.csv"
)


assert_cols(
  C,
  C_REQUIRED,
  "03_07H02_count_heat_collinearity.csv"
)


assert_cols(
  H,
  H_REQUIRED,
  "05_07H02_hard_checks.csv"
)


# ============================================================
# 9. Frozen 07H02 hard-check audit
# ============================================================

EXPECTED_07H02_HARD_NAMES <- c(
  "07H01_hard_checks_all_pass",
  "07H01_yearly_files_24",
  "07F_quasi_primary_rows_426657",
  "07F_quasi_primary_monitors_133",
  "07F_filter_formal_rows_107407",
  "07F_filter_formal_monitors_133",
  "filter_period_universe_174200",
  "model_specs_720",
  "effect_rows_960",
  "duplicate_effect_keys_zero",
  "count_heat_correlation_rows_240",
  "fit_error_rows_zero"
)


missing_hard_names <- setdiff(
  EXPECTED_07H02_HARD_NAMES,
  H$check
)


if (
  length(
    missing_hard_names
  ) >
    0L
) {

  stopf(
    "07H02 hard-check表缺少正式检查项：",
    paste(
      missing_hard_names,
      collapse = ", "
    )
  )
}


H_FORMAL <- H[
  check %in%
    EXPECTED_07H02_HARD_NAMES
]


if (
  H_FORMAL[
    status !=
      "PASS",
    .N
  ] >
    0L
) {

  print(
    H_FORMAL[
      status !=
        "PASS"
    ]
  )


  stopf(
    "07H02存在未PASS的正式hard check，禁止继续07H03。"
  )
}


# ============================================================
# 10. Design-value audit
# ============================================================

E[
  ,
  radius_km :=
    as.integer(
      radius_km
    )
]


C[
  ,
  radius_km :=
    as.integer(
      radius_km
    )
]


assert_exact_set(
  E$branch,
  BRANCHES,
  "branch"
)


assert_exact_set(
  E$fuel,
  FUELS,
  "fuel"
)


assert_exact_set(
  E$direction,
  DIRECTIONS,
  "direction"
)


assert_exact_set(
  E$radius_km,
  RADII,
  "radius_km"
)


assert_exact_set(
  E$model_type,
  MODEL_TYPES,
  "model_type"
)


VALID_TERM_COMBO <- E[
  !(
    (
      model_type ==
        "count_only" &
      term ==
        "count"
    ) |
    (
      model_type ==
        "heat_only" &
      term ==
        "heat"
    ) |
    (
      model_type ==
        "joint_count_heat" &
      term %in%
        c(
          "count",
          "heat"
        )
    )
  ),
  .N
]


if (
  VALID_TERM_COMBO >
    0L
) {

  stopf(
    "发现非法的model_type × term组合：",
    VALID_TERM_COMBO
  )
}


# ============================================================
# 11. Technical-status audit
# ============================================================

if (
  nrow(S) !=
    EXPECTED_MODEL_SPECS
) {

  stopf(
    "07H02 model-status行数不是720；实际=",
    nrow(S)
  )
}


if (
  uniqueN(
    S$spec_id
  ) !=
    EXPECTED_MODEL_SPECS
) {

  stopf(
    "07H02 model-status唯一spec_id不是720；实际=",
    uniqueN(
      S$spec_id
    )
  )
}


N_SINGULAR <- S[
  singular_any %in%
    TRUE,
  .N
]


N_CONVERGENCE <- S[
  convergence_message_any %in%
    TRUE,
  .N
]


N_WARNING <- S[
  warnings_any %in%
    TRUE,
  .N
]


if (
  N_SINGULAR >
    0L |
  N_CONVERGENCE >
    0L |
  N_WARNING >
    0L
) {

  stopf(
    "07H02技术状态异常：",
    " singular=",
    N_SINGULAR,
    "; convergence=",
    N_CONVERGENCE,
    "; warnings=",
    N_WARNING
  )
}


# ============================================================
# 12. Frozen result dimensions
# ============================================================

if (
  nrow(E) !=
    EXPECTED_EFFECT_ROWS
) {

  stopf(
    "effect rows不是960；实际=",
    nrow(E)
  )
}


if (
  uniqueN(
    E$spec_id
  ) !=
    EXPECTED_MODEL_SPECS
) {

  stopf(
    "effect表中的唯一model specs不是720；实际=",
    uniqueN(
      E$spec_id
    )
  )
}


N_PASS <- E[
  status ==
    "PASS",
  .N
]


N_SKIP <- E[
  status ==
    "SKIP",
  .N
]


N_FAIL <- E[
  status ==
    "FAIL",
  .N
]


if (
  N_PASS !=
    EXPECTED_PASS
) {

  stopf(
    "PASS effect rows不是846；实际=",
    N_PASS
  )
}


if (
  N_SKIP !=
    EXPECTED_SKIP
) {

  stopf(
    "SKIP effect rows不是114；实际=",
    N_SKIP
  )
}


if (
  N_FAIL !=
    EXPECTED_FAIL
) {

  stopf(
    "FAIL effect rows不是0；实际=",
    N_FAIL
  )
}


E_DUP <- E[
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
  N >
    1L,
  .N
]


if (
  E_DUP >
    0L
) {

  stopf(
    "07H02 effect结果存在重复主键；重复数=",
    E_DUP
  )
}


# ============================================================
# 13. Count–Heat collinearity audit
# ============================================================

if (
  nrow(C) !=
    EXPECTED_COR_ROWS
) {

  stopf(
    "Count–Heat correlation rows不是240；实际=",
    nrow(C)
  )
}


C_DUP <- C[
  ,
  .N,
  by = .(
    branch,
    fuel,
    direction,
    radius_km
  )
][
  N >
    1L,
  .N
]


if (
  C_DUP >
    0L
) {

  stopf(
    "Count–Heat correlation表存在重复主键；重复数=",
    C_DUP
  )
}


C[
  ,
  `:=`(

    pearson_r =
      safe_num(
        pearson_r
      ),

    spearman_rho =
      safe_num(
        spearman_rho
      ),

    two_predictor_vif =
      safe_num(
        two_predictor_vif
      )
  )
]


C[
  ,
  collinearity_flag_07H03 :=
    fifelse(

      is.infinite(
        two_predictor_vif
      ),

      "PERFECT_COLLINEARITY_VIF_INF",

      fifelse(

        is.finite(
          two_predictor_vif
        ) &
          two_predictor_vif >=
            10,

        "HIGH_VIF_GE10",

        fifelse(

          is.finite(
            two_predictor_vif
          ) &
            two_predictor_vif >=
              5,

          "MODERATE_VIF_5_TO_10",

          fifelse(

            is.finite(
              two_predictor_vif
            ),

            "NO_HIGH_VIF_FLAG",

            "VIF_NOT_ESTIMABLE"
          )
        )
      )
    )
]


N_LABEL_MISMATCH <- C[
  as.character(
    collinearity_flag
  ) !=
    collinearity_flag_07H03,
  .N
]


if (
  N_LABEL_MISMATCH >
    0L
) {

  print(
    C[
      as.character(
        collinearity_flag
      ) !=
        collinearity_flag_07H03,
      .(
        branch,
        fuel,
        direction,
        radius_km,
        two_predictor_vif,
        collinearity_flag,
        collinearity_flag_07H03
      )
    ]
  )


  stopf(
    "07H02共线性标签与07H03独立重算不一致；",
    "mismatch=",
    N_LABEL_MISMATCH
  )
}


N_INF <- C[
  is.infinite(
    two_predictor_vif
  ),
  .N
]


N_GE10 <- C[
  is.finite(
    two_predictor_vif
  ) &
    two_predictor_vif >=
      10,
  .N
]


N_5_10 <- C[
  is.finite(
    two_predictor_vif
  ) &
    two_predictor_vif >=
      5 &
    two_predictor_vif <
      10,
  .N
]


if (
  N_INF !=
    EXPECTED_VIF_INF
) {

  stopf(
    "VIF=Inf数量不是6；实际=",
    N_INF
  )
}


if (
  N_GE10 !=
    EXPECTED_VIF_GE10
) {

  stopf(
    "有限VIF>=10数量不是14；实际=",
    N_GE10
  )
}


if (
  N_5_10 !=
    EXPECTED_VIF_5_10
) {

  stopf(
    "5<=VIF<10数量不是56；实际=",
    N_5_10
  )
}


cat(
  "\n============================================================\n",
  "07H02 FROZEN INPUT AUDIT: PASS\n",
  "============================================================\n",
  "Effect rows = ",
  nrow(E),
  "\n",
  "Model specs = ",
  uniqueN(E$spec_id),
  "\n",
  "PASS = ",
  N_PASS,
  "\n",
  "SKIP = ",
  N_SKIP,
  "\n",
  "FAIL = ",
  N_FAIL,
  "\n",
  "Singular specs = ",
  N_SINGULAR,
  "\n",
  "Convergence-message specs = ",
  N_CONVERGENCE,
  "\n",
  "Warning specs = ",
  N_WARNING,
  "\n",
  "Correlation rows = ",
  nrow(C),
  "\n",
  "VIF=Inf = ",
  N_INF,
  "\n",
  "Finite VIF>=10 = ",
  N_GE10,
  "\n",
  "5<=VIF<10 = ",
  N_5_10,
  "\n",
  "Collinearity label mismatch = ",
  N_LABEL_MISMATCH,
  "\n",
  sep = ""
)


# ============================================================
# 14. Normalize numeric effect fields
# ============================================================

NUM_E <- c(
  "n_rows",
  "n_monitors",
  "n_unique_exposure",
  "scale_value",
  "raw_increment",
  "wind_speed_ms",
  "estimate_mBq_m3",
  "se_mBq_m3",
  "ci95_low_mBq_m3",
  "ci95_high_mBq_m3",
  "p_wald",
  "sd_raw_exposure",
  "estimate_per_1SD_mBq_m3"
)


for (
  cc in intersect(
    NUM_E,
    names(E)
  )
) {

  E[
    ,
    (cc) :=
      safe_num(
        get(cc)
      )
  ]
}


NUM_C_OPTIONAL <- c(
  "n_model_base",
  "n_count_finite",
  "n_heat_finite",
  "n_both_finite",
  "n_positive_count",
  "n_positive_heat",
  "count_zero_pct",
  "heat_zero_pct"
)


for (
  cc in intersect(
    NUM_C_OPTIONAL,
    names(C)
  )
) {

  C[
    ,
    (cc) :=
      safe_num(
        get(cc)
      )
  ]
}


# ============================================================
# 15. Canonical effect-component labels
# ============================================================

E[
  ,
  result_component :=
    component_label(
      model_type,
      term
    )
]


if (
  E[
    is.na(
      result_component
    ),
    .N
  ] >
    0L
) {

  stopf(
    "存在无法识别的model_type × term组合。"
  )
}


E[
  ,
  `:=`(

    result_component_cn =
      component_cn(
        result_component
      ),

    effect_direction =
      effect_direction_fun(
        estimate_mBq_m3
      ),

    ci_status =
      ci_class(
        status,
        ci95_low_mBq_m3,
        ci95_high_mBq_m3
      ),

    p_lt_0_05 =
      status ==
        "PASS" &
      is.finite(
        p_wald
      ) &
      p_wald <
        0.05
  )
]


# ============================================================
# ============================================================
# 16. Merge Count–Heat correlation metadata
#
# IMPORTANT:
# Pearson/Spearman/VIF同时为NA，并不代表merge失败。
# 对暴露变化不足的稀疏fuel组合，这些统计量本来就可能为NA。
# 因此必须使用显式主键匹配标记，而不能根据统计量NA判断。
# ============================================================

C_JOIN_SOURCE <- copy(C)


C_JOIN_SOURCE[
  ,
  .cor_key_present := TRUE
]


COR_KEEP <- c(
  "branch",
  "fuel",
  "direction",
  "radius_km",
  "n_both_finite",
  "count_zero_pct",
  "heat_zero_pct",
  "pearson_r",
  "spearman_rho",
  "two_predictor_vif",
  "collinearity_flag_07H03",
  ".cor_key_present"
)


COR_KEEP <- intersect(
  COR_KEEP,
  names(C_JOIN_SOURCE)
)


COR_JOIN <- C_JOIN_SOURCE[
  ,
  ..COR_KEEP
]


# ------------------------------------------------------------
# Correlation-key completeness before merge
# ------------------------------------------------------------

if (
  nrow(COR_JOIN) !=
    EXPECTED_COR_ROWS
) {

  stopf(
    "COR_JOIN不是240行；实际=",
    nrow(COR_JOIN)
  )
}


COR_JOIN_DUP <- COR_JOIN[
  ,
  .N,
  by = .(
    branch,
    fuel,
    direction,
    radius_km
  )
][
  N >
    1L,
  .N
]


if (
  COR_JOIN_DUP >
    0L
) {

  stopf(
    "COR_JOIN存在重复主键；重复数=",
    COR_JOIN_DUP
  )
}


# ------------------------------------------------------------
# Merge
# ------------------------------------------------------------

E <- merge(
  E,
  COR_JOIN,

  by = c(
    "branch",
    "fuel",
    "direction",
    "radius_km"
  ),

  all.x = TRUE,
  sort = FALSE
)


if (
  nrow(E) !=
    EXPECTED_EFFECT_ROWS
) {

  stopf(
    "合并共线性信息后effect row数发生变化；实际=",
    nrow(E)
  )
}


# ------------------------------------------------------------
# True merge-key audit
#
# Do NOT use Pearson/Spearman/VIF missingness here.
# Some sparse combinations legitimately have all three = NA.
# ------------------------------------------------------------

N_COR_UNMATCHED <- E[
  is.na(
    .cor_key_present
  ),
  .N
]


if (
  N_COR_UNMATCHED >
    0L
) {

  cat(
    "\n未匹配到Count–Heat correlation的真实主键：\n"
  )


  print(
    unique(
      E[
        is.na(
          .cor_key_present
        ),
        .(
          branch,
          fuel,
          direction,
          radius_km
        )
      ]
    ),
    nrows = Inf
  )


  stopf(
    "存在effect结果没有匹配到Count–Heat共线性主键；真实未匹配行数=",
    N_COR_UNMATCHED
  )
}


cat(
  "Count–Heat correlation key merge: PASS\n",
  "True unmatched effect rows = ",
  N_COR_UNMATCHED,
  "\n",
  sep = ""
)


# 主键审计完成后删除临时标记，不进入正式输出。
E[
  ,
  .cor_key_present :=
    NULL
]


# 合法情况：
# 某些稀疏fuel组合即使成功匹配，
# pearson_r / spearman_rho / two_predictor_vif也可能全部为NA。
# 这些组合保留，由VIF_NOT_ESTIMABLE和SKIP状态共同解释。


E[
  ,
  joint_interpretability :=
    joint_qc_class(
      model_type,
      two_predictor_vif
    )
]


setorder(
  E,
  branch,
  fuel,
  direction,
  radius_km,
  model_type,
  term
)


# 17. Master result table
# ============================================================

MASTER_COLS <- c(

  "branch",
  "fuel",
  "direction",
  "radius_km",

  "model_type",
  "term",

  "result_component",
  "result_component_cn",

  "status",

  "n_rows",
  "n_monitors",
  "n_unique_exposure",

  "effect_label",
  "raw_increment",
  "raw_unit",

  "wind_speed_ms",

  "estimate_mBq_m3",
  "se_mBq_m3",
  "ci95_low_mBq_m3",
  "ci95_high_mBq_m3",
  "p_wald",

  "effect_direction",
  "ci_status",
  "p_lt_0_05",

  "sd_raw_exposure",
  "estimate_per_1SD_mBq_m3",

  "pearson_r",
  "spearman_rho",
  "two_predictor_vif",
  "collinearity_flag_07H03",
  "joint_interpretability",

  "count_zero_pct",
  "heat_zero_pct",
  "n_both_finite",

  "singular",
  "convergence_message",
  "warnings",

  "note",
  "spec_id"
)


MASTER_COLS <- intersect(
  MASTER_COLS,
  names(E)
)


MASTER <- E[
  ,
  ..MASTER_COLS
]


MASTER_DUP <- MASTER[
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
  N >
    1L,
  .N
]


if (
  MASTER_DUP >
    0L
) {

  stopf(
    "07H03 MASTER存在重复主键。"
  )
}


fwrite(
  MASTER,
  OUT_MASTER
)


# ============================================================
# 18. Upwind main analysis
# ============================================================

UP <- MASTER[
  direction ==
    "upwind"
]


setorder(
  UP,
  branch,
  fuel,
  radius_km,
  model_type,
  term
)


fwrite(
  UP,
  OUT_UPWIND
)


# ============================================================
# 19. Downwind negative-control analysis
# ============================================================

DOWN <- MASTER[
  direction ==
    "downwind"
]


setorder(
  DOWN,
  branch,
  fuel,
  radius_km,
  model_type,
  term
)


fwrite(
  DOWN,
  OUT_DOWNWIND
)


# ============================================================
# 20. SKIP audit
# ============================================================

SKIP <- MASTER[
  status ==
    "SKIP"
]


SKIP[
  ,
  skip_reason_group :=
    fifelse(

      grepl(
        "Insufficient rows/monitors/joint exposure variation",
        note,
        fixed = TRUE
      ),

      "JOINT_INSUFFICIENT_VARIATION",

      fifelse(

        grepl(
          "Insufficient rows/monitors/exposure variation",
          note,
          fixed = TRUE
        ),

        "SINGLE_EXPOSURE_INSUFFICIENT_VARIATION",

        fifelse(

          grepl(
            "rank",
            note,
            ignore.case = TRUE
          ),

          "MODEL_MATRIX_RANK_LIMITATION",

          "OTHER"
        )
      )
    )
]


setorder(
  SKIP,
  branch,
  fuel,
  direction,
  radius_km,
  model_type,
  term
)


fwrite(
  SKIP,
  OUT_SKIP
)


# ============================================================
# 21. Count–Heat collinearity table
# ============================================================

COL <- copy(C)


setorder(
  COL,
  -two_predictor_vif,
  branch,
  fuel,
  direction,
  radius_km
)


fwrite(
  COL,
  OUT_COLLINEARITY
)


# ============================================================
# 22. Model coverage
# ============================================================

COVERAGE <- MASTER[
  ,
  .(

    n_effect_rows =
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

    n_ci_includes_zero =
      sum(
        ci_status ==
          "CI_INCLUDES_ZERO",
        na.rm = TRUE
      )
  ),

  by = .(
    branch,
    fuel,
    direction,
    result_component
  )
]


setorder(
  COVERAGE,
  branch,
  fuel,
  direction,
  result_component
)


fwrite(
  COVERAGE,
  OUT_COVERAGE
)


# ============================================================
# 23. Radius-pattern summary
#
# 不选择最佳半径。
# 只描述所有预设累计半径上的结果模式。
# ============================================================

RADIUS_SUMMARY <- MASTER[
  ,
  .(

    n_radii_total =
      uniqueN(
        radius_km
      ),

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

    n_positive =
      sum(
        status ==
          "PASS" &
        estimate_mBq_m3 >
          0,
        na.rm = TRUE
      ),

    n_negative =
      sum(
        status ==
          "PASS" &
        estimate_mBq_m3 <
          0,
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

    n_p_lt_0_05 =
      sum(
        p_lt_0_05 %in%
          TRUE
      ),

    positive_radii_km =
      collapse_radii(
        radius_km,
        status ==
          "PASS" &
        estimate_mBq_m3 >
          0
      ),

    negative_radii_km =
      collapse_radii(
        radius_km,
        status ==
          "PASS" &
        estimate_mBq_m3 <
          0
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

    p_lt_0_05_radii_km =
      collapse_radii(
        radius_km,
        p_lt_0_05
      ),

    all_estimable_effects_same_sign =
      {

        zz <- sign(
          estimate_mBq_m3[
            status ==
              "PASS" &
            is.finite(
              estimate_mBq_m3
            )
          ]
        )


        zz <- zz[
          zz !=
            0
        ]


        if (
          length(zz) ==
            0L
        ) {

          NA

        } else {

          uniqueN(
            zz
          ) ==
            1L
        }
      }
  ),

  by = .(
    branch,
    fuel,
    direction,
    result_component
  )
]


setorder(
  RADIUS_SUMMARY,
  branch,
  fuel,
  direction,
  result_component
)


fwrite(
  RADIUS_SUMMARY,
  OUT_RADIUS
)


# ============================================================
# 24. Upwind vs downwind paired comparison
#
# 仅做描述性方向负对照比较。
# 不创建winner，不选择最佳方向。
# ============================================================

PAIR_BASE <- MASTER[
  status ==
    "PASS",
  .(
    branch,
    fuel,
    direction,
    radius_km,
    result_component,

    estimate =
      estimate_mBq_m3,

    low =
      ci95_low_mBq_m3,

    high =
      ci95_high_mBq_m3,

    p =
      p_wald,

    ci_status,

    joint_interpretability
  )
]


UP_PAIR <- PAIR_BASE[
  direction ==
    "upwind"
]


DOWN_PAIR <- PAIR_BASE[
  direction ==
    "downwind"
]


setnames(
  UP_PAIR,
  c(
    "estimate",
    "low",
    "high",
    "p",
    "ci_status",
    "joint_interpretability"
  ),
  c(
    "upwind_estimate",
    "upwind_ci_low",
    "upwind_ci_high",
    "upwind_p",
    "upwind_ci_status",
    "upwind_joint_interpretability"
  )
)


UP_PAIR[
  ,
  direction :=
    NULL
]


setnames(
  DOWN_PAIR,
  c(
    "estimate",
    "low",
    "high",
    "p",
    "ci_status",
    "joint_interpretability"
  ),
  c(
    "downwind_estimate",
    "downwind_ci_low",
    "downwind_ci_high",
    "downwind_p",
    "downwind_ci_status",
    "downwind_joint_interpretability"
  )
)


DOWN_PAIR[
  ,
  direction :=
    NULL
]


DIR_COMPARE <- merge(
  UP_PAIR,
  DOWN_PAIR,

  by = c(
    "branch",
    "fuel",
    "radius_km",
    "result_component"
  ),

  all = TRUE,
  sort = FALSE
)


DIR_COMPARE[
  ,
  `:=`(

    both_available =
      is.finite(
        upwind_estimate
      ) &
      is.finite(
        downwind_estimate
      ),

    same_effect_sign =
      fifelse(
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

    upwind_ci_excludes_zero =
      upwind_ci_status %in%
        c(
          "CI_ABOVE_ZERO",
          "CI_BELOW_ZERO"
        ),

    downwind_ci_excludes_zero =
      downwind_ci_status %in%
        c(
          "CI_ABOVE_ZERO",
          "CI_BELOW_ZERO"
        )
  )
]


setorder(
  DIR_COMPARE,
  branch,
  fuel,
  radius_km,
  result_component
)


fwrite(
  DIR_COMPARE,
  OUT_DIRECTION
)


# ============================================================
# 25. Quasi-daily vs filter-period paired comparison
# ============================================================

Q_PAIR <- MASTER[
  branch ==
    "quasi_daily" &
  status ==
    "PASS",
  .(
    fuel,
    direction,
    radius_km,
    result_component,

    quasi_estimate =
      estimate_mBq_m3,

    quasi_ci_low =
      ci95_low_mBq_m3,

    quasi_ci_high =
      ci95_high_mBq_m3,

    quasi_p =
      p_wald,

    quasi_ci_status =
      ci_status,

    quasi_joint_interpretability =
      joint_interpretability
  )
]


F_PAIR <- MASTER[
  branch ==
    "filter_period" &
  status ==
    "PASS",
  .(
    fuel,
    direction,
    radius_km,
    result_component,

    filter_estimate =
      estimate_mBq_m3,

    filter_ci_low =
      ci95_low_mBq_m3,

    filter_ci_high =
      ci95_high_mBq_m3,

    filter_p =
      p_wald,

    filter_ci_status =
      ci_status,

    filter_joint_interpretability =
      joint_interpretability
  )
]


BRANCH_COMPARE <- merge(
  Q_PAIR,
  F_PAIR,

  by = c(
    "fuel",
    "direction",
    "radius_km",
    "result_component"
  ),

  all = TRUE,
  sort = FALSE
)


BRANCH_COMPARE[
  ,
  `:=`(

    both_available =
      is.finite(
        quasi_estimate
      ) &
      is.finite(
        filter_estimate
      ),

    same_effect_sign =
      fifelse(
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

    quasi_ci_excludes_zero =
      quasi_ci_status %in%
        c(
          "CI_ABOVE_ZERO",
          "CI_BELOW_ZERO"
        ),

    filter_ci_excludes_zero =
      filter_ci_status %in%
        c(
          "CI_ABOVE_ZERO",
          "CI_BELOW_ZERO"
        )
  )
]


setorder(
  BRANCH_COMPARE,
  fuel,
  direction,
  radius_km,
  result_component
)


fwrite(
  BRANCH_COMPARE,
  OUT_BRANCH
)


# ============================================================
# 26. Upwind key-result index
#
# 不是“显著性筛选表”。
# 六个fuel、两个branch、全部半径均保留。
# ============================================================

KEY <- MASTER[
  direction ==
    "upwind",
  .(
    branch,
    fuel,
    radius_km,

    result_component,
    result_component_cn,

    status,

    estimate_mBq_m3,
    ci95_low_mBq_m3,
    ci95_high_mBq_m3,
    p_wald,

    effect_direction,
    ci_status,

    joint_interpretability,

    pearson_r,
    two_predictor_vif,

    note
  )
]


FUEL_COMPLETENESS <- MASTER[
  ,
  .(
    n_skip =
      sum(
        status ==
          "SKIP"
      )
  ),
  by = fuel
]


FUEL_COMPLETENESS[
  ,
  fuel_data_completeness :=
    fifelse(
      n_skip ==
        0L,
      "ALL_07H02_EFFECT_ROWS_ESTIMABLE",
      "PARTIAL_07H02_EFFECT_ROWS_ESTIMABLE"
    )
]


KEY <- merge(
  KEY,
  FUEL_COMPLETENESS[
    ,
    .(
      fuel,
      fuel_data_completeness
    )
  ],
  by = "fuel",
  all.x = TRUE,
  sort = FALSE
)


setorder(
  KEY,
  fuel,
  branch,
  radius_km,
  result_component
)


fwrite(
  KEY,
  OUT_KEY
)


# ============================================================
# 27. 07H03 hard checks
# ============================================================

HARD <- data.table(

  check = c(

    "07H02_hard_checks_all_pass",

    "effect_rows_960",

    "model_specs_720",

    "pass_effect_rows_846",

    "skip_effect_rows_114",

    "fail_effect_rows_zero",

    "correlation_rows_240",

    "duplicate_master_keys_zero",

    "VIF_inf_count_6",

    "VIF_ge10_finite_count_14",

    "VIF_5_to_10_count_56",

    "upwind_and_downwind_present",

    "quasi_and_filter_present",

    "six_fuels_present",

    "ten_radii_present"
  ),

  observed = c(

    as.character(
      H_FORMAL[
        status !=
          "PASS",
        .N
      ]
    ),

    as.character(
      nrow(
        MASTER
      )
    ),

    as.character(
      uniqueN(
        MASTER$spec_id
      )
    ),

    as.character(
      N_PASS
    ),

    as.character(
      N_SKIP
    ),

    as.character(
      N_FAIL
    ),

    as.character(
      nrow(C)
    ),

    as.character(
      MASTER_DUP
    ),

    as.character(
      N_INF
    ),

    as.character(
      N_GE10
    ),

    as.character(
      N_5_10
    ),

    paste(
      sort(
        unique(
          MASTER$direction
        )
      ),
      collapse = "|"
    ),

    paste(
      sort(
        unique(
          MASTER$branch
        )
      ),
      collapse = "|"
    ),

    as.character(
      uniqueN(
        MASTER$fuel
      )
    ),

    as.character(
      uniqueN(
        MASTER$radius_km
      )
    )
  ),

  expected = c(

    "0",

    "960",

    "720",

    "846",

    "114",

    "0",

    "240",

    "0",

    "6",

    "14",

    "56",

    "downwind|upwind",

    "filter_period|quasi_daily",

    "6",

    "10"
  )
)


HARD[
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
  HARD,
  OUT_HARD
)


N_HARD_FAIL <- HARD[
  status ==
    "FAIL",
  .N
]


# ============================================================
# 28. Console summaries
# ============================================================

cat(
  "\n============================================================\n",
  "A. Effect status\n",
  "============================================================\n"
)


print(
  MASTER[
    ,
    .N,
    by = status
  ][
    order(status)
  ]
)


cat(
  "\n============================================================\n",
  "B. Fuel × branch coverage\n",
  "============================================================\n"
)


print(
  MASTER[
    ,
    .(

      PASS =
        sum(
          status ==
            "PASS"
        ),

      SKIP =
        sum(
          status ==
            "SKIP"
        ),

      FAIL =
        sum(
          status ==
            "FAIL"
        )
    ),

    by = .(
      branch,
      fuel
    )
  ],
  nrows = Inf
)


cat(
  "\n============================================================\n",
  "C. Count–Heat collinearity\n",
  "============================================================\n",
  "VIF=Inf: ",
  N_INF,
  "\n",
  "Finite VIF>=10: ",
  N_GE10,
  "\n",
  "5<=VIF<10: ",
  N_5_10,
  "\n",
  sep = ""
)


cat(
  "\n============================================================\n",
  "D. Upwind radius-pattern summary\n",
  "============================================================\n"
)


print(
  RADIUS_SUMMARY[
    direction ==
      "upwind"
  ],
  nrows = Inf
)


cat(
  "\n============================================================\n",
  "E. Hard checks\n",
  "============================================================\n"
)


print(
  HARD,
  nrows = Inf
)


# ============================================================
# 29. Markdown report
# ============================================================

COMPLETE_FUELS <- FUEL_COMPLETENESS[
  n_skip ==
    0L,
  fuel
]


PARTIAL_FUELS <- FUEL_COMPLETENESS[
  n_skip >
    0L,
  fuel
]


REPORT <- c(

  "# 07H03 燃料特异 Count + Heat 模型结果审计与汇总",

  "",

  paste0(
    "- Generated: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),

  paste0(
    "- Input directory: ",
    H02
  ),

  "",

  "## 1. 冻结输入状态",

  "",

  paste0(
    "- 07H02 effect rows: ",
    nrow(MASTER)
  ),

  paste0(
    "- Model specs: ",
    uniqueN(
      MASTER$spec_id
    )
  ),

  paste0(
    "- PASS effect rows: ",
    N_PASS
  ),

  paste0(
    "- SKIP effect rows: ",
    N_SKIP
  ),

  paste0(
    "- FAIL effect rows: ",
    N_FAIL
  ),

  paste0(
    "- Singular model specs: ",
    N_SINGULAR
  ),

  paste0(
    "- Convergence-message specs: ",
    N_CONVERGENCE
  ),

  paste0(
    "- Warning specs: ",
    N_WARNING
  ),

  "",

  "## 2. 正式分析框架",

  "",

  "- Upwind为预设主方向。",
  "- Downwind为方向性负对照。",
  "- 同时保留quasi-daily和filter-period两个冻结分支。",
  "- 累计半径固定为20、25、30、35、40、45、50、100、150、200 km。",
  "- 不根据P值选择最佳半径。",
  "- Count-only、Heat-only、Joint Count和Joint Heat全部保留。",
  "- Joint Count表示控制Heat后的plant-count条件关联。",
  "- Joint Heat表示控制Count后的Heat条件关联。",

  "",

  "## 3. Count–Heat共线性",

  "",

  paste0(
    "- Perfect collinearity (VIF=Inf): ",
    N_INF,
    " combinations."
  ),

  paste0(
    "- High collinearity (finite VIF>=10): ",
    N_GE10,
    " combinations."
  ),

  paste0(
    "- Moderate collinearity (5<=VIF<10): ",
    N_5_10,
    " combinations."
  ),

  paste0(
    "- 07H02/07H03 collinearity-label mismatch: ",
    N_LABEL_MISMATCH,
    "."
  ),

  "- 高/完全共线性组合中的joint条件效应必须谨慎解释。",
  "- 不因为共线性而修改燃料分类、距离半径或暴露定义。",

  "",

  "## 4. Fuel模型完整性",

  "",

  paste0(
    "- Zero-SKIP fuels: ",
    paste(
      COMPLETE_FUELS,
      collapse = ", "
    )
  ),

  paste0(
    "- Partial-estimability fuels: ",
    paste(
      PARTIAL_FUELS,
      collapse = ", "
    )
  ),

  "- SKIP保留为正式可估计性结果，不强制补拟合。",

  "",

  "## 5. 输出",

  "",

  "- 01：全部模型结果主表。",
  "- 02：upwind主分析结果。",
  "- 03：downwind方向负对照。",
  "- 04：全部SKIP及原因。",
  "- 05：Count–Heat共线性审计。",
  "- 06：模型覆盖情况。",
  "- 07：全部预设累计半径上的结果模式。",
  "- 08：upwind与downwind成对比较。",
  "- 09：quasi-daily与filter-period成对比较。",
  "- 10：upwind重点结果索引；不是显著性筛选。",
  "- 11：07H03 hard checks。",

  "",

  if (
    N_HARD_FAIL ==
      0L
  ) {

    "**07H03 RESULT AUDIT: PASS**"

  } else {

    "**07H03 RESULT AUDIT: FAIL**"
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
# 30. Finish
# ============================================================

cat(
  "\n============================================================\n"
)


if (
  N_HARD_FAIL ==
    0L
) {

  cat(
    "07H03 RESULT AUDIT: PASS\n"
  )

} else {

  cat(
    "07H03 RESULT AUDIT: FAIL\n"
  )
}


cat(
  "============================================================\n",
  "Master effect rows = ",
  nrow(MASTER),
  "\n",
  "Model specs = ",
  uniqueN(
    MASTER$spec_id
  ),
  "\n",
  "PASS = ",
  N_PASS,
  "\n",
  "SKIP = ",
  N_SKIP,
  "\n",
  "FAIL = ",
  N_FAIL,
  "\n",
  "VIF=Inf = ",
  N_INF,
  "\n",
  "Finite VIF>=10 = ",
  N_GE10,
  "\n",
  "5<=VIF<10 = ",
  N_5_10,
  "\n",
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
    "07H03存在hard failure。"
  )
}

