#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999
)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("缺少R包 data.table。", call. = FALSE)
}

if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("缺少R包 lme4。", call. = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})


# ============================================================
# 07H05
# 距离加权敏感性统计模型
#
# 目的：
#   与冻结07H02主分析进行空间暴露定义敏感性比较。
#
# 主分析（冻结07H02）：
#   Count
#   Heat
#
# 本阶段敏感性暴露：
#   1. Σ1/d
#   2. ΣHeat/d
#
# 只拟合单暴露模型：
#   proximity only
#   heat-proximity only
#
# 暂不拟合：
#   Σ1/d + ΣHeat/d joint model
#
# 原因：
#   先审计两种weighted metric之间的相关性/VIF，
#   再决定joint模型是否具有可解释性。
#
# 不修改：
#   07A–07H04任何冻结/正式文件。
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


safe_int <- function(x) {
  suppressWarnings(
    as.integer(x)
  )
}


as_flag <- function(
  x,
  missing_value = FALSE
) {

  if (is.logical(x)) {

    z <- x

  } else {

    xx <- tolower(
      trimws(
        as.character(x)
      )
    )

    z <- rep(
      NA,
      length(xx)
    )

    z[
      xx %in%
        c(
          "true",
          "t",
          "1",
          "yes",
          "y"
        )
    ] <- TRUE

    z[
      xx %in%
        c(
          "false",
          "f",
          "0",
          "no",
          "n"
        )
    ] <- FALSE
  }

  z[
    is.na(z)
  ] <- missing_value

  z
}


fmt_int <- function(x) {
  format(
    x,
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


collapse_messages <- function(x) {

  if (
    is.null(x) ||
    length(x) == 0L
  ) {
    return("")
  }

  y <- unlist(
    x,
    recursive = TRUE,
    use.names = FALSE
  )

  y <- trimws(
    as.character(y)
  )

  y <- unique(
    y[
      !is.na(y) &
        nzchar(y)
    ]
  )

  if (length(y) == 0L) {
    return("")
  }

  paste(
    y,
    collapse = " | "
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


assert_dir <- function(
  dir,
  label
) {

  if (!dir.exists(dir)) {
    stopf(
      "找不到",
      label,
      "：\n",
      dir
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


strict_mean <- function(x) {

  x <- safe_num(x)

  if (
    length(x) == 0L ||
    any(
      !is.finite(x)
    )
  ) {
    return(
      NA_real_
    )
  }

  mean(x)
}


ci_excludes_zero <- function(
  low,
  high
) {

  is.finite(low) &
    is.finite(high) &
    (
      low > 0 |
        high < 0
    )
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


# ============================================================
# 2. Formal paths
# ============================================================

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

PP <- file.path(
  ROOT,
  "发电站数据处理"
)


F07_ROOT <- file.path(
  PP,
  "04_07F_正式统计分析数据构建"
)


F07_QDIR <- file.path(
  F07_ROOT,
  "quasi_daily_yearly"
)


FILTER_FILE <- file.path(
  F07_ROOT,
  "03_07F_filter_period_analysis_base_2001_2024.csv.gz"
)


H_ROOT <- file.path(
  PP,
  "06_07H_燃料特异CAMPD逐日暴露与模型"
)


H02_ROOT <- file.path(
  H_ROOT,
  "07H02_model"
)


H02_EFFECT_FILE <- file.path(
  H02_ROOT,
  "01_07H02_model_effects_all.csv"
)


H03_FREEZE_FILE <- file.path(
  H_ROOT,
  "07H03_模型结果审计与汇总",
  "07H03_FINAL_FREEZE.md"
)


H04_ROOT <- file.path(
  H_ROOT,
  "07H04_距离加权逐日暴露"
)


H04_DAILY_DIR <- file.path(
  H04_ROOT,
  "weighted_daily_yearly"
)


H04_HARD_FILE <- file.path(
  H04_ROOT,
  "qc",
  "05_build_hard_checks.csv"
)


OUT <- file.path(
  H_ROOT,
  "07H05_model"
)


FILTER_OUT_DIR <- file.path(
  OUT,
  "filter_period_weighted_by_fuel"
)


CHECKPOINT_DIR <- file.path(
  OUT,
  "checkpoints"
)


dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  FILTER_OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


dir.create(
  CHECKPOINT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


OUT_EFFECT <- file.path(
  OUT,
  "01_07H05_weighted_model_effects_all.csv"
)


OUT_STATUS <- file.path(
  OUT,
  "02_07H05_weighted_model_status_all.csv"
)


OUT_COR <- file.path(
  OUT,
  "03_07H05_weighted_metric_collinearity.csv"
)


OUT_SAMPLE <- file.path(
  OUT,
  "04_07H05_branch_sample_summary.csv"
)


OUT_COMPARE <- file.path(
  OUT,
  "05_07H05_vs_07H02_raw_sensitivity_comparison.csv"
)


OUT_COVERAGE <- file.path(
  OUT,
  "06_07H05_exposure_model_coverage.csv"
)


OUT_DISTANCE <- file.path(
  OUT,
  "07_07H05_distance_floor_audit.csv"
)


OUT_HARD <- file.path(
  OUT,
  "08_07H05_hard_checks.csv"
)


OUT_REPORT <- file.path(
  OUT,
  "09_07H05_距离加权敏感性模型报告.md"
)


OUT_SESSION <- file.path(
  OUT,
  "10_sessionInfo.txt"
)


# ============================================================
# 3. Constants
# ============================================================

YEARS <- 2001:2024


FUELS <- c(
  "coal",
  "gas",
  "oil",
  "other_fossil",
  "biomass",
  "solid_waste"
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


METRICS <- c(
  "proximity_invkm",
  "heat_proximity_mmbtu_per_km"
)


BRANCHES <- c(
  "quasi_daily",
  "filter_period"
)


WIND_EFFECT_MS <- 1


EXPECTED_Q_ROWS <- 769770L

EXPECTED_Q_PRIMARY <- 426657L

EXPECTED_Q_MONITORS <- 133L


EXPECTED_FP_ROWS <- 174200L

EXPECTED_FP_BASE <- 107407L

EXPECTED_FP_MONITORS <- 133L


EXPECTED_H04_MONITOR_DAYS <- 769600L

EXPECTED_MODEL_SPECS <- 480L

EXPECTED_COR_ROWS <- 240L

EXPECTED_COMPARISON_ROWS <- 480L


PIPELINE_VERSION <-
  "07H05_weighted_single_prespecified_20260919"


# ============================================================
# 4. Weighted exposure naming
# ============================================================

weighted_col <- function(
  fuel,
  direction,
  radius,
  metric
) {

  if (
    metric ==
      "proximity_invkm"
  ) {

    return(
      paste0(
        "campd_",
        fuel,
        "_",
        direction,
        "_proximity_invkm_",
        radius,
        "km_90deg"
      )
    )
  }


  if (
    metric ==
      "heat_proximity_mmbtu_per_km"
  ) {

    return(
      paste0(
        "campd_",
        fuel,
        "_",
        direction,
        "_heat_proximity_mmbtu_per_km_",
        radius,
        "km_90deg"
      )
    )
  }


  stopf(
    "未知weighted metric：",
    metric
  )
}


weighted_cols_for_fuel <- function(fuel) {

  out <- character()

  for (d in DIRECTIONS) {
    for (r in RADII) {
      for (m in METRICS) {
        out <- c(
          out,
          weighted_col(
            fuel,
            d,
            r,
            m
          )
        )
      }
    }
  }

  unique(out)
}


metric_increment <- function(metric) {

  if (
    metric ==
      "proximity_invkm"
  ) {
    return(0.1)
  }

  if (
    metric ==
      "heat_proximity_mmbtu_per_km"
  ) {
    return(10000)
  }

  NA_real_
}


metric_unit <- function(metric) {

  if (
    metric ==
      "proximity_invkm"
  ) {
    return("1/km")
  }

  if (
    metric ==
      "heat_proximity_mmbtu_per_km"
  ) {
    return("MMBtu/km")
  }

  NA_character_
}


metric_effect_label <- function(metric) {

  if (
    metric ==
      "proximity_invkm"
  ) {
    return(
      "Per +0.1 1/km at wind speed = 1 m/s"
    )
  }

  if (
    metric ==
      "heat_proximity_mmbtu_per_km"
  ) {
    return(
      "Per +10,000 MMBtu/km at wind speed = 1 m/s"
    )
  }

  NA_character_
}


# ============================================================
# 5. Prespecified 480 model specs
# ============================================================

SPECS <- CJ(
  branch = BRANCHES,
  fuel = FUELS,
  direction = DIRECTIONS,
  radius_km = RADII,
  metric = METRICS,
  unique = TRUE
)


SPECS[
  ,
  `:=`(
    model_type =
      "weighted_single",

    term =
      "exposure",

    raw_increment =
      vapply(
        metric,
        metric_increment,
        numeric(1)
      ),

    raw_unit =
      vapply(
        metric,
        metric_unit,
        character(1)
      ),

    effect_label =
      vapply(
        metric,
        metric_effect_label,
        character(1)
      )
  )
]


SPECS[
  ,
  spec_id :=
    paste(
      "07H05",
      branch,
      fuel,
      direction,
      paste0(
        radius_km,
        "km"
      ),
      metric,
      sep = "__"
    )
]


if (
  nrow(SPECS) !=
    EXPECTED_MODEL_SPECS
) {

  stopf(
    "07H05 model specs不是480；实际=",
    nrow(SPECS)
  )
}


if (
  uniqueN(
    SPECS$spec_id
  ) !=
    EXPECTED_MODEL_SPECS
) {

  stopf(
    "07H05 spec_id不是480个唯一值。"
  )
}


# ============================================================
# 6. Frozen upstream audit
# ============================================================

assert_file(
  FILTER_FILE,
  "07F filter-period正式表"
)


assert_file(
  H02_EFFECT_FILE,
  "07H02正式effect结果"
)


assert_file(
  H03_FREEZE_FILE,
  "07H03 FINAL FREEZE"
)


assert_file(
  H04_HARD_FILE,
  "07H04 hard checks"
)


assert_dir(
  F07_QDIR,
  "07F quasi-daily目录"
)


assert_dir(
  H04_DAILY_DIR,
  "07H04 weighted daily目录"
)


H03_FREEZE_TEXT <- readLines(
  H03_FREEZE_FILE,
  warn = FALSE,
  encoding = "UTF-8"
)


if (
  !any(
    grepl(
      "07H03 FINAL FREEZE: PASS",
      H03_FREEZE_TEXT,
      fixed = TRUE
    )
  )
) {

  stopf(
    "07H03尚未FINAL FREEZE PASS。"
  )
}


H04_HARD <- fread(
  H04_HARD_FILE,
  showProgress = FALSE
)


assert_cols(
  H04_HARD,
  c(
    "check",
    "status"
  ),
  "07H04 hard checks"
)


N_H04_HARD_FAIL <- H04_HARD[
  status !=
    "PASS",
  .N
]


if (
  N_H04_HARD_FAIL >
    0L
) {

  print(
    H04_HARD[
      status !=
        "PASS"
    ]
  )

  stopf(
    "07H04存在未PASS hard check，禁止进入07H05。"
  )
}


H04_FILES <- file.path(
  H04_DAILY_DIR,
  paste0(
    "07H04_fuel_daily_distance_weighted_exposure_",
    YEARS,
    ".csv.gz"
  )
)


if (
  !all(
    file.exists(
      H04_FILES
    )
  )
) {

  stopf(
    "07H04 yearly files不是24/24。"
  )
}


cat(
  "\n============================================================\n",
  "07H05 FROZEN INPUT AUDIT\n",
  "============================================================\n",
  "07H03 FINAL FREEZE: PASS\n",
  "07H04 hard checks: PASS\n",
  "07H04 yearly files: 24/24\n",
  "07H05 prespecified model specs: ",
  nrow(SPECS),
  "\n",
  sep = ""
)


# ============================================================
# 7. Read formal quasi-daily 07F model base
# ============================================================

Q_BASE_COLS <- c(

  "loc_num",
  "radnet_monitor_id",

  "collect_start",
  "collect_end",
  "date",

  "gross_beta_bq_m3",
  "gross_beta_mBq_m3",

  "c_prop",
  "surface_u_ppm_article",
  "sunspot_monthly",

  "pm_same_day",

  "inv_hpbl",
  "hpbl",

  "rhum",
  "soilm",
  "vel",

  "lat_diff_35",
  "year_model",
  "temp",

  "flag_conus",

  "n_same_day_pm_monitor_days",

  "flag_07f_primary_quasi_sample"
)


Q_PARTS <- vector(
  "list",
  length(YEARS)
)


for (
  ii in seq_along(YEARS)
) {

  yy <- YEARS[ii]

  f <- file.path(
    F07_QDIR,
    paste0(
      "07F_quasi_daily_analysis_base_",
      yy,
      ".csv.gz"
    )
  )


  assert_file(
    f,
    paste0(
      "07F quasi ",
      yy
    )
  )


  H <- names(
    fread(
      f,
      nrows = 0,
      showProgress = FALSE
    )
  )


  miss <- setdiff(
    Q_BASE_COLS,
    H
  )


  if (
    length(miss) >
      0L
  ) {

    stopf(
      yy,
      "年07F quasi缺少字段：",
      paste(
        miss,
        collapse = ", "
      )
    )
  }


  Q_PARTS[[ii]] <- fread(
    f,
    select = Q_BASE_COLS,
    showProgress = FALSE
  )
}


Q <- rbindlist(
  Q_PARTS,
  use.names = TRUE,
  fill = TRUE
)


rm(
  Q_PARTS
)

invisible(
  gc(
    verbose = FALSE
  )
)


Q[
  ,
  `:=`(
    loc_num =
      safe_int(
        loc_num
      ),

    date =
      as.IDate(
        as.character(
          date
        )
      ),

    flag_conus =
      as_flag(
        flag_conus,
        FALSE
      ),

    flag_07f_primary_quasi_sample =
      as_flag(
        flag_07f_primary_quasi_sample,
        FALSE
      )
  )
]


N_Q <- nrow(Q)


N_Q_PRIMARY <- Q[
  flag_07f_primary_quasi_sample %in%
    TRUE,
  .N
]


N_Q_MONITORS <- Q[
  flag_07f_primary_quasi_sample %in%
    TRUE,
  uniqueN(
    loc_num
  )
]


ELIGIBLE_LOCS <- sort(
  unique(
    Q[
      flag_07f_primary_quasi_sample %in%
        TRUE,
      loc_num
    ]
  )
)


if (
  N_Q !=
    EXPECTED_Q_ROWS
) {

  stopf(
    "07F quasi rows不是769770；实际=",
    N_Q
  )
}


if (
  N_Q_PRIMARY !=
    EXPECTED_Q_PRIMARY
) {

  stopf(
    "07F quasi primary rows不是426657；实际=",
    N_Q_PRIMARY
  )
}


if (
  N_Q_MONITORS !=
    EXPECTED_Q_MONITORS
) {

  stopf(
    "07F quasi primary monitors不是133；实际=",
    N_Q_MONITORS
  )
}


# ============================================================
# 8. Read formal filter-period table
# ============================================================

PRIMARY_FP <-
  "fpmean_campd_upwind_heat_observed_20km_90deg"


FP_BASE_COLS <- c(

  "loc_num",
  "radnet_monitor_id",

  "collect_start",
  "collect_end",

  "gross_beta_bq_m3",
  "gross_beta_mBq_m3",

  "c_prop",
  "surface_u_ppm_article",
  "sunspot_monthly",

  "pm_filter_period_official",

  "inv_hpbl",
  "hpbl",

  "rhum",
  "soilm",
  "vel",

  "lat_diff_35",
  "year_model",
  "temp",

  "flag_conus",

  "flag_07f_filter_period_primary_base",

  PRIMARY_FP
)


FP_HEADER <- names(
  fread(
    FILTER_FILE,
    nrows = 0,
    showProgress = FALSE
  )
)


FP_MISSING <- setdiff(
  FP_BASE_COLS,
  FP_HEADER
)


if (
  length(
    FP_MISSING
  ) >
    0L
) {

  stopf(
    "07F filter-period正式表缺少字段：",
    paste(
      FP_MISSING,
      collapse = ", "
    )
  )
}


FP <- fread(
  FILTER_FILE,
  select = FP_BASE_COLS,
  showProgress = FALSE
)


FP[
  ,
  `:=`(
    loc_num =
      safe_int(
        loc_num
      ),

    collect_start =
      as.IDate(
        as.character(
          collect_start
        )
      ),

    collect_end =
      as.IDate(
        as.character(
          collect_end
        )
      ),

    flag_conus =
      as_flag(
        flag_conus,
        FALSE
      ),

    flag_07f_filter_period_primary_base =
      as_flag(
        flag_07f_filter_period_primary_base,
        FALSE
      )
  )
]


FP[
  ,
  period_id :=
    .I
]


FP[
  ,
  n_calendar_days_expected :=
    as.integer(
      collect_end -
        collect_start
    )
]


if (
  FP[
    !is.finite(
      n_calendar_days_expected
    ) |
      n_calendar_days_expected <=
        0,
    .N
  ] >
    0L
) {

  stopf(
    "filter-period存在非正周期长度。"
  )
}


# Exact frozen 07G formal filter-period sample definition.
FP[
  ,
  flag_07g_filter_analysis_base :=

    flag_conus %in%
      TRUE &

    loc_num %in%
      ELIGIBLE_LOCS &

    is.finite(
      gross_beta_bq_m3
    ) &

    gross_beta_bq_m3 <
      0.04 &

    is.finite(
      pm_filter_period_official
    ) &

    pm_filter_period_official <
      35 &

    is.finite(
      vel
    ) &

    vel <
      6 &

    is.finite(
      hpbl
    ) &

    hpbl >
      150 &

    hpbl <
      3000 &

    is.finite(
      c_prop
    ) &

    is.finite(
      soilm
    ) &

    is.finite(
      get(
        PRIMARY_FP
      )
    )
]


N_FP <- nrow(FP)


N_FP_BASE <- FP[
  flag_07g_filter_analysis_base %in%
    TRUE,
  .N
]


N_FP_MONITORS <- FP[
  flag_07g_filter_analysis_base %in%
    TRUE,
  uniqueN(
    loc_num
  )
]


if (
  N_FP !=
    EXPECTED_FP_ROWS
) {

  stopf(
    "filter-period universe不是174200；实际=",
    N_FP
  )
}


if (
  N_FP_BASE !=
    EXPECTED_FP_BASE
) {

  stopf(
    "filter-period formal base不是107407；实际=",
    N_FP_BASE
  )
}


if (
  N_FP_MONITORS !=
    EXPECTED_FP_MONITORS
) {

  stopf(
    "filter-period formal monitors不是133；实际=",
    N_FP_MONITORS
  )
}


# ============================================================
# 9. Expand filter periods to [collect_start, collect_end)
# ============================================================

FP_KEYS <- FP[
  ,
  .(
    period_id,
    loc_num,
    collect_start,
    collect_end,
    n_calendar_days_expected
  )
]


FP_DAY_KEYS <- FP_KEYS[
  ,
  .(
    date =
      as.IDate(
        seq(
          as.Date(
            collect_start[1]
          ),
          as.Date(
            collect_end[1] -
              1L
          ),
          by = "day"
        )
      )
  ),
  by = .(
    period_id,
    loc_num
  )
]


EXPECTED_EXPANDED_DAYS <- sum(
  FP$n_calendar_days_expected
)


if (
  nrow(
    FP_DAY_KEYS
  ) !=
    EXPECTED_EXPANDED_DAYS
) {

  stopf(
    "filter-period日期展开行数异常；actual=",
    nrow(FP_DAY_KEYS),
    "; expected=",
    EXPECTED_EXPANDED_DAYS
  )
}


# ============================================================
# 10. Formal model formula
#
# Same structure as frozen 07G / 07H02.
# Internal exposure scaling is numerical only.
# Final effects are back-transformed to raw units.
# ============================================================

MODEL_FORMULA <- as.formula(
  paste0(
    "gross_beta_mBq_m3 ~ ",
    "exposure_scaled + ",
    "c_prop + ",
    "surface_u_ppm_article + ",
    "sunspot_monthly + ",
    "pm_model_centered + ",
    "inv_hpbl + ",
    "rhum + ",
    "soilm + ",
    "vel + ",
    "lat_diff_35 + ",
    "year_model + I(year_model^2) + ",
    "temp_centered + ",
    "I(temp_centered^2) + ",
    "I(temp_centered^3) + ",
    "sunspot_monthly:lat_diff_35 + ",
    "surface_u_ppm_article:soilm + ",
    "exposure_scaled:vel + ",
    "(1 | loc_num_factor)"
  )
)


MODEL_REQUIRED_BASE <- c(

  "gross_beta_mBq_m3",

  "c_prop",
  "surface_u_ppm_article",
  "sunspot_monthly",

  "inv_hpbl",
  "rhum",
  "soilm",
  "vel",

  "lat_diff_35",
  "year_model",
  "temp"
)


# ============================================================
# 11. Model-result helper
# ============================================================

blank_result <- function(
  meta,
  status,
  n_rows,
  n_monitors,
  n_unique,
  note,
  q25 = NA_real_,
  q75 = NA_real_,
  q95 = NA_real_,
  scale_value = NA_real_,
  sd_raw = NA_real_,
  zero_pct = NA_real_,
  singular = NA,
  convergence_message = "",
  warning_text = ""
) {

  cbind(
    copy(meta),

    data.table(

      status =
        status,

      n_rows =
        n_rows,

      n_monitors =
        n_monitors,

      n_unique_exposure =
        n_unique,

      q25_raw_exposure =
        q25,

      q75_raw_exposure =
        q75,

      q95_raw_exposure =
        q95,

      scale_value =
        scale_value,

      sd_raw_exposure =
        sd_raw,

      exposure_zero_pct =
        zero_pct,

      wind_speed_ms =
        WIND_EFFECT_MS,

      estimate_mBq_m3 =
        NA_real_,

      se_mBq_m3 =
        NA_real_,

      ci95_low_mBq_m3 =
        NA_real_,

      ci95_high_mBq_m3 =
        NA_real_,

      p_wald =
        NA_real_,

      estimate_per_1SD_mBq_m3 =
        NA_real_,

      se_per_1SD_mBq_m3 =
        NA_real_,

      ci95_low_per_1SD_mBq_m3 =
        NA_real_,

      ci95_high_per_1SD_mBq_m3 =
        NA_real_,

      singular =
        singular,

      convergence_message =
        convergence_message,

      warnings =
        warning_text,

      note =
        note
    )
  )
}


# ============================================================
# 12. Fit one single-exposure weighted model
# ============================================================

fit_weighted_model <- function(
  DT,
  exposure_col,
  pm_col,
  meta
) {

  required <- unique(
    c(
      MODEL_REQUIRED_BASE,
      exposure_col,
      pm_col,
      "loc_num"
    )
  )


  miss <- setdiff(
    required,
    names(DT)
  )


  if (
    length(miss) >
      0L
  ) {

    stopf(
      meta$spec_id,
      "缺少模型字段：",
      paste(
        miss,
        collapse = ", "
      )
    )
  }


  d <- copy(DT)


  d[
    ,
    `:=`(
      exposure_raw =
        safe_num(
          get(
            exposure_col
          )
        ),

      pm_model =
        safe_num(
          get(
            pm_col
          )
        ),

      temp_centered =
        safe_num(
          temp
        ) -
        273.15,

      loc_num_factor =
        factor(
          loc_num
        )
    )
  ]


  d[
    ,
    pm_model_centered :=
      pm_model -
      10
  ]


  model_vars <- c(

    "gross_beta_mBq_m3",

    "exposure_raw",

    "c_prop",
    "surface_u_ppm_article",
    "sunspot_monthly",

    "pm_model_centered",

    "inv_hpbl",
    "rhum",
    "soilm",
    "vel",

    "lat_diff_35",
    "year_model",
    "temp_centered",

    "loc_num"
  )


  keep <- complete.cases(
    d[
      ,
      ..model_vars
    ]
  )


  d <- d[
    keep
  ]


  n_rows <- nrow(d)

  n_monitors <- uniqueN(
    d$loc_num
  )


  x <- d$exposure_raw


  n_unique <- uniqueN(x)


  if (
    n_rows <
      1000L ||
    n_monitors <
      2L ||
    n_unique <
      2L
  ) {

    return(
      blank_result(
        meta =
          meta,

        status =
          "SKIP",

        n_rows =
          n_rows,

        n_monitors =
          n_monitors,

        n_unique =
          n_unique,

        note =
          "Insufficient rows/monitors/exposure variation"
      )
    )
  }


  q25 <- as.numeric(
    quantile(
      x,
      0.25,
      na.rm = TRUE,
      names = FALSE
    )
  )


  q75 <- as.numeric(
    quantile(
      x,
      0.75,
      na.rm = TRUE,
      names = FALSE
    )
  )


  q95 <- as.numeric(
    quantile(
      x,
      0.95,
      na.rm = TRUE,
      names = FALSE
    )
  )


  sd_raw <- sd(
    x,
    na.rm = TRUE
  )


  zero_pct <- 100 *
    mean(
      abs(x) <
        1e-15,
      na.rm = TRUE
    )


  # Same numerical-stability idea used in the prior model workflow:
  # P95 first; then IQR; then maximum absolute value.
  scale_value <- q95


  if (
    !is.finite(
      scale_value
    ) ||
    scale_value <=
      0
  ) {

    scale_value <-
      q75 -
      q25
  }


  if (
    !is.finite(
      scale_value
    ) ||
    scale_value <=
      0
  ) {

    scale_value <- max(
      abs(x),
      na.rm = TRUE
    )
  }


  if (
    !is.finite(
      scale_value
    ) ||
    scale_value <=
      0
  ) {

    return(
      blank_result(
        meta =
          meta,

        status =
          "SKIP",

        n_rows =
          n_rows,

        n_monitors =
          n_monitors,

        n_unique =
          n_unique,

        note =
          "Unable to define positive numerical exposure scale",

        q25 =
          q25,

        q75 =
          q75,

        q95 =
          q95,

        scale_value =
          scale_value,

        sd_raw =
          sd_raw,

        zero_pct =
          zero_pct
      )
    )
  }


  d[
    ,
    exposure_scaled :=
      exposure_raw /
      scale_value
  ]


  warning_text <- character()


  fit <- tryCatch(

    withCallingHandlers(

      lmer(
        MODEL_FORMULA,
        data = d,
        REML = TRUE,
        na.action = na.fail,
        control = lmerControl(
          autoscale = TRUE,
          optimizer = "bobyqa",
          optCtrl = list(
            maxfun = 200000
          ),
          check.conv.singular =
            "ignore",
          check.rankX =
            "message+drop.cols"
        )
      ),

      warning = function(w) {

        warning_text <<- c(
          warning_text,
          conditionMessage(w)
        )

        invokeRestart(
          "muffleWarning"
        )
      }
    ),

    error = function(e) {

      structure(
        list(
          message =
            conditionMessage(e)
        ),
        class =
          "fit_error"
      )
    }
  )


  if (
    inherits(
      fit,
      "fit_error"
    )
  ) {

    return(
      blank_result(
        meta =
          meta,

        status =
          "FAIL",

        n_rows =
          n_rows,

        n_monitors =
          n_monitors,

        n_unique =
          n_unique,

        note =
          "lmer fitting error",

        q25 =
          q25,

        q75 =
          q75,

        q95 =
          q95,

        scale_value =
          scale_value,

        sd_raw =
          sd_raw,

        zero_pct =
          zero_pct,

        convergence_message =
          fit$message,

        warning_text =
          collapse_messages(
            warning_text
          )
      )
    )
  }


  beta <- fixef(
    fit
  )


  V <- as.matrix(
    vcov(
      fit
    )
  )


  main_term <- "exposure_scaled"


  int_candidates <- c(
    "exposure_scaled:vel",
    "vel:exposure_scaled"
  )


  int_term <- int_candidates[
    int_candidates %in%
      names(beta)
  ]


  singular_flag <- isSingular(
    fit,
    tol = 1e-4
  )


  conv_text <- collapse_messages(
    fit@optinfo$conv$lme4$messages
  )


  warn_text <- collapse_messages(
    warning_text
  )


  if (
    !main_term %in%
      names(beta) ||
    length(int_term) !=
      1L
  ) {

    return(
      blank_result(
        meta =
          meta,

        status =
          "SKIP",

        n_rows =
          n_rows,

        n_monitors =
          n_monitors,

        n_unique =
          n_unique,

        note =
          "Exposure main/interaction term unavailable after model matrix rank handling",

        q25 =
          q25,

        q75 =
          q75,

        q95 =
          q95,

        scale_value =
          scale_value,

        sd_raw =
          sd_raw,

        zero_pct =
          zero_pct,

        singular =
          singular_flag,

        convergence_message =
          conv_text,

        warning_text =
          warn_text
      )
    )
  }


  int_term <- int_term[1]


  contrast_effect <- function(
    raw_increment
  ) {

    if (
      !is.finite(
        raw_increment
      )
    ) {

      return(
        c(
          estimate =
            NA_real_,
          se =
            NA_real_,
          low =
            NA_real_,
          high =
            NA_real_,
          p =
            NA_real_
        )
      )
    }


    delta_scaled <-
      raw_increment /
      scale_value


    cvec <- rep(
      0,
      length(beta)
    )


    names(cvec) <- names(
      beta
    )


    cvec[
      main_term
    ] <-
      delta_scaled


    cvec[
      int_term
    ] <-
      delta_scaled *
      WIND_EFFECT_MS


    est <- sum(
      cvec *
        beta
    )


    var_est <- as.numeric(
      t(cvec) %*%
        V %*%
        cvec
    )


    se <- if (
      is.finite(
        var_est
      ) &&
      var_est >=
        -1e-12
    ) {

      sqrt(
        max(
          var_est,
          0
        )
      )

    } else {

      NA_real_
    }


    low <- if (
      is.finite(se)
    ) {
      est -
        1.96 *
        se
    } else {
      NA_real_
    }


    high <- if (
      is.finite(se)
    ) {
      est +
        1.96 *
        se
    } else {
      NA_real_
    }


    p <- if (
      is.finite(est) &&
      is.finite(se) &&
      se >
        0
    ) {

      2 *
        pnorm(
          abs(
            est /
              se
          ),
          lower.tail = FALSE
        )

    } else {

      NA_real_
    }


    c(
      estimate =
        est,
      se =
        se,
      low =
        low,
      high =
        high,
      p =
        p
    )
  }


  E_MAIN <- contrast_effect(
    meta$raw_increment[1]
  )


  E_SD <- contrast_effect(
    sd_raw
  )


  cbind(

    copy(meta),

    data.table(

      status =
        "PASS",

      n_rows =
        n_rows,

      n_monitors =
        n_monitors,

      n_unique_exposure =
        n_unique,

      q25_raw_exposure =
        q25,

      q75_raw_exposure =
        q75,

      q95_raw_exposure =
        q95,

      scale_value =
        scale_value,

      sd_raw_exposure =
        sd_raw,

      exposure_zero_pct =
        zero_pct,

      wind_speed_ms =
        WIND_EFFECT_MS,

      estimate_mBq_m3 =
        unname(
          E_MAIN[
            "estimate"
          ]
        ),

      se_mBq_m3 =
        unname(
          E_MAIN[
            "se"
          ]
        ),

      ci95_low_mBq_m3 =
        unname(
          E_MAIN[
            "low"
          ]
        ),

      ci95_high_mBq_m3 =
        unname(
          E_MAIN[
            "high"
          ]
        ),

      p_wald =
        unname(
          E_MAIN[
            "p"
          ]
        ),

      estimate_per_1SD_mBq_m3 =
        unname(
          E_SD[
            "estimate"
          ]
        ),

      se_per_1SD_mBq_m3 =
        unname(
          E_SD[
            "se"
          ]
        ),

      ci95_low_per_1SD_mBq_m3 =
        unname(
          E_SD[
            "low"
          ]
        ),

      ci95_high_per_1SD_mBq_m3 =
        unname(
          E_SD[
            "high"
          ]
        ),

      singular =
        singular_flag,

      convergence_message =
        conv_text,

      warnings =
        warn_text,

      note =
        ""
    )
  )
}


# ============================================================
# 13. Weighted metric pair correlation / VIF
#
# This is an audit only.
# No weighted joint model is fitted here.
# ============================================================

weighted_pair_qc <- function(
  DT,
  proximity_col,
  heat_col,
  branch,
  fuel,
  direction,
  radius
) {

  x <- safe_num(
    DT[[
      proximity_col
    ]]
  )


  y <- safe_num(
    DT[[
      heat_col
    ]]
  )


  keep <-
    is.finite(x) &
    is.finite(y)


  xx <- x[
    keep
  ]

  yy <- y[
    keep
  ]


  n_both <- length(xx)


  n_unique_x <- uniqueN(xx)

  n_unique_y <- uniqueN(yy)


  pearson <- NA_real_

  spearman <- NA_real_

  vif <- NA_real_


  if (
    n_both >=
      3L &&
    n_unique_x >=
      2L &&
    n_unique_y >=
      2L
  ) {

    pearson <- suppressWarnings(
      cor(
        xx,
        yy,
        method = "pearson"
      )
    )


    spearman <- suppressWarnings(
      cor(
        xx,
        yy,
        method = "spearman"
      )
    )


    if (
      is.finite(
        pearson
      )
    ) {

      if (
        abs(
          pearson
        ) >=
          1 -
          1e-12
      ) {

        vif <- Inf

      } else {

        vif <-
          1 /
          (
            1 -
              pearson^2
          )
      }
    }
  }


  flag <- if (
    is.infinite(vif)
  ) {

    "PERFECT_COLLINEARITY_VIF_INF"

  } else if (
    is.finite(vif) &&
    vif >=
      10
  ) {

    "HIGH_VIF_GE10"

  } else if (
    is.finite(vif) &&
    vif >=
      5
  ) {

    "MODERATE_VIF_5_TO_10"

  } else if (
    is.finite(vif)
  ) {

    "NO_HIGH_VIF_FLAG"

  } else {

    "VIF_NOT_ESTIMABLE"
  }


  data.table(

    branch =
      branch,

    fuel =
      fuel,

    direction =
      direction,

    radius_km =
      radius,

    n_model_base =
      nrow(DT),

    n_both_finite =
      n_both,

    n_unique_proximity =
      n_unique_x,

    n_unique_heat_proximity =
      n_unique_y,

    proximity_zero_pct =
      if (
        length(x) >
          0L
      ) {
        100 *
          mean(
            is.finite(x) &
              abs(x) <
                1e-15
          )
      } else {
        NA_real_
      },

    heat_proximity_zero_pct =
      if (
        length(y) >
          0L
      ) {
        100 *
          mean(
            is.finite(y) &
              abs(y) <
                1e-15
          )
      } else {
        NA_real_
      },

    pearson_r =
      pearson,

    spearman_rho =
      spearman,

    two_predictor_vif =
      vif,

    collinearity_flag =
      flag
  )
}


# ============================================================
# 14. Read one fuel's 07H04 daily weighted exposures
# ============================================================

read_weighted_fuel <- function(fuel) {

  exposure_cols <- weighted_cols_for_fuel(
    fuel
  )


  required <- c(
    "loc_num",
    "date",
    "wind_valid",
    exposure_cols
  )


  parts <- vector(
    "list",
    length(YEARS)
  )


  for (
    ii in seq_along(YEARS)
  ) {

    yy <- YEARS[ii]


    f <- file.path(
      H04_DAILY_DIR,
      paste0(
        "07H04_fuel_daily_distance_weighted_exposure_",
        yy,
        ".csv.gz"
      )
    )


    H <- names(
      fread(
        f,
        nrows = 0,
        showProgress = FALSE
      )
    )


    miss <- setdiff(
      required,
      H
    )


    if (
      length(miss) >
        0L
    ) {

      stopf(
        yy,
        "年07H04缺少",
        fuel,
        "字段：",
        paste(
          miss,
          collapse = ", "
        )
      )
    }


    parts[[ii]] <- fread(
      f,
      select = required,
      showProgress = FALSE
    )
  }


  W <- rbindlist(
    parts,
    use.names = TRUE,
    fill = TRUE
  )


  W[
    ,
    `:=`(
      loc_num =
        safe_int(
          loc_num
        ),

      date =
        as.IDate(
          as.character(
            date
          )
        ),

      wind_valid =
        as_flag(
          wind_valid,
          FALSE
        ),

      h04_key_present =
        TRUE
    )
  ]


  dup <- W[
    ,
    .N,
    by = .(
      loc_num,
      date
    )
  ][
    N >
      1L,
    .N
  ]


  if (
    dup >
      0L
  ) {

    stopf(
      fuel,
      " 07H04 daily存在重复monitor-day。"
    )
  }


  if (
    nrow(W) !=
      EXPECTED_H04_MONITOR_DAYS
  ) {

    stopf(
      fuel,
      " 07H04 monitor-days不是769600；实际=",
      nrow(W)
    )
  }


  list(
    data =
      W,
    exposure_cols =
      exposure_cols
  )
}


# ============================================================
# 15. Aggregate one fuel to filter periods
#
# Exact rule:
# [collect_start, collect_end)
# arithmetic mean
# ANY daily exposure NA -> period exposure NA
# ============================================================

aggregate_filter_fuel <- function(
  W,
  fuel,
  exposure_cols
) {

  JOINED <- W[
    FP_DAY_KEYS,
    on = .(
      loc_num,
      date
    )
  ]


  N_MISSING_DAILY_KEYS <- JOINED[
    is.na(
      h04_key_present
    ),
    .N
  ]


  AGG <- JOINED[
    ,
    c(
      list(
        n_join_rows =
          .N,

        n_missing_h04_daily_keys =
          sum(
            is.na(
              h04_key_present
            )
          )
      ),

      lapply(
        .SD,
        strict_mean
      )
    ),

    by = period_id,

    .SDcols =
      exposure_cols
  ]


  FP_NAMES <- paste0(
    "fpmean_",
    exposure_cols
  )


  setnames(
    AGG,
    exposure_cols,
    FP_NAMES
  )


  if (
    nrow(AGG) !=
      N_FP
  ) {

    stopf(
      fuel,
      " filter-period weighted聚合不是174200行；实际=",
      nrow(AGG)
    )
  }


  CHECK <- merge(
    FP[
      ,
      .(
        period_id,
        n_calendar_days_expected
      )
    ],
    AGG[
      ,
      .(
        period_id,
        n_join_rows,
        n_missing_h04_daily_keys
      )
    ],
    by = "period_id",
    all.x = TRUE,
    sort = FALSE
  )


  N_DAY_COUNT_MISMATCH <- CHECK[
    is.na(
      n_join_rows
    ) |
      n_join_rows !=
        n_calendar_days_expected,
    .N
  ]


  if (
    N_DAY_COUNT_MISMATCH >
      0L
  ) {

    stopf(
      fuel,
      " filter-period日期展开/聚合日数不一致：",
      N_DAY_COUNT_MISMATCH
    )
  }


  FPF <- merge(
    FP,
    AGG[
      ,
      c(
        "period_id",
        FP_NAMES
      ),
      with = FALSE
    ],
    by = "period_id",
    all.x = TRUE,
    sort = FALSE
  )


  if (
    nrow(FPF) !=
      N_FP
  ) {

    stopf(
      fuel,
      " filter-period merge后行数改变。"
    )
  }


  save_cols <- c(
    "period_id",
    "loc_num",
    "collect_start",
    "collect_end",
    "flag_07g_filter_analysis_base",
    FP_NAMES
  )


  fwrite(
    FPF[
      ,
      ..save_cols
    ],
    file.path(
      FILTER_OUT_DIR,
      paste0(
        "07H05_filter_period_weighted_",
        fuel,
        ".csv.gz"
      )
    )
  )


  list(
    data =
      FPF,

    fp_cols =
      FP_NAMES,

    n_missing_daily_keys =
      N_MISSING_DAILY_KEYS,

    n_day_count_mismatch =
      N_DAY_COUNT_MISMATCH
  )
}


# ============================================================
# 16. Checkpoint logic
# ============================================================

checkpoint_file <- function(fuel) {

  file.path(
    CHECKPOINT_DIR,
    paste0(
      "07H05_",
      fuel,
      "_checkpoint.rds"
    )
  )
}


valid_checkpoint <- function(
  x,
  fuel
) {

  if (
    !is.list(x) ||
    is.null(
      x$pipeline_version
    ) ||
    x$pipeline_version !=
      PIPELINE_VERSION ||
    is.null(
      x$fuel
    ) ||
    x$fuel !=
      fuel ||
    is.null(
      x$results
    ) ||
    is.null(
      x$correlation
    )
  ) {

    return(FALSE)
  }


  if (
    nrow(
      x$results
    ) !=
      80L ||
    nrow(
      x$correlation
    ) !=
      40L
  ) {

    return(FALSE)
  }


  if (
    x$results[
      status ==
        "FAIL",
      .N
    ] >
      0L
  ) {

    return(FALSE)
  }


  TRUE
}


# ============================================================
# 17. Fit all six fuels
# ============================================================

RESULT_LIST <- list()

COR_LIST <- list()


TOTAL_Q_KEY_MISS <- 0L

TOTAL_FP_DAILY_KEY_MISS <- 0L

TOTAL_FP_DAY_COUNT_MISMATCH <- 0L


for (
  fuel in FUELS
) {

  cat(
    "\n============================================================\n",
    "07H05 FUEL: ",
    fuel,
    "\n",
    "============================================================\n",
    sep = ""
  )


  WF <- read_weighted_fuel(
    fuel
  )


  W <- WF$data

  EXP_COLS <- WF$exposure_cols


  # ----------------------------------------------------------
  # Quasi-daily join
  # ----------------------------------------------------------

  QF <- merge(
    Q,
    W[
      ,
      c(
        "loc_num",
        "date",
        "h04_key_present",
        EXP_COLS
      ),
      with = FALSE
    ],
    by = c(
      "loc_num",
      "date"
    ),
    all.x = TRUE,
    sort = FALSE
  )


  if (
    nrow(QF) !=
      N_Q
  ) {

    stopf(
      fuel,
      " quasi merge后行数不是769770；实际=",
      nrow(QF)
    )
  }


  N_Q_KEY_MISS <- QF[
    is.na(
      h04_key_present
    ),
    .N
  ]


  TOTAL_Q_KEY_MISS <-
    TOTAL_Q_KEY_MISS +
    N_Q_KEY_MISS


  if (
    N_Q_KEY_MISS >
      0L
  ) {

    stopf(
      fuel,
      " quasi存在未匹配07H04 monitor-day：",
      N_Q_KEY_MISS
    )
  }


  QF[
    ,
    h04_key_present :=
      NULL
  ]


  QB <- QF[
    flag_07f_primary_quasi_sample %in%
      TRUE
  ]


  if (
    nrow(QB) !=
      EXPECTED_Q_PRIMARY
  ) {

    stopf(
      fuel,
      " quasi primary merge后不是426657行。"
    )
  }


  # ----------------------------------------------------------
  # Strict filter-period aggregation
  # ----------------------------------------------------------

  FA <- aggregate_filter_fuel(
    W =
      W,

    fuel =
      fuel,

    exposure_cols =
      EXP_COLS
  )


  FPF <- FA$data


  TOTAL_FP_DAILY_KEY_MISS <-
    TOTAL_FP_DAILY_KEY_MISS +
    FA$n_missing_daily_keys


  TOTAL_FP_DAY_COUNT_MISMATCH <-
    TOTAL_FP_DAY_COUNT_MISMATCH +
    FA$n_day_count_mismatch


  if (
    FA$n_missing_daily_keys >
      0L
  ) {

    stopf(
      fuel,
      " filter-period日期存在未匹配07H04 daily key：",
      FA$n_missing_daily_keys
    )
  }


  FB <- FPF[
    flag_07g_filter_analysis_base %in%
      TRUE
  ]


  if (
    nrow(FB) !=
      EXPECTED_FP_BASE
  ) {

    stopf(
      fuel,
      " filter-period formal base merge后不是107407行。"
    )
  }


  # ----------------------------------------------------------
  # Reuse valid checkpoint only after all input joins passed.
  # ----------------------------------------------------------

  CP_FILE <- checkpoint_file(
    fuel
  )


  USE_CP <- FALSE


  if (
    file.exists(
      CP_FILE
    )
  ) {

    CP <- tryCatch(
      readRDS(
        CP_FILE
      ),
      error = function(e) {
        NULL
      }
    )


    if (
      !is.null(CP) &&
      valid_checkpoint(
        CP,
        fuel
      )
    ) {

      USE_CP <- TRUE

      cat(
        "有效checkpoint：复用 ",
        fuel,
        " 已完成模型。\n",
        sep = ""
      )
    }
  }


  if (
    USE_CP
  ) {

    RESULT_FUEL <-
      CP$results

    COR_FUEL <-
      CP$correlation

  } else {

    RESULT_FUEL_LIST <- list()

    COR_FUEL_LIST <- list()

    kk <- 0L

    cc <- 0L


    # --------------------------------------------------------
    # quasi-daily
    # --------------------------------------------------------

    for (
      direction in DIRECTIONS
    ) {

      for (
        radius in RADII
      ) {

        pcol <- weighted_col(
          fuel,
          direction,
          radius,
          "proximity_invkm"
        )


        hcol <- weighted_col(
          fuel,
          direction,
          radius,
          "heat_proximity_mmbtu_per_km"
        )


        cc <- cc + 1L

        COR_FUEL_LIST[[cc]] <-
          weighted_pair_qc(
            DT =
              QB,

            proximity_col =
              pcol,

            heat_col =
              hcol,

            branch =
              "quasi_daily",

            fuel =
              fuel,

            direction =
              direction,

            radius =
              radius
          )


        for (
          metric in METRICS
        ) {

          fuel_key <- fuel
          direction_key <- direction
          metric_key <- metric

          META <- SPECS[
            branch ==
              "quasi_daily" &
            fuel ==
              fuel_key &
            direction ==
              direction_key &
            radius_km ==
              radius &
            metric ==
              metric_key
          ]


          if (
            nrow(META) !=
              1L
          ) {

            stopf(
              "无法唯一定位quasi spec：",
              fuel,
              "/",
              direction,
              "/",
              radius,
              "/",
              metric
            )
          }


          EXP <- weighted_col(
            fuel,
            direction,
            radius,
            metric
          )


          cat(
            "[quasi] ",
            META$spec_id,
            "\n",
            sep = ""
          )


          kk <- kk + 1L

          RESULT_FUEL_LIST[[kk]] <-
            fit_weighted_model(
              DT =
                QB,

              exposure_col =
                EXP,

              pm_col =
                "pm_same_day",

              meta =
                META
            )
        }
      }
    }


    # --------------------------------------------------------
    # filter-period
    # --------------------------------------------------------

    for (
      direction in DIRECTIONS
    ) {

      for (
        radius in RADII
      ) {

        p_raw <- weighted_col(
          fuel,
          direction,
          radius,
          "proximity_invkm"
        )


        h_raw <- weighted_col(
          fuel,
          direction,
          radius,
          "heat_proximity_mmbtu_per_km"
        )


        pcol <- paste0(
          "fpmean_",
          p_raw
        )


        hcol <- paste0(
          "fpmean_",
          h_raw
        )


        cc <- cc + 1L

        COR_FUEL_LIST[[cc]] <-
          weighted_pair_qc(
            DT =
              FB,

            proximity_col =
              pcol,

            heat_col =
              hcol,

            branch =
              "filter_period",

            fuel =
              fuel,

            direction =
              direction,

            radius =
              radius
          )


        for (
          metric in METRICS
        ) {

          fuel_key <- fuel
          direction_key <- direction
          metric_key <- metric

          META <- SPECS[
            branch ==
              "filter_period" &
            fuel ==
              fuel_key &
            direction ==
              direction_key &
            radius_km ==
              radius &
            metric ==
              metric_key
          ]


          if (
            nrow(META) !=
              1L
          ) {

            stopf(
              "无法唯一定位filter spec：",
              fuel,
              "/",
              direction,
              "/",
              radius,
              "/",
              metric
            )
          }


          EXP <- paste0(
            "fpmean_",
            weighted_col(
              fuel,
              direction,
              radius,
              metric
            )
          )


          cat(
            "[filter] ",
            META$spec_id,
            "\n",
            sep = ""
          )


          kk <- kk + 1L

          RESULT_FUEL_LIST[[kk]] <-
            fit_weighted_model(
              DT =
                FB,

              exposure_col =
                EXP,

              pm_col =
                "pm_filter_period_official",

              meta =
                META
            )
        }
      }
    }


    RESULT_FUEL <- rbindlist(
      RESULT_FUEL_LIST,
      use.names = TRUE,
      fill = TRUE
    )


    COR_FUEL <- rbindlist(
      COR_FUEL_LIST,
      use.names = TRUE,
      fill = TRUE
    )


    if (
      nrow(
        RESULT_FUEL
      ) !=
        80L
    ) {

      stopf(
        fuel,
        " 模型结果不是80行；实际=",
        nrow(
          RESULT_FUEL
        )
      )
    }


    if (
      nrow(
        COR_FUEL
      ) !=
        40L
    ) {

      stopf(
        fuel,
        " weighted correlation不是40行；实际=",
        nrow(
          COR_FUEL
        )
      )
    }


    saveRDS(
      list(
        pipeline_version =
          PIPELINE_VERSION,

        fuel =
          fuel,

        results =
          RESULT_FUEL,

        correlation =
          COR_FUEL
      ),
      CP_FILE
    )
  }


  RESULT_LIST[[
    fuel
  ]] <-
    RESULT_FUEL


  COR_LIST[[
    fuel
  ]] <-
    COR_FUEL


  rm(
    WF,
    W,
    QF,
    QB,
    FA,
    FPF,
    FB,
    RESULT_FUEL,
    COR_FUEL
  )


  invisible(
    gc(
      verbose = FALSE
    )
  )
}


# ============================================================
# 18. Combine 480 model results
# ============================================================

RESULTS <- rbindlist(
  RESULT_LIST,
  use.names = TRUE,
  fill = TRUE
)


COR <- rbindlist(
  COR_LIST,
  use.names = TRUE,
  fill = TRUE
)


setorder(
  RESULTS,
  branch,
  fuel,
  direction,
  radius_km,
  metric
)


setorder(
  COR,
  branch,
  fuel,
  direction,
  radius_km
)


N_DUP_SPEC <- RESULTS[
  ,
  .N,
  by = spec_id
][
  N >
    1L,
  .N
]


N_DUP_COR <- COR[
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


fwrite(
  RESULTS,
  OUT_EFFECT
)


STATUS_COLS <- c(

  "spec_id",

  "branch",
  "fuel",
  "direction",
  "radius_km",
  "metric",

  "model_type",

  "status",

  "n_rows",
  "n_monitors",
  "n_unique_exposure",

  "singular",
  "convergence_message",
  "warnings",
  "note"
)


STATUS <- RESULTS[
  ,
  ..STATUS_COLS
]


fwrite(
  STATUS,
  OUT_STATUS
)


fwrite(
  COR,
  OUT_COR
)


# ============================================================
# 19. Branch sample summary
# ============================================================

SAMPLE_SUMMARY <- data.table(

  branch = c(
    "quasi_daily",
    "filter_period"
  ),

  universe_rows = c(
    N_Q,
    N_FP
  ),

  formal_base_rows = c(
    N_Q_PRIMARY,
    N_FP_BASE
  ),

  formal_base_monitors = c(
    N_Q_MONITORS,
    N_FP_MONITORS
  ),

  PM_variable = c(
    "pm_same_day",
    "pm_filter_period_official"
  ),

  formal_base_definition = c(
    "Frozen 07F flag_07f_primary_quasi_sample",
    "Frozen 07G filter-period formal analysis-base reconstruction"
  )
)


fwrite(
  SAMPLE_SUMMARY,
  OUT_SAMPLE
)


# ============================================================
# 20. Compare weighted sensitivity models with frozen 07H02
#
# Count-only        vs Σ1/d
# Heat-only         vs ΣHeat/d
#
# No ranking / no best radius selection.
# ============================================================

RAW <- fread(
  H02_EFFECT_FILE,
  showProgress = FALSE
)


RAW_REQUIRED <- c(
  "branch",
  "fuel",
  "direction",
  "radius_km",
  "model_type",
  "term",
  "status",
  "estimate_mBq_m3",
  "ci95_low_mBq_m3",
  "ci95_high_mBq_m3"
)


assert_cols(
  RAW,
  RAW_REQUIRED,
  "07H02 effects"
)


OPTIONAL_RAW_NUM <- c(
  "estimate_per_1SD_mBq_m3",
  "sd_raw_exposure",
  "p_wald"
)


for (
  cc in OPTIONAL_RAW_NUM
) {

  if (
    !cc %in%
      names(RAW)
  ) {

    RAW[
      ,
      (cc) :=
        NA_real_
    ]
  }
}


RAW_SINGLE <- RAW[
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
  )
]


RAW_SINGLE[
  ,
  comparison_family :=
    fifelse(
      model_type ==
        "count_only",
      "Count_vs_inverse_distance_proximity",
      "Heat_vs_inverse_distance_heat_proximity"
    )
]


if (
  nrow(
    RAW_SINGLE
  ) !=
    EXPECTED_COMPARISON_ROWS
) {

  stopf(
    "07H02单暴露raw effect rows不是480；实际=",
    nrow(
      RAW_SINGLE
    )
  )
}


RAW_COMPARE <- RAW_SINGLE[
  ,
  .(

    branch,
    fuel,
    direction,
    radius_km,
    comparison_family,

    raw_model_type =
      model_type,

    raw_status =
      status,

    raw_estimate_mBq_m3 =
      estimate_mBq_m3,

    raw_ci95_low_mBq_m3 =
      ci95_low_mBq_m3,

    raw_ci95_high_mBq_m3 =
      ci95_high_mBq_m3,

    raw_p_wald =
      p_wald,

    raw_sd_exposure =
      sd_raw_exposure,

    raw_estimate_per_1SD_mBq_m3 =
      estimate_per_1SD_mBq_m3
  )
]


WEIGHT_COMPARE <- RESULTS[
  ,
  .(

    branch,
    fuel,
    direction,
    radius_km,

    comparison_family =
      fifelse(
        metric ==
          "proximity_invkm",
        "Count_vs_inverse_distance_proximity",
        "Heat_vs_inverse_distance_heat_proximity"
      ),

    weighted_metric =
      metric,

    weighted_status =
      status,

    weighted_raw_increment =
      raw_increment,

    weighted_raw_unit =
      raw_unit,

    weighted_estimate_mBq_m3 =
      estimate_mBq_m3,

    weighted_ci95_low_mBq_m3 =
      ci95_low_mBq_m3,

    weighted_ci95_high_mBq_m3 =
      ci95_high_mBq_m3,

    weighted_p_wald =
      p_wald,

    weighted_sd_exposure =
      sd_raw_exposure,

    weighted_estimate_per_1SD_mBq_m3 =
      estimate_per_1SD_mBq_m3,

    weighted_ci95_low_per_1SD_mBq_m3 =
      ci95_low_per_1SD_mBq_m3,

    weighted_ci95_high_per_1SD_mBq_m3 =
      ci95_high_per_1SD_mBq_m3
  )
]


COMPARE <- merge(
  WEIGHT_COMPARE,
  RAW_COMPARE,

  by = c(
    "branch",
    "fuel",
    "direction",
    "radius_km",
    "comparison_family"
  ),

  all.x = TRUE,
  sort = FALSE
)


COMPARE[
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

    same_sign_fixed_effect =
      fifelse(
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

    same_sign_1SD_effect =
      fifelse(
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
        raw_ci95_low_mBq_m3,
        raw_ci95_high_mBq_m3
      ),

    weighted_ci_excludes_zero =
      ci_excludes_zero(
        weighted_ci95_low_mBq_m3,
        weighted_ci95_high_mBq_m3
      )
  )
]


setorder(
  COMPARE,
  branch,
  fuel,
  direction,
  radius_km,
  comparison_family
)


N_RAW_UNMATCHED <- COMPARE[
  is.na(
    raw_status
  ),
  .N
]


fwrite(
  COMPARE,
  OUT_COMPARE
)


# ============================================================
# 21. Model coverage
# ============================================================

COVERAGE <- RESULTS[
  ,
  .(

    n_specs =
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

    n_ci_excludes_zero =
      sum(
        status ==
          "PASS" &
        ci_excludes_zero(
          ci95_low_mBq_m3,
          ci95_high_mBq_m3
        ),
        na.rm = TRUE
      )
  ),

  by = .(
    branch,
    fuel,
    direction,
    metric
  )
]


setorder(
  COVERAGE,
  branch,
  fuel,
  direction,
  metric
)


fwrite(
  COVERAGE,
  OUT_COVERAGE
)


# ============================================================
# 22. 1-km distance-floor audit
#
# Informational only.
# Radius inclusion still used original 07A distance.
# ============================================================

FORMAL_07A_ROOT <- file.path(
  PP,
  "00_既往发电站空间暴露与模型流程_冻结",
  "02_07A_RadNet发电站空间链接",
  "07A_radnet_powerplant_same_year_250km"
)


GEO_FILE <- file.path(
  FORMAL_07A_ROOT,
  "03_monitor_year_all_plant_link_250km_07A_2001_2024.csv"
)


assert_file(
  GEO_FILE,
  "正式07A geometry"
)


GEO_HEADER <- names(
  fread(
    GEO_FILE,
    nrows = 0,
    showProgress = FALSE
  )
)


DIST_CANDIDATES <- c(
  "distance_km",
  "distance"
)


DIST_COL <- DIST_CANDIDATES[
  DIST_CANDIDATES %in%
    GEO_HEADER
]


if (
  length(DIST_COL) ==
    0L
) {

  stopf(
    "正式07A geometry找不到distance字段。"
  )
}


DIST_COL <- DIST_COL[1]


GD <- fread(
  GEO_FILE,
  select = DIST_COL,
  showProgress = FALSE
)


setnames(
  GD,
  DIST_COL,
  "distance_km"
)


GD[
  ,
  distance_km :=
    safe_num(
      distance_km
    )
]


GD200 <- GD[
  is.finite(
    distance_km
  ) &
  distance_km >=
    0 &
  distance_km <=
    200 +
    1e-8
]


N_LT1 <- GD200[
  distance_km <
    1,
  .N
]


DIST_AUDIT <- data.table(

  distance_definition =
    "Frozen 07A WGS84 geodesic distance",

  maximum_radius_km =
    200,

  denominator_floor_km =
    1,

  n_links_0_200km =
    nrow(
      GD200
    ),

  n_links_distance_lt1km =
    N_LT1,

  pct_links_distance_lt1km =
    100 *
    N_LT1 /
    nrow(
      GD200
    ),

  minimum_distance_km =
    if (
      nrow(
        GD200
      ) >
        0L
    ) {
      min(
        GD200$distance_km,
        na.rm = TRUE
      )
    } else {
      NA_real_
    },

  note =
    "1-km floor affects inverse-distance denominator only; cumulative-radius membership uses original distance_km."
)


fwrite(
  DIST_AUDIT,
  OUT_DISTANCE
)


# ============================================================
# 23. Technical status
# ============================================================

N_PASS <- RESULTS[
  status ==
    "PASS",
  .N
]


N_SKIP <- RESULTS[
  status ==
    "SKIP",
  .N
]


N_FAIL <- RESULTS[
  status ==
    "FAIL",
  .N
]


N_SINGULAR <- RESULTS[
  status ==
    "PASS" &
  singular %in%
    TRUE,
  .N
]


N_CONVERGENCE <- RESULTS[
  status ==
    "PASS" &
  !is.na(
    convergence_message
  ) &
  nzchar(
    trimws(
      convergence_message
    )
  ),
  .N
]


N_WARNING <- RESULTS[
  status ==
    "PASS" &
  !is.na(
    warnings
  ) &
  nzchar(
    trimws(
      warnings
    )
  ),
  .N
]


N_COR_INF <- COR[
  is.infinite(
    two_predictor_vif
  ),
  .N
]


N_COR_GE10 <- COR[
  is.finite(
    two_predictor_vif
  ) &
  two_predictor_vif >=
    10,
  .N
]


N_COR_5_10 <- COR[
  is.finite(
    two_predictor_vif
  ) &
  two_predictor_vif >=
    5 &
  two_predictor_vif <
    10,
  .N
]


# ============================================================
# 24. Hard checks
# ============================================================

N_FILTER_FILES <- length(
  list.files(
    FILTER_OUT_DIR,
    pattern =
      "^07H05_filter_period_weighted_(coal|gas|oil|other_fossil|biomass|solid_waste)\\.csv\\.gz$"
  )
)


HARD <- data.table(

  check = c(

    "07H04_hard_checks_all_pass",

    "07H04_yearly_files_24",

    "07F_quasi_rows_769770",

    "07F_quasi_primary_rows_426657",

    "07F_quasi_primary_monitors_133",

    "07F_filter_rows_174200",

    "07G_filter_formal_rows_107407",

    "07G_filter_formal_monitors_133",

    "model_specs_480",

    "model_results_480",

    "duplicate_spec_ids_zero",

    "weighted_correlation_rows_240",

    "duplicate_correlation_keys_zero",

    "raw_vs_weighted_comparison_rows_480",

    "raw_comparison_unmatched_zero",

    "filter_weighted_fuel_files_6",

    "quasi_H04_key_unmatched_zero",

    "filter_H04_daily_key_unmatched_zero",

    "filter_period_day_count_mismatch_zero",

    "fit_error_rows_zero"
  ),

  observed = c(

    as.character(
      N_H04_HARD_FAIL
    ),

    as.character(
      length(
        H04_FILES
      )
    ),

    as.character(
      N_Q
    ),

    as.character(
      N_Q_PRIMARY
    ),

    as.character(
      N_Q_MONITORS
    ),

    as.character(
      N_FP
    ),

    as.character(
      N_FP_BASE
    ),

    as.character(
      N_FP_MONITORS
    ),

    as.character(
      nrow(
        SPECS
      )
    ),

    as.character(
      nrow(
        RESULTS
      )
    ),

    as.character(
      N_DUP_SPEC
    ),

    as.character(
      nrow(
        COR
      )
    ),

    as.character(
      N_DUP_COR
    ),

    as.character(
      nrow(
        COMPARE
      )
    ),

    as.character(
      N_RAW_UNMATCHED
    ),

    as.character(
      N_FILTER_FILES
    ),

    as.character(
      TOTAL_Q_KEY_MISS
    ),

    as.character(
      TOTAL_FP_DAILY_KEY_MISS
    ),

    as.character(
      TOTAL_FP_DAY_COUNT_MISMATCH
    ),

    as.character(
      N_FAIL
    )
  ),

  expected = c(

    "0",

    "24",

    "769770",

    "426657",

    "133",

    "174200",

    "107407",

    "133",

    "480",

    "480",

    "0",

    "240",

    "0",

    "480",

    "0",

    "6",

    "0",

    "0",

    "0",

    "0"
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
# 25. Console audit
# ============================================================

cat(
  "\n============================================================\n",
  "A. 07H05 MODEL STATUS\n",
  "============================================================\n",
  sep = ""
)


print(
  RESULTS[
    ,
    .N,
    by = status
  ][
    order(
      status
    )
  ]
)


cat(
  "\n============================================================\n",
  "B. BRANCH × FUEL COVERAGE\n",
  "============================================================\n",
  sep = ""
)


print(
  RESULTS[
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
  "C. WEIGHTED METRIC COLLINEARITY\n",
  "============================================================\n",
  "VIF=Inf = ",
  N_COR_INF,
  "\n",
  "Finite VIF>=10 = ",
  N_COR_GE10,
  "\n",
  "5<=VIF<10 = ",
  N_COR_5_10,
  "\n",
  sep = ""
)


cat(
  "\n============================================================\n",
  "D. MODEL TECHNICAL STATUS\n",
  "============================================================\n",
  "Singular PASS specs = ",
  N_SINGULAR,
  "\n",
  "Convergence-message PASS specs = ",
  N_CONVERGENCE,
  "\n",
  "Warning PASS specs = ",
  N_WARNING,
  "\n",
  sep = ""
)


cat(
  "\n============================================================\n",
  "E. DISTANCE FLOOR AUDIT\n",
  "============================================================\n",
  sep = ""
)


print(
  DIST_AUDIT,
  nrows = Inf
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
# 26. Report
# ============================================================

REPORT <- c(

  "# 07H05 距离加权敏感性统计模型",

  "",

  paste0(
    "- Generated: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),

  "",

  "## Analysis role",

  "",

  "- Sensitivity analysis only; does not replace frozen 07H02 raw Count/Heat analysis.",

  "- Weighted metric 1: Σ1/max(distance,1 km), unit 1/km.",

  "- Weighted metric 2: ΣHeat/max(distance,1 km), unit MMBtu/km.",

  "- Original 07A distance determines cumulative-radius membership.",

  "- The 1-km floor affects only the inverse-distance denominator.",

  "",

  "## Prespecified structure",

  "",

  "- Fuels: coal, gas, oil, other_fossil, biomass, solid_waste.",

  "- Directions: upwind and downwind.",

  "- Cumulative radii: 20,25,30,35,40,45,50,100,150,200 km.",

  "- Branches: quasi_daily and filter_period.",

  "- 2 weighted single-exposure metrics.",

  "- Total prespecified model specs: 480.",

  "- No significance-based radius selection.",

  "- No weighted joint model is fitted in 07H05.",

  "",

  "## Filter-period aggregation",

  "",

  "- Window: [collect_start, collect_end).",

  "- Aggregation: arithmetic mean of daily weighted exposure.",

  "- Any daily NA within the period causes period weighted exposure to remain NA.",

  "",

  "## Model",

  "",

  "- Same formal LME covariate architecture as frozen 07G/07H02.",

  "- Random intercept: monitor (loc_num).",

  "- REML = TRUE.",

  "- Optimizer = bobyqa; maxfun = 200000.",

  "- Exposure × wind-speed interaction retained.",

  "- Effects evaluated at wind speed = 1 m/s.",

  "- Internal scaling is numerical only; fixed-increment and 1-SD effects are returned on the raw exposure scale.",

  "",

  "## Effect increments",

  "",

  "- Σ1/d: +0.1 1/km.",

  "- ΣHeat/d: +10,000 MMBtu/km.",

  "- 1-SD effects are also reported for direct robustness comparison with frozen raw Count/Heat models.",

  "",

  "## Results status",

  "",

  paste0(
    "- PASS model specs: ",
    N_PASS,
    "."
  ),

  paste0(
    "- SKIP model specs: ",
    N_SKIP,
    "."
  ),

  paste0(
    "- FAIL model specs: ",
    N_FAIL,
    "."
  ),

  paste0(
    "- Singular PASS specs: ",
    N_SINGULAR,
    "."
  ),

  paste0(
    "- PASS specs with convergence message: ",
    N_CONVERGENCE,
    "."
  ),

  paste0(
    "- PASS specs with warning: ",
    N_WARNING,
    "."
  ),

  "",

  "## Weighted-metric collinearity audit",

  "",

  paste0(
    "- VIF=Inf: ",
    N_COR_INF,
    "."
  ),

  paste0(
    "- Finite VIF>=10: ",
    N_COR_GE10,
    "."
  ),

  paste0(
    "- 5<=VIF<10: ",
    N_COR_5_10,
    "."
  ),

  "- These diagnostics are used to decide whether a later weighted joint model is estimable; 07H05 itself does not fit that joint model.",

  "",

  "## Frozen raw-vs-weighted comparison",

  "",

  "- Count-only is paired with Σ1/d.",

  "- Heat-only is paired with ΣHeat/d.",

  "- Comparison preserves all prespecified radii and both directions/branches.",

  "- 1-SD estimates are the preferred magnitude comparison because raw units differ.",

  "",

  if (
    N_HARD_FAIL ==
      0L
  ) {

    "**07H05 DISTANCE-WEIGHTED SENSITIVITY MODEL: PASS**"

  } else {

    "**07H05 DISTANCE-WEIGHTED SENSITIVITY MODEL: FAIL**"
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
# 27. Final
# ============================================================

cat(
  "\n============================================================\n"
)


if (
  N_HARD_FAIL ==
    0L
) {

  cat(
    "07H05 DISTANCE-WEIGHTED SENSITIVITY MODEL: PASS\n"
  )

} else {

  cat(
    "07H05 DISTANCE-WEIGHTED SENSITIVITY MODEL: FAIL\n"
  )
}


cat(
  "============================================================\n",
  "Model specs = ",
  nrow(
    RESULTS
  ),
  " / 480\n",
  "PASS = ",
  N_PASS,
  "\n",
  "SKIP = ",
  N_SKIP,
  "\n",
  "FAIL = ",
  N_FAIL,
  "\n",
  "Weighted correlation rows = ",
  nrow(
    COR
  ),
  " / 240\n",
  "Raw-vs-weighted comparison rows = ",
  nrow(
    COMPARE
  ),
  " / 480\n",
  "VIF=Inf = ",
  N_COR_INF,
  "\n",
  "Finite VIF>=10 = ",
  N_COR_GE10,
  "\n",
  "5<=VIF<10 = ",
  N_COR_5_10,
  "\n",
  "Singular PASS specs = ",
  N_SINGULAR,
  "\n",
  "Convergence-message PASS specs = ",
  N_CONVERGENCE,
  "\n",
  "Warning PASS specs = ",
  N_WARNING,
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
    "07H05存在hard failure。checkpoint已保留，可修正后继续。"
  )
}


# 正式全部PASS以后删除临时checkpoint；
# 正式CSV结果与filter-period weighted exposure均保留。
if (
  dir.exists(
    CHECKPOINT_DIR
  )
) {

  unlink(
    CHECKPOINT_DIR,
    recursive = TRUE,
    force = TRUE
  )
}


cat(
  "07H05 checkpoints：正式PASS后已自动清理。\n"
)

