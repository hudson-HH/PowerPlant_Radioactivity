#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999,
  width = 220
)

required_packages <- c(
  "data.table",
  "fixest",
  "lme4"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(lme4)
})

N_CORES <- suppressWarnings(
  parallel::detectCores(
    logical = TRUE
  )
)

if (
  is.finite(N_CORES) &&
  N_CORES >= 2
) {
  fixest::setFixest_nthreads(
    max(
      1L,
      floor(N_CORES / 2)
    )
  )
}

options(
  fixest_notes = FALSE
)

CHECKPOINT_VERSION <- "07H10_20260919_v1"

FORCE_REBUILD <- FALSE

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

PP <- file.path(
  ROOT,
  "发电站数据处理"
)

HROOT <- file.path(
  PP,
  "06_07H_燃料特异CAMPD逐日暴露与模型"
)

F07 <- file.path(
  PP,
  "04_07F_正式统计分析数据构建"
)

QDIR <- file.path(
  F07,
  "quasi_daily_yearly"
)

FP_FILE <- file.path(
  F07,
  "03_07F_filter_period_analysis_base_2001_2024.csv.gz"
)

RAW_DAILY_DIR <- file.path(
  HROOT,
  "fuel_daily_yearly"
)

H04_DIR <- file.path(
  HROOT,
  "07H04_距离加权逐日暴露"
)

H07_DIR <- file.path(
  HROOT,
  "07H07_论文结果表与图形准备"
)

H07_FREEZE <- file.path(
  H07_DIR,
  "07H07_FINAL_FREEZE.md"
)

H07_MAIN <- file.path(
  H07_DIR,
  "01_论文结果表",
  "01_论文主结果_upwind单暴露_全部预设半径.csv"
)

H09B_DIR <- file.path(
  HROOT,
  "07H09B_监测站年份固定效应严格验证"
)

H09B_HARD <- file.path(
  H09B_DIR,
  "08_07H09B_hard_checks.csv"
)

OUT <- file.path(
  HROOT,
  "07H10_燃料特异站点年份内综合正式分析"
)

CHECK_DIR <- file.path(
  OUT,
  "00_checkpoint"
)

DIR_MAIN <- file.path(
  OUT,
  "01_主分析_监测站年份固定效应"
)

DIR_STRICT <- file.path(
  OUT,
  "02_严格时间固定效应"
)

DIR_WIND <- file.path(
  OUT,
  "03_连续风速分析"
)

DIR_DIRECTION <- file.path(
  OUT,
  "04_上风向下风向直接比较"
)

DIR_WEIGHTED <- file.path(
  OUT,
  "05_距离加权敏感性"
)

DIR_JOINT <- file.path(
  OUT,
  "06_数量热输入联合模型"
)

DIR_PM <- file.path(
  OUT,
  "07_细颗粒物敏感性"
)

DIR_LAG <- file.path(
  OUT,
  "08_准逐日滞后分析"
)

DIR_DECOMP <- file.path(
  OUT,
  "09_站点年份内外分解"
)

DIR_SUMMARY <- file.path(
  OUT,
  "10_综合汇总"
)

for (
  d in c(
    OUT,
    CHECK_DIR,
    DIR_MAIN,
    DIR_STRICT,
    DIR_WIND,
    DIR_DIRECTION,
    DIR_WEIGHTED,
    DIR_JOINT,
    DIR_PM,
    DIR_LAG,
    DIR_DECOMP,
    DIR_SUMMARY
  )
) {

  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

YEARS <- 2001:2024

RADII <- c(
  20,
  25,
  30,
  35,
  40,
  45,
  50,
  100,
  150,
  200
)

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

FUEL_CN <- c(
  coal = "煤",
  gas = "天然气",
  oil = "燃油",
  other_fossil = "其他化石燃料",
  biomass = "生物质",
  solid_waste = "固体废物"
)

BRANCHES <- c(
  "filter_period",
  "quasi_daily"
)

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

collapse_msg <- function(x) {

  x <- unique(
    as.character(x)
  )

  x <- x[
    !is.na(x) &
    nzchar(x)
  ]

  paste(
    x,
    collapse = " | "
  )
}

strict_mean <- function(x) {

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

first_finite <- function(x) {

  z <- x[
    is.finite(x)
  ]

  if (length(z) == 0L) {
    return(
      NA_real_
    )
  }

  z[1]
}

direction_text <- function(x) {

  fifelse(
    !is.finite(x),
    NA_character_,
    fifelse(
      x > 0,
      "正向",
      fifelse(
        x < 0,
        "负向",
        "零"
      )
    )
  )
}

ci_text <- function(
  low,
  high
) {

  fifelse(
    !is.finite(low) |
    !is.finite(high),
    NA_character_,
    fifelse(
      low > 0,
      "95%置信区间完全高于0",
      fifelse(
        high < 0,
        "95%置信区间完全低于0",
        "95%置信区间包含0"
      )
    )
  )
}

raw_col <- function(
  fuel,
  direction,
  metric,
  radius
) {

  if (metric == "count") {

    paste0(
      "campd_",
      fuel,
      "_",
      direction,
      "_n_plants_",
      radius,
      "km_90deg"
    )

  } else {

    paste0(
      "campd_",
      fuel,
      "_",
      direction,
      "_heat_observed_",
      radius,
      "km_90deg"
    )
  }
}

metric_cn <- function(metric) {

  if (metric == "count") {
    return(
      "电厂数量"
    )
  }

  if (metric == "heat") {
    return(
      "实际热输入"
    )
  }

  if (metric == "proximity") {
    return(
      "Σ1/d"
    )
  }

  if (metric == "heat_proximity") {
    return(
      "ΣHeat/d"
    )
  }

  metric
}

metric_scale <- function(metric) {

  if (metric == "count") {
    return(1)
  }

  if (metric == "heat") {
    return(10000)
  }

  if (metric == "proximity") {
    return(0.1)
  }

  if (metric == "heat_proximity") {
    return(10000)
  }

  1
}

metric_unit <- function(metric) {

  if (metric == "count") {
    return(
      "每增加1座活动电厂"
    )
  }

  if (metric == "heat") {
    return(
      "每增加10000 MMBtu实际热输入"
    )
  }

  if (metric == "proximity") {
    return(
      "每增加0.1 1/km"
    )
  }

  if (metric == "heat_proximity") {
    return(
      "每增加10000 MMBtu/km"
    )
  }

  ""
}

find_interaction <- function(
  coef_names,
  a,
  b
) {

  z <- intersect(
    c(
      paste0(a, ":", b),
      paste0(b, ":", a)
    ),
    coef_names
  )

  if (length(z) != 1L) {
    return(
      NA_character_
    )
  }

  z[1]
}

linear_combination <- function(
  beta,
  V,
  weights
) {

  if (
    length(weights) == 0L ||
    any(
      !names(weights) %in%
      names(beta)
    )
  ) {

    return(
      c(
        estimate = NA_real_,
        se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p = NA_real_
      )
    )
  }

  cvec <- setNames(
    rep(
      0,
      length(beta)
    ),
    names(beta)
  )

  cvec[
    names(weights)
  ] <- weights

  est <- sum(
    cvec *
    beta
  )

  vv <- as.numeric(
    t(cvec) %*%
      V %*%
      cvec
  )

  se <- if (
    is.finite(vv)
  ) {

    sqrt(
      max(
        vv,
        0
      )
    )

  } else {

    NA_real_
  }

  low <- if (
    is.finite(est) &&
    is.finite(se)
  ) {
    est -
      1.96 *
      se
  } else {
    NA_real_
  }

  high <- if (
    is.finite(est) &&
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
    se > 0
  ) {

    2 *
      pnorm(
        abs(
          est / se
        ),
        lower.tail = FALSE
      )

  } else {

    NA_real_
  }

  c(
    estimate = est,
    se = se,
    ci_low = low,
    ci_high = high,
    p = p
  )
}

effect_at_wind <- function(
  beta,
  V,
  main_term,
  interaction_term,
  wind
) {

  if (
    !main_term %in%
      names(beta) ||
    is.na(interaction_term) ||
    !interaction_term %in%
      names(beta)
  ) {

    return(
      c(
        estimate = NA_real_,
        se = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p = NA_real_
      )
    )
  }

  linear_combination(
    beta,
    V,
    setNames(
      c(
        1,
        wind
      ),
      c(
        main_term,
        interaction_term
      )
    )
  )
}

reference_winds <- function(x) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0L) {

    return(
      data.table(
        wind_reference =
          character(),
        wind_ms =
          numeric()
      )
    )
  }

  data.table(
    wind_reference = c(
      "参考论文1米每秒",
      "样本第25百分位",
      "样本中位数",
      "样本均值",
      "样本第75百分位"
    ),
    wind_ms = c(
      1,
      unname(
        quantile(
          x,
          0.25
        )
      ),
      median(x),
      mean(x),
      unname(
        quantile(
          x,
          0.75
        )
      )
    )
  )
}

make_lmer_control <- function() {

  z <- list(
    optimizer =
      "bobyqa",

    optCtrl =
      list(
        maxfun =
          200000
      ),

    check.conv.singular =
      "ignore",

    check.rankX =
      "message+drop.cols"
  )

  formal_names <- names(
    formals(
      lme4::lmerControl
    )
  )

  if (
    "autoscale" %in%
      formal_names
  ) {
    z$autoscale <-
      TRUE
  }

  do.call(
    lme4::lmerControl,
    z
  )
}

checkpoint_file <- function(
  radius,
  family
) {

  file.path(
    CHECK_DIR,
    paste0(
      "radius_",
      radius,
      "km__",
      family,
      ".rds"
    )
  )
}

checkpoint_run <- function(
  radius,
  family,
  fun
) {

  f <- checkpoint_file(
    radius,
    family
  )

  if (
    file.exists(f) &&
    !FORCE_REBUILD
  ) {

    obj <- tryCatch(
      readRDS(f),
      error =
        function(e)
          NULL
    )

    if (
      !is.null(obj) &&
      is.list(obj) &&
      identical(
        obj$checkpoint_version,
        CHECKPOINT_VERSION
      )
    ) {

      cat(
        "[checkpoint] ",
        radius,
        " km | ",
        family,
        "：已完成，直接读取。\n",
        sep = ""
      )

      return(
        obj$result
      )
    }
  }

  cat(
    "[运行] ",
    radius,
    " km | ",
    family,
    "\n",
    sep = ""
  )

  ans <- fun()

  saveRDS(
    list(
      checkpoint_version =
        CHECKPOINT_VERSION,
      radius_km =
        radius,
      family =
        family,
      result =
        ans
    ),
    f
  )

  ans
}

if (!file.exists(FP_FILE)) {
  stopf(
    "找不到07F filter-period正式表：",
    FP_FILE
  )
}

if (!file.exists(H07_FREEZE)) {
  stopf(
    "找不到07H07冻结标记。"
  )
}

if (!file.exists(H07_MAIN)) {
  stopf(
    "找不到07H07主结果表。"
  )
}

if (!file.exists(H09B_HARD)) {
  stopf(
    "找不到07H09B硬检查结果。"
  )
}

H09B_AUDIT <- fread(
  file =
    H09B_HARD,
  showProgress =
    FALSE
)

if (
  !"status" %in%
    names(H09B_AUDIT) ||
  H09B_AUDIT[
    status !=
      "PASS",
    .N
  ] >
    0L
) {

  stopf(
    "07H09B没有通过全部硬检查，禁止进入07H10。"
  )
}

RAW_FILES <- file.path(
  RAW_DAILY_DIR,
  paste0(
    "07H_fuel_daily_exposure_",
    YEARS,
    ".csv.gz"
  )
)

if (
  !all(
    file.exists(
      RAW_FILES
    )
  )
) {

  print(
    RAW_FILES[
      !file.exists(
        RAW_FILES
      )
    ]
  )

  stopf(
    "07H01逐日燃料暴露文件不是24/24。"
  )
}

if (!dir.exists(H04_DIR)) {

  stopf(
    "找不到07H04距离加权目录：",
    H04_DIR
  )
}

Q_BASE_COLS <- c(
  "loc_num",
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
  "flag_07f_primary_quasi_sample"
)

cat(
  "\n============================================================\n",
  "读取07F准逐日基础数据\n",
  "============================================================\n",
  sep = ""
)

Q_ALL_LIST <- vector(
  "list",
  length(YEARS)
)

for (
  ii in
  seq_along(YEARS)
) {

  yy <- YEARS[ii]

  f <- file.path(
    QDIR,
    paste0(
      "07F_quasi_daily_analysis_base_",
      yy,
      ".csv.gz"
    )
  )

  if (!file.exists(f)) {
    stopf(
      "缺少07F准逐日文件：",
      yy
    )
  }

  H <- names(
    fread(
      file = f,
      nrows = 0,
      showProgress = FALSE
    )
  )

  miss <- setdiff(
    Q_BASE_COLS,
    H
  )

  if (length(miss) > 0L) {

    stopf(
      yy,
      "年准逐日文件缺少字段：",
      paste(
        miss,
        collapse = ", "
      )
    )
  }

  Q_ALL_LIST[[ii]] <-
    fread(
      file = f,
      select =
        Q_BASE_COLS,
      showProgress =
        FALSE
    )
}

Q_ALL <- rbindlist(
  Q_ALL_LIST,
  use.names = TRUE
)

Q_ALL[
  ,
  `:=`(
    loc_num =
      as.character(
        loc_num
      ),
    date =
      as.IDate(
        date
      ),
    collect_start =
      as.IDate(
        collect_start
      ),
    collect_end =
      as.IDate(
        collect_end
      ),
    vel =
      safe_num(
        vel
      )
  )
]

WIND_CONFLICT <- Q_ALL[
  is.finite(vel),
  .(
    n_wind_values =
      uniqueN(
        round(
          vel,
          10
        )
      )
  ),
  by = .(
    loc_num,
    date
  )
][
  n_wind_values >
    1L,
  .N
]

if (
  WIND_CONFLICT >
    0L
) {

  stopf(
    "同一monitor-day存在不一致风速，不能构建滞后分析。"
  )
}

WIND_LOOKUP <- Q_ALL[
  ,
  .(
    exposure_day_wind =
      first_finite(
        vel
      )
  ),
  by = .(
    loc_num,
    date
  )
]

Q_BASE <- Q_ALL[
  flag_07f_primary_quasi_sample %in%
    TRUE
]

N_Q <- nrow(
  Q_BASE
)

ELIGIBLE_LOCS <- sort(
  unique(
    Q_BASE$loc_num
  )
)

rm(
  Q_ALL,
  Q_ALL_LIST
)

invisible(
  gc(
    verbose = FALSE
  )
)

cat(
  "正式准逐日样本 = ",
  format(
    N_Q,
    big.mark = ","
  ),
  "\n正式准逐日监测站 = ",
  length(
    ELIGIBLE_LOCS
  ),
  "\n",
  sep = ""
)

FP_REQUIRED <- c(
  "loc_num",
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
  "fpmean_campd_upwind_heat_observed_20km_90deg"
)

HFP <- names(
  fread(
    file =
      FP_FILE,
    nrows = 0,
    showProgress = FALSE
  )
)

FP_MISS <- setdiff(
  FP_REQUIRED,
  HFP
)

if (
  length(
    FP_MISS
  ) >
  0L
) {

  stopf(
    "filter-period正式表缺少：",
    paste(
      FP_MISS,
      collapse = ", "
    )
  )
}

FP_BASE <- fread(
  file =
    FP_FILE,
  select =
    FP_REQUIRED,
  showProgress =
    FALSE
)

FP_BASE[
  ,
  `:=`(
    loc_num =
      as.character(
        loc_num
      ),
    collect_start =
      as.IDate(
        collect_start
      ),
    collect_end =
      as.IDate(
        collect_end
      )
  )
]

N_FP_ALL <- nrow(
  FP_BASE
)

FP_BASE[
  ,
  flag_07h10_filter_base :=
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
      fpmean_campd_upwind_heat_observed_20km_90deg
    )
]

FP_BASE <- FP_BASE[
  flag_07h10_filter_base %in%
    TRUE
]

FP_BASE[
  ,
  period_id :=
    .I
]

FP_BASE[
  ,
  period_days :=
    as.integer(
      collect_end -
      collect_start
    )
]

if (
  FP_BASE[
    !is.finite(
      period_days
    ) |
    period_days <=
      0,
    .N
  ] >
  0L
) {

  stopf(
    "发现非法filter period。"
  )
}

FP_DAY_MAP <- FP_BASE[
  ,
  .(
    date =
      seq(
        collect_start,
        collect_end -
          1L,
        by = "day"
      )
  ),
  by = .(
    period_id,
    loc_num
  )
]

cat(
  "filter-period全部记录 = ",
  format(
    N_FP_ALL,
    big.mark = ","
  ),
  "\n正式filter-period样本 = ",
  format(
    nrow(FP_BASE),
    big.mark = ","
  ),
  "\n正式filter-period监测站 = ",
  uniqueN(
    FP_BASE$loc_num
  ),
  "\n",
  sep = ""
)

standardize_base <- function(
  D,
  branch
) {

  X <- copy(D)

  if (
    branch ==
      "quasi_daily"
  ) {

    X[
      ,
      `:=`(
        analysis_date =
          as.IDate(
            date
          ),
        pm_value =
          safe_num(
            pm_same_day
          ),
        time_fe =
          factor(
            as.character(
              date
            )
          )
      )
    ]

  } else {

    X[
      ,
      `:=`(
        analysis_date =
          as.IDate(
            collect_start
          ),
        pm_value =
          safe_num(
            pm_filter_period_official
          ),
        time_fe =
          interaction(
            collect_start,
            collect_end,
            drop = TRUE,
            lex.order = TRUE
          )
      )
    ]
  }

  X[
    ,
    `:=`(
      gross_beta_mBq_m3 =
        safe_num(
          gross_beta_mBq_m3
        ),
      c_prop =
        safe_num(
          c_prop
        ),
      surface_u_ppm_article =
        safe_num(
          surface_u_ppm_article
        ),
      sunspot_monthly =
        safe_num(
          sunspot_monthly
        ),
      inv_hpbl =
        safe_num(
          inv_hpbl
        ),
      rhum =
        safe_num(
          rhum
        ),
      soilm =
        safe_num(
          soilm
        ),
      vel =
        safe_num(
          vel
        ),
      lat_diff_35 =
        safe_num(
          lat_diff_35
        ),
      year_model =
        safe_num(
          year_model
        ),
      temp_centered =
        safe_num(
          temp
        ) -
        273.15
    )
  ]

  X[
    ,
    calendar_year :=
      as.integer(
        format(
          analysis_date,
          "%Y"
        )
      )
  ]

  X[
    ,
    monitor_year_fe :=
      interaction(
        loc_num,
        calendar_year,
        drop = TRUE,
        lex.order = TRUE
      )
  ]

  X
}

Q_BASE <- standardize_base(
  Q_BASE,
  "quasi_daily"
)

FP_BASE <- standardize_base(
  FP_BASE,
  "filter_period"
)

MODEL_BASE_COLS <- c(
  "gross_beta_mBq_m3",
  "c_prop",
  "surface_u_ppm_article",
  "sunspot_monthly",
  "pm_value",
  "inv_hpbl",
  "rhum",
  "soilm",
  "vel",
  "lat_diff_35",
  "year_model",
  "temp_centered",
  "loc_num",
  "calendar_year",
  "monitor_year_fe",
  "time_fe"
)

make_single_frame <- function(
  D,
  exposure_col,
  scale_value = 1
) {

  if (
    !exposure_col %in%
      names(D)
  ) {

    stopf(
      "缺少暴露字段：",
      exposure_col
    )
  }

  d <- D[
    ,
    c(
      MODEL_BASE_COLS,
      exposure_col
    ),
    with = FALSE
  ]

  d[
    ,
    exposure :=
      safe_num(
        get(
          exposure_col
        )
      ) /
      scale_value
  ]

  req <- c(
    MODEL_BASE_COLS,
    "exposure"
  )

  d <- d[
    complete.cases(
      d[
        ,
        ..req
      ]
    )
  ]

  d[
    ,
    pm_centered :=
      pm_value -
      10
  ]

  d
}

within_sd <- function(
  d,
  col
) {

  if (
    nrow(d) ==
      0L
  ) {
    return(
      NA_real_
    )
  }

  z <- d[
    ,
    get(col) -
      mean(
        get(col)
      ),
    by =
      monitor_year_fe
  ]$V1

  sd(
    z,
    na.rm = TRUE
  )
}

fe_rhs_covariates <- function(
  include_pm = TRUE,
  strict_time = FALSE,
  current_wind_name = "vel"
) {

  parts <- c(
    "c_prop"
  )

  if (
    !strict_time
  ) {
    parts <- c(
      parts,
      "sunspot_monthly"
    )
  }

  if (
    include_pm
  ) {
    parts <- c(
      parts,
      "pm_centered"
    )
  }

  parts <- c(
    parts,
    "inv_hpbl",
    "rhum",
    "soilm",
    current_wind_name,
    "temp_centered",
    "I(temp_centered^2)",
    "I(temp_centered^3)",
    "sunspot_monthly:lat_diff_35",
    "surface_u_ppm_article:soilm"
  )

  paste(
    unique(parts),
    collapse = " + "
  )
}

fit_single_fe <- function(
  D,
  exposure_col,
  scale_value,
  metadata,
  strict_time = FALSE,
  include_pm = TRUE,
  make_curve = FALSE,
  make_distribution_summary = FALSE
) {

  d <- make_single_frame(
    D,
    exposure_col,
    scale_value
  )

  meta <- copy(
    metadata
  )

  meta[
    ,
    `:=`(
      exposure_col =
        exposure_col,
      scale_value =
        scale_value,
      n_rows_pre =
        nrow(d),
      n_monitors =
        uniqueN(
          d$loc_num
        ),
      n_monitor_years =
        uniqueN(
          d$monitor_year_fe
        ),
      within_sd =
        within_sd(
          d,
          "exposure"
        )
    )
  ]

  if (
    nrow(d) <
      1000L ||
    uniqueN(
      d$loc_num
    ) <
      2L ||
    !is.finite(
      meta$within_sd
    ) ||
    meta$within_sd <=
      0
  ) {

    return(
      list(
        status =
          cbind(
            meta,
            data.table(
              status =
                "SKIP",
              nobs =
                NA_integer_,
              note =
                "同站同年内部暴露变异不足"
            )
          ),
        effects =
          data.table(),
        curve =
          data.table(),
        wind_summary =
          data.table()
      )
    )
  }

  rhs_cov <- fe_rhs_covariates(
    include_pm =
      include_pm,
    strict_time =
      strict_time,
    current_wind_name =
      "vel"
  )

  rhs <- paste0(
    "exposure + exposure:vel + ",
    rhs_cov
  )

  fe_part <- if (
    strict_time
  ) {
    "monitor_year_fe + time_fe"
  } else {
    "monitor_year_fe"
  }

  FORM <- as.formula(
    paste0(
      "gross_beta_mBq_m3 ~ ",
      rhs,
      " | ",
      fe_part
    )
  )

  fit <- tryCatch(
    feols(
      FORM,
      data = d,
      vcov = ~loc_num,
      notes = FALSE,
      warn = FALSE
    ),
    error =
      function(e) {
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
      list(
        status =
          cbind(
            meta,
            data.table(
              status =
                "FAIL",
              nobs =
                NA_integer_,
              note =
                fit$message
            )
          ),
        effects =
          data.table(),
        curve =
          data.table(),
        wind_summary =
          data.table()
      )
    )
  }

  beta <- coef(
    fit
  )

  V <- vcov(
    fit
  )

  int_term <- find_interaction(
    names(beta),
    "exposure",
    "vel"
  )

  if (
    !"exposure" %in%
      names(beta) ||
    is.na(
      int_term
    )
  ) {

    return(
      list(
        status =
          cbind(
            meta,
            data.table(
              status =
                "SKIP",
              nobs =
                nobs(fit),
              note =
                "暴露主效应或暴露×风速项被共线性删除"
            )
          ),
        effects =
          data.table(),
        curve =
          data.table(),
        wind_summary =
          data.table()
      )
    )
  }

  REF <- reference_winds(
    d$vel
  )

  EFFECTS <- rbindlist(
    lapply(
      seq_len(
        nrow(REF)
      ),
      function(ii) {

        eff <- effect_at_wind(
          beta,
          V,
          "exposure",
          int_term,
          REF$wind_ms[ii]
        )

        cbind(
          meta[
            ,
            setdiff(
              names(meta),
              c(
                "n_rows_pre",
                "n_monitors",
                "n_monitor_years",
                "within_sd"
              )
            ),
            with = FALSE
          ],
          data.table(
            wind_reference =
              REF$wind_reference[ii],
            wind_ms =
              REF$wind_ms[ii],
            estimate =
              unname(
                eff["estimate"]
              ),
            se =
              unname(
                eff["se"]
              ),
            ci95_low =
              unname(
                eff["ci_low"]
              ),
            ci95_high =
              unname(
                eff["ci_high"]
              ),
            p_wald =
              unname(
                eff["p"]
              ),
            direction =
              direction_text(
                unname(
                  eff["estimate"]
                )
              ),
            ci_status =
              ci_text(
                unname(
                  eff["ci_low"]
                ),
                unname(
                  eff["ci_high"]
                )
              )
          )
        )
      }
    ),
    use.names = TRUE,
    fill = TRUE
  )

  CURVE <- data.table()

  if (
    make_curve
  ) {

    q01 <- unname(
      quantile(
        d$vel,
        0.01
      )
    )

    q99 <- unname(
      quantile(
        d$vel,
        0.99
      )
    )

    grid <- seq(
      q01,
      q99,
      length.out =
        201
    )

    CURVE <- rbindlist(
      lapply(
        grid,
        function(w) {

          eff <- effect_at_wind(
            beta,
            V,
            "exposure",
            int_term,
            w
          )

          cbind(
            meta[
              ,
              setdiff(
                names(meta),
                c(
                  "n_rows_pre",
                  "n_monitors",
                  "n_monitor_years",
                  "within_sd"
                )
              ),
              with = FALSE
            ],
            data.table(
              wind_ms =
                w,
              estimate =
                unname(
                  eff["estimate"]
                ),
              ci95_low =
                unname(
                  eff["ci_low"]
                ),
              ci95_high =
                unname(
                  eff["ci_high"]
                ),
              direction =
                direction_text(
                  unname(
                    eff["estimate"]
                  )
                ),
              ci_status =
                ci_text(
                  unname(
                    eff["ci_low"]
                  ),
                  unname(
                    eff["ci_high"]
                  )
                )
            )
          )
        }
      ),
      use.names = TRUE,
      fill = TRUE
    )
  }

  WIND_SUM <- data.table()

  if (
    make_distribution_summary
  ) {

    b0 <- beta[["exposure"]]
    bi <- beta[[int_term]]

    w <- d$vel

    est_vec <-
      b0 +
      bi *
      w

    var_vec <-
      V[
        "exposure",
        "exposure"
      ] +
      w^2 *
      V[
        int_term,
        int_term
      ] +
      2 *
      w *
      V[
        "exposure",
        int_term
      ]

    se_vec <- sqrt(
      pmax(
        var_vec,
        0
      )
    )

    low_vec <-
      est_vec -
      1.96 *
      se_vec

    high_vec <-
      est_vec +
      1.96 *
      se_vec

    zero_cross <- if (
      is.finite(bi) &&
      abs(bi) >
        .Machine$double.eps
    ) {
      -b0 / bi
    } else {
      NA_real_
    }

    WIND_SUM <- cbind(
      meta[
        ,
        setdiff(
          names(meta),
          c(
            "n_rows_pre",
            "n_monitors",
            "n_monitor_years",
            "within_sd"
          )
        ),
        with = FALSE
      ],
      data.table(
        wind_min =
          min(
            w
          ),
        wind_q01 =
          unname(
            quantile(
              w,
              0.01
            )
          ),
        wind_q25 =
          unname(
            quantile(
              w,
              0.25
            )
          ),
        wind_median =
          median(w),
        wind_mean =
          mean(w),
        wind_q75 =
          unname(
            quantile(
              w,
              0.75
            )
          ),
        wind_q99 =
          unname(
            quantile(
              w,
              0.99
            )
          ),
        wind_max =
          max(w),
        zero_crossing_wind_ms =
          zero_cross,
        zero_crossing_inside_q01_q99 =
          is.finite(
            zero_cross
          ) &&
          zero_cross >=
            unname(
              quantile(
                w,
                0.01
              )
            ) &&
          zero_cross <=
            unname(
              quantile(
                w,
                0.99
              )
            ),
        pct_observations_point_effect_positive =
          100 *
          mean(
            est_vec >
              0
          ),
        pct_observations_point_effect_negative =
          100 *
          mean(
            est_vec <
              0
          ),
        pct_observations_CI_above_zero =
          100 *
          mean(
            low_vec >
              0
          ),
        pct_observations_CI_below_zero =
          100 *
          mean(
            high_vec <
              0
          ),
        pct_observations_CI_include_zero =
          100 *
          mean(
            low_vec <=
              0 &
            high_vec >=
              0
          )
      )
    )
  }

  STATUS <- cbind(
    meta,
    data.table(
      status =
        "PASS",
      nobs =
        nobs(fit),
      note =
        if (
          strict_time
        ) {
          "监测站×年份 + 共同时间固定效应；标准误按监测站聚类"
        } else {
          "监测站×年份固定效应；标准误按监测站聚类"
        }
    )
  )

  list(
    status =
      STATUS,
    effects =
      EFFECTS,
    curve =
      CURVE,
    wind_summary =
      WIND_SUM
  )
}

fit_direction_joint <- function(
  D,
  up_col,
  down_col,
  scale_value,
  metadata
) {

  req <- c(
    MODEL_BASE_COLS,
    up_col,
    down_col
  )

  d <- D[
    ,
    ..req
  ]

  d[
    ,
    `:=`(
      up_exposure =
        safe_num(
          get(
            up_col
          )
        ) /
        scale_value,
      down_exposure =
        safe_num(
          get(
            down_col
          )
        ) /
        scale_value
    )
  ]

  d <- d[
    complete.cases(
      d[
        ,
        c(
          MODEL_BASE_COLS,
          "up_exposure",
          "down_exposure"
        ),
        with = FALSE
      ]
    )
  ]

  d[
    ,
    pm_centered :=
      pm_value -
      10
  ]

  up_sd <- within_sd(
    d,
    "up_exposure"
  )

  down_sd <- within_sd(
    d,
    "down_exposure"
  )

  meta <- copy(
    metadata
  )

  meta[
    ,
    `:=`(
      up_col =
        up_col,
      down_col =
        down_col,
      n_rows_pre =
        nrow(d),
      up_within_sd =
        up_sd,
      down_within_sd =
        down_sd
    )
  ]

  if (
    nrow(d) <
      1000L ||
    !is.finite(
      up_sd
    ) ||
    !is.finite(
      down_sd
    ) ||
    up_sd <=
      0 ||
    down_sd <=
      0
  ) {

    return(
      list(
        status =
          cbind(
            meta,
            data.table(
              status =
                "SKIP",
              note =
                "上风向或下风向同站同年内部变异不足"
            )
          ),
        effects =
          data.table()
      )
    )
  }

  FORM <- as.formula(
    paste0(
      "gross_beta_mBq_m3 ~ ",
      "up_exposure + down_exposure + ",
      "up_exposure:vel + down_exposure:vel + ",
      fe_rhs_covariates(
        include_pm =
          TRUE,
        strict_time =
          FALSE,
        current_wind_name =
          "vel"
      ),
      " | monitor_year_fe"
    )
  )

  fit <- tryCatch(
    feols(
      FORM,
      data = d,
      vcov = ~loc_num,
      notes = FALSE,
      warn = FALSE
    ),
    error =
      function(e) {
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
      list(
        status =
          cbind(
            meta,
            data.table(
              status =
                "FAIL",
              note =
                fit$message
            )
          ),
        effects =
          data.table()
      )
    )
  }

  beta <- coef(
    fit
  )

  V <- vcov(
    fit
  )

  up_int <- find_interaction(
    names(beta),
    "up_exposure",
    "vel"
  )

  down_int <- find_interaction(
    names(beta),
    "down_exposure",
    "vel"
  )

  needed <- c(
    "up_exposure",
    "down_exposure"
  )

  if (
    any(
      !needed %in%
        names(beta)
    ) ||
    is.na(up_int) ||
    is.na(down_int)
  ) {

    return(
      list(
        status =
          cbind(
            meta,
            data.table(
              status =
                "SKIP",
              note =
                "上/下风向项被共线性删除"
            )
          ),
        effects =
          data.table()
      )
    )
  }

  REF <- reference_winds(
    d$vel
  )

  EFFECTS <- rbindlist(
    lapply(
      seq_len(
        nrow(REF)
      ),
      function(ii) {

        w <- REF$wind_ms[ii]

        up_eff <- effect_at_wind(
          beta,
          V,
          "up_exposure",
          up_int,
          w
        )

        down_eff <- effect_at_wind(
          beta,
          V,
          "down_exposure",
          down_int,
          w
        )

        diff_eff <- linear_combination(
          beta,
          V,
          setNames(
            c(
              1,
              w,
              -1,
              -w
            ),
            c(
              "up_exposure",
              up_int,
              "down_exposure",
              down_int
            )
          )
        )

        rbindlist(
          list(
            cbind(
              metadata,
              data.table(
                component =
                  "上风向效应",
                wind_reference =
                  REF$wind_reference[ii],
                wind_ms =
                  w,
                estimate =
                  unname(
                    up_eff["estimate"]
                  ),
                se =
                  unname(
                    up_eff["se"]
                  ),
                ci95_low =
                  unname(
                    up_eff["ci_low"]
                  ),
                ci95_high =
                  unname(
                    up_eff["ci_high"]
                  ),
                p_wald =
                  unname(
                    up_eff["p"]
                  ),
                direction =
                  direction_text(
                    unname(
                      up_eff["estimate"]
                    )
                  ),
                ci_status =
                  ci_text(
                    unname(
                      up_eff["ci_low"]
                    ),
                    unname(
                      up_eff["ci_high"]
                    )
                  )
              )
            ),
            cbind(
              metadata,
              data.table(
                component =
                  "下风向效应",
                wind_reference =
                  REF$wind_reference[ii],
                wind_ms =
                  w,
                estimate =
                  unname(
                    down_eff["estimate"]
                  ),
                se =
                  unname(
                    down_eff["se"]
                  ),
                ci95_low =
                  unname(
                    down_eff["ci_low"]
                  ),
                ci95_high =
                  unname(
                    down_eff["ci_high"]
                  ),
                p_wald =
                  unname(
                    down_eff["p"]
                  ),
                direction =
                  direction_text(
                    unname(
                      down_eff["estimate"]
                    )
                  ),
                ci_status =
                  ci_text(
                    unname(
                      down_eff["ci_low"]
                    ),
                    unname(
                      down_eff["ci_high"]
                    )
                  )
              )
            ),
            cbind(
              metadata,
              data.table(
                component =
                  "上风向减下风向",
                wind_reference =
                  REF$wind_reference[ii],
                wind_ms =
                  w,
                estimate =
                  unname(
                    diff_eff["estimate"]
                  ),
                se =
                  unname(
                    diff_eff["se"]
                  ),
                ci95_low =
                  unname(
                    diff_eff["ci_low"]
                  ),
                ci95_high =
                  unname(
                    diff_eff["ci_high"]
                  ),
                p_wald =
                  unname(
                    diff_eff["p"]
                  ),
                direction =
                  direction_text(
                    unname(
                      diff_eff["estimate"]
                    )
                  ),
                ci_status =
                  ci_text(
                    unname(
                      diff_eff["ci_low"]
                    ),
                    unname(
                      diff_eff["ci_high"]
                    )
                  )
              )
            )
          ),
          use.names = TRUE,
          fill = TRUE
        )
      }
    ),
    use.names = TRUE,
    fill = TRUE
  )

  list(
    status =
      cbind(
        meta,
        data.table(
          status =
            "PASS",
          note =
            "同一固定效应模型同时纳入上风向和下风向，并直接检验上风向减下风向"
        )
      ),
    effects =
      EFFECTS
  )
}

fit_count_heat_joint <- function(
  D,
  count_col,
  heat_col,
  metadata
) {

  req <- c(
    MODEL_BASE_COLS,
    count_col,
    heat_col
  )

  d <- D[
    ,
    ..req
  ]

  d[
    ,
    `:=`(
      count_exposure =
        safe_num(
          get(
            count_col
          )
        ),
      heat_exposure =
        safe_num(
          get(
            heat_col
          )
        ) /
        10000
    )
  ]

  d <- d[
    complete.cases(
      d[
        ,
        c(
          MODEL_BASE_COLS,
          "count_exposure",
          "heat_exposure"
        ),
        with = FALSE
      ]
    )
  ]

  d[
    ,
    pm_centered :=
      pm_value -
      10
  ]

  count_sd <- within_sd(
    d,
    "count_exposure"
  )

  heat_sd <- within_sd(
    d,
    "heat_exposure"
  )

  meta <- copy(
    metadata
  )

  meta[
    ,
    `:=`(
      count_col =
        count_col,
      heat_col =
        heat_col,
      n_rows_pre =
        nrow(d),
      count_within_sd =
        count_sd,
      heat_within_sd =
        heat_sd
    )
  ]

  if (
    nrow(d) <
      1000L ||
    !is.finite(
      count_sd
    ) ||
    !is.finite(
      heat_sd
    ) ||
    count_sd <=
      0 ||
    heat_sd <=
      0
  ) {

    return(
      list(
        status =
          cbind(
            meta,
            data.table(
              status =
                "SKIP",
              note =
                "Count或Heat同站同年内部变异不足"
            )
          ),
        effects =
          data.table(),
        collinearity =
          data.table()
      )
    )
  }

  WITHIN <- d[
    ,
    .(
      count_within =
        count_exposure -
        mean(
          count_exposure
        ),
      heat_within =
        heat_exposure -
        mean(
          heat_exposure
        )
    ),
    by =
      monitor_year_fe
  ]

  rr <- suppressWarnings(
    cor(
      WITHIN$count_within,
      WITHIN$heat_within,
      use = "complete.obs"
    )
  )

  vif <- if (
    is.finite(rr) &&
    abs(rr) <
      1
  ) {
    1 /
      (
        1 -
        rr^2
      )
  } else if (
    is.finite(rr)
  ) {
    Inf
  } else {
    NA_real_
  }

  COLL <- cbind(
    metadata,
    data.table(
      within_count_heat_correlation =
        rr,
      approximate_two_predictor_VIF =
        vif
    )
  )

  FORM <- as.formula(
    paste0(
      "gross_beta_mBq_m3 ~ ",
      "count_exposure + heat_exposure + ",
      "count_exposure:vel + heat_exposure:vel + ",
      fe_rhs_covariates(
        include_pm =
          TRUE,
        strict_time =
          FALSE,
        current_wind_name =
          "vel"
      ),
      " | monitor_year_fe"
    )
  )

  fit <- tryCatch(
    feols(
      FORM,
      data = d,
      vcov = ~loc_num,
      notes = FALSE,
      warn = FALSE
    ),
    error =
      function(e) {
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
      list(
        status =
          cbind(
            meta,
            data.table(
              status =
                "FAIL",
              note =
                fit$message
            )
          ),
        effects =
          data.table(),
        collinearity =
          COLL
      )
    )
  }

  beta <- coef(
    fit
  )

  V <- vcov(
    fit
  )

  count_int <- find_interaction(
    names(beta),
    "count_exposure",
    "vel"
  )

  heat_int <- find_interaction(
    names(beta),
    "heat_exposure",
    "vel"
  )

  if (
    !"count_exposure" %in%
      names(beta) ||
    !"heat_exposure" %in%
      names(beta) ||
    is.na(count_int) ||
    is.na(heat_int)
  ) {

    return(
      list(
        status =
          cbind(
            meta,
            data.table(
              status =
                "SKIP",
              note =
                "联合模型Count或Heat项被共线性删除"
            )
          ),
        effects =
          data.table(),
        collinearity =
          COLL
      )
    )
  }

  REF <- reference_winds(
    d$vel
  )

  EFFECTS <- rbindlist(
    lapply(
      seq_len(
        nrow(REF)
      ),
      function(ii) {

        w <- REF$wind_ms[ii]

        ec <- effect_at_wind(
          beta,
          V,
          "count_exposure",
          count_int,
          w
        )

        eh <- effect_at_wind(
          beta,
          V,
          "heat_exposure",
          heat_int,
          w
        )

        rbindlist(
          list(
            cbind(
              metadata,
              data.table(
                component =
                  "联合模型_Count",
                unit =
                  "总Heat相同时，每增加1座活动电厂",
                wind_reference =
                  REF$wind_reference[ii],
                wind_ms =
                  w,
                estimate =
                  unname(
                    ec["estimate"]
                  ),
                se =
                  unname(
                    ec["se"]
                  ),
                ci95_low =
                  unname(
                    ec["ci_low"]
                  ),
                ci95_high =
                  unname(
                    ec["ci_high"]
                  ),
                p_wald =
                  unname(
                    ec["p"]
                  ),
                direction =
                  direction_text(
                    unname(
                      ec["estimate"]
                    )
                  ),
                ci_status =
                  ci_text(
                    unname(
                      ec["ci_low"]
                    ),
                    unname(
                      ec["ci_high"]
                    )
                  )
              )
            ),
            cbind(
              metadata,
              data.table(
                component =
                  "联合模型_Heat",
                unit =
                  "Count相同时，每增加10000 MMBtu实际热输入",
                wind_reference =
                  REF$wind_reference[ii],
                wind_ms =
                  w,
                estimate =
                  unname(
                    eh["estimate"]
                  ),
                se =
                  unname(
                    eh["se"]
                  ),
                ci95_low =
                  unname(
                    eh["ci_low"]
                  ),
                ci95_high =
                  unname(
                    eh["ci_high"]
                  ),
                p_wald =
                  unname(
                    eh["p"]
                  ),
                direction =
                  direction_text(
                    unname(
                      eh["estimate"]
                    )
                  ),
                ci_status =
                  ci_text(
                    unname(
                      eh["ci_low"]
                    ),
                    unname(
                      eh["ci_high"]
                    )
                  )
              )
            )
          ),
          use.names = TRUE,
          fill = TRUE
        )
      }
    ),
    use.names = TRUE,
    fill = TRUE
  )

  list(
    status =
      cbind(
        meta,
        data.table(
          status =
            "PASS",
          note =
            "同燃料Count与Heat联合固定效应模型"
        )
      ),
    effects =
      EFFECTS,
    collinearity =
      COLL
  )
}

fit_mundlak <- function(
  D,
  exposure_col,
  scale_value,
  metadata
) {

  d <- make_single_frame(
    D,
    exposure_col,
    scale_value
  )

  if (
    nrow(d) <
      1000L
  ) {

    return(
      list(
        status =
          cbind(
            metadata,
            data.table(
              status =
                "SKIP",
              note =
                "样本不足"
            )
          ),
        effects =
          data.table()
      )
    )
  }

  d[
    ,
    exposure_between :=
      mean(
        exposure
      ),
    by =
      monitor_year_fe
  ]

  d[
    ,
    exposure_within :=
      exposure -
      exposure_between
  ]

  ws <- sd(
    d$exposure_within,
    na.rm = TRUE
  )

  bs <- sd(
    d$exposure_between,
    na.rm = TRUE
  )

  if (
    !is.finite(ws) ||
    ws <= 0 ||
    !is.finite(bs)
  ) {

    return(
      list(
        status =
          cbind(
            metadata,
            data.table(
              status =
                "SKIP",
              note =
                "站点年份内或站点年份间暴露变异不足"
            )
          ),
        effects =
          data.table()
      )
    )
  }

  d[
    ,
    loc_num_factor :=
      factor(
        loc_num
      )
  ]

  FORM <- as.formula(
    paste0(
      "gross_beta_mBq_m3 ~ ",
      "exposure_within + exposure_between + ",
      "c_prop + ",
      "surface_u_ppm_article + ",
      "sunspot_monthly + ",
      "pm_centered + ",
      "inv_hpbl + ",
      "rhum + ",
      "soilm + ",
      "vel + ",
      "lat_diff_35 + ",
      "year_model + I(year_model^2) + ",
      "temp_centered + I(temp_centered^2) + I(temp_centered^3) + ",
      "sunspot_monthly:lat_diff_35 + ",
      "surface_u_ppm_article:soilm + ",
      "exposure_within:vel + ",
      "exposure_between:vel + ",
      "(1 | loc_num_factor)"
    )
  )

  warning_text <- character()

  fit <- tryCatch(
    withCallingHandlers(
      lmer(
        FORM,
        data = d,
        REML = TRUE,
        na.action = na.fail,
        control =
          make_lmer_control()
      ),
      warning =
        function(w) {

          warning_text <<-
            c(
              warning_text,
              conditionMessage(w)
            )

          invokeRestart(
            "muffleWarning"
          )
        }
    ),
    error =
      function(e) {
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
      list(
        status =
          cbind(
            metadata,
            data.table(
              status =
                "FAIL",
              note =
                fit$message
            )
          ),
        effects =
          data.table()
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

  wi <- find_interaction(
    names(beta),
    "exposure_within",
    "vel"
  )

  bi <- find_interaction(
    names(beta),
    "exposure_between",
    "vel"
  )

  if (
    !"exposure_within" %in%
      names(beta) ||
    !"exposure_between" %in%
      names(beta) ||
    is.na(wi) ||
    is.na(bi)
  ) {

    return(
      list(
        status =
          cbind(
            metadata,
            data.table(
              status =
                "SKIP",
              note =
                "站点年份内/间项被秩缺陷删除"
            )
          ),
        effects =
          data.table()
      )
    )
  }

  REF <- reference_winds(
    d$vel
  )

  EFFECTS <- rbindlist(
    lapply(
      c(
        "站点年份内变化",
        "站点年份间平均差异"
      ),
      function(component) {

        if (
          component ==
            "站点年份内变化"
        ) {

          main_term <-
            "exposure_within"

          int_term <-
            wi

        } else {

          main_term <-
            "exposure_between"

          int_term <-
            bi
        }

        rbindlist(
          lapply(
            seq_len(
              nrow(REF)
            ),
            function(ii) {

              eff <- effect_at_wind(
                beta,
                V,
                main_term,
                int_term,
                REF$wind_ms[ii]
              )

              cbind(
                metadata,
                data.table(
                  component =
                    component,
                  wind_reference =
                    REF$wind_reference[ii],
                  wind_ms =
                    REF$wind_ms[ii],
                  estimate =
                    unname(
                      eff["estimate"]
                    ),
                  se =
                    unname(
                      eff["se"]
                    ),
                  ci95_low =
                    unname(
                      eff["ci_low"]
                    ),
                  ci95_high =
                    unname(
                      eff["ci_high"]
                    ),
                  p_wald =
                    unname(
                      eff["p"]
                    ),
                  direction =
                    direction_text(
                      unname(
                        eff["estimate"]
                      )
                    ),
                  ci_status =
                    ci_text(
                      unname(
                        eff["ci_low"]
                      ),
                      unname(
                        eff["ci_high"]
                      )
                    )
                )
              )
            }
          ),
          use.names = TRUE,
          fill = TRUE
        )
      }
    ),
    use.names = TRUE,
    fill = TRUE
  )

  list(
    status =
      cbind(
        metadata,
        data.table(
          status =
            "PASS",
          singular =
            as.integer(
              isSingular(
                fit,
                tol = 1e-4
              )
            ),
          warnings =
            collapse_msg(
              warning_text
            ),
          note =
            "站点年份内/间分解随机截距诊断"
        )
      ),
    effects =
      EFFECTS
  )
}

find_weighted_files <- function() {

  all_files <- list.files(
    H04_DIR,
    pattern =
      "\\.csv\\.gz$",
    recursive =
      TRUE,
    full.names =
      TRUE
  )

  out <- character(
    length(YEARS)
  )

  names(out) <- YEARS

  for (yy in YEARS) {

    cand <- all_files[
      grepl(
        as.character(yy),
        basename(
          all_files
        ),
        fixed = TRUE
      )
    ]

    if (
      length(cand) ==
        0L
    ) {

      stopf(
        "07H04中找不到",
        yy,
        "年的csv.gz。"
      )
    }

    scores <- vapply(
      cand,
      function(f) {

        H <- tryCatch(
          names(
            fread(
              file = f,
              nrows = 0,
              showProgress = FALSE
            )
          ),
          error =
            function(e)
              character()
        )

        if (
          !all(
            c(
              "loc_num",
              "date"
            ) %in%
              H
          )
        ) {
          return(-1)
        }

        sum(
          grepl(
            "proximity",
            H,
            ignore.case = TRUE
          )
        )
      },
      numeric(1)
    )

    if (
      max(
        scores
      ) <=
        0
    ) {

      stopf(
        "07H04中",
        yy,
        "年没有找到包含proximity暴露的逐日文件。"
      )
    }

    out[
      as.character(yy)
    ] <- cand[
      which.max(
        scores
      )
    ]
  }

  out
}

WEIGHTED_FILES <- find_weighted_files()

WEIGHTED_HEADER <- names(
  fread(
    file =
      WEIGHTED_FILES[1],
    nrows = 0,
    showProgress = FALSE
  )
)

find_weighted_col <- function(
  fuel,
  direction,
  radius,
  metric
) {

  H <- WEIGHTED_HEADER

  exact_suffix <- paste0(
    "_",
    radius,
    "km_90deg"
  )

  base <- H[
    grepl(
      paste0(
        "^campd_",
        fuel,
        "_",
        direction,
        "_"
      ),
      H
    ) &
    endsWith(
      H,
      exact_suffix
    ) &
    grepl(
      "proximity",
      H,
      ignore.case = TRUE
    )
  ]

  if (
    metric ==
      "proximity"
  ) {

    z <- base[
      !grepl(
        "heat",
        base,
        ignore.case = TRUE
      )
    ]

  } else if (
    metric ==
      "heat_proximity"
  ) {

    z <- base[
      grepl(
        "heat",
        base,
        ignore.case = TRUE
      )
    ]

  } else {

    stopf(
      "未知距离加权metric：",
      metric
    )
  }

  if (
    length(z) !=
      1L
  ) {

    stopf(
      "无法唯一识别07H04距离加权字段：fuel=",
      fuel,
      "; direction=",
      direction,
      "; radius=",
      radius,
      "; metric=",
      metric,
      "\n精确半径后缀：",
      exact_suffix,
      "\n候选字段：",
      paste(
        z,
        collapse = " | "
      ),
      "\n该燃料/方向下全部proximity字段：",
      paste(
        H[
          grepl(
            paste0(
              "^campd_",
              fuel,
              "_",
              direction,
              "_"
            ),
            H
          ) &
          grepl(
            "proximity",
            H,
            ignore.case = TRUE
          )
        ],
        collapse = " | "
      )
    )
  }

  z[1]
}

WEIGHTED_MAP <- rbindlist(
  lapply(
    RADII,
    function(r) {

      rbindlist(
        lapply(
          FUELS,
          function(fuel) {

            rbindlist(
              lapply(
                DIRECTIONS,
                function(direction) {

                  data.table(
                    radius_km =
                      r,
                    fuel =
                      fuel,
                    direction =
                      direction,
                    proximity_col =
                      find_weighted_col(
                        fuel,
                        direction,
                        r,
                        "proximity"
                      ),
                    heat_proximity_col =
                      find_weighted_col(
                        fuel,
                        direction,
                        r,
                        "heat_proximity"
                      )
                  )
                }
              )
            )
          }
        )
      )
    }
  )
)

fwrite(
  WEIGHTED_MAP,
  file.path(
    DIR_WEIGHTED,
    "00_07H04距离加权字段映射.csv"
  )
)

read_raw_radius <- function(
  radius
) {

  cols <- c(
    "loc_num",
    "date"
  )

  for (fuel in FUELS) {

    for (
      direction in
      DIRECTIONS
    ) {

      cols <- c(
        cols,
        raw_col(
          fuel,
          direction,
          "count",
          radius
        ),
        raw_col(
          fuel,
          direction,
          "heat",
          radius
        )
      )
    }
  }

  out <- rbindlist(
    lapply(
      RAW_FILES,
      function(f) {

        H <- names(
          fread(
            file = f,
            nrows = 0,
            showProgress = FALSE
          )
        )

        miss <- setdiff(
          cols,
          H
        )

        if (
          length(miss) >
            0L
        ) {

          stopf(
            basename(f),
            "缺少",
            radius,
            "km字段：",
            paste(
              miss,
              collapse = ", "
            )
          )
        }

        fread(
          file = f,
          select = cols,
          showProgress = FALSE
        )
      }
    ),
    use.names = TRUE
  )

  out[
    ,
    `:=`(
      loc_num =
        as.character(
          loc_num
        ),
      date =
        as.IDate(
          date
        )
    )
  ]

  if (
    out[
      ,
      anyDuplicated(
        paste(
          loc_num,
          date
        )
      )
    ] >
      0L
  ) {

    stopf(
      radius,
      "km原始逐日暴露存在重复monitor-day。"
    )
  }

  out
}

read_weighted_radius <- function(
  radius
) {

  M <- WEIGHTED_MAP[
    radius_km ==
      radius
  ]

  exposure_cols <- unique(
    c(
      M$proximity_col,
      M$heat_proximity_col
    )
  )

  cols <- c(
    "loc_num",
    "date",
    exposure_cols
  )

  out <- rbindlist(
    lapply(
      WEIGHTED_FILES,
      function(f) {

        H <- names(
          fread(
            file = f,
            nrows = 0,
            showProgress = FALSE
          )
        )

        miss <- setdiff(
          cols,
          H
        )

        if (
          length(miss) >
            0L
        ) {

          stopf(
            basename(f),
            "缺少07H04字段：",
            paste(
              miss,
              collapse = ", "
            )
          )
        }

        fread(
          file = f,
          select = cols,
          showProgress = FALSE
        )
      }
    ),
    use.names = TRUE
  )

  out[
    ,
    `:=`(
      loc_num =
        as.character(
          loc_num
        ),
      date =
        as.IDate(
          date
        )
    )
  ]

  if (
    out[
      ,
      anyDuplicated(
        paste(
          loc_num,
          date
        )
      )
    ] >
      0L
  ) {

    stopf(
      radius,
      "km距离加权逐日暴露存在重复monitor-day。"
    )
  }

  out
}

aggregate_to_filter <- function(
  DAILY,
  exposure_cols
) {

  DD <- merge(
    FP_DAY_MAP,
    DAILY[
      ,
      c(
        "loc_num",
        "date",
        exposure_cols
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

  AG <- DD[
    ,
    lapply(
      .SD,
      strict_mean
    ),
    by =
      period_id,
    .SDcols =
      exposure_cols
  ]

  setnames(
    AG,
    exposure_cols,
    paste0(
      "fpmean_",
      exposure_cols
    )
  )

  merge(
    FP_BASE,
    AG,
    by =
      "period_id",
    all.x =
      TRUE,
    sort =
      FALSE
  )
}

merge_quasi <- function(
  DAILY,
  exposure_cols
) {

  merge(
    Q_BASE,
    DAILY[
      ,
      c(
        "loc_num",
        "date",
        exposure_cols
      ),
      with = FALSE
    ],
    by = c(
      "loc_num",
      "date"
    ),
    all.x = TRUE,
    sort = FALSE,
    allow.cartesian = TRUE
  )
}

fit_lag_one <- function(
  QD,
  lag_col,
  lag_wind_col,
  scale_value,
  metadata
) {

  req <- c(
    MODEL_BASE_COLS,
    lag_col,
    lag_wind_col
  )

  d <- QD[
    ,
    ..req
  ]

  d[
    ,
    `:=`(
      lag_exposure =
        safe_num(
          get(
            lag_col
          )
        ) /
        scale_value,
      lag_wind =
        safe_num(
          get(
            lag_wind_col
          )
        )
    )
  ]

  d <- d[
    complete.cases(
      d[
        ,
        c(
          MODEL_BASE_COLS,
          "lag_exposure",
          "lag_wind"
        ),
        with = FALSE
      ]
    )
  ]

  d[
    ,
    pm_centered :=
      pm_value -
      10
  ]

  ws <- within_sd(
    d,
    "lag_exposure"
  )

  if (
    nrow(d) <
      1000L ||
    !is.finite(ws) ||
    ws <= 0
  ) {

    return(
      list(
        status =
          cbind(
            metadata,
            data.table(
              status =
                "SKIP",
              n_rows =
                nrow(d),
              within_sd =
                ws,
              note =
                "滞后暴露同站同年内部变异不足"
            )
          ),
        effects =
          data.table()
      )
    )
  }

  FORM <- as.formula(
    paste0(
      "gross_beta_mBq_m3 ~ ",
      "lag_exposure + lag_exposure:lag_wind + ",
      "c_prop + ",
      "sunspot_monthly + ",
      "pm_centered + ",
      "inv_hpbl + ",
      "rhum + ",
      "soilm + ",
      "vel + ",
      "temp_centered + ",
      "I(temp_centered^2) + ",
      "I(temp_centered^3) + ",
      "sunspot_monthly:lat_diff_35 + ",
      "surface_u_ppm_article:soilm ",
      "| monitor_year_fe"
    )
  )

  fit <- tryCatch(
    feols(
      FORM,
      data = d,
      vcov = ~loc_num,
      notes = FALSE,
      warn = FALSE
    ),
    error =
      function(e) {
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
      list(
        status =
          cbind(
            metadata,
            data.table(
              status =
                "FAIL",
              n_rows =
                nrow(d),
              within_sd =
                ws,
              note =
                fit$message
            )
          ),
        effects =
          data.table()
      )
    )
  }

  beta <- coef(
    fit
  )

  V <- vcov(
    fit
  )

  int_term <- find_interaction(
    names(beta),
    "lag_exposure",
    "lag_wind"
  )

  if (
    !"lag_exposure" %in%
      names(beta) ||
    is.na(int_term)
  ) {

    return(
      list(
        status =
          cbind(
            metadata,
            data.table(
              status =
                "SKIP",
              n_rows =
                nrow(d),
              within_sd =
                ws,
              note =
                "滞后暴露主效应或滞后风速交互被删除"
            )
          ),
        effects =
          data.table()
      )
    )
  }

  REF <- reference_winds(
    d$lag_wind
  )

  EFFECTS <- rbindlist(
    lapply(
      seq_len(
        nrow(REF)
      ),
      function(ii) {

        eff <- effect_at_wind(
          beta,
          V,
          "lag_exposure",
          int_term,
          REF$wind_ms[ii]
        )

        cbind(
          metadata,
          data.table(
            wind_reference =
              REF$wind_reference[ii],
            wind_ms =
              REF$wind_ms[ii],
            estimate =
              unname(
                eff["estimate"]
              ),
            se =
              unname(
                eff["se"]
              ),
            ci95_low =
              unname(
                eff["ci_low"]
              ),
            ci95_high =
              unname(
                eff["ci_high"]
              ),
            p_wald =
              unname(
                eff["p"]
              ),
            direction =
              direction_text(
                unname(
                  eff["estimate"]
                )
              ),
            ci_status =
              ci_text(
                unname(
                  eff["ci_low"]
                ),
                unname(
                  eff["ci_high"]
                )
              )
          )
        )
      }
    ),
    use.names = TRUE,
    fill = TRUE
  )

  list(
    status =
      cbind(
        metadata,
        data.table(
          status =
            "PASS",
          n_rows =
            nrow(d),
          within_sd =
            ws,
          note =
            "准逐日滞后暴露；方向由暴露日风向确定；交互使用暴露日风速"
        )
      ),
    effects =
      EFFECTS
  )
}

make_lag_daily <- function(
  RAW,
  radius
) {

  D <- merge(
    RAW,
    WIND_LOOKUP,
    by = c(
      "loc_num",
      "date"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  setorder(
    D,
    loc_num,
    date
  )

  exposure_cols <- setdiff(
    names(RAW),
    c(
      "loc_num",
      "date"
    )
  )

  lag_map <- list()

  idx <- 0L

  for (col in exposure_cols) {

    for (
      lag_n in
      0:3
    ) {

      exp_name <- paste0(
        col,
        "__lag",
        lag_n
      )

      wind_name <- paste0(
        col,
        "__lagwind",
        lag_n
      )

      if (
        lag_n ==
          0L
      ) {

        D[
          ,
          (exp_name) :=
            get(col)
        ]

        D[
          ,
          (wind_name) :=
            exposure_day_wind
        ]

      } else {

        D[
          ,
          (exp_name) :=
            shift(
              get(col),
              n =
                lag_n
            ),
          by =
            loc_num
        ]

        D[
          ,
          (wind_name) :=
            shift(
              exposure_day_wind,
              n =
                lag_n
            ),
          by =
            loc_num
        ]

        date_prev_name <- paste0(
          "__date_lag_",
          lag_n
        )

        if (
          !date_prev_name %in%
            names(D)
        ) {

          D[
            ,
            (date_prev_name) :=
              shift(
                date,
                n =
                  lag_n
              ),
            by =
              loc_num
          ]
        }

        bad_gap <- (
          as.integer(
            D$date -
            D[[
              date_prev_name
            ]]
          ) !=
            lag_n
        )

        bad_gap[
          is.na(
            bad_gap
          )
        ] <- TRUE

        set(
          D,
          which(
            bad_gap
          ),
          exp_name,
          NA_real_
        )

        set(
          D,
          which(
            bad_gap
          ),
          wind_name,
          NA_real_
        )
      }
    }

    lag0 <-
      paste0(
        col,
        "__lag0"
      )

    lag1 <-
      paste0(
        col,
        "__lag1"
      )

    lag2 <-
      paste0(
        col,
        "__lag2"
      )

    w0 <-
      paste0(
        col,
        "__lagwind0"
      )

    w1 <-
      paste0(
        col,
        "__lagwind1"
      )

    w2 <-
      paste0(
        col,
        "__lagwind2"
      )

    mean01 <- paste0(
      col,
      "__mean01"
    )

    mean02 <- paste0(
      col,
      "__mean02"
    )

    meanw01 <- paste0(
      col,
      "__meanwind01"
    )

    meanw02 <- paste0(
      col,
      "__meanwind02"
    )

    D[
      ,
      (mean01) :=
        fifelse(
          is.finite(
            get(lag0)
          ) &
          is.finite(
            get(lag1)
          ),
          (
            get(lag0) +
            get(lag1)
          ) /
            2,
          NA_real_
        )
    ]

    D[
      ,
      (mean02) :=
        fifelse(
          is.finite(
            get(lag0)
          ) &
          is.finite(
            get(lag1)
          ) &
          is.finite(
            get(lag2)
          ),
          (
            get(lag0) +
            get(lag1) +
            get(lag2)
          ) /
            3,
          NA_real_
        )
    ]

    D[
      ,
      (meanw01) :=
        fifelse(
          is.finite(
            get(w0)
          ) &
          is.finite(
            get(w1)
          ),
          (
            get(w0) +
            get(w1)
          ) /
            2,
          NA_real_
        )
    ]

    D[
      ,
      (meanw02) :=
        fifelse(
          is.finite(
            get(w0)
          ) &
          is.finite(
            get(w1)
          ) &
          is.finite(
            get(w2)
          ),
          (
            get(w0) +
            get(w1) +
            get(w2)
          ) /
            3,
          NA_real_
        )
    ]

    lag_map[[col]] <- data.table(
      original_col =
        col,
      lag_spec = c(
        "lag0",
        "lag1",
        "lag2",
        "lag3",
        "mean_0_1",
        "mean_0_2"
      ),
      exposure_col = c(
        paste0(
          col,
          "__lag0"
        ),
        paste0(
          col,
          "__lag1"
        ),
        paste0(
          col,
          "__lag2"
        ),
        paste0(
          col,
          "__lag3"
        ),
        mean01,
        mean02
      ),
      wind_col = c(
        paste0(
          col,
          "__lagwind0"
        ),
        paste0(
          col,
          "__lagwind1"
        ),
        paste0(
          col,
          "__lagwind2"
        ),
        paste0(
          col,
          "__lagwind3"
        ),
        meanw01,
        meanw02
      )
    )
  }

  date_lag_cols <- grep(
    "^__date_lag_",
    names(D),
    value = TRUE
  )

  if (
    length(
      date_lag_cols
    ) >
      0L
  ) {

    D[
      ,
      (date_lag_cols) :=
        NULL
    ]
  }

  list(
    data = D,
    map =
      rbindlist(
        lag_map,
        use.names = TRUE
      )
  )
}

process_radius <- function(
  radius
) {

  cat(
    "\n\n############################################################\n",
    "07H10：开始 ",
    radius,
    " km\n",
    "############################################################\n",
    sep = ""
  )

  RAW <- read_raw_radius(
    radius
  )

  RAW_EXPOSURE_COLS <- setdiff(
    names(RAW),
    c(
      "loc_num",
      "date"
    )
  )

  Q_RAW <- merge_quasi(
    RAW,
    RAW_EXPOSURE_COLS
  )

  FP_RAW <- aggregate_to_filter(
    RAW,
    RAW_EXPOSURE_COLS
  )

  MAIN <- checkpoint_run(
    radius,
    "01_main_FE",
    function() {

      statuses <- list()
      effects <- list()
      curves <- list()
      wind_summaries <- list()

      kk <- 0L

      for (branch in BRANCHES) {

        D <- if (
          branch ==
            "quasi_daily"
        ) {
          Q_RAW
        } else {
          FP_RAW
        }

        for (fuel in FUELS) {

          for (
            metric in
            c(
              "count",
              "heat"
            )
          ) {

            kk <- kk +
              1L

            col <- raw_col(
              fuel,
              "upwind",
              metric,
              radius
            )

            if (
              branch ==
                "filter_period"
            ) {
              col <- paste0(
                "fpmean_",
                col
              )
            }

            meta <- data.table(
              radius_km =
                radius,
              branch =
                branch,
              fuel =
                fuel,
              fuel_label =
                unname(
                  FUEL_CN[fuel]
                ),
              metric =
                metric,
              metric_label =
                metric_cn(metric),
              unit =
                metric_unit(metric),
              model_family =
                "主分析_监测站年份固定效应"
            )

            cat(
              "[主分析 ",
              radius,
              "km ",
              kk,
              "/24] ",
              branch,
              " | ",
              fuel,
              " | ",
              metric,
              "\n",
              sep = ""
            )

            z <- fit_single_fe(
              D,
              col,
              metric_scale(
                metric
              ),
              meta,
              strict_time =
                FALSE,
              include_pm =
                TRUE,
              make_curve =
                TRUE,
              make_distribution_summary =
                TRUE
            )

            statuses[[kk]] <-
              z$status

            effects[[kk]] <-
              z$effects

            curves[[kk]] <-
              z$curve

            wind_summaries[[kk]] <-
              z$wind_summary
          }
        }
      }

      list(
        status =
          rbindlist(
            statuses,
            fill = TRUE
          ),
        effects =
          rbindlist(
            effects,
            fill = TRUE
          ),
        curve =
          rbindlist(
            curves,
            fill = TRUE
          ),
        wind_summary =
          rbindlist(
            wind_summaries,
            fill = TRUE
          )
      )
    }
  )

  STRICT <- checkpoint_run(
    radius,
    "02_strict_time_FE",
    function() {

      statuses <- list()
      effects <- list()
      curves <- list()
      wind_summaries <- list()

      kk <- 0L

      for (branch in BRANCHES) {

        D <- if (
          branch ==
            "quasi_daily"
        ) {
          Q_RAW
        } else {
          FP_RAW
        }

        for (fuel in FUELS) {

          for (
            metric in
            c(
              "count",
              "heat"
            )
          ) {

            kk <- kk +
              1L

            col <- raw_col(
              fuel,
              "upwind",
              metric,
              radius
            )

            if (
              branch ==
                "filter_period"
            ) {
              col <- paste0(
                "fpmean_",
                col
              )
            }

            meta <- data.table(
              radius_km =
                radius,
              branch =
                branch,
              fuel =
                fuel,
              fuel_label =
                unname(
                  FUEL_CN[fuel]
                ),
              metric =
                metric,
              metric_label =
                metric_cn(metric),
              unit =
                metric_unit(metric),
              model_family =
                "严格时间固定效应"
            )

            cat(
              "[严格FE ",
              radius,
              "km ",
              kk,
              "/24] ",
              branch,
              " | ",
              fuel,
              " | ",
              metric,
              "\n",
              sep = ""
            )

            z <- fit_single_fe(
              D,
              col,
              metric_scale(
                metric
              ),
              meta,
              strict_time =
                TRUE,
              include_pm =
                TRUE,
              make_curve =
                TRUE,
              make_distribution_summary =
                TRUE
            )

            statuses[[kk]] <-
              z$status

            effects[[kk]] <-
              z$effects

            curves[[kk]] <-
              z$curve

            wind_summaries[[kk]] <-
              z$wind_summary
          }
        }
      }

      list(
        status =
          rbindlist(
            statuses,
            fill = TRUE
          ),
        effects =
          rbindlist(
            effects,
            fill = TRUE
          ),
        curve =
          rbindlist(
            curves,
            fill = TRUE
          ),
        wind_summary =
          rbindlist(
            wind_summaries,
            fill = TRUE
          )
      )
    }
  )

  DIRECTION <- checkpoint_run(
    radius,
    "03_direction_joint",
    function() {

      statuses <- list()
      effects <- list()

      kk <- 0L

      for (branch in BRANCHES) {

        D <- if (
          branch ==
            "quasi_daily"
        ) {
          Q_RAW
        } else {
          FP_RAW
        }

        for (fuel in FUELS) {

          for (
            metric in
            c(
              "count",
              "heat"
            )
          ) {

            kk <- kk +
              1L

            up <- raw_col(
              fuel,
              "upwind",
              metric,
              radius
            )

            down <- raw_col(
              fuel,
              "downwind",
              metric,
              radius
            )

            if (
              branch ==
                "filter_period"
            ) {

              up <- paste0(
                "fpmean_",
                up
              )

              down <- paste0(
                "fpmean_",
                down
              )
            }

            meta <- data.table(
              radius_km =
                radius,
              branch =
                branch,
              fuel =
                fuel,
              fuel_label =
                unname(
                  FUEL_CN[fuel]
                ),
              metric =
                metric,
              metric_label =
                metric_cn(metric),
              unit =
                metric_unit(metric)
            )

            cat(
              "[方向性 ",
              radius,
              "km ",
              kk,
              "/24] ",
              branch,
              " | ",
              fuel,
              " | ",
              metric,
              "\n",
              sep = ""
            )

            z <- fit_direction_joint(
              D,
              up,
              down,
              metric_scale(metric),
              meta
            )

            statuses[[kk]] <-
              z$status

            effects[[kk]] <-
              z$effects
          }
        }
      }

      list(
        status =
          rbindlist(
            statuses,
            fill = TRUE
          ),
        effects =
          rbindlist(
            effects,
            fill = TRUE
          )
      )
    }
  )

  JOINT <- checkpoint_run(
    radius,
    "04_count_heat_joint",
    function() {

      statuses <- list()
      effects <- list()
      collinearity <- list()

      kk <- 0L

      for (branch in BRANCHES) {

        D <- if (
          branch ==
            "quasi_daily"
        ) {
          Q_RAW
        } else {
          FP_RAW
        }

        for (fuel in FUELS) {

          kk <- kk +
            1L

          cc <- raw_col(
            fuel,
            "upwind",
            "count",
            radius
          )

          hh <- raw_col(
            fuel,
            "upwind",
            "heat",
            radius
          )

          if (
            branch ==
              "filter_period"
          ) {

            cc <- paste0(
              "fpmean_",
              cc
            )

            hh <- paste0(
              "fpmean_",
              hh
            )
          }

          meta <- data.table(
            radius_km =
              radius,
            branch =
              branch,
            fuel =
              fuel,
            fuel_label =
              unname(
                FUEL_CN[fuel]
              )
          )

          cat(
            "[Count+Heat ",
            radius,
            "km ",
            kk,
            "/12] ",
            branch,
            " | ",
            fuel,
            "\n",
            sep = ""
          )

          z <- fit_count_heat_joint(
            D,
            cc,
            hh,
            meta
          )

          statuses[[kk]] <-
            z$status

          effects[[kk]] <-
            z$effects

          collinearity[[kk]] <-
            z$collinearity
        }
      }

      list(
        status =
          rbindlist(
            statuses,
            fill = TRUE
          ),
        effects =
          rbindlist(
            effects,
            fill = TRUE
          ),
        collinearity =
          rbindlist(
            collinearity,
            fill = TRUE
          )
      )
    }
  )

  PM_SENS <- checkpoint_run(
    radius,
    "05_PM_sensitivity",
    function() {

      statuses <- list()
      effects <- list()

      kk <- 0L

      for (branch in BRANCHES) {

        D <- if (
          branch ==
            "quasi_daily"
        ) {
          Q_RAW
        } else {
          FP_RAW
        }

        for (fuel in FUELS) {

          for (
            metric in
            c(
              "count",
              "heat"
            )
          ) {

            kk <- kk +
              1L

            col <- raw_col(
              fuel,
              "upwind",
              metric,
              radius
            )

            if (
              branch ==
                "filter_period"
            ) {
              col <- paste0(
                "fpmean_",
                col
              )
            }

            meta <- data.table(
              radius_km =
                radius,
              branch =
                branch,
              fuel =
                fuel,
              fuel_label =
                unname(
                  FUEL_CN[fuel]
                ),
              metric =
                metric,
              metric_label =
                metric_cn(metric),
              unit =
                metric_unit(metric),
              pm_variant =
                "不调整PM2.5_但保持与主模型完全相同行"
            )

            cat(
              "[PM敏感性 ",
              radius,
              "km ",
              kk,
              "/24] ",
              branch,
              " | ",
              fuel,
              " | ",
              metric,
              "\n",
              sep = ""
            )

            z <- fit_single_fe(
              D,
              col,
              metric_scale(metric),
              meta,
              strict_time =
                FALSE,
              include_pm =
                FALSE,
              make_curve =
                FALSE,
              make_distribution_summary =
                FALSE
            )

            statuses[[kk]] <-
              z$status

            effects[[kk]] <-
              z$effects
          }
        }
      }

      list(
        status =
          rbindlist(
            statuses,
            fill = TRUE
          ),
        effects =
          rbindlist(
            effects,
            fill = TRUE
          )
      )
    }
  )

  DECOMP <- checkpoint_run(
    radius,
    "06_within_between",
    function() {

      statuses <- list()
      effects <- list()

      kk <- 0L

      for (branch in BRANCHES) {

        D <- if (
          branch ==
            "quasi_daily"
        ) {
          Q_RAW
        } else {
          FP_RAW
        }

        for (fuel in FUELS) {

          for (
            metric in
            c(
              "count",
              "heat"
            )
          ) {

            kk <- kk +
              1L

            col <- raw_col(
              fuel,
              "upwind",
              metric,
              radius
            )

            if (
              branch ==
                "filter_period"
            ) {

              col <- paste0(
                "fpmean_",
                col
              )
            }

            meta <- data.table(
              radius_km =
                radius,
              branch =
                branch,
              fuel =
                fuel,
              fuel_label =
                unname(
                  FUEL_CN[fuel]
                ),
              metric =
                metric,
              metric_label =
                metric_cn(metric),
              unit =
                metric_unit(metric)
            )

            cat(
              "[内外分解 ",
              radius,
              "km ",
              kk,
              "/24] ",
              branch,
              " | ",
              fuel,
              " | ",
              metric,
              "\n",
              sep = ""
            )

            z <- fit_mundlak(
              D,
              col,
              metric_scale(metric),
              meta
            )

            statuses[[kk]] <-
              z$status

            effects[[kk]] <-
              z$effects
          }
        }
      }

      list(
        status =
          rbindlist(
            statuses,
            fill = TRUE
          ),
        effects =
          rbindlist(
            effects,
            fill = TRUE
          )
      )
    }
  )

  WEIGHTED_DAILY <- read_weighted_radius(
    radius
  )

  WMAP <- WEIGHTED_MAP[
    radius_km ==
      radius
  ]

  WEIGHTED_COLS <- unique(
    c(
      WMAP$proximity_col,
      WMAP$heat_proximity_col
    )
  )

  Q_WEIGHTED <- merge_quasi(
    WEIGHTED_DAILY,
    WEIGHTED_COLS
  )

  FP_WEIGHTED <- aggregate_to_filter(
    WEIGHTED_DAILY,
    WEIGHTED_COLS
  )

  WEIGHTED <- checkpoint_run(
    radius,
    "07_weighted",
    function() {

      statuses <- list()
      effects <- list()

      kk <- 0L

      for (branch in BRANCHES) {

        D <- if (
          branch ==
            "quasi_daily"
        ) {
          Q_WEIGHTED
        } else {
          FP_WEIGHTED
        }

        for (fuel in FUELS) {

          for (
            direction in
            DIRECTIONS
          ) {

            fuel_value <- fuel
            direction_value <- direction

            rowmap <- WMAP[
              fuel == fuel_value &
              direction == direction_value
            ]

            if (nrow(rowmap) != 1L) {
              stopf(
                "07H04距离加权字段映射不是唯一一行：radius=",
                radius,
                "; fuel=",
                fuel_value,
                "; direction=",
                direction_value,
                "; nrow=",
                nrow(rowmap)
              )
            }

            for (
              metric in
              c(
                "proximity",
                "heat_proximity"
              )
            ) {

              kk <- kk +
                1L

              col <- if (
                metric ==
                  "proximity"
              ) {
                rowmap$proximity_col[1]
              } else {
                rowmap$heat_proximity_col[1]
              }

              if (
                branch ==
                  "filter_period"
              ) {

                col <- paste0(
                  "fpmean_",
                  col
                )
              }

              meta <- data.table(
                radius_km =
                  radius,
                branch =
                  branch,
                fuel =
                  fuel,
                fuel_label =
                  unname(
                    FUEL_CN[fuel]
                  ),
                direction =
                  direction,
                metric =
                  metric,
                metric_label =
                  metric_cn(metric),
                unit =
                  metric_unit(metric),
                model_family =
                  "距离加权敏感性"
              )

              cat(
                "[距离加权 ",
                radius,
                "km ",
                kk,
                "/48] ",
                branch,
                " | ",
                fuel,
                " | ",
                direction,
                " | ",
                metric,
                "\n",
                sep = ""
              )

              z <- fit_single_fe(
                D,
                col,
                metric_scale(metric),
                meta,
                strict_time =
                  FALSE,
                include_pm =
                  TRUE,
                make_curve =
                  FALSE,
                make_distribution_summary =
                  FALSE
              )

              statuses[[kk]] <-
                z$status

              effects[[kk]] <-
                z$effects
            }
          }
        }
      }

      list(
        status =
          rbindlist(
            statuses,
            fill = TRUE
          ),
        effects =
          rbindlist(
            effects,
            fill = TRUE
          )
      )
    }
  )

  LAG <- checkpoint_run(
    radius,
    "08_lag",
    function() {

      L <- make_lag_daily(
        RAW[
          ,
          c(
            "loc_num",
            "date",
            unlist(
              lapply(
                FUELS,
                function(fuel) {

                  c(
                    raw_col(
                      fuel,
                      "upwind",
                      "count",
                      radius
                    ),
                    raw_col(
                      fuel,
                      "upwind",
                      "heat",
                      radius
                    )
                  )
                }
              )
            )
          ),
          with = FALSE
        ],
        radius
      )

      LAG_DAILY <- L$data

      LAG_MAP <- L$map

      lag_cols_needed <- unique(
        c(
          LAG_MAP$exposure_col,
          LAG_MAP$wind_col
        )
      )

      Q_LAG <- merge(
        Q_BASE,
        LAG_DAILY[
          ,
          c(
            "loc_num",
            "date",
            lag_cols_needed
          ),
          with = FALSE
        ],
        by = c(
          "loc_num",
          "date"
        ),
        all.x = TRUE,
        sort = FALSE,
        allow.cartesian = TRUE
      )

      statuses <- list()
      effects <- list()

      kk <- 0L

      for (fuel in FUELS) {

        for (
          metric in
          c(
            "count",
            "heat"
          )
        ) {

          original <- raw_col(
            fuel,
            "upwind",
            metric,
            radius
          )

          LM <- LAG_MAP[
            original_col ==
              original
          ]

          for (
            ii in
            seq_len(
              nrow(LM)
            )
          ) {

            kk <- kk +
              1L

            meta <- data.table(
              radius_km =
                radius,
              branch =
                "quasi_daily",
              fuel =
                fuel,
              fuel_label =
                unname(
                  FUEL_CN[fuel]
                ),
              metric =
                metric,
              metric_label =
                metric_cn(metric),
              unit =
                metric_unit(metric),
              lag_spec =
                LM$lag_spec[ii]
            )

            cat(
              "[滞后 ",
              radius,
              "km ",
              kk,
              "/72] ",
              fuel,
              " | ",
              metric,
              " | ",
              LM$lag_spec[ii],
              "\n",
              sep = ""
            )

            z <- fit_lag_one(
              Q_LAG,
              LM$exposure_col[ii],
              LM$wind_col[ii],
              metric_scale(metric),
              meta
            )

            statuses[[kk]] <-
              z$status

            effects[[kk]] <-
              z$effects
          }
        }
      }

      list(
        status =
          rbindlist(
            statuses,
            fill = TRUE
          ),
        effects =
          rbindlist(
            effects,
            fill = TRUE
          )
      )
    }
  )

  rm(
    RAW,
    Q_RAW,
    FP_RAW,
    WEIGHTED_DAILY,
    Q_WEIGHTED,
    FP_WEIGHTED
  )

  invisible(
    gc(
      verbose = FALSE
    )
  )

  list(
    main =
      MAIN,
    strict =
      STRICT,
    direction =
      DIRECTION,
    joint =
      JOINT,
    pm =
      PM_SENS,
    decomposition =
      DECOMP,
    weighted =
      WEIGHTED,
    lag =
      LAG
  )
}

ALL_RADIUS_RESULTS <- vector(
  "list",
  length(RADII)
)

names(
  ALL_RADIUS_RESULTS
) <- as.character(
  RADII
)

for (
  rr in
  RADII
) {

  ALL_RADIUS_RESULTS[[
    as.character(
      rr
    )
  ]] <- process_radius(
    rr
  )
}

collect_table <- function(
  family,
  table_name
) {

  rbindlist(
    lapply(
      ALL_RADIUS_RESULTS,
      function(x) {

        z <- x[[family]][[
          table_name
        ]]

        if (
          is.null(z)
        ) {
          return(
            data.table()
          )
        }

        z
      }
    ),
    use.names = TRUE,
    fill = TRUE
  )
}

MAIN_STATUS <- collect_table(
  "main",
  "status"
)

MAIN_EFFECTS <- collect_table(
  "main",
  "effects"
)

MAIN_CURVE <- collect_table(
  "main",
  "curve"
)

MAIN_WIND <- collect_table(
  "main",
  "wind_summary"
)

STRICT_STATUS <- collect_table(
  "strict",
  "status"
)

STRICT_EFFECTS <- collect_table(
  "strict",
  "effects"
)

STRICT_CURVE <- collect_table(
  "strict",
  "curve"
)

STRICT_WIND <- collect_table(
  "strict",
  "wind_summary"
)

DIRECTION_STATUS <- collect_table(
  "direction",
  "status"
)

DIRECTION_EFFECTS <- collect_table(
  "direction",
  "effects"
)

JOINT_STATUS <- collect_table(
  "joint",
  "status"
)

JOINT_EFFECTS <- collect_table(
  "joint",
  "effects"
)

JOINT_COLL <- collect_table(
  "joint",
  "collinearity"
)

PM_STATUS <- collect_table(
  "pm",
  "status"
)

PM_EFFECTS <- collect_table(
  "pm",
  "effects"
)

DECOMP_STATUS <- collect_table(
  "decomposition",
  "status"
)

DECOMP_EFFECTS <- collect_table(
  "decomposition",
  "effects"
)

WEIGHTED_STATUS <- collect_table(
  "weighted",
  "status"
)

WEIGHTED_EFFECTS <- collect_table(
  "weighted",
  "effects"
)

LAG_STATUS <- collect_table(
  "lag",
  "status"
)

LAG_EFFECTS <- collect_table(
  "lag",
  "effects"
)

fwrite(
  MAIN_STATUS,
  file.path(
    DIR_MAIN,
    "01_全部半径主模型状态.csv"
  )
)

fwrite(
  MAIN_EFFECTS,
  file.path(
    DIR_MAIN,
    "02_全部半径主模型代表风速效应.csv"
  )
)

fwrite(
  MAIN_CURVE,
  file.path(
    DIR_WIND,
    "01_主模型连续风速效应曲线.csv"
  )
)

fwrite(
  MAIN_WIND,
  file.path(
    DIR_WIND,
    "02_主模型实际风速分布效应审计.csv"
  )
)

fwrite(
  STRICT_STATUS,
  file.path(
    DIR_STRICT,
    "01_全部半径严格时间固定效应模型状态.csv"
  )
)

fwrite(
  STRICT_EFFECTS,
  file.path(
    DIR_STRICT,
    "02_全部半径严格时间固定效应结果.csv"
  )
)

fwrite(
  STRICT_CURVE,
  file.path(
    DIR_WIND,
    "03_严格时间固定效应连续风速曲线.csv"
  )
)

fwrite(
  STRICT_WIND,
  file.path(
    DIR_WIND,
    "04_严格时间固定效应实际风速分布审计.csv"
  )
)

fwrite(
  DIRECTION_STATUS,
  file.path(
    DIR_DIRECTION,
    "01_上风向下风向联合模型状态.csv"
  )
)

fwrite(
  DIRECTION_EFFECTS,
  file.path(
    DIR_DIRECTION,
    "02_上风向下风向及直接差值检验.csv"
  )
)

fwrite(
  WEIGHTED_STATUS,
  file.path(
    DIR_WEIGHTED,
    "01_距离加权固定效应模型状态.csv"
  )
)

fwrite(
  WEIGHTED_EFFECTS,
  file.path(
    DIR_WEIGHTED,
    "02_距离加权全部结果.csv"
  )
)

fwrite(
  JOINT_STATUS,
  file.path(
    DIR_JOINT,
    "01_Count_Heat联合模型状态.csv"
  )
)

fwrite(
  JOINT_EFFECTS,
  file.path(
    DIR_JOINT,
    "02_Count_Heat联合模型效应.csv"
  )
)

fwrite(
  JOINT_COLL,
  file.path(
    DIR_JOINT,
    "03_Count_Heat站点年份内相关性与VIF.csv"
  )
)

fwrite(
  PM_STATUS,
  file.path(
    DIR_PM,
    "01_不调整PM2.5模型状态.csv"
  )
)

fwrite(
  PM_EFFECTS,
  file.path(
    DIR_PM,
    "02_不调整PM2.5同样本结果.csv"
  )
)

fwrite(
  DECOMP_STATUS,
  file.path(
    DIR_DECOMP,
    "01_站点年份内外分解模型状态.csv"
  )
)

fwrite(
  DECOMP_EFFECTS,
  file.path(
    DIR_DECOMP,
    "02_站点年份内外分解全部结果.csv"
  )
)

fwrite(
  LAG_STATUS,
  file.path(
    DIR_LAG,
    "01_准逐日滞后模型状态.csv"
  )
)

fwrite(
  LAG_EFFECTS,
  file.path(
    DIR_LAG,
    "02_准逐日滞后全部结果.csv"
  )
)

MAIN_MEDIAN <- MAIN_EFFECTS[
  wind_reference ==
    "样本中位数"
]

STRICT_MEDIAN <- STRICT_EFFECTS[
  wind_reference ==
    "样本中位数"
]

PM_MEDIAN <- PM_EFFECTS[
  wind_reference ==
    "样本中位数"
]

WEIGHTED_UP_MEDIAN <- WEIGHTED_EFFECTS[
  direction ==
    "upwind" &
  wind_reference ==
    "样本中位数"
]

DIRECTION_MEDIAN <- DIRECTION_EFFECTS[
  component ==
    "上风向减下风向" &
  wind_reference ==
    "样本中位数"
]

DECOMP_WITHIN_1 <- DECOMP_EFFECTS[
  component ==
    "站点年份内变化" &
  wind_reference ==
    "参考论文1米每秒"
]

DECOMP_BETWEEN_1 <- DECOMP_EFFECTS[
  component ==
    "站点年份间平均差异" &
  wind_reference ==
    "参考论文1米每秒"
]

OLD <- fread(
  file =
    H07_MAIN,
  showProgress =
    FALSE
)

OLD_KEY <- OLD[
  component %in%
    c(
      "raw_count_only",
      "raw_heat_only"
    ),
  .(
    branch,
    fuel,
    radius_km,
    metric =
      fifelse(
        component ==
          "raw_count_only",
        "count",
        "heat"
      ),
    old_status =
      status,
    old_estimate_wind1 =
      estimate_fixed,
    old_ci_low_wind1 =
      ci_low_fixed,
    old_ci_high_wind1 =
      ci_high_fixed
  )
]

WITHIN_BETWEEN_COMPARE <- merge(
  OLD_KEY,
  DECOMP_WITHIN_1[
    ,
    .(
      branch,
      fuel,
      radius_km,
      metric,
      within_estimate_wind1 =
        estimate,
      within_ci_low_wind1 =
        ci95_low,
      within_ci_high_wind1 =
        ci95_high
    )
  ],
  by = c(
    "branch",
    "fuel",
    "radius_km",
    "metric"
  ),
  all = TRUE,
  sort = FALSE
)

WITHIN_BETWEEN_COMPARE <- merge(
  WITHIN_BETWEEN_COMPARE,
  DECOMP_BETWEEN_1[
    ,
    .(
      branch,
      fuel,
      radius_km,
      metric,
      between_estimate_wind1 =
        estimate,
      between_ci_low_wind1 =
        ci95_low,
      between_ci_high_wind1 =
        ci95_high
    )
  ],
  by = c(
    "branch",
    "fuel",
    "radius_km",
    "metric"
  ),
  all = TRUE,
  sort = FALSE
)

WITHIN_BETWEEN_COMPARE[
  ,
  sign_pattern :=
    fifelse(
      is.finite(
        old_estimate_wind1
      ) &
      is.finite(
        within_estimate_wind1
      ) &
      is.finite(
        between_estimate_wind1
      ),
      paste0(
        "原模型",
        fifelse(
          old_estimate_wind1 >
            0,
          "正",
          "负"
        ),
        "_站内",
        fifelse(
          within_estimate_wind1 >
            0,
          "正",
          "负"
        ),
        "_站点年份间",
        fifelse(
          between_estimate_wind1 >
            0,
          "正",
          "负"
        )
      ),
      "无法完整比较"
    )
]

fwrite(
  WITHIN_BETWEEN_COMPARE,
  file.path(
    DIR_DECOMP,
    "03_冻结原模型与站点年份内外方向解释.csv"
  )
)

SUMMARY <- MAIN_MEDIAN[
  ,
  .(
    branch,
    fuel,
    fuel_label,
    radius_km,
    metric,
    metric_label,
    unit,
    main_estimate =
      estimate,
    main_ci_low =
      ci95_low,
    main_ci_high =
      ci95_high,
    main_direction =
      direction,
    main_ci_status =
      ci_status
  )
]

SUMMARY <- merge(
  SUMMARY,
  STRICT_MEDIAN[
    ,
    .(
      branch,
      fuel,
      radius_km,
      metric,
      strict_estimate =
        estimate,
      strict_ci_low =
        ci95_low,
      strict_ci_high =
        ci95_high,
      strict_direction =
        direction,
      strict_ci_status =
        ci_status
    )
  ],
  by = c(
    "branch",
    "fuel",
    "radius_km",
    "metric"
  ),
  all.x = TRUE,
  sort = FALSE
)

SUMMARY <- merge(
  SUMMARY,
  PM_MEDIAN[
    ,
    .(
      branch,
      fuel,
      radius_km,
      metric,
      no_pm_estimate =
        estimate,
      no_pm_ci_low =
        ci95_low,
      no_pm_ci_high =
        ci95_high,
      no_pm_direction =
        direction
    )
  ],
  by = c(
    "branch",
    "fuel",
    "radius_km",
    "metric"
  ),
  all.x = TRUE,
  sort = FALSE
)

SUMMARY <- merge(
  SUMMARY,
  DIRECTION_MEDIAN[
    ,
    .(
      branch,
      fuel,
      radius_km,
      metric,
      up_minus_down_estimate =
        estimate,
      up_minus_down_ci_low =
        ci95_low,
      up_minus_down_ci_high =
        ci95_high,
      up_minus_down_direction =
        direction,
      up_minus_down_ci_status =
        ci_status
    )
  ],
  by = c(
    "branch",
    "fuel",
    "radius_km",
    "metric"
  ),
  all.x = TRUE,
  sort = FALSE
)

WEIGHTED_MATCH <- copy(
  WEIGHTED_UP_MEDIAN
)

WEIGHTED_MATCH[
  ,
  raw_metric :=
    fifelse(
      metric ==
        "proximity",
      "count",
      fifelse(
        metric ==
          "heat_proximity",
        "heat",
        NA_character_
      )
    )
]

SUMMARY <- merge(
  SUMMARY,
  WEIGHTED_MATCH[
    !is.na(
      raw_metric
    ),
    .(
      branch,
      fuel,
      radius_km,
      metric =
        raw_metric,
      weighted_estimate =
        estimate,
      weighted_ci_low =
        ci95_low,
      weighted_ci_high =
        ci95_high,
      weighted_direction =
        direction,
      weighted_ci_status =
        ci_status
    )
  ],
  by = c(
    "branch",
    "fuel",
    "radius_km",
    "metric"
  ),
  all.x = TRUE,
  sort = FALSE
)

WIND_PATTERN <- MAIN_CURVE[
  ,
  .(
    wind_pattern =
      fifelse(
        all(
          estimate >
            0,
          na.rm = TRUE
        ),
        "实际主要风速范围内始终正向",
        fifelse(
          all(
            estimate <
              0,
            na.rm = TRUE
          ),
          "实际主要风速范围内始终负向",
          "实际主要风速范围内发生方向转换"
        )
      ),
    curve_min_estimate =
      min(
        estimate,
        na.rm = TRUE
      ),
    curve_max_estimate =
      max(
        estimate,
        na.rm = TRUE
      )
  ),
  by = .(
    branch,
    fuel,
    radius_km,
    metric
  )
]

SUMMARY <- merge(
  SUMMARY,
  WIND_PATTERN,
  by = c(
    "branch",
    "fuel",
    "radius_km",
    "metric"
  ),
  all.x = TRUE,
  sort = FALSE
)

SUMMARY[
  ,
  main_vs_strict :=
    fifelse(
      is.finite(
        main_estimate
      ) &
      is.finite(
        strict_estimate
      ),
      fifelse(
        sign(
          main_estimate
        ) ==
        sign(
          strict_estimate
        ),
        fifelse(
          main_estimate >
            0,
          "主模型与严格时间控制均正向",
          "主模型与严格时间控制均负向"
        ),
        "严格时间控制后方向改变"
      ),
      "无法比较"
    )
]

SUMMARY[
  ,
  pm_sensitivity :=
    fifelse(
      is.finite(
        main_estimate
      ) &
      is.finite(
        no_pm_estimate
      ),
      fifelse(
        sign(
          main_estimate
        ) ==
        sign(
          no_pm_estimate
        ),
        "调整与不调整PM2.5方向一致",
        "PM2.5调整后方向改变"
      ),
      "无法比较"
    )
]

SUMMARY[
  ,
  weighted_sensitivity :=
    fifelse(
      is.finite(
        main_estimate
      ) &
      is.finite(
        weighted_estimate
      ),
      fifelse(
        sign(
          main_estimate
        ) ==
        sign(
          weighted_estimate
        ),
        "普通累计暴露与距离加权方向一致",
        "距离加权后方向改变"
      ),
      "无法比较"
    )
]

SUMMARY[
  ,
  directionality_interpretation :=
    fifelse(
      is.finite(
        up_minus_down_ci_low
      ) &
      up_minus_down_ci_low >
        0,
      "上风向效应明确高于下风向",
      fifelse(
        is.finite(
          up_minus_down_ci_high
        ) &
        up_minus_down_ci_high <
          0,
        "下风向效应明确高于上风向",
        "上风向与下风向差值置信区间包含0"
      )
    )
]

setorder(
  SUMMARY,
  branch,
  fuel,
  metric,
  radius_km
)

fwrite(
  SUMMARY,
  file.path(
    DIR_SUMMARY,
    "01_全部半径综合结果主表.csv"
  )
)

LAG_MEDIAN <- LAG_EFFECTS[
  wind_reference ==
    "样本中位数"
]

LAG_SUMMARY <- LAG_MEDIAN[
  ,
  .(
    n_lag_specs =
      .N,
    n_positive =
      sum(
        estimate >
          0,
        na.rm = TRUE
      ),
    n_negative =
      sum(
        estimate <
          0,
        na.rm = TRUE
      ),
    n_CI_above_zero =
      sum(
        ci95_low >
          0,
        na.rm = TRUE
      ),
    n_CI_below_zero =
      sum(
        ci95_high <
          0,
        na.rm = TRUE
      ),
    lag_direction_sequence =
      paste(
        paste0(
          lag_spec,
          ":",
          direction
        ),
        collapse = " | "
      )
  ),
  by = .(
    fuel,
    fuel_label,
    radius_km,
    metric,
    metric_label
  )
]

fwrite(
  LAG_SUMMARY,
  file.path(
    DIR_LAG,
    "03_准逐日各滞后方向汇总_不按显著性选最佳滞后.csv"
  )
)

STATUS_SUMMARY <- rbindlist(
  list(
    MAIN_STATUS[
      ,
      .(
        family =
          "主分析",
        status
      )
    ],
    STRICT_STATUS[
      ,
      .(
        family =
          "严格时间固定效应",
        status
      )
    ],
    DIRECTION_STATUS[
      ,
      .(
        family =
          "上风向下风向联合模型",
        status
      )
    ],
    WEIGHTED_STATUS[
      ,
      .(
        family =
          "距离加权",
        status
      )
    ],
    JOINT_STATUS[
      ,
      .(
        family =
          "Count_Heat联合",
        status
      )
    ],
    PM_STATUS[
      ,
      .(
        family =
          "PM敏感性",
        status
      )
    ],
    DECOMP_STATUS[
      ,
      .(
        family =
          "站点年份内外分解",
        status
      )
    ],
    LAG_STATUS[
      ,
      .(
        family =
          "准逐日滞后",
        status
      )
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)

STATUS_COUNTS <- STATUS_SUMMARY[
  ,
  .N,
  by = .(
    family,
    status
  )
]

fwrite(
  STATUS_COUNTS,
  file.path(
    DIR_SUMMARY,
    "02_全部模型状态汇总.csv"
  )
)

EXPECTED <- data.table(
  family = c(
    "主分析",
    "严格时间固定效应",
    "上风向下风向联合模型",
    "距离加权",
    "Count_Heat联合",
    "PM敏感性",
    "站点年份内外分解",
    "准逐日滞后"
  ),
  expected_specs = c(
    240L,
    240L,
    240L,
    480L,
    120L,
    240L,
    240L,
    720L
  )
)

OBSERVED <- data.table(
  family = EXPECTED$family,
  observed_specs = c(
    nrow(
      MAIN_STATUS
    ),
    nrow(
      STRICT_STATUS
    ),
    nrow(
      DIRECTION_STATUS
    ),
    nrow(
      WEIGHTED_STATUS
    ),
    nrow(
      JOINT_STATUS
    ),
    nrow(
      PM_STATUS
    ),
    nrow(
      DECOMP_STATUS
    ),
    nrow(
      LAG_STATUS
    )
  ),
  n_FAIL = c(
    MAIN_STATUS[
      status ==
        "FAIL",
      .N
    ],
    STRICT_STATUS[
      status ==
        "FAIL",
      .N
    ],
    DIRECTION_STATUS[
      status ==
        "FAIL",
      .N
    ],
    WEIGHTED_STATUS[
      status ==
        "FAIL",
      .N
    ],
    JOINT_STATUS[
      status ==
        "FAIL",
      .N
    ],
    PM_STATUS[
      status ==
        "FAIL",
      .N
    ],
    DECOMP_STATUS[
      status ==
        "FAIL",
      .N
    ],
    LAG_STATUS[
      status ==
        "FAIL",
      .N
    ]
  )
)

FAMILY_AUDIT <- merge(
  EXPECTED,
  OBSERVED,
  by =
    "family",
  all =
    TRUE
)

FAMILY_AUDIT[
  ,
  status :=
    fifelse(
      expected_specs ==
        observed_specs &
      n_FAIL ==
        0L,
      "PASS",
      "FAIL"
    )
]

fwrite(
  FAMILY_AUDIT,
  file.path(
    DIR_SUMMARY,
    "03_分析家族完整性审计.csv"
  )
)

N_CHECKPOINTS <- length(
  list.files(
    CHECK_DIR,
    pattern =
      "^radius_.*\\.rds$"
  )
)

HARD <- data.table(
  check = c(
    "07H09B_upstream_all_PASS",
    "07H01_daily_files_24",
    "quasi_primary_rows_426657",
    "quasi_primary_monitors_133",
    "filter_all_rows_174200",
    "filter_primary_rows_107407",
    "filter_primary_monitors_133",
    "radii_10",
    "main_specs_240",
    "strict_specs_240",
    "direction_specs_240",
    "weighted_specs_480",
    "joint_specs_120",
    "PM_specs_240",
    "decomposition_specs_240",
    "lag_specs_720",
    "all_model_families_FAIL_zero",
    "main_continuous_wind_curve_exists",
    "strict_continuous_wind_curve_exists",
    "integrated_summary_nonempty"
  ),
  observed = c(
    as.character(
      H09B_AUDIT[
        status !=
          "PASS",
        .N
      ]
    ),
    as.character(
      sum(
        file.exists(
          RAW_FILES
        )
      )
    ),
    as.character(
      N_Q
    ),
    as.character(
      length(
        ELIGIBLE_LOCS
      )
    ),
    as.character(
      N_FP_ALL
    ),
    as.character(
      nrow(
        FP_BASE
      )
    ),
    as.character(
      uniqueN(
        FP_BASE$loc_num
      )
    ),
    as.character(
      length(
        RADII
      )
    ),
    as.character(
      nrow(
        MAIN_STATUS
      )
    ),
    as.character(
      nrow(
        STRICT_STATUS
      )
    ),
    as.character(
      nrow(
        DIRECTION_STATUS
      )
    ),
    as.character(
      nrow(
        WEIGHTED_STATUS
      )
    ),
    as.character(
      nrow(
        JOINT_STATUS
      )
    ),
    as.character(
      nrow(
        PM_STATUS
      )
    ),
    as.character(
      nrow(
        DECOMP_STATUS
      )
    ),
    as.character(
      nrow(
        LAG_STATUS
      )
    ),
    as.character(
      sum(
        FAMILY_AUDIT$n_FAIL
      )
    ),
    as.character(
      as.integer(
        nrow(
          MAIN_CURVE
        ) >
          0L
      )
    ),
    as.character(
      as.integer(
        nrow(
          STRICT_CURVE
        ) >
          0L
      )
    ),
    as.character(
      as.integer(
        nrow(
          SUMMARY
        ) >
          0L
      )
    )
  ),
  expected = c(
    "0",
    "24",
    "426657",
    "133",
    "174200",
    "107407",
    "133",
    "10",
    "240",
    "240",
    "240",
    "480",
    "120",
    "240",
    "240",
    "720",
    "0",
    "1",
    "1",
    "1"
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
  file.path(
    DIR_SUMMARY,
    "04_07H10_hard_checks.csv"
  )
)

MAIN_MEDIAN_SUM <- MAIN_MEDIAN[
  ,
  .(
    n_PASS_effects =
      .N,
    n_positive =
      sum(
        estimate >
          0,
        na.rm = TRUE
      ),
    n_negative =
      sum(
        estimate <
          0,
        na.rm = TRUE
      ),
    n_CI_above_zero =
      sum(
        ci95_low >
          0,
        na.rm = TRUE
      ),
    n_CI_below_zero =
      sum(
        ci95_high <
          0,
        na.rm = TRUE
      ),
    n_CI_include_zero =
      sum(
        ci95_low <=
          0 &
        ci95_high >=
          0,
        na.rm = TRUE
      )
  ),
  by = .(
    branch,
    fuel,
    fuel_label,
    metric,
    metric_label
  )
]

fwrite(
  MAIN_MEDIAN_SUM,
  file.path(
    DIR_SUMMARY,
    "05_主分析全部半径方向汇总.csv"
  )
)

RADIUS_PATTERN <- SUMMARY[
  ,
  .(
    n_radii =
      .N,
    n_main_positive =
      sum(
        main_estimate >
          0,
        na.rm = TRUE
      ),
    n_main_negative =
      sum(
        main_estimate <
          0,
        na.rm = TRUE
      ),
    n_strict_positive =
      sum(
        strict_estimate >
          0,
        na.rm = TRUE
      ),
    n_strict_negative =
      sum(
        strict_estimate <
          0,
        na.rm = TRUE
      ),
    n_main_strict_same_direction =
      sum(
        sign(
          main_estimate
        ) ==
        sign(
          strict_estimate
        ),
        na.rm = TRUE
      ),
    n_directionality_CI_above_zero =
      sum(
        up_minus_down_ci_low >
          0,
        na.rm = TRUE
      ),
    n_weighted_same_direction =
      sum(
        sign(
          main_estimate
        ) ==
        sign(
          weighted_estimate
        ),
        na.rm = TRUE
      ),
    n_PM_same_direction =
      sum(
        sign(
          main_estimate
        ) ==
        sign(
          no_pm_estimate
        ),
        na.rm = TRUE
      ),
    wind_patterns =
      paste(
        unique(
          wind_pattern
        ),
        collapse = " | "
      )
  ),
  by = .(
    branch,
    fuel,
    fuel_label,
    metric,
    metric_label
  )
]

fwrite(
  RADIUS_PATTERN,
  file.path(
    DIR_SUMMARY,
    "06_不按显著性选择半径的整体模式汇总.csv"
  )
)

cat(
  "\n============================================================\n",
  "A. 07H10全部模型状态\n",
  "============================================================\n",
  sep = ""
)

print(
  STATUS_COUNTS,
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "B. 主分析：样本中位风速下全部半径方向汇总\n",
  "============================================================\n",
  sep = ""
)

print(
  MAIN_MEDIAN_SUM,
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "C. 不按显著性选择半径的10半径综合模式\n",
  "============================================================\n",
  sep = ""
)

print(
  RADIUS_PATTERN,
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "D. 分析家族完整性审计\n",
  "============================================================\n",
  sep = ""
)

print(
  FAMILY_AUDIT,
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "E. HARD CHECKS\n",
  "============================================================\n",
  sep = ""
)

print(
  HARD,
  nrows = Inf
)

HARD_FAIL <- HARD[
  status ==
    "FAIL",
  .N
]

REPORT <- c(
  "# 07H10 燃料特异站点年份内综合正式分析",
  "",
  paste0(
    "- 生成时间：",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  "",
  "## 核心研究问题",
  "",
  "- 使用2001–2024全部数据。",
  "- 主要识别来源是同一个RadNet监测站、同一年内部随时间发生的暴露变化。",
  "- 每一个monitor×year拥有独立固定基线。",
  "- 因此不同城市、不同监测站以及同一监测站跨年份的长期基础水平不直接用于识别主要暴露系数。",
  "",
  "## 正式暴露",
  "",
  "- 六种燃烧燃料：coal、gas、oil、other_fossil、biomass、solid_waste。",
  "- 累计半径：20、25、30、35、40、45、50、100、150、200 km。",
  "- 方向：90°上风向和180°相反方向的下风向对照。",
  "- 普通暴露：活动电厂数量与实际CAMPD Heat Input。",
  "- 距离加权敏感性：Σ1/d与ΣHeat/d。",
  "- 距离加权仍是空间权重，不解释为真实大气浓度或放射性剂量。",
  "",
  "## 正式模型层级",
  "",
  "- 主分析：monitor×year固定效应，标准误按监测站聚类。",
  "- 严格时间敏感性：monitor×year + 共同时间固定效应。",
  "- quasi-daily共同时间固定效应为具体日期。",
  "- filter-period共同时间固定效应为collect_start + collect_end完整采样区间。",
  "- 风速始终为连续变量。",
  "- 主暴露与连续风速保留交互。",
  "- 1 m/s仅用于与参考论文比较，不作为唯一科学结论。",
  "- 同时报告实际样本第25百分位、中位数、均值和第75百分位。",
  "- 连续风速曲线覆盖实际样本第1–99百分位。",
  "",
  "## 方向性检验",
  "",
  "- 上风向与下风向暴露在同一个固定效应模型内同时估计。",
  "- 直接计算上风向效应减去下风向效应及95%置信区间。",
  "- 不使用“一个显著、另一个不显著”替代正式差值检验。",
  "",
  "## 补充分析",
  "",
  "- Count + Heat联合模型。",
  "- 站点年份内Count–Heat相关性及近似VIF。",
  "- 不调整PM2.5的完全同样本敏感性。",
  "- 原始累计暴露与Σ1/d、ΣHeat/d距离加权敏感性。",
  "- 站点年份内与站点年份间分解，用于解释旧随机截距模型为什么可能变号。",
  "- quasi-daily lag0、lag1、lag2、lag3、0–1日平均、0–2日平均。",
  "- 滞后分析的风向来自暴露日；暴露×风速交互使用暴露日风速。",
  "",
  "## 半径解释",
  "",
  "- 10个半径全部为预先规定的累计半径。",
  "- 不根据P值或显著性挑选最佳半径。",
  "- 输出整体半径模式，而不是用10个高度相关半径当作10次独立重复。",
  "",
  "## 运行状态",
  "",
  paste0(
    "- checkpoint文件数：",
    N_CHECKPOINTS
  ),
  paste0(
    "- 主分析规格数：",
    nrow(
      MAIN_STATUS
    )
  ),
  paste0(
    "- 严格时间规格数：",
    nrow(
      STRICT_STATUS
    )
  ),
  paste0(
    "- 上/下风向联合规格数：",
    nrow(
      DIRECTION_STATUS
    )
  ),
  paste0(
    "- 距离加权规格数：",
    nrow(
      WEIGHTED_STATUS
    )
  ),
  paste0(
    "- Count+Heat联合规格数：",
    nrow(
      JOINT_STATUS
    )
  ),
  paste0(
    "- PM敏感性规格数：",
    nrow(
      PM_STATUS
    )
  ),
  paste0(
    "- 站点年份内外分解规格数：",
    nrow(
      DECOMP_STATUS
    )
  ),
  paste0(
    "- 准逐日滞后规格数：",
    nrow(
      LAG_STATUS
    )
  ),
  "",
  if (
    HARD_FAIL ==
      0L
  ) {
    "**07H10 COMPREHENSIVE WITHIN-MONITOR-YEAR ANALYSIS: PASS**"
  } else {
    "**07H10 COMPREHENSIVE WITHIN-MONITOR-YEAR ANALYSIS: FAIL**"
  }
)

writeLines(
  REPORT,
  file.path(
    OUT,
    "07H10_综合正式分析报告.md"
  ),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUT,
    "07H10_sessionInfo.txt"
  )
)

cat(
  "\n============================================================\n"
)

if (
  HARD_FAIL ==
    0L
) {

  cat(
    "07H10 COMPREHENSIVE WITHIN-MONITOR-YEAR ANALYSIS: PASS\n"
  )

} else {

  cat(
    "07H10 COMPREHENSIVE WITHIN-MONITOR-YEAR ANALYSIS: FAIL\n"
  )
}

cat(
  "============================================================\n",
  "Hard FAIL = ",
  HARD_FAIL,
  "\n",
  sep = ""
)

if (
  HARD_FAIL >
    0L
) {

  stopf(
    "07H10存在硬性检查失败。禁止冻结或解释最终结果。"
  )
}
