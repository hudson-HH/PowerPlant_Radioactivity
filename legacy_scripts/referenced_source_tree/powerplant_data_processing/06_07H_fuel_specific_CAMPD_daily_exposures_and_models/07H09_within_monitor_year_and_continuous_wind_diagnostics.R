#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999
)

for (pkg in c("data.table", "lme4")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0("缺少R包：", pkg),
      call. = FALSE
    )
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
})

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
    x[
      !is.na(x) &
      nzchar(x)
    ]
  )

  paste(
    x,
    collapse = " | "
  )
}

make_lmer_control <- function() {

  z <- list(
    optimizer = "bobyqa",
    optCtrl = list(
      maxfun = 200000
    ),
    check.conv.singular = "ignore",
    check.rankX = "message+drop.cols"
  )

  fm <- names(
    formals(
      lme4::lmerControl
    )
  )

  if ("autoscale" %in% fm) {
    z$autoscale <- TRUE
  }

  do.call(
    lme4::lmerControl,
    z
  )
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
    return(NA_character_)
  }

  z[1]
}

linear_effect <- function(
  beta,
  V,
  main_term,
  interaction_term,
  wind_ms
) {

  if (
    !main_term %in% names(beta) ||
    is.na(interaction_term) ||
    !interaction_term %in% names(beta)
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

  b1 <- beta[[main_term]]
  b2 <- beta[[interaction_term]]

  est <- b1 + wind_ms * b2

  vv <-
    V[main_term, main_term] +
    wind_ms^2 *
      V[interaction_term, interaction_term] +
    2 *
      wind_ms *
      V[main_term, interaction_term]

  if (!is.finite(vv)) {
    se <- NA_real_
  } else {
    se <- sqrt(
      max(
        vv,
        0
      )
    )
  }

  if (
    is.finite(est) &&
    is.finite(se)
  ) {

    low <- est - 1.96 * se
    high <- est + 1.96 * se

    p <- if (se > 0) {
      2 *
        pnorm(
          abs(est / se),
          lower.tail = FALSE
        )
    } else {
      NA_real_
    }

  } else {

    low <- NA_real_
    high <- NA_real_
    p <- NA_real_
  }

  c(
    estimate = est,
    se = se,
    ci_low = low,
    ci_high = high,
    p = p
  )
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

ci_text <- function(low, high) {

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

strict_mean <- function(x) {

  if (
    length(x) == 0L ||
    any(!is.finite(x))
  ) {
    return(NA_real_)
  }

  mean(x)
}

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

DAILY_DIR <- file.path(
  HROOT,
  "fuel_daily_yearly"
)

H07_DIR <- file.path(
  HROOT,
  "07H07_论文结果表与图形准备"
)

H07_FREEZE <- file.path(
  H07_DIR,
  "07H07_FINAL_FREEZE.md"
)

H07_RESULT <- file.path(
  H07_DIR,
  "01_论文结果表",
  "01_论文主结果_upwind单暴露_全部预设半径.csv"
)

OUT <- file.path(
  HROOT,
  "07H09_站点年份内效应与连续风速诊断"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

for (
  f in c(
    FP_FILE,
    H07_FREEZE,
    H07_RESULT
  )
) {

  if (!file.exists(f)) {
    stopf(
      "缺少必要输入：\n",
      f
    )
  }
}

freeze_text <- readLines(
  H07_FREEZE,
  warn = FALSE,
  encoding = "UTF-8"
)

if (
  !any(
    grepl(
      "07H07 FINAL FREEZE: PASS",
      freeze_text,
      fixed = TRUE
    )
  )
) {
  stopf(
    "07H07没有确认FINAL FREEZE PASS。"
  )
}

YEARS <- 2001:2024

FUELS <- c(
  "coal",
  "gas",
  "oil",
  "other_fossil",
  "biomass",
  "solid_waste"
)

FUEL_CN <- c(
  coal = "煤",
  gas = "天然气",
  oil = "燃油",
  other_fossil = "其他化石燃料",
  biomass = "生物质",
  solid_waste = "固体废物"
)

COUNT_COLS <- paste0(
  "campd_",
  FUELS,
  "_upwind_n_plants_20km_90deg"
)

HEAT_COLS <- paste0(
  "campd_",
  FUELS,
  "_upwind_heat_observed_20km_90deg"
)

FUEL_EXPOSURES <- c(
  COUNT_COLS,
  HEAT_COLS
)

DAILY_FILES <- file.path(
  DAILY_DIR,
  paste0(
    "07H_fuel_daily_exposure_",
    YEARS,
    ".csv.gz"
  )
)

if (!all(file.exists(DAILY_FILES))) {

  print(
    DAILY_FILES[
      !file.exists(DAILY_FILES)
    ]
  )

  stopf(
    "07H01逐日燃料暴露文件不是24/24。"
  )
}

cat(
  "\n读取07H01燃料特异逐日暴露...\n"
)

DAILY <- rbindlist(
  lapply(
    DAILY_FILES,
    function(f) {

      H <- names(
        fread(
          file = f,
          nrows = 0,
          showProgress = FALSE
        )
      )

      req <- c(
        "loc_num",
        "date",
        FUEL_EXPOSURES
      )

      miss <- setdiff(
        req,
        H
      )

      if (length(miss) > 0L) {
        stopf(
          basename(f),
          " 缺少字段：",
          paste(
            miss,
            collapse = ", "
          )
        )
      }

      fread(
        file = f,
        select = req,
        showProgress = FALSE
      )
    }
  ),
  use.names = TRUE
)

DAILY[
  ,
  `:=`(
    loc_num =
      as.character(loc_num),

    date =
      as.IDate(date),

    daily_match =
      1L
  )
]

DAILY_DUP <- DAILY[
  ,
  .N,
  by = .(
    loc_num,
    date
  )
][
  N > 1L,
  .N
]

if (DAILY_DUP > 0L) {
  stopf(
    "07H01逐日燃料暴露存在重复monitor-day。"
  )
}

cat(
  "逐日暴露monitor-days = ",
  format(
    nrow(DAILY),
    big.mark = ","
  ),
  "\n",
  sep = ""
)

Q_BASE <- c(
  "loc_num",
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

Q_LIST <- vector(
  "list",
  length(YEARS)
)

cat(
  "\n读取07F准逐日正式数据...\n"
)

for (ii in seq_along(YEARS)) {

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
    Q_BASE,
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

  Q_LIST[[ii]] <- fread(
    file = f,
    select = Q_BASE,
    showProgress = FALSE
  )
}

Q <- rbindlist(
  Q_LIST,
  use.names = TRUE
)

Q <- Q[
  flag_07f_primary_quasi_sample %in%
    TRUE
]

Q[
  ,
  `:=`(
    loc_num =
      as.character(loc_num),

    date =
      as.IDate(date)
  )
]

N_Q_PREJOIN <- nrow(Q)

ELIGIBLE_LOCS <- sort(
  unique(
    Q$loc_num
  )
)

cat(
  "准逐日正式样本 = ",
  format(
    N_Q_PREJOIN,
    big.mark = ","
  ),
  "\n",
  "准逐日正式监测站 = ",
  length(
    ELIGIBLE_LOCS
  ),
  "\n",
  sep = ""
)

Q <- merge(
  Q,
  DAILY,
  by = c(
    "loc_num",
    "date"
  ),
  all.x = TRUE,
  sort = FALSE,
  allow.cartesian = TRUE
)

Q_UNMATCHED <- Q[
  is.na(daily_match),
  .N
]

Q[
  ,
  daily_match :=
    NULL
]

if (Q_UNMATCHED > 0L) {
  stopf(
    "准逐日正式样本有 ",
    Q_UNMATCHED,
    " 行无法连接07H01逐日暴露。"
  )
}

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

FP_HEADER <- names(
  fread(
    file = FP_FILE,
    nrows = 0,
    showProgress = FALSE
  )
)

FP_MISS <- setdiff(
  FP_REQUIRED,
  FP_HEADER
)

if (length(FP_MISS) > 0L) {
  stopf(
    "07F filter-period文件缺少字段：",
    paste(
      FP_MISS,
      collapse = ", "
    )
  )
}

cat(
  "\n读取07F filter-period数据...\n"
)

FP <- fread(
  file = FP_FILE,
  select = FP_REQUIRED,
  showProgress = FALSE
)

FP[
  ,
  `:=`(
    loc_num =
      as.character(loc_num),

    collect_start =
      as.IDate(collect_start),

    collect_end =
      as.IDate(collect_end)
  )
]

N_FP_ALL <- nrow(FP)

FP[
  ,
  flag_07g_filter_analysis_base :=
    flag_conus %in% TRUE &
    loc_num %in% ELIGIBLE_LOCS &
    is.finite(gross_beta_bq_m3) &
    gross_beta_bq_m3 < 0.04 &
    is.finite(pm_filter_period_official) &
    pm_filter_period_official < 35 &
    is.finite(vel) &
    vel < 6 &
    is.finite(hpbl) &
    hpbl > 150 &
    hpbl < 3000 &
    is.finite(c_prop) &
    is.finite(soilm) &
    is.finite(
      fpmean_campd_upwind_heat_observed_20km_90deg
    )
]

FP <- FP[
  flag_07g_filter_analysis_base %in%
    TRUE
]

FP[
  ,
  period_id :=
    .I
]

FP[
  ,
  n_period_days :=
    as.integer(
      collect_end -
      collect_start
    )
]

if (
  FP[
    !is.finite(n_period_days) |
    n_period_days <= 0,
    .N
  ] > 0L
) {
  stopf(
    "filter-period存在非法时间区间。"
  )
}

cat(
  "filter-period全部行 = ",
  format(
    N_FP_ALL,
    big.mark = ","
  ),
  "\n",
  "正式分析样本 = ",
  format(
    nrow(FP),
    big.mark = ","
  ),
  "\n",
  "正式监测站 = ",
  uniqueN(FP$loc_num),
  "\n",
  sep = ""
)

cat(
  "\n按冻结规则重构filter-period燃料暴露...\n"
)

FP_DAYS <- FP[
  ,
  .(
    date =
      seq(
        collect_start,
        collect_end - 1L,
        by = "day"
      )
  ),
  by = .(
    period_id,
    loc_num
  )
]

FP_DAYS <- merge(
  FP_DAYS,
  DAILY,
  by = c(
    "loc_num",
    "date"
  ),
  all.x = TRUE,
  sort = FALSE
)

FP_DAY_UNMATCHED <- FP_DAYS[
  is.na(daily_match),
  .N
]

if (FP_DAY_UNMATCHED > 0L) {
  stopf(
    "filter-period展开后有 ",
    FP_DAY_UNMATCHED,
    " 个monitor-day无法连接07H01暴露。"
  )
}

FP_EXP <- FP_DAYS[
  ,
  lapply(
    .SD,
    strict_mean
  ),
  by = period_id,
  .SDcols = FUEL_EXPOSURES
]

setnames(
  FP_EXP,
  FUEL_EXPOSURES,
  paste0(
    "fpmean_",
    FUEL_EXPOSURES
  )
)

FP <- merge(
  FP,
  FP_EXP,
  by = "period_id",
  all.x = TRUE,
  sort = FALSE
)

rm(
  FP_DAYS,
  FP_EXP
)

gc()

prepare_data <- function(
  DT,
  branch,
  fuel,
  metric
) {

  d <- copy(DT)

  if (metric == "电厂数量") {

    raw_col <- paste0(
      "campd_",
      fuel,
      "_upwind_n_plants_20km_90deg"
    )

    increment <- 1
    unit_label <- "每增加1座活动电厂"

  } else {

    raw_col <- paste0(
      "campd_",
      fuel,
      "_upwind_heat_observed_20km_90deg"
    )

    increment <- 10000
    unit_label <- "每增加10000 MMBtu实际热输入"
  }

  exp_col <- if (
    branch == "quasi_daily"
  ) {
    raw_col
  } else {
    paste0(
      "fpmean_",
      raw_col
    )
  }

  if (!exp_col %in% names(d)) {
    stopf(
      "缺少暴露字段：",
      exp_col
    )
  }

  if (branch == "quasi_daily") {

    d[
      ,
      `:=`(
        pm_model_centered =
          safe_num(pm_same_day) -
          10,

        analysis_date =
          as.IDate(date)
      )
    ]

  } else {

    d[
      ,
      `:=`(
        pm_model_centered =
          safe_num(
            pm_filter_period_official
          ) -
          10,

        analysis_date =
          as.IDate(
            collect_start
          )
      )
    ]
  }

  d[
    ,
    `:=`(
      gross_beta_mBq_m3 =
        safe_num(
          gross_beta_mBq_m3
        ),

      exposure =
        safe_num(
          get(exp_col)
        ) /
        increment,

      c_prop =
        safe_num(c_prop),

      surface_u_ppm_article =
        safe_num(
          surface_u_ppm_article
        ),

      sunspot_monthly =
        safe_num(
          sunspot_monthly
        ),

      inv_hpbl =
        safe_num(inv_hpbl),

      rhum =
        safe_num(rhum),

      soilm =
        safe_num(soilm),

      vel =
        safe_num(vel),

      lat_diff_35 =
        safe_num(lat_diff_35),

      year_model =
        safe_num(year_model),

      temp_centered =
        safe_num(temp) -
        273.15
    )
  ]

  vars <- c(
    "gross_beta_mBq_m3",
    "exposure",
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
    "loc_num",
    "analysis_date"
  )

  d <- d[
    complete.cases(
      d[
        ,
        ..vars
      ]
    )
  ]

  d[
    ,
    calendar_year :=
      as.integer(
        format(
          analysis_date,
          "%Y"
        )
      )
  ]

  d[
    ,
    monitor_year :=
      paste(
        loc_num,
        calendar_year,
        sep = "_"
      )
  ]

  d[
    ,
    exposure_between :=
      mean(
        exposure
      ),
    by = monitor_year
  ]

  d[
    ,
    exposure_within :=
      exposure -
      exposure_between
  ]

  d[
    ,
    loc_num_factor :=
      factor(
        loc_num
      )
  ]

  list(
    data = d,
    exposure_col = exp_col,
    increment = increment,
    unit_label = unit_label
  )
}

fit_one <- function(
  DT,
  branch,
  fuel,
  metric
) {

  P <- prepare_data(
    DT,
    branch,
    fuel,
    metric
  )

  d <- P$data

  meta <- data.table(
    branch = branch,
    fuel = fuel,
    fuel_label =
      unname(
        FUEL_CN[fuel]
      ),
    metric = metric,
    exposure_col =
      P$exposure_col,
    unit_label =
      P$unit_label,
    n_rows = nrow(d),
    n_monitors =
      uniqueN(
        d$loc_num
      ),
    n_monitor_years =
      uniqueN(
        d$monitor_year
      ),
    exposure_sd =
      sd(
        d$exposure,
        na.rm = TRUE
      ),
    within_sd =
      sd(
        d$exposure_within,
        na.rm = TRUE
      ),
    between_sd =
      sd(
        d$exposure_between,
        na.rm = TRUE
      ),
    wind_min =
      min(
        d$vel,
        na.rm = TRUE
      ),
    wind_q01 =
      unname(
        quantile(
          d$vel,
          0.01,
          na.rm = TRUE
        )
      ),
    wind_q25 =
      unname(
        quantile(
          d$vel,
          0.25,
          na.rm = TRUE
        )
      ),
    wind_median =
      unname(
        quantile(
          d$vel,
          0.50,
          na.rm = TRUE
        )
      ),
    wind_mean =
      mean(
        d$vel,
        na.rm = TRUE
      ),
    wind_q75 =
      unname(
        quantile(
          d$vel,
          0.75,
          na.rm = TRUE
        )
      ),
    wind_q99 =
      unname(
        quantile(
          d$vel,
          0.99,
          na.rm = TRUE
        )
      ),
    wind_max =
      max(
        d$vel,
        na.rm = TRUE
      )
  )

  if (
    nrow(d) < 1000L ||
    uniqueN(d$loc_num) < 2L ||
    !is.finite(
      sd(
        d$exposure,
        na.rm = TRUE
      )
    ) ||
    sd(
      d$exposure,
      na.rm = TRUE
    ) <= 0
  ) {

    return(
      list(
        status =
          cbind(
            meta,
            data.table(
              paper_model_status =
                "SKIP",
              decomposition_status =
                "SKIP",
              paper_note =
                "暴露变异不足",
              decomposition_note =
                "暴露变异不足"
            )
          ),
        reference =
          data.table(),
        decomposition =
          data.table(),
        curves =
          data.table(),
        transition =
          data.table()
      )
    )
  }

  FORM_PAPER <- as.formula(
    paste0(
      "gross_beta_mBq_m3 ~ ",
      "exposure + ",
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
      "temp_centered + I(temp_centered^2) + I(temp_centered^3) + ",
      "sunspot_monthly:lat_diff_35 + ",
      "surface_u_ppm_article:soilm + ",
      "exposure:vel + ",
      "(1 | loc_num_factor)"
    )
  )

  FORM_DECOMP <- as.formula(
    paste0(
      "gross_beta_mBq_m3 ~ ",
      "exposure_within + ",
      "exposure_between + ",
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
      "temp_centered + I(temp_centered^2) + I(temp_centered^3) + ",
      "sunspot_monthly:lat_diff_35 + ",
      "surface_u_ppm_article:soilm + ",
      "exposure_within:vel + ",
      "exposure_between:vel + ",
      "(1 | loc_num_factor)"
    )
  )

  warning_paper <- character()

  fit_paper <- tryCatch(

    withCallingHandlers(

      lmer(
        FORM_PAPER,
        data = d,
        REML = TRUE,
        na.action = na.fail,
        control =
          make_lmer_control()
      ),

      warning = function(w) {

        warning_paper <<- c(
          warning_paper,
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

  warning_decomp <- character()

  fit_decomp <- tryCatch(

    withCallingHandlers(

      lmer(
        FORM_DECOMP,
        data = d,
        REML = TRUE,
        na.action = na.fail,
        control =
          make_lmer_control()
      ),

      warning = function(w) {

        warning_decomp <<- c(
          warning_decomp,
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

  paper_status <- if (
    inherits(
      fit_paper,
      "fit_error"
    )
  ) {
    "FAIL"
  } else {
    "PASS"
  }

  decomp_status <- if (
    inherits(
      fit_decomp,
      "fit_error"
    )
  ) {
    "FAIL"
  } else {
    "PASS"
  }

  STATUS <- cbind(
    meta,
    data.table(
      paper_model_status =
        paper_status,

      decomposition_status =
        decomp_status,

      paper_singular =
        if (
          paper_status ==
            "PASS"
        ) {
          as.integer(
            isSingular(
              fit_paper,
              tol = 1e-4
            )
          )
        } else {
          NA_integer_
        },

      decomposition_singular =
        if (
          decomp_status ==
            "PASS"
        ) {
          as.integer(
            isSingular(
              fit_decomp,
              tol = 1e-4
            )
          )
        } else {
          NA_integer_
        },

      paper_convergence =
        if (
          paper_status ==
            "PASS"
        ) {
          collapse_msg(
            fit_paper@optinfo$conv$lme4$messages
          )
        } else {
          fit_paper$message
        },

      decomposition_convergence =
        if (
          decomp_status ==
            "PASS"
        ) {
          collapse_msg(
            fit_decomp@optinfo$conv$lme4$messages
          )
        } else {
          fit_decomp$message
        },

      paper_warnings =
        collapse_msg(
          warning_paper
        ),

      decomposition_warnings =
        collapse_msg(
          warning_decomp
        )
    )
  )

  REF <- data.table()
  DEC <- data.table()
  CURVE <- data.table()
  TRANS <- data.table()

  if (
    paper_status ==
      "PASS"
  ) {

    bp <- fixef(
      fit_paper
    )

    Vp <- as.matrix(
      vcov(
        fit_paper
      )
    )

    intp <- find_interaction(
      names(bp),
      "exposure",
      "vel"
    )

    reference_winds <- data.table(
      wind_reference = c(
        "参考论文1米每秒",
        "样本第25百分位",
        "样本中位数",
        "样本第75百分位"
      ),

      wind_ms = c(
        1,
        meta$wind_q25,
        meta$wind_median,
        meta$wind_q75
      )
    )

    REF <- rbindlist(
      lapply(
        seq_len(
          nrow(
            reference_winds
          )
        ),
        function(ii) {

          eff <- linear_effect(
            bp,
            Vp,
            "exposure",
            intp,
            reference_winds$wind_ms[ii]
          )

          data.table(
            branch = branch,
            fuel = fuel,
            fuel_label =
              unname(
                FUEL_CN[fuel]
              ),
            metric = metric,
            unit_label =
              P$unit_label,
            model =
              "参考论文式随机截距模型",
            component =
              "总体暴露",
            wind_reference =
              reference_winds$wind_reference[ii],
            wind_ms =
              reference_winds$wind_ms[ii],
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
        }
      )
    )
  }

  if (
    decomp_status ==
      "PASS"
  ) {

    bd <- fixef(
      fit_decomp
    )

    Vd <- as.matrix(
      vcov(
        fit_decomp
      )
    )

    int_within <- find_interaction(
      names(bd),
      "exposure_within",
      "vel"
    )

    int_between <- find_interaction(
      names(bd),
      "exposure_between",
      "vel"
    )

    reference_winds <- data.table(
      wind_reference = c(
        "参考论文1米每秒",
        "样本第25百分位",
        "样本中位数",
        "样本第75百分位"
      ),

      wind_ms = c(
        1,
        meta$wind_q25,
        meta$wind_median,
        meta$wind_q75
      )
    )

    DEC <- rbindlist(
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
              int_within

          } else {

            main_term <-
              "exposure_between"

            int_term <-
              int_between
          }

          rbindlist(
            lapply(
              seq_len(
                nrow(
                  reference_winds
                )
              ),
              function(ii) {

                eff <- linear_effect(
                  bd,
                  Vd,
                  main_term,
                  int_term,
                  reference_winds$wind_ms[ii]
                )

                data.table(
                  branch = branch,
                  fuel = fuel,
                  fuel_label =
                    unname(
                      FUEL_CN[fuel]
                    ),
                  metric = metric,
                  unit_label =
                    P$unit_label,
                  model =
                    "站点年份内外分解模型",
                  component =
                    component,
                  wind_reference =
                    reference_winds$wind_reference[ii],
                  wind_ms =
                    reference_winds$wind_ms[ii],
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
              }
            )
          )
        }
      )
    )

    WIND_GRID <- seq(
      from =
        meta$wind_q01,
      to =
        meta$wind_q99,
      length.out =
        201
    )

    CURVE <- rbindlist(
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
              int_within

          } else {

            main_term <-
              "exposure_between"

            int_term <-
              int_between
          }

          rbindlist(
            lapply(
              WIND_GRID,
              function(w) {

                eff <- linear_effect(
                  bd,
                  Vd,
                  main_term,
                  int_term,
                  w
                )

                data.table(
                  branch = branch,
                  fuel = fuel,
                  fuel_label =
                    unname(
                      FUEL_CN[fuel]
                    ),
                  metric = metric,
                  unit_label =
                    P$unit_label,
                  component =
                    component,
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
              }
            )
          )
        }
      )
    )

    make_transition <- function(
      component,
      main_term,
      int_term
    ) {

      if (
        !main_term %in% names(bd) ||
        is.na(int_term) ||
        !int_term %in% names(bd)
      ) {

        return(
          data.table(
            branch = branch,
            fuel = fuel,
            fuel_label =
              unname(
                FUEL_CN[fuel]
              ),
            metric = metric,
            component = component,
            main_beta = NA_real_,
            interaction_beta = NA_real_,
            zero_crossing_wind_ms =
              NA_real_,
            zero_crossing_inside_q01_q99 =
              NA,
            interpretation =
              "无法计算"
          )
        )
      }

      b0 <- bd[[main_term]]
      bi <- bd[[int_term]]

      zc <- if (
        is.finite(bi) &&
        abs(bi) >
          .Machine$double.eps
      ) {
        -b0 / bi
      } else {
        NA_real_
      }

      inside <-
        is.finite(zc) &&
        zc >=
          meta$wind_q01 &&
        zc <=
          meta$wind_q99

      txt <- if (
        !is.finite(zc)
      ) {
        "实际风速范围内没有可定义的效应过零点"
      } else if (
        inside
      ) {
        "效应过零点位于样本主要风速范围内"
      } else {
        "效应过零点位于样本主要风速范围外"
      }

      data.table(
        branch = branch,
        fuel = fuel,
        fuel_label =
          unname(
            FUEL_CN[fuel]
          ),
        metric = metric,
        component = component,
        main_beta = b0,
        interaction_beta = bi,
        zero_crossing_wind_ms =
          zc,
        zero_crossing_inside_q01_q99 =
          inside,
        interpretation =
          txt
      )
    }

    TRANS <- rbindlist(
      list(
        make_transition(
          "站点年份内变化",
          "exposure_within",
          int_within
        ),

        make_transition(
          "站点年份间平均差异",
          "exposure_between",
          int_between
        )
      ),
      use.names = TRUE
    )
  }

  list(
    status = STATUS,
    reference = REF,
    decomposition = DEC,
    curves = CURVE,
    transition = TRANS
  )
}

STATUS_LIST <- list()
REF_LIST <- list()
DEC_LIST <- list()
CURVE_LIST <- list()
TRANS_LIST <- list()

kk <- 0L

cat(
  "\n开始20公里、90°上风向诊断模型...\n"
)

for (
  branch in c(
    "filter_period",
    "quasi_daily"
  )
) {

  DT <- if (
    branch ==
      "filter_period"
  ) {
    FP
  } else {
    Q
  }

  for (fuel in FUELS) {

    for (
      metric in c(
        "电厂数量",
        "实际热输入"
      )
    ) {

      kk <- kk + 1L

      cat(
        "[",
        kk,
        "/24] ",
        branch,
        " | ",
        unname(
          FUEL_CN[fuel]
        ),
        " | ",
        metric,
        "\n",
        sep = ""
      )

      z <- fit_one(
        DT,
        branch,
        fuel,
        metric
      )

      STATUS_LIST[[kk]] <-
        z$status

      REF_LIST[[kk]] <-
        z$reference

      DEC_LIST[[kk]] <-
        z$decomposition

      CURVE_LIST[[kk]] <-
        z$curves

      TRANS_LIST[[kk]] <-
        z$transition
    }
  }
}

STATUS <- rbindlist(
  STATUS_LIST,
  use.names = TRUE,
  fill = TRUE
)

REF <- rbindlist(
  REF_LIST,
  use.names = TRUE,
  fill = TRUE
)

DEC <- rbindlist(
  DEC_LIST,
  use.names = TRUE,
  fill = TRUE
)

CURVE <- rbindlist(
  CURVE_LIST,
  use.names = TRUE,
  fill = TRUE
)

TRANS <- rbindlist(
  TRANS_LIST,
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  STATUS,
  file.path(
    OUT,
    "01_模型状态与实际风速分布.csv"
  )
)

fwrite(
  REF,
  file.path(
    OUT,
    "02_参考论文式模型_代表风速效应.csv"
  )
)

fwrite(
  DEC,
  file.path(
    OUT,
    "03_站点年份内外分解_代表风速效应.csv"
  )
)

fwrite(
  CURVE,
  file.path(
    OUT,
    "04_站点年份内外效应_连续风速曲线.csv"
  )
)

fwrite(
  TRANS,
  file.path(
    OUT,
    "05_效应过零风速审计.csv"
  )
)

OLD <- fread(
  file = H07_RESULT,
  showProgress = FALSE
)

OLD20 <- OLD[
  radius_km ==
    20 &
  component %in%
    c(
      "raw_count_only",
      "raw_heat_only"
    )
]

OLD20[
  ,
  metric :=
    fifelse(
      component ==
        "raw_count_only",
      "电厂数量",
      "实际热输入"
    )
]

NEW1 <- REF[
  wind_reference ==
    "参考论文1米每秒",
  .(
    branch,
    fuel,
    metric,
    new_estimate =
      estimate,
    new_ci_low =
      ci95_low,
    new_ci_high =
      ci95_high
  )
]

COMPARE <- merge(
  OLD20[
    ,
    .(
      branch,
      fuel,
      metric,
      old_status =
        status,
      old_estimate =
        estimate_fixed,
      old_ci_low =
        ci_low_fixed,
      old_ci_high =
        ci_high_fixed
    )
  ],
  NEW1,
  by = c(
    "branch",
    "fuel",
    "metric"
  ),
  all = TRUE,
  sort = FALSE
)

COMPARE[
  ,
  `:=`(
    estimate_abs_difference =
      abs(
        old_estimate -
        new_estimate
      ),

    ci_low_abs_difference =
      abs(
        old_ci_low -
        new_ci_low
      ),

    ci_high_abs_difference =
      abs(
        old_ci_high -
        new_ci_high
      )
  )
]

fwrite(
  COMPARE,
  file.path(
    OUT,
    "06_冻结07H与07H09参考论文式模型_1米每秒复现审计.csv"
  )
)

KEY <- merge(
  REF[
    wind_reference ==
      "参考论文1米每秒",
    .(
      branch,
      fuel,
      fuel_label,
      metric,
      paper_estimate_wind1 =
        estimate,
      paper_ci_low_wind1 =
        ci95_low,
      paper_ci_high_wind1 =
        ci95_high,
      paper_direction_wind1 =
        direction
    )
  ],
  DEC[
    wind_reference ==
      "参考论文1米每秒" &
    component ==
      "站点年份内变化",
    .(
      branch,
      fuel,
      metric,
      within_estimate_wind1 =
        estimate,
      within_ci_low_wind1 =
        ci95_low,
      within_ci_high_wind1 =
        ci95_high,
      within_direction_wind1 =
        direction
    )
  ],
  by = c(
    "branch",
    "fuel",
    "metric"
  ),
  all = TRUE,
  sort = FALSE
)

KEY <- merge(
  KEY,
  DEC[
    wind_reference ==
      "参考论文1米每秒" &
    component ==
      "站点年份间平均差异",
    .(
      branch,
      fuel,
      metric,
      between_estimate_wind1 =
        estimate,
      between_ci_low_wind1 =
        ci95_low,
      between_ci_high_wind1 =
        ci95_high,
      between_direction_wind1 =
        direction
    )
  ],
  by = c(
    "branch",
    "fuel",
    "metric"
  ),
  all = TRUE,
  sort = FALSE
)

KEY[
  ,
  paper_vs_within :=
    fifelse(
      is.finite(
        paper_estimate_wind1
      ) &
      is.finite(
        within_estimate_wind1
      ),
      fifelse(
        sign(
          paper_estimate_wind1
        ) ==
        sign(
          within_estimate_wind1
        ),
        "方向相同",
        "方向改变"
      ),
      "无法比较"
    )
]

fwrite(
  KEY,
  file.path(
    OUT,
    "07_最关键结果_原模型与站点年份内外方向比较.csv"
  )
)

MEDIAN_WITHIN <- DEC[
  wind_reference ==
    "样本中位数" &
  component ==
    "站点年份内变化"
]

fwrite(
  MEDIAN_WITHIN,
  file.path(
    OUT,
    "08_站点年份内效应_实际样本中位风速.csv"
  )
)

N_PAPER_FAIL <- STATUS[
  paper_model_status ==
    "FAIL",
  .N
]

N_DEC_FAIL <- STATUS[
  decomposition_status ==
    "FAIL",
  .N
]

N_PAPER_PASS <- STATUS[
  paper_model_status ==
    "PASS",
  .N
]

N_DEC_PASS <- STATUS[
  decomposition_status ==
    "PASS",
  .N
]

MAX_EST_DIFF <- if (
  COMPARE[
    old_status ==
      "PASS" &
    is.finite(
      new_estimate
    ),
    .N
  ] >
  0L
) {

  max(
    COMPARE[
      old_status ==
        "PASS" &
      is.finite(
        new_estimate
      ),
      estimate_abs_difference
    ],
    na.rm = TRUE
  )

} else {

  NA_real_
}

MAX_CI_DIFF <- if (
  nrow(COMPARE) >
    0L
) {

  max(
    c(
      COMPARE$ci_low_abs_difference,
      COMPARE$ci_high_abs_difference
    ),
    na.rm = TRUE
  )

} else {

  NA_real_
}

HARD <- data.table(
  check = c(
    "07H01_daily_files_24",
    "daily_duplicate_keys_zero",
    "quasi_primary_rows_426657",
    "quasi_primary_monitors_133",
    "quasi_daily_join_unmatched_zero",
    "filter_all_rows_174200",
    "filter_formal_rows_107407",
    "filter_formal_monitors_133",
    "filter_daily_join_unmatched_zero",
    "model_specs_24",
    "paper_model_FAIL_zero",
    "decomposition_model_FAIL_zero",
    "reference_model_has_results",
    "decomposition_model_has_results",
    "continuous_curve_has_results",
    "frozen_07H_comparison_rows_24"
  ),

  observed = c(
    as.character(
      sum(
        file.exists(
          DAILY_FILES
        )
      )
    ),
    as.character(
      DAILY_DUP
    ),
    as.character(
      N_Q_PREJOIN
    ),
    as.character(
      length(
        ELIGIBLE_LOCS
      )
    ),
    as.character(
      Q_UNMATCHED
    ),
    as.character(
      N_FP_ALL
    ),
    as.character(
      nrow(FP)
    ),
    as.character(
      uniqueN(
        FP$loc_num
      )
    ),
    as.character(
      FP_DAY_UNMATCHED
    ),
    as.character(
      nrow(STATUS)
    ),
    as.character(
      N_PAPER_FAIL
    ),
    as.character(
      N_DEC_FAIL
    ),
    as.character(
      as.integer(
        nrow(REF) >
          0L
      )
    ),
    as.character(
      as.integer(
        nrow(DEC) >
          0L
      )
    ),
    as.character(
      as.integer(
        nrow(CURVE) >
          0L
      )
    ),
    as.character(
      nrow(COMPARE)
    )
  ),

  expected = c(
    "24",
    "0",
    "426657",
    "133",
    "0",
    "174200",
    "107407",
    "133",
    "0",
    "24",
    "0",
    "0",
    "1",
    "1",
    "1",
    "24"
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
    OUT,
    "09_07H09_hard_checks.csv"
  )
)

cat(
  "\n============================================================\n",
  "A. 07H09模型状态\n",
  "============================================================\n",
  sep = ""
)

print(
  STATUS[
    ,
    .(
      branch,
      fuel,
      metric,
      paper_model_status,
      decomposition_status,
      n_rows,
      n_monitors,
      n_monitor_years,
      wind_q25,
      wind_median,
      wind_q75
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "B. 冻结07H结果与本轮论文式模型复现审计\n",
  "============================================================\n",
  sep = ""
)

cat(
  "最大效应绝对差 = ",
  format(
    MAX_EST_DIFF,
    scientific = TRUE,
    digits = 8
  ),
  "\n",
  "最大95%CI端点绝对差 = ",
  format(
    MAX_CI_DIFF,
    scientific = TRUE,
    digits = 8
  ),
  "\n",
  sep = ""
)

print(
  COMPARE[
    ,
    .(
      branch,
      fuel,
      metric,
      old_estimate,
      new_estimate,
      estimate_abs_difference
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "C. 最关键结果：1米/秒时 原模型 vs 同站同年内 vs 站点年份间\n",
  "============================================================\n",
  sep = ""
)

print(
  KEY[
    ,
    .(
      branch,
      fuel,
      metric,
      paper_estimate_wind1,
      within_estimate_wind1,
      between_estimate_wind1,
      paper_vs_within
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "D. 最关键结果：实际样本中位风速下的同站同年内效应\n",
  "============================================================\n",
  sep = ""
)

print(
  MEDIAN_WITHIN[
    ,
    .(
      branch,
      fuel,
      metric,
      wind_ms,
      estimate,
      ci95_low,
      ci95_high,
      direction,
      ci_status
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "E. 连续风速效应是否发生正负方向转换\n",
  "============================================================\n",
  sep = ""
)

print(
  TRANS[
    component ==
      "站点年份内变化",
    .(
      branch,
      fuel,
      metric,
      main_beta,
      interaction_beta,
      zero_crossing_wind_ms,
      zero_crossing_inside_q01_q99,
      interpretation
    )
  ],
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

HARD_FAIL <- HARD[
  status ==
    "FAIL",
  .N
]

REPORT <- c(
  "# 07H09 站点年份内效应与连续风速诊断",
  "",
  paste0(
    "- 生成时间：",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  "",
  "## 分析定位",
  "",
  "- 不修改07H01至07H08任何冻结结果。",
  "- 固定20 km、90°上风向。",
  "- 第一层重新拟合参考论文式随机截距模型。",
  "- 第二层将燃料暴露拆成同一监测站同一年内部变化和站点年份平均差异。",
  "- 第三层保留连续风速与暴露交互，并在实际样本第1至第99百分位风速范围生成连续效应曲线。",
  "- 1米/秒仅用于与参考论文和冻结07H结果直接比较。",
  "- 样本第25百分位、中位数和第75百分位用于描述本研究真实风速范围。",
  "",
  "## 主要解释对象",
  "",
  "- 站点年份内变化效应：同一监测站、同一年内，某天或某滤膜期的上风向暴露高于该站该年平均水平时，总贝塔如何变化。",
  "- 站点年份间平均差异：包含较强的地区和长期结构信息，仅用于解释原模型结果来源。",
  "",
  paste0(
    "- 论文式模型PASS：",
    N_PAPER_PASS,
    "/24"
  ),
  paste0(
    "- 分解模型PASS：",
    N_DEC_PASS,
    "/24"
  ),
  paste0(
    "- 最大冻结07H效应复现差：",
    format(
      MAX_EST_DIFF,
      scientific = TRUE,
      digits = 8
    )
  ),
  "",
  if (HARD_FAIL == 0L) {
    "**07H09 DIAGNOSTIC: PASS**"
  } else {
    "**07H09 DIAGNOSTIC: FAIL**"
  }
)

writeLines(
  REPORT,
  file.path(
    OUT,
    "10_07H09_诊断说明.md"
  ),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUT,
    "11_sessionInfo.txt"
  )
)

cat(
  "\n============================================================\n"
)

if (HARD_FAIL == 0L) {
  cat(
    "07H09 DIAGNOSTIC: PASS\n"
  )
} else {
  cat(
    "07H09 DIAGNOSTIC: FAIL\n"
  )
}

cat(
  "============================================================\n",
  "论文式模型PASS = ",
  N_PAPER_PASS,
  "/24\n",
  "站点年份内外分解模型PASS = ",
  N_DEC_PASS,
  "/24\n",
  "Hard FAIL = ",
  HARD_FAIL,
  "\n",
  sep = ""
)

if (HARD_FAIL > 0L) {
  stopf(
    "07H09存在硬性检查失败，先不要解释结果。"
  )
}
