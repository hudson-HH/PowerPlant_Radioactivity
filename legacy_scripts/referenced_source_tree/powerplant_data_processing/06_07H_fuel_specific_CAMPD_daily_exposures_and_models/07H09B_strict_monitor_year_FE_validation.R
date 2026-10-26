#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999
)

required_packages <- c(
  "data.table",
  "fixest"
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

strict_mean <- function(x) {

  if (
    length(x) == 0L ||
    any(!is.finite(x))
  ) {
    return(NA_real_)
  }

  mean(x)
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

  est <-
    beta[[main_term]] +
    wind_ms *
      beta[[interaction_term]]

  vv <-
    V[main_term, main_term] +
    wind_ms^2 *
      V[interaction_term, interaction_term] +
    2 *
      wind_ms *
      V[main_term, interaction_term]

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
    est - 1.96 * se
  } else {
    NA_real_
  }

  high <- if (
    is.finite(est) &&
    is.finite(se)
  ) {
    est + 1.96 * se
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
        abs(est / se),
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

H09 <- file.path(
  HROOT,
  "07H09_站点年份内效应与连续风速诊断"
)

H09_WITHIN <- file.path(
  H09,
  "03_站点年份内外分解_代表风速效应.csv"
)

OUT <- file.path(
  HROOT,
  "07H09B_监测站年份固定效应严格验证"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

for (
  f in c(
    FP_FILE,
    H09_WITHIN
  )
) {

  if (!file.exists(f)) {
    stopf(
      "缺少输入：\n",
      f
    )
  }
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

EXPOSURES <- c(
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
    "07H01逐日暴露文件不是24/24。"
  )
}

cat(
  "\n读取07H01逐日燃料暴露...\n"
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
        EXPOSURES
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
    "07H01逐日暴露存在重复monitor-day。"
  )
}

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
  "\n读取07F准逐日正式样本...\n"
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

ELIGIBLE_LOCS <- sort(
  unique(
    Q$loc_num
  )
)

N_Q <- nrow(Q)

Q <- merge(
  Q,
  DAILY,
  by = c(
    "loc_num",
    "date"
  ),
  all.x = TRUE,
  sort = FALSE
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
    "准逐日样本暴露连接缺失。"
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

HFP <- names(
  fread(
    file = FP_FILE,
    nrows = 0,
    showProgress = FALSE
  )
)

miss_fp <- setdiff(
  FP_REQUIRED,
  HFP
)

if (length(miss_fp) > 0L) {
  stopf(
    "filter-period缺少字段：",
    paste(
      miss_fp,
      collapse = ", "
    )
  )
}

cat(
  "\n读取07F滤膜采样期正式样本...\n"
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
  flag_formal :=
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
  flag_formal %in%
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
    "发现非法滤膜采样期。"
  )
}

cat(
  "filter-period正式样本 = ",
  format(
    nrow(FP),
    big.mark = ","
  ),
  "\n",
  sep = ""
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

FP_UNMATCHED <- FP_DAYS[
  is.na(daily_match),
  .N
]

if (FP_UNMATCHED > 0L) {
  stopf(
    "滤膜采样期逐日暴露连接缺失。"
  )
}

FP_EXP <- FP_DAYS[
  ,
  lapply(
    .SD,
    strict_mean
  ),
  by = period_id,
  .SDcols = EXPOSURES
]

setnames(
  FP_EXP,
  EXPOSURES,
  paste0(
    "fpmean_",
    EXPOSURES
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

    unit_label <-
      "每增加1座活动电厂"

  } else {

    raw_col <- paste0(
      "campd_",
      fuel,
      "_upwind_heat_observed_20km_90deg"
    )

    increment <- 10000

    unit_label <-
      "每增加10000 MMBtu实际热输入"
  }

  exp_col <- if (
    branch ==
      "quasi_daily"
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

  if (
    branch ==
      "quasi_daily"
  ) {

    d[
      ,
      `:=`(
        analysis_date =
          as.IDate(date),

        pm_model_centered =
          safe_num(pm_same_day) -
          10
      )
    ]

  } else {

    d[
      ,
      `:=`(
        analysis_date =
          as.IDate(
            collect_start
          ),

        pm_model_centered =
          safe_num(
            pm_filter_period_official
          ) -
          10
      )
    ]
  }

  d[
    ,
    `:=`(
      exposure =
        safe_num(
          get(exp_col)
        ) /
        increment,

      gross_beta_mBq_m3 =
        safe_num(
          gross_beta_mBq_m3
        ),

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
    monitor_year_fe :=
      interaction(
        loc_num,
        calendar_year,
        drop = TRUE,
        lex.order = TRUE
      )
  ]

  d[
    ,
    date_fe :=
      factor(
        as.character(
          analysis_date
        )
      )
  ]

  list(
    data = d,
    exposure_col = exp_col,
    unit_label = unit_label
  )
}

fit_fe <- function(
  DT,
  branch,
  fuel,
  metric,
  fe_structure
) {

  P <- prepare_data(
    DT,
    branch,
    fuel,
    metric
  )

  d <- P$data

  n_rows_pre <- nrow(d)

  n_monitor_years <-
    uniqueN(
      d$monitor_year_fe
    )

  exposure_sd <-
    sd(
      d$exposure,
      na.rm = TRUE
    )

  within_var <- d[
    ,
    var(
      exposure,
      na.rm = TRUE
    ),
    by = monitor_year_fe
  ]

  n_monitor_years_with_variation <-
    within_var[
      is.finite(V1) &
      V1 > 0,
      .N
    ]

  meta <- data.table(
    branch = branch,
    fuel = fuel,
    fuel_label =
      unname(
        FUEL_CN[fuel]
      ),
    metric = metric,
    unit_label =
      P$unit_label,
    fe_structure =
      fe_structure,
    n_rows_pre =
      n_rows_pre,
    n_monitors =
      uniqueN(
        d$loc_num
      ),
    n_monitor_years =
      n_monitor_years,
    n_monitor_years_with_exposure_variation =
      n_monitor_years_with_variation,
    exposure_sd =
      exposure_sd,
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
      )
  )

  if (
    n_rows_pre < 1000L ||
    uniqueN(d$loc_num) < 2L ||
    !is.finite(exposure_sd) ||
    exposure_sd <= 0 ||
    n_monitor_years_with_variation <
      2L
  ) {

    return(
      list(
        status =
          cbind(
            meta,
            data.table(
              status = "SKIP",
              nobs = NA_integer_,
              note =
                "同站同年内部暴露变异不足"
            )
          ),
        effects =
          data.table(),
        curve =
          data.table(),
        transition =
          data.table()
      )
    )
  }

  rhs <- if (
    fe_structure ==
      "monitor_year_FE"
  ) {

    paste0(
      "exposure + exposure:vel + ",
      "c_prop + ",
      "sunspot_monthly + ",
      "pm_model_centered + ",
      "inv_hpbl + ",
      "rhum + ",
      "soilm + ",
      "vel + ",
      "temp_centered + ",
      "I(temp_centered^2) + ",
      "I(temp_centered^3) + ",
      "sunspot_monthly:lat_diff_35 + ",
      "surface_u_ppm_article:soilm"
    )

  } else {

    paste0(
      "exposure + exposure:vel + ",
      "c_prop + ",
      "pm_model_centered + ",
      "inv_hpbl + ",
      "rhum + ",
      "soilm + ",
      "vel + ",
      "temp_centered + ",
      "I(temp_centered^2) + ",
      "I(temp_centered^3) + ",
      "sunspot_monthly:lat_diff_35 + ",
      "surface_u_ppm_article:soilm"
    )
  }

  FORM <- if (
    fe_structure ==
      "monitor_year_FE"
  ) {

    as.formula(
      paste0(
        "gross_beta_mBq_m3 ~ ",
        rhs,
        " | monitor_year_fe"
      )
    )

  } else {

    as.formula(
      paste0(
        "gross_beta_mBq_m3 ~ ",
        rhs,
        " | monitor_year_fe + date_fe"
      )
    )
  }

  warning_text <- character()

  fit <- tryCatch(

    withCallingHandlers(

      feols(
        FORM,
        data = d,
        vcov = ~loc_num,
        notes = FALSE,
        warn = FALSE
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
      list(
        status =
          cbind(
            meta,
            data.table(
              status = "FAIL",
              nobs = NA_integer_,
              warnings =
                paste(
                  unique(
                    warning_text
                  ),
                  collapse = " | "
                ),
              note =
                fit$message
            )
          ),
        effects =
          data.table(),
        curve =
          data.table(),
        transition =
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
    is.na(int_term)
  ) {

    return(
      list(
        status =
          cbind(
            meta,
            data.table(
              status = "SKIP",
              nobs =
                nobs(fit),
              warnings =
                paste(
                  unique(
                    warning_text
                  ),
                  collapse = " | "
                ),
              note =
                "暴露主效应或暴露×风速交互被共线性删除"
            )
          ),
        effects =
          data.table(),
        curve =
          data.table(),
        transition =
          data.table()
      )
    )
  }

  REF_WIND <- data.table(
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

  EFFECTS <- rbindlist(
    lapply(
      seq_len(
        nrow(
          REF_WIND
        )
      ),
      function(ii) {

        eff <- linear_effect(
          beta,
          V,
          "exposure",
          int_term,
          REF_WIND$wind_ms[ii]
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
          fe_structure =
            fe_structure,
          wind_reference =
            REF_WIND$wind_reference[ii],
          wind_ms =
            REF_WIND$wind_ms[ii],
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
      WIND_GRID,
      function(w) {

        eff <- linear_effect(
          beta,
          V,
          "exposure",
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
          fe_structure =
            fe_structure,
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

  b0 <- beta[["exposure"]]
  bi <- beta[[int_term]]

  zero_cross <- if (
    is.finite(bi) &&
    abs(bi) >
      .Machine$double.eps
  ) {
    -b0 / bi
  } else {
    NA_real_
  }

  zero_inside <-
    is.finite(zero_cross) &&
    zero_cross >=
      meta$wind_q01 &&
    zero_cross <=
      meta$wind_q99

  TRANS <- data.table(
    branch = branch,
    fuel = fuel,
    fuel_label =
      unname(
        FUEL_CN[fuel]
      ),
    metric = metric,
    fe_structure =
      fe_structure,
    exposure_beta =
      b0,
    exposure_wind_interaction_beta =
      bi,
    zero_crossing_wind_ms =
      zero_cross,
    zero_crossing_inside_q01_q99 =
      zero_inside,
    interpretation =
      fifelse(
        !is.finite(
          zero_cross
        ),
        "没有可定义的过零风速",
        fifelse(
          zero_inside,
          "过零点位于样本主要风速范围内",
          "过零点位于样本主要风速范围外"
        )
      )
  )

  STATUS <- cbind(
    meta,
    data.table(
      status = "PASS",
      nobs =
        nobs(fit),
      warnings =
        paste(
          unique(
            warning_text
          ),
          collapse = " | "
        ),
      note =
        "固定效应模型；标准误按监测站聚类"
    )
  )

  list(
    status =
      STATUS,
    effects =
      EFFECTS,
    curve =
      CURVE,
    transition =
      TRANS
  )
}

STATUS_LIST <- list()
EFFECT_LIST <- list()
CURVE_LIST <- list()
TRANS_LIST <- list()

kk <- 0L

cat(
  "\n开始监测站×年份固定效应严格验证...\n"
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

      for (
        fe_structure in c(
          "monitor_year_FE",
          "monitor_year_plus_date_FE"
        )
      ) {

        kk <- kk + 1L

        cat(
          "[",
          kk,
          "/48] ",
          branch,
          " | ",
          unname(
            FUEL_CN[fuel]
          ),
          " | ",
          metric,
          " | ",
          fe_structure,
          "\n",
          sep = ""
        )

        z <- fit_fe(
          DT,
          branch,
          fuel,
          metric,
          fe_structure
        )

        STATUS_LIST[[kk]] <-
          z$status

        EFFECT_LIST[[kk]] <-
          z$effects

        CURVE_LIST[[kk]] <-
          z$curve

        TRANS_LIST[[kk]] <-
          z$transition
      }
    }
  }
}

STATUS <- rbindlist(
  STATUS_LIST,
  use.names = TRUE,
  fill = TRUE
)

EFFECTS <- rbindlist(
  EFFECT_LIST,
  use.names = TRUE,
  fill = TRUE
)

CURVES <- rbindlist(
  CURVE_LIST,
  use.names = TRUE,
  fill = TRUE
)

TRANS <- rbindlist(
  TRANS_LIST,
  use.names = TRUE,
  fill = TRUE
)

setorder(
  STATUS,
  branch,
  fuel,
  metric,
  fe_structure
)

setorder(
  EFFECTS,
  branch,
  fuel,
  metric,
  fe_structure,
  wind_ms
)

fwrite(
  STATUS,
  file.path(
    OUT,
    "01_固定效应模型状态.csv"
  )
)

fwrite(
  EFFECTS,
  file.path(
    OUT,
    "02_固定效应模型_代表风速结果.csv"
  )
)

fwrite(
  CURVES,
  file.path(
    OUT,
    "03_固定效应模型_连续风速效应曲线.csv"
  )
)

fwrite(
  TRANS,
  file.path(
    OUT,
    "04_固定效应模型_过零风速.csv"
  )
)

H09_WITHIN_DT <- fread(
  file = H09_WITHIN,
  showProgress = FALSE
)

H09_WITHIN_KEY <- H09_WITHIN_DT[
  component ==
    "站点年份内变化" &
  wind_reference %in%
    c(
      "参考论文1米每秒",
      "样本中位数"
    ),
  .(
    branch,
    fuel,
    metric,
    wind_reference,
    h09_within_estimate =
      estimate,
    h09_within_ci_low =
      ci95_low,
    h09_within_ci_high =
      ci95_high,
    h09_within_direction =
      direction
  )
]

FE_KEY <- EFFECTS[
  fe_structure ==
    "monitor_year_FE" &
  wind_reference %in%
    c(
      "参考论文1米每秒",
      "样本中位数"
    ),
  .(
    branch,
    fuel,
    metric,
    wind_reference,
    fe_estimate =
      estimate,
    fe_ci_low =
      ci95_low,
    fe_ci_high =
      ci95_high,
    fe_direction =
      direction,
    fe_ci_status =
      ci_status
  )
]

COMPARE <- merge(
  H09_WITHIN_KEY,
  FE_KEY,
  by = c(
    "branch",
    "fuel",
    "metric",
    "wind_reference"
  ),
  all = TRUE,
  sort = FALSE
)

COMPARE[
  ,
  direction_compare :=
    fifelse(
      is.finite(
        h09_within_estimate
      ) &
      is.finite(
        fe_estimate
      ),
      fifelse(
        sign(
          h09_within_estimate
        ) ==
        sign(
          fe_estimate
        ),
        "方向一致",
        "方向改变"
      ),
      "无法比较"
    )
]

fwrite(
  COMPARE,
  file.path(
    OUT,
    "05_07H09分解模型与固定效应模型对照.csv"
  )
)

PRIMARY <- EFFECTS[
  fe_structure ==
    "monitor_year_FE" &
  wind_reference ==
    "样本中位数"
]

STRICT <- EFFECTS[
  fe_structure ==
    "monitor_year_plus_date_FE" &
  wind_reference ==
    "样本中位数"
]

fwrite(
  PRIMARY,
  file.path(
    OUT,
    "06_核心结果_监测站年份固定效应_样本中位风速.csv"
  )
)

fwrite(
  STRICT,
  file.path(
    OUT,
    "07_最严格结果_监测站年份加日期固定效应_样本中位风速.csv"
  )
)

N_FAIL <- STATUS[
  status ==
    "FAIL",
  .N
]

N_PASS <- STATUS[
  status ==
    "PASS",
  .N
]

N_SKIP <- STATUS[
  status ==
    "SKIP",
  .N
]

N_PRIMARY_PASS <- STATUS[
  fe_structure ==
    "monitor_year_FE" &
  status ==
    "PASS",
  .N
]

N_STRICT_PASS <- STATUS[
  fe_structure ==
    "monitor_year_plus_date_FE" &
  status ==
    "PASS",
  .N
]

HARD <- data.table(
  check = c(
    "daily_files_24",
    "daily_duplicate_keys_zero",
    "quasi_rows_426657",
    "quasi_monitors_133",
    "quasi_unmatched_zero",
    "filter_all_rows_174200",
    "filter_formal_rows_107407",
    "filter_monitors_133",
    "filter_unmatched_zero",
    "total_FE_specs_48",
    "FE_FAIL_zero",
    "status_partition_48",
    "monitor_year_FE_results_exist",
    "monitor_year_plus_date_FE_results_exist",
    "continuous_curve_results_exist"
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
      N_Q
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
      FP_UNMATCHED
    ),
    as.character(
      nrow(STATUS)
    ),
    as.character(
      N_FAIL
    ),
    as.character(
      N_PASS +
      N_SKIP +
      N_FAIL
    ),
    as.character(
      as.integer(
        N_PRIMARY_PASS >
          0L
      )
    ),
    as.character(
      as.integer(
        N_STRICT_PASS >
          0L
      )
    ),
    as.character(
      as.integer(
        nrow(CURVES) >
          0L
      )
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
    "48",
    "0",
    "48",
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
    OUT,
    "08_07H09B_hard_checks.csv"
  )
)

cat(
  "\n============================================================\n",
  "A. 固定效应模型状态\n",
  "============================================================\n",
  sep = ""
)

print(
  STATUS[
    ,
    .N,
    by = .(
      fe_structure,
      status
    )
  ]
)

cat(
  "\n============================================================\n",
  "B. 07H09站点年份内分解 vs 监测站年份固定效应\n",
  "============================================================\n",
  sep = ""
)

print(
  COMPARE[
    wind_reference ==
      "样本中位数",
    .(
      branch,
      fuel,
      metric,
      h09_within_estimate,
      fe_estimate,
      h09_within_direction,
      fe_direction,
      direction_compare
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "C. 核心验证：监测站×年份固定效应，样本中位风速\n",
  "============================================================\n",
  sep = ""
)

print(
  PRIMARY[
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
  "D. 最严格验证：监测站×年份 + 日期固定效应，样本中位风速\n",
  "============================================================\n",
  sep = ""
)

print(
  STRICT[
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
  "E. 固定效应模型连续风速过零点\n",
  "============================================================\n",
  sep = ""
)

print(
  TRANS[
    ,
    .(
      branch,
      fuel,
      metric,
      fe_structure,
      exposure_beta,
      exposure_wind_interaction_beta,
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
  "# 07H09B 监测站年份固定效应严格验证",
  "",
  paste0(
    "- 生成时间：",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  "",
  "## 目的",
  "",
  "- 固定20 km、90°上风向。",
  "- 以监测站×年份固定效应直接控制每个监测站每一年的独立基线。",
  "- 因而发电站效应主要由同一监测站同一年内部的暴露变化识别。",
  "- 进一步加入日期固定效应，吸收全国共同日期冲击。",
  "- 风速始终为连续变量，并保留暴露×风速交互。",
  "- 连续效应曲线覆盖各模型实际风速第1至第99百分位。",
  "- 标准误按RadNet监测站聚类。",
  "",
  "## 结果状态",
  "",
  paste0(
    "- 固定效应模型PASS：",
    N_PASS,
    "/48"
  ),
  paste0(
    "- 固定效应模型SKIP：",
    N_SKIP,
    "/48"
  ),
  paste0(
    "- 固定效应模型FAIL：",
    N_FAIL,
    "/48"
  ),
  "",
  if (HARD_FAIL == 0L) {
    "**07H09B STRICT FE VALIDATION: PASS**"
  } else {
    "**07H09B STRICT FE VALIDATION: FAIL**"
  }
)

writeLines(
  REPORT,
  file.path(
    OUT,
    "09_07H09B_严格验证说明.md"
  ),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUT,
    "10_sessionInfo.txt"
  )
)

cat(
  "\n============================================================\n"
)

if (HARD_FAIL == 0L) {

  cat(
    "07H09B STRICT FE VALIDATION: PASS\n"
  )

} else {

  cat(
    "07H09B STRICT FE VALIDATION: FAIL\n"
  )
}

cat(
  "============================================================\n",
  "PASS = ",
  N_PASS,
  "/48\n",
  "SKIP = ",
  N_SKIP,
  "/48\n",
  "FAIL = ",
  N_FAIL,
  "/48\n",
  "Hard FAIL = ",
  HARD_FAIL,
  "\n",
  sep = ""
)

if (HARD_FAIL > 0L) {
  stopf(
    "07H09B存在硬性检查失败，先不要解释结果。"
  )
}
