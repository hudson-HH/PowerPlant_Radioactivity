#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999,
  width = 220
)

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(ggplot2)
})

# ============================================================
# 07H15
#
# 两种站点内比较：
#
# A. monitor FE
#    同一monitor整个2001-2024内部比较
#
# B. monitor × year FE
#    同一monitor同一年内部比较
#
# 共同条件：
# - filter-period
# - actual daily Count / Heat
# - 20-50 km
# - 360° vs 90° upwind
# - same complete-case sample
# - Exposure × wind
# - SE clustered by monitor
# - 1 m/s / Q25 / median / mean / Q75
# - continuous Q01-Q99 wind curve
#
# 不修改任何冻结上游结果。
# ============================================================


# ============================================================
# 1. Paths
# ============================================================

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

UP_DAILY_DIR <- file.path(
  HROOT,
  "fuel_daily_yearly"
)

ND_DAILY_DIR <- file.path(
  HROOT,
  "07H14_模型阶梯诊断_论文式随机截距到站点年份固定效应",
  "non_directional_daily_yearly"
)

OUT <- file.path(
  HROOT,
  "07H15A_扩展燃料两种站点固定效应与多风速Figure4_5诊断"
)

FIG_DIR <- file.path(
  OUT,
  "figures"
)

AUDIT_DIR <- file.path(
  OUT,
  "audit"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  AUDIT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Parameters
# ============================================================

YEARS <- 2001:2024

FUELS <- c(
  "biomass",
  "coal",
  "gas",
  "oil",
  "other_fossil",
  "solid_waste"
)

METRICS <- c(
  "count",
  "heat"
)

RADII <- c(
  20L,
  25L,
  30L,
  35L,
  40L,
  45L,
  50L
)

HEAT_SCALE <- 10000

EXPECTED_FP_ROWS <- 107407L
EXPECTED_FP_MONITORS <- 133L

FUEL_LABEL <- c(
  biomass = "Biomass",
  coal = "Coal",
  gas = "Natural gas",
  oil = "Oil",
  other_fossil = "Other fossil",
  solid_waste = "Solid waste"
)

FE_LABEL <- c(
  monitor_FE =
    "Monitor FE: within monitor across 2001–2024",

  monitor_year_FE =
    "Monitor × year FE: within monitor and year"
)


# ============================================================
# 3. Helpers
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

clean_chr <- function(x) {

  y <- trimws(
    as.character(x)
  )

  y[
    y %in% c(
      "",
      "NA",
      "N/A",
      "NULL",
      "."
    )
  ] <- NA_character_

  sub(
    "\\.0+$",
    "",
    y
  )
}

as_flag <- function(x) {

  if (is.logical(x)) {
    y <- x

  } else if (
    is.numeric(x) ||
    is.integer(x)
  ) {
    y <- x != 0

  } else {

    y <- tolower(
      trimws(
        as.character(x)
      )
    ) %in% c(
      "true",
      "t",
      "1",
      "yes",
      "y"
    )
  }

  y[
    is.na(y)
  ] <- FALSE

  y
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

nd_col <- function(
  fuel,
  metric,
  radius
) {

  paste0(
    "nd_",
    fuel,
    "_",
    metric,
    "_",
    radius,
    "km"
  )
}

up_col <- function(
  fuel,
  metric,
  radius
) {

  if (
    metric == "count"
  ) {

    paste0(
      "campd_",
      fuel,
      "_upwind_n_plants_",
      radius,
      "km_90deg"
    )

  } else {

    paste0(
      "campd_",
      fuel,
      "_upwind_heat_observed_",
      radius,
      "km_90deg"
    )
  }
}

metric_scale <- function(metric) {

  if (
    metric == "heat"
  ) {
    HEAT_SCALE
  } else {
    1
  }
}

metric_unit <- function(metric) {

  if (
    metric == "heat"
  ) {
    "per 10,000 MMBtu/day"
  } else {
    "per additional active plant"
  }
}

find_interaction <- function(coef_names) {

  hit <- intersect(
    c(
      "exposure:vel",
      "vel:exposure"
    ),
    coef_names
  )

  if (
    length(hit) != 1L
  ) {
    return(NA_character_)
  }

  hit[1]
}

linear_effect <- function(
  beta,
  V,
  wind_ms
) {

  int_name <- find_interaction(
    names(beta)
  )

  if (
    !"exposure" %in% names(beta) ||
    is.na(int_name)
  ) {

    return(
      c(
        estimate = NA_real_,
        se = NA_real_,
        ci95_low = NA_real_,
        ci95_high = NA_real_,
        p_wald = NA_real_
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
    "exposure"
  ] <- 1

  cvec[
    int_name
  ] <- wind_ms

  est <- sum(
    cvec *
      beta
  )

  vv <- as.numeric(
    t(cvec) %*%
      V %*%
      cvec
  )

  se <- sqrt(
    max(
      vv,
      0
    )
  )

  low <- est -
    1.96 *
    se

  high <- est +
    1.96 *
    se

  p <- if (
    is.finite(se) &&
    se > 0
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
    estimate = est,
    se = se,
    ci95_low = low,
    ci95_high = high,
    p_wald = p
  )
}


# ============================================================
# 4. Input checks
# ============================================================

cat(
  "============================================================\n",
  "07H15 INPUT AUDIT\n",
  "============================================================\n",
  sep = ""
)

if (
  !file.exists(FP_FILE)
) {
  stopf(
    "缺少07F filter-period文件：",
    FP_FILE
  )
}

for (
  yy in YEARS
) {

  qf <- file.path(
    QDIR,
    paste0(
      "07F_quasi_daily_analysis_base_",
      yy,
      ".csv.gz"
    )
  )

  uf <- file.path(
    UP_DAILY_DIR,
    paste0(
      "07H_fuel_daily_exposure_",
      yy,
      ".csv.gz"
    )
  )

  nf <- file.path(
    ND_DAILY_DIR,
    paste0(
      "07H14_nondirectional_fuel_daily_exposure_",
      yy,
      ".csv.gz"
    )
  )

  if (!file.exists(qf)) {
    stopf(
      "缺少07F quasi文件：",
      yy
    )
  }

  if (!file.exists(uf)) {
    stopf(
      "缺少07H01 upwind文件：",
      yy
    )
  }

  if (!file.exists(nf)) {
    stopf(
      "缺少07H14 360°文件：",
      yy
    )
  }
}

cat(
  "07F：PASS\n",
  "07H01 upwind daily：24/24 PASS\n",
  "07H14 360° daily：24/24 PASS\n",
  sep = ""
)


# ============================================================
# 5. Get formal eligible monitor set from 07F quasi branch
# ============================================================

ELIGIBLE_LIST <- vector(
  "list",
  length(YEARS)
)

for (
  ii in seq_along(YEARS)
) {

  yy <- YEARS[ii]

  qf <- file.path(
    QDIR,
    paste0(
      "07F_quasi_daily_analysis_base_",
      yy,
      ".csv.gz"
    )
  )

  z <- fread(
    qf,
    select = c(
      "loc_num",
      "flag_07f_primary_quasi_sample"
    ),
    showProgress = FALSE
  )

  z[
    ,
    flag_07f_primary_quasi_sample :=
      as_flag(
        flag_07f_primary_quasi_sample
      )
  ]

  ELIGIBLE_LIST[[ii]] <-
    unique(
      z[
        flag_07f_primary_quasi_sample %in%
          TRUE,
        .(
          loc_num =
            clean_chr(
              loc_num
            )
        )
      ]
    )
}

ELIGIBLE_LOCS <- sort(
  unique(
    rbindlist(
      ELIGIBLE_LIST
    )$loc_num
  )
)

rm(
  ELIGIBLE_LIST
)

cat(
  "正式eligible monitors = ",
  length(
    ELIGIBLE_LOCS
  ),
  "\n",
  sep = ""
)


# ============================================================
# 6. Read formal filter-period base
# ============================================================

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
    FP_FILE,
    nrows = 0,
    showProgress = FALSE
  )
)

MISS_FP <- setdiff(
  FP_REQUIRED,
  HFP
)

if (
  length(MISS_FP) > 0L
) {

  stopf(
    "filter-period缺字段：",
    paste(
      MISS_FP,
      collapse = ", "
    )
  )
}

FP <- fread(
  FP_FILE,
  select = FP_REQUIRED,
  showProgress = TRUE
)

FP[
  ,
  `:=`(
    loc_num =
      clean_chr(
        loc_num
      ),

    collect_start =
      as.IDate(
        collect_start
      ),

    collect_end =
      as.IDate(
        collect_end
      ),

    flag_conus =
      as_flag(
        flag_conus
      ),

    gross_beta_bq_m3 =
      safe_num(
        gross_beta_bq_m3
      ),

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

    pm_filter_period_official =
      safe_num(
        pm_filter_period_official
      ),

    inv_hpbl =
      safe_num(
        inv_hpbl
      ),

    hpbl =
      safe_num(
        hpbl
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
      273.15,

    fpmean_campd_upwind_heat_observed_20km_90deg =
      safe_num(
        fpmean_campd_upwind_heat_observed_20km_90deg
      )
  )
]

FP[
  ,
  flag_formal :=
    flag_conus %in% TRUE &

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

FP[
  ,
  `:=`(
    pm_model_centered =
      pm_filter_period_official -
      10,

    calendar_year =
      as.integer(
        format(
          collect_start,
          "%Y"
        )
      )
  )
]

FP[
  ,
  monitor_year_fe :=
    interaction(
      loc_num,
      calendar_year,
      drop = TRUE,
      lex.order = TRUE
    )
]

if (
  nrow(FP) !=
    EXPECTED_FP_ROWS
) {

  stopf(
    "filter-period正式样本行数=",
    nrow(FP),
    "；预期=",
    EXPECTED_FP_ROWS
  )
}

if (
  uniqueN(
    FP$loc_num
  ) !=
    EXPECTED_FP_MONITORS
) {

  stopf(
    "filter-period monitor数错误。"
  )
}

cat(
  "filter-period = ",
  format(
    nrow(FP),
    big.mark = ","
  ),
  " rows / ",
  uniqueN(
    FP$loc_num
  ),
  " monitors\n",
  sep = ""
)


# ============================================================
# 7. Expand filter-period into actual dates
# ============================================================

FP_DAYS <- FP[
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
    loc_num,
    n_period_days
  )
]

FP_DAYS[
  ,
  date :=
    as.IDate(
      date
    )
]


# ============================================================
# 8. Model formulas
# ============================================================

# ------------------------------------------------------------
# A. monitor FE
#
# Uses within-monitor variation over full 2001-2024.
#
# year_model + year_model^2 retained to control secular trend.
#
# surface U and latitude main effects are time-invariant
# within monitor and therefore absorbed by monitor FE.
# ------------------------------------------------------------

FORM_MONITOR <- as.formula(
  paste0(
    "gross_beta_mBq_m3 ~ ",
    "exposure + ",
    "exposure:vel + ",
    "c_prop + ",
    "sunspot_monthly + ",
    "pm_model_centered + ",
    "inv_hpbl + ",
    "rhum + ",
    "soilm + ",
    "vel + ",
    "year_model + I(year_model^2) + ",
    "temp_centered + ",
    "I(temp_centered^2) + ",
    "I(temp_centered^3) + ",
    "sunspot_monthly:lat_diff_35 + ",
    "surface_u_ppm_article:soilm ",
    "| loc_num"
  )
)


# ------------------------------------------------------------
# B. monitor × year FE
#
# Uses only within-monitor-year variation.
#
# Year main trends are absorbed by monitor×year FE.
# ------------------------------------------------------------

FORM_MONITOR_YEAR <- as.formula(
  paste0(
    "gross_beta_mBq_m3 ~ ",
    "exposure + ",
    "exposure:vel + ",
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
    "surface_u_ppm_article:soilm ",
    "| monitor_year_fe"
  )
)


# ============================================================
# 9. Fit helper
# ============================================================

fit_one <- function(
  d,
  exposure_name,
  fe_structure
) {

  z <- copy(d)

  z[
    ,
    exposure :=
      get(
        exposure_name
      )
  ]

  if (
    nrow(z) < 1000L ||
    uniqueN(
      z$loc_num
    ) < 2L ||
    !is.finite(
      sd(
        z$exposure
      )
    ) ||
    sd(
      z$exposure
    ) <= 0
  ) {

    return(
      list(
        status = "SKIP",
        note = "Insufficient exposure variation",
        fit = NULL
      )
    )
  }

  FORM <- if (
    fe_structure ==
      "monitor_FE"
  ) {
    FORM_MONITOR
  } else {
    FORM_MONITOR_YEAR
  }

  fit <- tryCatch(

    feols(
      FORM,
      data = z,
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
        status = "FAIL",
        note = fit$message,
        fit = NULL
      )
    )
  }

  beta <- coef(
    fit
  )

  V <- as.matrix(
    vcov(
      fit
    )
  )

  int_name <- find_interaction(
    names(beta)
  )

  if (
    !"exposure" %in%
      names(beta) ||
    is.na(int_name)
  ) {

    return(
      list(
        status = "SKIP",
        note = "Exposure main/interaction unavailable",
        fit = NULL
      )
    )
  }

  list(
    status = "PASS",
    note = "",
    fit = fit,
    beta = beta,
    V = V
  )
}


# ============================================================
# 10. Results containers
# ============================================================

STATUS_LIST <- list()
EFFECT_LIST <- list()
CURVE_LIST <- list()
ZERO_LIST <- list()

k_status <- 0L
k_effect <- 0L
k_curve <- 0L
k_zero <- 0L


# ============================================================
# 11. Main loops
# ============================================================

for (
  fuel in FUELS
) {

  for (
    metric in METRICS
  ) {

    cat(
      "\n============================================================\n",
      fuel,
      " / ",
      metric,
      "\n",
      "============================================================\n",
      sep = ""
    )

    nd_names <- vapply(
      RADII,
      function(r) {
        nd_col(
          fuel,
          metric,
          r
        )
      },
      character(1)
    )

    up_names <- vapply(
      RADII,
      function(r) {
        up_col(
          fuel,
          metric,
          r
        )
      },
      character(1)
    )

    DAILY_LIST <- vector(
      "list",
      length(YEARS)
    )

    for (
      ii in seq_along(YEARS)
    ) {

      yy <- YEARS[ii]

      nd_file <- file.path(
        ND_DAILY_DIR,
        paste0(
          "07H14_nondirectional_fuel_daily_exposure_",
          yy,
          ".csv.gz"
        )
      )

      up_file <- file.path(
        UP_DAILY_DIR,
        paste0(
          "07H_fuel_daily_exposure_",
          yy,
          ".csv.gz"
        )
      )

      HN <- names(
        fread(
          nd_file,
          nrows = 0,
          showProgress = FALSE
        )
      )

      HU <- names(
        fread(
          up_file,
          nrows = 0,
          showProgress = FALSE
        )
      )

      miss_nd <- setdiff(
        c(
          "loc_num",
          "date",
          nd_names
        ),
        HN
      )

      miss_up <- setdiff(
        c(
          "loc_num",
          "date",
          up_names
        ),
        HU
      )

      if (
        length(miss_nd) > 0L
      ) {
        stopf(
          yy,
          "年360°文件缺字段：",
          paste(
            miss_nd,
            collapse = ", "
          )
        )
      }

      if (
        length(miss_up) > 0L
      ) {
        stopf(
          yy,
          "年upwind文件缺字段：",
          paste(
            miss_up,
            collapse = ", "
          )
        )
      }

      A <- fread(
        nd_file,
        select = c(
          "loc_num",
          "date",
          nd_names
        ),
        showProgress = FALSE
      )

      B <- fread(
        up_file,
        select = c(
          "loc_num",
          "date",
          up_names
        ),
        showProgress = FALSE
      )

      A[
        ,
        `:=`(
          loc_num =
            clean_chr(
              loc_num
            ),

          date =
            as.IDate(
              date
            )
        )
      ]

      B[
        ,
        `:=`(
          loc_num =
            clean_chr(
              loc_num
            ),

          date =
            as.IDate(
              date
            )
        )
      ]

      DAILY_LIST[[ii]] <- merge(
        A,
        B,
        by = c(
          "loc_num",
          "date"
        ),
        all.x = TRUE,
        sort = FALSE
      )
    }

    DAILY <- rbindlist(
      DAILY_LIST,
      use.names = TRUE,
      fill = TRUE
    )

    rm(
      DAILY_LIST
    )

    invisible(
      gc()
    )


    # --------------------------------------------------------
    # Strict-complete filter-period means
    # --------------------------------------------------------

    exposure_cols <- c(
      nd_names,
      up_names
    )

    J <- DAILY[
      FP_DAYS,
      on = .(
        loc_num,
        date
      )
    ]

    if (
      J[
        ,
        .N,
        by = .(
          period_id,
          n_period_days
        )
      ][
        N != n_period_days,
        .N
      ] >
      0L
    ) {

      stopf(
        "period-day数量错误：",
        fuel,
        " / ",
        metric
      )
    }

    FP_EXP <- J[
      ,
      lapply(
        .SD,
        strict_mean
      ),
      by = period_id,
      .SDcols = exposure_cols
    ]

    X <- merge(
      FP,
      FP_EXP,
      by = "period_id",
      all.x = TRUE,
      sort = FALSE
    )

    rm(
      J,
      FP_EXP,
      DAILY
    )

    invisible(
      gc()
    )


    # --------------------------------------------------------
    # Each radius
    # --------------------------------------------------------

    for (
      radius in RADII
    ) {

      nd_name <- nd_col(
        fuel,
        metric,
        radius
      )

      up_name <- up_col(
        fuel,
        metric,
        radius
      )

      sc <- metric_scale(
        metric
      )

      vars <- c(
        "gross_beta_mBq_m3",
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
        "monitor_year_fe",
        nd_name,
        up_name
      )

      d <- X[
        ,
        ..vars
      ]

      setnames(
        d,
        c(
          nd_name,
          up_name
        ),
        c(
          "exposure_nd_raw",
          "exposure_up_raw"
        )
      )

      d[
        ,
        `:=`(
          exposure_nd =
            safe_num(
              exposure_nd_raw
            ) /
            sc,

          exposure_up =
            safe_num(
              exposure_up_raw
            ) /
            sc
        )
      ]

      complete_vars <- c(
        "gross_beta_mBq_m3",
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
        "monitor_year_fe",
        "exposure_nd",
        "exposure_up"
      )

      d <- d[
        complete.cases(
          d[
            ,
            ..complete_vars
          ]
        )
      ]

      if (
        nrow(d) <
          1000L
      ) {

        stopf(
          "共同样本过小：",
          fuel,
          " / ",
          metric,
          " / ",
          radius,
          " km"
        )
      }


      # ------------------------------------------------------
      # Prespecified wind references
      # ------------------------------------------------------

      q25 <- unname(
        quantile(
          d$vel,
          0.25,
          na.rm = TRUE
        )
      )

      med <- unname(
        quantile(
          d$vel,
          0.50,
          na.rm = TRUE
        )
      )

      mn <- mean(
        d$vel,
        na.rm = TRUE
      )

      q75 <- unname(
        quantile(
          d$vel,
          0.75,
          na.rm = TRUE
        )
      )

      q01 <- unname(
        quantile(
          d$vel,
          0.01,
          na.rm = TRUE
        )
      )

      q99 <- unname(
        quantile(
          d$vel,
          0.99,
          na.rm = TRUE
        )
      )

      WIND_REFS <- data.table(
        wind_reference = c(
          "paper_1m_s",
          "sample_q25",
          "sample_median",
          "sample_mean",
          "sample_q75"
        ),

        wind_ms = c(
          1,
          q25,
          med,
          mn,
          q75
        )
      )

      WIND_GRID <- seq(
        q01,
        q99,
        length.out = 101
      )


      # ------------------------------------------------------
      # 2 FE structures × 2 geometries
      # ------------------------------------------------------

      for (
        fe_structure in c(
          "monitor_FE",
          "monitor_year_FE"
        )
      ) {

        for (
          geometry in c(
            "non_directional_360",
            "upwind_90"
          )
        ) {

          exposure_name <- if (
            geometry ==
              "non_directional_360"
          ) {
            "exposure_nd"
          } else {
            "exposure_up"
          }

          cat(
            "[",
            fe_structure,
            "] ",
            geometry,
            " | ",
            fuel,
            " | ",
            metric,
            " | ",
            radius,
            " km | N=",
            nrow(d),
            "\n",
            sep = ""
          )

          rr <- fit_one(
            d,
            exposure_name,
            fe_structure
          )

          k_status <- k_status +
            1L

          STATUS_LIST[[k_status]] <-
            data.table(
              fe_structure =
                fe_structure,

              fe_label =
                unname(
                  FE_LABEL[
                    fe_structure
                  ]
                ),

              geometry =
                geometry,

              fuel =
                fuel,

              metric =
                metric,

              radius_km =
                radius,

              status =
                rr$status,

              n_rows =
                nrow(d),

              n_monitors =
                uniqueN(
                  d$loc_num
                ),

              n_monitor_years =
                uniqueN(
                  d$monitor_year_fe
                ),

              wind_q01 =
                q01,

              wind_q25 =
                q25,

              wind_median =
                med,

              wind_mean =
                mn,

              wind_q75 =
                q75,

              wind_q99 =
                q99,

              note =
                rr$note
            )

          if (
            rr$status !=
              "PASS"
          ) {
            next
          }


          # --------------------------------------------------
          # Representative wind effects
          # --------------------------------------------------

          for (
            jj in seq_len(
              nrow(
                WIND_REFS
              )
            )
          ) {

            w <- WIND_REFS$wind_ms[
              jj
            ]

            ee <- linear_effect(
              rr$beta,
              rr$V,
              w
            )

            k_effect <- k_effect +
              1L

            EFFECT_LIST[[k_effect]] <-
              data.table(
                fe_structure =
                  fe_structure,

                fe_label =
                  unname(
                    FE_LABEL[
                      fe_structure
                    ]
                  ),

                geometry =
                  geometry,

                fuel =
                  fuel,

                metric =
                  metric,

                radius_km =
                  radius,

                wind_reference =
                  WIND_REFS$wind_reference[
                    jj
                  ],

                wind_ms =
                  w,

                estimate =
                  unname(
                    ee[
                      "estimate"
                    ]
                  ),

                se =
                  unname(
                    ee[
                      "se"
                    ]
                  ),

                ci95_low =
                  unname(
                    ee[
                      "ci95_low"
                    ]
                  ),

                ci95_high =
                  unname(
                    ee[
                      "ci95_high"
                    ]
                  ),

                p_wald =
                  unname(
                    ee[
                      "p_wald"
                    ]
                  ),

                direction =
                  fifelse(
                    ee[
                      "estimate"
                    ] >
                      0,
                    "positive",
                    "negative"
                  ),

                ci_status =
                  fifelse(
                    ee[
                      "ci95_low"
                    ] >
                      0,
                    "CI_above_0",
                    fifelse(
                      ee[
                        "ci95_high"
                      ] <
                        0,
                      "CI_below_0",
                      "CI_crosses_0"
                    )
                  ),

                unit =
                  metric_unit(
                    metric
                  )
              )
          }


          # --------------------------------------------------
          # Continuous wind curve
          # --------------------------------------------------

          CURVE_ONE <- rbindlist(
            lapply(
              WIND_GRID,
              function(w) {

                ee <- linear_effect(
                  rr$beta,
                  rr$V,
                  w
                )

                data.table(
                  fe_structure =
                    fe_structure,

                  geometry =
                    geometry,

                  fuel =
                    fuel,

                  metric =
                    metric,

                  radius_km =
                    radius,

                  wind_ms =
                    w,

                  estimate =
                    unname(
                      ee[
                        "estimate"
                      ]
                    ),

                  ci95_low =
                    unname(
                      ee[
                        "ci95_low"
                      ]
                    ),

                  ci95_high =
                    unname(
                      ee[
                        "ci95_high"
                      ]
                    )
                )
              }
            )
          )

          k_curve <- k_curve +
            1L

          CURVE_LIST[[k_curve]] <-
            CURVE_ONE


          # --------------------------------------------------
          # Zero-crossing wind
          # --------------------------------------------------

          int_name <- find_interaction(
            names(
              rr$beta
            )
          )

          b0 <- unname(
            rr$beta[
              "exposure"
            ]
          )

          b1 <- unname(
            rr$beta[
              int_name
            ]
          )

          zero_w <- if (
            is.finite(b1) &&
            abs(b1) >
              1e-12
          ) {

            -b0 /
              b1

          } else {

            NA_real_
          }

          k_zero <- k_zero +
            1L

          ZERO_LIST[[k_zero]] <-
            data.table(
              fe_structure =
                fe_structure,

              geometry =
                geometry,

              fuel =
                fuel,

              metric =
                metric,

              radius_km =
                radius,

              exposure_beta =
                b0,

              exposure_wind_beta =
                b1,

              zero_crossing_wind_ms =
                zero_w,

              zero_inside_q01_q99 =
                is.finite(
                  zero_w
                ) &&
                zero_w >=
                  q01 &&
                zero_w <=
                  q99
            )
        }
      }

      rm(
        d
      )

      invisible(
        gc()
      )
    }

    rm(
      X
    )

    invisible(
      gc()
    )
  }
}


# ============================================================
# 12. Combine and save
# ============================================================

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

ZERO <- rbindlist(
  ZERO_LIST,
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  STATUS,
  file.path(
    OUT,
    "01_model_status.csv"
  )
)

fwrite(
  EFFECTS,
  file.path(
    OUT,
    "02_representative_wind_effects.csv"
  )
)

fwrite(
  CURVES,
  file.path(
    OUT,
    "03_continuous_wind_effect_curves.csv"
  )
)

fwrite(
  ZERO,
  file.path(
    OUT,
    "04_zero_crossing_wind.csv"
  )
)


# ============================================================
# 13. Side-by-side tables
# ============================================================

MEDIAN_EFFECTS <- EFFECTS[
  wind_reference ==
    "sample_median"
]

SIDE_FE <- dcast(
  MEDIAN_EFFECTS,
  geometry +
    fuel +
    metric +
    radius_km ~
    fe_structure,
  value.var = c(
    "estimate",
    "ci95_low",
    "ci95_high"
  )
)

fwrite(
  SIDE_FE,
  file.path(
    OUT,
    "05_monitorFE_vs_monitorYearFE_at_median_wind.csv"
  )
)

SIDE_DIRECTION <- dcast(
  MEDIAN_EFFECTS,
  fe_structure +
    fuel +
    metric +
    radius_km ~
    geometry,
  value.var = c(
    "estimate",
    "ci95_low",
    "ci95_high"
  )
)

fwrite(
  SIDE_DIRECTION,
  file.path(
    OUT,
    "06_360_vs_upwind_at_median_wind.csv"
  )
)


# ============================================================
# 14. Direction summaries
# ============================================================

SUMMARY <- EFFECTS[
  ,
  .(
    n_radii =
      .N,

    n_positive =
      sum(
        estimate >
          0
      ),

    n_negative =
      sum(
        estimate <
          0
      ),

    n_ci_above_0 =
      sum(
        ci95_low >
          0
      ),

    n_ci_below_0 =
      sum(
        ci95_high <
          0
      ),

    n_ci_cross_0 =
      sum(
        ci95_low <=
          0 &
        ci95_high >=
          0
      )
  ),
  by = .(
    fe_structure,
    geometry,
    fuel,
    metric,
    wind_reference
  )
]

fwrite(
  SUMMARY,
  file.path(
    OUT,
    "07_direction_summary_20_50km.csv"
  )
)


# ============================================================
# 15. Plot helper
# ============================================================

plot_effect <- function(
  fe_value,
  metric_value,
  filename
) {

  D <- MEDIAN_EFFECTS[
    fe_structure ==
      fe_value &
    metric ==
      metric_value
  ]

  D[
    ,
    fuel_label :=
      factor(
        FUEL_LABEL[
          fuel
        ],
        levels =
          FUEL_LABEL[
            FUELS
          ]
      )
  ]

  D[
    ,
    geometry_label :=
      factor(
        fifelse(
          geometry ==
            "non_directional_360",
          "360° non-directional",
          "90° upwind"
        ),
        levels = c(
          "360° non-directional",
          "90° upwind"
        )
      )
  ]

  pd <- position_dodge(
    width = 1.4
  )

  p <- ggplot(
    D,
    aes(
      x = radius_km,
      y = estimate,
      group = geometry_label,
      shape = geometry_label,
      linetype = geometry_label
    )
  ) +

    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +

    geom_errorbar(
      aes(
        ymin = ci95_low,
        ymax = ci95_high
      ),
      width = 0.8,
      linewidth = 0.45,
      position = pd
    ) +

    geom_line(
      linewidth = 0.6,
      position = pd
    ) +

    geom_point(
      size = 1.8,
      position = pd
    ) +

    facet_wrap(
      ~fuel_label,
      scales = "free_y",
      ncol = 2
    ) +

    scale_x_continuous(
      breaks = RADII
    ) +

    labs(
      title =
        paste0(
          unname(
            FE_LABEL[
              fe_value
            ]
          ),
          " — ",
          ifelse(
            metric_value ==
              "count",
            "Active plant count",
            "Actual daily Heat Input"
          )
        ),

      subtitle =
        "Filter-period analysis; effects evaluated at the sample median wind speed; 360° and upwind use the same complete-case sample.",

      x =
        "Cumulative radius (km)",

      y =
        if (
          metric_value ==
            "count"
        ) {

          "Adjusted gross-beta change per additional active plant (mBq/m³)"

        } else {

          "Adjusted gross-beta change per 10,000 MMBtu/day (mBq/m³)"
        },

      shape =
        "Exposure geometry",

      linetype =
        "Exposure geometry"
    ) +

    theme_bw(
      base_size = 10
    ) +

    theme(
      panel.grid.minor =
        element_blank(),

      strip.background =
        element_rect(
          fill = "white"
        ),

      strip.text =
        element_text(
          face = "bold"
        ),

      plot.title =
        element_text(
          face = "bold"
        ),

      legend.position =
        "bottom"
    )

  ggsave(
    file.path(
      FIG_DIR,
      paste0(
        filename,
        ".png"
      )
    ),
    p,
    width = 11,
    height = 8,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  if (
    capabilities(
      "cairo"
    )
  ) {

    ggsave(
      file.path(
        FIG_DIR,
        paste0(
          filename,
          ".pdf"
        )
      ),
      p,
      width = 11,
      height = 8,
      units = "in",
      device = cairo_pdf
    )

  } else {

    ggsave(
      file.path(
        FIG_DIR,
        paste0(
          filename,
          ".pdf"
        )
      ),
      p,
      width = 11,
      height = 8,
      units = "in"
    )
  }
}


# ============================================================
# 16. Four principal figures
# ============================================================

plot_effect(
  "monitor_FE",
  "count",
  "01_monitorFE_Count_medianWind"
)

plot_effect(
  "monitor_FE",
  "heat",
  "02_monitorFE_Heat_medianWind"
)

plot_effect(
  "monitor_year_FE",
  "count",
  "03_monitorYearFE_Count_medianWind"
)

plot_effect(
  "monitor_year_FE",
  "heat",
  "04_monitorYearFE_Heat_medianWind"
)


# ============================================================
# 17. Hard checks
# ============================================================

EXPECTED_MODELS <-
  length(FUELS) *
  length(METRICS) *
  length(RADII) *
  2L *
  2L

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

EXPECTED_EFFECT_ROWS <-
  N_PASS *
  5L

CHECKS <- data.table(

  check = c(
    "formal_filter_rows",
    "formal_filter_monitors",
    "model_status_rows",
    "model_failures",
    "representative_effect_rows",
    "four_png",
    "four_pdf"
  ),

  observed = c(
    as.character(
      nrow(FP)
    ),

    as.character(
      uniqueN(
        FP$loc_num
      )
    ),

    as.character(
      nrow(STATUS)
    ),

    as.character(
      N_FAIL
    ),

    as.character(
      nrow(EFFECTS)
    ),

    as.character(
      length(
        list.files(
          FIG_DIR,
          pattern =
            "\\.png$"
        )
      )
    ),

    as.character(
      length(
        list.files(
          FIG_DIR,
          pattern =
            "\\.pdf$"
        )
      )
    )
  ),

  expected = c(
    as.character(
      EXPECTED_FP_ROWS
    ),

    as.character(
      EXPECTED_FP_MONITORS
    ),

    as.character(
      EXPECTED_MODELS
    ),

    "0",

    as.character(
      EXPECTED_EFFECT_ROWS
    ),

    "4",
    "4"
  )
)

CHECKS[
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
  CHECKS,
  file.path(
    AUDIT_DIR,
    "01_07H15_hard_checks.csv"
  )
)


# ============================================================
# 18. Terminal summaries
# ============================================================

cat(
  "\n============================================================\n",
  "07H15 MODEL STATUS\n",
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
  ][
    order(
      fe_structure,
      status
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "SAMPLE MEDIAN WIND SUMMARY\n",
  "============================================================\n",
  sep = ""
)

print(
  SUMMARY[
    wind_reference ==
      "sample_median",
    .(
      fe_structure,
      geometry,
      fuel,
      metric,
      n_radii,
      n_positive,
      n_negative,
      n_ci_above_0,
      n_ci_below_0,
      n_ci_cross_0
    )
  ][
    order(
      fe_structure,
      geometry,
      fuel,
      metric
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "1 m/s SUMMARY\n",
  "============================================================\n",
  sep = ""
)

print(
  SUMMARY[
    wind_reference ==
      "paper_1m_s",
    .(
      fe_structure,
      geometry,
      fuel,
      metric,
      n_radii,
      n_positive,
      n_negative,
      n_ci_above_0,
      n_ci_below_0,
      n_ci_cross_0
    )
  ][
    order(
      fe_structure,
      geometry,
      fuel,
      metric
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "07H15 HARD CHECKS\n",
  "============================================================\n",
  sep = ""
)

print(
  CHECKS,
  nrows = Inf
)

HARD_FAIL <- CHECKS[
  status ==
    "FAIL",
  .N
]

cat(
  "\n============================================================\n",
  ifelse(
    HARD_FAIL ==
      0L,
    "07H15 TWO-FE MULTI-WIND DIAGNOSIS: PASS\n",
    "07H15 TWO-FE MULTI-WIND DIAGNOSIS: FAIL\n"
  ),
  "Hard FAIL = ",
  HARD_FAIL,
  "\n",
  "Model FAIL = ",
  N_FAIL,
  "\n",
  "Model PASS = ",
  N_PASS,
  "\n",
  "Output = ",
  OUT,
  "\n",
  "============================================================\n",
  sep = ""
)

if (
  HARD_FAIL >
    0L
) {
  stopf(
    "07H15存在Hard FAIL。"
  )
}
