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
  "ggplot2",
  "patchwork"
)

for (p in required_packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(ggplot2)
  library(patchwork)
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


# ============================================================
# 07H16
#
# 最终正文 Figure 4 / Figure 5
#
# Figure 4
#   不同燃料类型的90°上风向实际运行暴露
#   与gross beta的关联
#
# Figure 5
#   上风向与下风向方向性直接比较
#   同一模型同时纳入upwind和downwind
#   正式绘制：
#   Upwind effect - Downwind effect
#
# 主文：
#   Monitor × year FE
#
# Supplementary：
#   Monitor FE + year + year²
#
# Fuel:
#   coal / gas / oil / biomass
#
# Metric:
#   Count
#   Heat Input
#
# Radius:
#   20–50 km, step 5 km
#
# Wind:
#   各模型complete-case样本的median wind
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

DAILY_DIR <- file.path(
  HROOT,
  "fuel_daily_yearly"
)

# ------------------------------------------------------------
# 已经验证过的07H15结果
# 用于交叉检查Figure 4
# ------------------------------------------------------------

H15_EFFECT_FILE <- file.path(
  HROOT,
  "07H15_两种站点固定效应与多风速Figure4_5诊断",
  "02_representative_wind_effects.csv"
)

# ------------------------------------------------------------
# 已经验证过的07H10方向直接比较
# 用于交叉检查Figure 5的Monitor×year FE
# ------------------------------------------------------------

H10_DIRECTION_FILE <- file.path(
  HROOT,
  "07H10_燃料特异站点年份内综合正式分析",
  "04_上风向下风向直接比较",
  "02_上风向下风向及直接差值检验.csv"
)

OUT <- file.path(
  HROOT,
  "07H16_Figure4_Figure5正文与MonitorFE补充图"
)

FIG_DIR <- file.path(
  OUT,
  "figures"
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


# ============================================================
# 2. Parameters
# ============================================================

YEARS <- 2001:2024

RADII <- c(
  20L,
  25L,
  30L,
  35L,
  40L,
  45L,
  50L
)

FUELS <- c(
  "coal",
  "gas",
  "oil",
  "biomass"
)

METRICS <- c(
  "count",
  "heat"
)

FUEL_LABEL <- c(
  coal = "Coal",
  gas = "Natural gas",
  oil = "Oil",
  biomass = "Biomass"
)

EXPECTED_FP_ROWS <- 107407L
EXPECTED_FP_MONITORS <- 133L

TOL <- 1e-7


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
    is.na(y) |
    y == "" |
    y %in% c(
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

as_bool <- function(x) {

  if (is.logical(x)) {
    ans <- x

  } else if (
    is.numeric(x) ||
    is.integer(x)
  ) {

    ans <- x != 0

  } else {

    z <- tolower(
      trimws(
        as.character(x)
      )
    )

    ans <- z %in% c(
      "true",
      "t",
      "1",
      "yes",
      "y"
    )
  }

  ans[is.na(ans)] <- FALSE

  ans
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

raw_col <- function(
  fuel,
  direction,
  metric,
  radius
) {

  if (
    metric == "count"
  ) {

    return(
      sprintf(
        "campd_%s_%s_n_plants_%dkm_90deg",
        fuel,
        direction,
        radius
      )
    )
  }

  sprintf(
    "campd_%s_%s_heat_observed_%dkm_90deg",
    fuel,
    direction,
    radius
  )
}

metric_scale <- function(metric) {

  if (
    metric == "heat"
  ) {
    return(
      10000
    )
  }

  1
}

find_interaction <- function(
  coef_names,
  a,
  b
) {

  z <- intersect(
    c(
      paste0(
        a,
        ":",
        b
      ),
      paste0(
        b,
        ":",
        a
      )
    ),
    coef_names
  )

  if (
    length(z) != 1L
  ) {
    return(
      NA_character_
    )
  }

  z[1]
}

effect_at_wind <- function(
  fit,
  main_term,
  interaction_term,
  wind_ms
) {

  b <- coef(
    fit
  )

  V <- as.matrix(
    vcov(
      fit
    )
  )

  if (
    !main_term %in%
      names(b) ||
    is.na(interaction_term) ||
    !interaction_term %in%
      names(b)
  ) {

    return(
      data.table(
        estimate = NA_real_,
        se = NA_real_,
        ci95_low = NA_real_,
        ci95_high = NA_real_,
        p_wald = NA_real_
      )
    )
  }

  est <-
    unname(
      b[main_term]
    ) +
    wind_ms *
    unname(
      b[interaction_term]
    )

  vv <-
    V[
      main_term,
      main_term
    ] +
    wind_ms^2 *
    V[
      interaction_term,
      interaction_term
    ] +
    2 *
    wind_ms *
    V[
      main_term,
      interaction_term
    ]

  se <- if (
    is.finite(vv)
  ) {

    sqrt(
      max(
        as.numeric(vv),
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
        abs(
          est /
          se
        ),
        lower.tail = FALSE
      )

  } else {

    NA_real_
  }

  data.table(
    estimate = est,
    se = se,
    ci95_low = low,
    ci95_high = high,
    p_wald = p
  )
}

linear_combination <- function(
  fit,
  terms,
  weights
) {

  b <- coef(
    fit
  )

  V <- as.matrix(
    vcov(
      fit
    )
  )

  if (
    any(
      !terms %in%
        names(b)
    )
  ) {

    return(
      data.table(
        estimate = NA_real_,
        se = NA_real_,
        ci95_low = NA_real_,
        ci95_high = NA_real_,
        p_wald = NA_real_
      )
    )
  }

  bb <- b[
    terms
  ]

  VV <- V[
    terms,
    terms,
    drop = FALSE
  ]

  est <- sum(
    weights *
    bb
  )

  vv <- as.numeric(
    t(weights) %*%
    VV %*%
    weights
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
        abs(
          est /
          se
        ),
        lower.tail = FALSE
      )

  } else {

    NA_real_
  }

  data.table(
    estimate = est,
    se = se,
    ci95_low = low,
    ci95_high = high,
    p_wald = p
  )
}


# ============================================================
# 4. Regression formula
#
# monitor_FE:
#   monitor FE
#   + linear/quadratic calendar-time trend
#
# monitor_year_FE:
#   monitor × calendar-year FE
# ============================================================

base_rhs <- function(
  fe_structure
) {

  rhs <- c(
    "c_prop",
    "sunspot_monthly",
    "pm_centered_10",
    "inv_hpbl",
    "rhum",
    "soilm",
    "vel",
    "temp_centered",
    "I(temp_centered^2)",
    "I(temp_centered^3)",
    "sunspot_monthly:lat_diff_35",
    "surface_u_ppm_article:soilm"
  )

  if (
    fe_structure ==
      "monitor_FE"
  ) {

    rhs <- c(
      rhs,
      "year_model",
      "I(year_model^2)"
    )
  }

  rhs
}

make_single_formula <- function(
  fe_structure
) {

  rhs <- c(
    "exposure",
    "exposure:vel",
    base_rhs(
      fe_structure
    )
  )

  fe_part <- if (
    fe_structure ==
      "monitor_FE"
  ) {
    "loc_num"
  } else {
    "monitor_year_fe"
  }

  as.formula(
    paste0(
      "gross_beta_mBq_m3 ~ ",
      paste(
        rhs,
        collapse = " + "
      ),
      " | ",
      fe_part
    )
  )
}

make_direction_formula <- function(
  fe_structure
) {

  rhs <- c(
    "up_exposure",
    "down_exposure",
    "up_exposure:vel",
    "down_exposure:vel",
    base_rhs(
      fe_structure
    )
  )

  fe_part <- if (
    fe_structure ==
      "monitor_FE"
  ) {
    "loc_num"
  } else {
    "monitor_year_fe"
  }

  as.formula(
    paste0(
      "gross_beta_mBq_m3 ~ ",
      paste(
        rhs,
        collapse = " + "
      ),
      " | ",
      fe_part
    )
  )
}


# ============================================================
# 5. Input audit
# ============================================================

cat(
  "============================================================\n",
  "07H16 INPUT AUDIT\n",
  "============================================================\n",
  sep = ""
)

REQ_FILES <- c(
  FP_FILE,
  H15_EFFECT_FILE,
  H10_DIRECTION_FILE
)

MISS_FILES <- REQ_FILES[
  !file.exists(
    REQ_FILES
  )
]

if (
  length(MISS_FILES) > 0L
) {
  stopf(
    "缺少输入：\n",
    paste(
      MISS_FILES,
      collapse = "\n"
    )
  )
}

DAILY_FILES <- file.path(
  DAILY_DIR,
  paste0(
    "07H_fuel_daily_exposure_",
    YEARS,
    ".csv.gz"
  )
)

if (
  sum(
    file.exists(
      DAILY_FILES
    )
  ) !=
    24L
) {
  stopf(
    "07H01 daily files不是24/24。"
  )
}

cat(
  "07F filter-period：PASS\n",
  "07H01 daily：24/24 PASS\n",
  "07H15 comparison source：PASS\n",
  "07H10 direction comparison source：PASS\n",
  sep = ""
)


# ============================================================
# 6. Recover formal 133-monitor universe
# ============================================================

Q_MON <- rbindlist(
  lapply(
    YEARS,
    function(yy) {

      f <- file.path(
        QDIR,
        paste0(
          "07F_quasi_daily_analysis_base_",
          yy,
          ".csv.gz"
        )
      )

      if (
        !file.exists(f)
      ) {
        stopf(
          "缺少07F quasi file：",
          f
        )
      }

      z <- fread(
        f,
        select = c(
          "loc_num",
          "flag_07f_primary_quasi_sample"
        ),
        showProgress = FALSE
      )

      z[
        ,
        loc_num :=
          clean_chr(
            loc_num
          )
      ]

      z[
        ,
        flag_07f_primary_quasi_sample :=
          as_bool(
            flag_07f_primary_quasi_sample
          )
      ]

      z
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

ELIGIBLE_LOCS <- sort(
  unique(
    Q_MON[
      flag_07f_primary_quasi_sample %in%
        TRUE,
      loc_num
    ]
  )
)

if (
  length(
    ELIGIBLE_LOCS
  ) !=
    EXPECTED_FP_MONITORS
) {
  stopf(
    "正式eligible monitors不是133；实际=",
    length(
      ELIGIBLE_LOCS
    )
  )
}

rm(Q_MON)
invisible(gc())

cat(
  "Formal eligible monitors = 133：PASS\n"
)


# ============================================================
# 7. Formal filter-period sample
# ============================================================

FP_COLS <- c(
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
    FP_FILE,
    nrows = 0,
    showProgress = FALSE
  )
)

FP_MISS <- setdiff(
  FP_COLS,
  FP_HEADER
)

if (
  length(FP_MISS) > 0L
) {
  stopf(
    "07F filter-period缺字段：",
    paste(
      FP_MISS,
      collapse = ", "
    )
  )
}

FP <- fread(
  FP_FILE,
  select = FP_COLS,
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
      as_bool(
        flag_conus
      )
  )
]

NUM_COLS <- c(
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
  "fpmean_campd_upwind_heat_observed_20km_90deg"
)

for (
  cc in NUM_COLS
) {

  set(
    FP,
    j = cc,
    value =
      safe_num(
        FP[[cc]]
      )
  )
}

# ------------------------------------------------------------
# 精确复现冻结07H10/07H15正式filter-period定义
# ------------------------------------------------------------

FP <- FP[
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

if (
  nrow(FP) !=
    EXPECTED_FP_ROWS
) {
  stopf(
    "正式filter-period不是107407；实际=",
    nrow(FP)
  )
}

if (
  uniqueN(
    FP$loc_num
  ) !=
    EXPECTED_FP_MONITORS
) {
  stopf(
    "正式filter-period不是133 monitors。"
  )
}

FP[
  ,
  `:=`(
    period_id = .I,

    period_days =
      as.integer(
        collect_end -
        collect_start
      ),

    pm_centered_10 =
      pm_filter_period_official -
      10,

    temp_centered =
      temp -
      273.15,

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
  FP[
    !is.finite(
      period_days
    ) |
    period_days <= 0,
    .N
  ] >
    0L
) {
  stopf(
    "存在非法filter period。"
  )
}

cat(
  "Formal filter-period = ",
  format(
    nrow(FP),
    big.mark = ","
  ),
  " / ",
  uniqueN(
    FP$loc_num
  ),
  " monitors：PASS\n",
  sep = ""
)


# ============================================================
# 8. Expand filter periods
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
    period_days
  )
]

FP_DAYS[
  ,
  date :=
    as.IDate(
      date
    )
]


MODEL_BASE_COLS <- c(
  "gross_beta_mBq_m3",
  "c_prop",
  "surface_u_ppm_article",
  "sunspot_monthly",
  "pm_centered_10",
  "inv_hpbl",
  "rhum",
  "soilm",
  "vel",
  "lat_diff_35",
  "year_model",
  "temp_centered",
  "loc_num",
  "monitor_year_fe"
)


# ============================================================
# 9. Result containers
# ============================================================

STATUS <- data.table()

UPWIND_EFFECTS <- data.table()

DIRECTION_EFFECTS <- data.table()

DAILY_AUDIT <- data.table()


# ============================================================
# 10. Radius loop
#     每次只读取当前半径，降低内存占用
# ============================================================

for (
  radius in RADII
) {

  cat(
    "\n============================================================\n",
    "RADIUS = ",
    radius,
    " km\n",
    "============================================================\n",
    sep = ""
  )

  EXPOSURE_COLS <- unlist(
    lapply(
      FUELS,
      function(fuel) {

        unlist(
          lapply(
            METRICS,
            function(metric) {

              c(
                raw_col(
                  fuel,
                  "upwind",
                  metric,
                  radius
                ),
                raw_col(
                  fuel,
                  "downwind",
                  metric,
                  radius
                )
              )
            }
          ),
          use.names = FALSE
        )
      }
    ),
    use.names = FALSE
  )

  # ----------------------------------------------------------
  # Read 24 years
  # ----------------------------------------------------------

  DAILY <- rbindlist(
    lapply(
      YEARS,
      function(yy) {

        f <- file.path(
          DAILY_DIR,
          paste0(
            "07H_fuel_daily_exposure_",
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

        req <- c(
          "loc_num",
          "date",
          EXPOSURE_COLS
        )

        miss <- setdiff(
          req,
          H
        )

        if (
          length(miss) > 0L
        ) {
          stopf(
            yy,
            "年07H01缺字段：",
            paste(
              miss,
              collapse = ", "
            )
          )
        }

        z <- fread(
          f,
          select = req,
          showProgress = FALSE
        )

        z[
          ,
          `:=`(
            loc_num =
              clean_chr(
                loc_num
              ),

            date =
              as.IDate(
                as.character(
                  date
                )
              ),

            daily_present =
              1L
          )
        ]

        for (
          cc in EXPOSURE_COLS
        ) {

          set(
            z,
            j = cc,
            value =
              safe_num(
                z[[cc]]
              )
          )
        }

        z
      }
    ),
    use.names = TRUE,
    fill = TRUE
  )

  if (
    DAILY[
      ,
      .N,
      by = .(
        loc_num,
        date
      )
    ][
      N > 1L,
      .N
    ] >
      0L
  ) {
    stopf(
      radius,
      " km daily exposure存在重复monitor-date。"
    )
  }

  # ----------------------------------------------------------
  # Join to formal filter-period dates
  # ----------------------------------------------------------

  J <- merge(
    FP_DAYS,
    DAILY,
    by = c(
      "loc_num",
      "date"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  unmatched <- J[
    is.na(
      daily_present
    ),
    .N
  ]

  DAILY_AUDIT <- rbind(
    DAILY_AUDIT,
    data.table(
      radius_km = radius,
      n_period_days = nrow(J),
      unmatched_daily_rows = unmatched
    ),
    fill = TRUE
  )

  if (
    unmatched != 0L
  ) {
    stopf(
      radius,
      " km存在",
      unmatched,
      "个formal period-day无法连接07H01 daily。"
    )
  }

  # ----------------------------------------------------------
  # Strict period mean
  # ----------------------------------------------------------

  FP_EXP <- J[
    ,
    lapply(
      .SD,
      strict_mean
    ),
    by = period_id,
    .SDcols =
      EXPOSURE_COLS
  ]

  D0 <- merge(
    FP,
    FP_EXP,
    by = "period_id",
    all.x = TRUE,
    sort = FALSE
  )

  if (
    nrow(D0) !=
      EXPECTED_FP_ROWS
  ) {
    stopf(
      radius,
      " km聚合后改变filter-period行数。"
    )
  }

  rm(
    DAILY,
    J,
    FP_EXP
  )

  invisible(
    gc()
  )


  # ==========================================================
  # Fuel × metric
  # ==========================================================

  for (
    fuel in FUELS
  ) {

    for (
      metric in METRICS
    ) {

      cat(
        fuel,
        " / ",
        metric,
        " / ",
        radius,
        " km\n",
        sep = ""
      )

      scale_value <- metric_scale(
        metric
      )

      up_col <- raw_col(
        fuel,
        "upwind",
        metric,
        radius
      )

      down_col <- raw_col(
        fuel,
        "downwind",
        metric,
        radius
      )


      # ======================================================
      # A. Figure 4 model
      #    Upwind single-exposure model
      # ======================================================

      DS <- D0[
        ,
        c(
          MODEL_BASE_COLS,
          up_col
        ),
        with = FALSE
      ]

      DS[
        ,
        exposure :=
          safe_num(
            get(
              up_col
            )
          ) /
          scale_value
      ]

      DS <- DS[
        complete.cases(
          DS[
            ,
            c(
              MODEL_BASE_COLS,
              "exposure"
            ),
            with = FALSE
          ]
        )
      ]

      median_wind_single <- median(
        DS$vel,
        na.rm = TRUE
      )


      # ======================================================
      # B. Figure 5 model
      #    Joint upwind + downwind
      # ======================================================

      DD <- D0[
        ,
        c(
          MODEL_BASE_COLS,
          up_col,
          down_col
        ),
        with = FALSE
      ]

      DD[
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

      DD <- DD[
        complete.cases(
          DD[
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

      median_wind_direction <- median(
        DD$vel,
        na.rm = TRUE
      )


      # ======================================================
      # FE structures
      # ======================================================

      for (
        fe_structure in c(
          "monitor_FE",
          "monitor_year_FE"
        )
      ) {

        # ----------------------------------------------------
        # Figure 4 single-upwind
        # ----------------------------------------------------

        FORM_SINGLE <- make_single_formula(
          fe_structure
        )

        fit_single <- tryCatch(
          feols(
            FORM_SINGLE,
            data = DS,
            vcov = ~loc_num,
            notes = FALSE,
            warn = FALSE
          ),
          error = function(e) {
            e
          }
        )

        single_status <- "PASS"
        single_note <- ""

        if (
          inherits(
            fit_single,
            "error"
          )
        ) {

          single_status <- "FAIL"

          single_note <- conditionMessage(
            fit_single
          )

        } else {

          bn <- names(
            coef(
              fit_single
            )
          )

          int_name <- find_interaction(
            bn,
            "exposure",
            "vel"
          )

          if (
            !"exposure" %in%
              bn ||
            is.na(
              int_name
            )
          ) {

            single_status <- "SKIP"

            single_note <-
              "exposure or exposure×wind dropped"
          }
        }

        STATUS <- rbind(
          STATUS,
          data.table(
            model_family =
              "Figure4_upwind_single",

            fe_structure =
              fe_structure,

            fuel =
              fuel,

            metric =
              metric,

            radius_km =
              radius,

            n_rows =
              nrow(
                DS
              ),

            n_monitors =
              uniqueN(
                DS$loc_num
              ),

            status =
              single_status,

            note =
              single_note
          ),
          fill = TRUE
        )

        if (
          single_status ==
            "PASS"
        ) {

          int_name <- find_interaction(
            names(
              coef(
                fit_single
              )
            ),
            "exposure",
            "vel"
          )

          ee <- effect_at_wind(
            fit_single,
            "exposure",
            int_name,
            median_wind_single
          )

          UPWIND_EFFECTS <- rbind(
            UPWIND_EFFECTS,
            data.table(
              fe_structure =
                fe_structure,

              fuel =
                fuel,

              metric =
                metric,

              radius_km =
                radius,

              wind_reference =
                "sample_median",

              wind_ms =
                median_wind_single,

              estimate =
                ee$estimate,

              se =
                ee$se,

              ci95_low =
                ee$ci95_low,

              ci95_high =
                ee$ci95_high,

              p_wald =
                ee$p_wald,

              n_rows =
                nrow(
                  DS
                )
            ),
            fill = TRUE
          )
        }


        # ----------------------------------------------------
        # Figure 5 joint direction model
        # ----------------------------------------------------

        FORM_DIRECTION <- make_direction_formula(
          fe_structure
        )

        fit_direction <- tryCatch(
          feols(
            FORM_DIRECTION,
            data = DD,
            vcov = ~loc_num,
            notes = FALSE,
            warn = FALSE
          ),
          error = function(e) {
            e
          }
        )

        direction_status <- "PASS"
        direction_note <- ""

        if (
          inherits(
            fit_direction,
            "error"
          )
        ) {

          direction_status <- "FAIL"

          direction_note <- conditionMessage(
            fit_direction
          )

        } else {

          bn <- names(
            coef(
              fit_direction
            )
          )

          up_int <- find_interaction(
            bn,
            "up_exposure",
            "vel"
          )

          down_int <- find_interaction(
            bn,
            "down_exposure",
            "vel"
          )

          if (
            !"up_exposure" %in%
              bn ||
            !"down_exposure" %in%
              bn ||
            is.na(
              up_int
            ) ||
            is.na(
              down_int
            )
          ) {

            direction_status <- "SKIP"

            direction_note <-
              "up/down exposure term or exposure×wind dropped"
          }
        }

        STATUS <- rbind(
          STATUS,
          data.table(
            model_family =
              "Figure5_direction_joint",

            fe_structure =
              fe_structure,

            fuel =
              fuel,

            metric =
              metric,

            radius_km =
              radius,

            n_rows =
              nrow(
                DD
              ),

            n_monitors =
              uniqueN(
                DD$loc_num
              ),

            status =
              direction_status,

            note =
              direction_note
          ),
          fill = TRUE
        )

        if (
          direction_status ==
            "PASS"
        ) {

          bn <- names(
            coef(
              fit_direction
            )
          )

          up_int <- find_interaction(
            bn,
            "up_exposure",
            "vel"
          )

          down_int <- find_interaction(
            bn,
            "down_exposure",
            "vel"
          )

          up_eff <- effect_at_wind(
            fit_direction,
            "up_exposure",
            up_int,
            median_wind_direction
          )

          down_eff <- effect_at_wind(
            fit_direction,
            "down_exposure",
            down_int,
            median_wind_direction
          )

          diff_eff <- linear_combination(
            fit_direction,
            c(
              "up_exposure",
              up_int,
              "down_exposure",
              down_int
            ),
            c(
              1,
              median_wind_direction,
              -1,
              -median_wind_direction
            )
          )

          DIRECTION_EFFECTS <- rbind(
            DIRECTION_EFFECTS,

            data.table(
              fe_structure =
                fe_structure,

              fuel =
                fuel,

              metric =
                metric,

              radius_km =
                radius,

              component =
                "upwind",

              wind_reference =
                "sample_median",

              wind_ms =
                median_wind_direction,

              estimate =
                up_eff$estimate,

              se =
                up_eff$se,

              ci95_low =
                up_eff$ci95_low,

              ci95_high =
                up_eff$ci95_high,

              p_wald =
                up_eff$p_wald,

              n_rows =
                nrow(
                  DD
                )
            ),

            data.table(
              fe_structure =
                fe_structure,

              fuel =
                fuel,

              metric =
                metric,

              radius_km =
                radius,

              component =
                "downwind",

              wind_reference =
                "sample_median",

              wind_ms =
                median_wind_direction,

              estimate =
                down_eff$estimate,

              se =
                down_eff$se,

              ci95_low =
                down_eff$ci95_low,

              ci95_high =
                down_eff$ci95_high,

              p_wald =
                down_eff$p_wald,

              n_rows =
                nrow(
                  DD
                )
            ),

            data.table(
              fe_structure =
                fe_structure,

              fuel =
                fuel,

              metric =
                metric,

              radius_km =
                radius,

              component =
                "upwind_minus_downwind",

              wind_reference =
                "sample_median",

              wind_ms =
                median_wind_direction,

              estimate =
                diff_eff$estimate,

              se =
                diff_eff$se,

              ci95_low =
                diff_eff$ci95_low,

              ci95_high =
                diff_eff$ci95_high,

              p_wald =
                diff_eff$p_wald,

              n_rows =
                nrow(
                  DD
                )
            ),

            fill = TRUE
          )
        }
      }
    }
  }

  rm(D0)
  invisible(gc())
}


# ============================================================
# 11. Save raw model outputs
# ============================================================

setorder(
  STATUS,
  model_family,
  fe_structure,
  fuel,
  metric,
  radius_km
)

setorder(
  UPWIND_EFFECTS,
  fe_structure,
  fuel,
  metric,
  radius_km
)

setorder(
  DIRECTION_EFFECTS,
  fe_structure,
  fuel,
  metric,
  radius_km,
  component
)

fwrite(
  STATUS,
  file.path(
    OUT,
    "01_model_status.csv"
  )
)

fwrite(
  UPWIND_EFFECTS,
  file.path(
    OUT,
    "02_Figure4_upwind_effects_all_FE.csv"
  )
)

fwrite(
  DIRECTION_EFFECTS,
  file.path(
    OUT,
    "03_Figure5_direction_joint_effects_all_FE.csv"
  )
)

fwrite(
  DAILY_AUDIT,
  file.path(
    OUT,
    "04_daily_join_audit.csv"
  )
)


# ============================================================
# 12. Cross-check Figure 4 against validated 07H15
# ============================================================

H15 <- fread(
  H15_EFFECT_FILE,
  showProgress = FALSE
)

OLD_UP <- H15[
  geometry ==
    "upwind_90" &
  wind_reference ==
    "sample_median" &
  fuel %in%
    FUELS &
  metric %in%
    METRICS &
  radius_km %in%
    RADII,
  .(
    fe_structure,
    fuel,
    metric,
    radius_km,

    old_estimate =
      estimate,

    old_ci95_low =
      ci95_low,

    old_ci95_high =
      ci95_high
  )
]

NEW_UP <- UPWIND_EFFECTS[
  ,
  .(
    fe_structure,
    fuel,
    metric,
    radius_km,

    new_estimate =
      estimate,

    new_ci95_low =
      ci95_low,

    new_ci95_high =
      ci95_high
  )
]

CHECK_H15 <- merge(
  NEW_UP,
  OLD_UP,
  by = c(
    "fe_structure",
    "fuel",
    "metric",
    "radius_km"
  ),
  all = TRUE,
  sort = FALSE
)

CHECK_H15[
  ,
  `:=`(
    estimate_abs_diff =
      abs(
        new_estimate -
        old_estimate
      ),

    ci_low_abs_diff =
      abs(
        new_ci95_low -
        old_ci95_low
      ),

    ci_high_abs_diff =
      abs(
        new_ci95_high -
        old_ci95_high
      )
  )
]

fwrite(
  CHECK_H15,
  file.path(
    OUT,
    "05_crosscheck_Figure4_vs_07H15.csv"
  )
)


# ============================================================
# 13. Cross-check Monitor×year Figure 5 against 07H10
# ============================================================

H10 <- fread(
  H10_DIRECTION_FILE,
  showProgress = FALSE
)

OLD_DIR <- H10[
  branch ==
    "filter_period" &
  component ==
    "上风向减下风向" &
  wind_reference ==
    "样本中位数" &
  fuel %in%
    FUELS &
  metric %in%
    METRICS &
  radius_km %in%
    RADII,
  .(
    fuel,
    metric,
    radius_km,

    old_estimate =
      estimate,

    old_ci95_low =
      ci95_low,

    old_ci95_high =
      ci95_high
  )
]

NEW_DIR <- DIRECTION_EFFECTS[
  fe_structure ==
    "monitor_year_FE" &
  component ==
    "upwind_minus_downwind",
  .(
    fuel,
    metric,
    radius_km,

    new_estimate =
      estimate,

    new_ci95_low =
      ci95_low,

    new_ci95_high =
      ci95_high
  )
]

CHECK_H10 <- merge(
  NEW_DIR,
  OLD_DIR,
  by = c(
    "fuel",
    "metric",
    "radius_km"
  ),
  all = TRUE,
  sort = FALSE
)

CHECK_H10[
  ,
  `:=`(
    estimate_abs_diff =
      abs(
        new_estimate -
        old_estimate
      ),

    ci_low_abs_diff =
      abs(
        new_ci95_low -
        old_ci95_low
      ),

    ci_high_abs_diff =
      abs(
        new_ci95_high -
        old_ci95_high
      )
  )
]

fwrite(
  CHECK_H10,
  file.path(
    OUT,
    "06_crosscheck_Figure5_vs_07H10.csv"
  )
)


# ============================================================
# 14. Plot helpers
# ============================================================

prepare_plot_data <- function(D) {

  X <- copy(D)

  X[
    ,
    fuel_label :=
      factor(
        FUEL_LABEL[
          fuel
        ],
        levels =
          unname(
            FUEL_LABEL[
              FUELS
            ]
          )
      )
  ]

  X
}

plot_upwind_metric <- function(
  D,
  metric_value,
  panel_title
) {

  X <- prepare_plot_data(
    D[
      metric ==
        metric_value
    ]
  )

  y_title <- if (
    metric_value ==
      "count"
  ) {

    "Adjusted gross-beta change per additional active plant (mBq/m³)"

  } else {

    "Adjusted gross-beta change per 10,000 MMBtu/day (mBq/m³)"
  }

  ggplot(
    X,
    aes(
      x = radius_km,
      y = estimate,
      group = 1
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
      width = 0.6,
      linewidth = 0.4
    ) +
    geom_line(
      linewidth = 0.75
    ) +
    geom_point(
      size = 2
    ) +
    facet_wrap(
      ~fuel_label,
      nrow = 1,
      scales = "free_y"
    ) +
    scale_x_continuous(
      breaks = RADII
    ) +
    labs(
      title =
        panel_title,

      x =
        "Cumulative radius (km)",

      y =
        y_title
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
        )
    )
}

plot_direction_metric <- function(
  D,
  metric_value,
  panel_title
) {

  X <- prepare_plot_data(
    D[
      metric ==
        metric_value &
      component ==
        "upwind_minus_downwind"
    ]
  )

  y_title <- if (
    metric_value ==
      "count"
  ) {

    "Upwind − downwind difference per additional active plant (mBq/m³)"

  } else {

    "Upwind − downwind difference per 10,000 MMBtu/day (mBq/m³)"
  }

  ggplot(
    X,
    aes(
      x = radius_km,
      y = estimate,
      group = 1
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
      width = 0.6,
      linewidth = 0.4
    ) +
    geom_line(
      linewidth = 0.75
    ) +
    geom_point(
      size = 2
    ) +
    facet_wrap(
      ~fuel_label,
      nrow = 1,
      scales = "free_y"
    ) +
    scale_x_continuous(
      breaks = RADII
    ) +
    labs(
      title =
        panel_title,

      x =
        "Cumulative radius (km)",

      y =
        y_title
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
        )
    )
}

save_both <- function(
  p,
  basename_out,
  width = 14,
  height = 8.5
) {

  ggsave(
    file.path(
      FIG_DIR,
      paste0(
        basename_out,
        ".png"
      )
    ),
    p,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  ggsave(
    file.path(
      FIG_DIR,
      paste0(
        basename_out,
        ".pdf"
      )
    ),
    p,
    width = width,
    height = height,
    units = "in"
  )
}


# ============================================================
# 15. MAIN Figure 4
#     Monitor × year FE
# ============================================================

F4 <- UPWIND_EFFECTS[
  fe_structure ==
    "monitor_year_FE"
]

F4A <- plot_upwind_metric(
  F4,
  "count",
  "A. Active plant count"
)

F4B <- plot_upwind_metric(
  F4,
  "heat",
  "B. Actual daily Heat Input"
)

FIG4 <- (
  F4A /
  F4B
) +
  plot_annotation(
    title =
      "Figure 4. Fuel-specific associations of 90° upwind actual operation with gross-beta particle radioactivity",

    subtitle =
      "Filter-period analysis; Monitor × year fixed effects; effects evaluated at the sample median wind speed.",

    theme =
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 14
          ),

        plot.subtitle =
          element_text(
            size = 10
          )
      )
  )

save_both(
  FIG4,
  "Figure4_upwind_actual_operation_monitorYearFE"
)


# ============================================================
# 16. MAIN Figure 5
#     Monitor × year FE
# ============================================================

F5 <- DIRECTION_EFFECTS[
  fe_structure ==
    "monitor_year_FE"
]

F5A <- plot_direction_metric(
  F5,
  "count",
  "A. Active plant count"
)

F5B <- plot_direction_metric(
  F5,
  "heat",
  "B. Actual daily Heat Input"
)

FIG5 <- (
  F5A /
  F5B
) +
  plot_annotation(
    title =
      "Figure 5. Direct directional comparison of upwind and downwind actual-operation associations",

    subtitle =
      "Filter-period analysis; Monitor × year fixed effects. Values are Upwind − Downwind from the same joint model; positive values indicate a stronger upwind association.",

    theme =
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 14
          ),

        plot.subtitle =
          element_text(
            size = 10
          )
      )
  )

save_both(
  FIG5,
  "Figure5_upwind_minus_downwind_monitorYearFE"
)


# ============================================================
# 17. Supplementary analogue of Figure 4
#     Monitor FE
# ============================================================

SF4 <- UPWIND_EFFECTS[
  fe_structure ==
    "monitor_FE"
]

SF4A <- plot_upwind_metric(
  SF4,
  "count",
  "A. Active plant count"
)

SF4B <- plot_upwind_metric(
  SF4,
  "heat",
  "B. Actual daily Heat Input"
)

SUPP_F4 <- (
  SF4A /
  SF4B
) +
  plot_annotation(
    title =
      "Supplementary Figure. Upwind actual-operation associations using Monitor FE",

    subtitle =
      "Filter-period analysis; Monitor fixed effects with linear and quadratic calendar-time trends; effects evaluated at the sample median wind speed.",

    theme =
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 14
          ),

        plot.subtitle =
          element_text(
            size = 10
          )
      )
  )

save_both(
  SUPP_F4,
  "Supplementary_Figure_upwind_actual_operation_monitorFE"
)


# ============================================================
# 18. Supplementary analogue of Figure 5
#     Monitor FE
# ============================================================

SF5 <- DIRECTION_EFFECTS[
  fe_structure ==
    "monitor_FE"
]

SF5A <- plot_direction_metric(
  SF5,
  "count",
  "A. Active plant count"
)

SF5B <- plot_direction_metric(
  SF5,
  "heat",
  "B. Actual daily Heat Input"
)

SUPP_F5 <- (
  SF5A /
  SF5B
) +
  plot_annotation(
    title =
      "Supplementary Figure. Direct upwind-versus-downwind comparison using Monitor FE",

    subtitle =
      "Monitor fixed effects with linear and quadratic calendar-time trends. Values are Upwind − Downwind from the same joint model.",

    theme =
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 14
          ),

        plot.subtitle =
          element_text(
            size = 10
          )
      )
  )

save_both(
  SUPP_F5,
  "Supplementary_Figure_upwind_minus_downwind_monitorFE"
)


# ============================================================
# 19. Figure-source CSVs
# ============================================================

fwrite(
  F4,
  file.path(
    OUT,
    "07_Figure4_source_monitorYearFE.csv"
  )
)

fwrite(
  F5,
  file.path(
    OUT,
    "08_Figure5_source_monitorYearFE.csv"
  )
)

fwrite(
  SF4,
  file.path(
    OUT,
    "09_Supplementary_Figure4_source_monitorFE.csv"
  )
)

fwrite(
  SF5,
  file.path(
    OUT,
    "10_Supplementary_Figure5_source_monitorFE.csv"
  )
)


# ============================================================
# 20. Hard checks
# ============================================================

N_FAIL <- STATUS[
  status ==
    "FAIL",
  .N
]

N_SKIP <- STATUS[
  status ==
    "SKIP",
  .N
]

N_PASS <- STATUS[
  status ==
    "PASS",
  .N
]

H15_MISSING <- CHECK_H15[
  !is.finite(
    new_estimate
  ) |
  !is.finite(
    old_estimate
  ),
  .N
]

H10_MISSING <- CHECK_H10[
  !is.finite(
    new_estimate
  ) |
  !is.finite(
    old_estimate
  ),
  .N
]

H15_MAX_DIFF <- if (
  H15_MISSING ==
    0L
) {

  max(
    c(
      CHECK_H15$estimate_abs_diff,
      CHECK_H15$ci_low_abs_diff,
      CHECK_H15$ci_high_abs_diff
    ),
    na.rm = TRUE
  )

} else {

  Inf
}

H10_MAX_DIFF <- if (
  H10_MISSING ==
    0L
) {

  max(
    c(
      CHECK_H10$estimate_abs_diff,
      CHECK_H10$ci_low_abs_diff,
      CHECK_H10$ci_high_abs_diff
    ),
    na.rm = TRUE
  )

} else {

  Inf
}

EXPECTED_FIGURES <- file.path(
  FIG_DIR,
  c(
    "Figure4_upwind_actual_operation_monitorYearFE.png",
    "Figure4_upwind_actual_operation_monitorYearFE.pdf",

    "Figure5_upwind_minus_downwind_monitorYearFE.png",
    "Figure5_upwind_minus_downwind_monitorYearFE.pdf",

    "Supplementary_Figure_upwind_actual_operation_monitorFE.png",
    "Supplementary_Figure_upwind_actual_operation_monitorFE.pdf",

    "Supplementary_Figure_upwind_minus_downwind_monitorFE.png",
    "Supplementary_Figure_upwind_minus_downwind_monitorFE.pdf"
  )
)

HARD <- data.table(
  check = c(
    "formal_filter_rows",
    "formal_filter_monitors",
    "daily_join_unmatched_zero",
    "model_status_rows",
    "model_FAIL_zero",
    "model_SKIP_zero",
    "model_PASS_224",
    "Figure4_effect_rows_112",
    "Figure5_joint_effect_rows_336",
    "H15_crosscheck_complete",
    "H15_crosscheck_within_tolerance",
    "H10_direction_crosscheck_complete",
    "H10_direction_crosscheck_within_tolerance",
    "four_png_four_pdf"
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
      sum(
        DAILY_AUDIT$unmatched_daily_rows
      )
    ),

    as.character(
      nrow(
        STATUS
      )
    ),

    as.character(
      N_FAIL
    ),

    as.character(
      N_SKIP
    ),

    as.character(
      N_PASS
    ),

    as.character(
      nrow(
        UPWIND_EFFECTS
      )
    ),

    as.character(
      nrow(
        DIRECTION_EFFECTS
      )
    ),

    as.character(
      H15_MISSING
    ),

    format(
      H15_MAX_DIFF,
      scientific = TRUE
    ),

    as.character(
      H10_MISSING
    ),

    format(
      H10_MAX_DIFF,
      scientific = TRUE
    ),

    as.character(
      sum(
        file.exists(
          EXPECTED_FIGURES
        )
      )
    )
  ),

  expected = c(
    "107407",
    "133",
    "0",
    "224",
    "0",
    "0",
    "224",
    "112",
    "336",
    "0",
    paste0(
      "<=",
      TOL
    ),
    "0",
    paste0(
      "<=",
      TOL
    ),
    "8"
  )
)

HARD[
  ,
  status := c(
    ifelse(
      nrow(FP) ==
        107407L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      uniqueN(
        FP$loc_num
      ) ==
        133L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      sum(
        DAILY_AUDIT$unmatched_daily_rows
      ) ==
        0L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      nrow(
        STATUS
      ) ==
        224L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      N_FAIL ==
        0L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      N_SKIP ==
        0L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      N_PASS ==
        224L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      nrow(
        UPWIND_EFFECTS
      ) ==
        112L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      nrow(
        DIRECTION_EFFECTS
      ) ==
        336L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      H15_MISSING ==
        0L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      is.finite(H15_MAX_DIFF),
      "INFO",
      "FAIL"
    ),

    ifelse(
      H10_MISSING ==
        0L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      is.finite(
        H10_MAX_DIFF
      ) &&
      H10_MAX_DIFF <=
        TOL,
      "PASS",
      "FAIL"
    ),

    ifelse(
      sum(
        file.exists(
          EXPECTED_FIGURES
        )
      ) ==
        8L,
      "PASS",
      "FAIL"
    )
  )
]

fwrite(
  HARD,
  file.path(
    OUT,
    "11_hard_checks.csv"
  )
)


# ============================================================
# 21. Terminal summary
# ============================================================

cat(
  "\n============================================================\n",
  "07H16 MODEL STATUS\n",
  "============================================================\n",
  sep = ""
)

print(
  STATUS[
    ,
    .N,
    by = .(
      model_family,
      fe_structure,
      status
    )
  ][
    order(
      model_family,
      fe_structure,
      status
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "FIGURE 4: MONITOR × YEAR FE SUMMARY\n",
  "============================================================\n",
  sep = ""
)

print(
  F4[
    ,
    .(
      n_radii = .N,

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

      n_ci_above_0 =
        sum(
          ci95_low >
            0,
          na.rm = TRUE
        ),

      n_ci_below_0 =
        sum(
          ci95_high <
            0,
          na.rm = TRUE
        ),

      n_ci_cross_0 =
        sum(
          ci95_low <=
            0 &
          ci95_high >=
            0,
          na.rm = TRUE
        )
    ),
    by = .(
      fuel,
      metric
    )
  ][
    order(
      fuel,
      metric
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "FIGURE 5: UPWIND − DOWNWIND SUMMARY\n",
  "============================================================\n",
  sep = ""
)

print(
  F5[
    component ==
      "upwind_minus_downwind",
    .(
      n_radii = .N,

      n_upwind_gt_downwind =
        sum(
          estimate >
            0,
          na.rm = TRUE
        ),

      n_downwind_gt_upwind =
        sum(
          estimate <
            0,
          na.rm = TRUE
        ),

      n_diff_ci_above_0 =
        sum(
          ci95_low >
            0,
          na.rm = TRUE
        ),

      n_diff_ci_below_0 =
        sum(
          ci95_high <
            0,
          na.rm = TRUE
        ),

      n_diff_ci_cross_0 =
        sum(
          ci95_low <=
            0 &
          ci95_high >=
            0,
          na.rm = TRUE
        )
    ),
    by = .(
      fuel,
      metric
    )
  ][
    order(
      fuel,
      metric
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "CROSS-CHECK\n",
  "============================================================\n",
  "07H15 max absolute difference = ",
  format(
    H15_MAX_DIFF,
    scientific = TRUE
  ),
  "\n",
  "07H10 direction max absolute difference = ",
  format(
    H10_MAX_DIFF,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "\n============================================================\n",
  "HARD CHECKS\n",
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

cat(
  "\n============================================================\n",
  ifelse(
    HARD_FAIL ==
      0L,
    "07H16 FIGURE 4 / FIGURE 5 FINAL BUILD: PASS\n",
    "07H16 FIGURE 4 / FIGURE 5 FINAL BUILD: FAIL\n"
  ),
  "Hard FAIL = ",
  HARD_FAIL,
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
    "07H16存在Hard FAIL，请先检查11_hard_checks.csv。"
  )
}
