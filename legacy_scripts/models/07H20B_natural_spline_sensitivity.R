#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999,
  width = 240
)

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(ggplot2)
})


# ============================================================
# 07H20B
# Natural-spline covariate functional-form sensitivity
#
# Purpose:
#   Refit the complete Figure 4 analysis (20–50 km)
#   using natural cubic splines for continuous background
#   covariates.
#
# Main exposure structure remains unchanged:
#
#   exposure
#   exposure × wind
#   Monitor × year FE
#   monitor-clustered SE
#
# Natural spline:
#   df = 3
#
# Applied uniformly to:
#   PM2.5
#   relative humidity
#   soil moisture
#   wind speed
#   continental air-mass proportion
#   monthly sunspot number
#   inverse PBLH
#   temperature
#
# Existing interactions retained:
#   sunspot × latitude
#   surface uranium × soil moisture
#
# IMPORTANT:
#   This is a sensitivity analysis.
#   It does NOT replace the frozen main model.
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

FIG4_SOURCE <- file.path(
  HROOT,
  "07H16_Figure4_Figure5正文与MonitorFE补充图",
  "07_Figure4_source_monitorYearFE.csv"
)

OUT <- file.path(
  HROOT,
  "07H20B_自然样条函数形式敏感性分析"
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

EXPECTED_ROWS <- 107407L
EXPECTED_MONITORS <- 133L

NS_DF <- 3L

FUEL_LABEL <- c(
  coal = "Coal",
  gas = "Natural gas",
  oil = "Oil",
  biomass = "Biomass"
)

METRIC_LABEL <- c(
  count = "Active plant count",
  heat = "Actual daily Heat Input"
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

  if (
    is.logical(x)
  ) {

    y <- x

  } else if (
    is.numeric(x) ||
    is.integer(x)
  ) {

    y <- x != 0

  } else {

    z <- tolower(
      trimws(
        as.character(x)
      )
    )

    y <- z %in% c(
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

  if (
    length(z) != 1L
  ) {

    return(
      NA_character_
    )
  }

  z
}


effect_at_wind <- function(
  fit,
  wind_ms
) {

  beta <- coef(
    fit
  )

  V <- as.matrix(
    vcov(
      fit
    )
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
        status = "FAIL",
        note =
          "Exposure or exposure×wind coefficient missing."
      )
    )
  }

  b1 <- unname(
    beta[[
      "exposure"
    ]]
  )

  b2 <- unname(
    beta[[
      int_term
    ]]
  )

  est <-
    b1 +
    wind_ms *
    b2

  vv <-
    V[
      "exposure",
      "exposure"
    ] +
    wind_ms^2 *
    V[
      int_term,
      int_term
    ] +
    2 *
    wind_ms *
    V[
      "exposure",
      int_term
    ]

  vv <- unname(
    vv
  )

  if (
    !is.finite(vv) ||
    vv < -1e-12
  ) {

    return(
      list(
        status = "FAIL",
        note =
          "Invalid variance for representative-wind exposure effect."
      )
    )
  }

  se <- sqrt(
    max(
      vv,
      0
    )
  )

  list(
    status = "PASS",

    estimate =
      unname(
        est
      ),

    se =
      unname(
        se
      ),

    ci95_low =
      unname(
        est -
        1.96 *
        se
      ),

    ci95_high =
      unname(
        est +
        1.96 *
        se
      ),

    p_wald =
      if (
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
      },

    note = ""
  )
}


# ============================================================
# 4. Formal monitor universe
# ============================================================

cat(
  "\n============================================================\n",
  "A. FORMAL MONITOR UNIVERSE\n",
  "============================================================\n",
  sep = ""
)

Q_LIST <- vector(
  "list",
  length(
    YEARS
  )
)

for (
  ii in seq_along(
    YEARS
  )
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

  if (
    !file.exists(f)
  ) {

    stopf(
      "缺少07F quasi-daily文件：",
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
    `:=`(
      loc_num =
        clean_chr(
          loc_num
        ),

      flag_07f_primary_quasi_sample =
        as_bool(
          flag_07f_primary_quasi_sample
        )
    )
  ]

  Q_LIST[[ii]] <- z
}

Q <- rbindlist(
  Q_LIST,
  use.names = TRUE,
  fill = TRUE
)

ELIGIBLE_LOCS <- sort(
  unique(
    Q[
      flag_07f_primary_quasi_sample %in%
        TRUE,
      loc_num
    ]
  )
)

rm(
  Q,
  Q_LIST
)

if (
  length(
    ELIGIBLE_LOCS
  ) !=
    EXPECTED_MONITORS
) {

  stopf(
    "Formal monitor universe不是133；实际=",
    length(
      ELIGIBLE_LOCS
    )
  )
}

cat(
  "Formal eligible monitors = 133：PASS\n"
)


# ============================================================
# 5. Formal filter-period sample
# ============================================================

cat(
  "\n============================================================\n",
  "B. FORMAL FILTER-PERIOD SAMPLE\n",
  "============================================================\n",
  sep = ""
)

PRIMARY_FP <-
  "fpmean_campd_upwind_heat_observed_20km_90deg"

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

  PRIMARY_FP
)

FP_HEADER <- names(
  fread(
    FP_FILE,
    nrows = 0,
    showProgress = FALSE
  )
)

FP_MISS <- setdiff(
  FP_REQUIRED,
  FP_HEADER
)

if (
  length(
    FP_MISS
  ) >
    0L
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

NUMERIC_COLS <- c(
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

  PRIMARY_FP
)

for (
  cc in NUMERIC_COLS
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
    get(
      PRIMARY_FP
    )
  )
]

if (
  nrow(
    FP
  ) !=
    EXPECTED_ROWS
) {

  stopf(
    "Formal filter-period不是107407；实际=",
    nrow(
      FP
    )
  )
}

if (
  uniqueN(
    FP$loc_num
  ) !=
    EXPECTED_MONITORS
) {

  stopf(
    "Formal filter-period不是133 monitors。"
  )
}

FP[
  ,
  `:=`(
    period_id =
      .I,

    period_days =
      as.integer(
        collect_end -
        collect_start
      ),

    pm_model_centered =
      pm_filter_period_official -
      10,

    temp_centered =
      temp -
      273.15,

    temp2 =
      (
        temp -
        273.15
      )^2,

    temp3 =
      (
        temp -
        273.15
      )^3,

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

cat(
  "Formal filter-period = ",
  format(
    nrow(
      FP
    ),
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
# 6. Build fixed natural-spline bases
#
# IMPORTANT:
# Knots are defined ONCE using the full formal filter-period
# sample and then reused for every fuel/radius/model.
# ============================================================

cat(
  "\n============================================================\n",
  "C. BUILD FIXED NATURAL-SPLINE BASES\n",
  "============================================================\n",
  sep = ""
)

SPLINE_VARS <- data.table(
  id = c(
    "cprop",
    "sunspot",
    "pm25",
    "invhpbl",
    "rhum",
    "soilm",
    "wind",
    "temperature"
  ),

  variable = c(
    "c_prop",
    "sunspot_monthly",
    "pm_model_centered",
    "inv_hpbl",
    "rhum",
    "soilm",
    "vel",
    "temp_centered"
  ),

  label = c(
    "Continental air-mass proportion",
    "Monthly sunspot number",
    "PM2.5",
    "Inverse PBLH",
    "Relative humidity",
    "Soil moisture",
    "Wind speed",
    "Temperature"
  )
)

SPLINE_SPEC <- data.table()

SPLINE_TERMS <- character()

for (
  ii in seq_len(
    nrow(
      SPLINE_VARS
    )
  )
) {

  id <- SPLINE_VARS$id[ii]
  var <- SPLINE_VARS$variable[ii]
  lab <- SPLINE_VARS$label[ii]

  x <- FP[[var]]

  ok <- is.finite(
    x
  )

  if (
    sum(
      ok
    ) <
      1000L
  ) {

    stopf(
      "Spline变量有效样本不足：",
      var
    )
  }

  B0 <- splines::ns(
    x[
      ok
    ],
    df = NS_DF,
    intercept = FALSE
  )

  B <- matrix(
    NA_real_,
    nrow =
      nrow(
        FP
      ),
    ncol =
      ncol(
        B0
      )
  )

  B[
    ok,
  ] <- B0

  basis_names <- paste0(
    "ns_",
    id,
    "_",
    seq_len(
      ncol(
        B
      )
    )
  )

  for (
    jj in seq_along(
      basis_names
    )
  ) {

    set(
      FP,
      j =
        basis_names[jj],
      value =
        B[
          ,
          jj
        ]
    )
  }

  knots <- attr(
    B0,
    "knots"
  )

  boundary <- attr(
    B0,
    "Boundary.knots"
  )

  SPLINE_SPEC <- rbind(
    SPLINE_SPEC,

    data.table(
      id =
        id,

      covariate =
        lab,

      variable =
        var,

      df =
        NS_DF,

      n_finite =
        sum(
          ok
        ),

      internal_knots =
        paste(
          format(
            knots,
            digits = 12,
            trim = TRUE
          ),
          collapse = " | "
        ),

      boundary_knots =
        paste(
          format(
            boundary,
            digits = 12,
            trim = TRUE
          ),
          collapse = " | "
        ),

      basis_columns =
        paste(
          basis_names,
          collapse = " | "
        )
    ),

    fill = TRUE
  )

  SPLINE_TERMS <- c(
    SPLINE_TERMS,
    basis_names
  )

  cat(
    lab,
    ": df=",
    NS_DF,
    " | finite=",
    sum(
      ok
    ),
    "\n",
    sep = ""
  )
}

fwrite(
  SPLINE_SPEC,
  file.path(
    OUT,
    "01_natural_spline_basis_specification.csv"
  )
)


# ============================================================
# 7. Read all 20–50 km actual-operation exposures
# ============================================================

cat(
  "\n============================================================\n",
  "D. BUILD 20–50 KM FUEL-SPECIFIC EXPOSURES\n",
  "============================================================\n",
  sep = ""
)

COUNT_COLS <- unlist(
  lapply(
    FUELS,
    function(fuel) {

      paste0(
        "campd_",
        fuel,
        "_upwind_n_plants_",
        RADII,
        "km_90deg"
      )
    }
  ),
  use.names = FALSE
)

HEAT_COLS <- unlist(
  lapply(
    FUELS,
    function(fuel) {

      paste0(
        "campd_",
        fuel,
        "_upwind_heat_observed_",
        RADII,
        "km_90deg"
      )
    }
  ),
  use.names = FALSE
)

EXPOSURE_COLS <- c(
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

if (
  !all(
    file.exists(
      DAILY_FILES
    )
  )
) {

  print(
    DAILY_FILES[
      !file.exists(
        DAILY_FILES
      )
    ]
  )

  stopf(
    "07H01 daily exposure files不是24/24。"
  )
}

DAILY <- rbindlist(
  lapply(
    DAILY_FILES,
    function(f) {

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
        length(
          miss
        ) >
          0L
      ) {

        stopf(
          basename(
            f
          ),
          " 缺字段：",
          paste(
            miss,
            collapse = ", "
          )
        )
      }

      fread(
        f,
        select =
          req,
        showProgress =
          FALSE
      )
    }
  ),
  use.names =
    TRUE,
  fill =
    TRUE
)

DAILY[
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
      )
  )
]

for (
  cc in EXPOSURE_COLS
) {

  set(
    DAILY,
    j =
      cc,
    value =
      safe_num(
        DAILY[[cc]]
      )
  )
}

DAILY_DUP <- DAILY[
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
  DAILY_DUP >
    0L
) {

  stopf(
    "07H01 daily exposure存在重复monitor-day。"
  )
}

DAILY[
  ,
  daily_match :=
    TRUE
]


# ============================================================
# 8. Aggregate daily exposure to filter periods
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
    loc_num
  )
]

FP_DAYS[
  ,
  date :=
    as.IDate(
      date
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

N_UNMATCHED <- FP_DAYS[
  is.na(
    daily_match
  ),
  .N
]

if (
  N_UNMATCHED >
    0L
) {

  stopf(
    "Filter-period daily exposure join有缺失：",
    N_UNMATCHED
  )
}

FP_EXP <- FP_DAYS[
  ,
  lapply(
    .SD,
    strict_mean
  ),
  by =
    period_id,
  .SDcols =
    EXPOSURE_COLS
]

setnames(
  FP_EXP,
  EXPOSURE_COLS,
  paste0(
    "fpmean_",
    EXPOSURE_COLS
  )
)

FP <- merge(
  FP,
  FP_EXP,
  by =
    "period_id",
  all.x =
    TRUE,
  sort =
    FALSE
)

rm(
  FP_DAYS,
  FP_EXP,
  DAILY
)

invisible(
  gc()
)


# ============================================================
# 9. Current main-model RHS
# ============================================================

CURRENT_RHS <- paste(
  c(
    "exposure",
    "exposure:vel",

    "c_prop",
    "sunspot_monthly",
    "pm_model_centered",
    "inv_hpbl",
    "rhum",
    "soilm",
    "vel",

    "temp_centered",
    "temp2",
    "temp3",

    "sunspot_monthly:lat_diff_35",
    "surface_u_ppm_article:soilm"
  ),
  collapse = " + "
)


# ============================================================
# 10. Natural-spline sensitivity RHS
#
# Raw main effects are replaced by the spline basis.
#
# The following raw interactions remain unchanged:
#   exposure × vel
#   sunspot × latitude
#   surface uranium × soil moisture
# ============================================================

NS_RHS <- paste(
  c(
    "exposure",
    "exposure:vel",

    SPLINE_TERMS,

    "sunspot_monthly:lat_diff_35",
    "surface_u_ppm_article:soilm"
  ),
  collapse = " + "
)


# ============================================================
# 11. Prepare one fuel × metric × radius sample
# ============================================================

prepare_model_data <- function(
  fuel,
  metric,
  radius
) {

  if (
    metric ==
      "count"
  ) {

    exposure_col <- paste0(
      "fpmean_campd_",
      fuel,
      "_upwind_n_plants_",
      radius,
      "km_90deg"
    )

    scale_value <- 1

  } else {

    exposure_col <- paste0(
      "fpmean_campd_",
      fuel,
      "_upwind_heat_observed_",
      radius,
      "km_90deg"
    )

    scale_value <- 10000
  }

  if (
    !exposure_col %in%
      names(
        FP
      )
  ) {

    stopf(
      "缺少暴露字段：",
      exposure_col
    )
  }

  d <- copy(
    FP
  )

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

  required <- c(
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
    "temp2",
    "temp3",

    SPLINE_TERMS,

    "loc_num",
    "monitor_year_fe"
  )

  d <- d[
    complete.cases(
      d[
        ,
        ..required
      ]
    )
  ]

  list(
    data =
      d,

    exposure_col =
      exposure_col,

    scale_value =
      scale_value
  )
}


# ============================================================
# 12. Fit one model
# ============================================================

fit_one <- function(
  d,
  model_type
) {

  rhs <- if (
    model_type ==
      "current"
  ) {

    CURRENT_RHS

  } else {

    NS_RHS
  }

  FORM <- as.formula(
    paste0(
      "gross_beta_mBq_m3 ~ ",
      rhs,
      " | monitor_year_fe"
    )
  )

  warning_text <- character()

  fit <- tryCatch(

    withCallingHandlers(

      feols(
        FORM,
        data =
          d,
        vcov =
          ~loc_num,
        notes =
          FALSE,
        warn =
          FALSE
      ),

      warning = function(w) {

        warning_text <<- c(
          warning_text,
          conditionMessage(
            w
          )
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
            conditionMessage(
              e
            )
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
          "FAIL",

        note =
          fit$message,

        warnings =
          paste(
            unique(
              warning_text
            ),
            collapse = " | "
          )
      )
    )
  }

  # Same representative-wind rule as Figure 4:
  # median wind of the model input sample.
  wind_ref <- median(
    d$vel,
    na.rm = TRUE
  )

  ee <- effect_at_wind(
    fit,
    wind_ref
  )

  if (
    ee$status !=
      "PASS"
  ) {

    return(
      list(
        status =
          "FAIL",

        note =
          ee$note,

        warnings =
          paste(
            unique(
              warning_text
            ),
            collapse = " | "
          )
      )
    )
  }

  list(
    status =
      "PASS",

    note =
      "",

    warnings =
      paste(
        unique(
          warning_text
        ),
        collapse = " | "
      ),

    n_input =
      nrow(
        d
      ),

    n_model =
      nobs(
        fit
      ),

    n_monitors =
      uniqueN(
        d$loc_num
      ),

    wind_ms =
      wind_ref,

    estimate =
      ee$estimate,

    se =
      ee$se,

    ci95_low =
      ee$ci95_low,

    ci95_high =
      ee$ci95_high,

    p_wald =
      ee$p_wald
  )
}


# ============================================================
# 13. Fit all 20–50 km models
# ============================================================

cat(
  "\n============================================================\n",
  "E. FIT CURRENT AND NATURAL-SPLINE MODELS\n",
  "============================================================\n",
  sep = ""
)

RESULTS <- data.table()

STATUS <- data.table()

for (
  radius in RADII
) {

  cat(
    "\nRadius = ",
    radius,
    " km\n",
    sep = ""
  )

  for (
    fuel in FUELS
  ) {

    for (
      metric in METRICS
    ) {

      cat(
        "  ",
        fuel,
        " / ",
        metric,
        "\n",
        sep = ""
      )

      P <- prepare_model_data(
        fuel,
        metric,
        radius
      )

      d <- P$data

      for (
        model_type in c(
          "current",
          "natural_spline"
        )
      ) {

        rr <- fit_one(
          d,
          model_type
        )

        STATUS <- rbind(
          STATUS,

          data.table(
            fuel =
              fuel,

            metric =
              metric,

            radius_km =
              radius,

            model =
              model_type,

            status =
              rr$status,

            note =
              rr$note,

            warnings =
              rr$warnings
          ),

          fill = TRUE
        )

        if (
          rr$status !=
            "PASS"
        ) {

          next
        }

        RESULTS <- rbind(
          RESULTS,

          data.table(
            fuel =
              fuel,

            metric =
              metric,

            radius_km =
              radius,

            model =
              model_type,

            spline_df =
              if (
                model_type ==
                  "natural_spline"
              ) {
                NS_DF
              } else {
                NA_integer_
              },

            N_input =
              rr$n_input,

            N_model =
              rr$n_model,

            N_monitors =
              rr$n_monitors,

            wind_reference =
              "sample_median",

            wind_ms =
              rr$wind_ms,

            estimate =
              rr$estimate,

            se =
              rr$se,

            ci95_low =
              rr$ci95_low,

            ci95_high =
              rr$ci95_high,

            p_wald =
              rr$p_wald
          ),

          fill = TRUE
        )
      }
    }
  }
}


# ============================================================
# 14. Fail-fast
# ============================================================

BAD <- STATUS[
  status !=
    "PASS"
]

if (
  nrow(
    BAD
  ) >
    0L
) {

  cat(
    "\nMODEL FAILURES:\n"
  )

  print(
    BAD,
    nrows = Inf
  )

  stopf(
    "Natural-spline analysis存在模型FAIL。"
  )
}


# ============================================================
# 15. Current vs natural spline comparison
# ============================================================

CURRENT <- RESULTS[
  model ==
    "current",
  .(
    fuel,
    metric,
    radius_km,

    N_input_current =
      N_input,

    N_model_current =
      N_model,

    N_monitors,

    wind_ms,

    current_estimate =
      estimate,

    current_se =
      se,

    current_ci95_low =
      ci95_low,

    current_ci95_high =
      ci95_high,

    current_p =
      p_wald
  )
]

NS <- RESULTS[
  model ==
    "natural_spline",
  .(
    fuel,
    metric,
    radius_km,

    N_input_ns =
      N_input,

    N_model_ns =
      N_model,

    ns_estimate =
      estimate,

    ns_se =
      se,

    ns_ci95_low =
      ci95_low,

    ns_ci95_high =
      ci95_high,

    ns_p =
      p_wald
  )
]

COMPARE <- merge(
  CURRENT,
  NS,
  by = c(
    "fuel",
    "metric",
    "radius_km"
  ),
  all = TRUE,
  sort = FALSE
)

COMPARE[
  ,
  `:=`(
    absolute_change =
      ns_estimate -
      current_estimate,

    percent_change =
      fifelse(
        abs(
          current_estimate
        ) >
          1e-12,

        100 *
        (
          ns_estimate -
          current_estimate
        ) /
        abs(
          current_estimate
        ),

        NA_real_
      ),

    sign_changed =
      sign(
        ns_estimate
      ) !=
      sign(
        current_estimate
      ),

    current_ci_above_zero =
      current_ci95_low >
      0,

    ns_ci_above_zero =
      ns_ci95_low >
      0,

    current_ci_below_zero =
      current_ci95_high <
      0,

    ns_ci_below_zero =
      ns_ci95_high <
      0,

    model_N_equal =
      N_model_current ==
      N_model_ns
  )
]


# ============================================================
# 16. Exact reproduction of current Figure 4
# ============================================================

cat(
  "\n============================================================\n",
  "F. CURRENT FIGURE 4 REPRODUCTION CHECK\n",
  "============================================================\n",
  sep = ""
)

if (
  !file.exists(
    FIG4_SOURCE
  )
) {

  stopf(
    "缺少Figure4 source：",
    FIG4_SOURCE
  )
}

SRC <- fread(
  FIG4_SOURCE,
  showProgress = FALSE
)

if (
  "fe_structure" %in%
    names(
      SRC
    )
) {

  SRC <- SRC[
    fe_structure ==
      "monitor_year_FE"
  ]
}

if (
  "direction" %in%
    names(
      SRC
    )
) {

  SRC <- SRC[
    direction ==
      "upwind"
  ]
}

if (
  "status" %in%
    names(
      SRC
    )
) {

  SRC <- SRC[
    status ==
      "PASS"
  ]
}

SRC <- SRC[
  fuel %in%
    FUELS &

  metric %in%
    METRICS &

  radius_km %in%
    RADII
]

if (
  "wind_reference" %in%
    names(
      SRC
    )
) {

  SRC <- SRC[
    grepl(
      "median|中位",
      wind_reference,
      ignore.case = TRUE
    )
  ]
}

if (
  nrow(
    SRC
  ) !=
    56L
) {

  cat(
    "\nFigure4 source rows after filtering:\n"
  )

  print(
    SRC
  )

  stopf(
    "Figure4 20–50 km source不是56行；实际=",
    nrow(
      SRC
    )
  )
}

CHECK <- merge(
  CURRENT,

  SRC[
    ,
    .(
      fuel,
      metric,
      radius_km,

      source_wind_ms =
        wind_ms,

      source_estimate =
        estimate,

      source_ci95_low =
        ci95_low,

      source_ci95_high =
        ci95_high
    )
  ],

  by = c(
    "fuel",
    "metric",
    "radius_km"
  ),

  all.x = TRUE,
  sort = FALSE
)

CHECK[
  ,
  `:=`(
    diff_estimate =
      abs(
        current_estimate -
        source_estimate
      ),

    diff_low =
      abs(
        current_ci95_low -
        source_ci95_low
      ),

    diff_high =
      abs(
        current_ci95_high -
        source_ci95_high
      ),

    diff_wind =
      abs(
        wind_ms -
        source_wind_ms
      )
  )
]

if (
  any(
    !is.finite(
      CHECK$diff_estimate
    )
  ) ||
  any(
    !is.finite(
      CHECK$diff_low
    )
  ) ||
  any(
    !is.finite(
      CHECK$diff_high
    )
  )
) {

  print(
    CHECK,
    nrows = Inf
  )

  stopf(
    "Figure4 reproduction check存在NA/Inf。"
  )
}

MAX_DIFF <- max(
  c(
    CHECK$diff_estimate,
    CHECK$diff_low,
    CHECK$diff_high
  )
)

cat(
  "Current-model reproduction max abs difference = ",
  format(
    MAX_DIFF,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

if (
  !is.finite(
    MAX_DIFF
  ) ||
  MAX_DIFF >
    1e-7
) {

  print(
    CHECK[
      order(
        -diff_estimate
      )
    ],
    nrows = Inf
  )

  stopf(
    "Current model did not reproduce Figure 4."
  )
}

cat(
  "Current Figure 4 reproduction：PASS\n"
)


# ============================================================
# 17. Biomass summary
# ============================================================

BIOMASS <- COMPARE[
  fuel ==
    "biomass"
][
  order(
    metric,
    radius_km
  )
]

BIOMASS_SUMMARY <- BIOMASS[
  ,
  .(
    n_radii =
      .N,

    n_current_positive =
      sum(
        current_estimate >
          0
      ),

    n_ns_positive =
      sum(
        ns_estimate >
          0
      ),

    n_current_ci_above_zero =
      sum(
        current_ci95_low >
          0
      ),

    n_ns_ci_above_zero =
      sum(
        ns_ci95_low >
          0
      ),

    n_sign_changed =
      sum(
        sign_changed
      ),

    median_percent_change =
      median(
        percent_change,
        na.rm = TRUE
      ),

    min_percent_change =
      min(
        percent_change,
        na.rm = TRUE
      ),

    max_percent_change =
      max(
        percent_change,
        na.rm = TRUE
      ),

    max_absolute_percent_change =
      max(
        abs(
          percent_change
        ),
        na.rm = TRUE
      )
  ),
  by =
    metric
]


# ============================================================
# 18. Save tables
# ============================================================

fwrite(
  STATUS,
  file.path(
    OUT,
    "02_model_status.csv"
  )
)

fwrite(
  RESULTS,
  file.path(
    OUT,
    "03_current_and_natural_spline_effects.csv"
  )
)

fwrite(
  COMPARE,
  file.path(
    OUT,
    "04_current_vs_natural_spline_effect_change.csv"
  )
)

fwrite(
  BIOMASS,
  file.path(
    OUT,
    "05_biomass_current_vs_natural_spline.csv"
  )
)

fwrite(
  BIOMASS_SUMMARY,
  file.path(
    OUT,
    "06_biomass_natural_spline_summary.csv"
  )
)

fwrite(
  CHECK,
  file.path(
    OUT,
    "07_Figure4_reproduction_check.csv"
  )
)


# ============================================================
# 19. Plot: all fuels
# ============================================================

PLOT <- copy(
  RESULTS
)

PLOT[
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

PLOT[
  ,
  metric_label :=
    factor(
      METRIC_LABEL[
        metric
      ],
      levels =
        unname(
          METRIC_LABEL[
            METRICS
          ]
        )
    )
]

PLOT[
  ,
  model_label :=
    factor(
      fifelse(
        model ==
          "current",
        "Current main model",
        "Natural spline sensitivity"
      ),
      levels = c(
        "Current main model",
        "Natural spline sensitivity"
      )
    )
]

pd <- position_dodge(
  width = 1.4
)

p_all <- ggplot(
  PLOT,
  aes(
    x =
      radius_km,
    y =
      estimate,
    group =
      model_label,
    linetype =
      model_label,
    shape =
      model_label
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35
  ) +

  geom_errorbar(
    aes(
      ymin =
        ci95_low,
      ymax =
        ci95_high
    ),
    position =
      pd,
    width =
      1.2,
    linewidth =
      0.4
  ) +

  geom_line(
    position =
      pd,
    linewidth =
      0.65
  ) +

  geom_point(
    position =
      pd,
    size =
      1.8
  ) +

  facet_grid(
    metric_label ~
      fuel_label,
    scales =
      "free_y"
  ) +

  scale_x_continuous(
    breaks =
      RADII
  ) +

  labs(
    title =
      "Fuel-specific associations: current vs natural-spline covariate models",

    subtitle =
      paste0(
        "20–50 km; Monitor × year FE; ",
        "natural cubic splines with df=3 for continuous background covariates"
      ),

    x =
      "Cumulative radius (km)",

    y =
      "Adjusted gross-beta change (mBq/m³)",

    linetype =
      NULL,

    shape =
      NULL
  ) +

  theme_bw(
    base_size = 10
  ) +

  theme(
    panel.grid.minor =
      element_blank(),

    strip.text =
      element_text(
        face = "bold"
      ),

    legend.position =
      "bottom"
  )

ggsave(
  file.path(
    FIG_DIR,
    "01_all_fuels_current_vs_natural_spline.png"
  ),
  p_all,
  width = 13,
  height = 7.5,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 20. Plot: biomass only
# ============================================================

P_BIO <- PLOT[
  fuel ==
    "biomass"
]

p_bio <- ggplot(
  P_BIO,
  aes(
    x =
      radius_km,
    y =
      estimate,
    group =
      model_label,
    linetype =
      model_label,
    shape =
      model_label
  )
) +

  geom_hline(
    yintercept =
      0,
    linetype =
      "dashed",
    linewidth =
      0.35
  ) +

  geom_errorbar(
    aes(
      ymin =
        ci95_low,
      ymax =
        ci95_high
    ),
    position =
      pd,
    width =
      1.2,
    linewidth =
      0.5
  ) +

  geom_line(
    position =
      pd,
    linewidth =
      0.75
  ) +

  geom_point(
    position =
      pd,
    size =
      2
  ) +

  facet_wrap(
    ~metric_label,
    scales =
      "free_y",
    nrow =
      1
  ) +

  scale_x_continuous(
    breaks =
      RADII
  ) +

  labs(
    title =
      "Biomass associations under natural-spline covariate adjustment",

    subtitle =
      "Monitor × year FE; current model vs natural-spline sensitivity",

    x =
      "Cumulative radius (km)",

    y =
      "Adjusted gross-beta change (mBq/m³)",

    linetype =
      NULL,

    shape =
      NULL
  ) +

  theme_bw(
    base_size = 10
  ) +

  theme(
    panel.grid.minor =
      element_blank(),

    strip.text =
      element_text(
        face = "bold"
      ),

    legend.position =
      "bottom"
  )

ggsave(
  file.path(
    FIG_DIR,
    "02_biomass_current_vs_natural_spline.png"
  ),
  p_bio,
  width = 10,
  height = 5.5,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 21. Hard checks
# ============================================================

PAIR_N <- merge(
  RESULTS[
    model ==
      "current",
    .(
      fuel,
      metric,
      radius_km,
      current_N =
        N_model
    )
  ],

  RESULTS[
    model ==
      "natural_spline",
    .(
      fuel,
      metric,
      radius_km,
      ns_N =
        N_model
    )
  ],

  by = c(
    "fuel",
    "metric",
    "radius_km"
  ),

  all = TRUE
)

N_PAIR_MISMATCH <- PAIR_N[
  current_N !=
    ns_N |
  is.na(
    current_N
  ) |
  is.na(
    ns_N
  ),
  .N
]

HARD <- data.table(
  check = c(
    "formal_rows",
    "formal_monitors",
    "daily_files",
    "daily_join_unmatched",
    "current_models_PASS",
    "natural_spline_models_PASS",
    "current_vs_ns_model_N_equal",
    "Figure4_source_rows",
    "Figure4_reproduction"
  ),

  observed = c(
    as.character(
      nrow(
        FP
      )
    ),

    as.character(
      uniqueN(
        FP$loc_num
      )
    ),

    as.character(
      sum(
        file.exists(
          DAILY_FILES
        )
      )
    ),

    as.character(
      N_UNMATCHED
    ),

    as.character(
      STATUS[
        model ==
          "current" &
        status ==
          "PASS",
        .N
      ]
    ),

    as.character(
      STATUS[
        model ==
          "natural_spline" &
        status ==
          "PASS",
        .N
      ]
    ),

    as.character(
      N_PAIR_MISMATCH
    ),

    as.character(
      nrow(
        SRC
      )
    ),

    format(
      MAX_DIFF,
      scientific = TRUE
    )
  ),

  expected = c(
    "107407",
    "133",
    "24",
    "0",
    "56",
    "56",
    "0",
    "56",
    "<=1e-7"
  ),

  status = c(
    ifelse(
      nrow(
        FP
      ) ==
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
        file.exists(
          DAILY_FILES
        )
      ) ==
        24L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      N_UNMATCHED ==
        0L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      STATUS[
        model ==
          "current" &
        status ==
          "PASS",
        .N
      ] ==
        56L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      STATUS[
        model ==
          "natural_spline" &
        status ==
          "PASS",
        .N
      ] ==
        56L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      N_PAIR_MISMATCH ==
        0L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      nrow(
        SRC
      ) ==
        56L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      is.finite(
        MAX_DIFF
      ) &&
      MAX_DIFF <=
        1e-7,
      "PASS",
      "FAIL"
    )
  )
)

fwrite(
  HARD,
  file.path(
    OUT,
    "08_hard_checks.csv"
  )
)


# ============================================================
# 22. Console summary
# ============================================================

cat(
  "\n============================================================\n",
  "BIOMASS: CURRENT VS NATURAL SPLINE\n",
  "============================================================\n",
  sep = ""
)

print(
  BIOMASS,
  nrows = Inf
)


cat(
  "\n============================================================\n",
  "BIOMASS SUMMARY\n",
  "============================================================\n",
  sep = ""
)

print(
  BIOMASS_SUMMARY,
  nrows = Inf
)


cat(
  "\n============================================================\n",
  "ALL FUEL-SPECIFIC EFFECT CHANGES\n",
  "============================================================\n",
  sep = ""
)

print(
  COMPARE[
    order(
      match(
        fuel,
        FUELS
      ),
      match(
        metric,
        METRICS
      ),
      radius_km
    )
  ],
  nrows = Inf
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


N_FAIL <- HARD[
  status ==
    "FAIL",
  .N
]

cat(
  "\n============================================================\n",
  ifelse(
    N_FAIL ==
      0L,
    "07H20B NATURAL-SPLINE SENSITIVITY: PASS\n",
    "07H20B NATURAL-SPLINE SENSITIVITY: FAIL\n"
  ),
  "Hard FAIL = ",
  N_FAIL,
  "\n",
  "Output = ",
  OUT,
  "\n",
  "============================================================\n",
  sep = ""
)

if (
  N_FAIL >
    0L
) {

  stopf(
    "07H20B存在Hard FAIL。"
  )
}

