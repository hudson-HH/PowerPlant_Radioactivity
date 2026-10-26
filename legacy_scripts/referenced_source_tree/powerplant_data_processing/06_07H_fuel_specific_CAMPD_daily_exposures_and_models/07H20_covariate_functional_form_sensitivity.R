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
# 07H20
# 协变量函数形式敏感性分析
#
# 目的：
# 1. 检查主要连续环境协变量是否存在二次/三次非线性
# 2. 不改变正式主模型
# 3. 建立一个统一的 flexible-covariate sensitivity model
# 4. 检查 Figure 4 的20 km燃料特异效应是否明显改变
#
# 主分析保持：
#   Monitor × year FE
#   exposure × wind
#   clustered SE by monitor
#
# 注意：
#   wind非线性这里只检查wind主效应；
#   exposure × wind仍维持正式模型中的线性交互。
# ============================================================


# ============================================================
# 1. PATHS
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
  "07H20_协变量函数形式敏感性分析"
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
# 2. PARAMETERS
# ============================================================

YEARS <- 2001:2024

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
# 3. HELPERS
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

  z
}


# ============================================================
# clustered joint Wald test
#
# Terms use SIMPLE explicit variable names.
# No I(x^2) coefficient-name matching.
# ============================================================

wald_joint <- function(
  fit,
  terms
) {

  beta <- coef(
    fit
  )

  V <- as.matrix(
    vcov(
      fit
    )
  )

  if (
    !all(
      terms %in%
        names(beta)
    )
  ) {

    return(
      list(
        statistic = NA_real_,
        df = NA_integer_,
        p = NA_real_,
        note =
          paste0(
            "Missing coefficient(s): ",
            paste(
              setdiff(
                terms,
                names(beta)
              ),
              collapse = ", "
            )
          )
      )
    )
  }

  b <- unname(
    beta[
      terms
    ]
  )

  VV <- V[
    terms,
    terms,
    drop = FALSE
  ]

  if (
    any(
      !is.finite(b)
    ) ||
    any(
      !is.finite(VV)
    )
  ) {

    return(
      list(
        statistic = NA_real_,
        df = NA_integer_,
        p = NA_real_,
        note =
          "Non-finite coefficient or covariance."
      )
    )
  }

  VV <-
    0.5 *
    (
      VV +
      t(VV)
    )

  ee <- eigen(
    VV,
    symmetric = TRUE
  )

  vmax <- max(
    abs(
      ee$values
    )
  )

  if (
    !is.finite(vmax) ||
    vmax <= 0
  ) {

    return(
      list(
        statistic = NA_real_,
        df = NA_integer_,
        p = NA_real_,
        note =
          "Covariance matrix has zero numerical rank."
      )
    )
  }

  tol <-
    vmax *
    1e-10

  keep <-
    ee$values >
    tol

  rk <- sum(
    keep
  )

  if (
    rk == 0L
  ) {

    return(
      list(
        statistic = NA_real_,
        df = NA_integer_,
        p = NA_real_,
        note =
          "Covariance matrix has zero positive rank."
      )
    )
  }

  Q <- ee$vectors[
    ,
    keep,
    drop = FALSE
  ]

  lambda <- ee$values[
    keep
  ]

  Vinv <-
    Q %*%
    diag(
      1 /
      lambda,
      nrow =
        rk
    ) %*%
    t(Q)

  stat <- as.numeric(
    t(b) %*%
    Vinv %*%
    b
  )

  p <- pchisq(
    stat,
    df =
      rk,
    lower.tail =
      FALSE
  )

  list(
    statistic =
      stat,
    df =
      rk,
    p =
      p,
    note =
      ""
  )
}


# ============================================================
# exposure effect at a specified wind speed
# ============================================================

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
        status =
          "FAIL",
        note =
          "exposure or exposure×wind coefficient missing"
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
        status =
          "FAIL",
        note =
          "invalid variance for exposure effect"
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
    status =
      "PASS",

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

    p =
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
          lower.tail =
            FALSE
        )
      } else {
        NA_real_
      },

    note =
      ""
  )
}


# ============================================================
# 4. FORMAL MONITOR UNIVERSE
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

  Q_LIST[[ii]] <-
    z
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
    "Formal eligible monitors不是133；实际=",
    length(
      ELIGIBLE_LOCS
    )
  )
}

cat(
  "Formal eligible monitors = 133：PASS\n"
)


# ============================================================
# 5. FORMAL FILTER-PERIOD SAMPLE
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
  select =
    FP_REQUIRED,
  showProgress =
    TRUE
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
    j =
      cc,
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
      drop =
        TRUE,
      lex.order =
        TRUE
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
# 6. READ 20-KM FUEL-SPECIFIC DAILY EXPOSURES
# ============================================================

cat(
  "\n============================================================\n",
  "C. BUILD 20-KM FUEL-SPECIFIC EXPOSURES\n",
  "============================================================\n",
  sep = ""
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
        EXPOSURES
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
  cc in EXPOSURES
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
# 7. AGGREGATE DAILY EXPOSURE TO FILTER PERIOD
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
  all.x =
    TRUE,
  sort =
    FALSE
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
    "filter-period逐日暴露连接缺失：",
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
    EXPOSURES
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
# 8. CURRENT COVARIATE MODEL
# ============================================================

BASE_TERMS <- c(
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
)

BASE_RHS <- paste(
  BASE_TERMS,
  collapse =
    " + "
)

EXPOSURE_RHS <- paste(
  c(
    "exposure",
    "exposure:vel",
    BASE_TERMS
  ),
  collapse =
    " + "
)


# ============================================================
# 9. COMMON DIAGNOSTIC SAMPLE
# ============================================================

DIAG_VARS <- c(
  "gross_beta_mBq_m3",
  "loc_num",
  "monitor_year_fe",
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
  "temp3"
)

DIAG <- FP[
  complete.cases(
    FP[
      ,
      ..DIAG_VARS
    ]
  )
]

cat(
  "\nFunctional-form diagnostic sample = ",
  format(
    nrow(
      DIAG
    ),
    big.mark = ","
  ),
  " / ",
  uniqueN(
    DIAG$loc_num
  ),
  " monitors\n",
  sep = ""
)


# ============================================================
# 10. CANDIDATES
# ============================================================

CANDIDATES <- data.table(
  id = c(
    "pm25",
    "rhum",
    "soilm",
    "wind",
    "cprop",
    "sunspot",
    "invhpbl",
    "temperature"
  ),

  variable = c(
    "pm_model_centered",
    "rhum",
    "soilm",
    "vel",
    "c_prop",
    "sunspot_monthly",
    "inv_hpbl",
    "temp_centered"
  ),

  label = c(
    "PM2.5",
    "Relative humidity",
    "Soil moisture",
    "Wind speed",
    "Continental air-mass proportion",
    "Monthly sunspot number",
    "Inverse PBLH",
    "Temperature"
  )
)


# ============================================================
# 11. FUNCTIONAL-FORM TESTS
# ============================================================

cat(
  "\n============================================================\n",
  "D. COVARIATE FUNCTIONAL-FORM DIAGNOSTICS\n",
  "============================================================\n",
  sep = ""
)

NONLINEAR <- data.table()

for (
  ii in seq_len(
    nrow(
      CANDIDATES
    )
  )
) {

  id <- CANDIDATES$id[ii]
  var <- CANDIDATES$variable[ii]
  lab <- CANDIDATES$label[ii]

  cat(
    "Testing ",
    lab,
    "...\n",
    sep = ""
  )

  d <- copy(
    DIAG
  )


  # ----------------------------------------------------------
  # Temperature:
  # directly test the current quadratic and cubic terms
  # ----------------------------------------------------------

  if (
    id ==
      "temperature"
  ) {

    FORM <- as.formula(
      paste0(
        "gross_beta_mBq_m3 ~ ",
        BASE_RHS,
        " | monitor_year_fe"
      )
    )

    fit <- feols(
      FORM,
      data =
        d,
      vcov =
        ~loc_num,
      notes =
        FALSE,
      warn =
        FALSE
    )

    ww <- wald_joint(
      fit,
      c(
        "temp2",
        "temp3"
      )
    )

    term2 <-
      "temp2"

    term3 <-
      "temp3"

  } else {

    mu <- mean(
      d[[var]],
      na.rm = TRUE
    )

    ss <- sd(
      d[[var]],
      na.rm = TRUE
    )

    if (
      !is.finite(
        ss
      ) ||
      ss <=
        0
    ) {

      stopf(
        lab,
        " SD无效。"
      )
    }

    z <- (
      d[[var]] -
      mu
    ) /
      ss

    d[
      ,
      nl2 :=
        z^2
    ]

    d[
      ,
      nl3 :=
        z^3
    ]

    FORM <- as.formula(
      paste0(
        "gross_beta_mBq_m3 ~ ",
        BASE_RHS,
        " + nl2 + nl3",
        " | monitor_year_fe"
      )
    )

    fit <- feols(
      FORM,
      data =
        d,
      vcov =
        ~loc_num,
      notes =
        FALSE,
      warn =
        FALSE
    )

    ww <- wald_joint(
      fit,
      c(
        "nl2",
        "nl3"
      )
    )

    term2 <-
      "nl2"

    term3 <-
      "nl3"
  }


  NONLINEAR <- rbind(
    NONLINEAR,

    data.table(
      id =
        id,

      covariate =
        lab,

      variable =
        var,

      N =
        nobs(
          fit
        ),

      N_monitors =
        uniqueN(
          d$loc_num
        ),

      quadratic_term =
        term2,

      cubic_term =
        term3,

      Wald_chisq =
        ww$statistic,

      df =
        ww$df,

      P_nonlinearity =
        ww$p,

      note =
        ww$note
    ),

    fill =
      TRUE
  )
}


# ------------------------------------------------------------
# Hard stop BEFORE plotting
# ------------------------------------------------------------

BAD_NL <- NONLINEAR[
  !is.finite(
    P_nonlinearity
  )
]

if (
  nrow(
    BAD_NL
  ) >
    0L
) {

  cat(
    "\nNONLINEARITY TEST FAILURES:\n"
  )

  print(
    BAD_NL,
    nrows =
      Inf
  )

  stopf(
    "存在无法计算的非线性联合检验。"
  )
}


NONLINEAR[
  ,
  P_BH :=
    p.adjust(
      P_nonlinearity,
      method =
        "BH"
    )
]

NONLINEAR[
  ,
  conclusion :=
    fifelse(
      P_BH <
        0.05,
      "Evidence of nonlinearity",
      "No strong evidence of nonlinearity"
    )
]

fwrite(
  NONLINEAR,
  file.path(
    OUT,
    "01_covariate_nonlinearity_tests.csv"
  )
)


# ============================================================
# 12. NONLINEARITY TEST FIGURE
# ============================================================

PPLOT <- copy(
  NONLINEAR
)

PPLOT[
  ,
  neg_log10_BH :=
    -log10(
      pmax(
        P_BH,
        1e-300
      )
    )
]

PPLOT[
  ,
  covariate :=
    factor(
      covariate,
      levels =
        rev(
          covariate[
            order(
              P_BH
            )
          ]
        )
    )
]

p1 <- ggplot(
  PPLOT,
  aes(
    x =
      neg_log10_BH,
    y =
      covariate
  )
) +

  geom_vline(
    xintercept =
      -log10(
        0.05
      ),
    linetype =
      "dashed",
    linewidth =
      0.45
  ) +

  geom_col(
    width =
      0.65
  ) +

  labs(
    title =
      "Covariate nonlinearity diagnostics",

    subtitle =
      "Joint clustered Wald test of quadratic and cubic terms; BH-adjusted P values",

    x =
      "-log10(BH-adjusted P)",

    y =
      NULL
  ) +

  theme_bw(
    base_size =
      10
  ) +

  theme(
    panel.grid.minor =
      element_blank()
  )

ggsave(
  file.path(
    FIG_DIR,
    "01_covariate_nonlinearity_tests.png"
  ),
  p1,
  width =
    8.5,
  height =
    6.5,
  dpi =
    300,
  bg =
    "white"
)


# ============================================================
# 13. PREPARE FUEL-SPECIFIC MODEL DATA
# ============================================================

prepare_exposure_data <- function(
  fuel,
  metric
) {

  if (
    metric ==
      "count"
  ) {

    exp_col <- paste0(
      "fpmean_campd_",
      fuel,
      "_upwind_n_plants_20km_90deg"
    )

    scale_value <- 1

  } else {

    exp_col <- paste0(
      "fpmean_campd_",
      fuel,
      "_upwind_heat_observed_20km_90deg"
    )

    scale_value <- 10000
  }

  if (
    !exp_col %in%
      names(
        FP
      )
  ) {

    stopf(
      "缺少暴露字段：",
      exp_col
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
          exp_col
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

  d
}


# ============================================================
# 14. FLEXIBLE COVARIATE TERMS
# ============================================================

FLEX_VARS <- data.table(
  id = c(
    "pm25",
    "rhum",
    "soilm",
    "wind",
    "cprop",
    "sunspot",
    "invhpbl"
  ),

  variable = c(
    "pm_model_centered",
    "rhum",
    "soilm",
    "vel",
    "c_prop",
    "sunspot_monthly",
    "inv_hpbl"
  )
)


add_flexible_terms <- function(
  d
) {

  x <- copy(
    d
  )

  new_terms <- character()

  for (
    ii in seq_len(
      nrow(
        FLEX_VARS
      )
    )
  ) {

    id <- FLEX_VARS$id[ii]
    var <- FLEX_VARS$variable[ii]

    mu <- mean(
      x[[var]],
      na.rm =
        TRUE
    )

    ss <- sd(
      x[[var]],
      na.rm =
        TRUE
    )

    if (
      !is.finite(
        ss
      ) ||
      ss <=
        0
    ) {

      stopf(
        "Flexible model SD无效：",
        var
      )
    }

    z <- (
      x[[var]] -
      mu
    ) /
      ss

    z2_name <- paste0(
      id,
      "_z2"
    )

    z3_name <- paste0(
      id,
      "_z3"
    )

    set(
      x,
      j =
        z2_name,
      value =
        z^2
    )

    set(
      x,
      j =
        z3_name,
      value =
        z^3
    )

    new_terms <- c(
      new_terms,
      z2_name,
      z3_name
    )
  }

  list(
    data =
      x,

    terms =
      new_terms
  )
}


# ============================================================
# 15. FIT ONE EXPOSURE MODEL
# ============================================================

fit_one <- function(
  d,
  flexible =
    FALSE
) {

  if (
    flexible
  ) {

    zz <- add_flexible_terms(
      d
    )

    dd <- zz$data

    rhs <- paste(
      c(
        EXPOSURE_RHS,
        zz$terms
      ),
      collapse =
        " + "
    )

  } else {

    dd <- copy(
      d
    )

    rhs <-
      EXPOSURE_RHS
  }

  FORM <- as.formula(
    paste0(
      "gross_beta_mBq_m3 ~ ",
      rhs,
      " | monitor_year_fe"
    )
  )

  fit <- tryCatch(

    feols(
      FORM,
      data =
        dd,
      vcov =
        ~loc_num,
      notes =
        FALSE,
      warn =
        FALSE
    ),

    error = function(e) {
      e
    }
  )

  if (
    inherits(
      fit,
      "error"
    )
  ) {

    return(
      list(
        status =
          "FAIL",

        note =
          conditionMessage(
            fit
          )
      )
    )
  }

  w <- median(
    dd$vel,
    na.rm =
      TRUE
  )

  ee <- effect_at_wind(
    fit,
    w
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
          ee$note
      )
    )
  }

  list(
    status =
      "PASS",

    note =
      "",

    n =
      nrow(
        dd
      ),

    n_monitors =
      uniqueN(
        dd$loc_num
      ),

    wind_ms =
      w,

    estimate =
      ee$estimate,

    se =
      ee$se,

    ci95_low =
      ee$ci95_low,

    ci95_high =
      ee$ci95_high,

    p =
      ee$p
  )
}


# ============================================================
# 16. CURRENT VS FLEXIBLE MODELS
# ============================================================

cat(
  "\n============================================================\n",
  "E. FIGURE 4 20-KM BASELINE VS FLEXIBLE COVARIATE MODEL\n",
  "============================================================\n",
  sep = ""
)

RESULTS <- data.table()

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
      "\n",
      sep = ""
    )

    d <- prepare_exposure_data(
      fuel,
      metric
    )

    rr0 <- fit_one(
      d,
      flexible =
        FALSE
    )

    rr1 <- fit_one(
      d,
      flexible =
        TRUE
    )

    for (
      mm in c(
        "current",
        "flexible"
      )
    ) {

      rr <- if (
        mm ==
          "current"
      ) {
        rr0
      } else {
        rr1
      }

      RESULTS <- rbind(
        RESULTS,

        data.table(
          fuel =
            fuel,

          metric =
            metric,

          model =
            mm,

          status =
            rr$status,

          N =
            if (
              rr$status ==
                "PASS"
            ) {
              rr$n
            } else {
              NA_integer_
            },

          N_monitors =
            if (
              rr$status ==
                "PASS"
            ) {
              rr$n_monitors
            } else {
              NA_integer_
            },

          wind_ms =
            if (
              rr$status ==
                "PASS"
            ) {
              rr$wind_ms
            } else {
              NA_real_
            },

          estimate =
            if (
              rr$status ==
                "PASS"
            ) {
              rr$estimate
            } else {
              NA_real_
            },

          se =
            if (
              rr$status ==
                "PASS"
            ) {
              rr$se
            } else {
              NA_real_
            },

          ci95_low =
            if (
              rr$status ==
                "PASS"
            ) {
              rr$ci95_low
            } else {
              NA_real_
            },

          ci95_high =
            if (
              rr$status ==
                "PASS"
            ) {
              rr$ci95_high
            } else {
              NA_real_
            },

          p_wald =
            if (
              rr$status ==
                "PASS"
            ) {
              rr$p
            } else {
              NA_real_
            },

          note =
            rr$note
        ),

        fill =
          TRUE
      )
    }
  }
}


# ============================================================
# 17. FAIL-FAST MODEL CHECK
# ============================================================

BAD <- RESULTS[
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
    BAD[
      ,
      .(
        fuel,
        metric,
        model,
        status,
        note
      )
    ],
    nrows =
      Inf
  )

  stopf(
    "存在未成功拟合的current/flexible模型。"
  )
}


# ============================================================
# 18. CURRENT VS FLEXIBLE COMPARISON
# ============================================================

CURRENT <- RESULTS[
  model ==
    "current",
  .(
    fuel,
    metric,
    N,
    N_monitors,
    wind_ms,
    current_estimate =
      estimate,
    current_ci95_low =
      ci95_low,
    current_ci95_high =
      ci95_high,
    current_p =
      p_wald
  )
]

FLEX <- RESULTS[
  model ==
    "flexible",
  .(
    fuel,
    metric,
    flexible_estimate =
      estimate,
    flexible_ci95_low =
      ci95_low,
    flexible_ci95_high =
      ci95_high,
    flexible_p =
      p_wald
  )
]

COMPARE <- merge(
  CURRENT,
  FLEX,
  by = c(
    "fuel",
    "metric"
  ),
  all =
    TRUE,
  sort =
    FALSE
)

COMPARE[
  ,
  `:=`(
    absolute_change =
      flexible_estimate -
      current_estimate,

    percent_change =
      fifelse(
        abs(
          current_estimate
        ) >
          1e-12,

        100 *
        (
          flexible_estimate -
          current_estimate
        ) /
        abs(
          current_estimate
        ),

        NA_real_
      ),

    sign_changed =
      sign(
        flexible_estimate
      ) !=
      sign(
        current_estimate
      )
  )
]

fwrite(
  RESULTS,
  file.path(
    OUT,
    "02_current_vs_flexible_model_results.csv"
  )
)

fwrite(
  COMPARE,
  file.path(
    OUT,
    "03_current_vs_flexible_effect_change.csv"
  )
)


# ============================================================
# 19. EXACT FIGURE 4 REPRODUCTION CHECK
# ============================================================

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
  showProgress =
    FALSE
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

SRC <- SRC[
  radius_km ==
    20 &
  fuel %in%
    FUELS &
  metric %in%
    METRICS &
  grepl(
    "median|中位",
    wind_reference,
    ignore.case =
      TRUE
  )
]

if (
  nrow(
    SRC
  ) !=
    8L
) {

  stopf(
    "Figure4 20-km median-wind source不是8行；实际=",
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
    "metric"
  ),

  all.x =
    TRUE,

  sort =
    FALSE
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
    CHECK
  )

  stopf(
    "Figure4 reproduction check中存在NA/Inf。"
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
  "\nCurrent-model reproduction max abs difference = ",
  format(
    MAX_DIFF,
    scientific =
      TRUE
  ),
  "\n",
  sep = ""
)

if (
  MAX_DIFF >
    1e-7
) {

  print(
    CHECK,
    nrows =
      Inf
  )

  stopf(
    "Current model did not reproduce Figure 4."
  )
}

cat(
  "Current Figure 4 reproduction：PASS\n"
)


# ============================================================
# 20. EFFECT COMPARISON FIGURE
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
        "Flexible covariates"
      ),
      levels = c(
        "Current main model",
        "Flexible covariates"
      )
    )
]

p2 <- ggplot(
  PLOT,
  aes(
    x =
      model_label,
    y =
      estimate
  )
) +

  geom_hline(
    yintercept =
      0,
    linetype =
      "dashed",
    linewidth =
      0.4
  ) +

  geom_errorbar(
    aes(
      ymin =
        ci95_low,
      ymax =
        ci95_high
    ),
    width =
      0.15,
    linewidth =
      0.45
  ) +

  geom_point(
    size =
      2
  ) +

  facet_grid(
    metric_label ~
      fuel_label,
    scales =
      "free_y"
  ) +

  labs(
    title =
      "20-km fuel-specific associations: current vs flexible covariate forms",

    subtitle =
      "Monitor × year FE; effects evaluated at sample median wind",

    x =
      NULL,

    y =
      "Adjusted gross-beta change (mBq/m³)"
  ) +

  theme_bw(
    base_size =
      10
  ) +

  theme(
    panel.grid.minor =
      element_blank(),

    strip.text =
      element_text(
        face =
          "bold"
      ),

    axis.text.x =
      element_text(
        angle =
          25,
        hjust =
          1
      )
  )

ggsave(
  file.path(
    FIG_DIR,
    "02_Figure4_20km_current_vs_flexible_covariates.png"
  ),
  p2,
  width =
    13,
  height =
    7.5,
  dpi =
    300,
  bg =
    "white"
)


# ============================================================
# 21. HARD CHECKS
# ============================================================

HARD <- data.table(
  check = c(
    "formal_rows",
    "formal_monitors",
    "daily_files",
    "daily_join_unmatched",
    "nonlinearity_tests_finite",
    "exposure_models_PASS",
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
      NONLINEAR[
        is.finite(
          P_nonlinearity
        ),
        .N
      ]
    ),

    as.character(
      RESULTS[
        status ==
          "PASS",
        .N
      ]
    ),

    format(
      MAX_DIFF,
      scientific =
        TRUE
    )
  ),

  expected = c(
    "107407",
    "133",
    "24",
    "0",
    "8",
    "16",
    "<=1e-7"
  )
)

HARD[
  ,
  status := c(
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
      NONLINEAR[
        is.finite(
          P_nonlinearity
        ),
        .N
      ] ==
        8L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      RESULTS[
        status ==
          "PASS",
        .N
      ] ==
        16L,
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
]

fwrite(
  HARD,
  file.path(
    OUT,
    "04_hard_checks.csv"
  )
)


# ============================================================
# 22. CONSOLE OUTPUT
# ============================================================

cat(
  "\n============================================================\n",
  "COVARIATE NONLINEARITY TESTS\n",
  "============================================================\n",
  sep = ""
)

print(
  NONLINEAR[
    order(
      P_BH
    )
  ],
  nrows =
    Inf
)


cat(
  "\n============================================================\n",
  "BIOMASS: CURRENT VS FLEXIBLE COVARIATES\n",
  "============================================================\n",
  sep = ""
)

print(
  COMPARE[
    fuel ==
      "biomass"
  ][
    order(
      metric
    )
  ],
  nrows =
    Inf
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
      )
    )
  ],
  nrows =
    Inf
)


cat(
  "\n============================================================\n",
  "HARD CHECKS\n",
  "============================================================\n",
  sep = ""
)

print(
  HARD,
  nrows =
    Inf
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
    "07H20 COVARIATE FUNCTIONAL-FORM SENSITIVITY: PASS\n",
    "07H20 COVARIATE FUNCTIONAL-FORM SENSITIVITY: FAIL\n"
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
    "07H20存在Hard FAIL。"
  )
}

