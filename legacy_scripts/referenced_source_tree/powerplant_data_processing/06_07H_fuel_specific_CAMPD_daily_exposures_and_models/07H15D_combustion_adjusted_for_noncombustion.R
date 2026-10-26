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
# 07H15D
# 燃烧组调整非燃烧 / 其他能源影响
#
# Primary combustion exposures:
#   coal / gas / oil / biomass
#
# Metric:
#   Count = actual-daily active CAMPD plant count
#   Heat  = actual-daily CAMPD Heat Input
#
# Non-combustion adjustment:
#   wind / solar / hydro / nuclear / geothermal
#
# 所有非燃烧调整变量：
#   同半径
#   同为90° upwind
#   Count
#
# Model ladder:
#
# M0_noncomb_matched
#   combustion exposure
#   + combustion exposure × wind speed
#
# M1_noncomb_adjusted
#   M0
#   + Wind Count + Wind Count × wind speed
#   + Solar Count + Solar Count × wind speed
#   + Hydro Count + Hydro Count × wind speed
#   + Nuclear Count + Nuclear Count × wind speed
#   + Geothermal Count + Geothermal Count × wind speed
#
# M0与M1使用完全相同complete-case样本。
#
# FE:
#   monitor_FE
#   monitor_year_FE
#
# representative winds:
#   1 m/s
#   sample Q25
#   sample median
#   sample mean
#   sample Q75
#
# Count effect:
#   per +1 active plant
#
# Heat effect:
#   per +10,000 MMBtu/day
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

F07_ROOT <- file.path(
  PP,
  "04_07F_正式统计分析数据构建"
)

FILTER_FILE <- file.path(
  F07_ROOT,
  "03_07F_filter_period_analysis_base_2001_2024.csv.gz"
)

QDIR <- file.path(
  F07_ROOT,
  "quasi_daily_yearly"
)

# 正式燃烧fuel-specific actual-daily exposure
COMB_DAILY_DIR <- file.path(
  HROOT,
  "fuel_daily_yearly"
)

# 已验证的非燃烧20–50 km方向性Count
NONCOMB_ROOT <- file.path(
  HROOT,
  "07H15C0_非燃烧能源90度上风向Count_20_50km"
)

NONCOMB_DAILY_DIR <- file.path(
  NONCOMB_ROOT,
  "daily_yearly"
)

NONCOMB_AUDIT_FILE <- file.path(
  NONCOMB_ROOT,
  "01_07H15C0_build_audit.csv"
)

OUT <- file.path(
  HROOT,
  "07H15D_燃烧组调整非燃烧组影响"
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

COMB_FUELS <- c(
  "coal",
  "gas",
  "oil",
  "biomass"
)

METRICS <- c(
  "count",
  "heat"
)

NONCOMB_FUELS <- c(
  "wind",
  "solar",
  "hydro",
  "nuclear",
  "geothermal"
)

COMB_LABEL <- c(
  coal = "Coal",
  gas = "Natural gas",
  oil = "Oil",
  biomass = "Biomass"
)

NONCOMB_LABEL <- c(
  wind = "Wind",
  solar = "Solar",
  hydro = "Hydro",
  nuclear = "Nuclear",
  geothermal = "Geothermal"
)

METRIC_LABEL <- c(
  count = "Active plant Count",
  heat = "Actual daily Heat Input"
)

FE_LABEL <- c(
  monitor_FE = "Monitor FE",
  monitor_year_FE = "Monitor × year FE"
)

MODEL_LABEL <- c(
  M0_noncomb_matched =
    "Unadjusted · matched sample",

  M1_noncomb_adjusted =
    "Adjusted for non-combustion energy"
)

EXPECTED_FILTER_ROWS <- 107407L
EXPECTED_FILTER_MONITORS <- 133L

HEAT_SCALE <- 10000


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


comb_col <- function(
  fuel,
  metric,
  radius
) {

  if (
    metric == "count"
  ) {

    return(
      sprintf(
        "campd_%s_upwind_n_plants_%dkm_90deg",
        fuel,
        radius
      )
    )
  }

  sprintf(
    "campd_%s_upwind_heat_observed_%dkm_90deg",
    fuel,
    radius
  )
}


noncomb_col <- function(
  fuel,
  radius
) {

  sprintf(
    "%s_upwind_n_plants_%dkm",
    fuel,
    radius
  )
}


safe_cor <- function(
  x,
  y
) {

  ok <-
    is.finite(x) &
    is.finite(y)

  if (
    sum(ok) <
      3L
  ) {

    return(
      NA_real_
    )
  }

  if (
    !is.finite(
      sd(
        x[ok]
      )
    ) ||
      !is.finite(
        sd(
          y[ok]
        )
      ) ||
      sd(
        x[ok]
      ) <=
        0 ||
      sd(
        y[ok]
      ) <=
        0
  ) {

    return(
      NA_real_
    )
  }

  cor(
    x[ok],
    y[ok]
  )
}


effect_from_fit <- function(
  fit,
  wind_ms
) {

  b <- coef(
    fit
  )

  bn <- names(
    b
  )

  int_name <- intersect(
    c(
      "exposure:vel",
      "vel:exposure"
    ),
    bn
  )

  if (
    !"exposure" %in%
      bn ||
      length(
        int_name
      ) !=
      1L
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

  int_name <- int_name[1]

  V <- as.matrix(
    vcov(
      fit
    )
  )

  est <-
    unname(
      b[
        "exposure"
      ]
    ) +
    wind_ms *
    unname(
      b[
        int_name
      ]
    )

  vv <-
    V[
      "exposure",
      "exposure"
    ] +
    wind_ms^2 *
    V[
      int_name,
      int_name
    ] +
    2 *
    wind_ms *
    V[
      "exposure",
      int_name
    ]

  se <- sqrt(
    max(
      as.numeric(
        vv
      ),
      0
    )
  )

  low <-
    est -
    1.96 *
    se

  high <-
    est +
    1.96 *
    se

  p <- if (
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

  data.table(
    estimate =
      est,

    se =
      se,

    ci95_low =
      low,

    ci95_high =
      high,

    p_wald =
      p
  )
}


# ============================================================
# 4. Model formula
# ============================================================

make_formula <- function(
  fe_structure,
  adjusted
) {

  rhs <- c(
    "exposure",
    "exposure:vel"
  )

  if (
    adjusted
  ) {

    for (
      fuel in NONCOMB_FUELS
    ) {

      a <- paste0(
        "adj_",
        fuel
      )

      rhs <- c(
        rhs,
        a,
        paste0(
          a,
          ":vel"
        )
      )
    }
  }

  rhs <- c(
    rhs,

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

    fe_part <- "loc_num"

  } else {

    fe_part <- "monitor_year_fe"
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
# 5. Verify retained model terms
# ============================================================

check_model_terms <- function(
  fit,
  adjusted
) {

  bn <- names(
    coef(
      fit
    )
  )

  main_interaction <- intersect(
    c(
      "exposure:vel",
      "vel:exposure"
    ),
    bn
  )

  if (
    !"exposure" %in%
      bn ||
      length(
        main_interaction
      ) !=
      1L
  ) {

    return(
      list(
        ok = FALSE,
        note =
          "primary exposure or primary exposure×wind dropped"
      )
    )
  }

  if (
    !adjusted
  ) {

    return(
      list(
        ok = TRUE,
        note = ""
      )
    )
  }

  missing_terms <- character()

  for (
    fuel in NONCOMB_FUELS
  ) {

    a <- paste0(
      "adj_",
      fuel
    )

    if (
      !a %in%
        bn
    ) {

      missing_terms <- c(
        missing_terms,
        a
      )
    }

    ints <- c(
      paste0(
        a,
        ":vel"
      ),
      paste0(
        "vel:",
        a
      )
    )

    if (
      !any(
        ints %in%
          bn
      )
    ) {

      missing_terms <- c(
        missing_terms,
        paste0(
          a,
          "×wind"
        )
      )
    }
  }

  if (
    length(
      missing_terms
    ) >
      0L
  ) {

    return(
      list(
        ok = FALSE,

        note =
          paste0(
            "non-combustion adjustment terms dropped: ",
            paste(
              missing_terms,
              collapse = ", "
            )
          )
      )
    )
  }

  list(
    ok = TRUE,
    note = ""
  )
}


# ============================================================
# 6. Fit one model
# ============================================================

fit_one <- function(
  d,
  fe_structure,
  adjusted
) {

  if (
    nrow(d) <
      1000L
  ) {

    return(
      list(
        status = "SKIP",
        note = "N < 1000",
        fit = NULL
      )
    )
  }

  if (
    !is.finite(
      sd(
        d$exposure
      )
    ) ||
      sd(
        d$exposure
      ) <=
      0
  ) {

    return(
      list(
        status = "SKIP",
        note = "primary combustion exposure has no variation",
        fit = NULL
      )
    )
  }

  fml <- make_formula(
    fe_structure,
    adjusted
  )

  fit <- tryCatch(
    feols(
      fml,
      data = d,
      vcov = ~loc_num,
      notes = FALSE,
      warn = FALSE
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
        status = "FAIL",
        note =
          conditionMessage(
            fit
          ),
        fit = NULL
      )
    )
  }

  chk <- check_model_terms(
    fit,
    adjusted
  )

  if (
    !chk$ok
  ) {

    return(
      list(
        status = "SKIP",
        note = chk$note,
        fit = fit
      )
    )
  }

  list(
    status = "PASS",
    note = "",
    fit = fit
  )
}


# ============================================================
# 7. Input audit
# ============================================================

cat(
  "============================================================\n",
  "07H15D INPUT AUDIT\n",
  "============================================================\n",
  sep = ""
)

if (
  !file.exists(
    FILTER_FILE
  )
) {

  stopf(
    "缺少07F filter-period：",
    FILTER_FILE
  )
}

if (
  !file.exists(
    NONCOMB_AUDIT_FILE
  )
) {

  stopf(
    "缺少07H15C0审计文件。"
  )
}


COMB_FILES <- file.path(
  COMB_DAILY_DIR,
  paste0(
    "07H_fuel_daily_exposure_",
    YEARS,
    ".csv.gz"
  )
)


NONCOMB_FILES <- file.path(
  NONCOMB_DAILY_DIR,
  paste0(
    "07H15C_noncombustion_upwind_count_20_50km_",
    YEARS,
    ".csv.gz"
  )
)


if (
  sum(
    file.exists(
      COMB_FILES
    )
  ) !=
    24L
) {

  stopf(
    "07H01 combustion daily不是24/24。"
  )
}


if (
  sum(
    file.exists(
      NONCOMB_FILES
    )
  ) !=
    24L
) {

  stopf(
    "07H15C0 non-combustion daily不是24/24。"
  )
}


A0 <- fread(
  NONCOMB_AUDIT_FILE,
  showProgress = FALSE
)


if (
  nrow(A0) !=
    24L ||
      sum(
        A0$mismatch_20km,
        na.rm = TRUE
      ) !=
      0L ||
      any(
        A0$status !=
          "PASS"
      )
) {

  stopf(
    "07H15C0没有通过20 km exact reproduction审计。"
  )
}


cat(
  "07H01 combustion daily：24/24 PASS\n",
  "07H15C0 non-combustion daily：24/24 PASS\n",
  "07H15C0 20 km mismatch：0 PASS\n",
  sep = ""
)


# ============================================================
# 8. Recover formal 133 monitor universe
# ============================================================

Q_LOCS <- rbindlist(
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
          "缺少07F quasi：",
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

      unique(
        z[
          flag_07f_primary_quasi_sample %in%
            TRUE,
          .(
            loc_num
          )
        ]
      )
    }
  ),
  use.names = TRUE,
  fill = TRUE
)


ELIGIBLE_LOCS <- sort(
  unique(
    Q_LOCS$loc_num
  )
)


if (
  length(
    ELIGIBLE_LOCS
  ) !=
    133L
) {

  stopf(
    "正式eligible monitor不是133；实际=",
    length(
      ELIGIBLE_LOCS
    )
  )
}


cat(
  "Formal eligible monitors = 133：PASS\n"
)


# ============================================================
# 9. Exact frozen formal filter-period base
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


FP_HEADER <- names(
  fread(
    FILTER_FILE,
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
  FILTER_FILE,
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


FP_NUMERIC <- c(
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
  cc in FP_NUMERIC
) {

  set(
    FP,
    j = cc,
    value =
      safe_num(
        FP[
          ,
          get(cc)
        ]
      )
  )
}


# Frozen 07G / 07H10 formal filter-period definition
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
    EXPECTED_FILTER_ROWS
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
    EXPECTED_FILTER_MONITORS
) {

  stopf(
    "正式filter-period不是133站；实际=",
    uniqueN(
      FP$loc_num
    )
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
      period_days <=
      0,
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
# 10. Expand [collect_start, collect_end)
# ============================================================

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
  "calendar_year",
  "monitor_year_fe"
)


# ============================================================
# 11. Output containers
# ============================================================

STATUS <- data.table()

EFFECTS <- data.table()

SUPPORT <- data.table()

CORRELATION <- data.table()

DAILY_AUDIT <- data.table()


# ============================================================
# 12. Radius loop
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


  # ----------------------------------------------------------
  # combustion columns
  # ----------------------------------------------------------

  COMB_COUNT_COLS <- vapply(
    COMB_FUELS,
    function(fuel) {

      comb_col(
        fuel,
        "count",
        radius
      )
    },
    character(1),
    USE.NAMES = FALSE
  )


  COMB_HEAT_COLS <- vapply(
    COMB_FUELS,
    function(fuel) {

      comb_col(
        fuel,
        "heat",
        radius
      )
    },
    character(1),
    USE.NAMES = FALSE
  )


  # ----------------------------------------------------------
  # non-combustion adjustment columns
  # ----------------------------------------------------------

  NONCOMB_COLS <- vapply(
    NONCOMB_FUELS,
    function(fuel) {

      noncomb_col(
        fuel,
        radius
      )
    },
    character(1),
    USE.NAMES = FALSE
  )


  # ==========================================================
  # 12A. Read combustion daily
  # ==========================================================

  COMB_DAILY <- rbindlist(
    lapply(
      COMB_FILES,
      function(f) {

        H <- names(
          fread(
            f,
            nrows = 0,
            showProgress = FALSE
          )
        )

        need <- c(
          "loc_num",
          "date",
          COMB_COUNT_COLS,
          COMB_HEAT_COLS
        )

        miss <- setdiff(
          need,
          H
        )

        if (
          length(
            miss
          ) >
            0L
        ) {

          stopf(
            basename(f),
            " 缺燃烧字段：",
            paste(
              miss,
              collapse = ", "
            )
          )
        }

        fread(
          f,
          select = need,
          showProgress = FALSE
        )
      }
    ),
    use.names = TRUE,
    fill = TRUE
  )


  COMB_DAILY[
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


  if (
    COMB_DAILY[
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
      " km combustion daily存在重复monitor-day。"
    )
  }


  # ==========================================================
  # 12B. Read non-combustion daily
  # ==========================================================

  NONCOMB_DAILY <- rbindlist(
    lapply(
      NONCOMB_FILES,
      function(f) {

        H <- names(
          fread(
            f,
            nrows = 0,
            showProgress = FALSE
          )
        )

        need <- c(
          "loc_num",
          "date",
          NONCOMB_COLS
        )

        miss <- setdiff(
          need,
          H
        )

        if (
          length(
            miss
          ) >
            0L
        ) {

          stopf(
            basename(f),
            " 缺非燃烧字段：",
            paste(
              miss,
              collapse = ", "
            )
          )
        }

        fread(
          f,
          select = need,
          showProgress = FALSE
        )
      }
    ),
    use.names = TRUE,
    fill = TRUE
  )


  NONCOMB_DAILY[
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

      noncomb_key_match =
        TRUE
    )
  ]


  if (
    NONCOMB_DAILY[
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
      " km non-combustion daily存在重复monitor-day。"
    )
  }


  # ==========================================================
  # 12C. Exact daily key merge
  # ==========================================================

  DAILY <- merge(
    COMB_DAILY,
    NONCOMB_DAILY,
    by = c(
      "loc_num",
      "date"
    ),
    all.x = TRUE,
    sort = FALSE
  )


  UNMATCHED_NONCOMB <- DAILY[
    is.na(
      noncomb_key_match
    ),
    .N
  ]


  if (
    UNMATCHED_NONCOMB >
      0L
  ) {

    stopf(
      radius,
      " km有 ",
      UNMATCHED_NONCOMB,
      " 个燃烧monitor-day无法连接非燃烧暴露。"
    )
  }


  DAILY[
    ,
    noncomb_key_match :=
      NULL
  ]


  DAILY_AUDIT <- rbind(
    DAILY_AUDIT,

    data.table(
      radius_km =
        radius,

      combustion_daily_rows =
        nrow(
          COMB_DAILY
        ),

      noncombustion_daily_rows =
        nrow(
          NONCOMB_DAILY
        ),

      merged_daily_rows =
        nrow(
          DAILY
        ),

      unmatched_noncombustion_keys =
        UNMATCHED_NONCOMB
    ),

    fill = TRUE
  )


  rm(
    COMB_DAILY,
    NONCOMB_DAILY
  )


  invisible(
    gc(
      verbose = FALSE
    )
  )


  # ==========================================================
  # 12D. Restrict to formal filter-period day map
  # ==========================================================

  DAILY[
    ,
    daily_match :=
      TRUE
  ]


  DD <- merge(
    FP_DAYS,
    DAILY,
    by = c(
      "loc_num",
      "date"
    ),
    all.x = TRUE,
    sort = FALSE
  )


  UNMATCHED_FP_DAYS <- DD[
    is.na(
      daily_match
    ),
    .N
  ]


  if (
    UNMATCHED_FP_DAYS >
      0L
  ) {

    stopf(
      radius,
      " km正式filter-period有 ",
      UNMATCHED_FP_DAYS,
      " 个calendar-day缺少daily exposure。"
    )
  }


  DD[
    ,
    daily_match :=
      NULL
  ]


  ALL_EXPOSURE_COLS <- c(
    COMB_COUNT_COLS,
    COMB_HEAT_COLS,
    NONCOMB_COLS
  )


  # ==========================================================
  # 12E. Strict filter-period mean
  # ==========================================================

  AG <- DD[
    ,
    lapply(
      .SD,
      strict_mean
    ),
    by =
      period_id,
    .SDcols =
      ALL_EXPOSURE_COLS
  ]


  PERIOD <- merge(
    FP,
    AG,
    by =
      "period_id",
    all.x =
      TRUE,
    sort =
      FALSE
  )


  if (
    nrow(
      PERIOD
    ) !=
      EXPECTED_FILTER_ROWS
  ) {

    stopf(
      radius,
      " km聚合后filter-period行数改变。"
    )
  }


  rm(
    DAILY,
    DD,
    AG
  )


  invisible(
    gc(
      verbose = FALSE
    )
  )


  # ==========================================================
  # 13. Primary combustion exposure loop
  # ==========================================================

  for (
    fuel in COMB_FUELS
  ) {

    for (
      metric in METRICS
    ) {

      cat(
        "\n",
        fuel,
        " | ",
        metric,
        " | ",
        radius,
        " km\n",
        sep = ""
      )


      primary_col <- comb_col(
        fuel,
        metric,
        radius
      )


      d0 <- PERIOD[
        ,
        c(
          MODEL_BASE_COLS,
          primary_col,
          NONCOMB_COLS
        ),
        with = FALSE
      ]


      # --------------------------------------------------------
      # Primary combustion exposure
      # --------------------------------------------------------

      primary_raw <- safe_num(
        d0[
          ,
          get(
            primary_col
          )
        ]
      )


      if (
        metric ==
          "heat"
      ) {

        primary_raw <-
          primary_raw /
          HEAT_SCALE
      }


      set(
        d0,
        j = "exposure",
        value = primary_raw
      )


      # --------------------------------------------------------
      # Non-combustion adjustment variables
      # --------------------------------------------------------

      for (
        ii in seq_along(
          NONCOMB_FUELS
        )
      ) {

        nf <- NONCOMB_FUELS[ii]

        set(
          d0,
          j = paste0(
            "adj_",
            nf
          ),
          value =
            safe_num(
              d0[
                ,
                get(
                  NONCOMB_COLS[ii]
                )
              ]
            )
        )
      }


      ADJ_NAMES <- paste0(
        "adj_",
        NONCOMB_FUELS
      )


      # ========================================================
      # EXACT matched complete-case sample
      # ========================================================

      NEEDED <- c(
        MODEL_BASE_COLS,
        "exposure",
        ADJ_NAMES
      )


      d <- d0[
        complete.cases(
          d0[
            ,
            ..NEEDED
          ]
        )
      ]


      # --------------------------------------------------------
      # Primary exposure support
      # --------------------------------------------------------

      WV <- d[
        ,
        .(
          exposure_var =
            if (
              .N >
                1L
            ) {

              var(
                exposure
              )

            } else {

              NA_real_
            }
        ),
        by =
          monitor_year_fe
      ]


      SUPPORT <- rbind(
        SUPPORT,

        data.table(
          fuel =
            fuel,

          metric =
            metric,

          radius_km =
            radius,

          n_rows =
            nrow(d),

          n_monitors =
            uniqueN(
              d$loc_num
            ),

          n_exposed_rows =
            sum(
              d$exposure >
                0,
              na.rm = TRUE
            ),

          n_exposed_monitors =
            uniqueN(
              d[
                exposure >
                  0,
                loc_num
              ]
            ),

          n_monitor_years_with_variation =
            WV[
              is.finite(
                exposure_var
              ) &
                exposure_var >
                0,
              .N
            ]
        ),

        fill = TRUE
      )


      # ========================================================
      # Correlations:
      # combustion exposure vs non-combustion adjusters
      # ========================================================

      for (
        nf in NONCOMB_FUELS
      ) {

        av <- paste0(
          "adj_",
          nf
        )


        raw_r <- safe_cor(
          d$exposure,
          d[
            ,
            get(av)
          ]
        )


        WT <- d[
          ,
          .(
            x_within =
              exposure -
              mean(
                exposure
              ),

            y_within =
              get(av) -
              mean(
                get(av)
              )
          ),
          by =
            monitor_year_fe
        ]


        within_r <- safe_cor(
          WT$x_within,
          WT$y_within
        )


        CORRELATION <- rbind(
          CORRELATION,

          data.table(
            combustion_fuel =
              fuel,

            metric =
              metric,

            radius_km =
              radius,

            noncombustion_fuel =
              nf,

            raw_correlation =
              raw_r,

            within_monitor_year_correlation =
              within_r,

            n_rows =
              nrow(d)
          ),

          fill = TRUE
        )
      }


      # ========================================================
      # Representative wind values
      # ========================================================

      q25 <- as.numeric(
        quantile(
          d$vel,
          0.25,
          na.rm = TRUE,
          names = FALSE
        )
      )


      med <- as.numeric(
        quantile(
          d$vel,
          0.50,
          na.rm = TRUE,
          names = FALSE
        )
      )


      mn <- mean(
        d$vel,
        na.rm = TRUE
      )


      q75 <- as.numeric(
        quantile(
          d$vel,
          0.75,
          na.rm = TRUE,
          names = FALSE
        )
      )


      WINDS <- data.table(
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


      # ========================================================
      # Fit M0 / M1 under two FE structures
      # ========================================================

      for (
        fe_structure in c(
          "monitor_FE",
          "monitor_year_FE"
        )
      ) {

        for (
          model_id in c(
            "M0_noncomb_matched",
            "M1_noncomb_adjusted"
          )
        ) {

          adjusted <-
            model_id ==
            "M1_noncomb_adjusted"


          cat(
            "[",
            fe_structure,
            "] ",
            model_id,
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
            fe_structure,
            adjusted
          )


          STATUS <- rbind(
            STATUS,

            data.table(
              fe_structure =
                fe_structure,

              model_id =
                model_id,

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

              note =
                rr$note
            ),

            fill = TRUE
          )


          if (
            rr$status !=
              "PASS"
          ) {

            next
          }


          for (
            wi in seq_len(
              nrow(
                WINDS
              )
            )
          ) {

            ee <- effect_from_fit(
              rr$fit,
              WINDS$wind_ms[wi]
            )


            EFFECTS <- rbind(
              EFFECTS,

              data.table(
                fe_structure =
                  fe_structure,

                model_id =
                  model_id,

                fuel =
                  fuel,

                metric =
                  metric,

                radius_km =
                  radius,

                wind_reference =
                  WINDS$wind_reference[wi],

                wind_ms =
                  WINDS$wind_ms[wi],

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
              ),

              fill = TRUE
            )
          }
        }
      }


      rm(
        d0,
        d
      )
    }
  }


  rm(
    PERIOD
  )


  invisible(
    gc(
      verbose = FALSE
    )
  )
}


# ============================================================
# 14. Save raw outputs
# ============================================================

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
    "02_effects_all_winds.csv"
  )
)


fwrite(
  SUPPORT,
  file.path(
    OUT,
    "03_combustion_exposure_support.csv"
  )
)


fwrite(
  CORRELATION,
  file.path(
    OUT,
    "04_combustion_noncombustion_correlations.csv"
  )
)


fwrite(
  DAILY_AUDIT,
  file.path(
    OUT,
    "05_daily_join_audit.csv"
  )
)


# ============================================================
# 15. Median wind before/after table
# ============================================================

MED <- EFFECTS[
  wind_reference ==
    "sample_median"
]


BASE <- MED[
  model_id ==
    "M0_noncomb_matched",
  .(
    fe_structure,
    fuel,
    metric,
    radius_km,
    wind_ms,

    baseline_estimate =
      estimate,

    baseline_ci95_low =
      ci95_low,

    baseline_ci95_high =
      ci95_high,

    baseline_p =
      p_wald
  )
]


ADJ <- MED[
  model_id ==
    "M1_noncomb_adjusted",
  .(
    fe_structure,
    fuel,
    metric,
    radius_km,

    adjusted_estimate =
      estimate,

    adjusted_ci95_low =
      ci95_low,

    adjusted_ci95_high =
      ci95_high,

    adjusted_p =
      p_wald
  )
]


COMPARE <- merge(
  BASE,
  ADJ,
  by = c(
    "fe_structure",
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
    estimate_change =
      adjusted_estimate -
      baseline_estimate,

    absolute_estimate_change =
      abs(
        adjusted_estimate -
          baseline_estimate
      ),

    adjusted_to_baseline_ratio =
      fifelse(
        is.finite(
          baseline_estimate
        ) &
          abs(
            baseline_estimate
          ) >
          1e-12,

        adjusted_estimate /
          baseline_estimate,

        NA_real_
      ),

    sign_changed =
      is.finite(
        baseline_estimate
      ) &
        is.finite(
          adjusted_estimate
        ) &
        sign(
          baseline_estimate
        ) !=
        sign(
          adjusted_estimate
        )
  )
]


fwrite(
  COMPARE,
  file.path(
    OUT,
    "06_median_wind_before_after_noncombustion_adjustment.csv"
  )
)


# ============================================================
# 16. Median-wind direction summary
# ============================================================

MED_SUMMARY <- MED[
  ,
  .(
    n_radii =
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
    fe_structure,
    model_id,
    fuel,
    metric
  )
]


fwrite(
  MED_SUMMARY,
  file.path(
    OUT,
    "07_median_wind_summary.csv"
  )
)


# ============================================================
# 17. Plot preparation
# ============================================================

decorate <- function(D) {

  D[
    ,
    fuel_label :=
      factor(
        COMB_LABEL[
          fuel
        ],
        levels =
          unname(
            COMB_LABEL[
              COMB_FUELS
            ]
          )
      )
  ]


  D[
    ,
    fe_label :=
      factor(
        FE_LABEL[
          fe_structure
        ],
        levels =
          unname(
            FE_LABEL[
              c(
                "monitor_FE",
                "monitor_year_FE"
              )
            ]
          )
      )
  ]


  D[
    ,
    model_label :=
      factor(
        MODEL_LABEL[
          model_id
        ],
        levels =
          unname(
            MODEL_LABEL[
              c(
                "M0_noncomb_matched",
                "M1_noncomb_adjusted"
              )
            ]
          )
      )
  ]


  D
}


theme_my <- theme_bw(
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

    legend.position =
      "bottom"
  )


# ============================================================
# 18. Figure S01: combustion Count
# ============================================================

D_COUNT <- decorate(
  copy(
    MED[
      metric ==
        "count"
    ]
  )
)


p1 <- ggplot(
  D_COUNT,
  aes(
    x = radius_km,
    y = estimate,
    group = fuel_label,
    color = fuel_label,
    linetype = fuel_label,
    shape = fuel_label
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
    width = 0.7,
    linewidth = 0.35
  ) +
  geom_line(
    linewidth = 0.7
  ) +
  geom_point(
    size = 1.8
  ) +
  facet_grid(
    fe_label ~ model_label,
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = RADII
  ) +
  labs(
    title =
      "Combustion Count associations before and after non-combustion adjustment",

    subtitle =
      "Matched samples; effects evaluated at sample median wind.",

    x =
      "Cumulative radius (km)",

    y =
      "Gross-beta change per additional active upwind plant (mBq/m³)",

    color =
      "Fuel",

    linetype =
      "Fuel",

    shape =
      "Fuel"
  ) +
  theme_my


ggsave(
  file.path(
    FIG_DIR,
    "S01_combustion_Count_before_after_noncombustion_adjustment.png"
  ),
  p1,
  width = 13,
  height = 8.5,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 19. Figure S02: combustion Heat
# ============================================================

D_HEAT <- decorate(
  copy(
    MED[
      metric ==
        "heat"
    ]
  )
)


p2 <- ggplot(
  D_HEAT,
  aes(
    x = radius_km,
    y = estimate,
    group = fuel_label,
    color = fuel_label,
    linetype = fuel_label,
    shape = fuel_label
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
    width = 0.7,
    linewidth = 0.35
  ) +
  geom_line(
    linewidth = 0.7
  ) +
  geom_point(
    size = 1.8
  ) +
  facet_grid(
    fe_label ~ model_label,
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = RADII
  ) +
  labs(
    title =
      "Combustion Heat associations before and after non-combustion adjustment",

    subtitle =
      "Matched samples; effects evaluated at sample median wind.",

    x =
      "Cumulative radius (km)",

    y =
      "Gross-beta change per 10,000 MMBtu/day upwind Heat Input (mBq/m³)",

    color =
      "Fuel",

    linetype =
      "Fuel",

    shape =
      "Fuel"
  ) +
  theme_my


ggsave(
  file.path(
    FIG_DIR,
    "S02_combustion_Heat_before_after_noncombustion_adjustment.png"
  ),
  p2,
  width = 13,
  height = 8.5,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 20. Figure S03:
# monitor×year FE effect shift
# ============================================================

SHIFT <- copy(
  COMPARE[
    fe_structure ==
      "monitor_year_FE"
  ]
)


SHIFT[
  ,
  fuel_label :=
    factor(
      COMB_LABEL[
        fuel
      ],
      levels =
        unname(
          COMB_LABEL[
            COMB_FUELS
          ]
        )
    )
]


SHIFT[
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


p3 <- ggplot(
  SHIFT,
  aes(
    x = radius_km,
    y = estimate_change,
    group = fuel_label,
    color = fuel_label,
    linetype = fuel_label,
    shape = fuel_label
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_line(
    linewidth = 0.7
  ) +
  geom_point(
    size = 1.8
  ) +
  facet_wrap(
    ~metric_label,
    ncol = 1,
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = RADII
  ) +
  labs(
    title =
      "Change in combustion effect after non-combustion adjustment",

    subtitle =
      "Monitor × year FE; adjusted estimate minus matched unadjusted estimate.",

    x =
      "Cumulative radius (km)",

    y =
      "Adjusted − unadjusted estimate",

    color =
      "Fuel",

    linetype =
      "Fuel",

    shape =
      "Fuel"
  ) +
  theme_my


ggsave(
  file.path(
    FIG_DIR,
    "S03_monitorYearFE_effect_shift_after_noncombustion_adjustment.png"
  ),
  p3,
  width = 11,
  height = 8.5,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 21. Figure S04:
# adjusted models only
# ============================================================

FINAL <- decorate(
  copy(
    MED[
      model_id ==
        "M1_noncomb_adjusted"
    ]
  )
)


FINAL[
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


p4 <- ggplot(
  FINAL,
  aes(
    x = radius_km,
    y = estimate,
    group = fuel_label,
    color = fuel_label,
    linetype = fuel_label,
    shape = fuel_label
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
    width = 0.7,
    linewidth = 0.35
  ) +
  geom_line(
    linewidth = 0.7
  ) +
  geom_point(
    size = 1.8
  ) +
  facet_grid(
    fe_label ~ metric_label,
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = RADII
  ) +
  labs(
    title =
      "Combustion associations after adjustment for non-combustion energy",

    subtitle =
      "Wind, solar, hydro, nuclear and geothermal upwind Count simultaneously adjusted.",

    x =
      "Cumulative radius (km)",

    y =
      "Adjusted gross-beta effect",

    color =
      "Fuel",

    linetype =
      "Fuel",

    shape =
      "Fuel"
  ) +
  theme_my


ggsave(
  file.path(
    FIG_DIR,
    "S04_adjusted_combustion_results_only.png"
  ),
  p4,
  width = 13,
  height = 8.5,
  dpi = 300,
  bg = "white"
)


# ============================================================
# 22. Hard checks
# ============================================================

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


N_FAIL <- STATUS[
  status ==
    "FAIL",
  .N
]


# 4 fuels × 2 metrics × 7 radii × 2 FE × 2 models
EXPECTED_MODEL_ROWS <-
  length(
    COMB_FUELS
  ) *
  length(
    METRICS
  ) *
  length(
    RADII
  ) *
  2L *
  2L


# M0 / M1 within one FE + fuel + metric + radius
# must use exactly same sample size
PAIR_AUDIT <- STATUS[
  ,
  .(
    n_unique_sample_sizes =
      uniqueN(
        n_rows
      )
  ),
  by = .(
    fe_structure,
    fuel,
    metric,
    radius_km
  )
]


N_SAMPLE_MISMATCH <- PAIR_AUDIT[
  n_unique_sample_sizes !=
    1L,
  .N
]


EXPECTED_PLOTS <- file.path(
  FIG_DIR,
  c(
    "S01_combustion_Count_before_after_noncombustion_adjustment.png",
    "S02_combustion_Heat_before_after_noncombustion_adjustment.png",
    "S03_monitorYearFE_effect_shift_after_noncombustion_adjustment.png",
    "S04_adjusted_combustion_results_only.png"
  )
)


HARD <- data.table(
  check = c(
    "H01_combustion_daily_files_24",
    "H15C0_noncombustion_daily_files_24",
    "H15C0_20km_mismatch_zero",
    "formal_filter_rows_107407",
    "formal_filter_monitors_133",
    "daily_join_unmatched_zero",
    "model_status_rows_224",
    "model_FAIL_zero",
    "model_SKIP_zero",
    "matched_pair_sample_size_equal",
    "effect_rows_equal_PASS_x5",
    "four_png"
  ),

  observed = c(
    as.character(
      sum(
        file.exists(
          COMB_FILES
        )
      )
    ),

    as.character(
      sum(
        file.exists(
          NONCOMB_FILES
        )
      )
    ),

    as.character(
      sum(
        A0$mismatch_20km,
        na.rm = TRUE
      )
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
      sum(
        DAILY_AUDIT$unmatched_noncombustion_keys
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
      N_SAMPLE_MISMATCH
    ),

    as.character(
      nrow(
        EFFECTS
      )
    ),

    as.character(
      sum(
        file.exists(
          EXPECTED_PLOTS
        )
      )
    )
  ),

  expected = c(
    "24",
    "24",
    "0",
    "107407",
    "133",
    "0",
    as.character(
      EXPECTED_MODEL_ROWS
    ),
    "0",
    "0",
    "0",
    as.character(
      N_PASS *
        5L
    ),
    "4"
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
    "08_hard_checks.csv"
  )
)


# ============================================================
# 23. Terminal summaries
# ============================================================

cat(
  "\n============================================================\n",
  "MODEL STATUS\n",
  "============================================================\n",
  sep = ""
)


print(
  STATUS[
    ,
    .N,
    by = .(
      model_id,
      fe_structure,
      status
    )
  ][
    order(
      model_id,
      fe_structure,
      status
    )
  ],
  nrows = Inf
)


cat(
  "\n============================================================\n",
  "MEDIAN-WIND SUMMARY\n",
  "============================================================\n",
  sep = ""
)


print(
  MED_SUMMARY[
    order(
      fe_structure,
      metric,
      model_id,
      fuel
    )
  ],
  nrows = Inf
)


cat(
  "\n============================================================\n",
  "MONITOR×YEAR FE：BEFORE / AFTER NON-COMBUSTION ADJUSTMENT\n",
  "============================================================\n",
  sep = ""
)


print(
  COMPARE[
    fe_structure ==
      "monitor_year_FE",
    .(
      fuel,
      metric,
      radius_km,
      wind_ms,
      baseline_estimate,
      baseline_ci95_low,
      baseline_ci95_high,
      adjusted_estimate,
      adjusted_ci95_low,
      adjusted_ci95_high,
      estimate_change,
      adjusted_to_baseline_ratio,
      sign_changed
    )
  ][
    order(
      metric,
      fuel,
      radius_km
    )
  ],
  nrows = Inf
)


# ============================================================
# Correlation summary
# ============================================================

COR_SUM <- CORRELATION[
  ,
  .(
    min_within_r =
      min(
        within_monitor_year_correlation,
        na.rm = TRUE
      ),

    max_within_r =
      max(
        within_monitor_year_correlation,
        na.rm = TRUE
      ),

    max_abs_within_r =
      max(
        abs(
          within_monitor_year_correlation
        ),
        na.rm = TRUE
      )
  ),
  by = .(
    combustion_fuel,
    metric,
    noncombustion_fuel
  )
]


cat(
  "\n============================================================\n",
  "WITHIN MONITOR×YEAR CORRELATION RANGE\n",
  "============================================================\n",
  sep = ""
)


print(
  COR_SUM[
    order(
      metric,
      combustion_fuel,
      noncombustion_fuel
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

    "07H15D COMBUSTION ADJUSTED FOR NON-COMBUSTION: PASS\n",

    "07H15D COMBUSTION ADJUSTED FOR NON-COMBUSTION: FAIL\n"
  ),

  "PASS models = ",
  N_PASS,
  "\n",

  "SKIP models = ",
  N_SKIP,
  "\n",

  "FAIL models = ",
  N_FAIL,
  "\n",

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
    "07H15D存在Hard FAIL。"
  )
}

