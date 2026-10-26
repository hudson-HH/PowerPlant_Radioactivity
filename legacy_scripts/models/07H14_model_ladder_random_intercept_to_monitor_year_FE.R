#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1, scipen = 999, width = 220)

suppressPackageStartupMessages({
  library(data.table)
  library(lme4)
  library(fixest)
  library(ggplot2)
})

# ============================================================
# 07H14 模型阶梯诊断
#
# fuels:
#   biomass / coal / gas / oil
#
# radii:
#   20,25,30,35,40,45,50 km
#
# exposures:
#   daily active plant Count
#   actual daily Heat Input
#
# geometry:
#   360° non-directional
#   90° upwind
#
# model ladder:
#   M1 quasi-daily + monitor random intercept
#   M2 filter-period + monitor random intercept
#   M3 filter-period + monitor fixed effect
#   M4 filter-period + monitor×year fixed effect
#
# 所有模型均：
#   Exposure × wind speed
#   输出 1 m/s 与样本中位风速效应
#
# 本阶段只诊断，不修改任何冻结上游。
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

FILTER_FILE <- file.path(
  F07,
  "03_07F_filter_period_analysis_base_2001_2024.csv.gz"
)

RAW_DAILY_DIR <- file.path(
  HROOT,
  "fuel_daily_yearly"
)

GEO_FILE <- file.path(
  PP,
  "00_既往发电站空间暴露与模型流程_冻结",
  "02_07A_RadNet发电站空间链接",
  "07A_radnet_powerplant_same_year_250km",
  "03_monitor_year_all_plant_link_250km_07A_2001_2024.csv"
)

PF_DIR <- file.path(
  PP,
  "01_发电站月度逐日源强增强流程_冻结",
  "07B2_CAMPD_daily_emissions",
  "05_plant_day_primary_fuel"
)

HARM_FILE <- file.path(
  PP,
  "01_发电站月度逐日源强增强流程_冻结",
  "07B3_EIA_CAMPD_harmonization",
  "05_EIA_CAMPD_plant_month_harmonization_2001_2024.csv.gz"
)

OUT <- file.path(
  HROOT,
  "07H14_模型阶梯诊断_论文式随机截距到站点年份固定效应"
)

ND_DAILY_DIR <- file.path(
  OUT,
  "non_directional_daily_yearly"
)

CHECK_DIR <- file.path(
  OUT,
  "checkpoints"
)

FIG_DIR <- file.path(
  OUT,
  "figures"
)

AUDIT_DIR <- file.path(
  OUT,
  "audit"
)

for (
  d in c(
    OUT,
    ND_DAILY_DIR,
    CHECK_DIR,
    FIG_DIR,
    AUDIT_DIR
  )
) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ============================================================
# 2. Parameters
# ============================================================

YEARS <- 2001:2024

FUELS <- c(
  "biomass",
  "coal",
  "gas",
  "oil"
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

METRICS <- c(
  "count",
  "heat"
)

MAX_RADIUS <- max(
  RADII
)

HEAT_SCALE <- 10000

EXPECTED_Q_ROWS <- 426657L
EXPECTED_Q_MONITORS <- 133L

EXPECTED_FP_ROWS <- 107407L
EXPECTED_FP_MONITORS <- 133L

CHECKPOINT_VERSION <- "07H14_v1_20260922"

FORCE_REBUILD_NONDIR <- FALSE
FORCE_REFIT <- FALSE

FUEL_LABEL <- c(
  biomass = "Biomass",
  coal = "Coal",
  gas = "Natural gas",
  oil = "Oil"
)

MODEL_LABEL <- c(
  M1 = "M1 Quasi-daily + monitor random intercept",
  M2 = "M2 Filter-period + monitor random intercept",
  M3 = "M3 Filter-period + monitor fixed effect",
  M4 = "M4 Filter-period + monitor×year fixed effect"
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

safe_int <- function(x) {
  suppressWarnings(
    as.integer(x)
  )
}

clean_chr <- function(x) {

  y <- trimws(
    as.character(x)
  )

  y[
    y %in%
      c(
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

as_flag <- function(
  x,
  default = FALSE
) {

  if (is.logical(x)) {

    y <- x

  } else if (
    is.numeric(x) ||
    is.integer(x)
  ) {

    y <- x != 0

  } else {

    xx <- tolower(
      trimws(
        as.character(x)
      )
    )

    y <- xx %in%
      c(
        "true",
        "t",
        "1",
        "yes",
        "y"
      )
  }

  y[
    is.na(y)
  ] <- default

  y
}

pick_col <- function(
  nms,
  candidates,
  label,
  required = TRUE
) {

  hit <- candidates[
    candidates %in%
      nms
  ]

  if (
    length(hit) ==
      0L
  ) {

    if (required) {

      stopf(
        "找不到",
        label,
        "字段。候选：",
        paste(
          candidates,
          collapse = " | "
        )
      )

    } else {

      return(
        NA_character_
      )
    }
  }

  hit[1]
}

sum_keep_all_na <- function(x) {

  if (
    length(x) ==
      0L ||
    all(
      !is.finite(x)
    )
  ) {

    return(
      NA_real_
    )
  }

  sum(
    x[
      is.finite(x)
    ]
  )
}

strict_mean <- function(
  x,
  expected_n
) {

  if (
    length(x) !=
      expected_n ||
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

collapse_messages <- function(x) {

  if (
    is.null(x) ||
    length(x) ==
      0L
  ) {

    return("")
  }

  z <- unique(
    trimws(
      as.character(x)
    )
  )

  z <- z[
    nzchar(z)
  ]

  paste(
    z,
    collapse = " | "
  )
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
    metric ==
      "count"
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
    metric ==
      "heat"
  ) {
    HEAT_SCALE
  } else {
    1
  }
}

metric_unit <- function(metric) {

  if (
    metric ==
      "heat"
  ) {

    "per 10,000 MMBtu/day"

  } else {

    "per additional active plant"
  }
}

find_interaction <- function(
  coef_names,
  a = "exposure",
  b = "vel"
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
    length(z) !=
      1L
  ) {

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
    length(weights) ==
      0L ||
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
  wind
) {

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
        "exposure",
        int_term
      )
    )
  )
}

make_lmer_control <- function() {

  z <- list(
    optimizer = "bobyqa",

    optCtrl = list(
      maxfun = 200000
    ),

    check.conv.singular =
      "ignore"
  )

  fn <- names(
    formals(
      lme4::lmerControl
    )
  )

  if (
    "check.scaleX" %in%
      fn
  ) {

    z$check.scaleX <-
      "ignore"
  }

  do.call(
    lme4::lmerControl,
    z
  )
}


# ============================================================
# 4. Input audit
# ============================================================

cat(
  "============================================================\n"
)

cat(
  "07H14 INPUT AUDIT\n"
)

cat(
  "============================================================\n"
)

for (
  f in c(
    GEO_FILE,
    HARM_FILE,
    FILTER_FILE
  )
) {

  if (
    !file.exists(f)
  ) {

    stopf(
      "缺少输入：\n",
      f
    )
  }
}

for (
  d in c(
    QDIR,
    RAW_DAILY_DIR,
    PF_DIR
  )
) {

  if (
    !dir.exists(d)
  ) {

    stopf(
      "缺少目录：\n",
      d
    )
  }
}

UP_FILES <- file.path(
  RAW_DAILY_DIR,
  paste0(
    "07H_fuel_daily_exposure_",
    YEARS,
    ".csv.gz"
  )
)

if (
  sum(
    file.exists(
      UP_FILES
    )
  ) !=
    length(YEARS)
) {

  stopf(
    "07H01逐日燃料暴露不是24/24。"
  )
}

PF_EXPECTED <- unlist(
  lapply(
    YEARS,
    function(yy) {

      file.path(
        PF_DIR,
        paste0(
          "CAMPD_plant_day_primary_fuel_",
          yy,
          "_q",
          1:4,
          ".csv.gz"
        )
      )
    }
  )
)

if (
  sum(
    file.exists(
      PF_EXPECTED
    )
  ) !=
    96L
) {

  stopf(
    "07B2 primary-fuel季度文件不是96/96。"
  )
}

cat(
  "07A geometry：PASS\n",
  "07B2 primary-fuel：96/96 PASS\n",
  "07B3 harmonization：PASS\n",
  "07F quasi/filter：PASS\n",
  "07H01 upwind daily：24/24 PASS\n",
  sep = ""
)


# ============================================================
# 5. Monitor crosswalk + geometry
# ============================================================

cat(
  "\n============================================================\n",
  "读取07F monitor crosswalk + 07A geometry\n",
  "============================================================\n",
  sep = ""
)

cross_parts <- vector(
  "list",
  length(YEARS)
)

for (
  ii in
  seq_along(YEARS)
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

  H <- names(
    fread(
      qf,
      nrows = 0,
      showProgress = FALSE
    )
  )

  miss <- setdiff(
    c(
      "loc_num",
      "radnet_monitor_id"
    ),
    H
  )

  if (
    length(miss) >
      0L
  ) {

    stopf(
      yy,
      "年07F缺monitor crosswalk字段：",
      paste(
        miss,
        collapse = ", "
      )
    )
  }

  z <- fread(
    qf,
    select = c(
      "loc_num",
      "radnet_monitor_id"
    ),
    showProgress = FALSE
  )

  cross_parts[[ii]] <-
    unique(
      z[
        ,
        .(
          loc_num =
            clean_chr(
              loc_num
            ),

          radnet_monitor_id =
            as.character(
              radnet_monitor_id
            )
        )
      ]
    )
}

CROSS <- unique(
  rbindlist(
    cross_parts,
    use.names = TRUE,
    fill = TRUE
  )
)

CROSS <- CROSS[
  !is.na(loc_num) &
  !is.na(
    radnet_monitor_id
  )
]

if (
  CROSS[
    ,
    .N,
    by = loc_num
  ][
    N >
      1L,
    .N
  ] >
    0L ||
  CROSS[
    ,
    .N,
    by = radnet_monitor_id
  ][
    N >
      1L,
    .N
  ] >
    0L
) {

  stopf(
    "loc_num <-> radnet_monitor_id不是一一对应。"
  )
}


# ------------------------------------------------------------
# Geometry
# ------------------------------------------------------------

GH <- names(
  fread(
    GEO_FILE,
    nrows = 0,
    showProgress = FALSE
  )
)

year_col <- pick_col(
  GH,
  c(
    "year",
    "calendar_year"
  ),
  "year"
)

plant_col <- pick_col(
  GH,
  c(
    "plant_id",
    "orispl",
    "orispl_code"
  ),
  "plant_id"
)

monitor_col <- pick_col(
  GH,
  c(
    "radnet_monitor_id"
  ),
  "radnet_monitor_id"
)

loc_col <- pick_col(
  GH,
  c(
    "loc_num"
  ),
  "loc_num",
  required = FALSE
)

dist_col <- pick_col(
  GH,
  c(
    "distance_km",
    "distance"
  ),
  "distance_km"
)

need_geo <- unique(
  na.omit(
    c(
      year_col,
      plant_col,
      monitor_col,
      loc_col,
      dist_col
    )
  )
)

G0 <- fread(
  GEO_FILE,
  select = need_geo,
  showProgress = TRUE
)

if (
  year_col !=
    "year"
) {

  setnames(
    G0,
    year_col,
    "year"
  )
}

if (
  plant_col !=
    "plant_id"
) {

  setnames(
    G0,
    plant_col,
    "plant_id"
  )
}

if (
  monitor_col !=
    "radnet_monitor_id"
) {

  setnames(
    G0,
    monitor_col,
    "radnet_monitor_id"
  )
}

if (
  dist_col !=
    "distance_km"
) {

  setnames(
    G0,
    dist_col,
    "distance_km"
  )
}

G0[
  ,
  `:=`(
    year =
      safe_int(
        year
      ),

    plant_id =
      clean_chr(
        plant_id
      ),

    radnet_monitor_id =
      as.character(
        radnet_monitor_id
      ),

    distance_km =
      safe_num(
        distance_km
      )
  )
]

if (
  !is.na(
    loc_col
  )
) {

  if (
    loc_col !=
      "loc_num"
  ) {

    setnames(
      G0,
      loc_col,
      "loc_num"
    )
  }

  G0[
    ,
    loc_num :=
      clean_chr(
        loc_num
      )
  ]

} else {

  G0 <- merge(
    G0,
    CROSS,
    by =
      "radnet_monitor_id",
    all.x =
      TRUE,
    sort =
      FALSE
  )
}

if (
  G0[
    is.na(
      loc_num
    ),
    .N
  ] >
    0L
) {

  stopf(
    "07A geometry存在无法映射到loc_num的RadNet ID。"
  )
}

GEO <- unique(
  G0[
    year %in%
      YEARS &
    !is.na(loc_num) &
    !is.na(plant_id) &
    plant_id !=
      "99999" &
    is.finite(
      distance_km
    ) &
    distance_km >=
      0 &
    distance_km <=
      MAX_RADIUS +
      1e-8,
    .(
      loc_num,
      year,
      plant_id,
      distance_km
    )
  ]
)

rm(
  G0
)

invisible(
  gc()
)

if (
  GEO[
    ,
    .N,
    by =
      .(
        loc_num,
        year,
        plant_id
      )
  ][
    N >
      1L,
    .N
  ] >
    0L
) {

  stopf(
    "07A geometry存在重复loc_num+year+plant_id。"
  )
}

cat(
  "<=50 km geometry rows = ",
  format(
    nrow(GEO),
    big.mark = ","
  ),
  "\n",
  sep = ""
)


# ============================================================
# 6. Preferred CAMPD plant-month + fuel mapping
# ============================================================

HARM <- fread(
  HARM_FILE,
  select = c(
    "plant_id",
    "year",
    "month",
    "preferred_source_resolution"
  ),
  showProgress = FALSE
)

HARM[
  ,
  `:=`(
    plant_id =
      clean_chr(
        plant_id
      ),

    year =
      safe_int(
        year
      ),

    month =
      safe_int(
        month
      ),

    preferred_source_resolution =
      clean_chr(
        preferred_source_resolution
      )
  )
]

if (
  !"CAMPD_DAILY_ACTUAL" %in%
    unique(
      HARM$preferred_source_resolution
    )
) {

  stopf(
    "07B3中不存在CAMPD_DAILY_ACTUAL。"
  )
}

ELIGIBLE_PM <- unique(
  HARM[
    preferred_source_resolution ==
      "CAMPD_DAILY_ACTUAL",
    .(
      plant_id,
      year,
      month
    )
  ]
)

rm(
  HARM
)

invisible(
  gc()
)


# ------------------------------------------------------------
# Same CAMPD fuel mapping as frozen 07H01
# ------------------------------------------------------------

CAMPD_COMPONENT_MAP <- c(

  COAL =
    "coal",

  `COAL REFUSE` =
    "coal",

  `PIPELINE NATURAL GAS` =
    "gas",

  `NATURAL GAS` =
    "gas",

  `PROCESS GAS` =
    "gas",

  `OTHER GAS` =
    "gas",

  `RESIDUAL OIL` =
    "oil",

  `DIESEL OIL` =
    "oil",

  `OTHER OIL` =
    "oil",

  `PETROLEUM COKE` =
    "oil",

  WOOD =
    "biomass",

  `TIRE DERIVED FUEL` =
    "solid_waste",

  `TIRE-DERIVED FUEL` =
    "solid_waste",

  `OTHER SOLID FUEL` =
    "other_fossil"
)

classify_campd_fuel_label <- function(x) {

  if (
    is.na(x) ||
    !nzchar(
      trimws(
        as.character(x)
      )
    )
  ) {

    return(
      data.table(
        fuel_group =
          "unknown",

        mapping_method =
          "blank_unknown"
      )
    )
  }

  parts <- trimws(
    strsplit(
      as.character(x),
      ",",
      fixed = TRUE
    )[[1]]
  )

  parts <- parts[
    nzchar(parts)
  ]

  keys <- toupper(
    parts
  )

  groups <- unname(
    CAMPD_COMPONENT_MAP[
      keys
    ]
  )

  if (
    any(
      is.na(groups)
    )
  ) {

    return(
      data.table(
        fuel_group =
          "unknown",

        mapping_method =
          "unmapped_component"
      )
    )
  }

  ug <- sort(
    unique(
      groups
    )
  )

  if (
    length(ug) ==
      1L
  ) {

    return(
      data.table(
        fuel_group =
          ug[1],

        mapping_method =
          ifelse(
            length(parts) ==
              1L,
            "single_component",
            "same_group_multi_component"
          )
      )
    )
  }

  data.table(
    fuel_group =
      "mixed",

    mapping_method =
      "cross_group_mixed"
  )
}


# ------------------------------------------------------------
# Read one year actual CAMPD activity
# ------------------------------------------------------------

read_activity_year <- function(
  yy
) {

  files <- file.path(
    PF_DIR,
    paste0(
      "CAMPD_plant_day_primary_fuel_",
      yy,
      "_q",
      1:4,
      ".csv.gz"
    )
  )

  parts <- lapply(
    files,
    function(f) {

      H <- names(
        fread(
          f,
          nrows = 0,
          showProgress = FALSE
        )
      )

      req <- c(
        "year",
        "month",
        "date",
        "plant_id",
        "primary_fuel_info",
        "daily_heat_input_mmbtu_observed"
      )

      miss <- setdiff(
        req,
        H
      )

      if (
        length(miss) >
          0L
      ) {

        stopf(
          "07B2 primary-fuel缺字段：",
          paste(
            miss,
            collapse = ", "
          )
        )
      }

      fread(
        f,
        select = req,
        showProgress = FALSE
      )
    }
  )

  D <- rbindlist(
    parts,
    use.names = TRUE,
    fill = TRUE
  )

  D[
    ,
    `:=`(
      year =
        safe_int(
          year
        ),

      month =
        safe_int(
          month
        ),

      date =
        as.IDate(
          date
        ),

      plant_id =
        clean_chr(
          plant_id
        ),

      primary_fuel_info =
        clean_chr(
          primary_fuel_info
        ),

      heat =
        safe_num(
          daily_heat_input_mmbtu_observed
        )
    )
  ]

  EL <- ELIGIBLE_PM[
    year ==
      yy
  ]

  D <- D[
    EL,
    on =
      .(
        plant_id,
        year,
        month
      ),
    nomatch =
      0
  ]

  relevant <- unique(
    GEO[
      year ==
        yy,
      .(
        plant_id,
        year
      )
    ]
  )

  D <- D[
    relevant,
    on =
      .(
        plant_id,
        year
      ),
    nomatch =
      0
  ]

  RAW <- unique(
    D[
      ,
      .(
        primary_fuel_info
      )
    ]
  )

  MAP <- rbindlist(
    lapply(
      RAW$primary_fuel_info,
      classify_campd_fuel_label
    ),
    use.names = TRUE,
    fill = TRUE
  )

  RAW <- cbind(
    RAW,
    MAP
  )

  D <- merge(
    D,
    RAW,
    by =
      "primary_fuel_info",
    all.x =
      TRUE,
    sort =
      FALSE
  )

  if (
    D[
      mapping_method ==
        "unmapped_component",
      .N
    ] >
      0L
  ) {

    stopf(
      yy,
      "年存在未映射CAMPD燃料标签。"
    )
  }

  A <- D[
    fuel_group %in%
      FUELS,
    .(
      heat =
        sum_keep_all_na(
          heat
        )
    ),
    by =
      .(
        plant_id,
        date,
        fuel_group
      )
  ]

  if (
    A[
      ,
      .N,
      by =
        .(
          plant_id,
          date,
          fuel_group
        )
    ][
      N >
        1L,
      .N
    ] >
      0L
  ) {

    stopf(
      yy,
      "年actual fuel activity仍有重复。"
    )
  }

  A
}


# ============================================================
# 7. Formal quasi-daily sample
# ============================================================

cat(
  "\n============================================================\n",
  "读取正式quasi-daily样本\n",
  "============================================================\n",
  sep = ""
)

Q_REQUIRED <- c(
  "loc_num",
  "date",
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
  "temp",
  "flag_07f_primary_quasi_sample"
)

Q_LIST <- vector(
  "list",
  length(YEARS)
)

for (
  ii in
  seq_along(YEARS)
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

  H <- names(
    fread(
      qf,
      nrows = 0,
      showProgress = FALSE
    )
  )

  miss <- setdiff(
    Q_REQUIRED,
    H
  )

  if (
    length(miss) >
      0L
  ) {

    stopf(
      yy,
      "年quasi-daily缺字段：",
      paste(
        miss,
        collapse = ", "
      )
    )
  }

  z <- fread(
    qf,
    select = Q_REQUIRED,
    showProgress = FALSE
  )

  z[
    ,
    flag_07f_primary_quasi_sample :=
      as_flag(
        flag_07f_primary_quasi_sample
      )
  ]

  Q_LIST[[ii]] <- z[
    flag_07f_primary_quasi_sample %in%
      TRUE
  ]
}

Q_BASE <- rbindlist(
  Q_LIST,
  use.names = TRUE,
  fill = TRUE
)

rm(
  Q_LIST
)

invisible(
  gc()
)

Q_BASE[
  ,
  `:=`(
    loc_num =
      clean_chr(
        loc_num
      ),

    date =
      as.IDate(
        date
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

    pm_centered =
      safe_num(
        pm_same_day
      ) -
      10,

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

    temp_c =
      safe_num(
        temp
      ) -
      273.15
  )
]

Q_BASE[
  ,
  calendar_year :=
    safe_int(
      format(
        date,
        "%Y"
      )
    )
]

Q_BASE[
  ,
  year_c :=
    calendar_year -
    2012.5
]

Q_BASE[
  ,
  monitor_year_fe :=
    interaction(
      loc_num,
      calendar_year,
      drop = TRUE,
      lex.order = TRUE
    )
]

Q_BASE[
  ,
  row_id :=
    .I
]

if (
  nrow(Q_BASE) !=
    EXPECTED_Q_ROWS
) {

  stopf(
    "quasi-daily正式样本行数=",
    nrow(Q_BASE),
    "，预期=",
    EXPECTED_Q_ROWS
  )
}

if (
  uniqueN(
    Q_BASE$loc_num
  ) !=
    EXPECTED_Q_MONITORS
) {

  stopf(
    "quasi-daily监测站数错误。"
  )
}

ELIGIBLE_LOCS <- sort(
  unique(
    Q_BASE$loc_num
  )
)

cat(
  "quasi-daily = ",
  format(
    nrow(Q_BASE),
    big.mark = ","
  ),
  " rows / ",
  uniqueN(
    Q_BASE$loc_num
  ),
  " monitors\n",
  sep = ""
)


# ============================================================
# 8. Formal filter-period sample
# ============================================================

cat(
  "\n============================================================\n",
  "读取正式filter-period样本\n",
  "============================================================\n",
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
  "hpbl",
  "inv_hpbl",
  "rhum",
  "soilm",
  "vel",
  "lat_diff_35",
  "temp",
  "flag_conus",
  "fpmean_campd_upwind_heat_observed_20km_90deg"
)

HFP <- names(
  fread(
    FILTER_FILE,
    nrows = 0,
    showProgress = FALSE
  )
)

miss <- setdiff(
  FP_REQUIRED,
  HFP
)

if (
  length(miss) >
    0L
) {

  stopf(
    "filter-period缺字段：",
    paste(
      miss,
      collapse = ", "
    )
  )
}

F <- fread(
  FILTER_FILE,
  select = FP_REQUIRED,
  showProgress = TRUE
)

if (
  nrow(F) !=
    174200L
) {

  stopf(
    "filter-period原始行数不是174200。"
  )
}

F[
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

    hpbl =
      safe_num(
        hpbl
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

    temp_c =
      safe_num(
        temp
      ) -
      273.15,

    flag_conus =
      as_flag(
        flag_conus
      ),

    fpmean_campd_upwind_heat_observed_20km_90deg =
      safe_num(
        fpmean_campd_upwind_heat_observed_20km_90deg
      )
  )
]

F[
  ,
  formal_filter :=
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

FP_BASE <- F[
  formal_filter %in%
    TRUE
]

rm(
  F
)

invisible(
  gc()
)

FP_BASE[
  ,
  `:=`(
    period_id =
      .I,

    n_period_days =
      as.integer(
        collect_end -
          collect_start
      ),

    pm_centered =
      pm_filter_period_official -
      10
  )
]

FP_BASE[
  ,
  calendar_year :=
    safe_int(
      format(
        collect_start,
        "%Y"
      )
    )
]

FP_BASE[
  ,
  year_c :=
    calendar_year -
    2012.5
]

FP_BASE[
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
  FP_BASE[
    !is.finite(
      n_period_days
    ) |
    n_period_days <=
      0L,
    .N
  ] >
    0L
) {

  stopf(
    "存在非法filter-period长度。"
  )
}

if (
  nrow(FP_BASE) !=
    EXPECTED_FP_ROWS
) {

  stopf(
    "filter-period正式样本行数=",
    nrow(FP_BASE),
    "，预期=",
    EXPECTED_FP_ROWS
  )
}

if (
  uniqueN(
    FP_BASE$loc_num
  ) !=
    EXPECTED_FP_MONITORS
) {

  stopf(
    "filter-period监测站数错误。"
  )
}

cat(
  "filter-period = ",
  format(
    nrow(FP_BASE),
    big.mark = ","
  ),
  " rows / ",
  uniqueN(
    FP_BASE$loc_num
  ),
  " monitors\n",
  sep = ""
)


# ============================================================
# 9. Expand filter-period to individual dates
# ============================================================

EXPANDED <- FP_BASE[
  ,
  .(
    date =
      seq(
        collect_start,
        collect_end -
          1L,
        by =
          "1 day"
      )
  ),
  by =
    .(
      period_id,
      loc_num,
      n_period_days
    )
]

EXPANDED[
  ,
  date :=
    as.IDate(
      date
    )
]

EXPANDED[
  ,
  date_year :=
    safe_int(
      format(
        date,
        "%Y"
      )
    )
]

if (
  EXPANDED[
    ,
    .N,
    by =
      .(
        period_id,
        n_period_days
      )
  ][
    N !=
      n_period_days,
    .N
  ] >
    0L
) {

  stopf(
    "filter-period日期展开数量错误。"
  )
}

CROSS_YEAR_PERIODS <- FP_BASE[
  safe_int(
    format(
      collect_start,
      "%Y"
    )
  ) !=
    safe_int(
      format(
        collect_end -
          1L,
        "%Y"
      )
    ),
  .N
]

cat(
  "跨自然年filter periods = ",
  CROSS_YEAR_PERIODS,
  "（H14按每个实际date所属年份匹配geometry）\n",
  sep = ""
)


# ============================================================
# 10. Build 360° actual-daily exposures
# ============================================================

Q_DOMAIN <- unique(
  Q_BASE[
    ,
    .(
      loc_num,
      date
    )
  ]
)

Q_DOMAIN[
  ,
  date_year :=
    safe_int(
      format(
        date,
        "%Y"
      )
    )
]

FP_DOMAIN <- unique(
  EXPANDED[
    ,
    .(
      loc_num,
      date,
      date_year
    )
  ]
)

ALL_DOMAIN <- unique(
  rbindlist(
    list(
      Q_DOMAIN,
      FP_DOMAIN
    ),
    use.names = TRUE,
    fill = TRUE
  )
)

ND_ALL_COLS <- unlist(
  lapply(
    FUELS,
    function(g) {

      unlist(
        lapply(
          RADII,
          function(r) {

            c(
              nd_col(
                g,
                "count",
                r
              ),

              nd_col(
                g,
                "heat",
                r
              )
            )
          }
        )
      )
    }
  )
)

build_nd_year <- function(
  yy
) {

  out_file <- file.path(
    ND_DAILY_DIR,
    paste0(
      "07H14_nondirectional_fuel_daily_exposure_",
      yy,
      ".csv.gz"
    )
  )

  expected_cols <- c(
    "loc_num",
    "date",
    ND_ALL_COLS
  )

  if (
    file.exists(
      out_file
    ) &&
    !FORCE_REBUILD_NONDIR
  ) {

    H <- names(
      fread(
        out_file,
        nrows = 0,
        showProgress = FALSE
      )
    )

    if (
      all(
        expected_cols %in%
          H
      )
    ) {

      cat(
        "[360° cache] ",
        yy,
        "：读取现有文件\n",
        sep = ""
      )

      return(
        data.table(
          year = yy,
          status = "REUSED",
          file = out_file
        )
      )
    }
  }

  cat(
    "[360° build] ",
    yy,
    "\n",
    sep = ""
  )

  DOMAIN <- ALL_DOMAIN[
    date_year ==
      yy,
    .(
      loc_num,
      date
    )
  ]

  if (
    nrow(DOMAIN) ==
      0L
  ) {

    stopf(
      yy,
      "年没有分析domain。"
    )
  }

  GY <- GEO[
    year ==
      yy &
    loc_num %in%
      unique(
        DOMAIN$loc_num
      )
  ]

  AY <- read_activity_year(
    yy
  )

  if (
    nrow(GY) >
      0L &&
    nrow(AY) >
      0L
  ) {

    X <- merge(
      GY,
      AY,
      by =
        "plant_id",
      all =
        FALSE,
      allow.cartesian =
        TRUE,
      sort =
        FALSE
    )

    X <- X[
      DOMAIN,
      on =
        .(
          loc_num,
          date
        ),
      nomatch =
        0
    ]

  } else {

    X <- data.table()
  }

  OUTY <- copy(
    DOMAIN
  )

  for (
    g in
    FUELS
  ) {

    XG <- if (
      nrow(X) >
        0L
    ) {

      X[
        fuel_group ==
          g
      ]

    } else {

      data.table()
    }

    for (
      r in
      RADII
    ) {

      cn <- nd_col(
        g,
        "count",
        r
      )

      hn <- nd_col(
        g,
        "heat",
        r
      )

      OUTY[
        ,
        (cn) :=
          0
      ]

      OUTY[
        ,
        (hn) :=
          0
      ]

      if (
        nrow(XG) >
          0L
      ) {

        XR <- XG[
          distance_km <=
            r +
            1e-8
        ]

        if (
          nrow(XR) >
            0L
        ) {

          CC <- XR[
            ,
            .(
              value =
                uniqueN(
                  plant_id
                )
            ),
            by =
              .(
                loc_num,
                date
              )
          ]

          HH <- XR[
            ,
            .(
              value =
                sum_keep_all_na(
                  heat
                )
            ),
            by =
              .(
                loc_num,
                date
              )
          ]

          OUTY[
            CC,
            on =
              .(
                loc_num,
                date
              ),
            (cn) :=
              i.value
          ]

          OUTY[
            HH,
            on =
              .(
                loc_num,
                date
              ),
            (hn) :=
              i.value
          ]
        }
      }
    }
  }

  fwrite(
    OUTY,
    out_file,
    compress =
      "gzip"
  )

  data.table(
    year =
      yy,

    status =
      "BUILT",

    n_domain =
      nrow(DOMAIN),

    n_geometry =
      nrow(GY),

    n_activity =
      nrow(AY),

    n_linked =
      nrow(X),

    file =
      out_file
  )
}

BUILD_AUDIT <- rbindlist(
  lapply(
    YEARS,
    build_nd_year
  ),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  BUILD_AUDIT,
  file.path(
    AUDIT_DIR,
    "01_nondirectional_daily_build_audit.csv"
  )
)

ND_FILES <- file.path(
  ND_DAILY_DIR,
  paste0(
    "07H14_nondirectional_fuel_daily_exposure_",
    YEARS,
    ".csv.gz"
  )
)

if (
  sum(
    file.exists(
      ND_FILES
    )
  ) !=
    24L
) {

  stopf(
    "H14 non-directional daily files不是24/24。"
  )
}


# ============================================================
# 11. Model functions
# ============================================================

BASE_VARS <- c(
  "gross_beta_mBq_m3",
  "c_prop",
  "surface_u_ppm_article",
  "sunspot_monthly",
  "pm_centered",
  "inv_hpbl",
  "rhum",
  "soilm",
  "vel",
  "lat_diff_35",
  "year_c",
  "temp_c",
  "loc_num",
  "monitor_year_fe"
)

make_common_data <- function(
  D,
  nd_name,
  up_name,
  scale_value
) {

  req <- c(
    BASE_VARS,
    nd_name,
    up_name
  )

  miss <- setdiff(
    req,
    names(D)
  )

  if (
    length(miss) >
      0L
  ) {

    stopf(
      "模型数据缺字段：",
      paste(
        miss,
        collapse = ", "
      )
    )
  }

  d <- D[
    ,
    ..req
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
        scale_value,

      exposure_up =
        safe_num(
          exposure_up_raw
        ) /
        scale_value
    )
  ]

  needed <- c(
    BASE_VARS,
    "exposure_nd",
    "exposure_up"
  )

  d <- d[
    complete.cases(
      d[
        ,
        ..needed
      ]
    )
  ]

  d
}

within_sd_group <- function(
  d,
  col,
  group_col
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
      .(
        grp =
          get(
            group_col
          )
      )
  ]$V1

  sd(
    z,
    na.rm = TRUE
  )
}


# ------------------------------------------------------------
# Fit one model
# ------------------------------------------------------------

fit_model <- function(
  d,
  exposure_name,
  model_stage
) {

  z <- copy(
    d
  )

  z[
    ,
    exposure :=
      get(
        exposure_name
      )
  ]

  if (
    nrow(z) <
      1000L ||
    uniqueN(
      z$loc_num
    ) <
      2L ||
    sd(
      z$exposure
    ) <=
      0
  ) {

    return(
      list(
        status =
          "SKIP",

        note =
          "Insufficient rows/monitors/exposure variation",

        fit =
          NULL
      )
    )
  }

  warning_text <- character()

  # ----------------------------------------------------------
  # M1 / M2
  # Paper-like random-intercept LME
  # ----------------------------------------------------------

  if (
    model_stage %in%
      c(
        "M1",
        "M2"
      )
  ) {

    FORM <-
      gross_beta_mBq_m3 ~
        exposure +
        exposure:vel +
        c_prop +
        surface_u_ppm_article +
        sunspot_monthly +
        pm_centered +
        inv_hpbl +
        rhum +
        soilm +
        vel +
        lat_diff_35 +
        year_c +
        I(year_c^2) +
        temp_c +
        I(temp_c^2) +
        I(temp_c^3) +
        sunspot_monthly:lat_diff_35 +
        surface_u_ppm_article:soilm +
        (1 | loc_num)

    fit <- tryCatch(

      withCallingHandlers(

        lme4::lmer(
          FORM,
          data =
            z,
          REML =
            TRUE,
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

  # ----------------------------------------------------------
  # M3
  # monitor fixed effect
  #
  # surface U and latitude main effects are absorbed by monitor FE.
  # Their interactions with time-varying covariates remain.
  # ----------------------------------------------------------

  } else if (
    model_stage ==
      "M3"
  ) {

    FORM <-
      gross_beta_mBq_m3 ~
        exposure +
        exposure:vel +
        c_prop +
        sunspot_monthly +
        pm_centered +
        inv_hpbl +
        rhum +
        soilm +
        vel +
        year_c +
        I(year_c^2) +
        temp_c +
        I(temp_c^2) +
        I(temp_c^3) +
        sunspot_monthly:lat_diff_35 +
        surface_u_ppm_article:soilm |
        loc_num

    fit <- tryCatch(

      withCallingHandlers(

        fixest::feols(
          FORM,
          data =
            z,
          vcov =
            ~loc_num,
          notes =
            FALSE,
          warn =
            FALSE
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

  # ----------------------------------------------------------
  # M4
  # monitor × year fixed effect
  #
  # year terms are absorbed by monitor×year FE.
  # ----------------------------------------------------------

  } else if (
    model_stage ==
      "M4"
  ) {

    FORM <-
      gross_beta_mBq_m3 ~
        exposure +
        exposure:vel +
        c_prop +
        sunspot_monthly +
        pm_centered +
        inv_hpbl +
        rhum +
        soilm +
        vel +
        temp_c +
        I(temp_c^2) +
        I(temp_c^3) +
        sunspot_monthly:lat_diff_35 +
        surface_u_ppm_article:soilm |
        monitor_year_fe

    fit <- tryCatch(

      withCallingHandlers(

        fixest::feols(
          FORM,
          data =
            z,
          vcov =
            ~loc_num,
          notes =
            FALSE,
          warn =
            FALSE
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

  } else {

    stopf(
      "未知model_stage：",
      model_stage
    )
  }

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

        fit =
          NULL
      )
    )
  }

  beta <- if (
    inherits(
      fit,
      "merMod"
    )
  ) {

    lme4::fixef(
      fit
    )

  } else {

    coef(
      fit
    )
  }

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
          "SKIP",

        note =
          paste(
            "Exposure main/interaction unavailable after rank handling",
            collapse_messages(
              warning_text
            )
          ),

        fit =
          NULL
      )
    )
  }

  list(
    status =
      "PASS",

    note =
      collapse_messages(
        warning_text
      ),

    fit =
      fit,

    beta =
      beta,

    V =
      V
  )
}


# ------------------------------------------------------------
# Model checkpoint
# ------------------------------------------------------------

checkpoint_path <- function(
  model_stage,
  fuel,
  metric,
  radius,
  geometry
) {

  file.path(
    CHECK_DIR,
    paste0(
      model_stage,
      "__",
      fuel,
      "__",
      metric,
      "__",
      radius,
      "km__",
      geometry,
      ".rds"
    )
  )
}

run_one_fit <- function(
  d,
  exposure_name,
  model_stage,
  fuel,
  metric,
  radius,
  geometry,
  branch
) {

  cp <- checkpoint_path(
    model_stage,
    fuel,
    metric,
    radius,
    geometry
  )

  if (
    file.exists(cp) &&
    !FORCE_REFIT
  ) {

    obj <- tryCatch(
      readRDS(cp),
      error =
        function(e)
          NULL
    )

    if (
      !is.null(obj) &&
      identical(
        obj$version,
        CHECKPOINT_VERSION
      )
    ) {

      cat(
        "  [checkpoint] ",
        model_stage,
        " | ",
        geometry,
        "\n",
        sep = ""
      )

      return(
        obj$result
      )
    }
  }

  ans <- fit_model(
    d,
    exposure_name,
    model_stage
  )

  median_wind <- median(
    d$vel,
    na.rm = TRUE
  )

  STATUS <- data.table(

    model_stage =
      model_stage,

    model_label =
      unname(MODEL_LABEL[model_stage]),

    branch =
      branch,

    geometry =
      geometry,

    fuel =
      fuel,

    metric =
      metric,

    radius_km =
      radius,

    status =
      ans$status,

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

    median_wind_ms =
      median_wind,

    note =
      ans$note
  )

  EFFECTS <- data.table()

  if (
    ans$status ==
      "PASS"
  ) {

    refs <- data.table(

      wind_reference =
        c(
          "paper_1m_s",
          "sample_median"
        ),

      wind_ms =
        c(
          1,
          median_wind
        )
    )

    EFFECTS <- rbindlist(
      lapply(
        seq_len(
          nrow(refs)
        ),
        function(i) {

          e <- effect_at_wind(
            ans$beta,
            ans$V,
            refs$wind_ms[i]
          )

          data.table(

            model_stage =
              model_stage,

            model_label =
              unname(MODEL_LABEL[model_stage]),

            branch =
              branch,

            geometry =
              geometry,

            fuel =
              fuel,

            metric =
              metric,

            radius_km =
              radius,

            wind_reference =
              refs$wind_reference[i],

            wind_ms =
              refs$wind_ms[i],

            estimate =
              unname(
                e[
                  "estimate"
                ]
              ),

            se =
              unname(
                e[
                  "se"
                ]
              ),

            ci95_low =
              unname(
                e[
                  "ci_low"
                ]
              ),

            ci95_high =
              unname(
                e[
                  "ci_high"
                ]
              ),

            p_wald =
              unname(
                e[
                  "p"
                ]
              ),

            direction =
              fifelse(
                e[
                  "estimate"
                ] >
                  0,
                "positive",
                fifelse(
                  e[
                    "estimate"
                  ] <
                    0,
                  "negative",
                  "zero"
                )
              ),

            ci_status =
              fifelse(
                e[
                  "ci_low"
                ] >
                  0,
                "CI_above_0",
                fifelse(
                  e[
                    "ci_high"
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
      )
    )
  }

  result <- list(
    status =
      STATUS,

    effects =
      EFFECTS
  )

  saveRDS(
    list(
      version =
        CHECKPOINT_VERSION,

      result =
        result
    ),
    cp
  )

  result
}


# ============================================================
# 12. Exposure diagnostics + model ladder
# ============================================================

STATUS_LIST <- list()
EFFECT_LIST <- list()
EXPOSURE_DIAG_LIST <- list()

kk_status <- 0L
kk_effect <- 0L
kk_diag <- 0L


for (
  fuel in
  FUELS
) {

  for (
    metric in
    METRICS
  ) {

    cat(
      "\n============================================================\n",
      "Fuel/metric: ",
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

    daily_parts <- vector(
      "list",
      length(YEARS)
    )

    for (
      ii in
      seq_along(YEARS)
    ) {

      yy <- YEARS[ii]

      ndf <- file.path(
        ND_DAILY_DIR,
        paste0(
          "07H14_nondirectional_fuel_daily_exposure_",
          yy,
          ".csv.gz"
        )
      )

      upf <- file.path(
        RAW_DAILY_DIR,
        paste0(
          "07H_fuel_daily_exposure_",
          yy,
          ".csv.gz"
        )
      )

      HN <- names(
        fread(
          ndf,
          nrows = 0,
          showProgress = FALSE
        )
      )

      HU <- names(
        fread(
          upf,
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
        length(miss_nd) >
          0L
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
        length(miss_up) >
          0L
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
        ndf,
        select = c(
          "loc_num",
          "date",
          nd_names
        ),
        showProgress = FALSE
      )

      B <- fread(
        upf,
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

      if (
        A[
          ,
          .N,
          by =
            .(
              loc_num,
              date
            )
        ][
          N >
            1L,
          .N
        ] >
          0L
      ) {

        stopf(
          yy,
          "年360° daily key重复。"
        )
      }

      if (
        B[
          ,
          .N,
          by =
            .(
              loc_num,
              date
            )
        ][
          N >
            1L,
          .N
        ] >
          0L
      ) {

        stopf(
          yy,
          "年upwind daily key重复。"
        )
      }

      daily_parts[[ii]] <- merge(
        A,
        B,
        by =
          c(
            "loc_num",
            "date"
          ),
        all.x =
          TRUE,
        sort =
          FALSE
      )
    }

    E_DAILY <- rbindlist(
      daily_parts,
      use.names = TRUE,
      fill = TRUE
    )

    rm(
      daily_parts
    )

    invisible(
      gc()
    )

    if (
      E_DAILY[
        ,
        .N,
        by =
          .(
            loc_num,
            date
          )
      ][
        N >
          1L,
        .N
      ] >
        0L
    ) {

      stopf(
        "合并后的daily exposure key重复：",
        fuel,
        " / ",
        metric
      )
    }


    # --------------------------------------------------------
    # quasi-daily
    # --------------------------------------------------------

    QX <- merge(
      Q_BASE,
      E_DAILY,
      by =
        c(
          "loc_num",
          "date"
        ),
      all.x =
        TRUE,
      sort =
        FALSE
    )

    if (
      nrow(QX) !=
        nrow(Q_BASE)
    ) {

      stopf(
        "quasi-daily暴露合并改变行数：",
        fuel,
        "/",
        metric
      )
    }


    # --------------------------------------------------------
    # filter-period strict-complete daily mean
    # --------------------------------------------------------

    exposure_cols <- c(
      nd_names,
      up_names
    )

    J <- E_DAILY[
      EXPANDED,
      on =
        .(
          loc_num,
          date
        )
    ]

    FPX <- J[
      ,
      lapply(
        .SD,
        function(x) {

          expected_n <- first(
            n_period_days
          )

          strict_mean(
            x,
            expected_n
          )
        }
      ),
      by =
        period_id,
      .SDcols =
        exposure_cols
    ]

    FX <- merge(
      FP_BASE,
      FPX,
      by =
        "period_id",
      all.x =
        TRUE,
      sort =
        FALSE
    )

    if (
      nrow(FX) !=
        nrow(FP_BASE)
    ) {

      stopf(
        "filter-period暴露合并改变行数：",
        fuel,
        "/",
        metric
      )
    }


    # ========================================================
    # Each radius
    # ========================================================

    for (
      r in
      RADII
    ) {

      ncn <- nd_col(
        fuel,
        metric,
        r
      )

      ucn <- up_col(
        fuel,
        metric,
        r
      )

      sc <- metric_scale(
        metric
      )

      # Same complete-case sample for 360° and upwind.
      dQ <- make_common_data(
        QX,
        ncn,
        ucn,
        sc
      )

      dF <- make_common_data(
        FX,
        ncn,
        ucn,
        sc
      )


      # ------------------------------------------------------
      # Exposure diagnostics
      # ------------------------------------------------------

      kk_diag <- kk_diag +
        1L

      EXPOSURE_DIAG_LIST[[kk_diag]] <- rbindlist(
        list(

          data.table(

            branch =
              "quasi_daily",

            fuel =
              fuel,

            metric =
              metric,

            radius_km =
              r,

            n_rows =
              nrow(dQ),

            n_monitors =
              uniqueN(
                dQ$loc_num
              ),

            nd_mean =
              mean(
                dQ$exposure_nd
              ),

            nd_sd =
              sd(
                dQ$exposure_nd
              ),

            nd_zero_pct =
              mean(
                dQ$exposure_nd ==
                  0
              ) *
              100,

            up_mean =
              mean(
                dQ$exposure_up
              ),

            up_sd =
              sd(
                dQ$exposure_up
              ),

            up_zero_pct =
              mean(
                dQ$exposure_up ==
                  0
              ) *
              100,

            nd_up_cor =
              suppressWarnings(
                cor(
                  dQ$exposure_nd,
                  dQ$exposure_up
                )
              ),

            nd_within_monitor_year_sd =
              within_sd_group(
                dQ,
                "exposure_nd",
                "monitor_year_fe"
              ),

            up_within_monitor_year_sd =
              within_sd_group(
                dQ,
                "exposure_up",
                "monitor_year_fe"
              )
          ),

          data.table(

            branch =
              "filter_period",

            fuel =
              fuel,

            metric =
              metric,

            radius_km =
              r,

            n_rows =
              nrow(dF),

            n_monitors =
              uniqueN(
                dF$loc_num
              ),

            nd_mean =
              mean(
                dF$exposure_nd
              ),

            nd_sd =
              sd(
                dF$exposure_nd
              ),

            nd_zero_pct =
              mean(
                dF$exposure_nd ==
                  0
              ) *
              100,

            up_mean =
              mean(
                dF$exposure_up
              ),

            up_sd =
              sd(
                dF$exposure_up
              ),

            up_zero_pct =
              mean(
                dF$exposure_up ==
                  0
              ) *
              100,

            nd_up_cor =
              suppressWarnings(
                cor(
                  dF$exposure_nd,
                  dF$exposure_up
                )
              ),

            nd_within_monitor_year_sd =
              within_sd_group(
                dF,
                "exposure_nd",
                "monitor_year_fe"
              ),

            up_within_monitor_year_sd =
              within_sd_group(
                dF,
                "exposure_up",
                "monitor_year_fe"
              )
          )
        ),
        use.names = TRUE,
        fill = TRUE
      )


      # ------------------------------------------------------
      # Model ladder
      # ------------------------------------------------------

      specs <- list(

        list(
          stage =
            "M1",

          branch =
            "quasi_daily",

          data =
            dQ
        ),

        list(
          stage =
            "M2",

          branch =
            "filter_period",

          data =
            dF
        ),

        list(
          stage =
            "M3",

          branch =
            "filter_period",

          data =
            dF
        ),

        list(
          stage =
            "M4",

          branch =
            "filter_period",

          data =
            dF
        )
      )

      cat(
        "\n[",
        fuel,
        " | ",
        metric,
        " | ",
        r,
        " km] Q rows=",
        nrow(dQ),
        "; FP rows=",
        nrow(dF),
        "\n",
        sep = ""
      )

      for (
        sp in
        specs
      ) {

        for (
          geom in
          c(
            "non_directional_360",
            "upwind_90"
          )
        ) {

          exp_name <- if (
            geom ==
              "non_directional_360"
          ) {

            "exposure_nd"

          } else {

            "exposure_up"
          }

          rr <- run_one_fit(

            sp$data,

            exp_name,

            sp$stage,

            fuel,

            metric,

            r,

            geom,

            sp$branch
          )

          kk_status <- kk_status +
            1L

          STATUS_LIST[[kk_status]] <- rr$status

          if (
            nrow(
              rr$effects
            ) >
              0L
          ) {

            kk_effect <- kk_effect +
              1L

            EFFECT_LIST[[kk_effect]] <- rr$effects
          }
        }
      }

      rm(
        dQ,
        dF
      )

      invisible(
        gc()
      )
    }

    rm(
      E_DAILY,
      QX,
      J,
      FPX,
      FX
    )

    invisible(
      gc()
    )
  }
}


# ============================================================
# 13. Save full results
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

EXPOSURE_DIAG <- rbindlist(
  EXPOSURE_DIAG_LIST,
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  EXPOSURE_DIAG,
  file.path(
    OUT,
    "01_exposure_diagnostics.csv"
  )
)

fwrite(
  STATUS,
  file.path(
    OUT,
    "02_model_status.csv"
  )
)

fwrite(
  EFFECTS,
  file.path(
    OUT,
    "03_all_effects_1ms_and_median.csv"
  )
)


# ============================================================
# 14. Summaries
# ============================================================

SUMMARY <- EFFECTS[
  ,
  .(
    n_effects =
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
      )
  ),
  by =
    .(
      model_stage,
      model_label,
      branch,
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
    "04_direction_summary_across_20_50km.csv"
  )
)

PAIR <- dcast(
  EFFECTS,
  model_stage +
    model_label +
    branch +
    fuel +
    metric +
    radius_km +
    wind_reference +
    wind_ms ~
    geometry,
  value.var =
    c(
      "estimate",
      "ci95_low",
      "ci95_high"
    )
)

fwrite(
  PAIR,
  file.path(
    OUT,
    "05_360_vs_upwind_side_by_side.csv"
  )
)

ONE_MS <- EFFECTS[
  wind_reference ==
    "paper_1m_s"
]

MED <- EFFECTS[
  wind_reference ==
    "sample_median"
]

fwrite(
  ONE_MS,
  file.path(
    OUT,
    "06_effects_at_1m_s.csv"
  )
)

fwrite(
  MED,
  file.path(
    OUT,
    "07_effects_at_sample_median_wind.csv"
  )
)

STAGE_WIDE <- dcast(
  EFFECTS,
  geometry +
    fuel +
    metric +
    radius_km +
    wind_reference +
    wind_ms ~
    model_stage,
  value.var =
    "estimate"
)

fwrite(
  STAGE_WIDE,
  file.path(
    OUT,
    "08_model_ladder_estimate_wide.csv"
  )
)


# ------------------------------------------------------------
# Print most important first diagnostic: 1 m/s
# ------------------------------------------------------------

cat(
  "\n============================================================\n",
  "1 m/s方向摘要（20–50 km）\n",
  "============================================================\n",
  sep = ""
)

print(
  SUMMARY[
    wind_reference ==
      "paper_1m_s",
    .(
      model_stage,
      geometry,
      fuel,
      metric,
      n_effects,
      n_positive,
      n_negative,
      n_ci_above_0,
      n_ci_below_0
    )
  ][
    order(
      model_stage,
      geometry,
      fuel,
      metric
    )
  ],
  nrows = Inf
)


# ============================================================
# 15. Diagnostic plots
# ============================================================

plot_diag <- function(
  metric,
  wind_ref,
  filename,
  title_text
) {

  metric_value <- metric
  wind_value <- wind_ref

  D <- EFFECTS[
    metric ==
      metric_value &
    wind_reference ==
      wind_value
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
    model_label_plot :=
      factor(
        model_label,
        levels =
          unname(
            MODEL_LABEL
          )
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
        levels =
          c(
            "360° non-directional",
            "90° upwind"
          )
      )
  ]

  pd <- position_dodge(
    width =
      1.8
  )

  p <- ggplot(
    D,
    aes(
      x =
        radius_km,

      y =
        estimate,

      group =
        model_label_plot,

      shape =
        model_label_plot
    )
  ) +

    geom_hline(
      yintercept =
        0,
      linetype =
        "dashed",
      linewidth =
        0.45
    ) +

    geom_errorbar(
      aes(
        ymin =
          ci95_low,

        ymax =
          ci95_high
      ),
      width =
        1.1,
      linewidth =
        0.45,
      position =
        pd
    ) +

    geom_line(
      aes(
        linetype =
          model_label_plot
      ),
      linewidth =
        0.55,
      position =
        pd
    ) +

    geom_point(
      size =
        1.7,
      position =
        pd
    ) +

    facet_grid(
      geometry_label ~
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
        title_text,

      subtitle =
        paste0(
          "Same 360°/upwind complete-case sample within each fuel × metric × radius; ",
          "error bars are 95% CIs."
        ),

      x =
        "Cumulative radius (km)",

      y =
        if (
          metric ==
            "heat"
        ) {

          "Adjusted gross-beta change per 10,000 MMBtu/day (mBq/m³)"

        } else {

          "Adjusted gross-beta change per additional active plant (mBq/m³)"
        },

      shape =
        "Model ladder",

      linetype =
        "Model ladder"
    ) +

    theme_bw(
      base_size =
        9
    ) +

    theme(
      panel.grid.minor =
        element_blank(),

      strip.background =
        element_rect(
          fill =
            "white"
        ),

      strip.text =
        element_text(
          face =
            "bold"
        ),

      plot.title =
        element_text(
          face =
            "bold",
          size =
            12
        ),

      legend.position =
        "bottom",

      axis.text.x =
        element_text(
          angle =
            45,
          hjust =
            1
        )
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
    width =
      15,
    height =
      8.5,
    units =
      "in",
    dpi =
      300,
    bg =
      "white"
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
      width =
        15,
      height =
        8.5,
      units =
        "in",
      device =
        cairo_pdf
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
      width =
        15,
      height =
        8.5,
      units =
        "in"
    )
  }
}

plot_diag(
  "count",
  "paper_1m_s",
  "01_Count_1ms_model_ladder",
  "Count: model-ladder diagnosis at 1 m/s"
)

plot_diag(
  "heat",
  "paper_1m_s",
  "02_Heat_1ms_model_ladder",
  "Heat Input: model-ladder diagnosis at 1 m/s"
)

plot_diag(
  "count",
  "sample_median",
  "03_Count_median_wind_model_ladder",
  "Count: model-ladder diagnosis at sample median wind"
)

plot_diag(
  "heat",
  "sample_median",
  "04_Heat_median_wind_model_ladder",
  "Heat Input: model-ladder diagnosis at sample median wind"
)


# ============================================================
# 16. Hard checks
# ============================================================

EXPECTED_STATUS <-
  length(FUELS) *
  length(METRICS) *
  length(RADII) *
  2L *
  4L

N_FAIL_MODELS <- STATUS[
  status ==
    "FAIL",
  .N
]

N_SKIP_MODELS <- STATUS[
  status ==
    "SKIP",
  .N
]

EXPECTED_EFFECT_ROWS <-
  STATUS[
    status ==
      "PASS",
    .N
  ] *
  2L

CHECKS <- data.table(

  check =
    c(
      "formal_quasi_rows",
      "formal_quasi_monitors",
      "formal_filter_rows",
      "formal_filter_monitors",
      "upwind_daily_files",
      "nondirectional_daily_files",
      "model_status_rows",
      "model_failures",
      "effect_rows_match_PASS_x2",
      "four_diagnostic_png",
      "four_diagnostic_pdf"
    ),

  observed =
    c(
      as.character(
        nrow(
          Q_BASE
        )
      ),

      as.character(
        uniqueN(
          Q_BASE$loc_num
        )
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
        sum(
          file.exists(
            UP_FILES
          )
        )
      ),

      as.character(
        sum(
          file.exists(
            ND_FILES
          )
        )
      ),

      as.character(
        nrow(
          STATUS
        )
      ),

      as.character(
        N_FAIL_MODELS
      ),

      as.character(
        nrow(
          EFFECTS
        )
      ),

      as.character(
        sum(
          file.exists(
            file.path(
              FIG_DIR,
              c(
                "01_Count_1ms_model_ladder.png",
                "02_Heat_1ms_model_ladder.png",
                "03_Count_median_wind_model_ladder.png",
                "04_Heat_median_wind_model_ladder.png"
              )
            )
          )
        )
      ),

      as.character(
        sum(
          file.exists(
            file.path(
              FIG_DIR,
              c(
                "01_Count_1ms_model_ladder.pdf",
                "02_Heat_1ms_model_ladder.pdf",
                "03_Count_median_wind_model_ladder.pdf",
                "04_Heat_median_wind_model_ladder.pdf"
              )
            )
          )
        )
      )
    ),

  expected =
    c(
      as.character(
        EXPECTED_Q_ROWS
      ),

      as.character(
        EXPECTED_Q_MONITORS
      ),

      as.character(
        EXPECTED_FP_ROWS
      ),

      as.character(
        EXPECTED_FP_MONITORS
      ),

      "24",
      "24",

      as.character(
        EXPECTED_STATUS
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
    "02_07H14_hard_checks.csv"
  )
)


# ============================================================
# 17. Final terminal summary
# ============================================================

cat(
  "\n============================================================\n",
  "07H14 MODEL STATUS\n",
  "============================================================\n",
  sep = ""
)

print(
  STATUS[
    ,
    .N,
    by =
      .(
        model_stage,
        status
      )
  ][
    order(
      model_stage,
      status
    )
  ],
  nrows = Inf
)

cat(
  "\nSKIP models = ",
  N_SKIP_MODELS,
  "\n",
  sep = ""
)

cat(
  "\n============================================================\n",
  "07H14 HARD CHECKS\n",
  "============================================================\n",
  sep = ""
)

print(
  CHECKS,
  nrows = Inf
)

N_HARD_FAIL <- CHECKS[
  status ==
    "FAIL",
  .N
]

cat(
  "\n============================================================\n"
)

if (
  N_HARD_FAIL ==
    0L
) {

  cat(
    "07H14 MODEL-LADDER DIAGNOSIS: PASS\n"
  )

} else {

  cat(
    "07H14 MODEL-LADDER DIAGNOSIS: FAIL\n"
  )
}

cat(
  "Hard FAIL = ",
  N_HARD_FAIL,
  "\n",
  sep = ""
)

cat(
  "Model FAIL = ",
  N_FAIL_MODELS,
  "\n",
  sep = ""
)

cat(
  "Model SKIP = ",
  N_SKIP_MODELS,
  "\n",
  sep = ""
)

cat(
  "Effects = ",
  nrow(
    EFFECTS
  ),
  "\n",
  sep = ""
)

cat(
  "Output = ",
  OUT,
  "\n",
  sep = ""
)

cat(
  "============================================================\n"
)


# ============================================================
# 18. Report
# ============================================================

REPORT <- c(

  "# 07H14 模型阶梯诊断",

  "",

  "- 目的：定位不同燃料实际运行暴露的关联方向在哪个分析层级发生变化。",

  "- 核心燃料：Biomass / Coal / Natural gas / Oil。",

  "- 半径：20–50 km，每5 km。",

  "- 指标：daily active plant Count 与 actual daily Heat Input。",

  "- 几何：360° non-directional 与 90° upwind；每个fuel×metric×radius内使用共同complete-case样本。",

  "- M1：quasi-daily + monitor random intercept，最接近参考O&G论文统计结构。",

  "- M2：filter-period + monitor random intercept。",

  "- M3：filter-period + monitor fixed effect。",

  "- M4：filter-period + monitor×year fixed effect，接近冻结07H10主模型。",

  "- 四层均包含Exposure×wind interaction，并报告1 m/s与样本中位风速效应。",

  "- Heat效应单位：每10,000 MMBtu/day。Count效应单位：每增加1座活动电厂。",

  "- 360°模型中的Exposure×wind仅用于matched diagnostic，使其与upwind模型结构一致；不表示这是最终非方向主模型的唯一选择。",

  paste0(
    "- Model FAIL = ",
    N_FAIL_MODELS,
    "; SKIP = ",
    N_SKIP_MODELS,
    "."
  ),

  "",

  if (
    N_HARD_FAIL ==
      0L
  ) {

    "**07H14 MODEL-LADDER DIAGNOSIS: PASS**"

  } else {

    "**07H14 MODEL-LADDER DIAGNOSIS: FAIL**"
  }
)

writeLines(
  REPORT,
  file.path(
    OUT,
    "07H14_模型阶梯诊断报告.md"
  ),
  useBytes = TRUE
)

if (
  N_HARD_FAIL >
    0L
) {

  stopf(
    "07H14存在Hard FAIL，请先检查审计文件。"
  )
}
