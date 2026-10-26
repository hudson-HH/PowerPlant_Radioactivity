#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999
)


# ============================================================
# 07H02
# Fuel-specific CAMPD Count + Heat cumulative-radius models
#
# 正式框架
# ------------------------------------------------------------
# Fuel:
#   coal
#   gas
#   oil
#   other_fossil
#   biomass
#   solid_waste
#
# Radius:
#   20 / 25 / 30 / 35 / 40 / 45 / 50 / 100 / 150 / 200 km
#
# Direction:
#   upwind   = main directional analysis
#   downwind = directional negative control
#
# For every fuel × radius × direction × branch:
#
#   A. Count-only
#   B. Heat-only
#   C. Joint Count + Heat
#
# Branch:
#   1. quasi_daily
#   2. filter_period
#
# Effect reporting:
#   Count = +1 plant
#   Heat  = +10,000 MMBtu
#   1 SD  = supplementary standardized increment
#
# Joint interpretation:
#   Count conditional on Heat
#   Heat conditional on Count
#
# 不修改任何：
#   07A / 07B2 / 07B3 / 07C / 07D / 07E / 07F / 07G / 07H01
#
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

required_packages <- c(
  "data.table",
  "lme4"
)

for (
  pkg in required_packages
) {

  if (
    !requireNamespace(
      pkg,
      quietly = TRUE
    )
  ) {

    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}


suppressPackageStartupMessages({

  library(
    data.table
  )

  library(
    lme4
  )
})


# ============================================================
# 1. Paths
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

F07_FREEZE <- file.path(
  F07_ROOT,
  "07F_FINAL_FREEZE.md"
)

F07_QDIR <- file.path(
  F07_ROOT,
  "quasi_daily_yearly"
)

H_ROOT <- file.path(
  PP,
  "06_07H_燃料特异CAMPD逐日暴露与模型"
)

H01_DAILY_DIR <- file.path(
  H_ROOT,
  "fuel_daily_yearly"
)

H01_QC_DIR <- file.path(
  H_ROOT,
  "qc"
)

H01_HARD <- file.path(
  H01_QC_DIR,
  "04_build_hard_checks.csv"
)

OUT <- file.path(
  H_ROOT,
  "07H02_model"
)

CHECKPOINT_DIR <- file.path(
  OUT,
  "checkpoints"
)

FILTER_EXP_DIR <- file.path(
  OUT,
  "filter_period_exposure"
)

QC_DIR <- file.path(
  OUT,
  "qc"
)

for (
  d in c(
    OUT,
    CHECKPOINT_DIR,
    FILTER_EXP_DIR,
    QC_DIR
  )
) {

  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


OUT_EFFECTS <- file.path(
  OUT,
  "01_07H02_model_effects_all.csv"
)

OUT_STATUS <- file.path(
  OUT,
  "02_07H02_model_status_all.csv"
)

OUT_COR <- file.path(
  OUT,
  "03_07H02_count_heat_collinearity.csv"
)

OUT_SAMPLE <- file.path(
  OUT,
  "04_07H02_branch_sample_summary.csv"
)

OUT_HARD <- file.path(
  OUT,
  "05_07H02_hard_checks.csv"
)

OUT_FILTER_LOCATOR <- file.path(
  QC_DIR,
  "01_07F_filter_base_locator_audit.csv"
)

OUT_PRIMARY_FLAG <- file.path(
  QC_DIR,
  "02_07F_primary_quasi_flag_audit.csv"
)

OUT_REPORT <- file.path(
  OUT,
  "07H02_燃料特异Count_Heat累计半径模型报告.md"
)

OUT_SESSION <- file.path(
  OUT,
  "sessionInfo.txt"
)


# ============================================================
# 2. Constants
# ============================================================

MODEL_VERSION <-
  "07H02_count_heat_joint_v2026-09-18"

FORCE_RESTART <- FALSE

YEARS <- 2001:2024

FUELS <- c(
  "coal",
  "gas",
  "oil",
  "other_fossil",
  "biomass",
  "solid_waste"
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

DIRECTIONS <- c(
  "upwind",
  "downwind"
)

BRANCHES <- c(
  "filter_period",
  "quasi_daily"
)

MODEL_TYPES <- c(
  "count_only",
  "heat_only",
  "joint_count_heat"
)

EFFECT_WIND_MS <- 1

HEAT_INCREMENT_MMBTU <- 10000

EXPECTED_07F_ROWS <- 769770L

EXPECTED_UNIQUE_MONITOR_DAYS <- 769600L

EXPECTED_FILTER_PERIODS <- 174200L

EXPECTED_QUASI_PRIMARY_ROWS <- 426657L

EXPECTED_FILTER_FORMAL_ROWS <- 107407L

EXPECTED_MODEL_MONITORS <- 133L

EXPECTED_MODEL_SPECS <-
  length(FUELS) *
  length(RADII) *
  length(DIRECTIONS) *
  length(BRANCHES) *
  length(MODEL_TYPES)

EXPECTED_EFFECT_ROWS <-
  length(FUELS) *
  length(RADII) *
  length(DIRECTIONS) *
  length(BRANCHES) *
  4L

EXPECTED_COR_ROWS <-
  length(FUELS) *
  length(RADII) *
  length(DIRECTIONS) *
  length(BRANCHES)


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


fmt_int <- function(x) {

  format(
    x,
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


as_flag <- function(
  x,
  missing_value = FALSE
) {

  if (
    is.logical(x)
  ) {

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


  if (
    !is.na(
      missing_value
    )
  ) {

    z[
      is.na(z)
    ] <- missing_value
  }


  z
}


collapse_messages <- function(x) {

  if (
    is.null(x) ||
    length(x) ==
      0L
  ) {

    return("")
  }


  x <- unique(
    trimws(
      as.character(x)
    )
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


safe_quantile <- function(
  x,
  p
) {

  x <- safe_num(x)

  x <- x[
    is.finite(x)
  ]


  if (
    length(x) ==
      0L
  ) {

    return(
      NA_real_
    )
  }


  as.numeric(
    quantile(
      x,
      probs = p,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
  )
}


choose_scale <- function(x) {

  x <- safe_num(x)

  x <- x[
    is.finite(x)
  ]


  if (
    length(x) <
      2L ||
    uniqueN(x) <
      2L
  ) {

    return(
      list(
        ok = FALSE,
        value = NA_real_,
        method = "no_variation"
      )
    )
  }


  s <- safe_quantile(
    x,
    0.95
  )


  if (
    is.finite(s) &&
    abs(s) >
      1e-12
  ) {

    return(
      list(
        ok = TRUE,
        value = abs(s),
        method = "P95_all"
      )
    )
  }


  nz <- abs(
    x[
      abs(x) >
        1e-12
    ]
  )


  if (
    length(nz) ==
      0L
  ) {

    return(
      list(
        ok = FALSE,
        value = NA_real_,
        method = "all_zero"
      )
    )
  }


  s2 <- safe_quantile(
    nz,
    0.95
  )


  if (
    !is.finite(s2) ||
    s2 <=
      0
  ) {

    return(
      list(
        ok = FALSE,
        value = NA_real_,
        method = "invalid_scale"
      )
    )
  }


  list(
    ok = TRUE,
    value = s2,
    method = "P95_nonzero_fallback"
  )
}


find_interaction_term <- function(
  beta_names,
  a,
  b
) {

  candidates <- c(
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
  )


  hit <- intersect(
    candidates,
    beta_names
  )


  if (
    length(hit) !=
      1L
  ) {

    return(
      NA_character_
    )
  }


  hit[1]
}


make_control <- function() {

  args <- list(

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


  if (
    "autoscale" %in%
      names(
        formals(
          lmerControl
        )
      )
  ) {

    args$autoscale <- TRUE
  }


  do.call(
    lmerControl,
    args
  )
}


safe_cor <- function(
  x,
  y,
  method
) {

  ok <-
    is.finite(x) &
    is.finite(y)

  x <- x[ok]
  y <- y[ok]


  if (
    length(x) <
      3L ||
    uniqueN(x) <
      2L ||
    uniqueN(y) <
      2L
  ) {

    return(
      NA_real_
    )
  }


  suppressWarnings(
    cor(
      x,
      y,
      method = method
    )
  )
}


term_info <- function(term) {

  if (
    term ==
      "count"
  ) {

    return(
      list(
        raw_increment =
          1,

        raw_unit =
          "plant",

        effect_label =
          "+1 plant"
      )
    )
  }


  list(
    raw_increment =
      HEAT_INCREMENT_MMBTU,

    raw_unit =
      "MMBtu",

    effect_label =
      "+10,000 MMBtu"
  )
}


interpretation_text <- function(
  model_type,
  term
) {

  if (
    model_type ==
      "count_only"
  ) {

    return(
      "Plant-count association"
    )
  }


  if (
    model_type ==
      "heat_only"
  ) {

    return(
      "Heat-Input association"
    )
  }


  if (
    term ==
      "count"
  ) {

    return(
      "Plant-count association conditional on Heat Input"
    )
  }


  "Heat-Input association conditional on plant count"
}


# ============================================================
# 4. Verify 07H01 and 07F frozen upstream
# ============================================================

if (
  !file.exists(
    H01_HARD
  )
) {

  stopf(
    "找不到07H01 hard-check文件：\n",
    H01_HARD
  )
}


H01_CHECK <- fread(
  H01_HARD,
  showProgress = FALSE
)


if (
  !"status" %in%
    names(
      H01_CHECK
    ) ||
  H01_CHECK[
    status !=
      "PASS",
    .N
  ] >
    0L
) {

  print(
    H01_CHECK
  )

  stopf(
    "07H01尚未全部hard PASS，禁止建模。"
  )
}


cat(
  "07H01 hard checks：PASS\n"
)


H01_FILES <- file.path(
  H01_DAILY_DIR,
  paste0(
    "07H_fuel_daily_exposure_",
    YEARS,
    ".csv.gz"
  )
)


if (
  !all(
    file.exists(
      H01_FILES
    )
  )
) {

  stopf(
    "07H01 daily yearly files不是24/24。"
  )
}


if (
  !file.exists(
    F07_FREEZE
  )
) {

  stopf(
    "找不到07F_FINAL_FREEZE.md。"
  )
}


freeze_text <- readLines(
  F07_FREEZE,
  warn = FALSE,
  encoding = "UTF-8"
)


if (
  !any(
    grepl(
      "07F FINAL FREEZE: PASS",
      freeze_text,
      fixed = TRUE
    )
  )
) {

  stopf(
    "07F尚未FINAL FREEZE PASS。"
  )
}


cat(
  "07F FINAL FREEZE：PASS\n"
)


# ============================================================
# 5. Read frozen 07F quasi universe
#    and automatically identify the frozen primary-quasi flag
# ============================================================

Q_FILES <- file.path(
  F07_QDIR,
  paste0(
    "07F_quasi_daily_analysis_base_",
    YEARS,
    ".csv.gz"
  )
)


if (
  !all(
    file.exists(
      Q_FILES
    )
  )
) {

  stopf(
    "07F quasi yearly files不是24/24。"
  )
}


Q_HEADERS <- lapply(
  Q_FILES,
  function(f) {

    names(
      fread(
        f,
        nrows = 0,
        showProgress = FALSE
      )
    )
  }
)


COMMON_Q_HEADER <- Reduce(
  intersect,
  Q_HEADERS
)


Q_CORE <- c(
  "loc_num",
  "date",
  "collect_start",
  "collect_end",
  "gross_beta_mBq_m3",
  "c_prop",
  "surface_u_ppm_article",
  "sunspot_monthly",
  "pm_same_day",
  "inv_hpbl",
  "rhum",
  "soilm",
  "vel",
  "lat_diff_35",
  "year_model",
  "temp"
)


missing_q_core <- setdiff(
  Q_CORE,
  COMMON_Q_HEADER
)


if (
  length(
    missing_q_core
  ) >
    0L
) {

  stopf(
    "07F quasi yearly文件缺少模型字段：",
    paste(
      missing_q_core,
      collapse = ", "
    )
  )
}


FLAG_CANDIDATES <- grep(
  "flag|sample|primary|formal",
  COMMON_Q_HEADER,
  ignore.case = TRUE,
  value = TRUE
)


FLAG_CANDIDATES <- setdiff(
  FLAG_CANDIDATES,
  Q_CORE
)


if (
  length(
    FLAG_CANDIDATES
  ) ==
    0L
) {

  stopf(
    "07F quasi中找不到任何候选sample/flag字段。"
  )
}


Q_SELECT <- unique(
  c(
    Q_CORE,
    FLAG_CANDIDATES
  )
)


cat(
  "\n读取07F 24年quasi数据...\n"
)


Q_ALL <- rbindlist(
  lapply(
    Q_FILES,
    function(f) {

      fread(
        f,
        select = Q_SELECT,
        showProgress = FALSE
      )
    }
  ),
  use.names = TRUE,
  fill = TRUE
)


Q_ALL[
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
      )
  )
]


if (
  nrow(
    Q_ALL
  ) !=
    EXPECTED_07F_ROWS
) {

  stopf(
    "07F quasi总行数不是冻结值",
    EXPECTED_07F_ROWS,
    "，实际=",
    nrow(Q_ALL)
  )
}


dup_q <- Q_ALL[
  ,
  .N,
  by = .(
    loc_num,
    collect_start,
    collect_end,
    date
  )
][
  N >
    1L,
  .N
]


if (
  dup_q >
    0L
) {

  stopf(
    "07F quasi observation key存在重复。"
  )
}


FLAG_AUDIT <- rbindlist(
  lapply(
    FLAG_CANDIDATES,
    function(cc) {

      z <- as_flag(
        Q_ALL[[cc]],
        missing_value = FALSE
      )


      data.table(
        flag =
          cc,

        n_true =
          sum(
            z %in%
              TRUE
          )
      )
    }
  )
)


setorder(
  FLAG_AUDIT,
  -n_true,
  flag
)


fwrite(
  FLAG_AUDIT,
  OUT_PRIMARY_FLAG
)


PRIMARY_HITS <- FLAG_AUDIT[
  n_true ==
    EXPECTED_QUASI_PRIMARY_ROWS,
  flag
]


if (
  length(
    PRIMARY_HITS
  ) ==
    0L
) {

  cat(
    "\n候选flag计数：\n"
  )

  print(
    FLAG_AUDIT
  )

  stopf(
    "没有任何07F frozen flag恰好对应426,657条primary quasi记录。\n",
    "候选计数已保存：",
    OUT_PRIMARY_FLAG
  )
}


if (
  length(
    PRIMARY_HITS
  ) >
    1L
) {

  ref_flag <- as_flag(
    Q_ALL[[
      PRIMARY_HITS[1]
    ]],
    FALSE
  )


  same_flag <- vapply(
    PRIMARY_HITS[-1],
    function(cc) {

      identical(
        ref_flag,
        as_flag(
          Q_ALL[[cc]],
          FALSE
        )
      )
    },
    logical(1)
  )


  if (
    !all(
      same_flag
    )
  ) {

    cat(
      "\n多个不同flag均得到426,657行：\n"
    )

    print(
      PRIMARY_HITS
    )

    stopf(
      "无法唯一确定07F primary quasi flag。"
    )
  }
}


PRIMARY_FLAG <- PRIMARY_HITS[1]


cat(
  "07F primary quasi flag = ",
  PRIMARY_FLAG,
  "\n",
  sep = ""
)


Q_ALL[
  ,
  .primary_flag :=
    as_flag(
      get(
        PRIMARY_FLAG
      ),
      FALSE
    )
]


Q_PRIMARY <- Q_ALL[
  .primary_flag %in%
    TRUE,
  ..Q_CORE
]

# quasi-daily正式PM定义为same-day PM。
# 模型内部统一使用pm_c这个工作变量名，
# 这里只在07H02内临时改名；不修改任何冻结07F文件。
setnames(
  Q_PRIMARY,
  "pm_same_day",
  "pm_c"
)


if (
  nrow(
    Q_PRIMARY
  ) !=
    EXPECTED_QUASI_PRIMARY_ROWS
) {

  stopf(
    "primary quasi行数错误。"
  )
}


if (
  uniqueN(
    Q_PRIMARY$loc_num
  ) !=
    EXPECTED_MODEL_MONITORS
) {

  stopf(
    "primary quasi monitor数不是133，实际=",
    uniqueN(
      Q_PRIMARY$loc_num
    )
  )
}


PERIOD_DAY_KEYS <- Q_ALL[
  ,
  .(
    loc_num,
    date,
    collect_start,
    collect_end
  )
]


N_PERIODS <- uniqueN(
  PERIOD_DAY_KEYS,
  by = c(
    "loc_num",
    "collect_start",
    "collect_end"
  )
)


if (
  N_PERIODS !=
    EXPECTED_FILTER_PERIODS
) {

  stopf(
    "正式filter-period universe不是",
    EXPECTED_FILTER_PERIODS,
    "，实际=",
    N_PERIODS
  )
}


rm(
  Q_ALL
)

invisible(
  gc()
)


# ============================================================

# ============================================================
# 6. Reproduce frozen 07G filter-period formal analysis base
# ============================================================
#
# IMPORTANT:
#
# There is NO separate pre-existing 107,407-row CSV.
#
# Frozen 07G reads:
#
#   03_07F_filter_period_analysis_base_2001_2024.csv.gz
#
# which contains all 174,200 filter periods, and then constructs
# its formal filter-period analysis base internally.
#
# 07H02 reproduces that exact 07G filtering logic here.
#
# Official filter-period PM field:
#
#   pm_filter_period_official
#
# NOT pm_c.
#
# For shared downstream model code only, we create:
#
#   pm_c := pm_filter_period_official
#
# inside the in-memory 07H02 FILTER_BASE.
#
# Frozen 07F is never modified.
# ============================================================


FILTER_BASE_FILE <- file.path(
  F07_ROOT,
  "03_07F_filter_period_analysis_base_2001_2024.csv.gz"
)


if (
  !file.exists(
    FILTER_BASE_FILE
  )
) {

  stopf(
    "找不到07F正式filter-period表：\n",
    FILTER_BASE_FILE
  )
}


# ------------------------------------------------------------
# Exact overall primary Heat variable used by frozen 07G when
# constructing its filter-period formal analysis base.
# ------------------------------------------------------------

FROZEN_PRIMARY_FP <- paste0(
  "fpmean_",
  "campd_upwind_heat_observed_20km_90deg"
)


FILTER_SOURCE_REQUIRED <- c(

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

  "flag_07f_filter_period_primary_base",

  FROZEN_PRIMARY_FP
)


FILTER_HEADER <- names(
  fread(
    FILTER_BASE_FILE,
    nrows = 0,
    showProgress = FALSE
  )
)


FILTER_MISSING <- setdiff(
  FILTER_SOURCE_REQUIRED,
  FILTER_HEADER
)


if (
  length(
    FILTER_MISSING
  ) >
    0L
) {

  stopf(
    "07F正式filter-period表缺少字段：",
    paste(
      FILTER_MISSING,
      collapse = ", "
    )
  )
}


FP07 <- fread(
  FILTER_BASE_FILE,
  select = FILTER_SOURCE_REQUIRED,
  showProgress = FALSE
)


# ------------------------------------------------------------
# Type normalization
# ------------------------------------------------------------

FP07[
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


FILTER_NUMERIC <- c(

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

  FROZEN_PRIMARY_FP
)


for (
  cc in FILTER_NUMERIC
) {

  FP07[
    ,
    (cc) :=
      safe_num(
        get(
          cc
        )
      )
  ]
}


# ------------------------------------------------------------
# Frozen dimensions
# ------------------------------------------------------------

if (
  nrow(
    FP07
  ) !=
    EXPECTED_FILTER_PERIODS
) {

  stopf(
    "07F filter-period正式表不是",
    EXPECTED_FILTER_PERIODS,
    "行；实际=",
    nrow(
      FP07
    )
  )
}


dup_fp07 <- FP07[
  ,
  .N,
  by = .(
    loc_num,
    collect_start,
    collect_end
  )
][
  N >
    1L,
  .N
]


if (
  dup_fp07 >
    0L
) {

  stopf(
    "07F正式filter-period表存在重复period key：",
    dup_fp07
  )
}


# ============================================================
# Exact frozen 07G eligible-monitor universe
#
# Frozen quasi primary sample contains 426,657 rows and
# 133 monitors. Those monitors are the filter-period eligible
# monitor universe used by 07G.
# ============================================================

ELIGIBLE_LOCS <- sort(
  unique(
    Q_PRIMARY$loc_num
  )
)


if (
  length(
    ELIGIBLE_LOCS
  ) !=
    EXPECTED_MODEL_MONITORS
) {

  stopf(
    "07G eligible monitor universe不是",
    EXPECTED_MODEL_MONITORS,
    "站；实际=",
    length(
      ELIGIBLE_LOCS
    )
  )
}


# ============================================================
# EXACT frozen 07G filter-period analysis-base reconstruction
#
# Frozen 07G definition:
#
#   CONUS
#   + loc_num in eligible 133 monitors
#   + beta finite and < 0.04 Bq/m3
#   + official filter-period PM finite and < 35
#   + wind speed finite and < 6
#   + HPBL finite and 150 < HPBL < 3000
#   + c_prop finite
#   + soilm finite
#   + primary overall 20-km upwind Heat finite
#
# Do NOT replace the final condition with fuel-specific Heat.
# The common 07G base must remain exposure-independent.
# ============================================================

FP07[
  ,
  flag_07g_filter_analysis_base_reproduced :=

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
        FROZEN_PRIMARY_FP
      )
    )
]


N_07F_PRIMARY_FLAG <- FP07[
  flag_07f_filter_period_primary_base %in%
    TRUE,
  .N
]


N_07G_FILTER_BASE <- FP07[
  flag_07g_filter_analysis_base_reproduced %in%
    TRUE,
  .N
]


FILTER_BASE <- FP07[
  flag_07g_filter_analysis_base_reproduced %in%
    TRUE
]


# ------------------------------------------------------------
# Shared internal model variable name.
#
# This does NOT change the PM definition.
#
# filter-period:
#   pm_c = official filter-period PM
#
# quasi-daily:
#   pm_c = same-day PM
#
# The two branches therefore preserve the exact frozen 07G
# timing definitions.
# ------------------------------------------------------------

FILTER_BASE[
  ,
  pm_c :=
    pm_filter_period_official
]


# ------------------------------------------------------------
# Hard reproduction checks
# ------------------------------------------------------------

if (
  nrow(
    FILTER_BASE
  ) !=
    EXPECTED_FILTER_FORMAL_ROWS
) {

  stopf(
    "冻结07G filter-period正式样本复现失败：",
    "expected=",
    EXPECTED_FILTER_FORMAL_ROWS,
    "; actual=",
    nrow(
      FILTER_BASE
    )
  )
}


if (
  uniqueN(
    FILTER_BASE$loc_num
  ) !=
    EXPECTED_MODEL_MONITORS
) {

  stopf(
    "冻结07G filter-period正式样本monitor数不是",
    EXPECTED_MODEL_MONITORS,
    "；实际=",
    uniqueN(
      FILTER_BASE$loc_num
    )
  )
}


dup_filter <- FILTER_BASE[
  ,
  .N,
  by = .(
    loc_num,
    collect_start,
    collect_end
  )
][
  N >
    1L,
  .N
]


if (
  dup_filter >
    0L
) {

  stopf(
    "冻结07G filter-period正式样本存在重复period key。"
  )
}


# ------------------------------------------------------------
# Audit
# ------------------------------------------------------------

FILTER_LOCATOR_AUDIT <- data.table(

  item = c(

    "07F_filter_period_total_rows",

    "07F_filter_period_unique_monitors",

    "07F_flag_filter_period_primary_base_rows",

    "07G_eligible_monitor_count",

    "07G_filter_analysis_base_reproduced_rows",

    "07G_filter_analysis_base_reproduced_monitors",

    "duplicate_period_keys",

    "official_PM_field",

    "overall_primary_filter_exposure"
  ),

  value = c(

    as.character(
      nrow(
        FP07
      )
    ),

    as.character(
      uniqueN(
        FP07$loc_num
      )
    ),

    as.character(
      N_07F_PRIMARY_FLAG
    ),

    as.character(
      length(
        ELIGIBLE_LOCS
      )
    ),

    as.character(
      N_07G_FILTER_BASE
    ),

    as.character(
      uniqueN(
        FILTER_BASE$loc_num
      )
    ),

    as.character(
      dup_filter
    ),

    "pm_filter_period_official",

    FROZEN_PRIMARY_FP
  )
)


fwrite(
  FILTER_LOCATOR_AUDIT,
  OUT_FILTER_LOCATOR
)


FILTER_BASE_SOURCE <-
  "exact_frozen_07G_filter_analysis_base_reproduction"


cat(
  "\n============================================================\n",
  "07G filter-period formal base reproduction\n",
  "============================================================\n",
  "07F filter periods = ",
  fmt_int(
    nrow(
      FP07
    )
  ),
  "\n",
  "07F primary-base flag rows = ",
  fmt_int(
    N_07F_PRIMARY_FLAG
  ),
  "\n",
  "Eligible monitors = ",
  length(
    ELIGIBLE_LOCS
  ),
  "\n",
  "07G reproduced formal base = ",
  fmt_int(
    nrow(
      FILTER_BASE
    )
  ),
  " rows / ",
  uniqueN(
    FILTER_BASE$loc_num
  ),
  " monitors\n",
  "Official filter PM = pm_filter_period_official\n",
  "07H02 working PM variable = pm_c\n",
  "Filter formal base reproduction: PASS\n",
  sep = ""
)


rm(
  FP07
)

invisible(
  gc()
)

# 7. Model formulas
#
# Same architecture as frozen 07G:
#
# gross beta
# ~ exposure
# + c_prop
# + surface U
# + sunspots
# + PM2.5
# + inverse HPBL
# + RH
# + soil moisture
# + wind
# + latitude
# + year + year^2
# + temperature + temperature^2 + temperature^3
# + sunspots × latitude
# + surface U × soil moisture
# + exposure × wind
# + random monitor intercept
# ============================================================

FORM_SINGLE <- as.formula(
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
    "temp_centered + I(temp_centered^2) + I(temp_centered^3) + ",
    "sunspot_monthly:lat_diff_35 + ",
    "surface_u_ppm_article:soilm + ",
    "exposure_scaled:vel + ",
    "(1 | loc_num_factor)"
  )
)


FORM_JOINT <- as.formula(
  paste0(
    "gross_beta_mBq_m3 ~ ",
    "count_scaled + ",
    "heat_scaled + ",
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
    "count_scaled:vel + ",
    "heat_scaled:vel + ",
    "(1 | loc_num_factor)"
  )
)


MODEL_BASE_VARS <- c(
  "loc_num",
  "gross_beta_mBq_m3",
  "c_prop",
  "surface_u_ppm_article",
  "sunspot_monthly",
  "pm_c",
  "inv_hpbl",
  "rhum",
  "soilm",
  "vel",
  "lat_diff_35",
  "year_model",
  "temp"
)


# ============================================================
# 8. Effect extraction
# ============================================================

extract_effect <- function(
  fit,
  main_term,
  scale_value,
  raw_increment,
  sd_raw
) {

  beta <- fixef(
    fit
  )


  V <- as.matrix(
    vcov(
      fit
    )
  )


  int_term <- find_interaction_term(
    names(beta),
    main_term,
    "vel"
  )


  if (
    !main_term %in%
      names(beta) ||
    is.na(
      int_term
    ) ||
    !int_term %in%
      names(beta)
  ) {

    return(
      list(
        ok = FALSE
      )
    )
  }


  cvec <- rep(
    0,
    length(beta)
  )

  names(cvec) <- names(beta)


  cvec[
    main_term
  ] <- 1


  cvec[
    int_term
  ] <- EFFECT_WIND_MS


  est_scaled <- sum(
    cvec *
      beta
  )


  vv <- as.numeric(
    t(cvec) %*%
      V %*%
      cvec
  )


  if (
    !is.finite(vv) ||
    vv <
      0
  ) {

    return(
      list(
        ok = FALSE
      )
    )
  }


  se_scaled <- sqrt(
    vv
  )


  fac_raw <-
    raw_increment /
    scale_value


  est_raw <-
    est_scaled *
    fac_raw


  se_raw <-
    se_scaled *
    fac_raw


  ci_low <-
    est_raw -
    1.96 *
      se_raw


  ci_high <-
    est_raw +
    1.96 *
      se_raw


  p_wald <- if (
    is.finite(
      se_raw
    ) &&
    se_raw >
      0
  ) {

    2 *
      pnorm(
        abs(
          est_raw /
            se_raw
        ),
        lower.tail = FALSE
      )

  } else {

    NA_real_
  }


  effect_1sd <- if (
    is.finite(
      sd_raw
    )
  ) {

    est_scaled *
      (
        sd_raw /
          scale_value
      )

  } else {

    NA_real_
  }


  list(
    ok =
      TRUE,

    estimate =
      est_raw,

    se =
      se_raw,

    low =
      ci_low,

    high =
      ci_high,

    p =
      p_wald,

    effect_1sd =
      effect_1sd,

    main_term =
      main_term,

    interaction_term =
      int_term
  )
}


# ============================================================
# 9. Result-row helper
# ============================================================

result_row <- function(
  meta,
  term,
  status,
  n_rows,
  n_monitors,
  n_unique_exposure,
  scale_value = NA_real_,
  scale_method = NA_character_,
  sd_raw = NA_real_,
  estimate = NA_real_,
  se = NA_real_,
  low = NA_real_,
  high = NA_real_,
  p = NA_real_,
  effect_1sd = NA_real_,
  singular = NA,
  convergence_message = "",
  warnings = "",
  elapsed_sec = NA_real_,
  note = ""
) {

  ti <- term_info(
    term
  )


  data.table(
    model_version =
      MODEL_VERSION,

    spec_id =
      meta$spec_id,

    branch =
      meta$branch,

    fuel =
      meta$fuel,

    direction =
      meta$direction,

    radius_km =
      meta$radius_km,

    model_type =
      meta$model_type,

    term =
      term,

    interpretation =
      interpretation_text(
        meta$model_type,
        term
      ),

    status =
      status,

    n_rows =
      n_rows,

    n_monitors =
      n_monitors,

    n_unique_exposure =
      n_unique_exposure,

    scale_value =
      scale_value,

    scale_method =
      scale_method,

    raw_increment =
      ti$raw_increment,

    raw_unit =
      ti$raw_unit,

    effect_label =
      ti$effect_label,

    wind_speed_ms =
      EFFECT_WIND_MS,

    estimate_mBq_m3 =
      estimate,

    se_mBq_m3 =
      se,

    ci95_low_mBq_m3 =
      low,

    ci95_high_mBq_m3 =
      high,

    p_wald =
      p,

    sd_raw_exposure =
      sd_raw,

    estimate_per_1SD_mBq_m3 =
      effect_1sd,

    singular =
      singular,

    convergence_message =
      convergence_message,

    warnings =
      warnings,

    elapsed_sec =
      elapsed_sec,

    note =
      note
  )
}


# ============================================================
# 10. Fit single-exposure model
# ============================================================

fit_single <- function(
  DT,
  exposure_col,
  term,
  meta
) {

  keep_cols <- unique(
    c(
      MODEL_BASE_VARS,
      exposure_col
    )
  )


  d <- copy(
    DT[
      ,
      ..keep_cols
    ]
  )


  setnames(
    d,
    exposure_col,
    "exposure_raw"
  )


  numeric_cols <- setdiff(
    MODEL_BASE_VARS,
    "loc_num"
  )


  for (
    cc in numeric_cols
  ) {

    d[
      ,
      (cc) :=
        safe_num(
          get(
            cc
          )
        )
    ]
  }


  d[
    ,
    exposure_raw :=
      safe_num(
        exposure_raw
      )
  ]


  d[
    ,
    `:=`(
      temp_centered =
        temp -
          273.15,

      pm_model_centered =
        pm_c -
          10
    )
  ]


  needed <- c(
    "loc_num",
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
    "temp_centered"
  )


  d <- d[
    complete.cases(
      d[
        ,
        ..needed
      ]
    )
  ]


  nr <- nrow(
    d
  )

  nm <- uniqueN(
    d$loc_num
  )

  nu <- uniqueN(
    d$exposure_raw
  )


  if (
    nr <
      1000L ||
    nm <
      2L ||
    nu <
      2L
  ) {

    return(
      result_row(
        meta,
        term,
        "SKIP",
        nr,
        nm,
        nu,
        note =
          "Insufficient rows/monitors/exposure variation"
      )
    )
  }


  sc <- choose_scale(
    d$exposure_raw
  )


  if (
    !sc$ok
  ) {

    return(
      result_row(
        meta,
        term,
        "SKIP",
        nr,
        nm,
        nu,
        scale_method =
          sc$method,
        note =
          "Unable to construct numerical exposure scale"
      )
    )
  }


  sd_raw <- sd(
    d$exposure_raw
  )


  d[
    ,
    `:=`(
      exposure_scaled =
        exposure_raw /
          sc$value,

      loc_num_factor =
        factor(
          loc_num
        )
    )
  ]


  warning_text <- character()


  t0 <- proc.time()[[
    "elapsed"
  ]]


  fit <- tryCatch(

    withCallingHandlers(

      lmer(
        FORM_SINGLE,
        data = d,
        REML = TRUE,
        na.action = na.fail,
        control = make_control()
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


  elapsed <- proc.time()[[
    "elapsed"
  ]] -
    t0


  if (
    inherits(
      fit,
      "fit_error"
    )
  ) {

    return(
      result_row(
        meta,
        term,
        "FAIL",
        nr,
        nm,
        nu,
        scale_value =
          sc$value,
        scale_method =
          sc$method,
        sd_raw =
          sd_raw,
        warnings =
          collapse_messages(
            warning_text
          ),
        elapsed_sec =
          elapsed,
        note =
          paste0(
            "lmer error: ",
            fit$message
          )
      )
    )
  }


  ti <- term_info(
    term
  )


  eff <- extract_effect(
    fit =
      fit,

    main_term =
      "exposure_scaled",

    scale_value =
      sc$value,

    raw_increment =
      ti$raw_increment,

    sd_raw =
      sd_raw
  )


  conv <- collapse_messages(
    fit@optinfo$conv$lme4$messages
  )


  singular <- isSingular(
    fit,
    tol = 1e-4
  )


  if (
    !eff$ok
  ) {

    return(
      result_row(
        meta,
        term,
        "SKIP",
        nr,
        nm,
        nu,
        scale_value =
          sc$value,
        scale_method =
          sc$method,
        sd_raw =
          sd_raw,
        singular =
          singular,
        convergence_message =
          conv,
        warnings =
          collapse_messages(
            warning_text
          ),
        elapsed_sec =
          elapsed,
        note =
          "Exposure main/interaction term unavailable after model matrix rank handling"
      )
    )
  }


  result_row(
    meta,
    term,
    "PASS",
    nr,
    nm,
    nu,
    scale_value =
      sc$value,
    scale_method =
      sc$method,
    sd_raw =
      sd_raw,
    estimate =
      eff$estimate,
    se =
      eff$se,
    low =
      eff$low,
    high =
      eff$high,
    p =
      eff$p,
    effect_1sd =
      eff$effect_1sd,
    singular =
      singular,
    convergence_message =
      conv,
    warnings =
      collapse_messages(
        warning_text
      ),
    elapsed_sec =
      elapsed
  )
}


# ============================================================
# 11. Fit joint Count + Heat model
# ============================================================

fit_joint <- function(
  DT,
  count_col,
  heat_col,
  meta
) {

  keep_cols <- unique(
    c(
      MODEL_BASE_VARS,
      count_col,
      heat_col
    )
  )


  d <- copy(
    DT[
      ,
      ..keep_cols
    ]
  )


  setnames(
    d,
    c(
      count_col,
      heat_col
    ),
    c(
      "count_raw",
      "heat_raw"
    )
  )


  numeric_cols <- setdiff(
    MODEL_BASE_VARS,
    "loc_num"
  )


  for (
    cc in numeric_cols
  ) {

    d[
      ,
      (cc) :=
        safe_num(
          get(
            cc
          )
        )
    ]
  }


  d[
    ,
    `:=`(
      count_raw =
        safe_num(
          count_raw
        ),

      heat_raw =
        safe_num(
          heat_raw
        ),

      temp_centered =
        temp -
          273.15,

      pm_model_centered =
        pm_c -
          10
    )
  ]


  needed <- c(
    "loc_num",
    "gross_beta_mBq_m3",
    "count_raw",
    "heat_raw",
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
    "temp_centered"
  )


  d <- d[
    complete.cases(
      d[
        ,
        ..needed
      ]
    )
  ]


  nr <- nrow(
    d
  )

  nm <- uniqueN(
    d$loc_num
  )

  nuc <- uniqueN(
    d$count_raw
  )

  nuh <- uniqueN(
    d$heat_raw
  )


  if (
    nr <
      1000L ||
    nm <
      2L ||
    nuc <
      2L ||
    nuh <
      2L
  ) {

    return(
      rbindlist(
        list(
          result_row(
            meta,
            "count",
            "SKIP",
            nr,
            nm,
            nuc,
            note =
              "Insufficient rows/monitors/joint exposure variation"
          ),

          result_row(
            meta,
            "heat",
            "SKIP",
            nr,
            nm,
            nuh,
            note =
              "Insufficient rows/monitors/joint exposure variation"
          )
        ),
        use.names = TRUE,
        fill = TRUE
      )
    )
  }


  sc_count <- choose_scale(
    d$count_raw
  )

  sc_heat <- choose_scale(
    d$heat_raw
  )


  if (
    !sc_count$ok ||
    !sc_heat$ok
  ) {

    return(
      rbindlist(
        list(
          result_row(
            meta,
            "count",
            "SKIP",
            nr,
            nm,
            nuc,
            scale_value =
              sc_count$value,
            scale_method =
              sc_count$method,
            note =
              "Unable to construct joint numerical exposure scale"
          ),

          result_row(
            meta,
            "heat",
            "SKIP",
            nr,
            nm,
            nuh,
            scale_value =
              sc_heat$value,
            scale_method =
              sc_heat$method,
            note =
              "Unable to construct joint numerical exposure scale"
          )
        ),
        use.names = TRUE,
        fill = TRUE
      )
    )
  }


  sd_count <- sd(
    d$count_raw
  )

  sd_heat <- sd(
    d$heat_raw
  )


  d[
    ,
    `:=`(
      count_scaled =
        count_raw /
          sc_count$value,

      heat_scaled =
        heat_raw /
          sc_heat$value,

      loc_num_factor =
        factor(
          loc_num
        )
    )
  ]


  warning_text <- character()


  t0 <- proc.time()[[
    "elapsed"
  ]]


  fit <- tryCatch(

    withCallingHandlers(

      lmer(
        FORM_JOINT,
        data = d,
        REML = TRUE,
        na.action = na.fail,
        control = make_control()
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


  elapsed <- proc.time()[[
    "elapsed"
  ]] -
    t0


  if (
    inherits(
      fit,
      "fit_error"
    )
  ) {

    return(
      rbindlist(
        list(
          result_row(
            meta,
            "count",
            "FAIL",
            nr,
            nm,
            nuc,
            scale_value =
              sc_count$value,
            scale_method =
              sc_count$method,
            sd_raw =
              sd_count,
            warnings =
              collapse_messages(
                warning_text
              ),
            elapsed_sec =
              elapsed,
            note =
              paste0(
                "joint lmer error: ",
                fit$message
              )
          ),

          result_row(
            meta,
            "heat",
            "FAIL",
            nr,
            nm,
            nuh,
            scale_value =
              sc_heat$value,
            scale_method =
              sc_heat$method,
            sd_raw =
              sd_heat,
            warnings =
              collapse_messages(
                warning_text
              ),
            elapsed_sec =
              elapsed,
            note =
              paste0(
                "joint lmer error: ",
                fit$message
              )
          )
        ),
        use.names = TRUE,
        fill = TRUE
      )
    )
  }


  conv <- collapse_messages(
    fit@optinfo$conv$lme4$messages
  )


  singular <- isSingular(
    fit,
    tol = 1e-4
  )


  count_info <- term_info(
    "count"
  )

  heat_info <- term_info(
    "heat"
  )


  ec <- extract_effect(
    fit =
      fit,
    main_term =
      "count_scaled",
    scale_value =
      sc_count$value,
    raw_increment =
      count_info$raw_increment,
    sd_raw =
      sd_count
  )


  eh <- extract_effect(
    fit =
      fit,
    main_term =
      "heat_scaled",
    scale_value =
      sc_heat$value,
    raw_increment =
      heat_info$raw_increment,
    sd_raw =
      sd_heat
  )


  RC <- if (
    ec$ok
  ) {

    result_row(
      meta,
      "count",
      "PASS",
      nr,
      nm,
      nuc,
      scale_value =
        sc_count$value,
      scale_method =
        sc_count$method,
      sd_raw =
        sd_count,
      estimate =
        ec$estimate,
      se =
        ec$se,
      low =
        ec$low,
      high =
        ec$high,
      p =
        ec$p,
      effect_1sd =
        ec$effect_1sd,
      singular =
        singular,
      convergence_message =
        conv,
      warnings =
        collapse_messages(
          warning_text
        ),
      elapsed_sec =
        elapsed
    )

  } else {

    result_row(
      meta,
      "count",
      "SKIP",
      nr,
      nm,
      nuc,
      scale_value =
        sc_count$value,
      scale_method =
        sc_count$method,
      sd_raw =
        sd_count,
      singular =
        singular,
      convergence_message =
        conv,
      warnings =
        collapse_messages(
          warning_text
        ),
      elapsed_sec =
        elapsed,
      note =
        "Count term unavailable after rank handling"
    )
  }


  RH <- if (
    eh$ok
  ) {

    result_row(
      meta,
      "heat",
      "PASS",
      nr,
      nm,
      nuh,
      scale_value =
        sc_heat$value,
      scale_method =
        sc_heat$method,
      sd_raw =
        sd_heat,
      estimate =
        eh$estimate,
      se =
        eh$se,
      low =
        eh$low,
      high =
        eh$high,
      p =
        eh$p,
      effect_1sd =
        eh$effect_1sd,
      singular =
        singular,
      convergence_message =
        conv,
      warnings =
        collapse_messages(
          warning_text
        ),
      elapsed_sec =
        elapsed
    )

  } else {

    result_row(
      meta,
      "heat",
      "SKIP",
      nr,
      nm,
      nuh,
      scale_value =
        sc_heat$value,
      scale_method =
        sc_heat$method,
      sd_raw =
        sd_heat,
      singular =
        singular,
      convergence_message =
        conv,
      warnings =
        collapse_messages(
          warning_text
        ),
      elapsed_sec =
        elapsed,
      note =
        "Heat term unavailable after rank handling"
    )
  }


  rbindlist(
    list(
      RC,
      RH
    ),
    use.names = TRUE,
    fill = TRUE
  )
}


# ============================================================
# 12. Daily 07H01 fuel reader
# ============================================================

exposure_columns <- function(
  fuel
) {

  unlist(
    lapply(
      DIRECTIONS,
      function(dir) {

        unlist(
          lapply(
            RADII,
            function(r) {

              c(
                paste0(
                  "campd_",
                  fuel,
                  "_",
                  dir,
                  "_n_plants_",
                  r,
                  "km_90deg"
                ),

                paste0(
                  "campd_",
                  fuel,
                  "_",
                  dir,
                  "_heat_observed_",
                  r,
                  "km_90deg"
                )
              )
            }
          )
        )
      }
    )
  )
}


read_fuel_daily <- function(
  fuel
) {

  excols <- exposure_columns(
    fuel
  )


  parts <- lapply(
    H01_FILES,
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
        excols
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
          "07H01 daily缺少字段：",
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
      loc_num =
        safe_int(
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
    nrow(D) !=
      EXPECTED_UNIQUE_MONITOR_DAYS
  ) {

    stopf(
      fuel,
      " daily monitor-day数不是",
      EXPECTED_UNIQUE_MONITOR_DAYS,
      "，实际=",
      nrow(D)
    )
  }


  dup <- D[
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
      " daily exposure存在重复loc_num+date。"
    )
  }


  D[
    ,
    .joined_07h :=
      TRUE
  ]


  D
}


# ============================================================
# 13. Count-Heat correlation / collinearity
# ============================================================

build_correlation_row <- function(
  DT,
  count_col,
  heat_col,
  branch,
  fuel,
  direction,
  radius
) {

  x <- safe_num(
    DT[[
      count_col
    ]]
  )

  y <- safe_num(
    DT[[
      heat_col
    ]]
  )


  both <-
    is.finite(x) &
    is.finite(y)


  xb <- x[both]
  yb <- y[both]


  rp <- safe_cor(
    xb,
    yb,
    "pearson"
  )


  rs <- safe_cor(
    xb,
    yb,
    "spearman"
  )


  vif <- if (
    is.finite(rp) &&
    abs(rp) <
      1
  ) {

    1 /
      (
        1 -
          rp^2
      )

  } else if (
    is.finite(rp)
  ) {

    Inf

  } else {

    NA_real_
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
      length(x),

    n_count_finite =
      sum(
        is.finite(x)
      ),

    n_heat_finite =
      sum(
        is.finite(y)
      ),

    n_both_finite =
      sum(
        both
      ),

    n_positive_count =
      sum(
        is.finite(x) &
          x >
            0
      ),

    n_positive_heat =
      sum(
        is.finite(y) &
          y >
            0
      ),

    count_zero_pct =
      if (
        sum(
          is.finite(x)
        ) >
          0
      ) {

        100 *
          mean(
            x[
              is.finite(x)
            ] ==
              0
          )

      } else {

        NA_real_
      },

    heat_zero_pct =
      if (
        sum(
          is.finite(y)
        ) >
          0
      ) {

        100 *
          mean(
            y[
              is.finite(y)
            ] ==
              0
          )

      } else {

        NA_real_
      },

    pearson_r =
      rp,

    spearman_rho =
      rs,

    two_predictor_vif =
      vif,

    collinearity_flag =
      fifelse(
        is.infinite(vif),
        "PERFECT_COLLINEARITY_VIF_INF",
        fifelse(
          is.finite(vif) &
            vif >=
              10,
          "HIGH_VIF_GE10",
          fifelse(
            is.finite(vif) &
              vif >=
                5,
            "MODERATE_VIF_5_TO_10",
            fifelse(
              is.finite(vif),
              "NO_HIGH_VIF_FLAG",
              "VIF_NOT_ESTIMABLE"
            )
          )
        )
      )
  )
}


# ============================================================
# 14. Filter-period strict aggregation
#
# EXACT rule:
#
# collect_start <= date < collect_end
#
# For each exposure:
#   every expected day must exist
#   every day must be finite
#   otherwise period exposure = NA
#
# If complete:
#   arithmetic mean
# ============================================================

build_filter_exposure <- function(
  fuel,
  DAILY
) {

  excols <- exposure_columns(
    fuel
  )


  PD <- merge(
    PERIOD_DAY_KEYS,
    DAILY,
    by = c(
      "loc_num",
      "date"
    ),
    all.x = TRUE,
    sort = FALSE
  )


  if (
    nrow(PD) !=
      EXPECTED_07F_ROWS
  ) {

    stopf(
      fuel,
      " filter-period day join改变了07F行数。"
    )
  }


  missing_join <- PD[
    is.na(
      .joined_07h
    ),
    .N
  ]


  if (
    missing_join >
      0L
  ) {

    stopf(
      fuel,
      " 有",
      missing_join,
      "个07F period-day无法匹配07H01 daily exposure。"
    )
  }


  PD[
    ,
    .joined_07h :=
      NULL
  ]


  invalid_period_day <- PD[
    is.na(
      collect_start
    ) |
      is.na(
        collect_end
      ) |
      date <
        collect_start |
      date >=
        collect_end |
      collect_end <=
        collect_start,
    .N
  ]


  if (
    invalid_period_day >
      0L
  ) {

    stopf(
      fuel,
      "存在违反[collect_start,collect_end)的period-day。"
    )
  }


  FP <- PD[
    ,
    {

      nn <- .N

      expected_days <-
        as.integer(
          collect_end[1] -
            collect_start[1]
        )


      vals <- lapply(
        .SD,
        function(v) {

          vv <- safe_num(
            v
          )


          if (
            nn ==
              expected_days &&
            length(vv) ==
              expected_days &&
            all(
              is.finite(
                vv
              )
            )
          ) {

            mean(
              vv
            )

          } else {

            NA_real_
          }
        }
      )


      c(
        list(
          n_days =
            nn,

          expected_days =
            expected_days
        ),
        vals
      )
    },
    by = .(
      loc_num,
      collect_start,
      collect_end
    ),
    .SDcols = excols
  ]


  if (
    nrow(FP) !=
      EXPECTED_FILTER_PERIODS
  ) {

    stopf(
      fuel,
      " filter-period聚合后period数不是",
      EXPECTED_FILTER_PERIODS,
      "。"
    )
  }


  bad_n_days <- FP[
    n_days !=
      expected_days,
    .N
  ]


  if (
    bad_n_days >
      0L
  ) {

    stopf(
      fuel,
      " 有",
      bad_n_days,
      "个filter period没有完整展开所有日期。"
    )
  }


  out_file <- file.path(
    FILTER_EXP_DIR,
    paste0(
      "07H_filter_period_exposure_",
      fuel,
      ".csv.gz"
    )
  )


  fwrite(
    FP,
    out_file
  )


  FP[
    ,
    .filter_exp_joined :=
      TRUE
  ]


  FM <- merge(
    FILTER_BASE,
    FP,
    by = c(
      "loc_num",
      "collect_start",
      "collect_end"
    ),
    all.x = TRUE,
    sort = FALSE
  )


  if (
    nrow(FM) !=
      EXPECTED_FILTER_FORMAL_ROWS
  ) {

    stopf(
      fuel,
      " filter formal exposure merge改变了107,407行。"
    )
  }


  miss_formal <- FM[
    is.na(
      .filter_exp_joined
    ),
    .N
  ]


  if (
    miss_formal >
      0L
  ) {

    stopf(
      fuel,
      " formal filter base有",
      miss_formal,
      "个period无法匹配07H filter exposure。"
    )
  }


  FM[
    ,
    .filter_exp_joined :=
      NULL
  ]


  FM
}


# ============================================================
# 15. Checkpoint helpers
# ============================================================

checkpoint_path <- function(
  spec_id
) {

  clean <- gsub(
    "[^A-Za-z0-9_]+",
    "_",
    spec_id
  )


  file.path(
    CHECKPOINT_DIR,
    paste0(
      "cp__",
      clean,
      ".csv"
    )
  )
}


checkpoint_read <- function(
  spec_id
) {

  f <- checkpoint_path(
    spec_id
  )


  if (
    FORCE_RESTART ||
    !file.exists(f)
  ) {

    return(
      NULL
    )
  }


  z <- tryCatch(
    fread(
      f,
      showProgress = FALSE
    ),
    error = function(e) {
      NULL
    }
  )


  if (
    is.null(z) ||
    nrow(z) ==
      0L ||
    !"model_version" %in%
      names(z) ||
    any(
      z$model_version !=
        MODEL_VERSION
    ) ||
    any(
      z$spec_id !=
        spec_id
    ) ||
    any(
      z$status ==
        "FAIL"
    )
  ) {

    return(
      NULL
    )
  }


  z
}


checkpoint_write <- function(
  z
) {

  f <- checkpoint_path(
    z$spec_id[1]
  )


  fwrite(
    z,
    f
  )
}


# ============================================================
# 16. Main modeling loop
# ============================================================

RESULT_LIST <- list()

COR_LIST <- list()

result_i <- 0L

cor_i <- 0L

spec_counter <- 0L


SAMPLE_SUMMARY <- data.table(

  branch =
    c(
      "quasi_daily",
      "filter_period"
    ),

  n_rows =
    c(
      nrow(
        Q_PRIMARY
      ),
      nrow(
        FILTER_BASE
      )
    ),

  n_monitors =
    c(
      uniqueN(
        Q_PRIMARY$loc_num
      ),
      uniqueN(
        FILTER_BASE$loc_num
      )
    ),

  source =
    c(
      paste0(
        "07F primary flag: ",
        PRIMARY_FLAG
      ),
      paste0(
        FILTER_BASE_FILE,
        " [",
        FILTER_BASE_SOURCE,
        "]"
      )
    )
)


fwrite(
  SAMPLE_SUMMARY,
  OUT_SAMPLE
)


cat(
  "\n============================================================\n",
  "07H02 MODELING START\n",
  "============================================================\n",
  "Expected model specs = ",
  EXPECTED_MODEL_SPECS,
  "\nExpected effect rows = ",
  EXPECTED_EFFECT_ROWS,
  "\nQuasi base = ",
  fmt_int(
    nrow(
      Q_PRIMARY
    )
  ),
  " rows / ",
  uniqueN(
    Q_PRIMARY$loc_num
  ),
  " monitors\n",
  "Filter base = ",
  fmt_int(
    nrow(
      FILTER_BASE
    )
  ),
  " rows / ",
  uniqueN(
    FILTER_BASE$loc_num
  ),
  " monitors\n",
  sep = ""
)


for (
  fuel in FUELS
) {

  cat(
    "\n============================================================\n",
    "FUEL = ",
    fuel,
    "\n",
    "============================================================\n",
    sep = ""
  )


  DAILY <- read_fuel_daily(
    fuel
  )


  QMODEL <- merge(
    Q_PRIMARY,
    DAILY,
    by = c(
      "loc_num",
      "date"
    ),
    all.x = TRUE,
    sort = FALSE
  )


  if (
    nrow(QMODEL) !=
      EXPECTED_QUASI_PRIMARY_ROWS
  ) {

    stopf(
      fuel,
      " quasi exposure join改变了426,657行。"
    )
  }


  miss_q_join <- QMODEL[
    is.na(
      .joined_07h
    ),
    .N
  ]


  if (
    miss_q_join >
      0L
  ) {

    stopf(
      fuel,
      " quasi base有",
      miss_q_join,
      "行无法匹配07H01。"
    )
  }


  QMODEL[
    ,
    .joined_07h :=
      NULL
  ]


  FMODEL <- build_filter_exposure(
    fuel,
    DAILY
  )


  rm(
    DAILY
  )

  invisible(
    gc()
  )


  for (
    branch in BRANCHES
  ) {

    BDT <- if (
      branch ==
        "filter_period"
    ) {

      FMODEL

    } else {

      QMODEL
    }


    for (
      direction in DIRECTIONS
    ) {

      for (
        radius in RADII
      ) {

        count_col <- paste0(
          "campd_",
          fuel,
          "_",
          direction,
          "_n_plants_",
          radius,
          "km_90deg"
        )


        heat_col <- paste0(
          "campd_",
          fuel,
          "_",
          direction,
          "_heat_observed_",
          radius,
          "km_90deg"
        )


        # ----------------------------------------------------
        # Count-Heat correlation audit
        # ----------------------------------------------------

        cor_i <- cor_i + 1L


        COR_LIST[[cor_i]] <-
          build_correlation_row(
            DT =
              BDT,
            count_col =
              count_col,
            heat_col =
              heat_col,
            branch =
              branch,
            fuel =
              fuel,
            direction =
              direction,
            radius =
              radius
          )


        # ----------------------------------------------------
        # A / B / C models
        # ----------------------------------------------------

        for (
          model_type in MODEL_TYPES
        ) {

          spec_counter <-
            spec_counter +
            1L


          spec_id <- paste(
            branch,
            fuel,
            direction,
            paste0(
              radius,
              "km"
            ),
            model_type,
            sep = "__"
          )


          meta <- list(
            spec_id =
              spec_id,
            branch =
              branch,
            fuel =
              fuel,
            direction =
              direction,
            radius_km =
              radius,
            model_type =
              model_type
          )


          old <- checkpoint_read(
            spec_id
          )


          if (
            !is.null(old)
          ) {

            RES <- old

            cat(
              "[",
              spec_counter,
              "/",
              EXPECTED_MODEL_SPECS,
              "] REUSE ",
              spec_id,
              "\n",
              sep = ""
            )


          } else {

            cat(
              "[",
              spec_counter,
              "/",
              EXPECTED_MODEL_SPECS,
              "] FIT   ",
              spec_id,
              "\n",
              sep = ""
            )


            RES <- if (
              model_type ==
                "count_only"
            ) {

              fit_single(
                DT =
                  BDT,
                exposure_col =
                  count_col,
                term =
                  "count",
                meta =
                  meta
              )


            } else if (
              model_type ==
                "heat_only"
            ) {

              fit_single(
                DT =
                  BDT,
                exposure_col =
                  heat_col,
                term =
                  "heat",
                meta =
                  meta
              )


            } else {

              fit_joint(
                DT =
                  BDT,
                count_col =
                  count_col,
                heat_col =
                  heat_col,
                meta =
                  meta
              )
            }


            checkpoint_write(
              RES
            )
          }


          result_i <- result_i + 1L

          RESULT_LIST[[result_i]] <-
            RES
        }
      }
    }
  }


  COR_TEMP <- rbindlist(
    COR_LIST,
    use.names = TRUE,
    fill = TRUE
  )


  fwrite(
    COR_TEMP,
    OUT_COR
  )


  rm(
    QMODEL,
    FMODEL
  )

  invisible(
    gc()
  )
}


# ============================================================
# 17. Final combine
# ============================================================

RESULT <- rbindlist(
  RESULT_LIST,
  use.names = TRUE,
  fill = TRUE
)


COR_RESULT <- rbindlist(
  COR_LIST,
  use.names = TRUE,
  fill = TRUE
)


setorder(
  RESULT,
  branch,
  fuel,
  direction,
  radius_km,
  model_type,
  term
)


setorder(
  COR_RESULT,
  branch,
  fuel,
  direction,
  radius_km
)


fwrite(
  RESULT,
  OUT_EFFECTS
)


fwrite(
  COR_RESULT,
  OUT_COR
)


STATUS <- RESULT[
  ,
  .(
    n_effect_rows =
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

    status_summary =
      paste(
        sort(
          unique(
            status
          )
        ),
        collapse = "|"
      ),

    n_rows_min =
      min(
        n_rows,
        na.rm = TRUE
      ),

    n_rows_max =
      max(
        n_rows,
        na.rm = TRUE
      ),

    n_monitors_min =
      min(
        n_monitors,
        na.rm = TRUE
      ),

    singular_any =
      any(
        singular %in%
          TRUE
      ),

    convergence_message_any =
      any(
        !is.na(
          convergence_message
        ) &
          nzchar(
            convergence_message
          )
      ),

    warnings_any =
      any(
        !is.na(
          warnings
        ) &
          nzchar(
            warnings
          )
      )
  ),
  by = .(
    spec_id,
    branch,
    fuel,
    direction,
    radius_km,
    model_type
  )
]


setorder(
  STATUS,
  branch,
  fuel,
  direction,
  radius_km,
  model_type
)


fwrite(
  STATUS,
  OUT_STATUS
)


# ============================================================
# 18. Hard checks
# ============================================================

n_specs <- uniqueN(
  RESULT$spec_id
)


n_effect_rows <- nrow(
  RESULT
)


dup_effect_key <- RESULT[
  ,
  .N,
  by = .(
    spec_id,
    term
  )
][
  N >
    1L,
  .N
]


n_fail <- RESULT[
  status ==
    "FAIL",
  .N
]


n_skip <- RESULT[
  status ==
    "SKIP",
  .N
]


n_pass <- RESULT[
  status ==
    "PASS",
  .N
]


n_singular <- RESULT[
  status ==
    "PASS" &
    singular %in%
      TRUE,
  uniqueN(
    spec_id
  )
]


n_conv <- RESULT[
  status ==
    "PASS" &
    !is.na(
      convergence_message
    ) &
    nzchar(
      convergence_message
    ),
  uniqueN(
    spec_id
  )
]


n_warning <- RESULT[
  status ==
    "PASS" &
    !is.na(
      warnings
    ) &
    nzchar(
      warnings
    ),
  uniqueN(
    spec_id
  )
]


HARD <- data.table(

  check = c(
    "07H01_hard_checks_all_pass",
    "07H01_yearly_files_24",
    "07F_quasi_primary_rows_426657",
    "07F_quasi_primary_monitors_133",
    "07F_filter_formal_rows_107407",
    "07F_filter_formal_monitors_133",
    "filter_period_universe_174200",
    "model_specs_720",
    "effect_rows_960",
    "duplicate_effect_keys_zero",
    "count_heat_correlation_rows_240",
    "fit_error_rows_zero"
  ),

  observed = c(
    as.character(
      H01_CHECK[
        status !=
          "PASS",
        .N
      ]
    ),

    as.character(
      sum(
        file.exists(
          H01_FILES
        )
      )
    ),

    as.character(
      nrow(
        Q_PRIMARY
      )
    ),

    as.character(
      uniqueN(
        Q_PRIMARY$loc_num
      )
    ),

    as.character(
      nrow(
        FILTER_BASE
      )
    ),

    as.character(
      uniqueN(
        FILTER_BASE$loc_num
      )
    ),

    as.character(
      N_PERIODS
    ),

    as.character(
      n_specs
    ),

    as.character(
      n_effect_rows
    ),

    as.character(
      dup_effect_key
    ),

    as.character(
      nrow(
        COR_RESULT
      )
    ),

    as.character(
      n_fail
    )
  ),

  expected = c(
    "0 failures",
    "24",
    "426657",
    "133",
    "107407",
    "133",
    "174200",
    as.character(
      EXPECTED_MODEL_SPECS
    ),
    as.character(
      EXPECTED_EFFECT_ROWS
    ),
    "0",
    as.character(
      EXPECTED_COR_ROWS
    ),
    "0"
  )
)


HARD[
  ,
  status :=
    c(
      ifelse(
        H01_CHECK[
          status !=
            "PASS",
          .N
        ] ==
          0L,
        "PASS",
        "FAIL"
      ),

      ifelse(
        sum(
          file.exists(
            H01_FILES
          )
        ) ==
          24L,
        "PASS",
        "FAIL"
      ),

      ifelse(
        nrow(
          Q_PRIMARY
        ) ==
          EXPECTED_QUASI_PRIMARY_ROWS,
        "PASS",
        "FAIL"
      ),

      ifelse(
        uniqueN(
          Q_PRIMARY$loc_num
        ) ==
          EXPECTED_MODEL_MONITORS,
        "PASS",
        "FAIL"
      ),

      ifelse(
        nrow(
          FILTER_BASE
        ) ==
          EXPECTED_FILTER_FORMAL_ROWS,
        "PASS",
        "FAIL"
      ),

      ifelse(
        uniqueN(
          FILTER_BASE$loc_num
        ) ==
          EXPECTED_MODEL_MONITORS,
        "PASS",
        "FAIL"
      ),

      ifelse(
        N_PERIODS ==
          EXPECTED_FILTER_PERIODS,
        "PASS",
        "FAIL"
      ),

      ifelse(
        n_specs ==
          EXPECTED_MODEL_SPECS,
        "PASS",
        "FAIL"
      ),

      ifelse(
        n_effect_rows ==
          EXPECTED_EFFECT_ROWS,
        "PASS",
        "FAIL"
      ),

      ifelse(
        dup_effect_key ==
          0L,
        "PASS",
        "FAIL"
      ),

      ifelse(
        nrow(
          COR_RESULT
        ) ==
          EXPECTED_COR_ROWS,
        "PASS",
        "FAIL"
      ),

      ifelse(
        n_fail ==
          0L,
        "PASS",
        "FAIL"
      )
    )
]


fwrite(
  HARD,
  OUT_HARD
)


print(
  HARD
)


# ============================================================
# 19. Report
# ============================================================

HARD_FAIL <- HARD[
  status ==
    "FAIL",
  .N
]


REPORT <- c(

  "# 07H02 Fuel-specific Count + Heat cumulative-radius models",

  "",

  paste0(
    "- Generated: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),

  paste0(
    "- Model version: ",
    MODEL_VERSION
  ),

  "- Frozen upstream 07A–07G and 07H01 were read-only.",

  paste0(
    "- 07F quasi primary flag: ",
    PRIMARY_FLAG
  ),

  paste0(
    "- 07F filter formal base: ",
    FILTER_BASE_FILE
  ),

  paste0(
    "- Filter base locator mode: ",
    FILTER_BASE_SOURCE
  ),

  "",

  "## Prespecified framework",

  "",

  "- Fuels: coal, gas, oil, other_fossil, biomass, solid_waste.",
  "- Cumulative radii: 20, 25, 30, 35, 40, 45, 50, 100, 150, 200 km.",
  "- Directions: upwind main analysis and downwind directional negative control.",
  "- Models: Count-only, Heat-only, and joint Count + Heat.",
  "- Count and Heat originate from the same daily CAMPD fuel-specific source layer constructed by 07H01.",
  "- Count effect is reported per +1 plant.",
  "- Heat effect is reported per +10,000 MMBtu.",
  "- Effect at wind speed = 1 m/s.",
  "- One-SD effects are supplementary.",
  "- Joint Count coefficient is conditional on Heat Input.",
  "- Joint Heat coefficient is conditional on plant count.",

  "",

  "## Filter-period rule",

  "",

  "- Period definition: collect_start <= date < collect_end.",
  "- Exposure aggregation: strict arithmetic mean.",
  "- If any expected daily exposure is NA, the corresponding filter-period exposure is NA.",

  "",

  "## Model architecture",

  "",

  "- Same covariate architecture as frozen 07G.",
  "- Temperature is centered as temp - 273.15 before polynomial terms.",
  "- PM2.5 is centered as pm_c - 10.",
  "- Exposure-by-wind interaction is retained.",
  "- Monitor-specific random intercept.",
  "- REML = TRUE; bobyqa; maxfun = 200000.",
  "- Exposure P95 scaling is numerical only; reported effects are converted back to physical units.",

  "",

  "## Completion",

  "",

  paste0(
    "- Expected model specs: ",
    EXPECTED_MODEL_SPECS
  ),

  paste0(
    "- Completed specs: ",
    n_specs
  ),

  paste0(
    "- Effect rows: ",
    n_effect_rows
  ),

  paste0(
    "- PASS effect rows: ",
    n_pass
  ),

  paste0(
    "- SKIP effect rows: ",
    n_skip
  ),

  paste0(
    "- FAIL effect rows: ",
    n_fail
  ),

  paste0(
    "- Singular PASS model specs: ",
    n_singular
  ),

  paste0(
    "- PASS specs with convergence message: ",
    n_conv
  ),

  paste0(
    "- PASS specs with warnings: ",
    n_warning
  ),

  "",

  if (
    HARD_FAIL ==
      0L
  ) {

    "**07H02 MODEL BUILD: PASS**"

  } else {

    "**07H02 MODEL BUILD: FAIL**"
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


cat(
  "\n============================================================\n"
)


if (
  HARD_FAIL ==
    0L
) {

  cat(
    "07H02 MODEL BUILD: PASS\n"
  )

} else {

  cat(
    "07H02 MODEL BUILD: FAIL\n"
  )
}


cat(
  "============================================================\n",
  "Model specs = ",
  n_specs,
  " / ",
  EXPECTED_MODEL_SPECS,
  "\n",
  "Effect rows = ",
  n_effect_rows,
  " / ",
  EXPECTED_EFFECT_ROWS,
  "\n",
  "PASS effect rows = ",
  n_pass,
  "\n",
  "SKIP effect rows = ",
  n_skip,
  "\n",
  "FAIL effect rows = ",
  n_fail,
  "\n",
  "Singular PASS specs = ",
  n_singular,
  "\n",
  "Convergence-message specs = ",
  n_conv,
  "\n",
  "Warning specs = ",
  n_warning,
  "\n",
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
    "07H02存在hard failure。"
  )
}
