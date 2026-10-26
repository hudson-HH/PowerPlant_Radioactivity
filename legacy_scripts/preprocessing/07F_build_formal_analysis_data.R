
#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999
)

# ============================================================
# 07F 正式统计分析数据构建
#
# 只读冻结的07E，不修改07D/07E。
#
# 目标：
#
# A. quasi-daily branch
#    - 保留08/07E原quasi-daily observation structure
#    - 独立重新构造same-day PM2.5
#    - 旧pm_c只作为legacy filter-period PM保留
#    - 构建新的same-day PM分析标志
#
# B. filter-period branch
#    - 一张RadNet filter period只保留一条gross-beta outcome
#    - 将逐日CAMPD exposure在整个采样周期取日均值
#    - 主分析要求该暴露在整个filter period内完整
#    - 正式filter-period PM仍保留原来的filter-period定义
#
# 07F不拟合任何回归模型。
# 07F通过后再进入下一模型阶段。
# ============================================================


# ============================================================
# 0. Package
# ============================================================

if (
  !requireNamespace(
    "data.table",
    quietly = TRUE
  )
) {

  install.packages(
    "data.table",
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(data.table)
})


# ============================================================
# 1. Paths
# ============================================================

BASE <- "/Users/hudson/Desktop/USA 辐射表征研究"

PP_ROOT <- file.path(
  BASE,
  "发电站数据处理"
)

E07_ROOT <- file.path(
  PP_ROOT,
  "03_07E_CAMPD逐日真实源强与模型表连接"
)

E07_YEAR_DIR <- file.path(
  E07_ROOT,
  "yearly"
)

E07_FREEZE <- file.path(
  E07_ROOT,
  "07E_FINAL_FREEZE.md"
)

E07_HIER <- file.path(
  E07_ROOT,
  "07_07E_model_exposure_hierarchy.csv"
)

# ------------------------------------------------------------
# same-day PM独立重构输入
# ------------------------------------------------------------

PM_DAILY_RETIRED <- file.path(
  BASE,
  "混杂因素",
  "PM2.5",
  "PM2.5处理版",
  "02_radnet_pm25_daily_50km_文章方法_2001_2024.csv.gz"
)

PM_LINK_FILE <- file.path(
  BASE,
  "混杂因素",
  "PM2.5",
  "PM2.5处理版",
  "01_radnet_pm25_site_links_50km_文章方法.csv"
)

PM_SITE_DAY_FILE <- file.path(
  BASE,
  "混杂因素",
  "PM2.5",
  "outputs",
  "pm25_site_day_2001_2024.csv"
)

# ------------------------------------------------------------
# 07F输出
# ------------------------------------------------------------

OUT <- file.path(
  PP_ROOT,
  "04_07F_正式统计分析数据构建"
)

FINAL_QUASI_DIR <- file.path(
  OUT,
  "quasi_daily_yearly"
)

BUILD_QUASI_DIR <- file.path(
  OUT,
  ".quasi_daily_building"
)

OUT_PM_AUDIT <- file.path(
  OUT,
  "01_07F_same_day_PM独立重构审计.csv"
)

OUT_PM_MISMATCH <- file.path(
  OUT,
  "02_07F_same_day_PM不一致明细.csv.gz"
)

OUT_FILTER <- file.path(
  OUT,
  "03_07F_filter_period_analysis_base_2001_2024.csv.gz"
)

OUT_FILTER_COV <- file.path(
  OUT,
  "04_07F_filter_period_exposure_coverage.csv"
)

OUT_QC <- file.path(
  OUT,
  "05_07F_dataset_qc_summary.csv"
)

OUT_HARD <- file.path(
  OUT,
  "06_07F_hard_failure_checks.csv"
)

OUT_WARN <- file.path(
  OUT,
  "07_07F_warning_checks.csv"
)

OUT_REPORT <- file.path(
  OUT,
  "08_07F_正式统计分析数据构建报告.md"
)

OUT_SESSION <- file.path(
  OUT,
  "09_sessionInfo.txt"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

if (
  dir.exists(
    BUILD_QUASI_DIR
  )
) {

  unlink(
    BUILD_QUASI_DIR,
    recursive = TRUE,
    force = TRUE
  )
}

dir.create(
  BUILD_QUASI_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Fixed constants
# ============================================================

YEARS <- 2001:2024

EXPECTED_ROWS <- 769770L

EXPECTED_UNIQUE_MONITOR_DAYS <- 769600L

EXPECTED_MONITORS <- 181L

EXPECTED_EXTRA_QUASI_ROWS <- 170L

EXPECTED_PM_MONITOR_DAYS <- 572251L

PM_AUDIT_TOL <- 1e-10

PRIMARY_EXPOSURE <-
  "campd_upwind_heat_observed_20km_90deg"

PRIMARY_SPATIAL <-
  "campd_upwind_heat_proximity_20km_90deg"


# ============================================================
# 3. Helpers
# ============================================================

safe_num <- function(x) {

  suppressWarnings(
    as.numeric(
      x
    )
  )
}


safe_int <- function(x) {

  suppressWarnings(
    as.integer(
      x
    )
  )
}


first_finite <- function(x) {

  x <- safe_num(x)

  z <- x[
    is.finite(
      x
    )
  ]

  if (
    length(z) == 0L
  ) {
    return(
      NA_real_
    )
  }

  z[1]
}


mean_finite <- function(x) {

  x <- safe_num(x)

  z <- x[
    is.finite(
      x
    )
  ]

  if (
    length(z) == 0L
  ) {
    return(
      NA_real_
    )
  }

  mean(z)
}


strict_mean_finite <- function(x) {

  x <- safe_num(x)

  if (
    length(x) == 0L ||
    any(
      !is.finite(
        x
      )
    )
  ) {

    return(
      NA_real_
    )
  }

  mean(x)
}


first_character <- function(x) {

  z <- as.character(x)

  z <- z[
    !is.na(z) &
      nzchar(z)
  ]

  if (
    length(z) == 0L
  ) {
    return(
      NA_character_
    )
  }

  z[1]
}


first_logical <- function(x) {

  z <- x[
    !is.na(x)
  ]

  if (
    length(z) == 0L
  ) {
    return(
      NA
    )
  }

  as.logical(
    z[1]
  )
}


rename_if_present <- function(
  DT,
  old,
  new
) {

  if (
    old %in%
      names(DT)
  ) {

    if (
      new %in%
        names(DT)
    ) {

      stop(
        "重命名目标已经存在：",
        new
      )
    }

    setnames(
      DT,
      old,
      new
    )
  }
}


fmt_int <- function(x) {

  format(
    x,
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


# ============================================================
# 4. Basic upstream checks
# ============================================================

cat(
  "============================================================\n"
)

cat(
  "07F：正式统计分析数据构建\n"
)

cat(
  "============================================================\n\n"
)

REQUIRED_INPUTS <- c(
  E07_FREEZE,
  E07_HIER,
  PM_DAILY_RETIRED,
  PM_LINK_FILE,
  PM_SITE_DAY_FILE
)

if (
  !all(
    file.exists(
      REQUIRED_INPUTS
    )
  )
) {

  cat(
    "缺失输入：\n"
  )

  print(
    REQUIRED_INPUTS[
      !file.exists(
        REQUIRED_INPUTS
      )
    ]
  )

  stop(
    "07F正式输入不完整。"
  )
}


freeze_text <- readLines(
  E07_FREEZE,
  warn = FALSE,
  encoding = "UTF-8"
)

E07_PASS <- any(
  grepl(
    "07E FINAL FREEZE: PASS",
    freeze_text,
    fixed = TRUE
  )
)

if (
  !E07_PASS
) {

  stop(
    "07E没有FINAL FREEZE PASS，07F禁止继续。"
  )
}

cat(
  "07E FINAL FREEZE：PASS\n"
)


# ============================================================
# 5. Read fixed 07E exposure hierarchy
# ============================================================

HIER <- fread(
  E07_HIER,
  showProgress = FALSE
)

REQ_HIER <- c(
  "variable",
  "source_metric",
  "model_role"
)

MISS_HIER <- setdiff(
  REQ_HIER,
  names(
    HIER
  )
)

if (
  length(
    MISS_HIER
  ) > 0L
) {

  stop(
    "07E模型层级表缺少：",
    paste(
      MISS_HIER,
      collapse = ", "
    )
  )
}


PRIMARY_ROWS <- HIER[
  model_role ==
    "PRIMARY_MAIN"
]

if (
  nrow(
    PRIMARY_ROWS
  ) != 1L
) {

  stop(
    "07E模型层级不是唯一PRIMARY_MAIN。"
  )
}

if (
  PRIMARY_ROWS$variable[1] !=
    PRIMARY_EXPOSURE
) {

  stop(
    "PRIMARY_MAIN与冻结约定不一致。"
  )
}


# ------------------------------------------------------------
# quasi-daily保留所有非QC exposure
# ------------------------------------------------------------

MODEL_EXP <- unique(
  HIER[
    model_role !=
      "QC_ONLY",
    variable
  ]
)

# ------------------------------------------------------------
# filter-period不必复制全部320列。
#
# 保留：
# 1. 主heat
# 2. heat proximity
# 3. 20 km parallel proxies
# 4. 20 km structural / downwind controls
# 5. 各半径heat与plant-count sensitivity
#
# 原始全320列永远保存在冻结07E中。
# ------------------------------------------------------------

FILTER_EXP <- unique(
  HIER[
    model_role %in%
      c(
        "PRIMARY_MAIN",
        "PRIMARY_SPATIAL_SENSITIVITY",
        "PARALLEL_PROXY",
        "PARALLEL_PROXY_SPATIAL_SENSITIVITY",
        "STRUCTURAL_COMPARATOR",
        "DIRECTIONAL_NEGATIVE_CONTROL"
      ) |
      (
        model_role ==
          "RADIUS_SENSITIVITY" &
          source_metric %in%
            c(
              "heat",
              "plant_count"
            )
      ),
    variable
  ]
)


cat(
  "07E non-QC model exposure variables = ",
  length(
    MODEL_EXP
  ),
  "\n",
  sep = ""
)

cat(
  "07F filter-period exposure variables = ",
  length(
    FILTER_EXP
  ),
  "\n\n",
  sep = ""
)


# ============================================================
# 6. Identify 07E yearly files
# ============================================================

YEAR_FILES <- file.path(
  E07_YEAR_DIR,
  paste0(
    "07E_quasi_daily_CAMPD_actual_daily_model_",
    YEARS,
    ".csv.gz"
  )
)

if (
  sum(
    file.exists(
      YEAR_FILES
    )
  ) !=
    24L
) {

  stop(
    "07E yearly文件不是24个。"
  )
}


# ============================================================
# 7. First pass:
#    extract target monitor-day domain
# ============================================================

cat(
  "============================================================\n"
)

cat(
  "A. 建立07E monitor-day目标域\n"
)

cat(
  "============================================================\n"
)


KEY_LIST <- vector(
  "list",
  length(
    YEARS
  )
)

for (
  i in seq_along(
    YEARS
  )
) {

  y <- YEARS[i]

  H <- names(
    fread(
      YEAR_FILES[i],
      nrows = 0
    )
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

    stop(
      y,
      " 07E缺少loc_num/date。"
    )
  }

  K <- fread(
    YEAR_FILES[i],
    select = c(
      "loc_num",
      "date"
    ),
    showProgress = FALSE
  )

  K[
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

  KEY_LIST[[i]] <- unique(
    K[
      ,
      .(
        loc_num,
        date
      )
    ]
  )

  rm(K)

  invisible(
    gc(
      verbose = FALSE
    )
  )
}


TARGET_PM_KEYS <- unique(
  rbindlist(
    KEY_LIST
  )
)

setorder(
  TARGET_PM_KEYS,
  loc_num,
  date
)

N_TARGET_DAYS <- nrow(
  TARGET_PM_KEYS
)

cat(
  "07E unique monitor-days = ",
  fmt_int(
    N_TARGET_DAYS
  ),
  "\n\n",
  sep = ""
)


# ============================================================
# 8. Independent same-day PM reconstruction
# ============================================================

cat(
  "============================================================\n"
)

cat(
  "B. 独立重构same-day PM2.5\n"
)

cat(
  "============================================================\n"
)


# ------------------------------------------------------------
# 8.1 retired daily PM:
# validation reference only
# ------------------------------------------------------------

PM_OFFICIAL <- fread(
  PM_DAILY_RETIRED,
  showProgress = FALSE
)

REQ_PM_OFFICIAL <- c(
  "loc_num",
  "date",
  "pm_c"
)

MISS_PM_OFFICIAL <- setdiff(
  REQ_PM_OFFICIAL,
  names(
    PM_OFFICIAL
  )
)

if (
  length(
    MISS_PM_OFFICIAL
  ) > 0L
) {

  stop(
    "退休same-day PM文件缺少：",
    paste(
      MISS_PM_OFFICIAL,
      collapse = ", "
    )
  )
}


PM_OFFICIAL[
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

    pm_c =
      safe_num(
        pm_c
      )
  )
]


N_PM_OFFICIAL_DUP <- PM_OFFICIAL[
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

if (
  N_PM_OFFICIAL_DUP >
    0L
) {

  stop(
    "退休same-day PM文件存在重复loc_num+date。"
  )
}


# ------------------------------------------------------------
# 8.2 50-km RadNet ↔ PM site link
# ------------------------------------------------------------

PM_LINK <- fread(
  PM_LINK_FILE,
  showProgress = FALSE
)

REQ_LINK <- c(
  "loc_num",
  "pm25_site_id"
)

MISS_LINK <- setdiff(
  REQ_LINK,
  names(
    PM_LINK
  )
)

if (
  length(
    MISS_LINK
  ) > 0L
) {

  stop(
    "PM站点link缺少：",
    paste(
      MISS_LINK,
      collapse = ", "
    )
  )
}


PM_LINK[
  ,
  `:=`(
    loc_num =
      safe_int(
        loc_num
      ),

    pm25_site_id =
      as.character(
        pm25_site_id
      )
  )
]


PM_LINK <- unique(
  PM_LINK[
    !is.na(
      loc_num
    ) &
      !is.na(
        pm25_site_id
      ),
    .(
      loc_num,
      pm25_site_id
    )
  ]
)


# ------------------------------------------------------------
# 8.3 Raw processed PM site-day table
# ------------------------------------------------------------

RAW_PM <- fread(
  PM_SITE_DAY_FILE,
  select = c(
    "pm25_site_id",
    "date_local",
    "pm25_daily_mean",
    "n_pm25_valid"
  ),
  showProgress = FALSE
)


RAW_PM[
  ,
  `:=`(
    pm25_site_id =
      as.character(
        pm25_site_id
      ),

    date_local =
      as.IDate(
        as.character(
          date_local
        )
      ),

    pm25_daily_mean =
      safe_num(
        pm25_daily_mean
      ),

    n_pm25_valid =
      safe_num(
        n_pm25_valid
      )
  )
]


RAW_PM <- RAW_PM[
  is.finite(
    pm25_daily_mean
  ) &
    is.finite(
      n_pm25_valid
    ) &
    n_pm25_valid >
      0
]


LINKED_SITES <- unique(
  PM_LINK$pm25_site_id
)

RAW_PM <- RAW_PM[
  pm25_site_id %in%
    LINKED_SITES
]

setkey(
  RAW_PM,
  pm25_site_id
)


# ------------------------------------------------------------
# 8.4 Recalculate loc_num × date same-day PM
# ------------------------------------------------------------

LOCS <- sort(
  unique(
    PM_LINK$loc_num
  )
)

PM_RECALC_LIST <- vector(
  "list",
  length(
    LOCS
  )
)


for (
  i in seq_along(
    LOCS
  )
) {

  loc <- LOCS[i]

  sites <- PM_LINK[
    loc_num ==
      loc,
    pm25_site_id
  ]

  rr <- RAW_PM[
    J(
      sites
    ),
    nomatch = 0
  ]

  if (
    nrow(rr) ==
      0L
  ) {
    next
  }

  PM_RECALC_LIST[[i]] <- rr[
    ,
    .(
      pm_recalc =
        mean(
          pm25_daily_mean,
          na.rm = TRUE
        ),

      n_pm25_sites_same_day =
        uniqueN(
          pm25_site_id
        )
    ),
    by = .(
      date =
        date_local
    )
  ][
    ,
    loc_num :=
      loc
  ]
}


PM_RECALC <- rbindlist(
  PM_RECALC_LIST,
  use.names = TRUE,
  fill = TRUE
)

setcolorder(
  PM_RECALC,
  c(
    "loc_num",
    "date",
    "pm_recalc",
    "n_pm25_sites_same_day"
  )
)


# ------------------------------------------------------------
# Restrict BOTH tables to exact 07E RadNet monitor-day domain
# ------------------------------------------------------------

PM_RECALC <- PM_RECALC[
  TARGET_PM_KEYS,
  on = .(
    loc_num,
    date
  ),
  nomatch = 0
]


PM_OFFICIAL <- PM_OFFICIAL[
  TARGET_PM_KEYS,
  on = .(
    loc_num,
    date
  ),
  nomatch = 0
]


PM_COMPARE <- merge(
  PM_OFFICIAL[
    ,
    .(
      loc_num,
      date,
      pm_official =
        pm_c
    )
  ],

  PM_RECALC,

  by = c(
    "loc_num",
    "date"
  ),

  all = TRUE,
  sort = FALSE
)


PM_COMPARE[
  ,
  abs_diff :=
    abs(
      pm_official -
        pm_recalc
    )
]


PM_AUDIT <- data.table(

  item = c(
    "official_rows",
    "recalculated_rows",
    "union_rows",
    "missing_official",
    "missing_recalculated",
    "finite_pairs",
    "pairs_abs_diff_gt_tolerance",
    "max_abs_difference",
    "mean_abs_difference",
    "tolerance"
  ),

  value = c(

    nrow(
      PM_OFFICIAL
    ),

    nrow(
      PM_RECALC
    ),

    nrow(
      PM_COMPARE
    ),

    PM_COMPARE[
      !is.finite(
        pm_official
      ) &
        is.finite(
          pm_recalc
        ),
      .N
    ],

    PM_COMPARE[
      is.finite(
        pm_official
      ) &
        !is.finite(
          pm_recalc
        ),
      .N
    ],

    PM_COMPARE[
      is.finite(
        pm_official
      ) &
        is.finite(
          pm_recalc
        ),
      .N
    ],

    PM_COMPARE[
      is.finite(
        abs_diff
      ) &
        abs_diff >
          PM_AUDIT_TOL,
      .N
    ],

    PM_COMPARE[
      is.finite(
        abs_diff
      ),
      max(
        abs_diff
      )
    ],

    PM_COMPARE[
      is.finite(
        abs_diff
      ),
      mean(
        abs_diff
      )
    ],

    PM_AUDIT_TOL
  )
)


fwrite(
  PM_AUDIT,
  OUT_PM_AUDIT
)


PM_MISMATCH <- PM_COMPARE[
  (
    !is.finite(
      pm_official
    ) |
      !is.finite(
        pm_recalc
      )
  ) |
    (
      is.finite(
        abs_diff
      ) &
        abs_diff >
          PM_AUDIT_TOL
    )
]


fwrite(
  PM_MISMATCH,
  OUT_PM_MISMATCH
)


print(
  PM_AUDIT
)


PM_AUDIT_PASS <-
  nrow(
    PM_OFFICIAL
  ) ==
    EXPECTED_PM_MONITOR_DAYS &&
  nrow(
    PM_RECALC
  ) ==
    EXPECTED_PM_MONITOR_DAYS &&
  PM_AUDIT[
    item ==
      "missing_official",
    value
  ] ==
    0 &&
  PM_AUDIT[
    item ==
      "missing_recalculated",
    value
  ] ==
    0 &&
  PM_AUDIT[
    item ==
      "pairs_abs_diff_gt_tolerance",
    value
  ] ==
    0


if (
  !PM_AUDIT_PASS
) {

  unlink(
    BUILD_QUASI_DIR,
    recursive = TRUE,
    force = TRUE
  )

  stop(
    paste0(
      "same-day PM独立重构未通过hard validation。",
      "07F停止，不生成正式分析数据。"
    )
  )
}


cat(
  "\nsame-day PM独立重构：PASS\n\n"
)


# ------------------------------------------------------------
# From now on use ONLY independently reconstructed values
# ------------------------------------------------------------

PM_VERIFIED <- PM_RECALC[
  ,
  .(
    loc_num,
    date,

    pm_same_day =
      pm_recalc,

    n_pm25_sites_same_day
  )
]


PM_MONITOR_COUNT <- PM_VERIFIED[
  is.finite(
    pm_same_day
  ),
  .(
    n_same_day_pm_monitor_days =
      uniqueN(
        date
      )
  ),
  by = loc_num
]


# ============================================================
# 9. Build quasi-daily analysis files
# ============================================================

cat(
  "============================================================\n"
)

cat(
  "C. 构建quasi-daily + same-day PM分析数据\n"
)

cat(
  "============================================================\n"
)


CORE_WANTED <- c(

  "loc_num",
  "radnet_monitor_id",

  "station",
  "city_name",
  "state_abbr",

  "monitor_latitude",
  "monitor_longitude",

  "main_analysis_conus",
  "main_analysis_scope",
  "powerplant_monitor_50km",

  "collect_start",
  "collect_end",
  "date",
  "year",
  "month",
  "doy",

  "gross_beta_bq_m3",
  "gross_beta_mBq_m3",
  "log_gross_beta",

  "surface_u_ppm_article",
  "c_prop",

  "temp",
  "temp_C",

  "hpbl",
  "pblh_km",
  "inv_hpbl",

  "rhum",
  "soilm",
  "vel",

  "pm_c",

  "sunspot_monthly",

  "year_model",
  "year_model_sq",

  "temp_sq",
  "temp_cu",

  "vel_centered_3",
  "pm_centered_10",

  "lat_diff_35",

  "flag_conus",
  "flag_pm_available",
  "flag_monitor_gt500_reference",
  "flag_pm_lt35",
  "flag_vel_lt6",
  "flag_hpbl_150_3000",
  "flag_cprop_soilm_complete",
  "flag_beta_reference_range",
  "flag_core_environment_complete",
  "flag_primary_07c_exposure_available",
  "flag_07c_valid_day",
  "flag_reference_main_sample",

  "campd07d_wind_ok",
  "flag_07d_key_matched"
)


REQUIRED_QUASI <- c(

  "loc_num",
  "radnet_monitor_id",
  "collect_start",
  "collect_end",
  "date",

  "gross_beta_bq_m3",
  "gross_beta_mBq_m3",

  "pm_c",

  "flag_conus",
  "flag_monitor_gt500_reference",
  "flag_vel_lt6",
  "flag_hpbl_150_3000",
  "flag_cprop_soilm_complete",
  "flag_beta_reference_range",

  "campd07d_wind_ok",
  "flag_07d_key_matched",

  PRIMARY_EXPOSURE
)


FP_LIST <- vector(
  "list",
  length(
    YEARS
  )
)


TOTAL_ROWS <- 0L

TOTAL_UNIQUE_DAYS <- 0L

TOTAL_PERIOD_DAY_DUP <- 0L

TOTAL_UNMATCHED_07D <- 0L

TOTAL_PRIMARY_QUASI_SAMPLE <- 0L

ALL_MONITORS <- character()


for (
  i in seq_along(
    YEARS
  )
) {

  y <- YEARS[i]

  cat(
    "07F YEAR: ",
    y,
    "\n",
    sep = ""
  )

  H <- names(
    fread(
      YEAR_FILES[i],
      nrows = 0
    )
  )


  MISS_REQ <- setdiff(
    REQUIRED_QUASI,
    H
  )

  if (
    length(
      MISS_REQ
    ) > 0L
  ) {

    stop(
      y,
      " 07E缺少07F必要字段：",
      paste(
        MISS_REQ,
        collapse = ", "
      )
    )
  }


  MISS_EXP <- setdiff(
    MODEL_EXP,
    H
  )

  if (
    length(
      MISS_EXP
    ) > 0L
  ) {

    stop(
      y,
      " 缺少冻结model exposure字段。"
    )
  }


  SELECT_COLS <- unique(
    c(
      intersect(
        CORE_WANTED,
        H
      ),
      MODEL_EXP
    )
  )


  D <- fread(
    YEAR_FILES[i],
    select = SELECT_COLS,
    showProgress = FALSE
  )


  # ----------------------------------------------------------
  # Make legacy PM meaning explicit
  # ----------------------------------------------------------

  rename_if_present(
    D,
    "pm_c",
    "pm_filter_period_legacy"
  )

  rename_if_present(
    D,
    "pm_centered_10",
    "pm_filter_period_centered_10_legacy"
  )

  rename_if_present(
    D,
    "flag_pm_available",
    "flag_pm_available_legacy_filter_period"
  )

  rename_if_present(
    D,
    "flag_pm_lt35",
    "flag_pm_lt35_legacy_filter_period"
  )

  rename_if_present(
    D,
    "flag_reference_main_sample",
    "flag_reference_main_sample_legacy_filter_period"
  )

  rename_if_present(
    D,
    "flag_primary_07c_exposure_available",
    "flag_primary_07c_exposure_available_legacy"
  )

  rename_if_present(
    D,
    "flag_07c_valid_day",
    "flag_07c_valid_day_legacy"
  )


  D[
    ,
    `:=`(
      loc_num =
        safe_int(
          loc_num
        ),

      radnet_monitor_id =
        as.character(
          radnet_monitor_id
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

      date =
        as.IDate(
          as.character(
            date
          )
        )
    )
  ]


  D[
    ,
    .row_id_07f :=
      .I
  ]


  # ----------------------------------------------------------
  # Verified same-day PM
  # ----------------------------------------------------------

  D <- merge(
    D,
    PM_VERIFIED,
    by = c(
      "loc_num",
      "date"
    ),
    all.x = TRUE,
    sort = FALSE
  )


  D <- merge(
    D,
    PM_MONITOR_COUNT,
    by = "loc_num",
    all.x = TRUE,
    sort = FALSE
  )


  D[
    is.na(
      n_same_day_pm_monitor_days
    ),
    n_same_day_pm_monitor_days :=
      0L
  ]


  D[
    ,
    `:=`(

      pm_same_day_centered_10 =
        pm_same_day -
          10,

      flag_pm_same_day_available =
        is.finite(
          pm_same_day
        ),

      flag_pm_same_day_lt35 =
        is.finite(
          pm_same_day
        ) &
          pm_same_day <
            35,

      flag_primary_heat_daily_available =
        is.finite(
          get(
            PRIMARY_EXPOSURE
          )
        )
    )
  ]


  # ----------------------------------------------------------
  # New quasi-daily paper-style eligibility flag
  #
  # This does NOT fit a model.
  # It only defines the correctly timed same-day-PM sample.
  # ----------------------------------------------------------

  D[
    ,
    flag_07f_primary_quasi_sample :=
      flag_conus %in%
        TRUE &
      flag_pm_same_day_available %in%
        TRUE &
      n_same_day_pm_monitor_days >
        500 &
      flag_pm_same_day_lt35 %in%
        TRUE &
      flag_vel_lt6 %in%
        TRUE &
      flag_hpbl_150_3000 %in%
        TRUE &
      flag_cprop_soilm_complete %in%
        TRUE &
      flag_beta_reference_range %in%
        TRUE &
      campd07d_wind_ok %in%
        TRUE &
      flag_primary_heat_daily_available %in%
        TRUE
  ]


  # ----------------------------------------------------------
  # Restore original row order
  # ----------------------------------------------------------

  setorder(
    D,
    .row_id_07f
  )

  D[
    ,
    .row_id_07f :=
      NULL
  ]


  # ----------------------------------------------------------
  # QC
  # ----------------------------------------------------------

  DUP <- D[
    ,
    .N,
    by = .(
      loc_num,
      collect_start,
      collect_end,
      date
    )
  ][
    N > 1L
  ]


  TOTAL_PERIOD_DAY_DUP <-
    TOTAL_PERIOD_DAY_DUP +
    nrow(
      DUP
    )


  TOTAL_UNMATCHED_07D <-
    TOTAL_UNMATCHED_07D +
    D[
      is.na(
        flag_07d_key_matched
      ) |
        flag_07d_key_matched !=
          TRUE,
      .N
    ]


  TOTAL_PRIMARY_QUASI_SAMPLE <-
    TOTAL_PRIMARY_QUASI_SAMPLE +
    D[
      flag_07f_primary_quasi_sample %in%
        TRUE,
      .N
    ]


  TOTAL_ROWS <-
    TOTAL_ROWS +
    nrow(D)


  TOTAL_UNIQUE_DAYS <-
    TOTAL_UNIQUE_DAYS +
    uniqueN(
      D[
        ,
        .(
          radnet_monitor_id,
          date
        )
      ]
    )


  ALL_MONITORS <- union(
    ALL_MONITORS,
    unique(
      D$radnet_monitor_id
    )
  )


  # ----------------------------------------------------------
  # Save formal quasi-daily building output
  # ----------------------------------------------------------

  out_y <- file.path(
    BUILD_QUASI_DIR,
    paste0(
      "07F_quasi_daily_analysis_base_",
      y,
      ".csv.gz"
    )
  )


  fwrite(
    D,
    out_y
  )


  # ----------------------------------------------------------
  # Keep only fields needed for filter-period construction
  # ----------------------------------------------------------

  FP_BASE_WANTED <- c(

    "loc_num",
    "radnet_monitor_id",

    "station",
    "city_name",
    "state_abbr",

    "monitor_latitude",
    "monitor_longitude",

    "main_analysis_conus",
    "flag_conus",

    "collect_start",
    "collect_end",
    "date",

    "gross_beta_bq_m3",
    "gross_beta_mBq_m3",
    "log_gross_beta",

    "surface_u_ppm_article",
    "c_prop",

    "temp",
    "temp_C",

    "hpbl",
    "pblh_km",
    "inv_hpbl",

    "rhum",
    "soilm",
    "vel",

    "sunspot_monthly",
    "year_model",
    "lat_diff_35",

    "pm_filter_period_legacy",
    "pm_same_day",
    "n_pm25_sites_same_day",

    "campd07d_wind_ok"
  )


  FP_KEEP <- unique(
    c(
      intersect(
        FP_BASE_WANTED,
        names(D)
      ),
      FILTER_EXP
    )
  )


  FP_LIST[[i]] <- D[
    ,
    ..FP_KEEP
  ]


  cat(
    "  rows = ",
    fmt_int(
      nrow(D)
    ),
    "\n",
    sep = ""
  )

  cat(
    "  primary quasi sample = ",
    fmt_int(
      D[
        flag_07f_primary_quasi_sample %in%
          TRUE,
        .N
      ]
    ),
    "\n",
    sep = ""
  )


  rm(
    D,
    DUP
  )

  invisible(
    gc(
      verbose = FALSE
    )
  )
}


# ============================================================
# 10. Build one-record-per-filter-period dataset
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "D. 构建filter-period单记录分析数据\n"
)

cat(
  "============================================================\n"
)


Q <- rbindlist(
  FP_LIST,
  use.names = TRUE,
  fill = TRUE
)


PERIOD_KEYS <- c(
  "loc_num",
  "radnet_monitor_id",
  "collect_start",
  "collect_end"
)


# ------------------------------------------------------------
# 10.1 Period constancy audits
# ------------------------------------------------------------

PERIOD_CONSTANCY <- Q[
  ,
  .(

    n_beta_bq =
      uniqueN(
        gross_beta_bq_m3[
          is.finite(
            gross_beta_bq_m3
          )
        ]
      ),

    n_beta_mbq =
      uniqueN(
        gross_beta_mBq_m3[
          is.finite(
            gross_beta_mBq_m3
          )
        ]
      ),

    n_pm_filter_period =
      uniqueN(
        pm_filter_period_legacy[
          is.finite(
            pm_filter_period_legacy
          )
        ]
      )
  ),

  by = PERIOD_KEYS
]


BAD_BETA_PERIODS <- PERIOD_CONSTANCY[
  n_beta_bq >
    1L |
    n_beta_mbq >
      1L,
  .N
]


BAD_PM_PERIODS <- PERIOD_CONSTANCY[
  n_pm_filter_period >
    1L,
  .N
]


# ------------------------------------------------------------
# 10.2 Base period table
# ------------------------------------------------------------

FP_BASE <- Q[
  ,
  .(

    station =
      if (
        "station" %in%
          names(Q)
      ) {
        first_character(
          station
        )
      } else {
        NA_character_
      },

    city_name =
      if (
        "city_name" %in%
          names(Q)
      ) {
        first_character(
          city_name
        )
      } else {
        NA_character_
      },

    state_abbr =
      if (
        "state_abbr" %in%
          names(Q)
      ) {
        first_character(
          state_abbr
        )
      } else {
        NA_character_
      },

    monitor_latitude =
      if (
        "monitor_latitude" %in%
          names(Q)
      ) {
        first_finite(
          monitor_latitude
        )
      } else {
        NA_real_
      },

    monitor_longitude =
      if (
        "monitor_longitude" %in%
          names(Q)
      ) {
        first_finite(
          monitor_longitude
        )
      } else {
        NA_real_
      },

    main_analysis_conus =
      if (
        "main_analysis_conus" %in%
          names(Q)
      ) {
        first_logical(
          main_analysis_conus
        )
      } else {
        NA
      },

    flag_conus =
      first_logical(
        flag_conus
      ),

    gross_beta_bq_m3 =
      first_finite(
        gross_beta_bq_m3
      ),

    gross_beta_mBq_m3 =
      first_finite(
        gross_beta_mBq_m3
      ),

    log_gross_beta =
      if (
        "log_gross_beta" %in%
          names(Q)
      ) {
        first_finite(
          log_gross_beta
        )
      } else {
        NA_real_
      },

    # station/static or slowly varying field
    surface_u_ppm_article =
      if (
        "surface_u_ppm_article" %in%
          names(Q)
      ) {
        first_finite(
          surface_u_ppm_article
        )
      } else {
        NA_real_
      },

    # Daily/period environmental variables:
    # average across the collection period
    c_prop =
      if (
        "c_prop" %in%
          names(Q)
      ) {
        mean_finite(
          c_prop
        )
      } else {
        NA_real_
      },

    temp =
      if (
        "temp" %in%
          names(Q)
      ) {
        mean_finite(
          temp
        )
      } else {
        NA_real_
      },

    temp_C =
      if (
        "temp_C" %in%
          names(Q)
      ) {
        mean_finite(
          temp_C
        )
      } else {
        NA_real_
      },

    hpbl =
      if (
        "hpbl" %in%
          names(Q)
      ) {
        mean_finite(
          hpbl
        )
      } else {
        NA_real_
      },

    pblh_km =
      if (
        "pblh_km" %in%
          names(Q)
      ) {
        mean_finite(
          pblh_km
        )
      } else {
        NA_real_
      },

    inv_hpbl =
      if (
        "inv_hpbl" %in%
          names(Q)
      ) {
        mean_finite(
          inv_hpbl
        )
      } else {
        NA_real_
      },

    rhum =
      if (
        "rhum" %in%
          names(Q)
      ) {
        mean_finite(
          rhum
        )
      } else {
        NA_real_
      },

    soilm =
      if (
        "soilm" %in%
          names(Q)
      ) {
        mean_finite(
          soilm
        )
      } else {
        NA_real_
      },

    vel =
      if (
        "vel" %in%
          names(Q)
      ) {
        mean_finite(
          vel
        )
      } else {
        NA_real_
      },

    sunspot_monthly =
      if (
        "sunspot_monthly" %in%
          names(Q)
      ) {
        mean_finite(
          sunspot_monthly
        )
      } else {
        NA_real_
      },

    year_model =
      if (
        "year_model" %in%
          names(Q)
      ) {
        mean_finite(
          year_model
        )
      } else {
        NA_real_
      },

    lat_diff_35 =
      if (
        "lat_diff_35" %in%
          names(Q)
      ) {
        first_finite(
          lat_diff_35
        )
      } else {
        NA_real_
      },

    # This PM variable is correctly defined at filter-period level.
    pm_filter_period_official =
      first_finite(
        pm_filter_period_legacy
      ),

    # Diagnostic / sensitivity only:
    # arithmetic mean of independently reconstructed daily PM
    pm_same_day_mean_available =
      mean_finite(
        pm_same_day
      ),

    n_pm_same_day_days =
      sum(
        is.finite(
          pm_same_day
        )
      ),

    n_quasi_days =
      .N,

    n_valid_wind_days =
      sum(
        campd07d_wind_ok %in%
          TRUE
      ),

    first_date =
      min(
        date
      ),

    last_date =
      max(
        date
      )
  ),

  by = PERIOD_KEYS
]


FP_BASE[
  ,
  `:=`(

    n_calendar_days_expected =
      as.integer(
        collect_end -
          collect_start
      ),

    pm_filter_period_centered_10 =
      pm_filter_period_official -
        10,

    pm_same_day_period_mean_centered_10 =
      pm_same_day_mean_available -
        10
  )
]


FP_BASE[
  ,
  pm_same_day_period_coverage :=
    fifelse(
      n_quasi_days >
        0L,

      n_pm_same_day_days /
        n_quasi_days,

      NA_real_
    )
]


# ------------------------------------------------------------
# Recompute transformations from period-averaged values
# ------------------------------------------------------------

FP_BASE[
  ,
  `:=`(

    temp_sq =
      temp^2,

    temp_cu =
      temp^3,

    year_model_sq =
      year_model^2,

    vel_centered_3 =
      vel -
        3
  )
]


# ------------------------------------------------------------
# 10.3 Strict-complete filter-period exposure averages
#
# Main principle:
# gross-beta is a period-average concentration,
# therefore use average DAILY source activity over the period,
# not total accumulated heat input.
#
# If one or more daily values are NA, the strict-complete
# period exposure remains NA.
# ------------------------------------------------------------

FP_EXP <- Q[
  ,
  lapply(
    .SD,
    strict_mean_finite
  ),
  by = PERIOD_KEYS,
  .SDcols = FILTER_EXP
]


setnames(
  FP_EXP,
  FILTER_EXP,
  paste0(
    "fpmean_",
    FILTER_EXP
  )
)


FP <- merge(
  FP_BASE,
  FP_EXP,
  by = PERIOD_KEYS,
  all = TRUE,
  sort = FALSE
)


PRIMARY_FP <-
  paste0(
    "fpmean_",
    PRIMARY_EXPOSURE
  )


FP[
  ,
  `:=`(

    period_year =
      as.integer(
        format(
          collect_start,
          "%Y"
        )
      ),

    period_month =
      as.integer(
        format(
          collect_start,
          "%m"
        )
      ),

    flag_primary_heat_filter_period_complete =
      is.finite(
        get(
          PRIMARY_FP
        )
      ),

    flag_pm_filter_period_available =
      is.finite(
        pm_filter_period_official
      ),

    flag_beta_available =
      is.finite(
        gross_beta_bq_m3
      )
  )
]


FP[
  ,
  flag_07f_filter_period_primary_base :=
    flag_conus %in%
      TRUE &
    flag_beta_available %in%
      TRUE &
    flag_primary_heat_filter_period_complete %in%
      TRUE
]


setorder(
  FP,
  loc_num,
  collect_start,
  collect_end
)


# ============================================================
# 11. Filter-period QC
# ============================================================

N_FILTER_PERIODS <- nrow(
  FP
)


N_FILTER_DUP <- FP[
  ,
  .N,
  by = PERIOD_KEYS
][
  N > 1L,
  .N
]


N_DATE_EXPANSION_MISMATCH <- FP[
  n_quasi_days !=
    n_calendar_days_expected |
    first_date !=
      collect_start |
    last_date !=
      (
        collect_end -
          1
      ),
  .N
]


N_PRIMARY_FP_COMPLETE <- FP[
  flag_primary_heat_filter_period_complete %in%
    TRUE,
  .N
]


# ============================================================
# 12. Filter-period exposure coverage
# ============================================================

FP_COV <- rbindlist(
  lapply(
    FILTER_EXP,

    function(v) {

      vv <- paste0(
        "fpmean_",
        v
      )

      x <- safe_num(
        FP[[vv]]
      )

      data.table(

        variable =
          v,

        filter_period_variable =
          vv,

        n_filter_periods =
          length(x),

        n_complete =
          sum(
            is.finite(
              x
            )
          ),

        pct_complete =
          100 *
          sum(
            is.finite(
              x
            )
          ) /
          length(x),

        n_zero =
          sum(
            is.finite(
              x
            ) &
              x ==
                0
          ),

        n_positive =
          sum(
            is.finite(
              x
            ) &
              x >
                0
          ),

        n_negative =
          sum(
            is.finite(
              x
            ) &
              x <
                0
          )
      )
    }
  )
)


FP_COV <- merge(
  FP_COV,

  HIER[
    ,
    .(
      variable,
      source_metric,
      representation,
      direction,
      radius_km,
      model_role
    )
  ],

  by = "variable",
  all.x = TRUE,
  sort = FALSE
)


setorder(
  FP_COV,
  model_role,
  source_metric,
  direction,
  radius_km,
  variable
)


# ============================================================
# 13. Final hard checks
# ============================================================

HARD <- data.table(
  check = character(),
  observed = character(),
  expected = character(),
  status = character()
)


add_hard <- function(
  check,
  condition,
  observed,
  expected
) {

  HARD <<- rbind(
    HARD,

    data.table(
      check =
        check,

      observed =
        as.character(
          observed
        ),

      expected =
        as.character(
          expected
        ),

      status =
        ifelse(
          isTRUE(
            condition
          ),
          "PASS",
          "FAIL"
        )
    )
  )
}


add_hard(
  "07E_final_freeze_PASS",
  E07_PASS,
  E07_PASS,
  TRUE
)


add_hard(
  "07E_yearly_files_24",
  sum(
    file.exists(
      YEAR_FILES
    )
  ) ==
    24L,
  sum(
    file.exists(
      YEAR_FILES
    )
  ),
  24
)


add_hard(
  "07E_unique_monitor_day_domain_769600",
  N_TARGET_DAYS ==
    EXPECTED_UNIQUE_MONITOR_DAYS,
  N_TARGET_DAYS,
  EXPECTED_UNIQUE_MONITOR_DAYS
)


add_hard(
  "same_day_PM_official_rows_572251",
  nrow(
    PM_OFFICIAL
  ) ==
    EXPECTED_PM_MONITOR_DAYS,
  nrow(
    PM_OFFICIAL
  ),
  EXPECTED_PM_MONITOR_DAYS
)


add_hard(
  "same_day_PM_recalculated_rows_572251",
  nrow(
    PM_RECALC
  ) ==
    EXPECTED_PM_MONITOR_DAYS,
  nrow(
    PM_RECALC
  ),
  EXPECTED_PM_MONITOR_DAYS
)


add_hard(
  "same_day_PM_missing_official_zero",
  PM_AUDIT[
    item ==
      "missing_official",
    value
  ] ==
    0,
  PM_AUDIT[
    item ==
      "missing_official",
    value
  ],
  0
)


add_hard(
  "same_day_PM_missing_recalculated_zero",
  PM_AUDIT[
    item ==
      "missing_recalculated",
    value
  ] ==
    0,
  PM_AUDIT[
    item ==
      "missing_recalculated",
    value
  ],
  0
)


add_hard(
  "same_day_PM_pairs_over_tolerance_zero",
  PM_AUDIT[
    item ==
      "pairs_abs_diff_gt_tolerance",
    value
  ] ==
    0,
  PM_AUDIT[
    item ==
      "pairs_abs_diff_gt_tolerance",
    value
  ],
  0
)


add_hard(
  "quasi_daily_rows_769770",
  TOTAL_ROWS ==
    EXPECTED_ROWS,
  TOTAL_ROWS,
  EXPECTED_ROWS
)


add_hard(
  "quasi_daily_unique_monitor_days_769600",
  TOTAL_UNIQUE_DAYS ==
    EXPECTED_UNIQUE_MONITOR_DAYS,
  TOTAL_UNIQUE_DAYS,
  EXPECTED_UNIQUE_MONITOR_DAYS
)


add_hard(
  "quasi_daily_extra_rows_170",
  TOTAL_ROWS -
    TOTAL_UNIQUE_DAYS ==
    EXPECTED_EXTRA_QUASI_ROWS,
  TOTAL_ROWS -
    TOTAL_UNIQUE_DAYS,
  EXPECTED_EXTRA_QUASI_ROWS
)


add_hard(
  "quasi_daily_monitors_181",
  length(
    ALL_MONITORS
  ) ==
    EXPECTED_MONITORS,
  length(
    ALL_MONITORS
  ),
  EXPECTED_MONITORS
)


add_hard(
  "quasi_period_day_duplicate_keys_zero",
  TOTAL_PERIOD_DAY_DUP ==
    0L,
  TOTAL_PERIOD_DAY_DUP,
  0
)


add_hard(
  "quasi_unmatched_07D_rows_zero",
  TOTAL_UNMATCHED_07D ==
    0L,
  TOTAL_UNMATCHED_07D,
  0
)


add_hard(
  "primary_quasi_sample_has_rows",
  TOTAL_PRIMARY_QUASI_SAMPLE >
    0L,
  TOTAL_PRIMARY_QUASI_SAMPLE,
  ">0"
)


add_hard(
  "filter_period_duplicate_keys_zero",
  N_FILTER_DUP ==
    0L,
  N_FILTER_DUP,
  0
)


add_hard(
  "filter_period_gross_beta_constant",
  BAD_BETA_PERIODS ==
    0L,
  BAD_BETA_PERIODS,
  0
)


add_hard(
  "filter_period_legacy_PM_constant",
  BAD_PM_PERIODS ==
    0L,
  BAD_PM_PERIODS,
  0
)


add_hard(
  "filter_period_daily_expansion_complete",
  N_DATE_EXPANSION_MISMATCH ==
    0L,
  N_DATE_EXPANSION_MISMATCH,
  0
)


add_hard(
  "filter_period_primary_heat_complete_has_rows",
  N_PRIMARY_FP_COMPLETE >
    0L,
  N_PRIMARY_FP_COMPLETE,
  ">0"
)


# ============================================================
# 14. Warnings / interpretation guards
# ============================================================

WARN <- data.table(

  item = c(

    "legacy_filter_period_PM_in_quasi_output",

    "quasi_daily_overlap_rows",

    "retired_daily_PM_role",

    "filter_period_exposure_definition",

    "filter_period_partial_daily_exposure",

    "no_model_fitted_in_07F"
  ),

  status = c(
    "WARN",
    "WARN",
    "INFO",
    "INFO",
    "INFO",
    "INFO"
  ),

  note = c(

    paste0(
      "Old pm_c is retained only as ",
      "pm_filter_period_legacy. ",
      "Do not use it as same-day PM in quasi-daily models."
    ),

    paste0(
      "The 170 quasi-daily rows above unique monitor-days are ",
      "intentional overlapping filter-period dates, not duplicate errors."
    ),

    paste0(
      "The retired daily PM file is used only as an independent ",
      "validation target. Formal 07F same-day PM values come from ",
      "the independent site-link + site-day reconstruction."
    ),

    paste0(
      "Filter-period CAMPD exposure is the arithmetic mean of ",
      "daily source activity over [collect_start, collect_end). ",
      "This matches a period-average gross-beta concentration ",
      "better than summing source activity across periods of ",
      "different duration."
    ),

    paste0(
      "Formal fpmean exposure is strict-complete: ",
      "if any daily exposure value in a filter period is NA, ",
      "the period exposure is NA rather than averaging only ",
      "the observed subset of days."
    ),

    "07F constructs analysis data only. No regression model is fitted."
  )
)


# ============================================================
# 15. QC summary
# ============================================================

QC <- data.table(

  item = c(

    "quasi_daily_rows",

    "quasi_daily_unique_monitor_days",

    "quasi_daily_monitors",

    "quasi_daily_overlap_rows",

    "same_day_PM_unique_monitor_days",

    "primary_quasi_sample_rows",

    "filter_periods",

    "filter_period_primary_heat_complete",

    "filter_period_date_expansion_mismatch",

    "filter_period_beta_nonconstant",

    "filter_period_PM_nonconstant"
  ),

  value = c(

    TOTAL_ROWS,

    TOTAL_UNIQUE_DAYS,

    length(
      ALL_MONITORS
    ),

    TOTAL_ROWS -
      TOTAL_UNIQUE_DAYS,

    nrow(
      PM_VERIFIED
    ),

    TOTAL_PRIMARY_QUASI_SAMPLE,

    N_FILTER_PERIODS,

    N_PRIMARY_FP_COMPLETE,

    N_DATE_EXPANSION_MISMATCH,

    BAD_BETA_PERIODS,

    BAD_PM_PERIODS
  )
)


# ============================================================
# 16. Write QC before promotion
# ============================================================

fwrite(
  FP_COV,
  OUT_FILTER_COV
)

fwrite(
  QC,
  OUT_QC
)

fwrite(
  HARD,
  OUT_HARD
)

fwrite(
  WARN,
  OUT_WARN
)


cat(
  "\n============================================================\n"
)

cat(
  "07F HARD CHECKS\n"
)

cat(
  "============================================================\n"
)

print(
  HARD
)


N_HARD_FAIL <- HARD[
  status ==
    "FAIL",
  .N
]

cat(
  "\nHard FAIL = ",
  N_HARD_FAIL,
  "\n",
  sep = ""
)


# ============================================================
# 17. Hard failure:
#     do not replace formal quasi outputs
# ============================================================

if (
  N_HARD_FAIL >
    0L
) {

  unlink(
    BUILD_QUASI_DIR,
    recursive = TRUE,
    force = TRUE
  )

  cat(
    "\n============================================================\n"
  )

  cat(
    "07F TECHNICAL STATUS: FAIL\n"
  )

  cat(
    "============================================================\n"
  )

  stop(
    "07F存在hard failure；未替换正式quasi-daily输出。"
  )
}


# ============================================================
# 18. Promote quasi-daily outputs
# ============================================================

if (
  dir.exists(
    FINAL_QUASI_DIR
  )
) {

  unlink(
    FINAL_QUASI_DIR,
    recursive = TRUE,
    force = TRUE
  )
}


PROMOTED <- file.rename(
  BUILD_QUASI_DIR,
  FINAL_QUASI_DIR
)


if (
  !PROMOTED
) {

  stop(
    "07F quasi-daily building目录无法提升为正式输出。"
  )
}


# ============================================================
# 19. Write filter-period formal output
# ============================================================

fwrite(
  FP,
  OUT_FILTER
)


# ============================================================
# 20. Report
# ============================================================

PRIMARY_FP_ROW <- FP_COV[
  variable ==
    PRIMARY_EXPOSURE
]


REPORT <- c(

  "# 07F 正式统计分析数据构建报告",

  "",

  paste0(
    "- Generated: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),

  "- Upstream 07E FINAL FREEZE: PASS",

  paste0(
    "- Hard failures: ",
    N_HARD_FAIL
  ),

  "",

  "## Purpose",

  "",

  paste0(
    "07F does not modify frozen 07D/07E and does not fit ",
    "statistical models."
  ),

  paste0(
    "It creates two downstream datasets with explicit temporal ",
    "definitions."
  ),

  "",

  "## Branch A — quasi-daily",

  "",

  paste0(
    "- Rows: ",
    fmt_int(
      TOTAL_ROWS
    )
  ),

  paste0(
    "- Unique monitor-days: ",
    fmt_int(
      TOTAL_UNIQUE_DAYS
    )
  ),

  paste0(
    "- Monitors: ",
    length(
      ALL_MONITORS
    )
  ),

  paste0(
    "- Intentional overlap rows: ",
    TOTAL_ROWS -
      TOTAL_UNIQUE_DAYS
  ),

  paste0(
    "- Same-day PM monitor-days independently reconstructed: ",
    fmt_int(
      nrow(
        PM_VERIFIED
      )
    )
  ),

  paste0(
    "- Primary quasi-daily eligible rows: ",
    fmt_int(
      TOTAL_PRIMARY_QUASI_SAMPLE
    )
  ),

  "",

  paste0(
    "The inherited pm_c variable has been explicitly renamed ",
    "pm_filter_period_legacy and must not be used as the ",
    "same-day PM covariate."
  ),

  "",

  "## Branch B — filter-period",

  "",

  paste0(
    "- Unique filter periods: ",
    fmt_int(
      N_FILTER_PERIODS
    )
  ),

  paste0(
    "- Primary heat exposure with complete daily coverage: ",
    fmt_int(
      N_PRIMARY_FP_COMPLETE
    )
  ),

  if (
    nrow(
      PRIMARY_FP_ROW
    ) ==
      1L
  ) {

    paste0(
      "- Primary heat complete-period coverage: ",
      sprintf(
        "%.3f%%",
        PRIMARY_FP_ROW$pct_complete
      )
    )

  } else {

    "- Primary heat complete-period coverage: unavailable"
  },

  "",

  paste0(
    "Filter-period exposure is the arithmetic mean of daily ",
    "source activity over [collect_start, collect_end)."
  ),

  paste0(
    "The strict-complete definition returns NA if any day in ",
    "the period lacks that source metric."
  ),

  "",

  "## Fixed primary exposure",

  "",

  paste0(
    "- Quasi-daily primary: ",
    PRIMARY_EXPOSURE
  ),

  paste0(
    "- Filter-period primary: ",
    PRIMARY_FP
  ),

  "",

  "Heat input is interpreted as the primary combustion-activity proxy.",

  paste0(
    "NOx/SO2/CO2 are parallel conventional-emission proxies; ",
    "gross load is a generation-activity proxy."
  ),

  "None is a direct measurement of radioactive emissions.",

  "",

  "## same-day PM validation",

  "",

  paste0(
    "- Expected monitor-days: ",
    fmt_int(
      EXPECTED_PM_MONITOR_DAYS
    )
  ),

  paste0(
    "- Reconstructed monitor-days: ",
    fmt_int(
      nrow(
        PM_RECALC
      )
    )
  ),

  paste0(
    "- Pairs above tolerance ",
    PM_AUDIT_TOL,
    ": ",
    PM_AUDIT[
      item ==
        "pairs_abs_diff_gt_tolerance",
      value
    ]
  ),

  "",

  "## Status",

  "",

  "**07F TECHNICAL STATUS: PASS**",

  "",

  paste0(
    "Do not freeze 07F until the resulting quasi-daily and ",
    "filter-period sample sizes/coverage are reviewed."
  )
)


writeLines(
  REPORT,
  OUT_REPORT,
  useBytes = TRUE
)


writeLines(
  capture.output(
    sessionInfo()
  ),
  OUT_SESSION,
  useBytes = TRUE
)


# ============================================================
# 21. Final console output
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "07F TECHNICAL STATUS: PASS\n"
)

cat(
  "============================================================\n"
)


cat(
  "Quasi-daily rows = ",
  fmt_int(
    TOTAL_ROWS
  ),
  "\n",
  sep = ""
)


cat(
  "Unique monitor-days = ",
  fmt_int(
    TOTAL_UNIQUE_DAYS
  ),
  "\n",
  sep = ""
)


cat(
  "RadNet monitors = ",
  length(
    ALL_MONITORS
  ),
  "\n",
  sep = ""
)


cat(
  "Verified same-day PM monitor-days = ",
  fmt_int(
    nrow(
      PM_VERIFIED
    )
  ),
  "\n",
  sep = ""
)


cat(
  "Primary quasi-daily sample rows = ",
  fmt_int(
    TOTAL_PRIMARY_QUASI_SAMPLE
  ),
  "\n",
  sep = ""
)


cat(
  "Filter periods = ",
  fmt_int(
    N_FILTER_PERIODS
  ),
  "\n",
  sep = ""
)


cat(
  "Primary heat complete filter periods = ",
  fmt_int(
    N_PRIMARY_FP_COMPLETE
  ),
  "\n",
  sep = ""
)


cat(
  "Filter-period date expansion mismatches = ",
  N_DATE_EXPANSION_MISMATCH,
  "\n",
  sep = ""
)


cat(
  "Hard FAIL = ",
  N_HARD_FAIL,
  "\n",
  sep = ""
)


cat(
  "\n正式quasi-daily输出：\n",
  FINAL_QUASI_DIR,
  "\n",
  sep = ""
)


cat(
  "\n正式filter-period输出：\n",
  OUT_FILTER,
  "\n",
  sep = ""
)


cat(
  "\nQC：\n",
  OUT_PM_AUDIT,
  "\n",
  OUT_PM_MISMATCH,
  "\n",
  OUT_FILTER_COV,
  "\n",
  OUT_QC,
  "\n",
  OUT_HARD,
  "\n",
  OUT_WARN,
  "\n",
  OUT_REPORT,
  "\n",
  sep = ""
)

