#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999,
  width = 220
)

# ============================================================
# 07H13
# 非方向燃料特异实际运行分析 + Figure 4 / Figure 5
#
# Figure 4:
#   不考虑风向
#   圆形缓冲区内燃料特异实际逐日运行
#
# Figure 5:
#   已冻结H10的90°上风向燃料特异结果
#
# 不修改：
#   07A / 07B2 / 07F / 07H01 / 07H10 / 07H11
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

PKGS <- c(
  "data.table",
  "fixest",
  "ggplot2",
  "patchwork"
)

for (p in PKGS) {
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

OUT <- file.path(
  HROOT,
  "07H13_非方向燃料特异实际运行分析"
)

FIG_DIR <- file.path(
  HROOT,
  "07H12_论文正式图表总整理",
  "01_正文图"
)

SOURCE_DIR <- file.path(
  HROOT,
  "07H12_论文正式图表总整理",
  "05_绘图数据"
)

AUDIT_DIR <- file.path(
  OUT,
  "审计"
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
  SOURCE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  AUDIT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ------------------------------------------------------------
# 冻结07A：总体空间连接
# ------------------------------------------------------------

GEO_FILE <- file.path(
  PP,
  "00_既往发电站空间暴露与模型流程_冻结",
  "02_07A_RadNet发电站空间链接",
  "07A_radnet_powerplant_same_year_250km",
  "03_monitor_year_all_plant_link_250km_07A_2001_2024.csv"
)


# ------------------------------------------------------------
# 冻结07B2：CAMPD plant-day × primary fuel
# 当前07H01实际使用的正式目录
# ------------------------------------------------------------

PF_DIR <- file.path(
  PP,
  "01_发电站月度逐日源强增强流程_冻结",
  "07B2_CAMPD_daily_emissions",
  "05_plant_day_primary_fuel"
)


# ------------------------------------------------------------
# 冻结07B3：用于限制正式preferred CAMPD plant-month
# ------------------------------------------------------------

HARM_FILE <- file.path(
  PP,
  "01_发电站月度逐日源强增强流程_冻结",
  "07B3_EIA_CAMPD_harmonization",
  "05_EIA_CAMPD_plant_month_harmonization_2001_2024.csv.gz"
)


# ------------------------------------------------------------
# 冻结07F filter-period正式分析数据
# ------------------------------------------------------------

FILTER_FILE <- file.path(
  PP,
  "04_07F_正式统计分析数据构建",
  "03_07F_filter_period_analysis_base_2001_2024.csv.gz"
)


# ------------------------------------------------------------
# 冻结07H10：
# Figure 5直接读取，不重新估计
# ------------------------------------------------------------

H10_EFFECT_FILE <- file.path(
  HROOT,
  "07H10_燃料特异站点年份内综合正式分析",
  "01_主分析_监测站年份固定效应",
  "02_全部半径主模型代表风速效应.csv"
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

CORE_FUELS <- c(
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
  50L,
  100L,
  150L,
  200L
)

MAX_RADIUS <- max(
  RADII
)

HEAT_INCREMENT_MMBTU <- 10000

EXPECTED_FILTER_ROWS <- 107407L
EXPECTED_FILTER_MONITORS <- 133L


FUEL_LABEL_EN <- c(
  biomass = "Biomass",
  coal = "Coal",
  gas = "Natural gas",
  oil = "Oil",
  other_fossil = "Other fossil",
  solid_waste = "Solid waste"
)

FUEL_LABEL_CN <- c(
  biomass = "生物质",
  coal = "煤",
  gas = "天然气",
  oil = "燃油",
  other_fossil = "其他化石燃料",
  solid_waste = "固体废物"
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

  y <- sub(
    "\\.0+$",
    "",
    y
  )

  y
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
    length(hit) == 0L
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
    length(x) == 0L
  ) {
    return(
      NA_real_
    )
  }

  if (
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
      expected_n
  ) {
    return(
      NA_real_
    )
  }

  if (
    any(
      !is.finite(x)
    )
  ) {
    return(
      NA_real_
    )
  }

  mean(
    x
  )
}


collapse_messages <- function(x) {

  if (
    is.null(x) ||
    length(x) == 0L
  ) {
    return("")
  }

  x <- unique(
    trimws(
      as.character(x)
    )
  )

  x <- x[
    nzchar(x)
  ]

  paste(
    x,
    collapse = " | "
  )
}


# ============================================================
# 4. Input audit
# ============================================================

cat(
  "\n============================================================\n",
  "07H13 INPUT AUDIT\n",
  "============================================================\n",
  sep = ""
)

for (
  f in c(
    GEO_FILE,
    HARM_FILE,
    FILTER_FILE,
    H10_EFFECT_FILE
  )
) {

  if (!file.exists(f)) {
    stopf(
      "缺少正式输入：\n",
      f
    )
  }
}

if (!dir.exists(PF_DIR)) {
  stopf(
    "缺少正式CAMPD primary-fuel目录：\n",
    PF_DIR
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
    "正式CAMPD primary-fuel季度文件不是96/96。"
  )
}

cat(
  "07A geometry：PASS\n",
  "07B2 primary-fuel：96/96 PASS\n",
  "07B3 harmonization：PASS\n",
  "07F filter-period：PASS\n",
  "07H10 existing upwind effects：PASS\n",
  sep = ""
)


# ============================================================
# 5. Read and standardize 07A geometry
#
# 冻结07A总体geometry不强制要求直接存在loc_num。
# 与冻结07H01保持完全一致：
# 若07A无loc_num，则通过07F的
# radnet_monitor_id <-> loc_num一一对应表恢复。
# ============================================================

cat(
  "\n============================================================\n",
  "读取冻结07A geometry\n",
  "============================================================\n",
  sep = ""
)


# ------------------------------------------------------------
# 5A. 从冻结07F逐年准逐日文件建立monitor ID crosswalk
# ------------------------------------------------------------

F07_QDIR <- file.path(
  PP,
  "04_07F_正式统计分析数据构建",
  "quasi_daily_yearly"
)

if (!dir.exists(F07_QDIR)) {
  stopf(
    "找不到冻结07F quasi_daily_yearly目录：\n",
    F07_QDIR
  )
}


cross_parts <- vector(
  "list",
  length(YEARS)
)


for (ii in seq_along(YEARS)) {

  yy <- YEARS[ii]

  qf <- file.path(
    F07_QDIR,
    paste0(
      "07F_quasi_daily_analysis_base_",
      yy,
      ".csv.gz"
    )
  )

  if (!file.exists(qf)) {
    stopf(
      "缺少冻结07F逐年文件：",
      qf
    )
  }

  qh <- names(
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
    qh
  )

  if (length(miss) > 0L) {
    stopf(
      yy,
      "年07F缺少monitor crosswalk字段：",
      paste(miss, collapse = ", ")
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

  cross_parts[[ii]] <- unique(
    z[
      ,
      .(
        loc_num = safe_int(loc_num),
        radnet_monitor_id = as.character(radnet_monitor_id)
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
  !is.na(radnet_monitor_id)
]


dup_loc <- CROSS[
  ,
  .N,
  by = loc_num
][
  N > 1L
]

dup_id <- CROSS[
  ,
  .N,
  by = radnet_monitor_id
][
  N > 1L
]


if (nrow(dup_loc) > 0L || nrow(dup_id) > 0L) {
  stopf(
    "07F loc_num <-> radnet_monitor_id不是一一对应。"
  )
}


cat(
  "07F monitor crosswalk = ",
  format(nrow(CROSS), big.mark = ","),
  " unique monitors\n",
  sep = ""
)


# ------------------------------------------------------------
# 5B. 读取07A；loc_num允许不存在
# ------------------------------------------------------------

GH <- names(
  fread(
    GEO_FILE,
    nrows = 0,
    showProgress = FALSE
  )
)

cat(
  "07A原始字段数 = ",
  length(GH),
  "\n",
  sep = ""
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


if (year_col != "year") {
  setnames(
    G0,
    year_col,
    "year"
  )
}

if (plant_col != "plant_id") {
  setnames(
    G0,
    plant_col,
    "plant_id"
  )
}

if (monitor_col != "radnet_monitor_id") {
  setnames(
    G0,
    monitor_col,
    "radnet_monitor_id"
  )
}

if (dist_col != "distance_km") {
  setnames(
    G0,
    dist_col,
    "distance_km"
  )
}


G0[
  ,
  `:=`(
    year = safe_int(year),
    plant_id = clean_chr(plant_id),
    radnet_monitor_id = as.character(radnet_monitor_id),
    distance_km = safe_num(distance_km)
  )
]


if (!is.na(loc_col)) {

  if (loc_col != "loc_num") {
    setnames(
      G0,
      loc_col,
      "loc_num"
    )
  }

  G0[
    ,
    loc_num := safe_int(loc_num)
  ]

  cat(
    "07A loc_num来源 = 07A原始字段\n"
  )

} else {

  G0 <- merge(
    G0,
    CROSS,
    by = "radnet_monitor_id",
    all.x = TRUE,
    sort = FALSE
  )

  cat(
    "07A loc_num来源 = 07F radnet_monitor_id crosswalk\n"
  )
}


N_UNMAPPED <- G0[
  is.na(loc_num),
  .N
]

if (N_UNMAPPED > 0L) {

  bad <- unique(
    G0[
      is.na(loc_num),
      .(
        radnet_monitor_id
      )
    ]
  )

  fwrite(
    bad,
    file.path(
      AUDIT_DIR,
      "ERROR_07A_unmapped_radnet_monitor_id.csv"
    )
  )

  stopf(
    "07A geometry存在无法映射到loc_num的RadNet ID：",
    N_UNMAPPED,
    " rows。"
  )
}


GEO <- unique(
  G0[
    year %in% YEARS &
    !is.na(loc_num) &
    !is.na(plant_id) &
    plant_id != "99999" &
    is.finite(distance_km) &
    distance_km >= 0 &
    distance_km <= MAX_RADIUS + 1e-8,
    .(
      loc_num,
      year,
      plant_id,
      distance_km
    )
  ]
)


rm(G0)
invisible(gc())


DUP_GEO <- GEO[
  ,
  .N,
  by = .(
    loc_num,
    year,
    plant_id
  )
][
  N > 1L
]


if (nrow(DUP_GEO) > 0L) {
  stopf(
    "07A geometry存在重复loc_num+year+plant_id。"
  )
}


cat(
  "200 km以内geometry rows = ",
  format(
    nrow(GEO),
    big.mark = ","
  ),
  "\n",
  "geometry monitors = ",
  uniqueN(GEO$loc_num),
  "\n",
  "geometry years = ",
  min(GEO$year),
  "–",
  max(GEO$year),
  "\n",
  sep = ""
)


# 6. Read preferred CAMPD plant-month universe
# ============================================================

cat(
  "\n============================================================\n",
  "读取冻结07B3 preferred CAMPD source\n",
  "============================================================\n",
  sep = ""
)

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


PREFERRED_CAMPD <- sort(
  unique(
    HARM[
      !is.na(
        preferred_source_resolution
      ) &
      grepl(
        "CAMPD",
        preferred_source_resolution,
        ignore.case = TRUE
      ),
      preferred_source_resolution
    ]
  )
)

if (
  length(
    PREFERRED_CAMPD
  ) ==
    0L
) {
  stopf(
    "07B3中没有识别到preferred CAMPD source。"
  )
}

cat(
  "preferred CAMPD source:\n"
)

print(
  PREFERRED_CAMPD
)


ELIGIBLE_PM <- unique(
  HARM[
    preferred_source_resolution %in%
      PREFERRED_CAMPD,
    .(
      plant_id,
      year,
      month
    )
  ]
)

rm(HARM)
invisible(gc())


# ============================================================
# 7. CAMPD fuel mapping
#    与冻结07H01保持一致
# ============================================================

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


# ============================================================
# 8. Read formal filter-period sample
#
# 完全复用已成功运行的07H09 / 07H09B / 07H10正式规则
# ============================================================

cat(
  "\n============================================================\n",
  "恢复07F正式准逐日监测站集合\n",
  "============================================================\n",
  sep = ""
)


# ------------------------------------------------------------
# 8A. 从冻结07F逐年文件读取真正的formal quasi-daily flag
# 冻结目标：426,657 rows / 133 monitors
# ------------------------------------------------------------

Q_LOC_LIST <- vector(
  "list",
  length(YEARS)
)

Q_PRIMARY_ROWS <- 0L


for (ii in seq_along(YEARS)) {

  yy <- YEARS[ii]

  qf <- file.path(
    F07_QDIR,
    paste0(
      "07F_quasi_daily_analysis_base_",
      yy,
      ".csv.gz"
    )
  )

  qh <- names(
    fread(
      qf,
      nrows = 0,
      showProgress = FALSE
    )
  )

  q_required <- c(
    "loc_num",
    "flag_07f_primary_quasi_sample"
  )

  q_missing <- setdiff(
    q_required,
    qh
  )

  if (length(q_missing) > 0L) {
    stopf(
      yy,
      "年07F准逐日文件缺少正式样本字段：",
      paste(q_missing, collapse = ", ")
    )
  }

  q <- fread(
    qf,
    select = q_required,
    showProgress = FALSE
  )

  q[
    ,
    `:=`(
      loc_num = safe_int(loc_num),
      flag_07f_primary_quasi_sample =
        as_flag(
          flag_07f_primary_quasi_sample
        )
    )
  ]

  qp <- q[
    flag_07f_primary_quasi_sample %in% TRUE
  ]

  Q_PRIMARY_ROWS <-
    Q_PRIMARY_ROWS +
    nrow(qp)

  Q_LOC_LIST[ii] <- list(
    unique(
      qp[
        ,
        .(
          loc_num
        )
      ]
    )
  )
}


ELIGIBLE_LOCS <- sort(
  unique(
    rbindlist(
      Q_LOC_LIST,
      use.names = TRUE,
      fill = TRUE
    )$loc_num
  )
)

ELIGIBLE_LOCS <- ELIGIBLE_LOCS[
  !is.na(ELIGIBLE_LOCS)
]


cat(
  "正式quasi-daily rows = ",
  format(Q_PRIMARY_ROWS, big.mark = ","),
  "\n",
  "正式quasi-daily monitors = ",
  length(ELIGIBLE_LOCS),
  "\n",
  sep = ""
)


if (Q_PRIMARY_ROWS != 426657L) {
  stopf(
    "正式quasi-daily行数错误：",
    Q_PRIMARY_ROWS,
    "；冻结值=426657"
  )
}

if (length(ELIGIBLE_LOCS) != 133L) {
  stopf(
    "正式quasi-daily监测站数错误：",
    length(ELIGIBLE_LOCS),
    "；冻结值=133"
  )
}


# ------------------------------------------------------------
# 8B. 读取全部174,200条filter-period数据
# ------------------------------------------------------------

cat(
  "\n============================================================\n",
  "读取07F filter-period数据并按冻结正式规则筛选\n",
  "============================================================\n",
  sep = ""
)


F <- fread(
  FILTER_FILE,
  showProgress = TRUE
)


if (nrow(F) != 174200L) {
  stopf(
    "filter-period全部行数错误：",
    nrow(F),
    "；冻结值=174200"
  )
}


required_filter <- c(
  "loc_num",
  "collect_start",
  "collect_end",
  "gross_beta_bq_m3",
  "gross_beta_mBq_m3",
  "c_prop",
  "surface_u_ppm_article",
  "sunspot_monthly",
  "pm_filter_period_official",
  "pm_filter_period_centered_10",
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


MISS_F <- setdiff(
  required_filter,
  names(F)
)

if (length(MISS_F) > 0L) {
  stopf(
    "07F filter-period缺少字段：\n",
    paste(MISS_F, collapse = "\n")
  )
}


F[
  ,
  `:=`(
    loc_num = safe_int(loc_num),
    collect_start = as.IDate(collect_start),
    collect_end = as.IDate(collect_end),
    flag_conus = as_flag(flag_conus)
  )
]


# ------------------------------------------------------------
# 8C. 与07H09 / 07H09B / 07H10完全相同的filter-period筛选
# ------------------------------------------------------------

F[
  ,
  flag_07h13_formal_filter := 
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


F <- F[
  flag_07h13_formal_filter %in% TRUE
]


# ------------------------------------------------------------
# 8D. 冻结硬审计
# ------------------------------------------------------------

if (nrow(F) != EXPECTED_FILTER_ROWS) {
  stopf(
    "filter-period正式样本行数错误：",
    nrow(F),
    "；预期=",
    EXPECTED_FILTER_ROWS
  )
}


if (uniqueN(F$loc_num) != EXPECTED_FILTER_MONITORS) {
  stopf(
    "filter-period正式监测站数错误：",
    uniqueN(F$loc_num),
    "；预期=",
    EXPECTED_FILTER_MONITORS
  )
}


# ------------------------------------------------------------
# 8E. 后续模型所需变量
# ------------------------------------------------------------

if (!"year" %in% names(F)) {

  F[
    ,
    year := as.integer(
      format(
        collect_start,
        "%Y"
      )
    )
  ]

} else {

  F[
    ,
    year := safe_int(year)
  ]
}


F[
  ,
  `:=`(
    period_id = .I,

    n_period_days =
      as.integer(
        collect_end - collect_start
      ),

    temp_centered =
      safe_num(temp) - 273.15,

    monitor_year_fe =
      paste(
        loc_num,
        year,
        sep = "_"
      )
  )
]


if (
  F[
    n_period_days <= 0 |
    !is.finite(n_period_days),
    .N
  ] > 0L
) {
  stopf(
    "正式filter-period存在非正采样长度。"
  )
}


cat(
  "filter-period全部数据 = 174,200 rows / 181 monitors\n",
  "filter-period正式样本 = ",
  format(nrow(F), big.mark = ","),
  " rows / ",
  uniqueN(F$loc_num),
  " monitors\n",
  sep = ""
)


# 9. Expand filter periods to [collect_start, collect_end)
# ============================================================

cat(
  "\n展开filter-period日期...\n"
)

EXPANDED <- F[
  ,
  .(
    date =
      seq(
        collect_start,
        collect_end -
          1L,
        by = "1 day"
      )
  ),
  by =
    .(
      period_id,
      loc_num,
      year,
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


EXP_CHECK <- EXPANDED[
  ,
  .N,
  by =
    .(
      period_id,
      n_period_days
    )
]


if (
  EXP_CHECK[
    N !=
      n_period_days,
    .N
  ] >
    0L
) {
  stopf(
    "filter-period日期展开数量与采样长度不一致。"
  )
}

rm(EXP_CHECK)
invisible(gc())


# ============================================================
# 10. Period-level exposure storage
# ============================================================

FPX <- data.table(
  period_id =
    F$period_id
)

for (
  g in FUELS
) {

  for (
    r in RADII
  ) {

    ncol_name <- paste0(
      "all_",
      g,
      "_count_",
      r,
      "km"
    )

    hcol_name <- paste0(
      "all_",
      g,
      "_heat_",
      r,
      "km"
    )

    FPX[
      ,
      (ncol_name) :=
        NA_real_
    ]

    FPX[
      ,
      (hcol_name) :=
        NA_real_
    ]
  }
}


# ============================================================
# 11. Read one year of actual fuel-specific CAMPD activity
# ============================================================

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
          "07B2 primary-fuel缺少字段：",
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


  DUP <- A[
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
      1L
  ]


  if (
    nrow(
      DUP
    ) >
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
# 12. Build NON-DIRECTIONAL daily exposure and aggregate
#     directly to filter period
#
# No wind sector is applied here.
# ============================================================

cat(
  "\n============================================================\n",
  "构建非方向燃料特异filter-period暴露\n",
  "============================================================\n",
  sep = ""
)


YEAR_AUDIT <- list()


for (
  yy in YEARS
) {

  cat(
    "\n[",
    yy,
    "] 构建非方向暴露...\n",
    sep = ""
  )


  EY <- EXPANDED[
    year ==
      yy
  ]


  if (
    nrow(EY) ==
      0L
  ) {
    next
  }


  DOMAIN <- unique(
    EY[
      ,
      .(
        loc_num,
        date
      )
    ]
  )


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


  cat(
    "  monitor-days=",
    format(
      nrow(DOMAIN),
      big.mark = ","
    ),
    "; linked active rows=",
    format(
      nrow(X),
      big.mark = ","
    ),
    "\n",
    sep = ""
  )


  for (
    g in FUELS
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
      r in RADII
    ) {

      ncol_name <- paste0(
        "all_",
        g,
        "_count_",
        r,
        "km"
      )

      hcol_name <- paste0(
        "all_",
        g,
        "_heat_",
        r,
        "km"
      )


      DD <- copy(
        DOMAIN
      )

      DD[
        ,
        `:=`(
          count =
            0,

          heat =
            0
        )
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
              count =
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
              heat =
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


          DD[
            CC,
            on =
              .(
                loc_num,
                date
              ),
            count :=
              i.count
          ]


          DD[
            HH,
            on =
              .(
                loc_num,
                date
              ),
            heat :=
              i.heat
          ]
        }
      }


      # ------------------------------------------
      # Expand to filter-period days and use
      # strict-complete arithmetic mean.
      # ------------------------------------------

      J <- DD[
        EY,
        on =
          .(
            loc_num,
            date
          )
      ]


      AG <- J[
        ,
        .(
          count_mean =
            strict_mean(
              count,
              first(
                n_period_days
              )
            ),

          heat_mean =
            strict_mean(
              heat,
              first(
                n_period_days
              )
            )
        ),
        by =
          period_id
      ]


      FPX[
        AG,
        on =
          "period_id",
        (ncol_name) :=
          i.count_mean
      ]


      FPX[
        AG,
        on =
          "period_id",
        (hcol_name) :=
          i.heat_mean
      ]
    }
  }


  YEAR_AUDIT[[as.character(
        yy
      )]] <- data.table(
    year =
      yy,

    n_period_day_rows =
      nrow(EY),

    n_unique_monitor_days =
      nrow(DOMAIN),

    n_geometry_rows =
      nrow(GY),

    n_actual_fuel_rows =
      nrow(AY),

    n_linked_active_rows =
      nrow(X)
  )


  rm(
    EY,
    DOMAIN,
    GY,
    AY,
    X
  )

  invisible(
    gc()
  )
}


YEAR_AUDIT <- rbindlist(
  YEAR_AUDIT,
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  YEAR_AUDIT,
  file.path(
    AUDIT_DIR,
    "01_non_directional_build_by_year.csv"
  )
)


# ============================================================
# 13. Merge period exposures into formal filter sample
# ============================================================

DF <- merge(
  F,
  FPX,
  by =
    "period_id",
  all.x =
    TRUE,
  sort =
    FALSE
)


if (
  nrow(DF) !=
    nrow(F)
) {
  stopf(
    "filter-period暴露合并改变正式样本行数。"
  )
}


# ============================================================
# 14. Fit non-directional monitor-year FE models
#
# IMPORTANT:
# Wind speed is adjusted as a covariate.
# NO exposure × wind interaction.
# ============================================================

cat(
  "\n============================================================\n",
  "拟合非方向燃料特异模型\n",
  "============================================================\n",
  sep = ""
)


fit_one <- function(
  fuel,
  radius,
  metric
) {

  exposure_col <- if (
    metric ==
      "count"
  ) {

    paste0(
      "all_",
      fuel,
      "_count_",
      radius,
      "km"
    )

  } else {

    paste0(
      "all_",
      fuel,
      "_heat_",
      radius,
      "km"
    )
  }


  d <- copy(
    DF
  )


  d[
    ,
    exposure_raw :=
      safe_num(
        get(
          exposure_col
        )
      )
  ]


  if (metric == "heat") {

    d[
      ,
      exposure_model :=
        exposure_raw /
        HEAT_INCREMENT_MMBTU
    ]

  } else {

    d[
      ,
      exposure_model :=
        exposure_raw
    ]

  }


  needed <- c(
    "gross_beta_mBq_m3",
    "exposure_model",
    "c_prop",
    "sunspot_monthly",
    "pm_filter_period_centered_10",
    "inv_hpbl",
    "rhum",
    "soilm",
    "vel",
    "temp_centered",
    "lat_diff_35",
    "surface_u_ppm_article",
    "loc_num",
    "monitor_year_fe"
  )


  keep <- complete.cases(
    d[
      ,
      ..needed
    ]
  )


  d <- d[
    keep
  ]


  model_id <- paste(
    "nondirectional",
    fuel,
    metric,
    radius,
    sep = "_"
  )


  if (
    nrow(d) <
      1000L ||
    uniqueN(
      d$loc_num
    ) <
      2L ||
    uniqueN(
      d$exposure_model
    ) <
      2L
  ) {

    return(
      data.table(
        model_id =
          model_id,

        fuel =
          fuel,

        radius_km =
          radius,

        metric =
          metric,

        status =
          "SKIP",

        n_rows =
          nrow(d),

        n_monitors =
          uniqueN(
            d$loc_num
          ),

        estimate =
          NA_real_,

        se =
          NA_real_,

        ci95_low =
          NA_real_,

        ci95_high =
          NA_real_,

        note =
          "Insufficient rows/monitors/exposure variation"
      )
    )
  }


  FORM <- as.formula(
    paste0(
      "gross_beta_mBq_m3 ~ ",
      "exposure_model + ",
      "c_prop + ",
      "sunspot_monthly + ",
      "pm_filter_period_centered_10 + ",
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
      data.table(
        model_id =
          model_id,

        fuel =
          fuel,

        radius_km =
          radius,

        metric =
          metric,

        status =
          "FAIL",

        n_rows =
          nrow(d),

        n_monitors =
          uniqueN(
            d$loc_num
          ),

        estimate =
          NA_real_,

        se =
          NA_real_,

        ci95_low =
          NA_real_,

        ci95_high =
          NA_real_,

        note =
          fit$message
      )
    )
  }


  beta <- coef(
    fit
  )

  V <- vcov(
    fit
  )


  if (
    !"exposure_model" %in%
      names(beta)
  ) {

    return(
      data.table(
        model_id =
          model_id,

        fuel =
          fuel,

        radius_km =
          radius,

        metric =
          metric,

        status =
          "SKIP",

        n_rows =
          nrow(d),

        n_monitors =
          uniqueN(
            d$loc_num
          ),

        estimate =
          NA_real_,

        se =
          NA_real_,

        ci95_low =
          NA_real_,

        ci95_high =
          NA_real_,

        note =
          "Exposure term unavailable after rank handling"
      )
    )
  }


  est <- unname(
    beta[
      "exposure_model"
    ]
  )

  se <- sqrt(
    V[
      "exposure_model",
      "exposure_model"
    ]
  )


  data.table(
    model_id =
      model_id,

    fuel =
      fuel,

    radius_km =
      radius,

    metric =
      metric,

    status =
      "PASS",

    n_rows =
      nrow(d),

    n_monitors =
      uniqueN(
        d$loc_num
      ),

    estimate =
      est,

    se =
      se,

    ci95_low =
      est -
      1.96 *
      se,

    ci95_high =
      est +
      1.96 *
      se,

    note =
      collapse_messages(
        warning_text
      )
  )
}


MODEL_LIST <- list()

kk <- 0L

TOTAL_SPECS <- length(
  FUELS
) *
  length(
    RADII
  ) *
  2L


for (
  g in FUELS
) {

  for (
    metric in c(
      "count",
      "heat"
    )
  ) {

    for (
      r in RADII
    ) {

      kk <- kk +
        1L

      cat(
        "[",
        kk,
        "/",
        TOTAL_SPECS,
        "] ",
        g,
        " | ",
        metric,
        " | ",
        r,
        " km\n",
        sep = ""
      )

      MODEL_LIST[[kk]] <- fit_one(
        g,
        r,
        metric
      )
    }
  }
}


ND_EFFECTS <- rbindlist(
  MODEL_LIST,
  use.names = TRUE,
  fill = TRUE
)


ND_EFFECTS[
  ,
  `:=`(
    fuel_label =
      unname(
        FUEL_LABEL_EN[
          fuel
        ]
      ),

    metric_label =
      fifelse(
        metric ==
          "count",
        "Active plant count",
        "Actual daily Heat Input"
      ),

    effect_unit =
      fifelse(
        metric ==
          "count",
        "per 1 additional active plant",
        paste0(
          "per ",
          format(
            HEAT_INCREMENT_MMBTU,
            big.mark = ",",
            scientific = FALSE
          ),
          " MMBtu/day"
        )
      )
  )
]


fwrite(
  ND_EFFECTS,
  file.path(
    OUT,
    "01_non_directional_filter_period_effects.csv"
  )
)


# ============================================================
# 15. Figure 4
#     NON-DIRECTIONAL fuel-specific actual daily activity
# ============================================================

cat(
  "\n============================================================\n",
  "生成 Figure 4：不考虑风向\n",
  "============================================================\n",
  sep = ""
)


F4 <- ND_EFFECTS[
  status ==
    "PASS"
]


F4[
  ,
  fuel_label :=
    factor(
      fuel_label,
      levels =
        unname(
          FUEL_LABEL_EN[
            FUELS
          ]
        )
    )
]


theme_paper <- function(
  base_size = 10
) {

  theme_bw(
    base_size =
      base_size
  ) +

    theme(
      panel.grid.minor =
        element_blank(),

      panel.grid.major =
        element_line(
          color =
            "grey90",
          linewidth =
            0.3
        ),

      strip.background =
        element_rect(
          fill =
            "white",
          color =
            "grey40"
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
            base_size +
            2
        ),

      axis.title =
        element_text(
          size =
            base_size
        ),

      axis.text =
        element_text(
          size =
            base_size -
            1
        )
    )
}


F4_COUNT <- F4[
  metric ==
    "count"
]

F4_HEAT <- F4[
  metric ==
    "heat"
]


p4a <- ggplot(
  F4_COUNT,
  aes(
    x =
      radius_km,
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
      2.5,
    linewidth =
      0.55
  ) +

  geom_line(
    linewidth =
      0.65
  ) +

  geom_point(
    size =
      1.8
  ) +

  facet_wrap(
    ~fuel_label,
    scales =
      "free_y",
    ncol =
      3
  ) +

  scale_x_continuous(
    breaks =
      RADII
  ) +

  labs(
    title =
      "A. Active plant count",

    x =
      "Cumulative radius (km)",

    y =
      "Adjusted gross-beta change (mBq/m³)\nper additional active plant"
  ) +

  theme_paper(
    10
  ) +

  theme(
    axis.text.x =
      element_text(
        angle =
          45,
        hjust =
          1
      )
  )


p4b <- ggplot(
  F4_HEAT,
  aes(
    x =
      radius_km,
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
      2.5,
    linewidth =
      0.55
  ) +

  geom_line(
    linewidth =
      0.65
  ) +

  geom_point(
    size =
      1.8
  ) +

  facet_wrap(
    ~fuel_label,
    scales =
      "free_y",
    ncol =
      3
  ) +

  scale_x_continuous(
    breaks =
      RADII
  ) +

  labs(
    title =
      "B. Actual daily Heat Input",

    x =
      "Cumulative radius (km)",

    y =
      paste0(
        "Adjusted gross-beta change (mBq/m³)\nper ",
        format(
          HEAT_INCREMENT_MMBTU,
          big.mark = ",",
          scientific = FALSE
        ),
        " MMBtu/day"
      )
  ) +

  theme_paper(
    10
  ) +

  theme(
    axis.text.x =
      element_text(
        angle =
          45,
        hjust =
          1
      )
  )


FIG4 <- p4a /
  p4b +

  plot_annotation(
    title =
      "Figure 4. Fuel-specific associations with actual daily power-plant activity regardless of wind direction",

    subtitle =
      paste0(
        "Filter-period analysis; monitor-by-year fixed effects; ",
        "circular cumulative buffers without directional filtering."
      ),

    theme =
      theme(
        plot.title =
          element_text(
            face =
              "bold",
            size =
              14
          ),

        plot.subtitle =
          element_text(
            size =
              10
          )
      )
  )


F4_PNG <- file.path(
  FIG_DIR,
  "Figure4_fuel_specific_non_directional_actual_activity.png"
)

F4_PDF <- file.path(
  FIG_DIR,
  "Figure4_fuel_specific_non_directional_actual_activity.pdf"
)


ggsave(
  F4_PNG,
  FIG4,
  width =
    13,
  height =
    11.5,
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
    F4_PDF,
    FIG4,
    width =
      13,
    height =
      11.5,
    units =
      "in",
    device =
      cairo_pdf
  )

} else {

  ggsave(
    F4_PDF,
    FIG4,
    width =
      13,
    height =
      11.5,
    units =
      "in"
  )
}


fwrite(
  F4,
  file.path(
    SOURCE_DIR,
    "Figure4_fuel_specific_non_directional_source.csv"
  )
)


# ============================================================
# 16. Figure 5
#     Existing frozen H10 UPWIND analysis
#
# Read only; do not refit H10.
# ============================================================

cat(
  "\n============================================================\n",
  "生成 Figure 5：冻结H10上风向燃料特异结果\n",
  "============================================================\n",
  sep = ""
)


E0 <- fread(
  H10_EFFECT_FILE,
  showProgress = FALSE
)

EN <- names(
  E0
)


branch_name <- pick_col(
  EN,
  c(
    "branch"
  ),
  "branch"
)

fuel_name <- pick_col(
  EN,
  c(
    "fuel"
  ),
  "fuel"
)

metric_name <- pick_col(
  EN,
  c(
    "metric"
  ),
  "metric"
)

radius_name <- pick_col(
  EN,
  c(
    "radius_km",
    "radius"
  ),
  "radius"
)

estimate_name <- pick_col(
  EN,
  c(
    "estimate",
    "effect_estimate",
    "estimate_mBq_m3"
  ),
  "estimate"
)

low_name <- pick_col(
  EN,
  c(
    "ci95_low",
    "ci_low",
    "ci95_low_mBq_m3"
  ),
  "CI lower"
)

high_name <- pick_col(
  EN,
  c(
    "ci95_high",
    "ci_high",
    "ci95_high_mBq_m3"
  ),
  "CI upper"
)

status_name <- pick_col(
  EN,
  c(
    "status"
  ),
  "status",
  required =
    FALSE
)

direction_name <- pick_col(
  EN,
  c(
    "direction"
  ),
  "direction",
  required =
    FALSE
)

wind_ref_name <- pick_col(
  EN,
  c(
    "wind_reference",
    "wind_reference_label",
    "reference_wind"
  ),
  "wind reference",
  required =
    FALSE
)


E <- data.table(

  branch =
    as.character(
      E0[[branch_name]]
    ),

  fuel =
    as.character(
      E0[[fuel_name]]
    ),

  metric =
    as.character(
      E0[[metric_name]]
    ),

  radius_km =
    safe_num(
      E0[[radius_name]]
    ),

  estimate =
    safe_num(
      E0[[estimate_name]]
    ),

  ci95_low =
    safe_num(
      E0[[low_name]]
    ),

  ci95_high =
    safe_num(
      E0[[high_name]]
    ),

  status =
    if (
      is.na(
        status_name
      )
    ) {
      rep(
        "PASS",
        nrow(E0)
      )
    } else {
      as.character(
        E0[[status_name]]
      )
    },

  direction =
    if (
      is.na(
        direction_name
      )
    ) {
      rep(
        "upwind",
        nrow(E0)
      )
    } else {
      as.character(
        E0[[direction_name]]
      )
    },

  wind_reference =
    if (
      is.na(
        wind_ref_name
      )
    ) {
      rep(
        NA_character_,
        nrow(E0)
      )
    } else {
      as.character(
        E0[[wind_ref_name]]
      )
    }
)


F5 <- E[
  branch ==
    "filter_period" &
  fuel %in%
    FUELS &
  metric %in%
    c(
      "count",
      "heat"
    ) &
  radius_km %in%
    RADII &
  direction ==
    "upwind" &
  status ==
    "PASS"
]


# If representative-effect file contains multiple reference winds,
# retain sample-median wind only.
if (
  any(
    !is.na(
      F5$wind_reference
    )
  )
) {

  median_rows <- grepl(
    "median|中位",
    F5$wind_reference,
    ignore.case = TRUE
  )

  if (
    any(
      median_rows
    )
  ) {

    F5 <- F5[
      median_rows
    ]
  }
}


F5[
  ,
  fuel_label :=
    factor(
      unname(
        FUEL_LABEL_EN[
          fuel
        ]
      ),
      levels =
        unname(
          FUEL_LABEL_EN[
            FUELS
          ]
        )
    )
]


F5_COUNT <- F5[
  metric ==
    "count"
]

F5_HEAT <- F5[
  metric ==
    "heat"
]


p5a <- ggplot(
  F5_COUNT,
  aes(
    x =
      radius_km,
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
      2.5,
    linewidth =
      0.55
  ) +

  geom_line(
    linewidth =
      0.65
  ) +

  geom_point(
    size =
      1.8
  ) +

  facet_wrap(
    ~fuel_label,
    scales =
      "free_y",
    ncol =
      3
  ) +

  scale_x_continuous(
    breaks =
      RADII
  ) +

  labs(
    title =
      "A. Active plant count",

    x =
      "Cumulative radius (km)",

    y =
      "Adjusted gross-beta change (mBq/m³)"
  ) +

  theme_paper(
    10
  ) +

  theme(
    axis.text.x =
      element_text(
        angle =
          45,
        hjust =
          1
      )
  )


p5b <- ggplot(
  F5_HEAT,
  aes(
    x =
      radius_km,
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
      2.5,
    linewidth =
      0.55
  ) +

  geom_line(
    linewidth =
      0.65
  ) +

  geom_point(
    size =
      1.8
  ) +

  facet_wrap(
    ~fuel_label,
    scales =
      "free_y",
    ncol =
      3
  ) +

  scale_x_continuous(
    breaks =
      RADII
  ) +

  labs(
    title =
      "B. Actual daily Heat Input",

    x =
      "Cumulative radius (km)",

    y =
      "Adjusted gross-beta change (mBq/m³)"
  ) +

  theme_paper(
    10
  ) +

  theme(
    axis.text.x =
      element_text(
        angle =
          45,
        hjust =
          1
      )
  )


FIG5 <- p5a /
  p5b +

  plot_annotation(
    title =
      "Figure 5. Fuel-specific upwind associations across predefined cumulative radii",

    subtitle =
      paste0(
        "Filter-period analysis; monitor-by-year fixed effects; ",
        "90° upwind directional sector; effects evaluated at the sample median wind speed."
      ),

    theme =
      theme(
        plot.title =
          element_text(
            face =
              "bold",
            size =
              14
          ),

        plot.subtitle =
          element_text(
            size =
              10
          )
      )
  )


F5_PNG <- file.path(
  FIG_DIR,
  "Figure5_fuel_specific_upwind_associations.png"
)

F5_PDF <- file.path(
  FIG_DIR,
  "Figure5_fuel_specific_upwind_associations.pdf"
)


ggsave(
  F5_PNG,
  FIG5,
  width =
    13,
  height =
    11.5,
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
    F5_PDF,
    FIG5,
    width =
      13,
    height =
      11.5,
    units =
      "in",
    device =
      cairo_pdf
  )

} else {

  ggsave(
    F5_PDF,
    FIG5,
    width =
      13,
    height =
      11.5,
    units =
      "in"
  )
}


fwrite(
  F5,
  file.path(
    SOURCE_DIR,
    "Figure5_fuel_specific_upwind_source.csv"
  )
)


# ============================================================
# 17. Hard checks
# ============================================================

ND_FAIL <- ND_EFFECTS[
  status ==
    "FAIL",
  .N
]

CORE_PASS <- ND_EFFECTS[
  fuel %in%
    CORE_FUELS &
  status ==
    "PASS",
  .N
]

EXPECTED_CORE_PASS <- length(
  CORE_FUELS
) *
  length(
    RADII
  ) *
  2L


F5_CORE_ROWS <- F5[
  fuel %in%
    CORE_FUELS,
  .N
]


CHECKS <- data.table(

  check = c(
    "filter_period_rows_107407",
    "filter_period_monitors_133",
    "CAMPD_quarters_96",
    "nondirectional_model_FAIL_zero",
    "nondirectional_core_models_complete",
    "Figure4_source_nonempty",
    "Figure4_PNG_exists",
    "Figure4_PDF_exists",
    "Figure5_core_results_present",
    "Figure5_PNG_exists",
    "Figure5_PDF_exists"
  ),

  observed = c(

    as.character(
      nrow(F)
    ),

    as.character(
      uniqueN(
        F$loc_num
      )
    ),

    as.character(
      sum(
        file.exists(
          PF_EXPECTED
        )
      )
    ),

    as.character(
      ND_FAIL
    ),

    as.character(
      CORE_PASS
    ),

    as.character(
      nrow(F4)
    ),

    as.character(
      as.integer(
        file.exists(
          F4_PNG
        )
      )
    ),

    as.character(
      as.integer(
        file.exists(
          F4_PDF
        )
      )
    ),

    as.character(
      F5_CORE_ROWS
    ),

    as.character(
      as.integer(
        file.exists(
          F5_PNG
        )
      )
    ),

    as.character(
      as.integer(
        file.exists(
          F5_PDF
        )
      )
    )
  ),

  expected = c(
    "107407",
    "133",
    "96",
    "0",
    as.character(
      EXPECTED_CORE_PASS
    ),
    ">0",
    "1",
    "1",
    ">0",
    "1",
    "1"
  )
)


CHECKS[
  ,
  status := c(

    ifelse(
      nrow(F) ==
        107407L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      uniqueN(
        F$loc_num
      ) ==
        133L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      sum(
        file.exists(
          PF_EXPECTED
        )
      ) ==
        96L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      ND_FAIL ==
        0L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      CORE_PASS ==
        EXPECTED_CORE_PASS,
      "PASS",
      "FAIL"
    ),

    ifelse(
      nrow(F4) >
        0L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      file.exists(
        F4_PNG
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      file.exists(
        F4_PDF
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      F5_CORE_ROWS >
        0L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      file.exists(
        F5_PNG
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      file.exists(
        F5_PDF
      ),
      "PASS",
      "FAIL"
    )
  )
]


fwrite(
  CHECKS,
  file.path(
    AUDIT_DIR,
    "02_07H13_hard_checks.csv"
  )
)


cat(
  "\n============================================================\n",
  "07H13 HARD CHECKS\n",
  "============================================================\n",
  sep = ""
)

print(
  CHECKS,
  nrows = Inf
)


cat(
  "\n============================================================\n",
  "NON-DIRECTIONAL MODEL STATUS\n",
  "============================================================\n",
  sep = ""
)

print(
  ND_EFFECTS[
    ,
    .N,
    by =
      status
  ]
)


N_FAIL <- CHECKS[
  status ==
    "FAIL",
  .N
]


cat(
  "\n============================================================\n",
  ifelse(
    N_FAIL ==
      0L,
    "07H13 FIGURE 4/5 GENERATION: PASS\n",
    "07H13 FIGURE 4/5 GENERATION: FAIL\n"
  ),
  "Hard FAIL = ",
  N_FAIL,
  "\n\n",
  "Figure 4 PNG = ",
  F4_PNG,
  "\n",
  "Figure 5 PNG = ",
  F5_PNG,
  "\n",
  "Non-directional effects = ",
  file.path(
    OUT,
    "01_non_directional_filter_period_effects.csv"
  ),
  "\n",
  "Audit = ",
  file.path(
    AUDIT_DIR,
    "02_07H13_hard_checks.csv"
  ),
  "\n",
  "============================================================\n",
  sep = ""
)


if (
  N_FAIL >
    0L
) {
  stopf(
    "07H13存在Hard FAIL，请勿使用Figure 4/5。"
  )
}
