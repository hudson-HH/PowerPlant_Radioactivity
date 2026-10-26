#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999,
  width = 220
)

# ============================================================
# 0. Packages
# ============================================================

required_packages <- c(
  "data.table",
  "lme4",
  "ggplot2",
  "patchwork",
  "scales"
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
  library(lme4)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

cat(
  "\n============================================================\n",
  "Figure 3：扩展半径 × 风速 × 方向性分析\n",
  "============================================================\n",
  sep = ""
)

# ============================================================
# 1. Paths
# ============================================================

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

INPUT08 <- file.path(
  ROOT,
  "model_build",
  "08_最终模型分析数据集构建",
  "01_final_quasi_daily_model_dataset_2001_2024.csv.gz"
)

INPUT09_EFFECTS <- file.path(
  ROOT,
  "model_build",
  "09_正式LME模型分析",
  "02_所有暴露模型效应估计.csv"
)

INPUT11_SECTOR <- file.path(
  ROOT,
  "model_build",
  "11_敏感性验证",
  "03_sector_angle_sensitivity_models.csv"
)

STATION_FILE <- file.path(
  ROOT,
  "radnet命名标准化",
  "radnet_181站_CONUS主分析标记与PM覆盖.csv"
)

OVERALL_LINK_FILE <- file.path(
  ROOT,
  "发电站数据处理",
  "00_既往发电站空间暴露与模型流程_冻结",
  "02_07A_RadNet发电站空间链接",
  "07A_radnet_powerplant_same_year_250km",
  "03_monitor_year_all_plant_link_250km_07A_2001_2024.csv"
)

WIND_DIR <- file.path(
  ROOT,
  "发电站数据处理",
  "00_既往发电站空间暴露与模型流程_冻结",
  "07C_daily_upwind_powerplant_exposure",
  "yearly"
)

H12 <- file.path(
  ROOT,
  "发电站数据处理",
  "06_07H_燃料特异CAMPD逐日暴露与模型",
  "07H12_论文正式图表总整理"
)

FIG_DIR <- file.path(
  H12,
  "01_正文图"
)

SOURCE_DIR <- file.path(
  H12,
  "05_绘图数据"
)

AUDIT_DIR <- file.path(
  H12,
  "06_审计"
)

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(AUDIT_DIR, recursive = TRUE, showWarnings = FALSE)

for (f in c(
  INPUT08,
  INPUT09_EFFECTS,
  INPUT11_SECTOR,
  STATION_FILE,
  OVERALL_LINK_FILE
)) {
  if (!file.exists(f)) {
    stop(
      paste0("缺少输入文件：\n", f),
      call. = FALSE
    )
  }
}

# ============================================================
# 2. Formal parameters
# ============================================================

YEARS <- 2001:2024

RADII <- c(
  20L,
  25L,
  30L,
  35L,
  40L,
  45L,
  50L,
  75L,
  100L,
  125L,
  150L
)

FIXED_WINDS <- c(
  1,
  2,
  3,
  4,
  5
)

MAIN_WIND <- 1

SECTOR_ANGLE <- 90

EXPECTED_MAIN_ROWS <- 558915L
EXPECTED_MAIN_MONITORS <- 140L

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

clean_id <- function(x) {

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
  DT,
  candidates,
  label
) {

  hit <- candidates[
    candidates %in%
      names(DT)
  ]

  if (length(hit) == 0L) {
    stopf(
      "找不到",
      label,
      "字段。候选：",
      paste(
        candidates,
        collapse = " | "
      )
    )
  }

  hit[1]
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

# 与冻结07C完全一致的球面initial-bearing公式
initial_bearing_deg_07c <- function(
  lon1,
  lat1,
  lon2,
  lat2
) {

  lon1r <- lon1 * pi / 180
  lat1r <- lat1 * pi / 180
  lon2r <- lon2 * pi / 180
  lat2r <- lat2 * pi / 180

  dlon <- lon2r - lon1r

  yy <-
    sin(dlon) *
    cos(lat2r)

  xx <-
    cos(lat1r) *
    sin(lat2r) -
    sin(lat1r) *
    cos(lat2r) *
    cos(dlon)

  bearing <-
    atan2(
      yy,
      xx
    ) *
    180 /
    pi

  (
    bearing +
      360
  ) %%
    360
}

# ============================================================
# 4. Read frozen formal model dataset
# ============================================================

cat(
  "\n============================================================\n",
  "读取冻结08正式模型样本\n",
  "============================================================\n",
  sep = ""
)

D <- fread(
  INPUT08,
  showProgress = TRUE
)

required08 <- c(
  "loc_num",
  "date",
  "gross_beta_mBq_m3",
  "surface_u_ppm_article",
  "c_prop",
  "sunspot_monthly",
  "pm_centered_10",
  "inv_hpbl",
  "rhum",
  "soilm",
  "vel",
  "lat_diff_35",
  "year_model",
  "temp",
  "flag_reference_main_sample",
  "all_upwind_n_plants_20km",
  "all_upwind_n_plants_25km",
  "all_upwind_n_plants_30km",
  "all_upwind_n_plants_35km",
  "all_upwind_n_plants_40km",
  "all_upwind_n_plants_45km",
  "all_upwind_n_plants_50km"
)

miss08 <- setdiff(
  required08,
  names(D)
)

if (length(miss08) > 0L) {
  stopf(
    "08缺少字段：\n",
    paste(
      miss08,
      collapse = "\n"
    )
  )
}

D[
  ,
  loc_num :=
    safe_int(
      loc_num
    )
]

D[
  ,
  date :=
    as.IDate(
      date
    )
]

main <- D[
  flag_reference_main_sample %in%
    TRUE
]

rm(D)
invisible(gc())

if (
  nrow(main) !=
    EXPECTED_MAIN_ROWS
) {
  stopf(
    "主分析行数不是冻结值：",
    nrow(main)
  )
}

if (
  uniqueN(main$loc_num) !=
    EXPECTED_MAIN_MONITORS
) {
  stopf(
    "主分析监测站数不是冻结值：",
    uniqueN(main$loc_num)
  )
}

if (!"year" %in% names(main)) {
  main[
    ,
    year :=
      as.integer(
        format(
          date,
          "%Y"
        )
      )
  ]
}

main[
  ,
  `:=`(
    loc_num_factor =
      factor(
        loc_num
      ),

    temp_centered =
      temp -
      273.15
  )
]

cat(
  "正式样本：",
  format(nrow(main), big.mark = ","),
  " 行；",
  uniqueN(main$loc_num),
  " 个监测站\n",
  sep = ""
)

# ============================================================
# 5. Calculate data-derived representative wind speed
#
# 不根据效应大小选择。
# 在linear exposure × wind模型中，
# effect at mean(wind) = average marginal effect.
# ============================================================

daily_wind_distribution <- unique(
  main[
    ,
    .(
      loc_num,
      date,
      vel
    )
  ]
)

daily_wind_distribution <- daily_wind_distribution[
  is.finite(vel)
]

if (
  nrow(daily_wind_distribution) ==
    0L
) {
  stopf(
    "正式样本没有有效风速。"
  )
}

REPRESENTATIVE_WIND_RAW <- mean(
  daily_wind_distribution$vel,
  na.rm = TRUE
)

# 图中固定考察1–5 m/s，所以代表值限制在同一展示区间
REPRESENTATIVE_WIND <- min(
  5,
  max(
    1,
    REPRESENTATIVE_WIND_RAW
  )
)

cat(
  "\n正式样本平均风速 = ",
  sprintf("%.4f", REPRESENTATIVE_WIND_RAW),
  " m/s\n",
  "图中数据驱动代表风速 = ",
  sprintf("%.4f", REPRESENTATIVE_WIND),
  " m/s\n",
  sep = ""
)

WIND_LEVELS <- c(
  FIXED_WINDS,
  REPRESENTATIVE_WIND
)

# 避免浮点完全重复
WIND_LEVELS <- unique(
  round(
    WIND_LEVELS,
    6
  )
)

# ============================================================
# 6. Read station crosswalk
# ============================================================

stations <- fread(
  STATION_FILE,
  showProgress = FALSE
)

required_station <- c(
  "loc_num",
  "radnet_monitor_id"
)

miss_station <- setdiff(
  required_station,
  names(stations)
)

if (length(miss_station) > 0L) {
  stopf(
    "站点表缺少字段：",
    paste(
      miss_station,
      collapse = ", "
    )
  )
}

crosswalk <- unique(
  stations[
    ,
    .(
      loc_num =
        safe_int(
          loc_num
        ),

      radnet_monitor_id =
        as.character(
          radnet_monitor_id
        )
    )
  ]
)

crosswalk <- crosswalk[
  loc_num %in%
    unique(
      main$loc_num
    )
]

# ============================================================
# 7. Read frozen 07A spatial link
# ============================================================

cat(
  "\n============================================================\n",
  "读取冻结07A空间连接并扩展到150 km\n",
  "============================================================\n",
  sep = ""
)

overall <- fread(
  OVERALL_LINK_FILE,
  colClasses =
    c(
      plant_id =
        "character"
    ),
  showProgress = FALSE
)

year_col <- pick_col(
  overall,
  c(
    "year",
    "calendar_year"
  ),
  "年份"
)

plant_col <- pick_col(
  overall,
  c(
    "plant_id",
    "orispl",
    "orispl_code"
  ),
  "发电站ID"
)

dist_col <- pick_col(
  overall,
  c(
    "distance_km"
  ),
  "距离"
)

monitor_id_col <- pick_col(
  overall,
  c(
    "radnet_monitor_id"
  ),
  "RadNet ID"
)

monitor_lat_col <- pick_col(
  overall,
  c(
    "monitor_latitude",
    "radnet_latitude",
    "monitor_lat"
  ),
  "监测站纬度"
)

monitor_lon_col <- pick_col(
  overall,
  c(
    "monitor_longitude",
    "radnet_longitude",
    "monitor_lon"
  ),
  "监测站经度"
)

plant_lat_col <- pick_col(
  overall,
  c(
    "plant_latitude",
    "analysis_latitude_final",
    "latitude_final"
  ),
  "发电站纬度"
)

plant_lon_col <- pick_col(
  overall,
  c(
    "plant_longitude",
    "analysis_longitude_final",
    "longitude_final"
  ),
  "发电站经度"
)

links <- data.table(
  year =
    safe_int(
      overall[[year_col]]
    ),

  plant_id =
    clean_id(
      overall[[plant_col]]
    ),

  radnet_monitor_id =
    as.character(
      overall[[monitor_id_col]]
    ),

  distance_km =
    safe_num(
      overall[[dist_col]]
    ),

  monitor_latitude =
    safe_num(
      overall[[monitor_lat_col]]
    ),

  monitor_longitude =
    safe_num(
      overall[[monitor_lon_col]]
    ),

  plant_latitude =
    safe_num(
      overall[[plant_lat_col]]
    ),

  plant_longitude =
    safe_num(
      overall[[plant_lon_col]]
    )
)

rm(overall)
invisible(gc())

links <- merge(
  links,
  crosswalk,
  by =
    "radnet_monitor_id",
  all.x =
    TRUE,
  sort =
    FALSE
)

links <- links[
  year %in%
    YEARS &
  loc_num %in%
    unique(main$loc_num) &
  is.finite(distance_km) &
  distance_km <=
    max(RADII) +
    1e-8
]

bad_coords <- links[
  !is.finite(monitor_latitude) |
  !is.finite(monitor_longitude) |
  !is.finite(plant_latitude) |
  !is.finite(plant_longitude),
  .N
]

if (bad_coords > 0L) {
  stopf(
    "07A空间连接存在",
    bad_coords,
    "条无效坐标。"
  )
}

links[
  ,
  bearing_deg :=
    initial_bearing_deg_07c(
      monitor_longitude,
      monitor_latitude,
      plant_longitude,
      plant_latitude
    )
]

dup_links <- links[
  ,
  .N,
  by =
    .(
      loc_num,
      year,
      plant_id
    )
][
  N > 1
]

if (nrow(dup_links) > 0L) {
  stopf(
    "07A链接存在重复loc_num+year+plant_id。"
  )
}

links <- links[
  ,
  .(
    loc_num,
    year,
    plant_id,
    distance_km,
    bearing_deg
  )
]

cat(
  "150 km以内空间连接：",
  format(
    nrow(links),
    big.mark = ","
  ),
  " 行\n",
  sep = ""
)

# ============================================================
# 8. Read frozen 07C wind
# ============================================================

cat(
  "\n============================================================\n",
  "读取冻结07C逐日风场\n",
  "============================================================\n",
  sep = ""
)

main_days <- unique(
  main[
    ,
    .(
      loc_num,
      date,
      year
    )
  ]
)

if (
  main_days[
    ,
    .N,
    by =
      .(
        loc_num,
        date
      )
  ][
    N > 1,
    .N
  ] > 0L
) {
  stopf(
    "main_days出现重复loc_num+date。"
  )
}

read_wind_year <- function(
  yy,
  target_days
) {

  path <- file.path(
    WIND_DIR,
    paste0(
      "monitor_day_wind_",
      yy,
      ".csv.gz"
    )
  )

  if (!file.exists(path)) {
    stopf(
      "缺少",
      yy,
      "年正式07C风场：",
      path
    )
  }

  W <- fread(
    path,
    showProgress = FALSE
  )

  required <- c(
    "loc_num",
    "date",
    "wind_from_deg",
    "wind_available",
    "monitor_year_07A_eligible",
    "plant_source_geography_covered"
  )

  miss <- setdiff(
    required,
    names(W)
  )

  if (length(miss) > 0L) {
    stopf(
      yy,
      "年风场缺少字段：",
      paste(
        miss,
        collapse = ", "
      )
    )
  }

  W[
    ,
    `:=`(
      loc_num =
        safe_int(
          loc_num
        ),

      date =
        as.IDate(
          date
        ),

      wind_from_deg =
        safe_num(
          wind_from_deg
        ) %%
        360,

      wind_available =
        as_flag(
          wind_available,
          FALSE
        ),

      monitor_year_07A_eligible =
        as_flag(
          monitor_year_07A_eligible,
          FALSE
        ),

      plant_source_geography_covered =
        as_flag(
          plant_source_geography_covered,
          FALSE
        )
    )
  ]

  W <- W[
    target_days,
    on =
      .(
        loc_num,
        date
      ),
    nomatch =
      0
  ]

  W[
    ,
    .(
      loc_num,
      date,
      wind_from_deg,
      wind_available,
      monitor_year_07A_eligible,
      plant_source_geography_covered
    )
  ]
}

wind_list <- vector(
  "list",
  length(YEARS)
)

names(wind_list) <- YEARS

for (yy in YEARS) {

  cat(
    "读取风场 ",
    yy,
    "...\n",
    sep = ""
  )

  wind_list[[
    as.character(yy)
  ]] <- read_wind_year(
    yy,
    main_days[
      year ==
        yy
    ]
  )
}

wind_main <- rbindlist(
  wind_list,
  use.names = TRUE,
  fill = TRUE
)

rm(wind_list)
invisible(gc())

wind_main <- merge(
  main_days,
  wind_main,
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

wind_main[
  ,
  valid_day :=
    wind_available %in%
      TRUE &
    monitor_year_07A_eligible %in%
      TRUE &
    plant_source_geography_covered %in%
      TRUE &
    is.finite(
      wind_from_deg
    )
]

if (
  wind_main[
    valid_day !=
      TRUE,
    .N
  ] > 0L
) {
  stopf(
    "冻结正式主样本中存在无效07C日期：",
    wind_main[
      valid_day != TRUE,
      .N
    ]
  )
}

# ============================================================
# 9. Rebuild daily 90° upwind plant counts
#    at all 11 requested radii
# ============================================================

cat(
  "\n============================================================\n",
  "重新构造11个半径的90°逐日上风向电厂数量\n",
  "============================================================\n",
  sep = ""
)

daily_radius <- copy(
  wind_main[
    ,
    .(
      loc_num,
      date,
      year,
      wind_from_deg,
      valid_day
    )
  ]
)

for (rr in RADII) {

  nm <- paste0(
    "all_upwind_n_plants_",
    rr,
    "km_rebuilt"
  )

  daily_radius[
    ,
    (nm) :=
      NA_real_
  ]

  daily_radius[
    valid_day %in%
      TRUE,
    (nm) :=
      0
  ]
}

monitor_years <- unique(
  daily_radius[
    ,
    .(
      loc_num,
      year
    )
  ]
)

for (
  ii in seq_len(
    nrow(monitor_years)
  )
) {

  if (
    ii %%
      100 ==
      0
  ) {
    cat(
      "处理monitor-year ",
      ii,
      " / ",
      nrow(monitor_years),
      "\n",
      sep = ""
    )
  }

  ln <- monitor_years$loc_num[ii]
  yy <- monitor_years$year[ii]

  row_idx <- which(
    daily_radius$loc_num ==
      ln &
    daily_radius$year ==
      yy &
    daily_radius$valid_day %in%
      TRUE
  )

  if (length(row_idx) == 0L) {
    next
  }

  L <- links[
    loc_num ==
      ln &
    year ==
      yy
  ]

  if (nrow(L) == 0L) {
    next
  }

  wdir <- daily_radius$wind_from_deg[
    row_idx
  ]

  delta <- abs(
    (
      outer(
        wdir,
        L$bearing_deg,
        "-"
      ) +
      180
    ) %%
      360 -
      180
  )

  upwind90 <- delta <=
    SECTOR_ANGLE /
    2

  for (rr in RADII) {

    nm <- paste0(
      "all_upwind_n_plants_",
      rr,
      "km_rebuilt"
    )

    eligible <- L$distance_km <=
      rr +
      1e-8

    if (!any(eligible)) {

      daily_radius[
        row_idx,
        (nm) :=
          0
      ]

    } else {

      daily_radius[
        row_idx,
        (nm) :=
          rowSums(
            upwind90[
              ,
              eligible,
              drop = FALSE
            ]
          )
      ]
    }
  }
}

# ============================================================
# 10. HARD AUDIT:
#     rebuilt 20–50 km must exactly equal frozen 08
# ============================================================

cat(
  "\n============================================================\n",
  "硬核对：重新构造20–50 km vs 冻结08\n",
  "============================================================\n",
  sep = ""
)

AUDIT_RADII <- c(
  20L,
  25L,
  30L,
  35L,
  40L,
  45L,
  50L
)

audit_list <- list()

for (rr in AUDIT_RADII) {

  frozen_col <- paste0(
    "all_upwind_n_plants_",
    rr,
    "km"
  )

  rebuilt_col <- paste0(
    "all_upwind_n_plants_",
    rr,
    "km_rebuilt"
  )

  frozen <- unique(
    main[
      ,
      .(
        loc_num,
        date,
        frozen =
          get(
            frozen_col
          )
      )
    ]
  )

  rebuilt <- daily_radius[
    ,
    .(
      loc_num,
      date,
      rebuilt =
        get(
          rebuilt_col
        )
    )
  ]

  aa <- merge(
    frozen,
    rebuilt,
    by =
      c(
        "loc_num",
        "date"
      ),
    all =
      TRUE
  )

  aa[
    ,
    abs_diff :=
      abs(
        frozen -
        rebuilt
      )
  ]

  audit_list[[
    as.character(rr)
  ]] <- data.table(
    radius_km =
      rr,

    n_rows =
      nrow(aa),

    n_nonzero_difference =
      aa[
        is.finite(abs_diff) &
        abs_diff >
          0,
        .N
      ],

    max_abs_difference =
      if (
        aa[
          is.finite(abs_diff),
          .N
        ] >
          0
      ) {
        max(
          aa$abs_diff,
          na.rm = TRUE
        )
      } else {
        NA_real_
      }
  )
}

radius_audit <- rbindlist(
  audit_list
)

print(
  radius_audit,
  nrows = Inf
)

if (
  radius_audit[
    n_nonzero_difference >
      0 |
    !is.finite(
      max_abs_difference
    ) |
    max_abs_difference >
      0,
    .N
  ] > 0L
) {
  stopf(
    "20–50 km重新构造结果没有100%复现冻结08；",
    "75/100/125/150 km结果禁止继续使用。"
  )
}

cat(
  "20–50 km重新构造 vs 冻结08：PASS\n"
)

# ============================================================
# 11. Merge rebuilt radii into formal main sample
# ============================================================

radius_keep_cols <- c(
  "loc_num",
  "date",
  paste0(
    "all_upwind_n_plants_",
    RADII,
    "km_rebuilt"
  )
)

main2 <- merge(
  main,
  daily_radius[
    ,
    ..radius_keep_cols
  ],
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
  nrow(main2) !=
    nrow(main)
) {
  stopf(
    "合并扩展半径后主样本行数改变。"
  )
}

# ============================================================
# 12. Exact original model specification
# ============================================================

MODEL_FORMULA <- as.formula(
  paste0(
    "gross_beta_mBq_m3 ~ ",
    "exposure_scaled + ",
    "c_prop + ",
    "surface_u_ppm_article + ",
    "sunspot_monthly + ",
    "pm_centered_10 + ",
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

control_args <- list(
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
  control_args$autoscale <- TRUE
}

LME_CONTROL <- do.call(
  lmerControl,
  control_args
)

# ============================================================
# 13. Fit 11 radius models
#     and evaluate effects at 1/2/3/4/5 + representative wind
# ============================================================

cat(
  "\n============================================================\n",
  "拟合11个半径模型\n",
  "============================================================\n",
  sep = ""
)

fit_radius <- function(rr) {

  exposure_col <- paste0(
    "all_upwind_n_plants_",
    rr,
    "km_rebuilt"
  )

  x <- safe_num(
    main2[[
      exposure_col
    ]]
  )

  finite_x <- is.finite(
    x
  )

  q95 <- as.numeric(
    quantile(
      x[
        finite_x
      ],
      probs =
        0.95,
      na.rm =
        TRUE,
      names =
        FALSE,
      type =
        7
    )
  )

  if (
    !is.finite(q95) ||
    q95 <=
      0
  ) {
    stopf(
      rr,
      " km暴露Q95无效。"
    )
  }

  d <- copy(
    main2[
      finite_x
    ]
  )

  d[
    ,
    exposure_scaled :=
      safe_num(
        get(
          exposure_col
        )
      ) /
      q95
  ]

  warning_text <- character()

  fit <- withCallingHandlers(

    lmer(
      formula =
        MODEL_FORMULA,

      data =
        d,

      REML =
        TRUE,

      na.action =
        na.fail,

      control =
        LME_CONTROL
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
  )

  b <- fixef(
    fit
  )

  V <- as.matrix(
    vcov(
      fit
    )
  )

  exposure_term <- "exposure_scaled"

  int_term <- c(
    "exposure_scaled:vel",
    "vel:exposure_scaled"
  )

  int_term <- int_term[
    int_term %in%
      names(b)
  ][1]

  if (
    is.na(int_term) ||
    length(int_term) ==
      0L
  ) {
    stopf(
      rr,
      " km模型缺少暴露×风速交互。"
    )
  }

  out <- rbindlist(
    lapply(
      WIND_LEVELS,
      function(w) {

        estimate_scaled <-
          b[
            exposure_term
          ] +
          w *
          b[
            int_term
          ]

        variance_scaled <-
          V[
            exposure_term,
            exposure_term
          ] +
          (
            w ^
              2
          ) *
          V[
            int_term,
            int_term
          ] +
          2 *
          w *
          V[
            exposure_term,
            int_term
          ]

        if (
          !is.finite(
            variance_scaled
          ) ||
          variance_scaled <
            0
        ) {
          stopf(
            "方差计算失败：",
            rr,
            " km / wind=",
            w
          )
        }

        se_scaled <- sqrt(
          variance_scaled
        )

        estimate <-
          estimate_scaled /
          q95

        se <-
          se_scaled /
          q95

        data.table(
          radius_km =
            rr,

          wind_speed_ms =
            w,

          is_representative_wind =
            abs(
              w -
              REPRESENTATIVE_WIND
            ) <
              1e-6,

          estimate_mBq_m3 =
            estimate,

          se_mBq_m3 =
            se,

          ci95_low_mBq_m3 =
            estimate -
            1.96 *
            se,

          ci95_high_mBq_m3 =
            estimate +
            1.96 *
            se,

          q95_scale =
            q95,

          n_rows =
            nrow(d),

          n_monitors =
            uniqueN(
              d$loc_num
            ),

          warnings =
            collapse_messages(
              warning_text
            )
        )
      }
    )
  )

  rm(d, fit)
  invisible(gc())

  out
}

A_DATA <- rbindlist(
  lapply(
    RADII,
    function(rr) {

      cat(
        "模型半径：",
        rr,
        " km\n",
        sep = ""
      )

      fit_radius(
        rr
      )
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

# ============================================================
# 14. Validate original 20–50 km / wind=1 against frozen 09
# ============================================================

E09 <- fread(
  INPUT09_EFFECTS,
  showProgress = FALSE
)

frozen09 <- E09[
  exposure_group ==
    "all" &
  direction ==
    "upwind" &
  metric ==
    "count" &
  radius_km %in%
    AUDIT_RADII &
  wind_speed_ms ==
    1 &
  effect_type ==
    "per_raw_unit",
  .(
    radius_km,
    frozen_estimate =
      estimate_mBq_m3,
    frozen_low =
      ci95_low_mBq_m3,
    frozen_high =
      ci95_high_mBq_m3
  )
]

refit09 <- A_DATA[
  radius_km %in%
    AUDIT_RADII &
  abs(
    wind_speed_ms -
    1
  ) <
    1e-8,
  .(
    radius_km,
    refit_estimate =
      estimate_mBq_m3,
    refit_low =
      ci95_low_mBq_m3,
    refit_high =
      ci95_high_mBq_m3
  )
]

model_audit <- merge(
  frozen09,
  refit09,
  by =
    "radius_km",
  all =
    TRUE
)

model_audit[
  ,
  `:=`(
    diff_estimate =
      abs(
        frozen_estimate -
        refit_estimate
      ),

    diff_low =
      abs(
        frozen_low -
        refit_low
      ),

    diff_high =
      abs(
        frozen_high -
        refit_high
      )
  )
]

MODEL_MAX_DIFF <- max(
  c(
    model_audit$diff_estimate,
    model_audit$diff_low,
    model_audit$diff_high
  ),
  na.rm = TRUE
)

cat(
  "\n原模型复现最大绝对差 = ",
  format(
    MODEL_MAX_DIFF,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

if (
  !is.finite(MODEL_MAX_DIFF) ||
  MODEL_MAX_DIFF >
    1e-7
) {

  print(
    model_audit
  )

  stopf(
    "扩展模型没有复现冻结09主结果。"
  )
}

# ============================================================
# 15. Panel B:
#     frozen 20-km upwind / downwind at 1 m/s
# ============================================================

B_DATA <- E09[
  exposure_group ==
    "all" &
  radius_km ==
    20 &
  metric ==
    "count" &
  wind_speed_ms ==
    1 &
  effect_type ==
    "per_raw_unit" &
  direction %in%
    c(
      "downwind",
      "upwind"
    ),
  .(
    direction,
    estimate_mBq_m3,
    ci95_low_mBq_m3,
    ci95_high_mBq_m3
  )
]

if (nrow(B_DATA) != 2L) {
  stopf(
    "B面板冻结upwind/downwind结果不是2行。"
  )
}

B_DATA[
  ,
  direction_label_en :=
    factor(
      fifelse(
        direction ==
          "upwind",
        "Upwind",
        "Downwind"
      ),
      levels =
        c(
          "Downwind",
          "Upwind"
        )
    )
]

B_DATA[
  ,
  direction_label_cn :=
    factor(
      fifelse(
        direction ==
          "upwind",
        "上风向",
        "下风向"
      ),
      levels =
        c(
          "下风向",
          "上风向"
        )
    )
]

# ============================================================
# 16. Panel C:
#     frozen 60 / 90 / 120°
# ============================================================

S11 <- fread(
  INPUT11_SECTOR,
  showProgress = FALSE
)

required_sector <- c(
  "sector_angle_deg",
  "status",
  "estimate_per_1",
  "ci95_low_per_1",
  "ci95_high_per_1"
)

miss_sector <- setdiff(
  required_sector,
  names(S11)
)

if (length(miss_sector) > 0L) {
  stopf(
    "角度敏感性文件缺少字段：",
    paste(
      miss_sector,
      collapse = ", "
    )
  )
}

C_DATA <- S11[
  status ==
    "PASS" &
  sector_angle_deg %in%
    c(
      60,
      90,
      120
    ),
  .(
    sector_angle_deg,
    estimate_mBq_m3 =
      estimate_per_1,
    ci95_low_mBq_m3 =
      ci95_low_per_1,
    ci95_high_mBq_m3 =
      ci95_high_per_1
  )
]

setorder(
  C_DATA,
  sector_angle_deg
)

if (nrow(C_DATA) != 3L) {
  stopf(
    "C面板60/90/120°结果不是3行。"
  )
}

# ============================================================
# 17. Save source data and audit
# ============================================================

fwrite(
  A_DATA,
  file.path(
    SOURCE_DIR,
    "Figure3A_20_150km_wind1_5_and_representative.csv"
  )
)

fwrite(
  B_DATA,
  file.path(
    SOURCE_DIR,
    "Figure3B_20km_upwind_downwind_count.csv"
  )
)

fwrite(
  C_DATA,
  file.path(
    SOURCE_DIR,
    "Figure3C_sector_angle_60_90_120.csv"
  )
)

fwrite(
  radius_audit,
  file.path(
    AUDIT_DIR,
    "Figure3_rebuilt_20_50km_vs_frozen08_audit.csv"
  )
)

fwrite(
  model_audit,
  file.path(
    AUDIT_DIR,
    "Figure3_refit_20_50km_vs_frozen09_audit.csv"
  )
)

wind_audit <- data.table(
  item = c(
    "sample_mean_wind_ms",
    "representative_wind_used_ms",
    "fixed_wind_1",
    "fixed_wind_2",
    "fixed_wind_3",
    "fixed_wind_4",
    "fixed_wind_5"
  ),

  value = c(
    REPRESENTATIVE_WIND_RAW,
    REPRESENTATIVE_WIND,
    1,
    2,
    3,
    4,
    5
  )
)

fwrite(
  wind_audit,
  file.path(
    AUDIT_DIR,
    "Figure3_wind_definition_audit.csv"
  )
)

# ============================================================
# 18. Plot labels
# ============================================================

REP_LABEL_EN <- paste0(
  "Representative (",
  sprintf(
    "%.2f",
    REPRESENTATIVE_WIND
  ),
  " m/s)"
)

REP_LABEL_CN <- paste0(
  "数据代表风速（",
  sprintf(
    "%.2f",
    REPRESENTATIVE_WIND
  ),
  " 米/秒）"
)

A_DATA[
  ,
  wind_label_en :=
    fifelse(
      is_representative_wind &
      !wind_speed_ms %in%
        FIXED_WINDS,
      REP_LABEL_EN,
      paste0(
        sprintf(
          "%.0f",
          wind_speed_ms
        ),
        " m/s"
      )
    )
]

A_DATA[
  ,
  wind_label_cn :=
    fifelse(
      is_representative_wind &
      !wind_speed_ms %in%
        FIXED_WINDS,
      REP_LABEL_CN,
      paste0(
        sprintf(
          "%.0f",
          wind_speed_ms
        ),
        " 米/秒"
      )
    )
]

# 如果代表风速恰好等于1–5中的一个，
# 为了仍然画出单独“代表风速”线，复制该结果一份。
if (
  any(
    abs(
      FIXED_WINDS -
      REPRESENTATIVE_WIND
    ) <
      1e-6
  )
) {

  representative_copy <- A_DATA[
    abs(
      wind_speed_ms -
      REPRESENTATIVE_WIND
    ) <
      1e-6
  ]

  representative_copy[
    ,
    `:=`(
      wind_label_en =
        REP_LABEL_EN,

      wind_label_cn =
        REP_LABEL_CN
    )
  ]

  A_DATA <- rbind(
    A_DATA,
    representative_copy,
    fill = TRUE
  )
}

# 明确排序
wind_order_en <- c(
  "1 m/s",
  "2 m/s",
  "3 m/s",
  "4 m/s",
  "5 m/s",
  REP_LABEL_EN
)

wind_order_cn <- c(
  "1 米/秒",
  "2 米/秒",
  "3 米/秒",
  "4 米/秒",
  "5 米/秒",
  REP_LABEL_CN
)

A_DATA[
  ,
  wind_label_en :=
    factor(
      wind_label_en,
      levels =
        wind_order_en
    )
]

A_DATA[
  ,
  wind_label_cn :=
    factor(
      wind_label_cn,
      levels =
        wind_order_cn
    )
]

# ============================================================
# 19. Plot style
# ============================================================

theme_fig3 <- function(
  base_size = 10.5
) {

  theme_classic(
    base_size =
      base_size
  ) +
    theme(
      plot.title =
        element_text(
          face =
            "bold",
          size =
            base_size +
            1,
          hjust =
            0
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
            0.3
        ),

      legend.title =
        element_text(
          face =
            "bold",
          size =
            base_size -
            0.2
        ),

      legend.text =
        element_text(
          size =
            base_size -
            0.6
        ),

      panel.grid.major =
        element_line(
          color =
            "grey92",
          linewidth =
            0.3
        ),

      panel.grid.minor =
        element_blank(),

      plot.margin =
        margin(
          5,
          8,
          5,
          5
        )
    )
}

COLORS_EN <- c(
  "1 m/s" =
    "#174A7E",

  "2 m/s" =
    "#2A9D8F",

  "3 m/s" =
    "#7A9E2F",

  "4 m/s" =
    "#E9A23B",

  "5 m/s" =
    "#C95D3A"
)

COLORS_EN[
  REP_LABEL_EN
] <- "#000000"

COLORS_CN <- c(
  "1 米/秒" =
    "#174A7E",

  "2 米/秒" =
    "#2A9D8F",

  "3 米/秒" =
    "#7A9E2F",

  "4 米/秒" =
    "#E9A23B",

  "5 米/秒" =
    "#C95D3A"
)

COLORS_CN[
  REP_LABEL_CN
] <- "#000000"

LINETYPES_EN <- c(
  "1 m/s" =
    "solid",

  "2 m/s" =
    "solid",

  "3 m/s" =
    "solid",

  "4 m/s" =
    "solid",

  "5 m/s" =
    "solid"
)

LINETYPES_EN[
  REP_LABEL_EN
] <- "longdash"

LINETYPES_CN <- c(
  "1 米/秒" =
    "solid",

  "2 米/秒" =
    "solid",

  "3 米/秒" =
    "solid",

  "4 米/秒" =
    "solid",

  "5 米/秒" =
    "solid"
)

LINETYPES_CN[
  REP_LABEL_CN
] <- "longdash"

LINEWIDTHS_EN <- setNames(
  rep(
    0.7,
    length(
      wind_order_en
    )
  ),
  wind_order_en
)

LINEWIDTHS_EN[
  REP_LABEL_EN
] <- 1.2

LINEWIDTHS_CN <- setNames(
  rep(
    0.7,
    length(
      wind_order_cn
    )
  ),
  wind_order_cn
)

LINEWIDTHS_CN[
  REP_LABEL_CN
] <- 1.2

# 共用纵轴
all_y <- c(
  A_DATA$ci95_low_mBq_m3,
  A_DATA$ci95_high_mBq_m3,
  B_DATA$ci95_low_mBq_m3,
  B_DATA$ci95_high_mBq_m3,
  C_DATA$ci95_low_mBq_m3,
  C_DATA$ci95_high_mBq_m3,
  0
)

YMIN <- min(
  all_y,
  na.rm = TRUE
)

YMAX <- max(
  all_y,
  na.rm = TRUE
)

PAD <- (
  YMAX -
  YMIN
) *
  0.05

if (
  !is.finite(PAD) ||
  PAD <=
    0
) {
  PAD <- 0.0001
}

YMIN <- min(
  0,
  YMIN
) -
  PAD

YMAX <- YMAX +
  PAD

# ============================================================
# 20. English Panel A
# ============================================================

pA_en <- ggplot(
  A_DATA,
  aes(
    x =
      radius_km,
    y =
      estimate_mBq_m3,
    color =
      wind_label_en,
    linetype =
      wind_label_en,
    linewidth =
      wind_label_en,
    group =
      wind_label_en
  )
) +

  geom_hline(
    yintercept =
      0,
    linewidth =
      0.4,
    color =
      "black"
  ) +

  geom_line() +

  geom_point(
    size =
      1.9
  ) +

  # 原主分析1 m/s显示CI
  geom_errorbar(
    data =
      A_DATA[
        wind_label_en ==
          "1 m/s"
      ],
    aes(
      x =
        radius_km,
      ymin =
        ci95_low_mBq_m3,
      ymax =
        ci95_high_mBq_m3
    ),
    inherit.aes =
      FALSE,
    width =
      1.1,
    linewidth =
      0.45,
    color =
      "#174A7E"
  ) +

  scale_x_continuous(
    breaks =
      RADII
  ) +

  scale_y_continuous(
    limits =
      c(
        YMIN,
        YMAX
      ),
    labels =
      label_number(
        accuracy =
          0.0001
      )
  ) +

  scale_color_manual(
    values =
      COLORS_EN,
    drop =
      FALSE
  ) +

  scale_linetype_manual(
    values =
      LINETYPES_EN,
    drop =
      FALSE
  ) +

  scale_linewidth_manual(
    values =
      LINEWIDTHS_EN,
    drop =
      FALSE
  ) +

  guides(
    linewidth =
      "none",
    linetype =
      "none"
  ) +

  labs(
    title =
      "A  Association across buffer radii at different wind speeds",

    x =
      "Buffer radius (km)",

    y =
      "Change in gross beta per additional\nupwind power plant (mBq/m³)",

    color =
      "Wind speed"
  ) +

  theme_fig3() +

  theme(
    legend.position =
      "right"
  )

# ============================================================
# 21. Chinese Panel A
# ============================================================

pA_cn <- ggplot(
  A_DATA,
  aes(
    x =
      radius_km,
    y =
      estimate_mBq_m3,
    color =
      wind_label_cn,
    linetype =
      wind_label_cn,
    linewidth =
      wind_label_cn,
    group =
      wind_label_cn
  )
) +

  geom_hline(
    yintercept =
      0,
    linewidth =
      0.4,
    color =
      "black"
  ) +

  geom_line() +

  geom_point(
    size =
      1.9
  ) +

  geom_errorbar(
    data =
      A_DATA[
        wind_label_cn ==
          "1 米/秒"
      ],
    aes(
      x =
        radius_km,
      ymin =
        ci95_low_mBq_m3,
      ymax =
        ci95_high_mBq_m3
    ),
    inherit.aes =
      FALSE,
    width =
      1.1,
    linewidth =
      0.45,
    color =
      "#174A7E"
  ) +

  scale_x_continuous(
    breaks =
      RADII
  ) +

  scale_y_continuous(
    limits =
      c(
        YMIN,
        YMAX
      ),
    labels =
      label_number(
        accuracy =
          0.0001
      )
  ) +

  scale_color_manual(
    values =
      COLORS_CN,
    drop =
      FALSE
  ) +

  scale_linetype_manual(
    values =
      LINETYPES_CN,
    drop =
      FALSE
  ) +

  scale_linewidth_manual(
    values =
      LINEWIDTHS_CN,
    drop =
      FALSE
  ) +

  guides(
    linewidth =
      "none",
    linetype =
      "none"
  ) +

  labs(
    title =
      "A  不同风速下的缓冲半径—效应关系",

    x =
      "缓冲半径（公里）",

    y =
      "每增加1座上风向发电站的\n总β变化（毫贝可/立方米）",

    color =
      "风速"
  ) +

  theme_fig3() +

  theme(
    legend.position =
      "right"
  )

# ============================================================
# 22. Panel B
# ============================================================

pB_en <- ggplot(
  B_DATA,
  aes(
    x =
      direction_label_en,
    y =
      estimate_mBq_m3
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
        ci95_low_mBq_m3,
      ymax =
        ci95_high_mBq_m3
    ),
    width =
      0.15,
    linewidth =
      0.65
  ) +

  geom_point(
    size =
      3
  ) +

  scale_y_continuous(
    limits =
      c(
        YMIN,
        YMAX
      ),
    labels =
      label_number(
        accuracy =
          0.0001
      )
  ) +

  labs(
    title =
      "B  Upwind versus downwind at 20 km",

    subtitle =
      "Wind speed = 1 m/s",

    x =
      NULL,

    y =
      "Change in gross beta\n(mBq/m³) per plant"
  ) +

  theme_fig3()

pB_cn <- ggplot(
  B_DATA,
  aes(
    x =
      direction_label_cn,
    y =
      estimate_mBq_m3
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
        ci95_low_mBq_m3,
      ymax =
        ci95_high_mBq_m3
    ),
    width =
      0.15,
    linewidth =
      0.65
  ) +

  geom_point(
    size =
      3
  ) +

  scale_y_continuous(
    limits =
      c(
        YMIN,
        YMAX
      ),
    labels =
      label_number(
        accuracy =
          0.0001
      )
  ) +

  labs(
    title =
      "B  20公里上风向与下风向比较",

    subtitle =
      "风速 = 1 米/秒",

    x =
      NULL,

    y =
      "每增加1座发电站的\n总β变化（毫贝可/立方米）"
  ) +

  theme_fig3()

# ============================================================
# 23. Panel C
# ============================================================

pC_en <- ggplot(
  C_DATA,
  aes(
    x =
      sector_angle_deg,
    y =
      estimate_mBq_m3
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
        ci95_low_mBq_m3,
      ymax =
        ci95_high_mBq_m3
    ),
    width =
      5,
    linewidth =
      0.65
  ) +

  geom_point(
    aes(
      fill =
        factor(
          sector_angle_deg ==
            90
        )
    ),
    shape =
      21,
    size =
      3,
    color =
      "black"
  ) +

  scale_fill_manual(
    values =
      c(
        "FALSE" =
          "white",

        "TRUE" =
          "#174A7E"
      ),
    guide =
      "none"
  ) +

  scale_x_continuous(
    breaks =
      c(
        60,
        90,
        120
      )
  ) +

  scale_y_continuous(
    limits =
      c(
        YMIN,
        YMAX
      ),
    labels =
      label_number(
        accuracy =
          0.0001
      )
  ) +

  labs(
    title =
      "C  Directional-sector sensitivity at 20 km",

    subtitle =
      "Wind speed = 1 m/s; 90° is the main specification",

    x =
      "Directional-sector width (°)",

    y =
      "Change in gross beta\n(mBq/m³) per plant"
  ) +

  theme_fig3()

pC_cn <- ggplot(
  C_DATA,
  aes(
    x =
      sector_angle_deg,
    y =
      estimate_mBq_m3
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
        ci95_low_mBq_m3,
      ymax =
        ci95_high_mBq_m3
    ),
    width =
      5,
    linewidth =
      0.65
  ) +

  geom_point(
    aes(
      fill =
        factor(
          sector_angle_deg ==
            90
        )
    ),
    shape =
      21,
    size =
      3,
    color =
      "black"
  ) +

  scale_fill_manual(
    values =
      c(
        "FALSE" =
          "white",

        "TRUE" =
          "#174A7E"
      ),
    guide =
      "none"
  ) +

  scale_x_continuous(
    breaks =
      c(
        60,
        90,
        120
      )
  ) +

  scale_y_continuous(
    limits =
      c(
        YMIN,
        YMAX
      ),
    labels =
      label_number(
        accuracy =
          0.0001
      )
  ) +

  labs(
    title =
      "C  20公里方向扇区角度敏感性",

    subtitle =
      "风速 = 1 米/秒；90°为主分析设定",

    x =
      "方向扇区宽度（°）",

    y =
      "每增加1座发电站的\n总β变化（毫贝可/立方米）"
  ) +

  theme_fig3()

# ============================================================
# 24. Combine Figure 3
# ============================================================

FIG3_EN <- (
  pA_en |
  pB_en |
  pC_en
) +

  plot_layout(
    widths =
      c(
        2.15,
        1,
        1
      )
  ) +

  plot_annotation(
    title =
      "Figure 3. Association between upwind power-plant count and gross-beta particle radioactivity",

    subtitle =
      paste0(
        "Effects across 20–150 km and wind speeds of 1–5 m/s; ",
        "the data-derived representative wind speed is ",
        sprintf(
          "%.2f",
          REPRESENTATIVE_WIND
        ),
        " m/s"
      ),

    caption =
      paste0(
        "The representative wind speed is the mean wind speed in the formal analysis sample and was determined independently of the outcome. ",
        "Panel A shows model-implied effects per additional upwind plant under a 90° directional sector. ",
        "Panel B compares 20-km upwind and downwind estimates at 1 m/s. ",
        "Panel C compares 60°, 90°, and 120° directional sectors at 20 km and 1 m/s."
      ),

    theme =
      theme(
        plot.title =
          element_text(
            face =
              "bold",
            size =
              15
          ),

        plot.subtitle =
          element_text(
            size =
              10.5
          ),

        plot.caption =
          element_text(
            size =
              8.2,
            hjust =
              0,
            color =
              "grey30"
          )
      )
  )

FIG3_CN <- (
  pA_cn |
  pB_cn |
  pC_cn
) +

  plot_layout(
    widths =
      c(
        2.15,
        1,
        1
      )
  ) +

  plot_annotation(
    title =
      "图3. 上风向发电站数量与总β的关联及方向性敏感性",

    subtitle =
      paste0(
        "20–150公里、1–5米/秒风速；数据计算代表风速 = ",
        sprintf(
          "%.2f",
          REPRESENTATIVE_WIND
        ),
        " 米/秒"
      ),

    caption =
      paste0(
        "数据代表风速为正式分析样本实际风速的平均值，其确定不使用总β结果。图A为90°方向扇区下每增加1座上风向发电站的模型效应；",
        "图B为20公里、1米/秒时上风向与下风向比较；图C为20公里、1米/秒时60°、90°和120°方向扇区敏感性。"
      ),

    theme =
      theme(
        plot.title =
          element_text(
            face =
              "bold",
            size =
              15
          ),

        plot.subtitle =
          element_text(
            size =
              10.5
          ),

        plot.caption =
          element_text(
            size =
              8.2,
            hjust =
              0,
            color =
              "grey30"
          )
      )
  )

# ============================================================
# 25. Save
# ============================================================

OUT_EN_PNG <- file.path(
  FIG_DIR,
  "Figure3_overall_count_radius_wind_direction_EN.png"
)

OUT_EN_PDF <- file.path(
  FIG_DIR,
  "Figure3_overall_count_radius_wind_direction_EN.pdf"
)

OUT_CN_PNG <- file.path(
  FIG_DIR,
  "Figure3_overall_count_radius_wind_direction_CN.png"
)

OUT_CN_PDF <- file.path(
  FIG_DIR,
  "Figure3_overall_count_radius_wind_direction_CN.pdf"
)

ggsave(
  OUT_EN_PNG,
  FIG3_EN,
  width =
    15.5,
  height =
    6.6,
  units =
    "in",
  dpi =
    300,
  bg =
    "white"
)

ggsave(
  OUT_CN_PNG,
  FIG3_CN,
  width =
    15.5,
  height =
    6.6,
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
    OUT_EN_PDF,
    FIG3_EN,
    width =
      15.5,
    height =
      6.6,
    units =
      "in",
    device =
      cairo_pdf
  )

  ggsave(
    OUT_CN_PDF,
    FIG3_CN,
    width =
      15.5,
    height =
      6.6,
    units =
      "in",
    device =
      cairo_pdf
  )

} else {

  ggsave(
    OUT_EN_PDF,
    FIG3_EN,
    width =
      15.5,
    height =
      6.6,
    units =
      "in"
  )

  ggsave(
    OUT_CN_PDF,
    FIG3_CN,
    width =
      15.5,
    height =
      6.6,
    units =
      "in"
  )
}

# ============================================================
# 26. Final hard audit
# ============================================================

AUDIT <- data.table(
  check = c(
    "formal_rows_558915",
    "formal_monitors_140",
    "requested_radii_11",
    "rebuild_20_50_exact_match",
    "refit_20_50_wind1_matches_09",
    "wind_1_5_present",
    "representative_wind_calculated",
    "panel_B_two_results",
    "panel_C_three_results",
    "EN_PNG_exists",
    "EN_PDF_exists",
    "CN_PNG_exists",
    "CN_PDF_exists"
  ),

  observed = c(
    as.character(
      nrow(main)
    ),

    as.character(
      uniqueN(main$loc_num)
    ),

    as.character(
      length(RADII)
    ),

    as.character(
      sum(
        radius_audit$n_nonzero_difference
      )
    ),

    as.character(
      MODEL_MAX_DIFF
    ),

    as.character(
      all(
        FIXED_WINDS %in%
          unique(
            A_DATA$wind_speed_ms
          )
      )
    ),

    sprintf(
      "%.6f",
      REPRESENTATIVE_WIND
    ),

    as.character(
      nrow(B_DATA)
    ),

    as.character(
      nrow(C_DATA)
    ),

    as.character(
      as.integer(
        file.exists(
          OUT_EN_PNG
        )
      )
    ),

    as.character(
      as.integer(
        file.exists(
          OUT_EN_PDF
        )
      )
    ),

    as.character(
      as.integer(
        file.exists(
          OUT_CN_PNG
        )
      )
    ),

    as.character(
      as.integer(
        file.exists(
          OUT_CN_PDF
        )
      )
    )
  ),

  status = c(
    ifelse(
      nrow(main) ==
        EXPECTED_MAIN_ROWS,
      "PASS",
      "FAIL"
    ),

    ifelse(
      uniqueN(main$loc_num) ==
        EXPECTED_MAIN_MONITORS,
      "PASS",
      "FAIL"
    ),

    ifelse(
      identical(
        RADII,
        c(
          20L,
          25L,
          30L,
          35L,
          40L,
          45L,
          50L,
          75L,
          100L,
          125L,
          150L
        )
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      sum(
        radius_audit$n_nonzero_difference
      ) ==
        0,
      "PASS",
      "FAIL"
    ),

    ifelse(
      is.finite(
        MODEL_MAX_DIFF
      ) &&
      MODEL_MAX_DIFF <=
        1e-7,
      "PASS",
      "FAIL"
    ),

    ifelse(
      all(
        FIXED_WINDS %in%
          unique(
            A_DATA$wind_speed_ms
          )
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      is.finite(
        REPRESENTATIVE_WIND
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      nrow(B_DATA) ==
        2L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      nrow(C_DATA) ==
        3L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      file.exists(
        OUT_EN_PNG
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      file.exists(
        OUT_EN_PDF
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      file.exists(
        OUT_CN_PNG
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      file.exists(
        OUT_CN_PDF
      ),
      "PASS",
      "FAIL"
    )
  )
)

AUDIT_FILE <- file.path(
  AUDIT_DIR,
  "Figure3_extended_radius_wind_direction_audit.csv"
)

fwrite(
  AUDIT,
  AUDIT_FILE
)

cat(
  "\n============================================================\n",
  "FIGURE 3 FINAL AUDIT\n",
  "============================================================\n",
  sep = ""
)

print(
  AUDIT,
  nrows = Inf
)

N_FAIL <- AUDIT[
  status ==
    "FAIL",
  .N
]

cat(
  "\n数据计算代表风速 = ",
  sprintf(
    "%.4f",
    REPRESENTATIVE_WIND
  ),
  " m/s\n",
  sep = ""
)

cat(
  "半径 = ",
  paste(
    RADII,
    collapse = ", "
  ),
  " km\n",
  sep = ""
)

cat(
  "固定风速 = ",
  paste(
    FIXED_WINDS,
    collapse = ", "
  ),
  " m/s\n",
  sep = ""
)

cat(
  "\n============================================================\n",
  ifelse(
    N_FAIL == 0L,
    "FIGURE 3 GENERATION: PASS\n",
    "FIGURE 3 GENERATION: FAIL\n"
  ),
  "Hard FAIL = ",
  N_FAIL,
  "\n",
  "English PNG = ",
  OUT_EN_PNG,
  "\n",
  "Chinese PNG = ",
  OUT_CN_PNG,
  "\n",
  "Audit = ",
  AUDIT_FILE,
  "\n",
  "============================================================\n",
  sep = ""
)

if (N_FAIL > 0L) {
  stopf(
    "Figure 3存在Hard FAIL，请勿使用。"
  )
}

