#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999
)

suppressPackageStartupMessages({
  library(data.table)
})


# ============================================================
# 07H04
# Fuel-specific actual-daily inverse-distance-weighted exposure
#
# Sensitivity exposure only.
# NO statistical model fitting in this stage.
#
# Exposure 1:
#   plant proximity
#   sum_j 1 / max(distance_km, 1)
#   unit = 1/km
#
# Exposure 2:
#   Heat proximity
#   sum_j Heat_j / max(distance_km, 1)
#   unit = MMBtu/km
#
# Radius selection always uses original distance_km.
# The 1-km floor is used ONLY in the inverse-distance denominator.
#
# Same as frozen 07H01:
#   - CAMPD actual daily source
#   - six formal combustion fuel groups
#   - 90-degree upwind/downwind sectors
#   - cumulative radii
#   - mixed/unknown fuel labels excluded
#   - negative Heat values preserved
#   - invalid wind/geography -> NA
#   - valid wind + no plant -> 0
#
# Upstream 07A–07H03 are read-only.
# ============================================================


# ============================================================
# 1. Constants
# ============================================================

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

PP <- file.path(
  ROOT,
  "发电站数据处理"
)

OLD_FROZEN <- file.path(
  PP,
  "00_既往发电站空间暴露与模型流程_冻结"
)

ENH_FROZEN <- file.path(
  PP,
  "01_发电站月度逐日源强增强流程_冻结"
)

F07_ROOT <- file.path(
  PP,
  "04_07F_正式统计分析数据构建"
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

H01_HARD_FILE <- file.path(
  H_ROOT,
  "qc",
  "04_build_hard_checks.csv"
)

H03_FREEZE <- file.path(
  H_ROOT,
  "07H03_模型结果审计与汇总",
  "07H03_FINAL_FREEZE.md"
)

OUT <- file.path(
  H_ROOT,
  "07H04_距离加权逐日暴露"
)

DAILY_DIR <- file.path(
  OUT,
  "weighted_daily_yearly"
)

QC_DIR <- file.path(
  OUT,
  "qc"
)


dir.create(
  DAILY_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  QC_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


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
  20,
  25,
  30,
  35,
  40,
  45,
  50,
  100,
  150,
  200
)

DIRECTIONS <- c(
  "upwind",
  "downwind"
)

HALF_ANGLE <- 45

MAX_RADIUS <- 200

DISTANCE_FLOOR_KM <- 1

TOL <- 1e-12


# ============================================================
# 2. Helpers
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
    is.na(y) |
      y == "" |
      toupper(y) %in%
        c(
          "NA",
          "N/A",
          "NULL"
        )
  ] <- NA_character_

  y
}


as_flag <- function(
  x,
  missing_value = FALSE
) {

  if (is.logical(x)) {

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

  z[
    is.na(z)
  ] <- missing_value

  z
}


sum_keep_all_na <- function(x) {

  x <- safe_num(x)

  if (
    length(x) == 0L ||
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


angular_difference <- function(
  angle1,
  angle2
) {

  abs(
    (
      (
        angle1 -
          angle2 +
          180
      ) %%
        360
    ) -
      180
  )
}


initial_bearing_deg <- function(
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


pick_name <- function(
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
        "找不到字段：",
        label,
        "\n候选：",
        paste(
          candidates,
          collapse = ", "
        )
      )
    }

    return(
      NA_character_
    )
  }

  hit[1]
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
# 3. Frozen upstream checks
# ============================================================

if (
  !file.exists(
    H01_HARD_FILE
  )
) {

  stopf(
    "找不到07H01 hard checks：\n",
    H01_HARD_FILE
  )
}


H01_HARD <- fread(
  H01_HARD_FILE,
  showProgress = FALSE
)


if (
  !all(
    c(
      "check",
      "status"
    ) %in%
      names(H01_HARD)
  )
) {

  stopf(
    "07H01 hard-check文件结构异常。"
  )
}


if (
  H01_HARD[
    status != "PASS",
    .N
  ] >
    0L
) {

  stopf(
    "07H01尚未全部PASS。"
  )
}


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
    "07H01 daily exposure不是24/24。"
  )
}


if (
  !file.exists(
    H03_FREEZE
  )
) {

  stopf(
    "找不到07H03_FINAL_FREEZE.md。\n",
    "请先运行07H03_final_freeze_audit.R。"
  )
}


freeze_text <- readLines(
  H03_FREEZE,
  warn = FALSE,
  encoding = "UTF-8"
)


if (
  !any(
    grepl(
      "07H03 FINAL FREEZE: PASS",
      freeze_text,
      fixed = TRUE
    )
  )
) {

  stopf(
    "07H03尚未FINAL FREEZE PASS。"
  )
}


F07_FREEZE <- file.path(
  F07_ROOT,
  "07F_FINAL_FREEZE.md"
)


if (
  !file.exists(
    F07_FREEZE
  )
) {

  stopf(
    "找不到07F_FINAL_FREEZE.md。"
  )
}


f07_text <- readLines(
  F07_FREEZE,
  warn = FALSE,
  encoding = "UTF-8"
)


if (
  !any(
    grepl(
      "07F FINAL FREEZE: PASS",
      f07_text,
      fixed = TRUE
    )
  )
) {

  stopf(
    "07F尚未FINAL FREEZE PASS。"
  )
}


cat(
  "07F FINAL FREEZE：PASS\n",
  "07H01 hard checks：PASS\n",
  "07H03 FINAL FREEZE：PASS\n",
  sep = ""
)


# ============================================================
# 4. Formal 07A geometry
# ============================================================

FORMAL_07A_ROOT <- file.path(
  OLD_FROZEN,
  "02_07A_RadNet发电站空间链接",
  "07A_radnet_powerplant_same_year_250km"
)


GEO_ALL_FILE <- file.path(
  FORMAL_07A_ROOT,
  "03_monitor_year_all_plant_link_250km_07A_2001_2024.csv"
)


GEO_FUEL_FILE <- file.path(
  FORMAL_07A_ROOT,
  "04_monitor_year_fuel_group_link_250km_07A_2001_2024.csv"
)


if (
  !file.exists(
    GEO_ALL_FILE
  )
) {

  stopf(
    "找不到正式07A overall geometry：\n",
    GEO_ALL_FILE
  )
}


if (
  !file.exists(
    GEO_FUEL_FILE
  )
) {

  stopf(
    "找不到正式07A fuel geometry：\n",
    GEO_FUEL_FILE
  )
}


# ============================================================
# 5. Formal 07C wind
# ============================================================

WIND_DIR <- file.path(
  OLD_FROZEN,
  "07C_daily_upwind_powerplant_exposure",
  "yearly"
)


WIND_FILES <- file.path(
  WIND_DIR,
  paste0(
    "monitor_day_wind_",
    YEARS,
    ".csv.gz"
  )
)


if (
  !all(
    file.exists(
      WIND_FILES
    )
  )
) {

  stopf(
    "正式TOP_LEVEL 07C wind不是24/24。"
  )
}


# ============================================================
# 6. Formal 07B2 / 07B3
# ============================================================

PLANT_DIR <- file.path(
  ENH_FROZEN,
  "07B2_CAMPD_daily_emissions",
  "04_plant_day"
)


PF_DIR <- file.path(
  ENH_FROZEN,
  "07B2_CAMPD_daily_emissions",
  "05_plant_day_primary_fuel"
)


B07B3_HARM_FILE <- file.path(
  ENH_FROZEN,
  "07B3_EIA_CAMPD_harmonization",
  "05_EIA_CAMPD_plant_month_harmonization_2001_2024.csv.gz"
)


if (
  !dir.exists(
    PLANT_DIR
  )
) {

  stopf(
    "找不到正式07B2 04_plant_day：\n",
    PLANT_DIR
  )
}


if (
  !dir.exists(
    PF_DIR
  )
) {

  stopf(
    "找不到正式07B2 05_plant_day_primary_fuel：\n",
    PF_DIR
  )
}


if (
  !file.exists(
    B07B3_HARM_FILE
  )
) {

  stopf(
    "找不到正式07B3 harmonization：\n",
    B07B3_HARM_FILE
  )
}


expected_tags <- sort(
  unlist(
    lapply(
      YEARS,
      function(y) {
        paste0(
          y,
          "_q",
          1:4
        )
      }
    )
  )
)


check_quarter_dir <- function(
  dir,
  prefix,
  label
) {

  fs <- list.files(
    dir,
    pattern = paste0(
      "^",
      prefix,
      "[0-9]{4}_q[1-4]\\.csv\\.gz$"
    ),
    full.names = TRUE
  )

  tag <- basename(fs)

  tag <- sub(
    paste0(
      "^",
      prefix
    ),
    "",
    tag
  )

  tag <- sub(
    "\\.csv\\.gz$",
    "",
    tag
  )

  miss <- setdiff(
    expected_tags,
    tag
  )

  extra <- setdiff(
    tag,
    expected_tags
  )

  dup <- unique(
    tag[
      duplicated(tag)
    ]
  )

  if (
    length(fs) != 96L ||
    length(miss) > 0L ||
    length(extra) > 0L ||
    length(dup) > 0L
  ) {

    stopf(
      label,
      "季度文件不是完整2001–2024 × 4。"
    )
  }

  invisible(TRUE)
}


check_quarter_dir(
  PLANT_DIR,
  "CAMPD_plant_day_",
  "07B2 04 plant-day"
)


check_quarter_dir(
  PF_DIR,
  "CAMPD_plant_day_primary_fuel_",
  "07B2 05 primary-fuel"
)


# ============================================================
# 7. loc_num <-> RadNet crosswalk from frozen 07F
# ============================================================

cross_parts <- vector(
  "list",
  length(YEARS)
)


for (
  ii in seq_along(YEARS)
) {

  yy <- YEARS[ii]

  qf <- file.path(
    F07_QDIR,
    paste0(
      "07F_quasi_daily_analysis_base_",
      yy,
      ".csv.gz"
    )
  )

  if (
    !file.exists(qf)
  ) {

    stopf(
      "缺少07F year：",
      yy
    )
  }

  qh <- names(
    fread(
      qf,
      nrows = 0,
      showProgress = FALSE
    )
  )

  req <- c(
    "loc_num",
    "radnet_monitor_id"
  )

  miss <- setdiff(
    req,
    qh
  )

  if (
    length(miss) >
      0L
  ) {

    stopf(
      yy,
      " 07F缺少crosswalk字段：",
      paste(
        miss,
        collapse = ", "
      )
    )
  }

  z <- fread(
    qf,
    select = req,
    showProgress = FALSE
  )

  cross_parts[[ii]] <- unique(
    z[
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
}


CROSS <- unique(
  rbindlist(
    cross_parts
  )
)


if (
  CROSS[
    ,
    .N,
    by = loc_num
  ][
    N > 1L,
    .N
  ] >
    0L ||
  CROSS[
    ,
    .N,
    by = radnet_monitor_id
  ][
    N > 1L,
    .N
  ] >
    0L
) {

  stopf(
    "07F loc_num <-> radnet_monitor_id不是一一对应。"
  )
}


# ============================================================
# 8. Standardize formal 07A geometry
# ============================================================

standardize_geometry <- function(file) {

  H <- names(
    fread(
      file,
      nrows = 0,
      showProgress = FALSE
    )
  )

  year_col <- pick_name(
    H,
    c(
      "year",
      "calendar_year"
    ),
    "year"
  )

  plant_col <- pick_name(
    H,
    c(
      "plant_id",
      "orispl",
      "orispl_code"
    ),
    "plant_id"
  )

  monitor_col <- pick_name(
    H,
    c(
      "radnet_monitor_id"
    ),
    "radnet_monitor_id"
  )

  loc_col <- pick_name(
    H,
    c(
      "loc_num"
    ),
    "loc_num",
    FALSE
  )

  dist_col <- pick_name(
    H,
    c(
      "distance_km",
      "distance"
    ),
    "distance_km"
  )

  bearing_col <- pick_name(
    H,
    c(
      "plant_bearing_deg",
      "bearing_deg",
      "monitor_to_plant_bearing_deg"
    ),
    "bearing",
    FALSE
  )

  need <- unique(
    na.omit(
      c(
        year_col,
        plant_col,
        monitor_col,
        loc_col,
        dist_col,
        bearing_col
      )
    )
  )

  if (
    is.na(
      bearing_col
    )
  ) {

    monitor_lat_col <- pick_name(
      H,
      c(
        "monitor_latitude",
        "radnet_latitude",
        "monitor_lat"
      ),
      "monitor latitude"
    )

    monitor_lon_col <- pick_name(
      H,
      c(
        "monitor_longitude",
        "radnet_longitude",
        "monitor_lon"
      ),
      "monitor longitude"
    )

    plant_lat_col <- pick_name(
      H,
      c(
        "plant_latitude",
        "analysis_latitude_final",
        "latitude_final"
      ),
      "plant latitude"
    )

    plant_lon_col <- pick_name(
      H,
      c(
        "plant_longitude",
        "analysis_longitude_final",
        "longitude_final"
      ),
      "plant longitude"
    )

    need <- unique(
      c(
        need,
        monitor_lat_col,
        monitor_lon_col,
        plant_lat_col,
        plant_lon_col
      )
    )
  }

  D <- fread(
    file,
    select = need,
    showProgress = FALSE
  )

  if (
    year_col !=
      "year"
  ) {
    setnames(
      D,
      year_col,
      "year"
    )
  }

  if (
    plant_col !=
      "plant_id"
  ) {
    setnames(
      D,
      plant_col,
      "plant_id"
    )
  }

  if (
    monitor_col !=
      "radnet_monitor_id"
  ) {
    setnames(
      D,
      monitor_col,
      "radnet_monitor_id"
    )
  }

  if (
    dist_col !=
      "distance_km"
  ) {
    setnames(
      D,
      dist_col,
      "distance_km"
    )
  }

  D[
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
        D,
        loc_col,
        "loc_num"
      )
    }

    D[
      ,
      loc_num :=
        safe_int(
          loc_num
        )
    ]

  } else {

    D <- merge(
      D,
      CROSS,
      by = "radnet_monitor_id",
      all.x = TRUE,
      sort = FALSE
    )
  }

  if (
    any(
      is.na(
        D$loc_num
      )
    )
  ) {

    stopf(
      "07A geometry存在无法映射的RadNet ID。"
    )
  }

  if (
    is.na(
      bearing_col
    )
  ) {

    D[
      ,
      bearing_deg :=
        initial_bearing_deg(
          safe_num(
            get(
              monitor_lon_col
            )
          ),
          safe_num(
            get(
              monitor_lat_col
            )
          ),
          safe_num(
            get(
              plant_lon_col
            )
          ),
          safe_num(
            get(
              plant_lat_col
            )
          )
        )
    ]

  } else {

    D[
      ,
      bearing_deg :=
        safe_num(
          get(
            bearing_col
          )
        ) %%
        360
    ]
  }

  D <- D[
    year %in%
      YEARS &
    !is.na(
      plant_id
    ) &
    plant_id !=
      "99999" &
    is.finite(
      distance_km
    ) &
    distance_km >=
      0 &
    distance_km <=
      MAX_RADIUS +
      1e-8 &
    is.finite(
      bearing_deg
    )
  ]

  D <- unique(
    D[
      ,
      .(
        loc_num,
        radnet_monitor_id,
        year,
        plant_id,
        distance_km,
        bearing_deg
      )
    ]
  )

  dup <- D[
    ,
    .N,
    by = .(
      loc_num,
      year,
      plant_id
    )
  ][
    N >
      1L
  ]

  if (
    nrow(dup) >
      0L
  ) {

    stopf(
      "07A overall geometry存在重复loc_num + year + plant_id。"
    )
  }

  D[
    ,
    distance_weight_km :=
      pmax(
        distance_km,
        DISTANCE_FLOOR_KM
      )
  ]

  D[
    ,
    inv_distance_weight :=
      1 /
      distance_weight_km
  ]

  D
}


cat(
  "\n读取正式冻结07A geometry...\n"
)


GEO_ALL <- standardize_geometry(
  GEO_ALL_FILE
)


cat(
  "overall links <=200 km = ",
  fmt_int(
    nrow(
      GEO_ALL
    )
  ),
  "\n",
  "links with distance <1 km = ",
  fmt_int(
    GEO_ALL[
      distance_km <
        DISTANCE_FLOOR_KM,
      .N
    ]
  ),
  "\n",
  sep = ""
)


# ============================================================
# 9. CAMPD fuel crosswalk
# ============================================================

CAMPD_COMPONENT_MAP <- c(

  "COAL" =
    "coal",

  "COAL REFUSE" =
    "coal",

  "PIPELINE NATURAL GAS" =
    "gas",

  "NATURAL GAS" =
    "gas",

  "PROCESS GAS" =
    "gas",

  "OTHER GAS" =
    "gas",

  "RESIDUAL OIL" =
    "oil",

  "DIESEL OIL" =
    "oil",

  "OTHER OIL" =
    "oil",

  "PETROLEUM COKE" =
    "oil",

  "WOOD" =
    "biomass",

  "TIRE DERIVED FUEL" =
    "solid_waste",

  "TIRE-DERIVED FUEL" =
    "solid_waste",

  "OTHER SOLID FUEL" =
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
          "blank_unknown",

        components =
          NA_character_,

        component_groups =
          NA_character_
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
          "unmapped_component",

        components =
          paste(
            parts,
            collapse = "|"
          ),

        component_groups =
          paste(
            ifelse(
              is.na(groups),
              "UNMAPPED",
              groups
            ),
            collapse = "|"
          )
      )
    )
  }

  unique_groups <- sort(
    unique(
      groups
    )
  )

  if (
    length(
      unique_groups
    ) ==
      1L
  ) {

    method <- if (
      length(parts) ==
        1L
    ) {
      "single_component"
    } else {
      "same_group_multi_component"
    }

    return(
      data.table(
        fuel_group =
          unique_groups[1],

        mapping_method =
          method,

        components =
          paste(
            parts,
            collapse = "|"
          ),

        component_groups =
          paste(
            groups,
            collapse = "|"
          )
      )
    )
  }

  data.table(
    fuel_group =
      "mixed",

    mapping_method =
      "cross_group_mixed",

    components =
      paste(
        parts,
        collapse = "|"
      ),

    component_groups =
      paste(
        unique_groups,
        collapse = "|"
      )
  )
}


# ============================================================
# 10. 07B3 preferred CAMPD source
# ============================================================

HARM_SOURCE <- fread(
  B07B3_HARM_FILE,
  select = c(
    "plant_id",
    "year",
    "month",
    "source_resolution_status",
    "preferred_source_resolution"
  ),
  showProgress = FALSE
)


HARM_SOURCE[
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

    source_resolution_status =
      clean_chr(
        source_resolution_status
      ),

    preferred_source_resolution =
      clean_chr(
        preferred_source_resolution
      )
  )
]


PREFERRED_VALUES <- sort(
  unique(
    HARM_SOURCE[
      !is.na(
        preferred_source_resolution
      ),
      preferred_source_resolution
    ]
  )
)


cat(
  "\n07B3 preferred_source_resolution实际类别：\n"
)

print(
  PREFERRED_VALUES
)


if (
  !"CAMPD_DAILY_ACTUAL" %in%
    PREFERRED_VALUES
) {

  stopf(
    "07B3中没有CAMPD_DAILY_ACTUAL。"
  )
}


CAMPD_ELIGIBLE_PM <- unique(
  HARM_SOURCE[
    preferred_source_resolution ==
      "CAMPD_DAILY_ACTUAL",
    .(
      plant_id,
      year,
      month
    )
  ]
)


dup_pm <- CAMPD_ELIGIBLE_PM[
  ,
  .N,
  by = .(
    plant_id,
    year,
    month
  )
][
  N >
    1L
]


if (
  nrow(
    dup_pm
  ) >
    0L
) {

  stopf(
    "07B3 CAMPD_DAILY_ACTUAL plant-month主键重复。"
  )
}


# ============================================================
# 11. Read wind
# ============================================================

read_wind_year <- function(
  yy,
  target
) {

  f <- file.path(
    WIND_DIR,
    paste0(
      "monitor_day_wind_",
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
    H
  )

  if (
    length(miss) >
      0L
  ) {

    stopf(
      yy,
      " 正式07C wind缺少：",
      paste(
        miss,
        collapse = ", "
      )
    )
  }

  W <- fread(
    f,
    select = required,
    showProgress = FALSE
  )

  W[
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

      wind_from_deg =
        safe_num(
          wind_from_deg
        ),

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

  W <- merge(
    target,
    W,
    by = c(
      "loc_num",
      "date"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  W[
    ,
    wind_valid :=
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

  W
}


# ============================================================
# 12. Read yearly CAMPD actual fuel activity
# ============================================================

MAPPING_AUDIT_LIST <- list()

SOURCE_AUDIT_LIST <- list()


read_pf_year <- function(
  yy
) {

  pf_files <- file.path(
    PF_DIR,
    paste0(
      "CAMPD_plant_day_primary_fuel_",
      yy,
      "_q",
      1:4,
      ".csv.gz"
    )
  )

  plant_files <- file.path(
    PLANT_DIR,
    paste0(
      "CAMPD_plant_day_",
      yy,
      "_q",
      1:4,
      ".csv.gz"
    )
  )


  pf_parts <- lapply(
    pf_files,
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
          "07B2 05 primary-fuel缺少字段：",
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
    pf_parts,
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
          as.character(
            date
          )
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


  plant_parts <- lapply(
    plant_files,
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
          "07B2 04 plant-day缺少字段：",
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


  P <- rbindlist(
    plant_parts,
    use.names = TRUE,
    fill = TRUE
  )


  P[
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
          as.character(
            date
          )
        ),

      plant_id =
        clean_chr(
          plant_id
        ),

      heat =
        safe_num(
          daily_heat_input_mmbtu_observed
        )
    )
  ]


  eligible_pm <- CAMPD_ELIGIBLE_PM[
    year ==
      yy
  ]


  D <- D[
    eligible_pm,
    on = .(
      plant_id,
      year,
      month
    ),
    nomatch = 0
  ]


  P <- P[
    eligible_pm,
    on = .(
      plant_id,
      year,
      month
    ),
    nomatch = 0
  ]


  relevant_py <- unique(
    GEO_ALL[
      year ==
        yy,
      .(
        year,
        plant_id
      )
    ]
  )


  D <- D[
    relevant_py,
    on = .(
      year,
      plant_id
    ),
    nomatch = 0
  ]


  P <- P[
    relevant_py,
    on = .(
      year,
      plant_id
    ),
    nomatch = 0
  ]


  dup_plant_day <- P[
    ,
    .N,
    by = .(
      plant_id,
      date
    )
  ][
    N >
      1L
  ]


  if (
    nrow(
      dup_plant_day
    ) >
      0L
  ) {

    stopf(
      yy,
      " 07B2 04 plant-day存在重复主键。"
    )
  }


  D[
    ,
    fuel_label_key :=
      fifelse(
        is.na(
          primary_fuel_info
        ),
        "<BLANK>",
        primary_fuel_info
      )
  ]


  raw_dict <- unique(
    D[
      ,
      .(
        fuel_label_key,
        primary_fuel_info
      )
    ]
  )


  mapped_parts <- lapply(
    raw_dict$primary_fuel_info,
    classify_campd_fuel_label
  )


  mapped <- rbindlist(
    mapped_parts,
    use.names = TRUE,
    fill = TRUE
  )


  raw_dict <- cbind(
    raw_dict,
    mapped
  )


  D <- merge(
    D,
    raw_dict[
      ,
      .(
        fuel_label_key,
        fuel_group,
        mapping_method,
        components,
        component_groups
      )
    ],
    by = "fuel_label_key",
    all.x = TRUE,
    sort = FALSE
  )


  unmapped_nonblank <- D[
    mapping_method ==
      "unmapped_component",
    .N
  ]


  if (
    unmapped_nonblank >
      0L
  ) {

    bad <- D[
      mapping_method ==
        "unmapped_component"
    ]

    fwrite(
      bad,
      file.path(
        QC_DIR,
        paste0(
          "UNMAPPED_CAMPD_FUEL_COMPONENT_",
          yy,
          ".csv.gz"
        )
      )
    )

    stopf(
      yy,
      "存在未纳入正式crosswalk的非空CAMPD燃料标签：",
      unmapped_nonblank,
      " rows。"
    )
  }


  mapping_audit <- D[
    ,
    .(
      n_rows =
        .N,

      n_unique_plants =
        uniqueN(
          plant_id
        ),

      n_heat_nonmissing =
        sum(
          is.finite(
            heat
          )
        ),

      total_heat_mmbtu =
        sum(
          heat[
            is.finite(
              heat
            )
          ]
        )
    ),
    by = .(
      year,
      primary_fuel_info =
        fuel_label_key,
      fuel_group,
      mapping_method,
      components,
      component_groups
    )
  ]


  mapping_audit[
    ,
    included_in_formal_fuel_analysis :=
      fuel_group %in%
        FUELS
  ]


  MAPPING_AUDIT_LIST[[
    length(
      MAPPING_AUDIT_LIST
    ) +
      1L
  ]] <<-
    mapping_audit


  PKEY <- unique(
    P[
      ,
      .(
        plant_id,
        date
      )
    ]
  )


  DKEY <- unique(
    D[
      ,
      .(
        plant_id,
        date
      )
    ]
  )


  missing_in_05 <- PKEY[
    !DKEY,
    on = .(
      plant_id,
      date
    )
  ]


  extra_in_05 <- DKEY[
    !PKEY,
    on = .(
      plant_id,
      date
    )
  ]


  source_audit <- data.table(
    year =
      yy,

    n_04_plant_days =
      nrow(
        PKEY
      ),

    n_05_plant_days =
      nrow(
        DKEY
      ),

    n_04_missing_from_05 =
      nrow(
        missing_in_05
      ),

    n_05_not_in_04 =
      nrow(
        extra_in_05
      )
  )


  SOURCE_AUDIT_LIST[[
    length(
      SOURCE_AUDIT_LIST
    ) +
      1L
  ]] <<-
    source_audit


  FUEL_ACTIVITY <- D[
    fuel_group %in%
      FUELS,
    .(
      heat =
        sum_keep_all_na(
          heat
        )
    ),
    by = .(
      plant_id,
      date,
      fuel_group
    )
  ]


  dup_activity <- FUEL_ACTIVITY[
    ,
    .N,
    by = .(
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
      dup_activity
    ) >
      0L
  ) {

    stopf(
      yy,
      " fuel-specific plant-day聚合后仍有重复。"
    )
  }


  list(
    fuel_activity =
      FUEL_ACTIVITY,

    source_audit =
      source_audit
  )
}


# ============================================================
# 13. Build one monitor-year
# ============================================================

build_monitor_year <- function(
  loc,
  W,
  links_all,
  fuel_activity
) {

  T <- copy(
    W[
      loc_num ==
        loc
    ]
  )

  setorder(
    T,
    date
  )


  D <- T[
    ,
    .(
      loc_num,
      date,
      wind_valid
    )
  ]


  for (
    g in FUELS
  ) {

    for (
      dir in DIRECTIONS
    ) {

      for (
        r in RADII
      ) {

        pc <- paste0(
          "campd_",
          g,
          "_",
          dir,
          "_proximity_invkm_",
          r,
          "km_90deg"
        )

        hc <- paste0(
          "campd_",
          g,
          "_",
          dir,
          "_heat_proximity_mmbtu_per_km_",
          r,
          "km_90deg"
        )

        D[
          ,
          (pc) :=
            NA_real_
        ]

        D[
          ,
          (hc) :=
            NA_real_
        ]

        D[
          wind_valid %in%
            TRUE,
          (pc) :=
            0
        ]

        D[
          wind_valid %in%
            TRUE,
          (hc) :=
            0
        ]
      }
    }
  }


  LALL <- links_all[
    loc_num ==
      loc
  ]


  if (
    nrow(
      LALL
    ) >
      0L &&
    nrow(
      fuel_activity
    ) >
      0L
  ) {

    plant_ids <- unique(
      LALL$plant_id
    )

    FH <- fuel_activity[
      plant_id %in%
        plant_ids
    ]


    if (
      nrow(
        FH
      ) >
        0L
    ) {

      HH <- merge(
        FH,
        LALL[
          ,
          .(
            plant_id,
            distance_km,
            distance_weight_km,
            inv_distance_weight,
            bearing_deg
          )
        ],
        by = "plant_id",
        all = FALSE,
        allow.cartesian = TRUE,
        sort = FALSE
      )


      HH <- merge(
        HH,
        T[
          ,
          .(
            date,
            wind_from_deg,
            wind_valid
          )
        ],
        by = "date",
        all = FALSE,
        sort = FALSE
      )


      HH <- HH[
        wind_valid %in%
          TRUE
      ]


      HH[
        ,
        `:=`(
          upwind_flag =
            angular_difference(
              bearing_deg,
              wind_from_deg
            ) <=
            HALF_ANGLE,

          downwind_flag =
            angular_difference(
              bearing_deg,
              (
                wind_from_deg +
                  180
              ) %%
                360
            ) <=
            HALF_ANGLE
        )
      ]

    } else {

      HH <- data.table()
    }

  } else {

    HH <- data.table()
  }


  for (
    g in FUELS
  ) {

    for (
      dir in DIRECTIONS
    ) {

      flag_col <- paste0(
        dir,
        "_flag"
      )


      for (
        r in RADII
      ) {

        pc <- paste0(
          "campd_",
          g,
          "_",
          dir,
          "_proximity_invkm_",
          r,
          "km_90deg"
        )

        hc <- paste0(
          "campd_",
          g,
          "_",
          dir,
          "_heat_proximity_mmbtu_per_km_",
          r,
          "km_90deg"
        )


        if (
          nrow(
            HH
          ) >
            0L
        ) {

          PX <- HH[
            fuel_group ==
              g &
            distance_km <=
              r +
              1e-8 &
            get(
              flag_col
            ) %in%
              TRUE,
            .(
              value =
                sum(
                  inv_distance_weight
                )
            ),
            by = date
          ]


          if (
            nrow(
              PX
            ) >
              0L
          ) {

            D[
              PX,
              on = "date",
              (pc) :=
                i.value
            ]
          }


          HX <- HH[
            fuel_group ==
              g &
            distance_km <=
              r +
              1e-8 &
            get(
              flag_col
            ) %in%
              TRUE,
            .(
              value =
                sum_keep_all_na(
                  heat *
                    inv_distance_weight
                )
            ),
            by = date
          ]


          if (
            nrow(
              HX
            ) >
              0L
          ) {

            D[
              HX,
              on = "date",
              (hc) :=
                i.value
            ]
          }
        }
      }
    }
  }


  D
}


# ============================================================
# 14. Build 2001–2024
# ============================================================

YEAR_QC_LIST <- list()

SEMANTIC_QC_LIST <- list()


for (
  yy in YEARS
) {

  cat(
    "\n============================================================\n",
    "07H04 YEAR ",
    yy,
    "\n",
    "============================================================\n",
    sep = ""
  )


  h01_file <- file.path(
    H01_DAILY_DIR,
    paste0(
      "07H_fuel_daily_exposure_",
      yy,
      ".csv.gz"
    )
  )


  H1 <- fread(
    h01_file,
    showProgress = FALSE
  )


  required_h01_base <- c(
    "loc_num",
    "date",
    "wind_valid"
  )


  miss_h01_base <- setdiff(
    required_h01_base,
    names(H1)
  )


  if (
    length(
      miss_h01_base
    ) >
      0L
  ) {

    stopf(
      yy,
      " 07H01 daily缺少：",
      paste(
        miss_h01_base,
        collapse = ", "
      )
    )
  }


  H1[
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

      wind_valid =
        as_flag(
          wind_valid,
          FALSE
        )
    )
  ]


  dup_h01 <- H1[
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
    dup_h01 >
      0L
  ) {

    stopf(
      yy,
      " 07H01 daily存在重复monitor-day。"
    )
  }


  TARGET <- H1[
    ,
    .(
      loc_num,
      date
    )
  ]


  W <- read_wind_year(
    yy,
    TARGET
  )


  if (
    nrow(W) !=
      nrow(TARGET)
  ) {

    stopf(
      yy,
      " wind目标域行数不一致。"
    )
  }


  wind_mismatch <- merge(
    H1[
      ,
      .(
        loc_num,
        date,
        h01_wind_valid =
          wind_valid
      )
    ],
    W[
      ,
      .(
        loc_num,
        date,
        reconstructed_wind_valid =
          wind_valid
      )
    ],
    by = c(
      "loc_num",
      "date"
    ),
    all = TRUE,
    sort = FALSE
  )[
    is.na(
      h01_wind_valid
    ) |
    is.na(
      reconstructed_wind_valid
    ) |
    h01_wind_valid !=
      reconstructed_wind_valid,
    .N
  ]


  if (
    wind_mismatch >
      0L
  ) {

    stopf(
      yy,
      " 07H01 wind_valid与正式07C重构不一致：",
      wind_mismatch
    )
  }


  PF <- read_pf_year(
    yy
  )


  LA <- GEO_ALL[
    year ==
      yy
  ]


  locs <- sort(
    unique(
      TARGET$loc_num
    )
  )


  OUT_PARTS <- vector(
    "list",
    length(locs)
  )


  for (
    jj in seq_along(locs)
  ) {

    OUT_PARTS[[jj]] <- build_monitor_year(
      loc =
        locs[jj],

      W =
        W,

      links_all =
        LA,

      fuel_activity =
        PF$fuel_activity
    )
  }


  DY <- rbindlist(
    OUT_PARTS,
    use.names = TRUE,
    fill = TRUE
  )


  setorder(
    DY,
    loc_num,
    date
  )


  dup_day <- DY[
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


  exposure_cols <- grep(
    "^campd_(coal|gas|oil|other_fossil|biomass|solid_waste)_(upwind|downwind)_(proximity_invkm|heat_proximity_mmbtu_per_km)_",
    names(DY),
    value = TRUE
  )


  invalid_non_na <- 0L

  invalid_idx <- !(
    DY$wind_valid %in%
      TRUE
  )


  if (
    any(
      invalid_idx
    )
  ) {

    for (
      cc in exposure_cols
    ) {

      invalid_non_na <-
        invalid_non_na +
        sum(
          !is.na(
            DY[[cc]][invalid_idx]
          )
        )
    }
  }


  inf_cells <- 0L


  for (
    cc in exposure_cols
  ) {

    inf_cells <-
      inf_cells +
      sum(
        is.infinite(
          DY[[cc]]
        )
      )
  }


  proximity_monotonic_fail <- 0L


  for (
    g in FUELS
  ) {

    for (
      dir in DIRECTIONS
    ) {

      pc <- paste0(
        "campd_",
        g,
        "_",
        dir,
        "_proximity_invkm_",
        RADII,
        "km_90deg"
      )


      M <- as.matrix(
        DY[
          wind_valid %in%
            TRUE,
          ..pc
        ]
      )


      if (
        nrow(M) >
          0L
      ) {

        proximity_monotonic_fail <-
          proximity_monotonic_fail +
          sum(
            apply(
              M,
              1,
              function(x) {

                any(
                  diff(x) <
                    -TOL,
                  na.rm = TRUE
                )
              }
            )
          )
      }
    }
  }


  # ==========================================================
  # Exact semantic comparison against frozen 07H01
  # ==========================================================

  J <- merge(
    H1,
    DY,
    by = c(
      "loc_num",
      "date",
      "wind_valid"
    ),
    all = FALSE,
    sort = FALSE
  )


  if (
    nrow(J) !=
      nrow(H1) ||
    nrow(J) !=
      nrow(DY)
  ) {

    stopf(
      yy,
      " 07H04与07H01 monitor-day domain不完全一致。"
    )
  }


  zero_pattern_fail <- 0L
  positive_pattern_fail <- 0L
  heat_missingness_fail <- 0L


  for (
    g in FUELS
  ) {

    for (
      dir in DIRECTIONS
    ) {

      for (
        r in RADII
      ) {

        nc <- paste0(
          "campd_",
          g,
          "_",
          dir,
          "_n_plants_",
          r,
          "km_90deg"
        )

        rh <- paste0(
          "campd_",
          g,
          "_",
          dir,
          "_heat_observed_",
          r,
          "km_90deg"
        )

        pc <- paste0(
          "campd_",
          g,
          "_",
          dir,
          "_proximity_invkm_",
          r,
          "km_90deg"
        )

        hc <- paste0(
          "campd_",
          g,
          "_",
          dir,
          "_heat_proximity_mmbtu_per_km_",
          r,
          "km_90deg"
        )


        required_now <- c(
          nc,
          rh,
          pc,
          hc
        )


        miss_now <- setdiff(
          required_now,
          names(J)
        )


        if (
          length(
            miss_now
          ) >
            0L
        ) {

          stopf(
            yy,
            " 缺少语义对照字段：",
            paste(
              miss_now,
              collapse = ", "
            )
          )
        }


        zero_pattern_fail <-
          zero_pattern_fail +
          J[
            wind_valid %in%
              TRUE &
            get(
              nc
            ) ==
              0L &
            (
              !is.finite(
                get(
                  pc
                )
              ) |
              abs(
                get(
                  pc
                )
              ) >
                TOL |
              !is.finite(
                get(
                  hc
                )
              ) |
              abs(
                get(
                  hc
                )
              ) >
                TOL
            ),
            .N
          ]


        positive_pattern_fail <-
          positive_pattern_fail +
          J[
            wind_valid %in%
              TRUE &
            get(
              nc
            ) >
              0L &
            (
              !is.finite(
                get(
                  pc
                )
              ) |
              get(
                pc
              ) <=
                0
            ),
            .N
          ]


        heat_missingness_fail <-
          heat_missingness_fail +
          J[
            wind_valid %in%
              TRUE &
            (
              is.na(
                get(
                  rh
                )
              ) !=
              is.na(
                get(
                  hc
                )
              )
            ),
            .N
          ]
      }
    }
  }


  source_mismatch <-
    PF$source_audit$n_04_missing_from_05 +
    PF$source_audit$n_05_not_in_04


  YEAR_QC_LIST[[
    length(
      YEAR_QC_LIST
    ) +
      1L
  ]] <- data.table(

    year =
      yy,

    n_monitor_days =
      nrow(DY),

    n_exposure_columns =
      length(
        exposure_cols
      ),

    duplicate_monitor_days =
      dup_day,

    invalid_wind_nonNA_cells =
      invalid_non_na,

    infinite_exposure_cells =
      inf_cells,

    proximity_radius_monotonic_fail_rows =
      proximity_monotonic_fail,

    zero_pattern_fail =
      zero_pattern_fail,

    positive_pattern_fail =
      positive_pattern_fail,

    heat_missingness_fail =
      heat_missingness_fail,

    source_layer_mismatch =
      source_mismatch
  )


  SEMANTIC_QC_LIST[[
    length(
      SEMANTIC_QC_LIST
    ) +
      1L
  ]] <- data.table(

    year =
      yy,

    zero_pattern_fail =
      zero_pattern_fail,

    positive_pattern_fail =
      positive_pattern_fail,

    heat_missingness_fail =
      heat_missingness_fail,

    wind_valid_mismatch =
      wind_mismatch
  )


  if (
    length(
      exposure_cols
    ) !=
      240L ||
    dup_day >
      0L ||
    invalid_non_na >
      0L ||
    inf_cells >
      0L ||
    proximity_monotonic_fail >
      0L ||
    zero_pattern_fail >
      0L ||
    positive_pattern_fail >
      0L ||
    heat_missingness_fail >
      0L ||
    source_mismatch >
      0L
  ) {

    stopf(
      yy,
      " 07H04 weighted exposure硬QC失败。"
    )
  }


  out_file <- file.path(
    DAILY_DIR,
    paste0(
      "07H04_fuel_daily_distance_weighted_exposure_",
      yy,
      ".csv.gz"
    )
  )


  fwrite(
    DY,
    out_file
  )


  cat(
    "monitor-days = ",
    fmt_int(
      nrow(DY)
    ),
    "\n",
    "weighted exposure cols = ",
    length(
      exposure_cols
    ),
    "\n",
    "zero-pattern mismatch = ",
    zero_pattern_fail,
    "\n",
    "heat-missingness mismatch = ",
    heat_missingness_fail,
    "\n",
    "source-layer mismatch = ",
    source_mismatch,
    "\n",
    sep = ""
  )


  rm(
    H1,
    TARGET,
    W,
    PF,
    LA,
    OUT_PARTS,
    DY,
    J
  )

  invisible(
    gc(
      verbose = FALSE
    )
  )
}


# ============================================================
# 15. Final audits
# ============================================================

MAPPING_AUDIT <- rbindlist(
  MAPPING_AUDIT_LIST,
  use.names = TRUE,
  fill = TRUE
)


SOURCE_AUDIT <- rbindlist(
  SOURCE_AUDIT_LIST,
  use.names = TRUE,
  fill = TRUE
)


YEAR_QC <- rbindlist(
  YEAR_QC_LIST,
  use.names = TRUE,
  fill = TRUE
)


SEMANTIC_QC <- rbindlist(
  SEMANTIC_QC_LIST,
  use.names = TRUE,
  fill = TRUE
)


fwrite(
  MAPPING_AUDIT,
  file.path(
    QC_DIR,
    "01_fuel_mapping_audit.csv"
  )
)


fwrite(
  SOURCE_AUDIT,
  file.path(
    QC_DIR,
    "02_source_layer_audit.csv"
  )
)


fwrite(
  YEAR_QC,
  file.path(
    QC_DIR,
    "03_weighted_exposure_qc_by_year.csv"
  )
)


fwrite(
  SEMANTIC_QC,
  file.path(
    QC_DIR,
    "04_weighted_vs_07H01_semantics_audit.csv"
  )
)


N_YEAR_FILES <- length(
  list.files(
    DAILY_DIR,
    pattern =
      "^07H04_fuel_daily_distance_weighted_exposure_[0-9]{4}\\.csv\\.gz$"
  )
)


HARD <- data.table(

  check = c(

    "07H01_hard_checks_all_pass",

    "07H03_final_freeze_pass",

    "yearly_files_24",

    "all_years_240_weighted_exposure_columns",

    "duplicate_monitor_days_zero",

    "invalid_wind_nonNA_zero",

    "infinite_exposure_cells_zero",

    "proximity_radius_monotonic_fail_zero",

    "raw_count_zero_pattern_mismatch_zero",

    "raw_count_positive_pattern_mismatch_zero",

    "raw_heat_missingness_mismatch_zero",

    "source_layer_mismatch_zero"
  ),

  observed = c(

    as.character(
      H01_HARD[
        status != "PASS",
        .N
      ]
    ),

    "1",

    as.character(
      N_YEAR_FILES
    ),

    as.character(
      sum(
        YEAR_QC$n_exposure_columns !=
          240L
      )
    ),

    as.character(
      sum(
        YEAR_QC$duplicate_monitor_days
      )
    ),

    as.character(
      sum(
        YEAR_QC$invalid_wind_nonNA_cells
      )
    ),

    as.character(
      sum(
        YEAR_QC$infinite_exposure_cells
      )
    ),

    as.character(
      sum(
        YEAR_QC$proximity_radius_monotonic_fail_rows
      )
    ),

    as.character(
      sum(
        YEAR_QC$zero_pattern_fail
      )
    ),

    as.character(
      sum(
        YEAR_QC$positive_pattern_fail
      )
    ),

    as.character(
      sum(
        YEAR_QC$heat_missingness_fail
      )
    ),

    as.character(
      sum(
        YEAR_QC$source_layer_mismatch
      )
    )
  ),

  expected = c(
    "0",
    "1",
    "24",
    "0",
    "0",
    "0",
    "0",
    "0",
    "0",
    "0",
    "0",
    "0"
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
    QC_DIR,
    "05_build_hard_checks.csv"
  )
)


cat(
  "\n============================================================\n",
  "07H04 HARD CHECKS\n",
  "============================================================\n"
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


REPORT <- c(

  "# 07H04 Fuel-specific inverse-distance-weighted daily exposure",

  "",

  paste0(
    "- Generated: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),

  "- This is a sensitivity exposure definition and does not replace frozen 07H01 Count/Heat.",

  "- Source universe: EPA CAMPD actual daily plant-days, identical fuel/source logic to frozen 07H01.",

  "- Formal fuel groups: coal, gas, oil, other_fossil, biomass, solid_waste.",

  "- Direction: frozen 90-degree upwind/downwind definition.",

  "- Cumulative radii: 20, 25, 30, 35, 40, 45, 50, 100, 150, 200 km.",

  "- Geographic distance: frozen 07A distance_km.",

  "- Radius membership uses original distance_km.",

  "- Inverse-distance denominator uses max(distance_km, 1 km).",

  "- plant proximity = sum(1 / distance_weight_km), unit 1/km.",

  "- Heat proximity = sum(actual daily Heat Input / distance_weight_km), unit MMBtu/km.",

  "- Valid wind + zero active fuel-specific plants gives both weighted exposures = 0.",

  "- Active plants with all Heat Input missing gives Heat proximity = NA.",

  "- Negative CAMPD Heat values are preserved.",

  "- Cross-group mixed and unknown fuel labels are excluded exactly as in frozen 07H01.",

  "- 07H04 exposure zero/missingness semantics were cross-checked against frozen 07H01.",

  "",

  if (
    N_FAIL ==
      0L
  ) {
    "**07H04 DISTANCE-WEIGHTED EXPOSURE BUILD: PASS**"
  } else {
    "**07H04 DISTANCE-WEIGHTED EXPOSURE BUILD: FAIL**"
  }
)


writeLines(
  REPORT,
  file.path(
    OUT,
    "07H04_距离加权逐日暴露构建报告.md"
  ),
  useBytes = TRUE
)


capture.output(
  sessionInfo(),
  file = file.path(
    OUT,
    "sessionInfo.txt"
  )
)


cat(
  "\n============================================================\n"
)


if (
  N_FAIL ==
    0L
) {

  cat(
    "07H04 DISTANCE-WEIGHTED EXPOSURE BUILD: PASS\n"
  )

} else {

  cat(
    "07H04 DISTANCE-WEIGHTED EXPOSURE BUILD: FAIL\n"
  )
}


cat(
  "============================================================\n",
  "Yearly files = ",
  N_YEAR_FILES,
  " / 24\n",
  "Distance floor = ",
  DISTANCE_FLOOR_KM,
  " km\n",
  "Hard FAIL = ",
  N_FAIL,
  "\n",
  sep = ""
)


if (
  N_FAIL >
    0L
) {

  stopf(
    "07H04存在hard failure。"
  )
}

