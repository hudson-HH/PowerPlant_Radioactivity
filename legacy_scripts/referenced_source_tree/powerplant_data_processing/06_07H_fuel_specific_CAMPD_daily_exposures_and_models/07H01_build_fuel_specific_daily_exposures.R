
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
# 07H01
# Fuel-specific CAMPD actual-daily directional exposure
#
# 不修改：
# 07A / 07C / 07D / 07E / 07F / 07G
#
# 正式暴露：
# fuel × direction × cumulative radius
#
# fuel:
# coal / gas / oil / other_fossil / biomass / solid_waste
#
# metric:
# 1. n_plants
# 2. heat_observed
#
# count:
#   daily active fuel-specific CAMPD plants; same source layer as Heat
#
# Heat:
#   EPA CAMPD actual daily Heat Input
#
# direction:
#   formal frozen 07C definition
#   wind-from vs monitor->plant spherical initial bearing
#   90° sector = +/-45°
#
# ============================================================


# ============================================================
# 1. Paths
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

F07_FREEZE <- file.path(
  F07_ROOT,
  "07F_FINAL_FREEZE.md"
)

F07_QDIR <- file.path(
  F07_ROOT,
  "quasi_daily_yearly"
)

# 正式TOP_LEVEL 07C wind，禁止搜索旧NESTED_04
WIND_DIR <- file.path(
  OLD_FROZEN,
  "07C_daily_upwind_powerplant_exposure",
  "yearly"
)

FUEL_MASTER <- file.path(
  ROOT,
  "powerplant_cleaned",
  "06C_powerplant_analysis_roles_frozen",
  "01_plant_fuel_group_year_analysis_role_master_06C_2001_2024.csv"
)

OUT <- file.path(
  PP,
  "06_07H_燃料特异CAMPD逐日暴露与模型"
)

DAILY_DIR <- file.path(
  OUT,
  "fuel_daily_yearly"
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


# ============================================================
# 2. Constants
# ============================================================

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
  20, 25, 30, 35, 40, 45, 50,
  100, 150, 200
)

DIRECTIONS <- c(
  "upwind",
  "downwind"
)

HALF_ANGLE <- 45

MAX_RADIUS <- 200

ALIGN_TOL <- 1e-7


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
    is.na(y) |
      y == "" |
      toupper(y) %in% c(
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
      xx %in% c(
        "true",
        "t",
        "1",
        "yes",
        "y"
      )
    ] <- TRUE

    z[
      xx %in% c(
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


find_exact_file <- function(
  root,
  pattern,
  label
) {

  x <- list.files(
    root,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE
  )

  if (
    length(x) !=
      1L
  ) {

    cat(
      "\n",
      label,
      " candidates:\n",
      sep = ""
    )

    print(x)

    stopf(
      "无法唯一确定",
      label,
      "。"
    )
  }

  x
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
# 4. Upstream freeze gate
# ============================================================

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
# 5. Formal frozen inputs
#
# 07A已有正式目录和旧副本同时存在。
# 禁止递归搜索、禁止自动选择。
# 只读取正式02_07A_RadNet发电站空间链接。
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


if (!dir.exists(FORMAL_07A_ROOT)) {
  stopf(
    "找不到正式07A目录：\n",
    FORMAL_07A_ROOT
  )
}


if (!file.exists(GEO_ALL_FILE)) {
  stopf(
    "找不到正式07A overall geometry：\n",
    GEO_ALL_FILE
  )
}


if (!file.exists(GEO_FUEL_FILE)) {
  stopf(
    "找不到正式07A fuel geometry：\n",
    GEO_FUEL_FILE
  )
}


cat(
  "\n正式07A输入：PASS\n",
  "overall:\n",
  GEO_ALL_FILE,
  "\n\nfuel:\n",
  GEO_FUEL_FILE,
  "\n",
  sep = ""
)


if (
  !dir.exists(
    WIND_DIR
  )
) {

  stopf(
    "找不到正式TOP_LEVEL 07C wind：",
    WIND_DIR
  )
}


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
    "TOP_LEVEL 07C wind不是24/24。"
  )
}


if (
  !file.exists(
    FUEL_MASTER
  )
) {

  stopf(
    "找不到06C canonical fuel master。"
  )
}


# ------------------------------------------------------------
# ============================================================
# Formal frozen 07B2 CAMPD plant-day × primary fuel
#
# 已经通过磁盘搜索确认正式归档位置：
# 01_发电站月度逐日源强增强流程_冻结/
#   07B2_CAMPD_daily_emissions/
#     05_plant_day_primary_fuel/
#
# 禁止递归搜索、禁止自动猜测。
# ============================================================

PF_DIR <- file.path(
  PP,
  "01_发电站月度逐日源强增强流程_冻结",
  "07B2_CAMPD_daily_emissions",
  "05_plant_day_primary_fuel"
)


if (
  !dir.exists(
    PF_DIR
  )
) {

  stopf(
    "找不到正式冻结07B2 primary-fuel目录：\n",
    PF_DIR
  )
}


PF_FILES_CHECK <- list.files(
  PF_DIR,
  pattern =
    "^CAMPD_plant_day_primary_fuel_[0-9]{4}_q[1-4]\\.csv\\.gz$",
  full.names = TRUE
)


PF_BASENAME <- basename(
  PF_FILES_CHECK
)


PF_TAG <- sub(
  "^CAMPD_plant_day_primary_fuel_",
  "",
  PF_BASENAME
)


PF_TAG <- sub(
  "\\.csv\\.gz$",
  "",
  PF_TAG
)


EXPECTED_PF_TAG <- sort(
  unlist(
    lapply(
      2001:2024,
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


PF_TAG_SORTED <- sort(
  PF_TAG
)


MISSING_PF <- setdiff(
  EXPECTED_PF_TAG,
  PF_TAG_SORTED
)


EXTRA_PF <- setdiff(
  PF_TAG_SORTED,
  EXPECTED_PF_TAG
)


DUP_PF <- unique(
  PF_TAG_SORTED[
    duplicated(
      PF_TAG_SORTED
    )
  ]
)


cat(
  "\n============================================================\n",
  "正式冻结07B2 primary-fuel输入检查\n",
  "============================================================\n",
  "目录：\n",
  PF_DIR,
  "\n季度文件数 = ",
  length(
    PF_FILES_CHECK
  ),
  "\n缺失季度数 = ",
  length(
    MISSING_PF
  ),
  "\n额外季度数 = ",
  length(
    EXTRA_PF
  ),
  "\n重复季度数 = ",
  length(
    DUP_PF
  ),
  "\n",
  sep = ""
)


if (
  length(
    MISSING_PF
  ) > 0L
) {

  cat("\n缺失季度：\n")
  print(MISSING_PF)
}


if (
  length(
    EXTRA_PF
  ) > 0L
) {

  cat("\n额外季度：\n")
  print(EXTRA_PF)
}


if (
  length(
    DUP_PF
  ) > 0L
) {

  cat("\n重复季度：\n")
  print(DUP_PF)
}


if (
  length(
    PF_FILES_CHECK
  ) != 96L ||
  length(
    MISSING_PF
  ) != 0L ||
  length(
    EXTRA_PF
  ) != 0L ||
  length(
    DUP_PF
  ) != 0L
) {

  stopf(
    "正式冻结07B2 primary-fuel文件结构不是完整的2001–2024 × 4季度。"
  )
}


cat(
  "\n07B2 frozen primary-fuel：PASS\n",
  "2001–2024 × 4季度 = 96/96\n",
  sep = ""
)


cat(
  "\n07A overall geometry:\n",
  GEO_ALL_FILE,
  "\n\n07A fuel geometry:\n",
  GEO_FUEL_FILE,
  "\n\n07B2 plant-day x primary fuel:\n",
  PF_DIR,
  "\n\n正式07C wind:\n",
  WIND_DIR,
  "\n",
  sep = ""
)


# ============================================================
# 6. Build loc_num <-> RadNet ID crosswalk from frozen 07F
# ============================================================

cross_parts <- vector(
  "list",
  length(YEARS)
)


for (
  ii in seq_along(
    YEARS
  )
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
    !file.exists(
      qf
    )
  ) {

    stopf(
      "缺少07F year：",
      yy
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
    "07F loc_num <-> radnet_monitor_id不是一一对应。"
  )
}


# ============================================================
# 7. Geometry standardization
# ============================================================

standardize_geometry <- function(
  file,
  fuel_specific = FALSE
) {

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


  fuel_col <- if (
    fuel_specific
  ) {

    pick_name(
      H,
      c(
        "fuel_group",
        "plant_fuel_group"
      ),
      "fuel_group"
    )

  } else {

    NA_character_
  }


  need <- unique(
    na.omit(
      c(
        year_col,
        plant_col,
        monitor_col,
        loc_col,
        dist_col,
        bearing_col,
        fuel_col
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

  if (
    fuel_specific &&
    fuel_col !=
      "fuel_group"
  ) {
    setnames(
      D,
      fuel_col,
      "fuel_group"
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


  if (
    fuel_specific
  ) {

    D[
      ,
      fuel_group :=
        tolower(
          trimws(
            as.character(
              fuel_group
            )
          )
        )
    ]

    D <- D[
      fuel_group %in%
        FUELS
    ]

    D <- unique(
      D[
        ,
        .(
          loc_num,
          radnet_monitor_id,
          year,
          plant_id,
          fuel_group,
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
        plant_id,
        fuel_group
      )
    ][
      N >
        1L
    ]

  } else {

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
  }


  if (
    nrow(
      dup
    ) >
      0L
  ) {

    stopf(
      "07A geometry存在重复主键。"
    )
  }


  D
}


cat(
  "\n读取冻结07A geometry...\n"
)


GEO_ALL <- standardize_geometry(
  GEO_ALL_FILE,
  FALSE
)


GEO_FUEL <- standardize_geometry(
  GEO_FUEL_FILE,
  TRUE
)


cat(
  "overall links = ",
  fmt_int(
    nrow(
      GEO_ALL
    )
  ),
  "\n",
  "fuel links = ",
  fmt_int(
    nrow(
      GEO_FUEL
    )
  ),
  "\n",
  sep = ""
)


# ============================================================

# ============================================================
# 8. Canonical fuel labels
# ============================================================
#
# IMPORTANT:
#
# 07B2 preserves raw EPA CAMPD primaryFuelInfo.
#
# The daily CAMPD fuel label is the formal fuel source for
# BOTH:
#
#   1. fuel-specific active plant count
#   2. fuel-specific actual Heat Input
#
# 06C annual fuel groups are NOT used as a hard compatibility
# gate for daily CAMPD fuel.
#
# Mixed labels are handled without double counting:
#
#   same-group combination:
#     Natural Gas, Pipeline Natural Gas -> gas
#
#   cross-group combination:
#     Coal, Pipeline Natural Gas -> mixed
#
# Cross-group "mixed" rows are preserved in the mapping audit
# and all-source Heat reconstruction, but are NOT duplicated
# into coal and gas.
#
# Blank labels are preserved as unknown and are not guessed.
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

# 9. Read formal wind
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
      "正式07C wind缺少：",
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


  W[
    ,
    month :=
      as.integer(
        format(
          date,
          "%m"
        )
      )
  ]


  W
}


# ============================================================

# ============================================================

# ============================================================
# 10. Read formal 07B2 source layers + apply frozen 07B3 source resolution
# ============================================================
#
# IMPORTANT
#
# corrected 07D source hierarchy is defined at plant × year × month
# by frozen 07B3 preferred_source_resolution.
#
# Therefore 07H uses exactly the same eligible CAMPD plant-month
# universe.
#
# Two distinct 07B2 layers are intentionally used:
#
#   04_plant_day
#       -> overall CAMPD plant-day source
#       -> used ONLY for direct reconstruction of corrected 07D
#
#   05_plant_day_primary_fuel
#       -> plant-day × raw EPA primaryFuelInfo
#       -> used for formal fuel-specific Count + Heat
#
# Do NOT reconstruct the overall benchmark by summing 05.
# ============================================================


PLANT_DIR <- file.path(
  PP,
  "01_发电站月度逐日源强增强流程_冻结",
  "07B2_CAMPD_daily_emissions",
  "04_plant_day"
)


B07B3_HARM_FILE <- file.path(
  PP,
  "01_发电站月度逐日源强增强流程_冻结",
  "07B3_EIA_CAMPD_harmonization",
  "05_EIA_CAMPD_plant_month_harmonization_2001_2024.csv.gz"
)


if (!dir.exists(PLANT_DIR)) {

  stopf(
    "找不到正式冻结07B2 04_plant_day：\n",
    PLANT_DIR
  )
}


if (!file.exists(B07B3_HARM_FILE)) {

  stopf(
    "找不到正式冻结07B3 harmonization：\n",
    B07B3_HARM_FILE
  )
}


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


PREFERRED_CAMPD_VALUES <- sort(
  unique(
    HARM_SOURCE[
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
    PREFERRED_CAMPD_VALUES
  ) ==
    0L
) {

  stopf(
    "冻结07B3中没有识别到preferred CAMPD source。"
  )
}


cat(
  "\n冻结07B3 preferred CAMPD source：\n"
)

print(
  PREFERRED_CAMPD_VALUES
)


CAMPD_ELIGIBLE_PM <- unique(
  HARM_SOURCE[
    preferred_source_resolution %in%
      PREFERRED_CAMPD_VALUES,
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
    "07B3 preferred CAMPD plant-month主键重复。"
  )
}


MAPPING_AUDIT_LIST <- list()
MAPPING_AUDIT_I <- 0L

SOURCE_LAYER_AUDIT_LIST <- list()
SOURCE_LAYER_AUDIT_I <- 0L


read_pf_year <- function(yy) {

  # ==========================================================
  # A. Read 05 plant-day × primary fuel
  # ==========================================================

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


  if (
    !all(
      file.exists(
        pf_files
      )
    )
  ) {

    stopf(
      yy,
      "年07B2 05 primary-fuel不是4/4季度。"
    )
  }


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


  # ==========================================================
  # B. Read 04 overall plant-day
  # ==========================================================

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


  if (
    !all(
      file.exists(
        plant_files
      )
    )
  ) {

    stopf(
      yy,
      "年07B2 04 plant-day不是4/4季度。"
    )
  }


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


  # ==========================================================
  # C. Apply EXACT frozen 07B3 preferred CAMPD plant-month pool
  # ==========================================================

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


  # ==========================================================
  # D. Restrict to frozen 07A overall source geography
  # ==========================================================

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


  # ==========================================================
  # E. Hard duplicate check for 04 overall plant-day
  # ==========================================================

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

    fwrite(
      dup_plant_day,
      file.path(
        QC_DIR,
        paste0(
          "DUPLICATE_04_PLANT_DAY_",
          yy,
          ".csv"
        )
      )
    )


    stopf(
      yy,
      " 07B2 04 plant-day存在重复主键。"
    )
  }


  # ==========================================================
  # F. Map raw CAMPD fuel labels from 05
  # ==========================================================

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


  # ==========================================================
  # G. Mapping audit
  # ==========================================================

  AUDIT <- D[
    ,
    .(
      n_plant_days =
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


  AUDIT[
    ,
    included_in_formal_fuel_analysis :=
      fuel_group %in%
        FUELS
  ]


  MAPPING_AUDIT_I <<-
    MAPPING_AUDIT_I +
    1L


  MAPPING_AUDIT_LIST[[MAPPING_AUDIT_I]] <<-
    AUDIT


  # ==========================================================
  # H. 04 vs 05 plant-day source-layer coverage audit
  # ==========================================================

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
      ),

    n_04_heat_nonmissing =
      P[
        is.finite(
          heat
        ),
        .N
      ],

    n_04_heat_zero =
      P[
        is.finite(
          heat
        ) &
          heat ==
            0,
        .N
      ]
  )


  SOURCE_LAYER_AUDIT_I <<-
    SOURCE_LAYER_AUDIT_I +
    1L


  SOURCE_LAYER_AUDIT_LIST[[SOURCE_LAYER_AUDIT_I]] <<-
    source_audit


  cat(
    "07B2 source-layer ",
    yy,
    ": 04 plant-days=",
    source_audit$n_04_plant_days,
    "; 05 plant-days=",
    source_audit$n_05_plant_days,
    "; 04 missing from 05=",
    source_audit$n_04_missing_from_05,
    "\n",
    sep = ""
  )


  # ==========================================================
  # I. Overall all-source Heat:
  #    use 04 DIRECTLY, exactly matching corrected 07D source
  # ==========================================================

  ALL_HEAT <- P[
    ,
    .(
      plant_id,
      date,
      heat
    )
  ]


  # ==========================================================
  # J. Formal fuel-specific daily source:
  #    use 05
  #
  # Count and Heat come from the same plant×date×fuel table.
  # ==========================================================

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


  # ==========================================================
  # K. Informative fuel coverage
  # ==========================================================

  total_04_heat <- P[
    is.finite(
      heat
    ),
    sum(
      heat
    )
  ]


  formal_05_heat <- D[
    fuel_group %in%
      FUELS &
      is.finite(
        heat
      ),
    sum(
      heat
    )
  ]


  mixed_05_heat <- D[
    fuel_group ==
      "mixed" &
      is.finite(
        heat
      ),
    sum(
      heat
    )
  ]


  unknown_05_heat <- D[
    fuel_group ==
      "unknown" &
      is.finite(
        heat
      ),
    sum(
      heat
    )
  ]


  cat(
    "CAMPD fuel mapping ",
    yy,
    ": formal rows=",
    D[
      fuel_group %in%
        FUELS,
      .N
    ],
    "; mixed rows=",
    D[
      fuel_group ==
        "mixed",
      .N
    ],
    "; unknown rows=",
    D[
      fuel_group ==
        "unknown",
      .N
    ],
    "\n",
    sep = ""
  )


  cat(
    "Heat MMBtu: 04 overall=",
    format(
      total_04_heat,
      scientific = FALSE,
      trim = TRUE
    ),
    "; 05 formal fuel=",
    format(
      formal_05_heat,
      scientific = FALSE,
      trim = TRUE
    ),
    "; mixed excluded=",
    format(
      mixed_05_heat,
      scientific = FALSE,
      trim = TRUE
    ),
    "; unknown excluded=",
    format(
      unknown_05_heat,
      scientific = FALSE,
      trim = TRUE
    ),
    "\n",
    sep = ""
  )


  list(
    all_heat =
      ALL_HEAT,

    fuel_heat =
      FUEL_ACTIVITY
  )
}

# 11. Build one monitor-year
# ============================================================
#
# Fuel-specific count and Heat now use the SAME source:
#
#   CAMPD daily primaryFuelInfo
#          ↓
#   formal daily fuel crosswalk
#          ↓
#   plant × date × fuel_group
#          ↓
#   frozen 07A OVERALL geometry
#          ↓
#   frozen 07C daily wind
#
# Thus:
#
#   n_plants = number of daily CAMPD plants in that fuel group
#   Heat     = actual daily Heat Input from the same plants
#
# links_fuel is retained only for call compatibility and is not
# used for exposure construction.
# ============================================================


build_monitor_year <- function(
  loc,
  yy,
  W,
  links_all,
  links_fuel,
  all_heat,
  fuel_heat,
  frozen_total
) {

  # Retained only so the existing main loop need not be changed.
  invisible(
    links_fuel
  )


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


  # ----------------------------------------------------------
  # Pre-create 240 exposure fields:
  #
  # 6 fuels × 2 directions × 10 radii × 2 metrics
  # ----------------------------------------------------------

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


        hc <- paste0(
          "campd_",
          g,
          "_",
          dir,
          "_heat_observed_",
          r,
          "km_90deg"
        )


        D[
          ,
          (nc) :=
            NA_integer_
        ]


        # Within the observed CAMPD source layer:
        # valid wind + no matching daily fuel-specific plant = 0.
        D[
          wind_valid %in%
            TRUE,
          (nc) :=
            0L
        ]


        D[
          ,
          (hc) :=
            NA_real_
        ]
      }
    }
  }


  # ----------------------------------------------------------
  # Frozen 07A OVERALL monitor-plant geometry.
  #
  # Fuel identity is NOT taken from 07A/06C here.
  # ----------------------------------------------------------

  LALL <- links_all[
    loc_num ==
      loc
  ]


  # ----------------------------------------------------------
  # Join formal daily fuel activity to geometry + wind.
  #
  # HH is the common source for BOTH count and Heat.
  # ----------------------------------------------------------

  if (
    nrow(
      LALL
    ) >
      0L &&
    nrow(
      fuel_heat
    ) >
      0L
  ) {

    plant_ids <- unique(
      LALL$plant_id
    )


    FH <- fuel_heat[
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


  # ----------------------------------------------------------
  # Fuel-specific daily count + Heat.
  #
  # SAME HH source table for both metrics.
  # ----------------------------------------------------------

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

        nc <- paste0(
          "campd_",
          g,
          "_",
          dir,
          "_n_plants_",
          r,
          "km_90deg"
        )


        hc <- paste0(
          "campd_",
          g,
          "_",
          dir,
          "_heat_observed_",
          r,
          "km_90deg"
        )


        if (
          nrow(
            HH
          ) >
            0L
        ) {

          # --------------------------------------------------
          # Daily active fuel-specific plant count
          # --------------------------------------------------

          AA <- HH[
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
                uniqueN(
                  plant_id
                )
            ),
            by = date
          ]


          if (
            nrow(
              AA
            ) >
              0L
          ) {

            D[
              AA,
              on = "date",
              (nc) :=
                i.value
            ]
          }
        }


        # ----------------------------------------------------
        # No daily plant in this fuel/radius/direction:
        # observed-source-layer zero.
        # ----------------------------------------------------

        D[
          wind_valid %in%
            TRUE &
            get(
              nc
            ) ==
              0L,
          (hc) :=
            0
        ]


        # ----------------------------------------------------
        # Actual daily Heat from the exact same plant-day rows.
        #
        # If count >0 but every Heat is missing:
        # sum_keep_all_na returns NA.
        #
        # If at least one plant reports Heat:
        # sum available Heat values.
        #
        # Negative CAMPD values are preserved.
        # ----------------------------------------------------

        if (
          nrow(
            HH
          ) >
            0L
        ) {

          HA <- HH[
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
                  heat
                )
            ),
            by = date
          ]


          if (
            nrow(
              HA
            ) >
              0L
          ) {

            D[
              HA,
              on = "date",
              (hc) :=
                i.value
            ]
          }
        }
      }
    }
  }


  # ==========================================================
  # Direct source-pool alignment
  #
  # Reconstruct ALL-fuel 20-km Heat from raw 07B2 and compare
  # with frozen corrected 07D/07F.
  #
  # Mixed and unknown fuel labels remain INCLUDED here.
  #
  # Therefore fuel-specific exclusions cannot artificially
  # improve this alignment.
  # ==========================================================

  if (
    nrow(
      LALL
    ) >
      0L
  ) {

    AH <- all_heat[
      plant_id %in%
        unique(
          LALL$plant_id
        )
    ]


  } else {

    AH <- data.table()
  }


  if (
    nrow(
      AH
    ) >
      0L
  ) {

    X <- merge(
      AH,
      LALL[
        ,
        .(
          plant_id,
          distance_km,
          bearing_deg
        )
      ],
      by = "plant_id",
      all = FALSE,
      allow.cartesian = TRUE,
      sort = FALSE
    )


    X <- merge(
      X,
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


    X <- X[
      wind_valid %in%
        TRUE &
        distance_km <=
          20 +
          1e-8
    ]


    X[
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

    X <- data.table()
  }


  align <- list()


  for (
    dir in DIRECTIONS
  ) {

    frozen_n_col <- paste0(
      "campd_",
      dir,
      "_n_plants_20km_90deg"
    )


    frozen_h_col <- paste0(
      "campd_",
      dir,
      "_heat_observed_20km_90deg"
    )


    flag_col <- paste0(
      dir,
      "_flag"
    )


    AA <- if (
      nrow(
        X
      ) >
        0L
    ) {

      X[
        get(
          flag_col
        ) %in%
          TRUE,
        .(
          reconstructed =
            sum_keep_all_na(
              heat
            )
        ),
        by = date
      ]


    } else {

      data.table(
        date =
          as.IDate(
            character()
          ),

        reconstructed =
          numeric()
      )
    }


    Z <- frozen_total[
      loc_num ==
        loc,
      .(
        date,

        frozen_n =
          get(
            frozen_n_col
          ),

        frozen_heat =
          get(
            frozen_h_col
          )
      )
    ]


    Z[
      ,
      reconstructed :=
        NA_real_
    ]


    # Frozen total CAMPD count == 0 is the exact zero rule used
    # by corrected 07D.
    Z[
      is.finite(
        frozen_n
      ) &
        frozen_n ==
          0,
      reconstructed :=
        0
    ]


    if (
      nrow(
        AA
      ) >
        0L
    ) {

      Z[
        AA,
        on = "date",
        reconstructed :=
          i.reconstructed
      ]
    }


    Z[
      ,
      mismatch :=
        fifelse(

          is.na(
            frozen_heat
          ) &
            is.na(
              reconstructed
            ),

          FALSE,

          fifelse(

            is.finite(
              frozen_heat
            ) &
              is.finite(
                reconstructed
              ),

            abs(
              frozen_heat -
                reconstructed
            ) >
              ALIGN_TOL,

            TRUE
          )
        )
    ]


    Z[
      ,
      `:=`(
        loc_num =
          loc,

        direction =
          dir
      )
    ]


    align[[dir]] <- Z
  }


  list(

    exposure =
      D,

    alignment =
      rbindlist(
        align,
        use.names = TRUE,
        fill = TRUE
      )
  )
}

# 12. Main yearly loop
# ============================================================

ALIGN_SUM_LIST <- list()
YEAR_QC_LIST <- list()

align_i <- 0L
qc_i <- 0L


for (
  yy in YEARS
) {

  cat(
    "\n============================================================\n",
    "07H01 YEAR ",
    yy,
    "\n",
    "============================================================\n",
    sep = ""
  )


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


  frozen_cols <- c(
    "loc_num",
    "date",
    "radnet_monitor_id",
    "campd_upwind_n_plants_20km_90deg",
    "campd_upwind_heat_observed_20km_90deg",
    "campd_downwind_n_plants_20km_90deg",
    "campd_downwind_heat_observed_20km_90deg"
  )


  miss <- setdiff(
    frozen_cols,
    qh
  )


  if (
    length(miss) >
      0L
  ) {

    stopf(
      yy,
      " 07F缺少benchmark字段：",
      paste(
        miss,
        collapse = ", "
      )
    )
  }


  QT <- fread(
    qf,
    select = frozen_cols,
    showProgress = FALSE
  )


  QT[
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


  # 170 overlapping quasi rows are expected.
  # Collapse to unique monitor-day for exposure construction.
  value_cols <- setdiff(
    frozen_cols,
    c(
      "loc_num",
      "date",
      "radnet_monitor_id"
    )
  )


  # ==========================================================
  # Check duplicated quasi rows within the same monitor-day.
  #
  # 07F has 170 intentional extra quasi rows because one
  # monitor-day may belong to overlapping filter periods.
  # The frozen 07D exposure values must nevertheless be
  # identical within loc_num + date.
  #
  # Use explicit per-column checks instead of rowSums(.SD > 1),
  # which is version-sensitive in data.table.
  # ==========================================================

  consistency_check <- QT[
    ,
    lapply(
      .SD,
      function(x) {

        y <- x[
          !is.na(x)
        ]

        if (
          length(y) == 0L
        ) {

          0L

        } else {

          data.table::uniqueN(
            y
          )
        }
      }
    ),
    by = .(
      loc_num,
      date
    ),
    .SDcols = value_cols
  ]


  bad_consistency_row <- rep(
    FALSE,
    nrow(
      consistency_check
    )
  )


  for (
    cc in value_cols
  ) {

    bad_consistency_row <-
      bad_consistency_row |
      consistency_check[[cc]] > 1L
  }


  consistency_bad <- sum(
    bad_consistency_row
  )


  if (
    consistency_bad > 0L
  ) {

    bad_keys <- consistency_check[
      bad_consistency_row,
      .(
        loc_num,
        date
      )
    ]


    fwrite(
      bad_keys,
      file.path(
        QC_DIR,
        paste0(
          "FROZEN_07D_DUPLICATE_VALUE_MISMATCH_",
          yy,
          ".csv"
        )
      )
    )


    stopf(
      yy,
      " overlapping quasi rows的冻结07D值不一致：",
      consistency_bad,
      " monitor-days。"
    )
  }


  # ==========================================================
  # Collapse intentional duplicate quasi rows to one
  # monitor-day record.
  #
  # Because the consistency check above has already proved
  # frozen 07D values are identical within monitor-day,
  # taking the first record is exact and safe.
  # ==========================================================

  FROZEN_DAY <- unique(
    QT,
    by = c(
      "loc_num",
      "date"
    )
  )


  expected_unique_days <- uniqueN(
    QT,
    by = c(
      "loc_num",
      "date"
    )
  )


  if (
    nrow(FROZEN_DAY) !=
      expected_unique_days
  ) {

    stopf(
      yy,
      " FROZEN_DAY collapse行数异常。"
    )
  }


  TARGET <- FROZEN_DAY[
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
    nrow(
      W
    ) !=
      nrow(
        TARGET
      )
  ) {

    stopf(
      yy,
      " wind目标域行数不一致。"
    )
  }


  PF <- read_pf_year(
    yy
  )


  LA <- GEO_ALL[
    year ==
      yy
  ]

  LF <- GEO_FUEL[
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
    length(
      locs
    )
  )


  ALIGN_PARTS <- vector(
    "list",
    length(
      locs
    )
  )


  for (
    jj in seq_along(
      locs
    )
  ) {

    loc <- locs[jj]


    z <- build_monitor_year(

      loc =
        loc,

      yy =
        yy,

      W =
        W,

      links_all =
        LA,

      links_fuel =
        LF,

      all_heat =
        PF$all_heat,

      fuel_heat =
        PF$fuel_heat,

      frozen_total =
        FROZEN_DAY
    )


    OUT_PARTS[[jj]] <- z$exposure

    ALIGN_PARTS[[jj]] <- z$alignment
  }


  DY <- rbindlist(
    OUT_PARTS,
    use.names = TRUE,
    fill = TRUE
  )


  AL <- rbindlist(
    ALIGN_PARTS,
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


  exp_cols <- grep(
    "^campd_(coal|gas|oil|other_fossil|biomass|solid_waste)_",
    names(
      DY
    ),
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
      cc in exp_cols
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


  # Count radius monotonicity
  count_monotonic_fail <- 0L

  for (
    g in FUELS
  ) {

    for (
      dir in DIRECTIONS
    ) {

      cc <- paste0(
        "campd_",
        g,
        "_",
        dir,
        "_n_plants_",
        RADII,
        "km_90deg"
      )


      M <- as.matrix(
        DY[
          wind_valid %in%
            TRUE,
          ..cc
        ]
      )


      if (
        nrow(
          M
        ) >
          0L
      ) {

        count_monotonic_fail <-
          count_monotonic_fail +
          sum(
            apply(
              M,
              1,
              function(x) {
                any(
                  diff(x) <
                    0,
                  na.rm = TRUE
                )
              }
            )
          )
      }
    }
  }


  zero_heat_fail <- 0L

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

        hc <- paste0(
          "campd_",
          g,
          "_",
          dir,
          "_heat_observed_",
          r,
          "km_90deg"
        )


        zero_heat_fail <-
          zero_heat_fail +
          DY[
            wind_valid %in%
              TRUE &
              get(
                nc
              ) ==
                0L &
              (
                !is.finite(
                  get(
                    hc
                  )
                ) |
                get(
                  hc
                ) !=
                  0
              ),
            .N
          ]
      }
    }
  }


  alignment_fail <- AL[
    mismatch %in%
      TRUE,
    .N
  ]


  align_i <- align_i + 1L

  ALIGN_SUM_LIST[[align_i]] <- AL[
    ,
    .(
      n_rows =
        .N,

      n_frozen_finite =
        sum(
          is.finite(
            frozen_heat
          )
        ),

      n_reconstructed_finite =
        sum(
          is.finite(
            reconstructed
          )
        ),

      n_mismatch =
        sum(
          mismatch %in%
            TRUE
        ),

      max_abs_diff =
        if (
          any(
            is.finite(
              frozen_heat
            ) &
              is.finite(
                reconstructed
              )
          )
        ) {

          max(
            abs(
              frozen_heat -
                reconstructed
            ),
            na.rm = TRUE
          )

        } else {

          NA_real_
        }
    ),
    by = direction
  ][
    ,
    year :=
      yy
  ]


  qc_i <- qc_i + 1L

  YEAR_QC_LIST[[qc_i]] <- data.table(

    year =
      yy,

    n_monitor_days =
      nrow(
        DY
      ),

    n_exposure_columns =
      length(
        exp_cols
      ),

    duplicate_monitor_days =
      dup_day,

    invalid_wind_nonNA_cells =
      invalid_non_na,

    count_radius_monotonic_fail_rows =
      count_monotonic_fail,

    zero_plant_heat_rule_fail =
      zero_heat_fail,

    source_pool_alignment_fail =
      alignment_fail
  )


  if (
    alignment_fail >
      0L
  ) {

    bad <- AL[
      mismatch %in%
        TRUE
    ]


    fwrite(
      bad,
      file.path(
        QC_DIR,
        paste0(
          "SOURCE_POOL_ALIGNMENT_MISMATCH_",
          yy,
          ".csv.gz"
        )
      )
    )


    fwrite(
      rbindlist(
        ALIGN_SUM_LIST,
        use.names = TRUE,
        fill = TRUE
      ),
      file.path(
        QC_DIR,
        "02_source_pool_alignment_by_year.csv"
      )
    )


    stopf(
      yy,
      " fuel-specific source构建与冻结corrected 07D的20 km Heat不一致：",
      alignment_fail,
      " monitor-day-direction rows。",
      "\n已停止，禁止继续建模。"
    )
  }


  if (
    dup_day >
      0L ||
    invalid_non_na >
      0L ||
    count_monotonic_fail >
      0L ||
    zero_heat_fail >
      0L
  ) {

    stopf(
      yy,
      " fuel-specific exposure硬QC失败。"
    )
  }


  out_file <- file.path(
    DAILY_DIR,
    paste0(
      "07H_fuel_daily_exposure_",
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
      nrow(
        DY
      )
    ),
    "\nexposure cols = ",
    length(
      exp_cols
    ),
    "\nsource alignment mismatch = ",
    alignment_fail,
    "\n",
    sep = ""
  )


  rm(
    QT,
    FROZEN_DAY,
    TARGET,
    W,
    PF,
    LA,
    LF,
    DY,
    AL,
    OUT_PARTS,
    ALIGN_PARTS
  )

  invisible(
    gc(
      verbose = FALSE
    )
  )
}


# ============================================================
# 13. Final audit outputs
# ============================================================

MAPPING_AUDIT <- rbindlist(
  MAPPING_AUDIT_LIST,
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


ALIGN_SUM <- rbindlist(
  ALIGN_SUM_LIST,
  use.names = TRUE,
  fill = TRUE
)


fwrite(
  ALIGN_SUM,
  file.path(
    QC_DIR,
    "02_source_pool_alignment_by_year.csv"
  )
)


YEAR_QC <- rbindlist(
  YEAR_QC_LIST,
  use.names = TRUE,
  fill = TRUE
)


fwrite(
  YEAR_QC,
  file.path(
    QC_DIR,
    "03_exposure_qc_by_year.csv"
  )
)


HARD <- data.table(

  check = c(
    "yearly_files_24",
    "all_years_240_exposure_columns",
    "duplicate_monitor_days_zero",
    "invalid_wind_nonNA_zero",
    "count_radius_monotonic_fail_zero",
    "zero_plant_heat_rule_fail_zero",
    "source_pool_alignment_fail_zero"
  ),

  observed = c(

    length(
      list.files(
        DAILY_DIR,
        pattern =
          "^07H_fuel_daily_exposure_[0-9]{4}\\.csv\\.gz$"
      )
    ),

    sum(
      YEAR_QC$n_exposure_columns !=
        240L
    ),

    sum(
      YEAR_QC$duplicate_monitor_days
    ),

    sum(
      YEAR_QC$invalid_wind_nonNA_cells
    ),

    sum(
      YEAR_QC$count_radius_monotonic_fail_rows
    ),

    sum(
      YEAR_QC$zero_plant_heat_rule_fail
    ),

    sum(
      YEAR_QC$source_pool_alignment_fail
    )
  ),

  expected = c(
    "24",
    "0 bad years",
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
    c(
      ifelse(
        observed[1] ==
          24,
        "PASS",
        "FAIL"
      ),

      ifelse(
        observed[2] ==
          0,
        "PASS",
        "FAIL"
      ),

      ifelse(
        observed[3] ==
          0,
        "PASS",
        "FAIL"
      ),

      ifelse(
        observed[4] ==
          0,
        "PASS",
        "FAIL"
      ),

      ifelse(
        observed[5] ==
          0,
        "PASS",
        "FAIL"
      ),

      ifelse(
        observed[6] ==
          0,
        "PASS",
        "FAIL"
      ),

      ifelse(
        observed[7] ==
          0,
        "PASS",
        "FAIL"
      )
    )
]


fwrite(
  HARD,
  file.path(
    QC_DIR,
    "04_build_hard_checks.csv"
  )
)


print(
  HARD
)


N_FAIL <- HARD[
  status ==
    "FAIL",
  .N
]


cat(
  "\n============================================================\n"
)


if (
  N_FAIL ==
    0L
) {

  cat(
    "07H01 FUEL-SPECIFIC EXPOSURE BUILD: PASS\n"
  )

} else {

  cat(
    "07H01 FUEL-SPECIFIC EXPOSURE BUILD: FAIL\n"
  )
}


cat(
  "============================================================\n",
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
    "07H01存在hard failure。"
  )
}


REPORT <- c(

  "# 07H01 Fuel-specific CAMPD actual-daily exposure",

  "",

  "- Status: PASS",

  "- Frozen upstream files were read-only.",

  "- Fuel groups: coal, gas, oil, other_fossil, biomass, solid_waste.
- Cross-group CAMPD mixed-fuel labels are audited but excluded from formal fuel-specific exposure to prevent double counting.
- Blank CAMPD fuel labels are retained as unknown in the audit and are not guessed.",

  "- Cumulative radii: 20, 25, 30, 35, 40, 45, 50, 100, 150, 200 km.",

  "- Direction: 90-degree upwind and downwind sectors.",

  "- Plant count: daily active fuel-specific CAMPD plants from the same source layer as Heat Input.",

  "- Heat Input: actual daily EPA CAMPD Heat Input.",

  "- Invalid wind/geography days remain NA.",

  "- Valid wind + zero fuel-specific plants gives count=0 and Heat=0.",

  "- Count>0 but no observed CAMPD Heat remains Heat=NA.",

  "- Negative CAMPD Heat values are preserved.",

  "- Reconstructed all-fuel 20-km Heat was directly cross-checked against frozen corrected 07D/07F.",

  "",

  "**07H01 FUEL-SPECIFIC EXPOSURE BUILD: PASS**"
)


writeLines(
  REPORT,
  file.path(
    OUT,
    "07H01_燃料特异逐日暴露构建报告.md"
  ),
  useBytes = TRUE
)

