#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999,
  width = 220
)

suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 07H15A0
#
# Purpose:
#   基于已经验证的07H14 360°逐日文件，
#   补建：
#
#     other_fossil
#     solid_waste
#
#   的：
#     Count
#     actual daily Heat Input
#
#   radii:
#     20,25,30,35,40,45,50 km
#
# 原07H14结果不修改、不覆盖。
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

OLD_ND_DIR <- file.path(
  HROOT,
  "07H14_模型阶梯诊断_论文式随机截距到站点年份固定效应",
  "non_directional_daily_yearly"
)

OUT <- file.path(
  HROOT,
  "07H15A_扩展燃料360逐日暴露"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

GEO_FILE <- file.path(
  PP,
  "00_既往发电站空间暴露与模型流程_冻结",
  "02_07A_RadNet发电站空间链接",
  "07A_radnet_powerplant_same_year_250km",
  "03_monitor_year_all_plant_link_250km_07A_2001_2024.csv"
)

HARM_FILE <- file.path(
  PP,
  "01_发电站月度逐日源强增强流程_冻结",
  "07B3_EIA_CAMPD_harmonization",
  "05_EIA_CAMPD_plant_month_harmonization_2001_2024.csv.gz"
)

PF_DIR <- file.path(
  PP,
  "01_发电站月度逐日源强增强流程_冻结",
  "07B2_CAMPD_daily_emissions",
  "05_plant_day_primary_fuel"
)

QDIR <- file.path(
  PP,
  "04_07F_正式统计分析数据构建",
  "quasi_daily_yearly"
)


# ============================================================
# 2. Parameters
# ============================================================

YEARS <- 2001:2024

EXTRA_FUELS <- c(
  "other_fossil",
  "solid_waste"
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

MAX_RADIUS <- 50


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
        "字段。"
      )
    }

    return(
      NA_character_
    )
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


# ============================================================
# 4. Input audit
# ============================================================

cat(
  "============================================================\n",
  "07H15A0 INPUT AUDIT\n",
  "============================================================\n",
  sep = ""
)

if (
  !file.exists(GEO_FILE)
) {
  stopf(
    "缺少07A geometry。"
  )
}

if (
  !file.exists(HARM_FILE)
) {
  stopf(
    "缺少07B3 harmonization。"
  )
}

OLD_FILES <- file.path(
  OLD_ND_DIR,
  paste0(
    "07H14_nondirectional_fuel_daily_exposure_",
    YEARS,
    ".csv.gz"
  )
)

if (
  sum(
    file.exists(
      OLD_FILES
    )
  ) != 24L
) {
  stopf(
    "07H14 360°文件不是24/24。"
  )
}

PF_FILES <- unlist(
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
      PF_FILES
    )
  ) != 96L
) {
  stopf(
    "07B2季度primary-fuel文件不是96/96。"
  )
}

cat(
  "07A geometry：PASS\n",
  "07B2 primary fuel：96/96 PASS\n",
  "07B3 harmonization：PASS\n",
  "07H14 original 360°：24/24 PASS\n",
  sep = ""
)


# ============================================================
# 5. loc_num <-> radnet_monitor_id crosswalk
# ============================================================

CROSS_LIST <- vector(
  "list",
  length(YEARS)
)

for (
  ii in seq_along(YEARS)
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

  z <- fread(
    f,
    select = c(
      "loc_num",
      "radnet_monitor_id"
    ),
    showProgress = FALSE
  )

  CROSS_LIST[[ii]] <-
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
    CROSS_LIST,
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

rm(
  CROSS_LIST
)


# ============================================================
# 6. Geometry
# ============================================================

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
  "distance"
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

GEO <- fread(
  GEO_FILE,
  select = need_geo,
  showProgress = TRUE
)

if (
  year_col !=
    "year"
) {
  setnames(
    GEO,
    year_col,
    "year"
  )
}

if (
  plant_col !=
    "plant_id"
) {
  setnames(
    GEO,
    plant_col,
    "plant_id"
  )
}

if (
  monitor_col !=
    "radnet_monitor_id"
) {
  setnames(
    GEO,
    monitor_col,
    "radnet_monitor_id"
  )
}

if (
  dist_col !=
    "distance_km"
) {
  setnames(
    GEO,
    dist_col,
    "distance_km"
  )
}

GEO[
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
  is.na(loc_col)
) {

  GEO <- merge(
    GEO,
    CROSS,
    by = "radnet_monitor_id",
    all.x = TRUE,
    sort = FALSE
  )

} else {

  if (
    loc_col !=
      "loc_num"
  ) {
    setnames(
      GEO,
      loc_col,
      "loc_num"
    )
  }

  GEO[
    ,
    loc_num :=
      clean_chr(
        loc_num
      )
  ]
}

GEO <- unique(
  GEO[
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
# 7. CAMPD eligible plant-month
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


# ============================================================
# 8. Exact fuel mapping used by H14
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


classify_fuel <- function(x) {

  if (
    is.na(x) ||
    !nzchar(
      trimws(
        as.character(x)
      )
    )
  ) {
    return(
      "unknown"
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

  groups <- unname(
    CAMPD_COMPONENT_MAP[
      toupper(parts)
    ]
  )

  if (
    any(
      is.na(groups)
    )
  ) {
    return(
      "unknown"
    )
  }

  ug <- unique(
    groups
  )

  if (
    length(ug) ==
      1L
  ) {
    return(
      ug[1]
    )
  }

  "mixed"
}


# ============================================================
# 9. Read actual activity for one year
# ============================================================

read_activity_year <- function(yy) {

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

  L <- lapply(
    files,
    function(f) {

      fread(
        f,
        select = c(
          "year",
          "month",
          "date",
          "plant_id",
          "primary_fuel_info",
          "daily_heat_input_mmbtu_observed"
        ),
        showProgress = FALSE
      )
    }
  )

  D <- rbindlist(
    L,
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

  D <- D[
    ELIGIBLE_PM[
      year ==
        yy
    ],
    on = .(
      plant_id,
      year,
      month
    ),
    nomatch = 0
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
    on = .(
      plant_id,
      year
    ),
    nomatch = 0
  ]

  labels <- unique(
    D$primary_fuel_info
  )

  MAP <- data.table(
    primary_fuel_info =
      labels,

    fuel_group =
      vapply(
        labels,
        classify_fuel,
        character(1)
      )
  )

  D <- merge(
    D,
    MAP,
    by =
      "primary_fuel_info",
    all.x =
      TRUE,
    sort =
      FALSE
  )

  D[
    fuel_group %in%
      EXTRA_FUELS,
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
}


# ============================================================
# 10. Build extra 360° columns
# ============================================================

AUDIT_LIST <- vector(
  "list",
  length(YEARS)
)

for (
  ii in seq_along(YEARS)
) {

  yy <- YEARS[ii]

  cat(
    "\n============================================================\n",
    "[BUILD] ",
    yy,
    "\n",
    "============================================================\n",
    sep = ""
  )

  old_file <- file.path(
    OLD_ND_DIR,
    paste0(
      "07H14_nondirectional_fuel_daily_exposure_",
      yy,
      ".csv.gz"
    )
  )

  out_file <- file.path(
    OUT,
    paste0(
      "07H15A_nondirectional_fuel_daily_exposure_",
      yy,
      ".csv.gz"
    )
  )

  BASE <- fread(
    old_file,
    showProgress = FALSE
  )

  BASE[
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
    BASE[
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
    ] >
      0L
  ) {
    stopf(
      yy,
      "年H14 base存在重复monitor-day。"
    )
  }

  DOMAIN <- BASE[
    ,
    .(
      loc_num,
      date
    )
  ]

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
      on = .(
        loc_num,
        date
      ),
      nomatch = 0
    ]

  } else {

    X <- data.table()
  }


  # ----------------------------------------------------------
  # initialize + fill
  # ----------------------------------------------------------

  for (
    fuel in EXTRA_FUELS
  ) {

    XF <- if (
      nrow(X) >
        0L
    ) {
      X[
        fuel_group ==
          fuel
      ]
    } else {
      data.table()
    }

    for (
      radius in RADII
    ) {

      count_name <- nd_col(
        fuel,
        "count",
        radius
      )

      heat_name <- nd_col(
        fuel,
        "heat",
        radius
      )

      BASE[
        ,
        (count_name) :=
          0
      ]

      BASE[
        ,
        (heat_name) :=
          0
      ]

      if (
        nrow(XF) >
          0L
      ) {

        XR <- XF[
          distance_km <=
            radius +
            1e-8
        ]

        if (
          nrow(XR) >
            0L
        ) {

          C <- XR[
            ,
            .(
              value =
                uniqueN(
                  plant_id
                )
            ),
            by = .(
              loc_num,
              date
            )
          ]

          H <- XR[
            ,
            .(
              value =
                sum_keep_all_na(
                  heat
                )
            ),
            by = .(
              loc_num,
              date
            )
          ]

          BASE[
            C,
            on = .(
              loc_num,
              date
            ),
            (count_name) :=
              i.value
          ]

          BASE[
            H,
            on = .(
              loc_num,
              date
            ),
            (heat_name) :=
              i.value
          ]
        }
      }
    }
  }


  # ----------------------------------------------------------
  # Audit
  # ----------------------------------------------------------

  extra_cols <- unlist(
    lapply(
      EXTRA_FUELS,
      function(fuel) {
        unlist(
          lapply(
            RADII,
            function(radius) {
              c(
                nd_col(
                  fuel,
                  "count",
                  radius
                ),
                nd_col(
                  fuel,
                  "heat",
                  radius
                )
              )
            }
          )
        )
      }
    )
  )

  miss <- setdiff(
    extra_cols,
    names(BASE)
  )

  if (
    length(miss) >
      0L
  ) {
    stopf(
      yy,
      "年补建后仍缺字段：",
      paste(
        miss,
        collapse = ", "
      )
    )
  }

  fwrite(
    BASE,
    out_file,
    compress =
      "gzip"
  )

  AUDIT_LIST[[ii]] <-
    data.table(
      year =
        yy,

      domain_rows =
        nrow(DOMAIN),

      geometry_rows =
        nrow(GY),

      extra_activity_rows =
        nrow(AY),

      linked_rows =
        nrow(X),

      other_fossil_nonzero_count_50km =
        sum(
          BASE[[nd_col("other_fossil", "count", 50)]] >
            0,
          na.rm = TRUE
        ),

      solid_waste_nonzero_count_50km =
        sum(
          BASE[[nd_col("solid_waste", "count", 50)]] >
            0,
          na.rm = TRUE
        ),

      status =
        "PASS"
    )

  cat(
    "saved: ",
    out_file,
    "\n",
    sep = ""
  )

  rm(
    BASE,
    DOMAIN,
    GY,
    AY,
    X
  )

  invisible(
    gc()
  )
}


# ============================================================
# 11. Final audit
# ============================================================

AUDIT <- rbindlist(
  AUDIT_LIST,
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  AUDIT,
  file.path(
    OUT,
    "07H15A_扩展燃料360逐日暴露审计.csv"
  )
)

OUT_FILES <- file.path(
  OUT,
  paste0(
    "07H15A_nondirectional_fuel_daily_exposure_",
    YEARS,
    ".csv.gz"
  )
)

N_FILES <- sum(
  file.exists(
    OUT_FILES
  )
)

cat(
  "\n============================================================\n",
  "07H15A0 FINAL AUDIT\n",
  "============================================================\n",
  sep = ""
)

print(
  AUDIT,
  nrows = Inf
)

cat(
  "\nfiles = ",
  N_FILES,
  "/24\n",
  sep = ""
)

if (
  N_FILES !=
    24L
) {
  stopf(
    "扩展360°文件不是24/24。"
  )
}

cat(
  "07H15A0 EXPANDED NON-DIRECTIONAL EXPOSURE: PASS\n"
)
