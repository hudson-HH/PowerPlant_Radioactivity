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
# 07H15C0
#
# 扩展原冻结07C fuel-specific exposure:
#
# energy:
#   wind / solar / hydro / nuclear / geothermal
#
# direction:
#   90° upwind = wind-from ±45°
#
# radii:
#   20 / 25 / 30 / 35 / 40 / 45 / 50 km
#
# metric:
#   plant count ONLY
#
# source:
#   frozen 07A fuel-group plant landscape
#   +
#   frozen formal 07C daily wind
#
# 重要：
#   这不是CAMPD actual-daily operation。
#   它代表每个monitor-day上风向扇区内的该能源类型电厂数量。
#
# Hard validation:
#   新建20 km必须逐日逐站完全复现冻结07C的20 km结果。
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

FORMAL_07A_ROOT <- file.path(
  OLD_FROZEN,
  "02_07A_RadNet发电站空间链接",
  "07A_radnet_powerplant_same_year_250km"
)

GEO_FUEL_FILE <- file.path(
  FORMAL_07A_ROOT,
  "04_monitor_year_fuel_group_link_250km_07A_2001_2024.csv"
)

WIND_DIR <- file.path(
  OLD_FROZEN,
  "07C_daily_upwind_powerplant_exposure",
  "yearly"
)

OLD_FUEL_DIR <- file.path(
  OLD_FROZEN,
  "07C_daily_upwind_powerplant_exposure",
  "daily_fuel_main90"
)

CROSSWALK_FILE <- file.path(
  ROOT,
  "radnet命名标准化",
  "radnet_181站_loc_num_radnet_monitor_id对照表.csv"
)

OUT <- file.path(
  PP,
  "06_07H_燃料特异CAMPD逐日暴露与模型",
  "07H15C0_非燃烧能源90度上风向Count_20_50km"
)

DAILY_DIR <- file.path(
  OUT,
  "daily_yearly"
)

dir.create(
  DAILY_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Parameters
# ============================================================

YEARS <- 2001:2024

FUELS <- c(
  "wind",
  "solar",
  "hydro",
  "nuclear",
  "geothermal"
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

HALF_ANGLE <- 45


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
    toupper(y) %in%
      c(
        "NA",
        "N/A",
        "NULL"
      )
  ] <- NA_character_

  y
}

norm_logical <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  z <- tolower(
    trimws(
      as.character(x)
    )
  )

  ans <- rep(
    NA,
    length(z)
  )

  ans[
    z %in%
      c(
        "true",
        "t",
        "1",
        "yes",
        "y"
      )
  ] <- TRUE

  ans[
    z %in%
      c(
        "false",
        "f",
        "0",
        "no",
        "n"
      )
  ] <- FALSE

  ans
}

pick_col <- function(
  nms,
  candidates,
  label
) {

  hit <- candidates[
    candidates %in%
      nms
  ]

  if (
    length(hit) ==
      0L
  ) {
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

  hit[1]
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

  lon1r <- lon1 *
    pi /
    180

  lat1r <- lat1 *
    pi /
    180

  lon2r <- lon2 *
    pi /
    180

  lat2r <- lat2 *
    pi /
    180

  dlon <- lon2r -
    lon1r

  yy <- sin(dlon) *
    cos(lat2r)

  xx <- cos(lat1r) *
    sin(lat2r) -
    sin(lat1r) *
    cos(lat2r) *
    cos(dlon)

  bearing <- atan2(
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

count_col <- function(
  fuel,
  radius
) {

  sprintf(
    "%s_upwind_n_plants_%dkm",
    fuel,
    radius
  )
}


# ============================================================
# 4. Input audit
# ============================================================

cat(
  "============================================================\n",
  "07H15C0 INPUT AUDIT\n",
  "============================================================\n",
  sep = ""
)

if (
  !file.exists(
    GEO_FUEL_FILE
  )
) {
  stopf(
    "缺少正式07A fuel geometry：\n",
    GEO_FUEL_FILE
  )
}

if (
  !file.exists(
    CROSSWALK_FILE
  )
) {
  stopf(
    "缺少crosswalk：\n",
    CROSSWALK_FILE
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

OLD_FUEL_FILES <- file.path(
  OLD_FUEL_DIR,
  paste0(
    "daily_fuel_90deg_20km_",
    YEARS,
    ".csv.gz"
  )
)

if (
  sum(
    file.exists(
      WIND_FILES
    )
  ) !=
    24L
) {
  stopf(
    "正式07C wind不是24/24。"
  )
}

if (
  sum(
    file.exists(
      OLD_FUEL_FILES
    )
  ) !=
    24L
) {
  stopf(
    "冻结07C fuel 20km不是24/24，无法做20km复现审计。"
  )
}

cat(
  "07A fuel geometry：PASS\n",
  "Formal 07C wind：24/24 PASS\n",
  "Frozen 07C fuel 20 km：24/24 PASS\n",
  sep = ""
)


# ============================================================
# 5. Crosswalk
# ============================================================

CW <- fread(
  CROSSWALK_FILE,
  showProgress = FALSE
)

if (
  !all(
    c(
      "loc_num",
      "radnet_monitor_id"
    ) %in%
      names(CW)
  )
) {
  stopf(
    "crosswalk缺loc_num/radnet_monitor_id。"
  )
}

CW <- unique(
  CW[
    ,
    .(
      loc_num =
        safe_int(
          loc_num
        ),

      radnet_monitor_id =
        clean_chr(
          radnet_monitor_id
        )
    )
  ]
)

if (
  nrow(CW) !=
    181L ||
  uniqueN(
    CW$loc_num
  ) !=
    181L ||
  uniqueN(
    CW$radnet_monitor_id
  ) !=
    181L
) {
  stopf(
    "crosswalk不是181站一对一。"
  )
}


# ============================================================
# 6. Read frozen 07A fuel geometry
# ============================================================

H <- names(
  fread(
    GEO_FUEL_FILE,
    nrows = 0,
    showProgress = FALSE
  )
)

id_col <- pick_col(
  H,
  c(
    "radnet_monitor_id"
  ),
  "radnet_monitor_id"
)

year_col <- pick_col(
  H,
  c(
    "year"
  ),
  "year"
)

plant_col <- pick_col(
  H,
  c(
    "plant_id"
  ),
  "plant_id"
)

fuel_col <- pick_col(
  H,
  c(
    "fuel_group"
  ),
  "fuel_group"
)

dist_col <- pick_col(
  H,
  c(
    "distance_km"
  ),
  "distance_km"
)

mlat_col <- pick_col(
  H,
  c(
    "monitor_latitude"
  ),
  "monitor_latitude"
)

mlon_col <- pick_col(
  H,
  c(
    "monitor_longitude"
  ),
  "monitor_longitude"
)

plat_col <- pick_col(
  H,
  c(
    "plant_latitude",
    "latitude_final",
    "analysis_latitude_final"
  ),
  "plant_latitude"
)

plon_col <- pick_col(
  H,
  c(
    "plant_longitude",
    "longitude_final",
    "analysis_longitude_final"
  ),
  "plant_longitude"
)

prox_col <- pick_col(
  H,
  c(
    "proximity_eligible"
  ),
  "proximity_eligible"
)

other_col <- pick_col(
  H,
  c(
    "supplemental_other_eligible"
  ),
  "supplemental_other_eligible"
)

G <- fread(
  GEO_FUEL_FILE,
  select = unique(
    c(
      id_col,
      year_col,
      plant_col,
      fuel_col,
      dist_col,
      mlat_col,
      mlon_col,
      plat_col,
      plon_col,
      prox_col,
      other_col
    )
  ),
  colClasses = setNames(
    "character",
    plant_col
  ),
  showProgress = TRUE
)

G <- data.table(

  radnet_monitor_id =
    clean_chr(
      G[[id_col]]
    ),

  year =
    safe_int(
      G[[year_col]]
    ),

  plant_id =
    clean_chr(
      G[[plant_col]]
    ),

  fuel_group =
    tolower(
      clean_chr(
        G[[fuel_col]]
      )
    ),

  distance_km =
    safe_num(
      G[[dist_col]]
    ),

  monitor_latitude =
    safe_num(
      G[[mlat_col]]
    ),

  monitor_longitude =
    safe_num(
      G[[mlon_col]]
    ),

  plant_latitude =
    safe_num(
      G[[plat_col]]
    ),

  plant_longitude =
    safe_num(
      G[[plon_col]]
    ),

  proximity_eligible =
    norm_logical(
      G[[prox_col]]
    ),

  supplemental_other_eligible =
    norm_logical(
      G[[other_col]]
    )
)

G <- G[
  year %in%
    YEARS &
  fuel_group %in%
    FUELS &
  is.finite(
    distance_km
  ) &
  distance_km >=
    0 &
  distance_km <=
    50 +
    1e-8 &
  (
    proximity_eligible %in%
      TRUE |
    supplemental_other_eligible %in%
      TRUE
  )
]

G <- merge(
  G,
  CW,
  by =
    "radnet_monitor_id",
  all.x =
    TRUE,
  sort =
    FALSE
)

if (
  anyNA(
    G$loc_num
  )
) {
  stopf(
    "07A fuel geometry存在无法映射的monitor。"
  )
}

if (
  G[
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
      1L,
    .N
  ] >
    0L
) {
  stopf(
    "07A fuel geometry存在重复key。"
  )
}

G[
  ,
  plant_bearing_deg :=
    initial_bearing_deg(
      monitor_longitude,
      monitor_latitude,
      plant_longitude,
      plant_latitude
    )
]

if (
  G[
    !is.finite(
      plant_bearing_deg
    ),
    .N
  ] >
    0L
) {
  stopf(
    "有效fuel geometry存在无法计算bearing的记录。"
  )
}

cat(
  "<=50 km eligible non-combustion links = ",
  format(
    nrow(G),
    big.mark = ","
  ),
  "\n",
  sep = ""
)


# ============================================================
# 7. Build yearly daily exposures
# ============================================================

AUDIT_LIST <- list()

for (
  yy in YEARS
) {

  cat(
    "\n============================================================\n",
    "[BUILD] ",
    yy,
    "\n",
    "============================================================\n",
    sep = ""
  )

  WF <- file.path(
    WIND_DIR,
    paste0(
      "monitor_day_wind_",
      yy,
      ".csv.gz"
    )
  )

  W <- fread(
    WF,
    showProgress = FALSE
  )

  REQ_W <- c(
    "loc_num",
    "date",
    "wind_from_deg",
    "wind_available",
    "monitor_year_07A_eligible",
    "plant_source_geography_covered"
  )

  MISS_W <- setdiff(
    REQ_W,
    names(W)
  )

  if (
    length(MISS_W) >
      0L
  ) {
    stopf(
      yy,
      "年正式wind缺字段：",
      paste(
        MISS_W,
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
        ),

      wind_available =
        norm_logical(
          wind_available
        ),

      monitor_year_07A_eligible =
        norm_logical(
          monitor_year_07A_eligible
        ),

      plant_source_geography_covered =
        norm_logical(
          plant_source_geography_covered
        )
    )
  ]

  OUTY <- unique(
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
  )

  if (
    OUTY[
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
      "年wind存在重复monitor-day。"
    )
  }

  valid_day <-
    OUTY$wind_available %in%
      TRUE &
    OUTY$monitor_year_07A_eligible %in%
      TRUE &
    OUTY$plant_source_geography_covered %in%
      TRUE

  for (
    fuel in FUELS
  ) {

    for (
      radius in RADII
    ) {

      cn <- count_col(
        fuel,
        radius
      )

      OUTY[
        ,
        (cn) :=
          fifelse(
            valid_day,
            0L,
            NA_integer_
          )
      ]
    }
  }


  # ----------------------------------------------------------
  # monitor-by-monitor to keep matrix small
  # ----------------------------------------------------------

  LOCS <- unique(
    OUTY$loc_num
  )

  for (
    loc in LOCS
  ) {

    idx <- which(
      OUTY$loc_num ==
        loc
    )

    if (
      length(idx) ==
        0L
    ) {
      next
    }

    wind_vec <- OUTY$wind_from_deg[
      idx
    ]

    GL <- G[
      year ==
        yy &
      loc_num ==
        loc
    ]

    if (
      nrow(GL) ==
        0L
    ) {
      next
    }

    for (
      fuel in FUELS
    ) {

      GF <- GL[
        fuel_group ==
          fuel
      ]

      if (
        nrow(GF) ==
          0L
      ) {
        next
      }

      for (
        radius in RADII
      ) {

        GR <- GF[
          distance_km <=
            radius +
            1e-8
        ]

        if (
          nrow(GR) ==
            0L
        ) {
          next
        }

        UP <- outer(
          wind_vec,
          GR$plant_bearing_deg,
          FUN = function(w, p) {
            angular_difference(
              p,
              w
            ) <=
              HALF_ANGLE
          }
        )

        counts <- rowSums(
          UP,
          na.rm = FALSE
        )

        cn <- count_col(
          fuel,
          radius
        )

        valid_loc <- valid_day[
          idx
        ]

        OUTY[
          idx[
            valid_loc
          ],
          (cn) :=
            as.integer(
              counts[
                valid_loc
              ]
            )
        ]
      }
    }
  }


  # ==========================================================
  # 20 km exact reproduction audit
  # ==========================================================

  OLD_FILE <- file.path(
    OLD_FUEL_DIR,
    paste0(
      "daily_fuel_90deg_20km_",
      yy,
      ".csv.gz"
    )
  )

  OLD_COLS <- c(
    "loc_num",
    "date",
    unname(
      vapply(
        FUELS,
        function(f) {
          count_col(
            f,
            20
          )
        },
        character(1)
      )
    )
  )

  OLD <- fread(
    OLD_FILE,
    select = OLD_COLS,
    showProgress = FALSE
  )

  OLD[
    ,
    `:=`(
      loc_num =
        safe_int(
          loc_num
        ),

      date =
        as.IDate(
          date
        )
    )
  ]

  NEW20 <- OUTY[
    ,
    ..OLD_COLS
  ]

  if (
    nrow(OLD) !=
      nrow(NEW20)
  ) {
    stopf(
      yy,
      "年20km复现审计行数不一致。"
    )
  }

  CMP <- merge(
    OLD,
    NEW20,
    by = c(
      "loc_num",
      "date"
    ),
    suffixes = c(
      "_old",
      "_new"
    ),
    all = TRUE,
    sort = FALSE
  )

  if (
    nrow(CMP) !=
      nrow(OLD)
  ) {
    stopf(
      yy,
      "年20km复现审计key不一致。"
    )
  }

  mismatch_by_fuel <- integer(
    length(
      FUELS
    )
  )

  names(
    mismatch_by_fuel
  ) <- FUELS

  for (
    fuel in FUELS
  ) {

    base <- count_col(
      fuel,
      20
    )

    old_name <- paste0(
      base,
      "_old"
    )

    new_name <- paste0(
      base,
      "_new"
    )

    a <- CMP[[old_name]]
    b <- CMP[[new_name]]

    equal_flag <-
      (
        is.na(a) &
        is.na(b)
      ) |
      (
        !is.na(a) &
        !is.na(b) &
        a ==
          b
      )

    mismatch_by_fuel[
      fuel
    ] <- sum(
      !equal_flag
    )
  }

  total_mismatch <- sum(
    mismatch_by_fuel
  )

  if (
    total_mismatch >
      0L
  ) {

    print(
      mismatch_by_fuel
    )

    stopf(
      yy,
      "年新建20km没有精确复现冻结07C，停止。"
    )
  }


  # ----------------------------------------------------------
  # Write only keys + new exposure columns
  # ----------------------------------------------------------

  EXP_COLS <- unlist(
    lapply(
      FUELS,
      function(f) {

        vapply(
          RADII,
          function(r) {
            count_col(
              f,
              r
            )
          },
          character(1)
        )
      }
    )
  )

  FINAL <- OUTY[
    ,
    c(
      "loc_num",
      "date",
      EXP_COLS
    ),
    with = FALSE
  ]

  OUT_FILE <- file.path(
    DAILY_DIR,
    paste0(
      "07H15C_noncombustion_upwind_count_20_50km_",
      yy,
      ".csv.gz"
    )
  )

  fwrite(
    FINAL,
    OUT_FILE,
    compress = "gzip"
  )

  AUDIT_LIST[[length(AUDIT_LIST) + 1L]] <- data.table(
    year =
      yy,

    n_monitor_days =
      nrow(FINAL),

    n_valid_days =
      sum(
        valid_day
      ),

    mismatch_20km =
      total_mismatch,

    wind_nonzero_50km =
      sum(
        FINAL[[
          count_col(
            "wind",
            50
          )
        ]] >
          0,
        na.rm = TRUE
      ),

    solar_nonzero_50km =
      sum(
        FINAL[[
          count_col(
            "solar",
            50
          )
        ]] >
          0,
        na.rm = TRUE
      ),

    hydro_nonzero_50km =
      sum(
        FINAL[[
          count_col(
            "hydro",
            50
          )
        ]] >
          0,
        na.rm = TRUE
      ),

    nuclear_nonzero_50km =
      sum(
        FINAL[[
          count_col(
            "nuclear",
            50
          )
        ]] >
          0,
        na.rm = TRUE
      ),

    geothermal_nonzero_50km =
      sum(
        FINAL[[
          count_col(
            "geothermal",
            50
          )
        ]] >
          0,
        na.rm = TRUE
      ),

    status =
      "PASS"
  )

  cat(
    "20 km exact reproduction：PASS\n",
    "saved: ",
    OUT_FILE,
    "\n",
    sep = ""
  )
}


# ============================================================
# 8. Final audit
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
    "01_07H15C0_build_audit.csv"
  )
)

FILES <- file.path(
  DAILY_DIR,
  paste0(
    "07H15C_noncombustion_upwind_count_20_50km_",
    YEARS,
    ".csv.gz"
  )
)

cat(
  "\n============================================================\n",
  "07H15C0 FINAL AUDIT\n",
  "============================================================\n",
  sep = ""
)

print(
  AUDIT,
  nrows = Inf
)

cat(
  "\nfiles = ",
  sum(
    file.exists(
      FILES
    )
  ),
  "/24\n",
  sep = ""
)

cat(
  "20 km total mismatches = ",
  sum(
    AUDIT$mismatch_20km
  ),
  "\n",
  sep = ""
)

if (
  sum(
    file.exists(
      FILES
    )
  ) !=
    24L
) {
  stopf(
    "输出文件不是24/24。"
  )
}

if (
  sum(
    AUDIT$mismatch_20km
  ) !=
    0L
) {
  stopf(
    "20km复现存在差异。"
  )
}

cat(
  "07H15C0 NON-COMBUSTION 20–50 km DAILY EXPOSURE: PASS\n"
)
