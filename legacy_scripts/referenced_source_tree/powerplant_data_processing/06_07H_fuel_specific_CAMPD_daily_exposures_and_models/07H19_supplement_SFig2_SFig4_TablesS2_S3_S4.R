#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999,
  width = 240
)

required_packages <- c(
  "data.table",
  "fixest",
  "ggplot2",
  "patchwork",
  "sf",
  "maps"
)

for (p in required_packages) {
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
  library(sf)
  library(maps)
})

N_CORES <- suppressWarnings(
  parallel::detectCores(
    logical = TRUE
  )
)

if (
  is.finite(N_CORES) &&
  N_CORES >= 2
) {
  fixest::setFixest_nthreads(
    max(
      1L,
      floor(N_CORES / 2)
    )
  )
}


# ============================================================
# 07H19
#
# Supplementary Figure S2
#   Three descriptive geographic panels with
#   20-km and 50-km monitor-centered buffers
#
# Supplementary Figure S4
#   Leave-one-monitor-out influence analysis
#   Core fuels × Count / Heat × 20 km
#
# Supplementary Table S2
#   Source data for current Figure 4
#
# Supplementary Table S3
#   Temporal stability:
#   2001–2008 / 2009–2016 / 2017–2024
#
# Supplementary Table S4
#   Previously generated directional-sector sensitivity
#   + downwind negative-control/direct contrast
#
# Main fuel-specific models:
#   filter-period
#   90° upwind
#   monitor × year FE
#   cluster by monitor
#   sample-median wind
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

FP_FILE <- file.path(
  F07,
  "03_07F_filter_period_analysis_base_2001_2024.csv.gz"
)

DAILY_DIR <- file.path(
  HROOT,
  "fuel_daily_yearly"
)

H12 <- file.path(
  HROOT,
  "07H12_论文正式图表总整理"
)

H16 <- file.path(
  HROOT,
  "07H16_Figure4_Figure5正文与MonitorFE补充图"
)

MAP_DIR <- file.path(
  H12,
  "05_绘图数据"
)

PLANT_MAP_FILE <- file.path(
  MAP_DIR,
  "Figure1_powerplant_locations_by_type.csv"
)

MONITOR_MAP_FILE <- file.path(
  MAP_DIR,
  "Figure1_radnet_monitor_locations.csv"
)

FIG4_SOURCE <- file.path(
  H16,
  "07_Figure4_source_monitorYearFE.csv"
)

OUT <- file.path(
  HROOT,
  "07H19_补充材料_SFig2_SFig4_TableS2_S3_S4"
)

FIG_DIR <- file.path(
  OUT,
  "figures"
)

TABLE_DIR <- file.path(
  OUT,
  "tables"
)

SOURCE_DIR <- file.path(
  OUT,
  "source_data"
)

AUDIT_DIR <- file.path(
  OUT,
  "audit"
)

for (
  d in c(
    OUT,
    FIG_DIR,
    TABLE_DIR,
    SOURCE_DIR,
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

CORE_FUELS <- c(
  "coal",
  "gas",
  "oil",
  "biomass"
)

METRICS <- c(
  "count",
  "heat"
)

FUEL_LABEL <- c(
  coal = "Coal",
  gas = "Natural gas",
  oil = "Oil",
  biomass = "Biomass"
)

METRIC_LABEL <- c(
  count = "Active plant count",
  heat = "Actual daily Heat Input"
)

EXPECTED_FP_ROWS <- 107407L
EXPECTED_MONITORS <- 133L

RADIUS_MAIN <- 20L

TOL <- 1e-7


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

  if (is.logical(x)) {

    ans <- x

  } else if (
    is.numeric(x) ||
    is.integer(x)
  ) {

    ans <- x != 0

  } else {

    z <- tolower(
      trimws(
        as.character(x)
      )
    )

    ans <- z %in% c(
      "true",
      "t",
      "1",
      "yes",
      "y"
    )
  }

  ans[
    is.na(ans)
  ] <- FALSE

  ans
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

raw_col <- function(
  fuel,
  direction,
  metric,
  radius
) {

  if (
    metric == "count"
  ) {

    return(
      sprintf(
        "campd_%s_%s_n_plants_%dkm_90deg",
        fuel,
        direction,
        radius
      )
    )
  }

  sprintf(
    "campd_%s_%s_heat_observed_%dkm_90deg",
    fuel,
    direction,
    radius
  )
}

metric_scale <- function(metric) {

  if (
    metric == "heat"
  ) {
    return(
      10000
    )
  }

  1
}

find_interaction <- function(
  coef_names,
  a,
  b
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
    length(z) != 1L
  ) {
    return(
      NA_character_
    )
  }

  z[1]
}

first_existing <- function(paths) {

  hit <- paths[
    file.exists(
      paths
    )
  ]

  if (
    length(hit) == 0L
  ) {
    return(
      NA_character_
    )
  }

  hit[1]
}


# ============================================================
# 4. Restore formal 133-monitor universe
# ============================================================

cat(
  "============================================================\n",
  "A. FORMAL MONITOR UNIVERSE\n",
  "============================================================\n",
  sep = ""
)

Q_MON <- rbindlist(
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
          "缺少quasi-daily文件：",
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

      z
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

ELIGIBLE_LOCS <- sort(
  unique(
    Q_MON[
      flag_07f_primary_quasi_sample %in%
        TRUE,
      loc_num
    ]
  )
)

rm(Q_MON)
invisible(gc())

if (
  length(
    ELIGIBLE_LOCS
  ) !=
    EXPECTED_MONITORS
) {

  stopf(
    "正式monitor universe不是133；实际=",
    length(
      ELIGIBLE_LOCS
    )
  )
}

cat(
  "Formal eligible monitors = 133：PASS\n"
)


# ============================================================
# 5. Supplementary Figure S2
#
# Reference-paper analogue:
# three descriptive geographic panels
# with 20-km and 50-km monitor-centered buffers.
#
# These West/Central/East panels are DESCRIPTIVE ONLY,
# not statistical analysis strata.
# ============================================================

cat(
  "\n============================================================\n",
  "B. SUPPLEMENTARY FIGURE S2\n",
  "Regional exposure geometry\n",
  "============================================================\n",
  sep = ""
)

for (
  f in c(
    PLANT_MAP_FILE,
    MONITOR_MAP_FILE
  )
) {

  if (
    !file.exists(f)
  ) {
    stopf(
      "缺少Figure1地图源数据：",
      f
    )
  }
}

P_MAP <- fread(
  PLANT_MAP_FILE,
  showProgress = FALSE
)

M_MAP <- fread(
  MONITOR_MAP_FILE,
  showProgress = FALSE
)

REQ_P <- c(
  "plant_longitude",
  "plant_latitude"
)

REQ_M <- c(
  "monitor_longitude",
  "monitor_latitude"
)

MISS_P <- setdiff(
  REQ_P,
  names(P_MAP)
)

MISS_M <- setdiff(
  REQ_M,
  names(M_MAP)
)

if (
  length(MISS_P) > 0L
) {
  stopf(
    "发电站地图文件缺字段：",
    paste(
      MISS_P,
      collapse = ", "
    )
  )
}

if (
  length(MISS_M) > 0L
) {
  stopf(
    "监测站地图文件缺字段：",
    paste(
      MISS_M,
      collapse = ", "
    )
  )
}

P_MAP[
  ,
  `:=`(
    plant_longitude =
      safe_num(
        plant_longitude
      ),

    plant_latitude =
      safe_num(
        plant_latitude
      )
  )
]

M_MAP[
  ,
  `:=`(
    monitor_longitude =
      safe_num(
        monitor_longitude
      ),

    monitor_latitude =
      safe_num(
        monitor_latitude
      )
  )
]

P_MAP <- P_MAP[
  is.finite(
    plant_longitude
  ) &
  is.finite(
    plant_latitude
  )
]

M_MAP <- M_MAP[
  is.finite(
    monitor_longitude
  ) &
  is.finite(
    monitor_latitude
  )
]

# ------------------------------------------------------------
# If loc_num is present, restrict to the formal 133 monitors.
# Otherwise retain the already curated Figure 1 monitor layer.
# ------------------------------------------------------------

MAP_MONITOR_SCOPE <- "Figure1 curated RadNet monitor layer"

if (
  "loc_num" %in%
    names(M_MAP)
) {

  M_MAP[
    ,
    loc_num :=
      clean_chr(
        loc_num
      )
  ]

  M_FORMAL <- M_MAP[
    loc_num %in%
      ELIGIBLE_LOCS
  ]

  if (
    uniqueN(
      M_FORMAL$loc_num
    ) ==
      EXPECTED_MONITORS
  ) {

    M_MAP <- M_FORMAL

    MAP_MONITOR_SCOPE <-
      "formal 133-monitor universe"
  }
}

# ------------------------------------------------------------
# Three descriptive CONUS panels.
# These are NOT model strata.
# ------------------------------------------------------------

assign_region <- function(lon) {

  fifelse(
    lon < -105,
    "A. Western CONUS",
    fifelse(
      lon < -90,
      "B. Central CONUS",
      "C. Eastern CONUS"
    )
  )
}

P_MAP[
  ,
  region :=
    assign_region(
      plant_longitude
    )
]

M_MAP[
  ,
  region :=
    assign_region(
      monitor_longitude
    )
]

# ------------------------------------------------------------
# Convert to sf and construct 20 / 50 km buffers
# ------------------------------------------------------------

P_SF <- st_as_sf(
  P_MAP,
  coords = c(
    "plant_longitude",
    "plant_latitude"
  ),
  crs = 4326,
  remove = FALSE
)

M_SF <- st_as_sf(
  M_MAP,
  coords = c(
    "monitor_longitude",
    "monitor_latitude"
  ),
  crs = 4326,
  remove = FALSE
)

P_5070 <- st_transform(
  P_SF,
  5070
)

M_5070 <- st_transform(
  M_SF,
  5070
)

# only show plants within 50 km of at least one displayed monitor
WITHIN50 <- lengths(
  st_is_within_distance(
    P_5070,
    M_5070,
    dist = 50000
  )
) > 0L

P_SF_NEAR <- P_SF[
  WITHIN50,
]

BUF20 <- st_transform(
  st_buffer(
    M_5070,
    dist = 20000
  ),
  4326
)

BUF50 <- st_transform(
  st_buffer(
    M_5070,
    dist = 50000
  ),
  4326
)

US_MAP_RAW <- maps::map(
  "state",
  plot = FALSE,
  fill = TRUE
)

US_SF <- st_as_sf(
  US_MAP_RAW
)

if (is.na(st_crs(US_SF))) {
  st_crs(US_SF) <- 4326
}

REGION_INFO <- data.table(
  region = c(
    "A. Western CONUS",
    "B. Central CONUS",
    "C. Eastern CONUS"
  ),

  xmin = c(
    -125,
    -105,
    -90
  ),

  xmax = c(
    -105,
    -90,
    -66
  ),

  ymin = c(
    24,
    24,
    24
  ),

  ymax = c(
    50,
    50,
    50
  )
)

make_region_map <- function(ii) {

  rr <- REGION_INFO[ii]

  reg <- rr$region

  PP <- P_SF_NEAR[
    P_SF_NEAR$region ==
      reg,
  ]

  MM <- M_SF[
    M_SF$region ==
      reg,
  ]

  B20 <- BUF20[
    BUF20$region ==
      reg,
  ]

  B50 <- BUF50[
    BUF50$region ==
      reg,
  ]

  ggplot() +
    geom_sf(
      data = US_SF,
      fill = "white",
      color = "grey80",
      linewidth = 0.25
    ) +
    geom_sf(
      data = B50,
      fill = NA,
      color = "grey55",
      linewidth = 0.35,
      linetype = "dashed"
    ) +
    geom_sf(
      data = B20,
      fill = NA,
      color = "grey20",
      linewidth = 0.45
    ) +
    geom_sf(
      data = PP,
      size = 0.55,
      alpha = 0.55,
      color = "grey45"
    ) +
    geom_sf(
      data = MM,
      shape = 17,
      size = 1.8,
      color = "black"
    ) +
    coord_sf(
      xlim = c(
        rr$xmin,
        rr$xmax
      ),
      ylim = c(
        rr$ymin,
        rr$ymax
      ),
      expand = FALSE
    ) +
    labs(
      title =
        reg,

      subtitle =
        paste0(
          nrow(MM),
          " RadNet monitors; ",
          nrow(PP),
          " power plants within 50 km"
        )
    ) +
    theme_void(
      base_size = 10
    ) +
    theme(
      plot.title =
        element_text(
          face = "bold"
        ),

      plot.subtitle =
        element_text(
          size = 8.5
        )
    )
}

S2_A <- make_region_map(1)
S2_B <- make_region_map(2)
S2_C <- make_region_map(3)

SFIG2 <- (
  S2_A |
  S2_B |
  S2_C
) +
  plot_annotation(
    title =
      "Supplementary Figure S2. Regional distribution of RadNet monitors and nearby power plants",

    subtitle =
      paste0(
        "Inner solid circles indicate 20-km buffers; outer dashed circles indicate 50-km buffers. ",
        "Panels are descriptive geographic divisions and are not statistical analysis strata."
      )
  )

SFIG2_PNG <- file.path(
  FIG_DIR,
  "Supplementary_Figure_S2_regional_exposure_geometry.png"
)

SFIG2_PDF <- file.path(
  FIG_DIR,
  "Supplementary_Figure_S2_regional_exposure_geometry.pdf"
)

ggsave(
  SFIG2_PNG,
  SFIG2,
  width = 15,
  height = 6.3,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  SFIG2_PDF,
  SFIG2,
  width = 15,
  height = 6.3,
  units = "in"
)

MAP_AUDIT <- data.table(
  item = c(
    "monitor_scope",
    "n_displayed_monitors",
    "n_nearby_power_plants",
    "n_west_monitors",
    "n_central_monitors",
    "n_east_monitors"
  ),

  value = c(
    MAP_MONITOR_SCOPE,
    nrow(M_MAP),
    nrow(P_SF_NEAR),
    M_MAP[
      region ==
        "A. Western CONUS",
      .N
    ],
    M_MAP[
      region ==
        "B. Central CONUS",
      .N
    ],
    M_MAP[
      region ==
        "C. Eastern CONUS",
      .N
    ]
  )
)

fwrite(
  MAP_AUDIT,
  file.path(
    AUDIT_DIR,
    "01_SFigureS2_map_audit.csv"
  )
)


# ============================================================
# 6. Formal filter-period sample for S4 and Table S3
# ============================================================

cat(
  "\n============================================================\n",
  "C. BUILD FORMAL FILTER-PERIOD SAMPLE\n",
  "============================================================\n",
  sep = ""
)

FP_COLS <- c(
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
  "temp",
  "flag_conus",
  "fpmean_campd_upwind_heat_observed_20km_90deg"
)

FP_HEADER <- names(
  fread(
    FP_FILE,
    nrows = 0,
    showProgress = FALSE
  )
)

FP_MISS <- setdiff(
  FP_COLS,
  FP_HEADER
)

if (
  length(FP_MISS) > 0L
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
  FP_FILE,
  select = FP_COLS,
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

NUM_COLS <- c(
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
  "temp",
  "fpmean_campd_upwind_heat_observed_20km_90deg"
)

for (
  cc in NUM_COLS
) {

  set(
    FP,
    j = cc,
    value =
      safe_num(
        FP[[cc]]
      )
  )
}

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
    EXPECTED_FP_ROWS
) {
  stopf(
    "formal filter-period不是107407；实际=",
    nrow(FP)
  )
}

if (
  uniqueN(
    FP$loc_num
  ) !=
    EXPECTED_MONITORS
) {
  stopf(
    "formal filter-period不是133 monitors。"
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
# 7. Aggregate 20-km fuel-specific daily exposures
# ============================================================

cat(
  "\n============================================================\n",
  "D. AGGREGATE 20-KM ACTUAL-OPERATION EXPOSURES\n",
  "============================================================\n",
  sep = ""
)

EXPOSURE_COLS <- unlist(
  lapply(
    CORE_FUELS,
    function(fuel) {

      unlist(
        lapply(
          METRICS,
          function(metric) {

            raw_col(
              fuel,
              "upwind",
              metric,
              RADIUS_MAIN
            )
          }
        ),
        use.names = FALSE
      )
    }
  ),
  use.names = FALSE
)

DAILY <- rbindlist(
  lapply(
    YEARS,
    function(yy) {

      f <- file.path(
        DAILY_DIR,
        paste0(
          "07H_fuel_daily_exposure_",
          yy,
          ".csv.gz"
        )
      )

      if (
        !file.exists(f)
      ) {
        stopf(
          "缺少07H01 daily exposure：",
          f
        )
      }

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
        EXPOSURE_COLS
      )

      miss <- setdiff(
        req,
        H
      )

      if (
        length(miss) > 0L
      ) {

        stopf(
          yy,
          "年daily exposure缺字段：",
          paste(
            miss,
            collapse = ", "
          )
        )
      }

      z <- fread(
        f,
        select = req,
        showProgress = FALSE
      )

      z[
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

          daily_present =
            1L
        )
      ]

      for (
        cc in EXPOSURE_COLS
      ) {

        set(
          z,
          j = cc,
          value =
            safe_num(
              z[[cc]]
            )
        )
      }

      z
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

if (
  DAILY[
    ,
    .N,
    by = .(
      loc_num,
      date
    )
  ][
    N > 1L,
    .N
  ] > 0L
) {

  stopf(
    "daily exposure存在重复monitor-date。"
  )
}

FP_DAYS <- FP[
  ,
  .(
    date =
      seq(
        collect_start,
        collect_end -
          1L,
        by = "day"
      )
  ),
  by = .(
    period_id,
    loc_num
  )
]

FP_DAYS[
  ,
  date :=
    as.IDate(
      date
    )
]

J <- merge(
  FP_DAYS,
  DAILY,
  by = c(
    "loc_num",
    "date"
  ),
  all.x = TRUE,
  sort = FALSE
)

UNMATCHED_DAILY <- J[
  is.na(
    daily_present
  ),
  .N
]

if (
  UNMATCHED_DAILY != 0L
) {

  stopf(
    "存在",
    UNMATCHED_DAILY,
    "个filter-period day无法连接daily exposure。"
  )
}

FP_EXP <- J[
  ,
  lapply(
    .SD,
    strict_mean
  ),
  by = period_id,
  .SDcols =
    EXPOSURE_COLS
]

D0 <- merge(
  FP,
  FP_EXP,
  by = "period_id",
  all.x = TRUE,
  sort = FALSE
)

if (
  nrow(D0) !=
    EXPECTED_FP_ROWS
) {
  stopf(
    "20-km exposure聚合后行数异常。"
  )
}

rm(
  DAILY,
  FP_DAYS,
  J,
  FP_EXP
)

invisible(gc())


# ============================================================
# 8. Main-model helpers
# ============================================================

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
  "temp_centered",
  "loc_num",
  "calendar_year",
  "monitor_year_fe"
)

FORM_MAIN <- as.formula(
  paste0(
    "gross_beta_mBq_m3 ~ ",
    "exposure + ",
    "exposure:vel + ",
    "c_prop + ",
    "sunspot_monthly + ",
    "pm_centered_10 + ",
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

make_model_data <- function(
  D,
  fuel,
  metric
) {

  exp_col <- raw_col(
    fuel,
    "upwind",
    metric,
    RADIUS_MAIN
  )

  scl <- metric_scale(
    metric
  )

  d <- D[
    ,
    c(
      MODEL_BASE_COLS,
      exp_col
    ),
    with = FALSE
  ]

  d[
    ,
    exposure :=
      safe_num(
        get(
          exp_col
        )
      ) /
      scl
  ]

  req <- c(
    MODEL_BASE_COLS,
    "exposure"
  )

  d <- d[
    complete.cases(
      d[
        ,
        ..req
      ]
    )
  ]

  d
}

fit_effect <- function(
  d,
  wind_ref
) {

  if (
    nrow(d) == 0L ||
    uniqueN(
      d$loc_num
    ) < 2L
  ) {

    return(
      data.table(
        status = "SKIP",
        estimate = NA_real_,
        se = NA_real_,
        ci95_low = NA_real_,
        ci95_high = NA_real_,
        p_wald = NA_real_,
        note = "Insufficient rows/monitors"
      )
    )
  }

  fit <- tryCatch(
    feols(
      FORM_MAIN,
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
      data.table(
        status = "FAIL",
        estimate = NA_real_,
        se = NA_real_,
        ci95_low = NA_real_,
        ci95_high = NA_real_,
        p_wald = NA_real_,
        note = conditionMessage(
          fit
        )
      )
    )
  }

  b <- coef(
    fit
  )

  int_name <- find_interaction(
    names(b),
    "exposure",
    "vel"
  )

  if (
    !"exposure" %in%
      names(b) ||
    is.na(
      int_name
    )
  ) {

    return(
      data.table(
        status = "SKIP",
        estimate = NA_real_,
        se = NA_real_,
        ci95_low = NA_real_,
        ci95_high = NA_real_,
        p_wald = NA_real_,
        note = "Exposure or exposure×wind term dropped"
      )
    )
  }

  V <- as.matrix(
    vcov(
      fit
    )
  )

  est <-
    unname(
      b["exposure"]
    ) +
    wind_ref *
    unname(
      b[int_name]
    )

  vv <-
    V[
      "exposure",
      "exposure"
    ] +
    wind_ref^2 *
    V[
      int_name,
      int_name
    ] +
    2 *
    wind_ref *
    V[
      "exposure",
      int_name
    ]

  se <- sqrt(
    max(
      as.numeric(vv),
      0
    )
  )

  low <- est -
    1.96 *
    se

  high <- est +
    1.96 *
    se

  p <- if (
    is.finite(se) &&
    se > 0
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
    status = "PASS",
    estimate = est,
    se = se,
    ci95_low = low,
    ci95_high = high,
    p_wald = p,
    note = ""
  )
}


# ============================================================
# 9. Full-sample 20-km models
#
# Used for:
# - leave-one-out reference
# - hard validation against current Figure 4
# - full-period row in Table S3
# ============================================================

cat(
  "\n============================================================\n",
  "E. FULL-SAMPLE 20-KM REFERENCE MODELS\n",
  "============================================================\n",
  sep = ""
)

FULL <- data.table()

for (
  fuel in CORE_FUELS
) {

  for (
    metric in METRICS
  ) {

    d <- make_model_data(
      D0,
      fuel,
      metric
    )

    ref_wind <- median(
      d$vel,
      na.rm = TRUE
    )

    ee <- fit_effect(
      d,
      ref_wind
    )

    FULL <- rbind(
      FULL,
      data.table(
        fuel = fuel,
        metric = metric,
        radius_km = RADIUS_MAIN,
        wind_reference = "sample_median",
        wind_ms = ref_wind,
        n_rows = nrow(d),
        n_monitors = uniqueN(
          d$loc_num
        ),
        status = ee$status,
        estimate = ee$estimate,
        se = ee$se,
        ci95_low = ee$ci95_low,
        ci95_high = ee$ci95_high,
        p_wald = ee$p_wald,
        note = ee$note
      ),
      fill = TRUE
    )
  }
}

if (
  FULL[
    status !=
      "PASS",
    .N
  ] > 0L
) {

  print(
    FULL[
      status !=
        "PASS"
    ]
  )

  stopf(
    "20-km完整样本模型存在非PASS。"
  )
}

fwrite(
  FULL,
  file.path(
    SOURCE_DIR,
    "SFigureS4_full_sample_reference_20km.csv"
  )
)


# ============================================================
# 10. Hard validation against current Figure 4
# ============================================================

if (
  !file.exists(
    FIG4_SOURCE
  )
) {
  stopf(
    "找不到当前Figure4 source：",
    FIG4_SOURCE
  )
}

F4 <- fread(
  FIG4_SOURCE,
  showProgress = FALSE
)

F4_20 <- F4[
  fe_structure ==
    "monitor_year_FE" &
  radius_km ==
    20 &
  fuel %in%
    CORE_FUELS &
  metric %in%
    METRICS &
  wind_reference ==
    "sample_median",
  .(
    fuel,
    metric,
    old_wind_ms =
      wind_ms,
    old_estimate =
      estimate,
    old_ci_low =
      ci95_low,
    old_ci_high =
      ci95_high
  )
]

CHECK_FULL <- merge(
  FULL[
    ,
    .(
      fuel,
      metric,
      new_wind_ms =
        wind_ms,
      new_estimate =
        estimate,
      new_ci_low =
        ci95_low,
      new_ci_high =
        ci95_high
    )
  ],
  F4_20,
  by = c(
    "fuel",
    "metric"
  ),
  all = TRUE,
  sort = FALSE
)

CHECK_FULL[
  ,
  max_abs_diff :=
    pmax(
      abs(
        new_estimate -
        old_estimate
      ),
      abs(
        new_ci_low -
        old_ci_low
      ),
      abs(
        new_ci_high -
        old_ci_high
      ),
      na.rm = TRUE
    )
]

FULL_MAX_DIFF <- max(
  CHECK_FULL$max_abs_diff,
  na.rm = TRUE
)

fwrite(
  CHECK_FULL,
  file.path(
    AUDIT_DIR,
    "02_SFigureS4_full_model_vs_Figure4_audit.csv"
  )
)

if (
  !is.finite(
    FULL_MAX_DIFF
  ) ||
  FULL_MAX_DIFF >
    TOL
) {

  stopf(
    "20-km完整模型没有复现当前Figure4；max diff=",
    FULL_MAX_DIFF
  )
}

cat(
  "20-km models reproduce current Figure 4：PASS\n",
  "max abs difference = ",
  format(
    FULL_MAX_DIFF,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)


# ============================================================
# 11. Supplementary Figure S4
# Leave-one-monitor-out
#
# Wind evaluation:
# for each fuel × metric, use the FULL-sample median wind
# for every leave-one-out refit.
# This prevents changes in reference wind from masquerading
# as monitor influence.
# ============================================================

cat(
  "\n============================================================\n",
  "F. SUPPLEMENTARY FIGURE S4\n",
  "Leave-one-monitor-out\n",
  "============================================================\n",
  sep = ""
)

LOO <- data.table()

for (
  fuel in CORE_FUELS
) {

  for (
    metric in METRICS
  ) {

    d <- make_model_data(
      D0,
      fuel,
      metric
    )

    fuel_i <- fuel
    metric_i <- metric

    ref_wind <- FULL[
      which(
        FULL$fuel == fuel_i &
        FULL$metric == metric_i
      ),
      wind_ms
    ][1]

    if (
      length(ref_wind) != 1L ||
      !is.finite(ref_wind)
    ) {
      stopf(
        "无法唯一取得LOO完整样本参考风速：",
        fuel_i,
        " / ",
        metric_i
      )
    }

    MONS <- sort(
      unique(
        d$loc_num
      )
    )

    cat(
      "\nLOO: ",
      fuel,
      " / ",
      metric,
      " | monitors = ",
      length(MONS),
      "\n",
      sep = ""
    )

    counter <- 0L

    for (
      omitted in MONS
    ) {

      counter <- counter +
        1L

      if (
        counter == 1L ||
        counter %% 20L == 0L ||
        counter == length(MONS)
      ) {

        cat(
          "  ",
          counter,
          "/",
          length(MONS),
          "\n",
          sep = ""
        )
      }

      dd <- d[
        loc_num !=
          omitted
      ]

      ee <- fit_effect(
        dd,
        ref_wind
      )

      LOO <- rbind(
        LOO,
        data.table(
          fuel = fuel,
          metric = metric,
          radius_km = 20L,
          omitted_monitor = omitted,
          wind_reference =
            "full_sample_median",
          wind_ms = ref_wind,
          n_rows = nrow(dd),
          n_monitors =
            uniqueN(
              dd$loc_num
            ),
          status = ee$status,
          estimate = ee$estimate,
          se = ee$se,
          ci95_low = ee$ci95_low,
          ci95_high = ee$ci95_high,
          p_wald = ee$p_wald,
          note = ee$note
        ),
        fill = TRUE
      )
    }
  }
}

fwrite(
  LOO,
  file.path(
    SOURCE_DIR,
    "SFigureS4_leave_one_monitor_out_results.csv"
  )
)

LOO[
  ,
  omitted_index :=
    frank(
      omitted_monitor,
      ties.method = "dense"
    ),
  by = .(
    fuel,
    metric
  )
]

LOO[
  ,
  fuel_label :=
    factor(
      FUEL_LABEL[
        fuel
      ],
      levels =
        unname(
          FUEL_LABEL[
            CORE_FUELS
          ]
        )
    )
]

LOO[
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

FULL_PLOT <- copy(
  FULL
)

FULL_PLOT[
  ,
  fuel_label :=
    factor(
      FUEL_LABEL[
        fuel
      ],
      levels =
        unname(
          FUEL_LABEL[
            CORE_FUELS
          ]
        )
    )
]

FULL_PLOT[
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

SFIG4 <- ggplot(
  LOO[
    status ==
      "PASS"
  ],
  aes(
    x = omitted_index,
    y = estimate
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dotted",
    linewidth = 0.35
  ) +
  geom_hline(
    data = FULL_PLOT,
    aes(
      yintercept =
        ci95_low
    ),
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey55"
  ) +
  geom_hline(
    data = FULL_PLOT,
    aes(
      yintercept =
        ci95_high
    ),
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey55"
  ) +
  geom_hline(
    data = FULL_PLOT,
    aes(
      yintercept =
        estimate
    ),
    linewidth = 0.65
  ) +
  geom_point(
    size = 1.25,
    alpha = 0.70
  ) +
  facet_grid(
    metric_label ~
      fuel_label,
    scales = "free_y"
  ) +
  labs(
    title =
      "Supplementary Figure S4. Leave-one-monitor-out influence analysis",

    subtitle =
      paste0(
        "20-km 90° upwind actual-operation models with monitor × year fixed effects. ",
        "Each point is the estimate after omitting one RadNet monitor; ",
        "solid lines are full-sample estimates and dashed lines are full-sample 95% confidence limits."
      ),

    x =
      "RadNet monitor omitted (ordered monitor index)",

    y =
      "Adjusted change in gross beta (mBq/m³)"
  ) +
  theme_bw(
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

    plot.title =
      element_text(
        face = "bold"
      )
  )

SFIG4_PNG <- file.path(
  FIG_DIR,
  "Supplementary_Figure_S4_leave_one_monitor_out.png"
)

SFIG4_PDF <- file.path(
  FIG_DIR,
  "Supplementary_Figure_S4_leave_one_monitor_out.pdf"
)

ggsave(
  SFIG4_PNG,
  SFIG4,
  width = 15,
  height = 8.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  SFIG4_PDF,
  SFIG4,
  width = 15,
  height = 8.5,
  units = "in"
)


# ============================================================
# 12. Supplementary Table S2
#
# Source data for CURRENT Figure 4.
# No model refitting.
# ============================================================

cat(
  "\n============================================================\n",
  "G. SUPPLEMENTARY TABLE S2\n",
  "Source data for Figure 4\n",
  "============================================================\n",
  sep = ""
)

S2 <- F4[
  fe_structure ==
    "monitor_year_FE" &
  fuel %in%
    CORE_FUELS &
  metric %in%
    METRICS &
  radius_km %in%
    c(
      20,
      25,
      30,
      35,
      40,
      45,
      50
    ) &
  wind_reference ==
    "sample_median"
]

if (
  nrow(S2) !=
    56L
) {
  stopf(
    "Table S2预期56行，实际=",
    nrow(S2)
  )
}

S2_OUT <- S2[
  ,
  .(
    Fuel =
      FUEL_LABEL[
        fuel
      ],

    Exposure_metric =
      METRIC_LABEL[
        metric
      ],

    Radius_km =
      radius_km,

    Wind_reference =
      "Sample median",

    Wind_speed_m_s =
      wind_ms,

    Estimate_mBq_m3 =
      estimate,

    CI95_lower =
      ci95_low,

    CI95_upper =
      ci95_high,

    N =
      n_rows
  )
]

setorder(
  S2_OUT,
  Fuel,
  Exposure_metric,
  Radius_km
)

TABLE_S2_FILE <- file.path(
  TABLE_DIR,
  "Supplementary_Table_S2_Figure4_source_data.csv"
)

fwrite(
  S2_OUT,
  TABLE_S2_FILE
)


# ============================================================
# 13. Supplementary Table S3
#
# Temporal stability:
# 2001–2008
# 2009–2016
# 2017–2024
#
# These are non-overlapping equal-duration periods.
# Radius = 20 km.
# Each period uses its own model-sample median wind.
# ============================================================

cat(
  "\n============================================================\n",
  "H. SUPPLEMENTARY TABLE S3\n",
  "Temporal stability\n",
  "============================================================\n",
  sep = ""
)

PERIODS <- data.table(
  period = c(
    "2001–2008",
    "2009–2016",
    "2017–2024"
  ),

  start_year = c(
    2001L,
    2009L,
    2017L
  ),

  end_year = c(
    2008L,
    2016L,
    2024L
  )
)

S3 <- data.table()

# ------------------------------------------------------------
# Full-period reference rows
# ------------------------------------------------------------

for (
  ii in seq_len(
    nrow(FULL)
  )
) {

  rr <- FULL[ii]

  S3 <- rbind(
    S3,
    data.table(
      period = "2001–2024",
      start_year = 2001L,
      end_year = 2024L,
      fuel = rr$fuel,
      metric = rr$metric,
      radius_km = 20L,
      wind_reference = "sample_median",
      wind_ms = rr$wind_ms,
      n_rows = rr$n_rows,
      n_monitors = rr$n_monitors,
      status = rr$status,
      estimate = rr$estimate,
      se = rr$se,
      ci95_low = rr$ci95_low,
      ci95_high = rr$ci95_high,
      p_wald = rr$p_wald,
      note = "Full-period reference"
    ),
    fill = TRUE
  )
}

# ------------------------------------------------------------
# Three non-overlapping periods
# ------------------------------------------------------------

for (
  pp in seq_len(
    nrow(PERIODS)
  )
) {

  period_name <- PERIODS$period[pp]
  y0 <- PERIODS$start_year[pp]
  y1 <- PERIODS$end_year[pp]

  cat(
    "\nPeriod ",
    period_name,
    "\n",
    sep = ""
  )

  for (
    fuel in CORE_FUELS
  ) {

    for (
      metric in METRICS
    ) {

      d <- make_model_data(
        D0,
        fuel,
        metric
      )

      d <- d[
        calendar_year >=
          y0 &
        calendar_year <=
          y1
      ]

      ref_wind <- if (
        nrow(d) > 0L
      ) {
        median(
          d$vel,
          na.rm = TRUE
        )
      } else {
        NA_real_
      }

      ee <- if (
        is.finite(
          ref_wind
        )
      ) {

        fit_effect(
          d,
          ref_wind
        )

      } else {

        data.table(
          status = "SKIP",
          estimate = NA_real_,
          se = NA_real_,
          ci95_low = NA_real_,
          ci95_high = NA_real_,
          p_wald = NA_real_,
          note = "No finite wind reference"
        )
      }

      S3 <- rbind(
        S3,
        data.table(
          period =
            period_name,
          start_year =
            y0,
          end_year =
            y1,
          fuel =
            fuel,
          metric =
            metric,
          radius_km =
            20L,
          wind_reference =
            "period_sample_median",
          wind_ms =
            ref_wind,
          n_rows =
            nrow(d),
          n_monitors =
            uniqueN(
              d$loc_num
            ),
          status =
            ee$status,
          estimate =
            ee$estimate,
          se =
            ee$se,
          ci95_low =
            ee$ci95_low,
          ci95_high =
            ee$ci95_high,
          p_wald =
            ee$p_wald,
          note =
            ee$note
        ),
        fill = TRUE
      )
    }
  }
}

S3[
  ,
  `:=`(
    Fuel =
      FUEL_LABEL[
        fuel
      ],

    Exposure_metric =
      METRIC_LABEL[
        metric
      ]
  )
]

S3_OUT <- S3[
  ,
  .(
    Period =
      period,

    Fuel,

    Exposure_metric,

    Radius_km =
      radius_km,

    Wind_reference =
      wind_reference,

    Wind_speed_m_s =
      wind_ms,

    N =
      n_rows,

    N_monitors =
      n_monitors,

    Estimate_mBq_m3 =
      estimate,

    CI95_lower =
      ci95_low,

    CI95_upper =
      ci95_high,

    Status =
      status
  )
]

PERIOD_ORDER <- c(
  "2001–2024",
  "2001–2008",
  "2009–2016",
  "2017–2024"
)

S3_OUT[
  ,
  Period :=
    factor(
      Period,
      levels =
        PERIOD_ORDER
    )
]

setorder(
  S3_OUT,
  Period,
  Fuel,
  Exposure_metric
)

S3_OUT[
  ,
  Period :=
    as.character(
      Period
    )
]

TABLE_S3_FILE <- file.path(
  TABLE_DIR,
  "Supplementary_Table_S3_temporal_stability_20km.csv"
)

fwrite(
  S3_OUT,
  TABLE_S3_FILE
)


# ============================================================
# 14. Supplementary Table S4
#
# Re-present already generated Figure 3 directional results.
#
# Latest preferred sources:
# - representative-wind joint up/down model
# - representative-wind 60/90/120 sector model
#
# This belongs to the EARLIER overall/broad-class analysis,
# not the 133-monitor fuel-specific Figure 4–6 sample.
# ============================================================

cat(
  "\n============================================================\n",
  "I. SUPPLEMENTARY TABLE S4\n",
  "Directional-sector sensitivity and downwind control\n",
  "============================================================\n",
  sep = ""
)

B_FILE <- first_existing(
  c(
    file.path(
      MAP_DIR,
      "Figure3B_joint_direction_at_representative_wind.csv"
    ),

    file.path(
      MAP_DIR,
      "Figure3B_20km_upwind_downwind_count.csv"
    )
  )
)

C_FILE <- first_existing(
  c(
    file.path(
      MAP_DIR,
      "Figure3C_sector_angle_at_representative_wind.csv"
    ),

    file.path(
      MAP_DIR,
      "Figure3C_sector_angle_60_90_120.csv"
    )
  )
)

if (
  is.na(B_FILE) ||
  is.na(C_FILE)
) {

  stopf(
    "无法定位当前Figure3 B/C source CSV。"
  )
}

cat(
  "Figure3B source = ",
  B_FILE,
  "\n",
  "Figure3C source = ",
  C_FILE,
  "\n",
  sep = ""
)

B3 <- fread(
  B_FILE,
  showProgress = FALSE
)

C3 <- fread(
  C_FILE,
  showProgress = FALSE
)

# ------------------------------------------------------------
# Normalize B
# ------------------------------------------------------------

if (
  "component" %in%
    names(B3)
) {

  B3[
    ,
    component_std :=
      tolower(
        trimws(
          as.character(
            component
          )
        )
      )
  ]

} else if (
  "direction" %in%
    names(B3)
) {

  B3[
    ,
    component_std :=
      tolower(
        trimws(
          as.character(
            direction
          )
        )
      )
  ]

} else {

  stopf(
    "Figure3B source缺少component/direction字段。"
  )
}

if (
  !"estimate_mBq_m3" %in%
    names(B3)
) {
  stopf(
    "Figure3B source缺少estimate_mBq_m3。"
  )
}

if (
  !"ci95_low_mBq_m3" %in%
    names(B3) ||
  !"ci95_high_mBq_m3" %in%
    names(B3)
) {
  stopf(
    "Figure3B source缺少95%CI字段。"
  )
}

if (
  !"wind_speed_ms" %in%
    names(B3)
) {

  B3[
    ,
    wind_speed_ms :=
      NA_real_
  ]
}

# normalize possible difference labels
B3[
  grepl(
    "difference|upwind.*downwind|上风向.*下风向",
    component_std
  ),
  component_std :=
    "difference"
]

B3[
  grepl(
    "^upwind$|上风向",
    component_std
  ) &
  component_std !=
    "difference",
  component_std :=
    "upwind"
]

B3[
  grepl(
    "^downwind$|下风向",
    component_std
  ) &
  component_std !=
    "difference",
  component_std :=
    "downwind"
]

if (
  B3[
    component_std ==
      "difference",
    .N
  ] !=
    1L
) {

  stopf(
    "Figure3B中无法唯一找到Upwind−Downwind直接差值。"
  )
}

if (
  B3[
    component_std ==
      "downwind",
    .N
  ] !=
    1L
) {

  stopf(
    "Figure3B中无法唯一找到downwind结果。"
  )
}

# ------------------------------------------------------------
# Normalize C
# ------------------------------------------------------------

REQ_C <- c(
  "sector_angle_deg",
  "estimate_mBq_m3",
  "ci95_low_mBq_m3",
  "ci95_high_mBq_m3"
)

MISS_C <- setdiff(
  REQ_C,
  names(C3)
)

if (
  length(MISS_C) > 0L
) {

  stopf(
    "Figure3C缺字段：",
    paste(
      MISS_C,
      collapse = ", "
    )
  )
}

if (
  !"wind_speed_ms" %in%
    names(C3)
) {

  C3[
    ,
    wind_speed_ms :=
      NA_real_
  ]
}

C3 <- C3[
  sector_angle_deg %in%
    c(
      60,
      90,
      120
    )
]

if (
  nrow(C3) !=
    3L
) {
  stopf(
    "Figure3C的60/90/120°结果不是3行。"
  )
}

# ------------------------------------------------------------
# Optional hard consistency check:
# B upwind 90° vs C 90°
# ------------------------------------------------------------

S4_ANGLE_DIFF <- NA_real_

if (
  B3[
    component_std ==
      "upwind",
    .N
  ] ==
    1L
) {

  S4_ANGLE_DIFF <- abs(
    B3[
      component_std ==
        "upwind",
      estimate_mBq_m3
    ] -
    C3[
      sector_angle_deg ==
        90,
      estimate_mBq_m3
    ]
  )
}

S4_SECTOR <- C3[
  ,
  .(
    Analysis =
      "Directional-sector sensitivity",

    Wind_direction =
      "Upwind",

    Sector_width_deg =
      sector_angle_deg,

    Radius_km =
      20,

    Wind_speed_m_s =
      wind_speed_ms,

    Estimate_mBq_m3 =
      estimate_mBq_m3,

    CI95_lower =
      ci95_low_mBq_m3,

    CI95_upper =
      ci95_high_mBq_m3
  )
]

S4_DOWN <- B3[
  component_std ==
    "downwind",
  .(
    Analysis =
      "Directional negative control",

    Wind_direction =
      "Downwind",

    Sector_width_deg =
      90,

    Radius_km =
      20,

    Wind_speed_m_s =
      wind_speed_ms,

    Estimate_mBq_m3 =
      estimate_mBq_m3,

    CI95_lower =
      ci95_low_mBq_m3,

    CI95_upper =
      ci95_high_mBq_m3
  )
]

S4_DIFF <- B3[
  component_std ==
    "difference",
  .(
    Analysis =
      "Direct directional contrast",

    Wind_direction =
      "Upwind − downwind",

    Sector_width_deg =
      90,

    Radius_km =
      20,

    Wind_speed_m_s =
      wind_speed_ms,

    Estimate_mBq_m3 =
      estimate_mBq_m3,

    CI95_lower =
      ci95_low_mBq_m3,

    CI95_upper =
      ci95_high_mBq_m3
  )
]

S4_OUT <- rbind(
  S4_SECTOR,
  S4_DOWN,
  S4_DIFF,
  fill = TRUE
)

setorder(
  S4_OUT,
  Analysis,
  Sector_width_deg
)

TABLE_S4_FILE <- file.path(
  TABLE_DIR,
  "Supplementary_Table_S4_direction_sector_and_downwind_control.csv"
)

fwrite(
  S4_OUT,
  TABLE_S4_FILE
)


# ============================================================
# 15. Save LOO summary
# ============================================================

# 先把完整样本estimate显式merge到LOO结果，
# 避免在data.table by=.()内部再嵌套查FULL。

FULL_KEY <- FULL[
  ,
  .(
    fuel,
    metric,
    full_estimate =
      estimate,
    full_ci95_low =
      ci95_low,
    full_ci95_high =
      ci95_high,
    full_wind_ms =
      wind_ms
  )
]

if (
  FULL_KEY[
    ,
    .N,
    by = .(
      fuel,
      metric
    )
  ][
    N != 1L,
    .N
  ] > 0L
) {
  stopf(
    "FULL_KEY存在非唯一fuel×metric组合。"
  )
}

LOO_SUM <- merge(
  LOO,
  FULL_KEY,
  by = c(
    "fuel",
    "metric"
  ),
  all.x = TRUE,
  sort = FALSE
)

if (
  any(
    !is.finite(
      LOO_SUM$full_estimate
    )
  )
) {
  stopf(
    "LOO结果无法完整匹配FULL reference estimate。"
  )
}

LOO_SUMMARY <- LOO_SUM[
  ,
  .(
    n_leave_one_out =
      .N,

    n_PASS =
      sum(
        status ==
          "PASS"
      ),

    n_nonPASS =
      sum(
        status !=
          "PASS"
      ),

    full_estimate =
      full_estimate[1],

    full_ci95_low =
      full_ci95_low[1],

    full_ci95_high =
      full_ci95_high[1],

    wind_ms =
      full_wind_ms[1],

    loo_min =
      if (
        any(
          status ==
            "PASS" &
          is.finite(
            estimate
          )
        )
      ) {
        min(
          estimate[
            status ==
              "PASS" &
            is.finite(
              estimate
            )
          ]
        )
      } else {
        NA_real_
      },

    loo_max =
      if (
        any(
          status ==
            "PASS" &
          is.finite(
            estimate
          )
        )
      ) {
        max(
          estimate[
            status ==
              "PASS" &
            is.finite(
              estimate
            )
          ]
        )
      } else {
        NA_real_
      },

    all_same_sign_as_full =
      if (
        all(
          status ==
            "PASS"
        )
      ) {
        all(
          sign(
            estimate
          ) ==
          sign(
            full_estimate
          )
        )
      } else {
        FALSE
      },

    max_abs_shift =
      if (
        any(
          status ==
            "PASS" &
          is.finite(
            estimate
          )
        )
      ) {
        max(
          abs(
            estimate[
              status ==
                "PASS" &
              is.finite(
                estimate
              )
            ] -
            full_estimate[1]
          )
        )
      } else {
        NA_real_
      }
  ),
  by = .(
    fuel,
    metric
  )
]

LOO_SUMMARY[
  ,
  fuel_order :=
    match(
      fuel,
      CORE_FUELS
    )
]

LOO_SUMMARY[
  ,
  metric_order :=
    match(
      metric,
      METRICS
    )
]

setorder(
  LOO_SUMMARY,
  fuel_order,
  metric_order
)

LOO_SUMMARY[
  ,
  c(
    "fuel_order",
    "metric_order"
  ) := NULL
]

fwrite(
  LOO_SUMMARY,
  file.path(
    TABLE_DIR,
    "Supplementary_Figure_S4_leave_one_out_summary.csv"
  )
)


# 16. Final hard checks
# ============================================================

EXPECTED_LOO_ROWS <- sum(
  FULL$n_monitors
)

N_LOO_NONPASS <- LOO[
  status !=
    "PASS",
  .N
]

N_S3_NONPASS <- S3[
  status !=
    "PASS",
  .N
]

OUTPUT_FILES <- c(
  SFIG2_PNG,
  SFIG2_PDF,
  SFIG4_PNG,
  SFIG4_PDF,
  TABLE_S2_FILE,
  TABLE_S3_FILE,
  TABLE_S4_FILE
)

HARD <- data.table(
  check = c(
    "formal_filter_rows_107407",
    "formal_monitors_133",
    "daily_join_unmatched_zero",
    "full_models_8_PASS",
    "full_models_match_current_Figure4",
    "leave_one_out_expected_rows",
    "leave_one_out_all_PASS",
    "TableS2_rows_56",
    "TableS3_rows_32",
    "TableS3_all_PASS",
    "TableS4_rows_5",
    "TableS4_joint_vs_sector_difference_informational",
    "seven_primary_outputs_exist"
  ),

  observed = c(
    as.character(
      nrow(FP)
    ),

    as.character(
      uniqueN(
        FP$loc_num
      )
    ),

    as.character(
      UNMATCHED_DAILY
    ),

    as.character(
      FULL[
        status ==
          "PASS",
        .N
      ]
    ),

    format(
      FULL_MAX_DIFF,
      scientific = TRUE
    ),

    as.character(
      nrow(LOO)
    ),

    as.character(
      N_LOO_NONPASS
    ),

    as.character(
      nrow(S2_OUT)
    ),

    as.character(
      nrow(S3_OUT)
    ),

    as.character(
      N_S3_NONPASS
    ),

    as.character(
      nrow(S4_OUT)
    ),

    ifelse(
      is.finite(
        S4_ANGLE_DIFF
      ),
      format(
        S4_ANGLE_DIFF,
        scientific = TRUE
      ),
      "not_available"
    ),

    as.character(
      sum(
        file.exists(
          OUTPUT_FILES
        )
      )
    )
  ),

  expected = c(
    "107407",
    "133",
    "0",
    "8",
    paste0(
      "<=",
      TOL
    ),
    as.character(
      EXPECTED_LOO_ROWS
    ),
    "0",
    "56",
    "32",
    "0",
    "5",
    "informational_only",
    "7"
  )
)

HARD[
  ,
  status := c(
    ifelse(
      nrow(FP) ==
        EXPECTED_FP_ROWS,
      "PASS",
      "FAIL"
    ),

    ifelse(
      uniqueN(
        FP$loc_num
      ) ==
        EXPECTED_MONITORS,
      "PASS",
      "FAIL"
    ),

    ifelse(
      UNMATCHED_DAILY ==
        0L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      FULL[
        status ==
          "PASS",
        .N
      ] ==
        8L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      is.finite(
        FULL_MAX_DIFF
      ) &&
      FULL_MAX_DIFF <=
        TOL,
      "PASS",
      "FAIL"
    ),

    ifelse(
      nrow(LOO) ==
        EXPECTED_LOO_ROWS,
      "PASS",
      "FAIL"
    ),

    ifelse(
      N_LOO_NONPASS ==
        0L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      nrow(S2_OUT) ==
        56L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      nrow(S3_OUT) ==
        32L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      N_S3_NONPASS ==
        0L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      nrow(S4_OUT) ==
        5L,
      "PASS",
      "FAIL"
    ),

    "PASS",

    ifelse(
      sum(
        file.exists(
          OUTPUT_FILES
        )
      ) ==
        7L,
      "PASS",
      "FAIL"
    )
  )
]

fwrite(
  HARD,
  file.path(
    AUDIT_DIR,
    "03_07H19_hard_checks.csv"
  )
)


# ============================================================
# 17. Documentation
# ============================================================

NOTE <- c(
  "# 07H19 Supplementary figures and tables",
  "",
  "## Supplementary Figure S2",
  "",
  "- Descriptive regional exposure-geometry figure.",
  "- Three geographic panels: Western, Central, and Eastern CONUS.",
  "- Solid circles represent 20-km monitor-centered buffers.",
  "- Dashed circles represent 50-km monitor-centered buffers.",
  "- These geographic panels are descriptive only and are not statistical analysis strata.",
  "",
  "## Supplementary Figure S4",
  "",
  "- Leave-one-monitor-out influence analysis.",
  "- Core fuels: coal, natural gas, oil, biomass.",
  "- Metrics: active plant Count and actual daily Heat Input.",
  "- Radius: 20 km.",
  "- 90-degree upwind exposure.",
  "- Monitor-by-year fixed effects.",
  "- Standard errors clustered by monitor.",
  "- For each fuel/metric, all leave-one-out models are evaluated at the full-sample median wind speed.",
  "",
  "## Supplementary Table S2",
  "",
  "- Source data for the current Figure 4.",
  "- 4 fuels × 2 metrics × 7 radii = 56 rows.",
  "",
  "## Supplementary Table S3",
  "",
  "- Temporal stability analysis.",
  "- Full period: 2001–2024.",
  "- Non-overlapping 8-year periods: 2001–2008, 2009–2016, 2017–2024.",
  "- 20-km 90-degree upwind actual-operation models.",
  "- Each sub-period effect is evaluated at its own complete-case sample median wind speed.",
  "",
  "## Supplementary Table S4",
  "",
  "- Re-presents the previously generated overall Figure 3 directional analyses.",
  "- Includes 60/90/120-degree upwind sector sensitivity at 20 km.",
  "- Includes 90-degree downwind estimate.",
  "- Includes the direct Upwind-minus-Downwind contrast.",
  "- This table belongs to the earlier overall/broad-class analysis sample and should not be described as the 133-monitor fuel-specific sample."
)

writeLines(
  NOTE,
  file.path(
    OUT,
    "00_07H19说明.md"
  ),
  useBytes = TRUE
)


# ============================================================
# 18. Terminal summary
# ============================================================

cat(
  "\n============================================================\n",
  "SUPPLEMENTARY FIGURE S4 LOO SUMMARY\n",
  "============================================================\n",
  sep = ""
)

print(
  LOO_SUMMARY[
    order(
      fuel,
      metric
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "SUPPLEMENTARY TABLE S3 TEMPORAL STABILITY\n",
  "============================================================\n",
  sep = ""
)

print(
  S3_OUT,
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "SUPPLEMENTARY TABLE S4\n",
  "============================================================\n",
  sep = ""
)

print(
  S4_OUT,
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

N_FAIL <- HARD[
  status ==
    "FAIL",
  .N
]

cat(
  "\n============================================================\n",
  ifelse(
    N_FAIL ==
      0L,
    "07H19 SUPPLEMENTARY BUILD: PASS\n",
    "07H19 SUPPLEMENTARY BUILD: FAIL\n"
  ),
  "Hard FAIL = ",
  N_FAIL,
  "\n",
  "Output = ",
  OUT,
  "\n",
  "============================================================\n",
  sep = ""
)

if (
  N_FAIL >
    0L
) {
  stopf(
    "07H19存在Hard FAIL，请检查audit/03_07H19_hard_checks.csv。"
  )
}

# >>> 07H19 FINAL REVISED S2/S4 REDRAW START >>>

cat(
  "\n============================================================\n",
  "J. FINAL REVISED SUPPLEMENTARY FIGURES S2 / S4\n",
  "============================================================\n",
  sep = ""
)

required_final_packages <- c(
  "data.table",
  "ggplot2",
  "patchwork",
  "sf",
  "maps",
  "lme4"
)

for (p in required_final_packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(sf)
  library(maps)
  library(lme4)
})


# ============================================================
# J1. Paths
# ============================================================

ROOT_FINAL <- "/Users/hudson/Desktop/USA 辐射表征研究"

HROOT_FINAL <- file.path(
  ROOT_FINAL,
  "发电站数据处理",
  "06_07H_燃料特异CAMPD逐日暴露与模型"
)

H12_FINAL <- file.path(
  HROOT_FINAL,
  "07H12_论文正式图表总整理"
)

OUT_FINAL <- file.path(
  HROOT_FINAL,
  "07H19_补充材料_SFig2_SFig4_TableS2_S3_S4"
)

FIG_FINAL <- file.path(
  OUT_FINAL,
  "figures"
)

SOURCE_FINAL <- file.path(
  OUT_FINAL,
  "source_data"
)

TABLE_FINAL <- file.path(
  OUT_FINAL,
  "tables"
)

AUDIT_FINAL <- file.path(
  OUT_FINAL,
  "audit"
)

dir.create(
  FIG_FINAL,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  SOURCE_FINAL,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TABLE_FINAL,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  AUDIT_FINAL,
  recursive = TRUE,
  showWarnings = FALSE
)


safe_num_final <- function(x) {
  suppressWarnings(
    as.numeric(x)
  )
}

clean_chr_final <- function(x) {

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

as_bool_final <- function(x) {

  if (is.logical(x)) {

    ans <- x

  } else if (
    is.numeric(x) ||
    is.integer(x)
  ) {

    ans <- x != 0

  } else {

    z <- tolower(
      trimws(
        as.character(x)
      )
    )

    ans <- z %in% c(
      "true",
      "t",
      "1",
      "yes",
      "y"
    )
  }

  ans[
    is.na(ans)
  ] <- FALSE

  ans
}


# ============================================================
# J2. Rebuild formal 133-monitor universe for S2
# ============================================================

cat(
  "\n------------------------------------------------------------\n",
  "J2. Rebuilding Supplementary Figure S2\n",
  "------------------------------------------------------------\n",
  sep = ""
)

QDIR_FINAL <- file.path(
  ROOT_FINAL,
  "发电站数据处理",
  "04_07F_正式统计分析数据构建",
  "quasi_daily_yearly"
)

QMON_FINAL <- rbindlist(
  lapply(
    2001:2024,
    function(yy) {

      f <- file.path(
        QDIR_FINAL,
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
          "flag_07f_primary_quasi_sample"
        ),
        showProgress = FALSE
      )

      z[
        ,
        loc_num :=
          clean_chr_final(
            loc_num
          )
      ]

      z[
        ,
        flag_07f_primary_quasi_sample :=
          as_bool_final(
            flag_07f_primary_quasi_sample
          )
      ]

      z
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

ELIGIBLE_FINAL <- sort(
  unique(
    QMON_FINAL[
      flag_07f_primary_quasi_sample %in%
        TRUE,
      loc_num
    ]
  )
)

rm(QMON_FINAL)

if (
  length(
    ELIGIBLE_FINAL
  ) != 133L
) {

  stop(
    paste0(
      "S2正式monitor universe不是133；实际=",
      length(
        ELIGIBLE_FINAL
      )
    ),
    call. = FALSE
  )
}


# ============================================================
# J3. Revised S2
# ============================================================

MAP_DIR_FINAL <- file.path(
  H12_FINAL,
  "05_绘图数据"
)

PLANT_FILE_FINAL <- file.path(
  MAP_DIR_FINAL,
  "Figure1_powerplant_locations_by_type.csv"
)

MONITOR_FILE_FINAL <- file.path(
  MAP_DIR_FINAL,
  "Figure1_radnet_monitor_locations.csv"
)

P_FINAL <- fread(
  PLANT_FILE_FINAL,
  showProgress = FALSE
)

M_FINAL <- fread(
  MONITOR_FILE_FINAL,
  showProgress = FALSE
)

P_FINAL[
  ,
  `:=`(
    plant_longitude =
      safe_num_final(
        plant_longitude
      ),

    plant_latitude =
      safe_num_final(
        plant_latitude
      )
  )
]

M_FINAL[
  ,
  `:=`(
    monitor_longitude =
      safe_num_final(
        monitor_longitude
      ),

    monitor_latitude =
      safe_num_final(
        monitor_latitude
      )
  )
]

P_FINAL <- P_FINAL[
  is.finite(
    plant_longitude
  ) &
  is.finite(
    plant_latitude
  )
]

M_FINAL <- M_FINAL[
  is.finite(
    monitor_longitude
  ) &
  is.finite(
    monitor_latitude
  )
]

if (
  "loc_num" %in%
    names(M_FINAL)
) {

  M_FINAL[
    ,
    loc_num :=
      clean_chr_final(
        loc_num
      )
  ]

  M_FINAL <- M_FINAL[
    loc_num %in%
      ELIGIBLE_FINAL
  ]
}

if (
  nrow(
    M_FINAL
  ) !=
    133L
) {

  stop(
    paste0(
      "S2地图正式monitor不是133；实际=",
      nrow(
        M_FINAL
      )
    ),
    call. = FALSE
  )
}


assign_region_final <- function(lon) {

  fifelse(
    lon < -105,
    "A. Western CONUS",
    fifelse(
      lon < -90,
      "B. Central CONUS",
      "C. Eastern CONUS"
    )
  )
}

P_FINAL[
  ,
  region :=
    assign_region_final(
      plant_longitude
    )
]

M_FINAL[
  ,
  region :=
    assign_region_final(
      monitor_longitude
    )
]

P_SF_FINAL <- st_as_sf(
  P_FINAL,
  coords = c(
    "plant_longitude",
    "plant_latitude"
  ),
  crs = 4326,
  remove = FALSE
)

M_SF_FINAL <- st_as_sf(
  M_FINAL,
  coords = c(
    "monitor_longitude",
    "monitor_latitude"
  ),
  crs = 4326,
  remove = FALSE
)

P_5070_FINAL <- st_transform(
  P_SF_FINAL,
  5070
)

M_5070_FINAL <- st_transform(
  M_SF_FINAL,
  5070
)

WITHIN50_FINAL <- lengths(
  st_is_within_distance(
    P_5070_FINAL,
    M_5070_FINAL,
    dist = 50000
  )
) > 0L

P_NEAR_FINAL <- P_SF_FINAL[
  WITHIN50_FINAL,
]

B20_FINAL <- st_transform(
  st_buffer(
    M_5070_FINAL,
    dist = 20000
  ),
  4326
)

B50_FINAL <- st_transform(
  st_buffer(
    M_5070_FINAL,
    dist = 50000
  ),
  4326
)

US_RAW_FINAL <- maps::map(
  "state",
  plot = FALSE,
  fill = TRUE
)

US_FINAL <- st_as_sf(
  US_RAW_FINAL
)

if (
  is.na(
    st_crs(
      US_FINAL
    )
  )
) {

  US_FINAL <- st_set_crs(
    US_FINAL,
    4326
  )
}


REGION_FINAL <- data.table(
  region = c(
    "A. Western CONUS",
    "B. Central CONUS",
    "C. Eastern CONUS"
  ),

  xmin = c(
    -125,
    -105,
    -90
  ),

  xmax = c(
    -105,
    -90,
    -66
  ),

  ymin = c(
    24,
    24,
    24
  ),

  ymax = c(
    50,
    50,
    50
  )
)


REGION_COUNT_FINAL <- rbindlist(
  lapply(
    REGION_FINAL$region,
    function(reg) {

      data.table(
        region =
          reg,

        n_monitors =
          M_FINAL[
            region ==
              reg,
            .N
          ],

        n_plants =
          nrow(
            P_NEAR_FINAL[
              P_NEAR_FINAL$region ==
                reg,
            ]
          )
      )
    }
  )
)


make_region_final <- function(ii) {

  rr <- REGION_FINAL[ii]

  reg <- rr$region

  MM <- M_SF_FINAL[
    M_SF_FINAL$region ==
      reg,
  ]

  PP <- P_NEAR_FINAL[
    P_NEAR_FINAL$region ==
      reg,
  ]

  BB20 <- B20_FINAL[
    B20_FINAL$region ==
      reg,
  ]

  BB50 <- B50_FINAL[
    B50_FINAL$region ==
      reg,
  ]

  CC <- REGION_COUNT_FINAL[
    region ==
      reg
  ]

  ggplot() +

    geom_sf(
      data = US_FINAL,
      fill = "white",
      color = "grey82",
      linewidth = 0.25
    ) +

    geom_sf(
      data = BB50,
      fill = NA,
      color = "grey55",
      linewidth = 0.42,
      linetype = "dashed"
    ) +

    geom_sf(
      data = BB20,
      fill = NA,
      color = "grey20",
      linewidth = 0.50
    ) +

    geom_sf(
      data = PP,
      size = 0.52,
      alpha = 0.35,
      color = "grey45"
    ) +

    geom_sf(
      data = MM,
      shape = 17,
      size = 1.8,
      color = "black"
    ) +

    coord_sf(
      xlim = c(
        rr$xmin,
        rr$xmax
      ),
      ylim = c(
        rr$ymin,
        rr$ymax
      ),
      expand = FALSE
    ) +

    labs(
      title =
        reg,

      subtitle =
        paste0(
          CC$n_monitors,
          " RadNet monitors\n",
          format(
            CC$n_plants,
            big.mark = ","
          ),
          " power plants within 50 km"
        )
    ) +

    theme_void(
      base_size = 10
    ) +

    theme(
      plot.title =
        element_text(
          face = "bold",
          size = 11,
          hjust = 0
        ),

      plot.subtitle =
        element_text(
          size = 8.5,
          lineheight = 1.05
        ),

      plot.margin =
        margin(
          8,
          12,
          8,
          12
        ),

      panel.border =
        element_rect(
          color = "grey75",
          fill = NA,
          linewidth = 0.45
        )
    )
}


S2A_FINAL <- make_region_final(1)
S2B_FINAL <- make_region_final(2)
S2C_FINAL <- make_region_final(3)

S2_TOP_FINAL <- (
  S2A_FINAL |
  S2B_FINAL
)

S2_FINAL <- (
  S2_TOP_FINAL /
  S2C_FINAL
) +

  plot_layout(
    heights = c(
      1,
      1
    )
  ) +

  plot_annotation(
    title =
      "Supplementary Figure S2. Regional distribution of RadNet monitors and nearby power plants",

    subtitle =
      paste0(
        "Solid circles indicate 20-km monitor-centered buffers; ",
        "dashed circles indicate 50-km buffers. ",
        "Black triangles represent RadNet monitors and grey points represent nearby power plants. ",
        "Panels are descriptive geographic divisions and are not statistical analysis strata."
      ),

    theme =
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 14
          ),

        plot.subtitle =
          element_text(
            size = 9.5
          )
      )
  )

S2_PNG_FINAL <- file.path(
  FIG_FINAL,
  "Supplementary_Figure_S2_regional_exposure_geometry.png"
)

S2_PDF_FINAL <- file.path(
  FIG_FINAL,
  "Supplementary_Figure_S2_regional_exposure_geometry.pdf"
)

ggsave(
  S2_PNG_FINAL,
  S2_FINAL,
  width = 14,
  height = 10.5,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  S2_PDF_FINAL,
  S2_FINAL,
  width = 14,
  height = 10.5,
  units = "in"
)

fwrite(
  REGION_COUNT_FINAL,
  file.path(
    SOURCE_FINAL,
    "SFigureS2_regional_counts.csv"
  )
)

cat(
  "Revised S2 saved：PASS\n"
)


# ============================================================
# J4. Existing fuel-specific LOO
# ============================================================

cat(
  "\n------------------------------------------------------------\n",
  "J4. Reading existing fuel-specific LOO results\n",
  "------------------------------------------------------------\n",
  sep = ""
)

FUEL_LOO_FILE_FINAL <- file.path(
  SOURCE_FINAL,
  "SFigureS4_leave_one_monitor_out_results.csv"
)

FUEL_FULL_FILE_FINAL <- file.path(
  SOURCE_FINAL,
  "SFigureS4_full_sample_reference_20km.csv"
)

if (
  !file.exists(
    FUEL_LOO_FILE_FINAL
  ) ||
  !file.exists(
    FUEL_FULL_FILE_FINAL
  )
) {

  stop(
    "缺少现有fuel-specific LOO结果。",
    call. = FALSE
  )
}

FUEL_LOO_FINAL <- fread(
  FUEL_LOO_FILE_FINAL,
  showProgress = FALSE
)

FUEL_FULL_FINAL <- fread(
  FUEL_FULL_FILE_FINAL,
  showProgress = FALSE
)

if (
  nrow(
    FUEL_LOO_FINAL
  ) !=
    1064L ||
  FUEL_LOO_FINAL[
    status !=
      "PASS",
    .N
  ] !=
    0L
) {

  stop(
    "现有fuel-specific LOO不是完整1064个PASS结果。",
    call. = FALSE
  )
}

cat(
  "Fuel-specific LOO = 1064 PASS；不重新拟合。\n"
)


# ============================================================
# J5. Overall Figure 3 data
# ============================================================

cat(
  "\n------------------------------------------------------------\n",
  "J5. Overall plant-count leave-one-monitor-out\n",
  "------------------------------------------------------------\n",
  sep = ""
)

OVERALL_INPUT_FINAL <- file.path(
  ROOT_FINAL,
  "model_build",
  "08_最终模型分析数据集构建",
  "01_final_quasi_daily_model_dataset_2001_2024.csv.gz"
)

FIG3A_FILE_FINAL <- file.path(
  H12_FINAL,
  "05_绘图数据",
  "Figure3A_20_150km_wind1_5_and_representative.csv"
)

if (
  !file.exists(
    OVERALL_INPUT_FINAL
  )
) {

  stop(
    paste0(
      "缺少overall frozen dataset：",
      OVERALL_INPUT_FINAL
    ),
    call. = FALSE
  )
}

if (
  !file.exists(
    FIG3A_FILE_FINAL
  )
) {

  stop(
    paste0(
      "缺少Figure3A source：",
      FIG3A_FILE_FINAL
    ),
    call. = FALSE
  )
}


OVERALL_FINAL <- fread(
  OVERALL_INPUT_FINAL,
  showProgress = TRUE
)

REQ_OVERALL_FINAL <- c(
  "loc_num",
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
  "all_upwind_n_plants_20km"
)

MISS_OVERALL_FINAL <- setdiff(
  REQ_OVERALL_FINAL,
  names(
    OVERALL_FINAL
  )
)

if (
  length(
    MISS_OVERALL_FINAL
  ) >
    0L
) {

  stop(
    paste0(
      "overall dataset缺字段：",
      paste(
        MISS_OVERALL_FINAL,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

OVERALL_FINAL <- OVERALL_FINAL[
  flag_reference_main_sample %in%
    TRUE
]

OVERALL_FINAL[
  ,
  loc_num :=
    clean_chr_final(
      loc_num
    )
]

OVERALL_FINAL[
  ,
  `:=`(
    loc_num_factor =
      factor(
        loc_num
      ),

    temp_centered =
      safe_num_final(
        temp
      ) -
      273.15
  )
]

NUM_OVERALL_FINAL <- c(
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
  "temp_centered",
  "all_upwind_n_plants_20km"
)

for (
  cc in NUM_OVERALL_FINAL
) {

  set(
    OVERALL_FINAL,
    j = cc,
    value =
      safe_num_final(
        OVERALL_FINAL[[cc]]
      )
  )
}

if (
  nrow(
    OVERALL_FINAL
  ) !=
    558915L
) {

  stop(
    paste0(
      "overall样本不是558,915；实际=",
      nrow(
        OVERALL_FINAL
      )
    ),
    call. = FALSE
  )
}

if (
  uniqueN(
    OVERALL_FINAL$loc_num
  ) !=
    140L
) {

  stop(
    paste0(
      "overall样本不是140 monitors；实际=",
      uniqueN(
        OVERALL_FINAL$loc_num
      )
    ),
    call. = FALSE
  )
}


# ============================================================
# J6. Figure 3 representative wind
# ============================================================

FIG3A_FINAL <- fread(
  FIG3A_FILE_FINAL,
  showProgress = FALSE
)

FIG3_REQUIRED_FINAL <- c(
  "radius_km",
  "wind_speed_ms",
  "estimate_mBq_m3",
  "ci95_low_mBq_m3",
  "ci95_high_mBq_m3"
)

FIG3_MISS_FINAL <- setdiff(
  FIG3_REQUIRED_FINAL,
  names(
    FIG3A_FINAL
  )
)

if (
  length(
    FIG3_MISS_FINAL
  ) >
    0L
) {

  stop(
    paste0(
      "Figure3A source缺字段：",
      paste(
        FIG3_MISS_FINAL,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

OVERALL_MEAN_WIND_FINAL <- mean(
  OVERALL_FINAL$vel,
  na.rm = TRUE
)

if (
  "is_representative_wind" %in%
    names(
      FIG3A_FINAL
    ) &&
  any(
    as_bool_final(
      FIG3A_FINAL$is_representative_wind
    )
  )
) {

  FIG3_REP_FINAL <- FIG3A_FINAL[
    radius_km ==
      20 &
    as_bool_final(
      is_representative_wind
    )
  ]

} else {

  FIG3_REP_FINAL <- FIG3A_FINAL[
    radius_km ==
      20
  ][
    which.min(
      abs(
        wind_speed_ms -
          OVERALL_MEAN_WIND_FINAL
      )
    )
  ]
}

if (
  nrow(
    FIG3_REP_FINAL
  ) !=
    1L
) {

  stop(
    "无法唯一识别Figure3 20-km representative-wind result。",
    call. = FALSE
  )
}

OVERALL_WIND_FINAL <-
  FIG3_REP_FINAL$wind_speed_ms[1]

cat(
  "Figure 3 representative wind = ",
  sprintf(
    "%.6f",
    OVERALL_WIND_FINAL
  ),
  " m/s\n",
  sep = ""
)


# ============================================================
# J7. Overall model
# ============================================================

OVERALL_Q95_FINAL <- as.numeric(
  quantile(
    OVERALL_FINAL$all_upwind_n_plants_20km,
    probs = 0.95,
    na.rm = TRUE,
    names = FALSE
  )
)

if (
  !is.finite(
    OVERALL_Q95_FINAL
  ) ||
  OVERALL_Q95_FINAL <=
    0
) {

  stop(
    "overall exposure P95 scale无效。",
    call. = FALSE
  )
}

OVERALL_FORMULA_FINAL <- as.formula(
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


fit_overall_final <- function(
  D,
  wind_ref
) {

  d <- copy(D)

  d[
    ,
    exposure_scaled :=
      all_upwind_n_plants_20km /
      OVERALL_Q95_FINAL
  ]

  d[
    ,
    loc_num_factor :=
      droplevels(
        factor(
          loc_num
        )
      )
  ]

  VARS <- c(
    "gross_beta_mBq_m3",
    "exposure_scaled",
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
    "loc_num_factor"
  )

  d <- d[
    complete.cases(
      d[
        ,
        ..VARS
      ]
    )
  ]

  warnings_final <- character()

  fit <- tryCatch(

    withCallingHandlers(

      lmer(
        OVERALL_FORMULA_FINAL,
        data = d,
        REML = TRUE,
        na.action = na.fail,
        control =
          lmerControl(
            optimizer = "bobyqa",
            optCtrl =
              list(
                maxfun = 200000
              ),
            check.conv.singular =
              "ignore",
            check.rankX =
              "message+drop.cols"
          )
      ),

      warning = function(w) {

        warnings_final <<- c(
          warnings_final,
          conditionMessage(
            w
          )
        )

        invokeRestart(
          "muffleWarning"
        )
      }
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
      data.table(
        status = "FAIL",
        n_rows = nrow(d),
        n_monitors =
          uniqueN(
            d$loc_num
          ),
        estimate = NA_real_,
        se = NA_real_,
        ci95_low = NA_real_,
        ci95_high = NA_real_,
        note =
          conditionMessage(
            fit
          )
      )
    )
  }

  b <- lme4::fixef(
    fit
  )

  V <- as.matrix(
    vcov(
      fit
    )
  )

  int_name <- intersect(
    c(
      "exposure_scaled:vel",
      "vel:exposure_scaled"
    ),
    names(b)
  )

  if (
    !"exposure_scaled" %in%
      names(b) ||
    length(
      int_name
    ) !=
      1L
  ) {

    return(
      data.table(
        status = "FAIL",
        n_rows = nrow(d),
        n_monitors =
          uniqueN(
            d$loc_num
          ),
        estimate = NA_real_,
        se = NA_real_,
        ci95_low = NA_real_,
        ci95_high = NA_real_,
        note =
          "Exposure coefficient missing"
      )
    )
  }

  int_name <- int_name[1]

  est_scaled <-
    unname(
      b[
        "exposure_scaled"
      ]
    ) +
    wind_ref *
    unname(
      b[
        int_name
      ]
    )

  var_scaled <-
    V[
      "exposure_scaled",
      "exposure_scaled"
    ] +
    wind_ref^2 *
    V[
      int_name,
      int_name
    ] +
    2 *
    wind_ref *
    V[
      "exposure_scaled",
      int_name
    ]

  if (
    !is.finite(
      var_scaled
    ) ||
    var_scaled <
      0
  ) {

    return(
      data.table(
        status = "FAIL",
        n_rows = nrow(d),
        n_monitors =
          uniqueN(
            d$loc_num
          ),
        estimate = NA_real_,
        se = NA_real_,
        ci95_low = NA_real_,
        ci95_high = NA_real_,
        note =
          "Invalid contrast variance"
      )
    )
  }

  se_scaled <- sqrt(
    var_scaled
  )

  estimate <-
    est_scaled /
    OVERALL_Q95_FINAL

  se <-
    se_scaled /
    OVERALL_Q95_FINAL

  data.table(
    status = "PASS",

    n_rows =
      nrow(d),

    n_monitors =
      uniqueN(
        d$loc_num
      ),

    estimate =
      estimate,

    se =
      se,

    ci95_low =
      estimate -
      1.96 *
      se,

    ci95_high =
      estimate +
      1.96 *
      se,

    note =
      paste(
        unique(
          warnings_final
        ),
        collapse = " | "
      )
  )
}


# ============================================================
# J8. Full overall model validation
# ============================================================

cat(
  "Fitting overall full-sample model...\n"
)

OVERALL_FULL_FINAL <- fit_overall_final(
  OVERALL_FINAL,
  OVERALL_WIND_FINAL
)

if (
  OVERALL_FULL_FINAL$status !=
    "PASS"
) {

  print(
    OVERALL_FULL_FINAL
  )

  stop(
    "overall full model failed。",
    call. = FALSE
  )
}

OVERALL_DIFF_FINAL <- max(
  c(
    abs(
      OVERALL_FULL_FINAL$estimate -
        FIG3_REP_FINAL$estimate_mBq_m3[1]
    ),

    abs(
      OVERALL_FULL_FINAL$ci95_low -
        FIG3_REP_FINAL$ci95_low_mBq_m3[1]
    ),

    abs(
      OVERALL_FULL_FINAL$ci95_high -
        FIG3_REP_FINAL$ci95_high_mBq_m3[1]
    )
  )
)

cat(
  "Overall full model vs Figure 3 max abs difference = ",
  format(
    OVERALL_DIFF_FINAL,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

if (
  !is.finite(
    OVERALL_DIFF_FINAL
  ) ||
  OVERALL_DIFF_FINAL >
    1e-7
) {

  stop(
    "overall full model没有复现当前Figure 3。",
    call. = FALSE
  )
}


# ============================================================
# J9. Overall leave-one-monitor-out with checkpoint
# ============================================================

OVERALL_LOO_FILE_FINAL <- file.path(
  SOURCE_FINAL,
  "SFigureS4_overall_leave_one_monitor_out_results.csv"
)

if (
  file.exists(
    OVERALL_LOO_FILE_FINAL
  )
) {

  OVERALL_LOO_FINAL <- fread(
    OVERALL_LOO_FILE_FINAL,
    showProgress = FALSE
  )

} else {

  OVERALL_LOO_FINAL <- data.table()
}

if (
  nrow(
    OVERALL_LOO_FINAL
  ) >
    0L
) {

  OVERALL_LOO_FINAL[
    ,
    omitted_monitor :=
      clean_chr_final(
        omitted_monitor
      )
  ]

  OVERALL_LOO_FINAL <- unique(
    OVERALL_LOO_FINAL,
    by = "omitted_monitor",
    fromLast = TRUE
  )
}

OVERALL_MON_FINAL <- sort(
  unique(
    OVERALL_FINAL$loc_num
  )
)

DONE_FINAL <- if (
  nrow(
    OVERALL_LOO_FINAL
  ) >
    0L
) {

  OVERALL_LOO_FINAL[
    status ==
      "PASS",
    omitted_monitor
  ]

} else {

  character()
}

TODO_FINAL <- setdiff(
  OVERALL_MON_FINAL,
  DONE_FINAL
)

cat(
  "Overall LOO total = 140; completed = ",
  length(
    DONE_FINAL
  ),
  "; remaining = ",
  length(
    TODO_FINAL
  ),
  "\n",
  sep = ""
)

if (
  length(
    TODO_FINAL
  ) >
    0L
) {

  for (
    ii in seq_along(
      TODO_FINAL
    )
  ) {

    omitted <- TODO_FINAL[ii]

    cat(
      "Overall LOO ",
      ii,
      "/",
      length(
        TODO_FINAL
      ),
      " | omit ",
      omitted,
      "\n",
      sep = ""
    )

    dd <- OVERALL_FINAL[
      loc_num !=
        omitted
    ]

    ee <- fit_overall_final(
      dd,
      OVERALL_WIND_FINAL
    )

    rr <- data.table(
      omitted_monitor =
        omitted,

      wind_ms =
        OVERALL_WIND_FINAL,

      status =
        ee$status,

      n_rows =
        ee$n_rows,

      n_monitors =
        ee$n_monitors,

      estimate =
        ee$estimate,

      se =
        ee$se,

      ci95_low =
        ee$ci95_low,

      ci95_high =
        ee$ci95_high,

      note =
        ee$note
    )

    OVERALL_LOO_FINAL <- rbind(
      OVERALL_LOO_FINAL,
      rr,
      fill = TRUE
    )

    OVERALL_LOO_FINAL <- unique(
      OVERALL_LOO_FINAL,
      by = "omitted_monitor",
      fromLast = TRUE
    )

    fwrite(
      OVERALL_LOO_FINAL,
      OVERALL_LOO_FILE_FINAL
    )
  }
}

if (
  nrow(
    OVERALL_LOO_FINAL
  ) !=
    140L ||
  OVERALL_LOO_FINAL[
    status !=
      "PASS",
    .N
  ] !=
    0L
) {

  stop(
    "overall leave-one-monitor-out未完整得到140个PASS结果。",
    call. = FALSE
  )
}

setorder(
  OVERALL_LOO_FINAL,
  omitted_monitor
)


# ============================================================
# J10. Overall influence summary
# ============================================================

OVERALL_LOO_FINAL[
  ,
  abs_shift :=
    abs(
      estimate -
        OVERALL_FULL_FINAL$estimate
    )
]

OVERALL_INFLUENTIAL_FINAL <- OVERALL_LOO_FINAL[
  which.max(
    abs_shift
  )
]

OVERALL_SUMMARY_FINAL <- data.table(
  analysis =
    "Overall upwind plant count",

  full_n =
    OVERALL_FULL_FINAL$n_rows,

  full_monitors =
    OVERALL_FULL_FINAL$n_monitors,

  wind_ms =
    OVERALL_WIND_FINAL,

  full_estimate =
    OVERALL_FULL_FINAL$estimate,

  full_ci95_low =
    OVERALL_FULL_FINAL$ci95_low,

  full_ci95_high =
    OVERALL_FULL_FINAL$ci95_high,

  loo_min =
    min(
      OVERALL_LOO_FINAL$estimate
    ),

  loo_max =
    max(
      OVERALL_LOO_FINAL$estimate
    ),

  max_abs_shift =
    max(
      OVERALL_LOO_FINAL$abs_shift
    ),

  all_same_sign_as_full =
    all(
      sign(
        OVERALL_LOO_FINAL$estimate
      ) ==
      sign(
        OVERALL_FULL_FINAL$estimate
      )
    ),

  n_ci_above_zero =
    sum(
      OVERALL_LOO_FINAL$ci95_low >
        0
    ),

  n_ci_below_zero =
    sum(
      OVERALL_LOO_FINAL$ci95_high <
        0
    ),

  n_ci_cross_zero =
    sum(
      OVERALL_LOO_FINAL$ci95_low <=
        0 &
      OVERALL_LOO_FINAL$ci95_high >=
        0
    ),

  most_influential_monitor =
    OVERALL_INFLUENTIAL_FINAL$omitted_monitor,

  estimate_after_omission =
    OVERALL_INFLUENTIAL_FINAL$estimate
)

fwrite(
  OVERALL_SUMMARY_FINAL,
  file.path(
    TABLE_FINAL,
    "Supplementary_Figure_S4_overall_leave_one_out_summary.csv"
  )
)


# ============================================================
# J11. Reference-style overall S4 panel
# ============================================================

OVERALL_LOO_FINAL[
  ,
  monitor_factor :=
    factor(
      omitted_monitor,
      levels =
        rev(
          omitted_monitor
        )
    )
]

P_OVERALL_FINAL <- ggplot(
  OVERALL_LOO_FINAL,
  aes(
    x = estimate,
    y = monitor_factor
  )
) +

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.35
  ) +

  geom_segment(
    aes(
      x = ci95_low,
      xend = ci95_high,
      y = monitor_factor,
      yend = monitor_factor
    ),
    linewidth = 0.28,
    alpha = 0.65
  ) +

  geom_point(
    size = 0.9
  ) +

  geom_vline(
    xintercept =
      OVERALL_FULL_FINAL$estimate,
    linewidth = 0.8,
    color = "red"
  ) +

  geom_point(
    data =
      OVERALL_INFLUENTIAL_FINAL,
    aes(
      x = estimate,
      y =
        factor(
          omitted_monitor,
          levels =
            levels(
              OVERALL_LOO_FINAL$monitor_factor
            )
        )
    ),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    size = 2.7,
    stroke = 0.8
  ) +

  labs(
    title =
      "A. Overall upwind power-plant count",

    subtitle =
      paste0(
        "20 km, 90° upwind; 140-monitor overall analysis; ",
        "effect evaluated at ",
        sprintf(
          "%.2f",
          OVERALL_WIND_FINAL
        ),
        " m/s"
      ),

    x =
      "Change in gross beta per additional upwind power plant (mBq/m³)",

    y =
      "RadNet monitor omitted"
  ) +

  theme_bw(
    base_size = 9
  ) +

  theme(
    panel.grid.major.y =
      element_blank(),

    panel.grid.minor =
      element_blank(),

    axis.text.y =
      element_text(
        size = 3.2
      ),

    plot.title =
      element_text(
        face = "bold",
        size = 11
      ),

    plot.subtitle =
      element_text(
        size = 8.5
      )
  )


# ============================================================
# J12. Fuel-specific reference-style forest panels
# ============================================================

FUEL_NAMES_FINAL <- c(
  coal = "Coal",
  gas = "Natural gas",
  oil = "Oil",
  biomass = "Biomass"
)

METRIC_NAMES_FINAL <- c(
  count = "Active plant count",
  heat = "Actual daily Heat Input"
)

FULL_KEY_FINAL <- FUEL_FULL_FINAL[
  ,
  .(
    fuel,
    metric,
    full_estimate =
      estimate,
    full_ci_low =
      ci95_low,
    full_ci_high =
      ci95_high,
    full_wind_ms =
      wind_ms
  )
]

FUEL_PLOT_FINAL <- merge(
  FUEL_LOO_FINAL,
  FULL_KEY_FINAL,
  by = c(
    "fuel",
    "metric"
  ),
  all.x = TRUE,
  sort = FALSE
)

FUEL_PLOT_FINAL[
  ,
  omitted_index :=
    frank(
      omitted_monitor,
      ties.method = "dense"
    ),
  by = .(
    fuel,
    metric
  )
]

FUEL_PLOT_FINAL[
  ,
  abs_shift :=
    abs(
      estimate -
        full_estimate
    )
]

FUEL_INFLUENTIAL_FINAL <- FUEL_PLOT_FINAL[
  ,
  .SD[
    which.max(
      abs_shift
    )
  ],
  by = .(
    fuel,
    metric
  )
]


make_fuel_forest_final <- function(
  fuel_i,
  metric_i
) {

  Z <- FUEL_PLOT_FINAL[
    fuel ==
      fuel_i &
    metric ==
      metric_i
  ]

  R <- FULL_KEY_FINAL[
    fuel ==
      fuel_i &
    metric ==
      metric_i
  ]

  I <- FUEL_INFLUENTIAL_FINAL[
    fuel ==
      fuel_i &
    metric ==
      metric_i
  ]

  ggplot(
    Z,
    aes(
      x = estimate,
      y = omitted_index
    )
  ) +

    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.3
    ) +

    geom_segment(
      aes(
        x = ci95_low,
        xend = ci95_high,
        y = omitted_index,
        yend = omitted_index
      ),
      linewidth = 0.20,
      alpha = 0.45
    ) +

    geom_point(
      size = 0.65,
      alpha = 0.75
    ) +

    geom_vline(
      xintercept =
        R$full_estimate,
      linewidth = 0.65,
      color = "red"
    ) +

    geom_point(
      data = I,
      aes(
        x = estimate,
        y = omitted_index
      ),
      inherit.aes = FALSE,
      shape = 21,
      fill = "white",
      size = 2.2,
      stroke = 0.7
    ) +

    labs(
      title =
        paste0(
          FUEL_NAMES_FINAL[
            fuel_i
          ],
          "\n",
          METRIC_NAMES_FINAL[
            metric_i
          ]
        ),

      x = NULL,
      y = NULL
    ) +

    theme_bw(
      base_size = 8
    ) +

    theme(
      panel.grid.minor =
        element_blank(),

      panel.grid.major.y =
        element_blank(),

      plot.title =
        element_text(
          face = "bold",
          size = 8,
          hjust = 0.5,
          lineheight = 1
        ),

      axis.text =
        element_text(
          size = 6.5
        )
    )
}


FUEL_PANEL_LIST_FINAL <- list()

kk_final <- 0L

for (
  metric_i in c(
    "count",
    "heat"
  )
) {

  for (
    fuel_i in c(
      "coal",
      "gas",
      "oil",
      "biomass"
    )
  ) {

    kk_final <- kk_final + 1L

    FUEL_PANEL_LIST_FINAL[[kk_final]] <-
      make_fuel_forest_final(
        fuel_i,
        metric_i
      )
  }
}

P_FUEL_FINAL <- wrap_plots(
  FUEL_PANEL_LIST_FINAL,
  ncol = 4
)


# ============================================================
# J13. Final revised S4
# ============================================================

S4_FINAL <- (
  P_OVERALL_FINAL /
  P_FUEL_FINAL
) +

  plot_layout(
    heights = c(
      1.20,
      1.65
    )
  ) +

  plot_annotation(
    title =
      "Supplementary Figure S4. Influence of omitting one RadNet monitor on estimated power-plant associations",

    subtitle =
      paste0(
        "Each point and horizontal line represents the effect estimate and 95% confidence interval ",
        "after omitting all observations from one monitor. ",
        "Red vertical lines indicate the corresponding full-sample estimates; ",
        "open circles identify the omission producing the largest absolute change."
      ),

    caption =
      paste0(
        "Panel A evaluates overall 20-km 90° upwind power-plant count in the earlier 140-monitor random-intercept analysis. ",
        "Panels B–I evaluate 20-km fuel-specific actual-operation Count and Heat Input in the 133-monitor ",
        "Monitor × year fixed-effects analysis. Panel-specific x-axis scales are used for the fuel-specific analyses."
      ),

    theme =
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 14
          ),

        plot.subtitle =
          element_text(
            size = 9
          ),

        plot.caption =
          element_text(
            size = 8,
            hjust = 0
          )
      )
  )

S4_PNG_FINAL <- file.path(
  FIG_FINAL,
  "Supplementary_Figure_S4_leave_one_monitor_out.png"
)

S4_PDF_FINAL <- file.path(
  FIG_FINAL,
  "Supplementary_Figure_S4_leave_one_monitor_out.pdf"
)

ggsave(
  S4_PNG_FINAL,
  S4_FINAL,
  width = 16,
  height = 18,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  S4_PDF_FINAL,
  S4_FINAL,
  width = 16,
  height = 18,
  units = "in"
)


# ============================================================
# J14. Combined influence summary
# ============================================================

FUEL_SUMMARY_FINAL <- FUEL_PLOT_FINAL[
  ,
  .(
    full_estimate =
      full_estimate[1],

    full_ci95_low =
      full_ci_low[1],

    full_ci95_high =
      full_ci_high[1],

    wind_ms =
      full_wind_ms[1],

    loo_min =
      min(
        estimate
      ),

    loo_max =
      max(
        estimate
      ),

    max_abs_shift =
      max(
        abs_shift
      ),

    all_same_sign_as_full =
      all(
        sign(
          estimate
        ) ==
        sign(
          full_estimate[1]
        )
      ),

    n_ci_above_zero =
      sum(
        ci95_low >
          0
      ),

    n_ci_below_zero =
      sum(
        ci95_high <
          0
      ),

    n_ci_cross_zero =
      sum(
        ci95_low <=
          0 &
        ci95_high >=
          0
      ),

    most_influential_monitor =
      omitted_monitor[
        which.max(
          abs_shift
        )
      ],

    estimate_after_omission =
      estimate[
        which.max(
          abs_shift
        )
      ]
  ),
  by = .(
    fuel,
    metric
  )
]

fwrite(
  FUEL_SUMMARY_FINAL,
  file.path(
    TABLE_FINAL,
    "Supplementary_Figure_S4_fuel_specific_leave_one_out_summary.csv"
  )
)


# ============================================================
# J15. Final redraw audit
# ============================================================

FINAL_AUDIT <- data.table(
  check = c(
    "S2_formal_monitors_133",
    "fuel_LOO_rows_1064",
    "fuel_LOO_all_PASS",
    "overall_sample_rows_558915",
    "overall_monitors_140",
    "overall_full_matches_Figure3",
    "overall_LOO_rows_140",
    "overall_LOO_all_PASS",
    "revised_S2_PNG_exists",
    "revised_S2_PDF_exists",
    "revised_S4_PNG_exists",
    "revised_S4_PDF_exists"
  ),

  observed = c(
    as.character(
      nrow(
        M_FINAL
      )
    ),

    as.character(
      nrow(
        FUEL_LOO_FINAL
      )
    ),

    as.character(
      FUEL_LOO_FINAL[
        status !=
          "PASS",
        .N
      ]
    ),

    as.character(
      nrow(
        OVERALL_FINAL
      )
    ),

    as.character(
      uniqueN(
        OVERALL_FINAL$loc_num
      )
    ),

    format(
      OVERALL_DIFF_FINAL,
      scientific = TRUE
    ),

    as.character(
      nrow(
        OVERALL_LOO_FINAL
      )
    ),

    as.character(
      OVERALL_LOO_FINAL[
        status !=
          "PASS",
        .N
      ]
    ),

    as.character(
      file.exists(
        S2_PNG_FINAL
      )
    ),

    as.character(
      file.exists(
        S2_PDF_FINAL
      )
    ),

    as.character(
      file.exists(
        S4_PNG_FINAL
      )
    ),

    as.character(
      file.exists(
        S4_PDF_FINAL
      )
    )
  ),

  expected = c(
    "133",
    "1064",
    "0",
    "558915",
    "140",
    "<=1e-7",
    "140",
    "0",
    "TRUE",
    "TRUE",
    "TRUE",
    "TRUE"
  )
)

FINAL_AUDIT[
  ,
  status := c(
    ifelse(
      nrow(
        M_FINAL
      ) ==
        133L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      nrow(
        FUEL_LOO_FINAL
      ) ==
        1064L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      FUEL_LOO_FINAL[
        status !=
          "PASS",
        .N
      ] ==
        0L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      nrow(
        OVERALL_FINAL
      ) ==
        558915L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      uniqueN(
        OVERALL_FINAL$loc_num
      ) ==
        140L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      is.finite(
        OVERALL_DIFF_FINAL
      ) &&
      OVERALL_DIFF_FINAL <=
        1e-7,
      "PASS",
      "FAIL"
    ),

    ifelse(
      nrow(
        OVERALL_LOO_FINAL
      ) ==
        140L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      OVERALL_LOO_FINAL[
        status !=
          "PASS",
        .N
      ] ==
        0L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      file.exists(
        S2_PNG_FINAL
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      file.exists(
        S2_PDF_FINAL
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      file.exists(
        S4_PNG_FINAL
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      file.exists(
        S4_PDF_FINAL
      ),
      "PASS",
      "FAIL"
    )
  )
]

fwrite(
  FINAL_AUDIT,
  file.path(
    AUDIT_FINAL,
    "05_final_revised_S2_S4_audit.csv"
  )
)

cat(
  "\n============================================================\n",
  "OVERALL LEAVE-ONE-MONITOR-OUT SUMMARY\n",
  "============================================================\n",
  sep = ""
)

print(
  OVERALL_SUMMARY_FINAL,
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "FUEL-SPECIFIC MOST INFLUENTIAL OMISSIONS\n",
  "============================================================\n",
  sep = ""
)

print(
  FUEL_INFLUENTIAL_FINAL[
    ,
    .(
      fuel,
      metric,
      omitted_monitor,
      full_estimate,
      estimate_after_omission =
        estimate,
      abs_shift,
      ci95_low,
      ci95_high
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "FINAL REVISED S2 / S4 AUDIT\n",
  "============================================================\n",
  sep = ""
)

print(
  FINAL_AUDIT,
  nrows = Inf
)

N_FINAL_FAIL <- FINAL_AUDIT[
  status ==
    "FAIL",
  .N
]

cat(
  "\n============================================================\n",
  ifelse(
    N_FINAL_FAIL ==
      0L,
    "07H19 FINAL REVISED S2/S4: PASS\n",
    "07H19 FINAL REVISED S2/S4: FAIL\n"
  ),
  "Hard FAIL = ",
  N_FINAL_FAIL,
  "\n",
  "S2 = ",
  S2_PNG_FINAL,
  "\n",
  "S4 = ",
  S4_PNG_FINAL,
  "\n",
  "============================================================\n",
  sep = ""
)

if (
  N_FINAL_FAIL >
    0L
) {

  stop(
    "最终S2/S4审计存在FAIL。",
    call. = FALSE
  )
}

# >>> 07H19 FINAL REVISED S2/S4 REDRAW END >>>

# >>> 07H19 PUBLICATION REDRAW S2 S4 START >>>

cat(
  "\n============================================================\n",
  "07H19 PUBLICATION REDRAW — S2 / S4\n",
  "No statistical models will be refitted\n",
  "============================================================\n",
  sep = ""
)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(sf)
  library(maps)
  library(scales)
})

ROOT_PUB <- "/Users/hudson/Desktop/USA 辐射表征研究"

HROOT_PUB <- file.path(
  ROOT_PUB,
  "发电站数据处理",
  "06_07H_燃料特异CAMPD逐日暴露与模型"
)

H12_PUB <- file.path(
  HROOT_PUB,
  "07H12_论文正式图表总整理"
)

OUT_PUB <- file.path(
  HROOT_PUB,
  "07H19_补充材料_SFig2_SFig4_TableS2_S3_S4"
)

FIG_PUB <- file.path(
  OUT_PUB,
  "figures"
)

SOURCE_PUB <- file.path(
  OUT_PUB,
  "source_data"
)

TABLE_PUB <- file.path(
  OUT_PUB,
  "tables"
)

AUDIT_PUB <- file.path(
  OUT_PUB,
  "audit"
)

dir.create(
  FIG_PUB,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  SOURCE_PUB,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  TABLE_PUB,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  AUDIT_PUB,
  recursive = TRUE,
  showWarnings = FALSE
)

clean_chr_pub <- function(x) {

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

safe_num_pub <- function(x) {

  suppressWarnings(
    as.numeric(x)
  )
}

as_bool_pub <- function(x) {

  if (
    is.logical(x)
  ) {

    out <- x

  } else if (
    is.numeric(x) ||
      is.integer(x)
  ) {

    out <- x != 0

  } else {

    z <- tolower(
      trimws(
        as.character(x)
      )
    )

    out <- z %in% c(
      "true",
      "t",
      "1",
      "yes",
      "y"
    )
  }

  out[
    is.na(out)
  ] <- FALSE

  out
}


# ============================================================
# A. S2
# ============================================================

cat(
  "\n============================================================\n",
  "A. FINAL PUBLICATION S2\n",
  "============================================================\n",
  sep = ""
)

QDIR_PUB <- file.path(
  ROOT_PUB,
  "发电站数据处理",
  "04_07F_正式统计分析数据构建",
  "quasi_daily_yearly"
)

QMON_PUB <- rbindlist(
  lapply(
    2001:2024,
    function(yy) {

      f <- file.path(
        QDIR_PUB,
        paste0(
          "07F_quasi_daily_analysis_base_",
          yy,
          ".csv.gz"
        )
      )

      if (
        !file.exists(f)
      ) {

        stop(
          paste0(
            "缺少07F quasi-daily：",
            f
          ),
          call. = FALSE
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
          clean_chr_pub(
            loc_num
          )
      ]

      z[
        ,
        flag_07f_primary_quasi_sample :=
          as_bool_pub(
            flag_07f_primary_quasi_sample
          )
      ]

      z
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

ELIGIBLE_PUB <- sort(
  unique(
    QMON_PUB[
      flag_07f_primary_quasi_sample %in%
        TRUE,
      loc_num
    ]
  )
)

rm(
  QMON_PUB
)

if (
  length(
    ELIGIBLE_PUB
  ) !=
    133L
) {

  stop(
    paste0(
      "Formal monitor universe不是133；实际=",
      length(
        ELIGIBLE_PUB
      )
    ),
    call. = FALSE
  )
}


MAP_DIR_PUB <- file.path(
  H12_PUB,
  "05_绘图数据"
)

PLANT_FILE_PUB <- file.path(
  MAP_DIR_PUB,
  "Figure1_powerplant_locations_by_type.csv"
)

MONITOR_FILE_PUB <- file.path(
  MAP_DIR_PUB,
  "Figure1_radnet_monitor_locations.csv"
)

P_PUB <- fread(
  PLANT_FILE_PUB,
  showProgress = FALSE
)

M_PUB <- fread(
  MONITOR_FILE_PUB,
  showProgress = FALSE
)

P_PUB[
  ,
  `:=`(
    plant_longitude =
      safe_num_pub(
        plant_longitude
      ),

    plant_latitude =
      safe_num_pub(
        plant_latitude
      )
  )
]

M_PUB[
  ,
  `:=`(
    monitor_longitude =
      safe_num_pub(
        monitor_longitude
      ),

    monitor_latitude =
      safe_num_pub(
        monitor_latitude
      )
  )
]

P_PUB <- P_PUB[
  is.finite(
    plant_longitude
  ) &
    is.finite(
      plant_latitude
    )
]

M_PUB <- M_PUB[
  is.finite(
    monitor_longitude
  ) &
    is.finite(
      monitor_latitude
    )
]

if (
  !"loc_num" %in%
    names(
      M_PUB
    )
) {

  stop(
    "Figure1 monitor source缺少loc_num。",
    call. = FALSE
  )
}

M_PUB[
  ,
  loc_num :=
    clean_chr_pub(
      loc_num
    )
]

M_PUB <- M_PUB[
  loc_num %in%
    ELIGIBLE_PUB
]

if (
  uniqueN(
    M_PUB$loc_num
  ) !=
    133L
) {

  stop(
    paste0(
      "S2地图monitor不是133；实际=",
      uniqueN(
        M_PUB$loc_num
      )
    ),
    call. = FALSE
  )
}


assign_region_pub <- function(lon) {

  fifelse(
    lon < -105,
    "A. Western CONUS",
    fifelse(
      lon < -90,
      "B. Central CONUS",
      "C. Eastern CONUS"
    )
  )
}

M_PUB[
  ,
  region :=
    assign_region_pub(
      monitor_longitude
    )
]


P_SF_PUB <- st_as_sf(
  P_PUB,
  coords = c(
    "plant_longitude",
    "plant_latitude"
  ),
  crs = 4326,
  remove = FALSE
)

M_SF_PUB <- st_as_sf(
  M_PUB,
  coords = c(
    "monitor_longitude",
    "monitor_latitude"
  ),
  crs = 4326,
  remove = FALSE
)

P_5070_PUB <- st_transform(
  P_SF_PUB,
  5070
)

M_5070_PUB <- st_transform(
  M_SF_PUB,
  5070
)

B20_PUB <- st_transform(
  st_buffer(
    M_5070_PUB,
    dist = 20000
  ),
  4326
)

B50_PUB <- st_transform(
  st_buffer(
    M_5070_PUB,
    dist = 50000
  ),
  4326
)

US_RAW_PUB <- maps::map(
  "state",
  plot = FALSE,
  fill = TRUE
)

US_SF_PUB <- st_as_sf(
  US_RAW_PUB
)

if (
  is.na(
    st_crs(
      US_SF_PUB
    )
  )
) {

  US_SF_PUB <- st_set_crs(
    US_SF_PUB,
    4326
  )
}


REGIONS_PUB <- c(
  "A. Western CONUS",
  "B. Central CONUS",
  "C. Eastern CONUS"
)

REGIONAL_OBJECTS_PUB <- list()

REGIONAL_COUNTS_PUB <- data.table()

for (
  reg in REGIONS_PUB
) {

  MM_5070 <- M_5070_PUB[
    M_5070_PUB$region ==
      reg,
  ]

  MM <- M_SF_PUB[
    M_SF_PUB$region ==
      reg,
  ]

  BB20 <- B20_PUB[
    B20_PUB$region ==
      reg,
  ]

  BB50 <- B50_PUB[
    B50_PUB$region ==
      reg,
  ]

  if (
    nrow(
      MM
    ) ==
      0L
  ) {

    stop(
      paste0(
        "区域无monitor：",
        reg
      ),
      call. = FALSE
    )
  }

  PLANT_MATCH <- lengths(
    st_is_within_distance(
      P_5070_PUB,
      MM_5070,
      dist = 50000
    )
  ) > 0L

  PP <- P_SF_PUB[
    PLANT_MATCH,
  ]

  REGIONAL_COUNTS_PUB <- rbind(
    REGIONAL_COUNTS_PUB,

    data.table(
      region =
        reg,

      n_monitors =
        nrow(
          MM
        ),

      n_unique_plants_within_50km =
        nrow(
          PP
        )
    )
  )

  if (
    nrow(PP) > 0L
  ) {

    GEOM_FOR_BBOX <- c(
      st_geometry(MM),
      st_geometry(PP)
    )

  } else {

    GEOM_FOR_BBOX <- st_geometry(MM)
  }

  BB <- st_bbox(
    GEOM_FOR_BBOX
  )

  XSPAN <- as.numeric(
    BB["xmax"] -
      BB["xmin"]
  )

  YSPAN <- as.numeric(
    BB["ymax"] -
      BB["ymin"]
  )

  XPAD <- max(
    0.9,
    0.07 *
      XSPAN
  )

  YPAD <- max(
    0.7,
    0.07 *
      YSPAN
  )

  XLIM <- c(
    max(
      -125,
      as.numeric(
        BB["xmin"]
      ) -
        XPAD
    ),

    min(
      -66,
      as.numeric(
        BB["xmax"]
      ) +
        XPAD
    )
  )

  YLIM <- c(
    max(
      24,
      as.numeric(
        BB["ymin"]
      ) -
        YPAD
    ),

    min(
      50,
      as.numeric(
        BB["ymax"]
      ) +
        YPAD
    )
  )

  REGIONAL_OBJECTS_PUB[[reg]] <- list(
    monitors =
      MM,

    plants =
      PP,

    buffer20 =
      BB20,

    buffer50 =
      BB50,

    xlim =
      XLIM,

    ylim =
      YLIM
  )
}


fwrite(
  REGIONAL_COUNTS_PUB,
  file.path(
    SOURCE_PUB,
    "SFigureS2_regional_counts.csv"
  )
)


make_S2_panel_pub <- function(reg) {

  OBJ <- REGIONAL_OBJECTS_PUB[[reg]]

  CC <- REGIONAL_COUNTS_PUB[
    region ==
      reg
  ]

  XLIM <- OBJ$xlim
  YLIM <- OBJ$ylim

  DX <- diff(
    XLIM
  )

  DY <- diff(
    YLIM
  )

  LABEL_DATA <- data.frame(
    x =
      XLIM[1] +
      0.025 *
      DX,

    y =
      YLIM[2] -
      0.025 *
      DY,

    label =
      paste0(
        CC$n_monitors,
        " RadNet monitors\n",
        format(
          CC$n_unique_plants_within_50km,
          big.mark = ","
        ),
        " unique power plants within 50 km"
      )
  )

  ggplot() +

    geom_sf(
      data = US_SF_PUB,
      fill = "white",
      color = "grey82",
      linewidth = 0.25
    ) +

    geom_sf(
      data = OBJ$buffer50,
      fill = NA,
      color = "grey58",
      linewidth = 0.38,
      linetype = "dashed"
    ) +

    geom_sf(
      data = OBJ$buffer20,
      fill = NA,
      color = "grey20",
      linewidth = 0.48
    ) +

    geom_sf(
      data = OBJ$plants,
      color = "grey48",
      size = 0.48,
      alpha = 0.34
    ) +

    geom_sf(
      data = OBJ$monitors,
      shape = 17,
      color = "black",
      size = 1.7
    ) +

    geom_label(
      data =
        LABEL_DATA,

      aes(
        x = x,
        y = y,
        label = label
      ),

      inherit.aes = FALSE,

      hjust = 0,
      vjust = 1,

      size = 2.65,

      lineheight = 1.03,

      fill = "white",

      color = "black",

      label.size = 0.15,

      label.padding =
        unit(
          0.11,
          "lines"
        )
    ) +

    coord_sf(
      xlim = XLIM,
      ylim = YLIM,
      expand = FALSE
    ) +

    labs(
      title =
        reg
    ) +

    theme_void(
      base_size = 10
    ) +

    theme(
      plot.title =
        element_text(
          face = "bold",
          size = 11,
          hjust = 0
        ),

      panel.border =
        element_rect(
          color = "grey75",
          fill = NA,
          linewidth = 0.45
        ),

      plot.margin =
        margin(
          5,
          6,
          5,
          6
        )
    )
}


S2_A_PUB <- make_S2_panel_pub(
  "A. Western CONUS"
)

S2_B_PUB <- make_S2_panel_pub(
  "B. Central CONUS"
)

S2_C_PUB <- make_S2_panel_pub(
  "C. Eastern CONUS"
)

S2_PUBLICATION <- (
  S2_A_PUB |
    S2_B_PUB |
    S2_C_PUB
) +

  plot_layout(
    widths = c(
      1,
      1,
      1
    )
  ) +

  plot_annotation(
    title =
      "Supplementary Figure S2. Regional distribution of RadNet monitors and nearby power plants",

    subtitle =
      paste0(
        "Solid circles indicate 20-km monitor-centered buffers; dashed circles indicate 50-km buffers.\n",
        "Black triangles represent RadNet monitors; grey points represent unique nearby power plants. ",
        "Regional divisions are descriptive and are not statistical analysis strata."
      ),

    caption =
      paste0(
        "Regional plant counts represent unique power plants located within 50 km of at least one RadNet monitor shown in that panel. ",
        "A plant may contribute to more than one descriptive panel if it lies within 50 km of monitors on opposite sides of a regional boundary."
      ),

    theme =
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 14
          ),

        plot.subtitle =
          element_text(
            size = 9.5,
            lineheight = 1.08
          ),

        plot.caption =
          element_text(
            size = 7.5,
            hjust = 0
          )
      )
  )


S2_PNG_PUB <- file.path(
  FIG_PUB,
  "Supplementary_Figure_S2_regional_exposure_geometry.png"
)

S2_PDF_PUB <- file.path(
  FIG_PUB,
  "Supplementary_Figure_S2_regional_exposure_geometry.pdf"
)

ggsave(
  S2_PNG_PUB,
  S2_PUBLICATION,
  width = 18,
  height = 7.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  S2_PDF_PUB,
  S2_PUBLICATION,
  width = 18,
  height = 7.2,
  units = "in",
  bg = "white"
)

cat(
  "S2 publication redraw：PASS\n"
)


# ============================================================
# B. S4
# ============================================================

cat(
  "\n============================================================\n",
  "B. FINAL PUBLICATION S4\n",
  "============================================================\n",
  sep = ""
)

FUEL_LOO_FILE_PUB <- file.path(
  SOURCE_PUB,
  "SFigureS4_leave_one_monitor_out_results.csv"
)

FUEL_FULL_FILE_PUB <- file.path(
  SOURCE_PUB,
  "SFigureS4_full_sample_reference_20km.csv"
)

OVERALL_LOO_FILE_PUB <- file.path(
  SOURCE_PUB,
  "SFigureS4_overall_leave_one_monitor_out_results.csv"
)

FIG3A_FILE_PUB <- file.path(
  H12_PUB,
  "05_绘图数据",
  "Figure3A_20_150km_wind1_5_and_representative.csv"
)

REQ_S4_FILES_PUB <- c(
  FUEL_LOO_FILE_PUB,
  FUEL_FULL_FILE_PUB,
  OVERALL_LOO_FILE_PUB,
  FIG3A_FILE_PUB
)

MISS_S4_FILES_PUB <- REQ_S4_FILES_PUB[
  !file.exists(
    REQ_S4_FILES_PUB
  )
]

if (
  length(
    MISS_S4_FILES_PUB
  ) >
    0L
) {

  stop(
    paste0(
      "S4缺少既有结果文件：\n",
      paste(
        MISS_S4_FILES_PUB,
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}


FUEL_LOO_PUB <- fread(
  FUEL_LOO_FILE_PUB,
  showProgress = FALSE
)

FUEL_FULL_PUB <- fread(
  FUEL_FULL_FILE_PUB,
  showProgress = FALSE
)

OVERALL_LOO_PUB <- fread(
  OVERALL_LOO_FILE_PUB,
  showProgress = FALSE
)

FIG3A_PUB <- fread(
  FIG3A_FILE_PUB,
  showProgress = FALSE
)


if (
  nrow(
    FUEL_LOO_PUB
  ) !=
    1064L ||
  FUEL_LOO_PUB[
    status !=
      "PASS",
    .N
  ] !=
    0L
) {

  stop(
    "Fuel-specific LOO不是1064个PASS结果。",
    call. = FALSE
  )
}


if (
  nrow(
    OVERALL_LOO_PUB
  ) !=
    140L ||
  OVERALL_LOO_PUB[
    status !=
      "PASS",
    .N
  ] !=
    0L
) {

  stop(
    "Overall LOO不是140个PASS结果。",
    call. = FALSE
  )
}


FUEL_LOO_PUB[
  ,
  omitted_monitor :=
    clean_chr_pub(
      omitted_monitor
    )
]

OVERALL_LOO_PUB[
  ,
  omitted_monitor :=
    clean_chr_pub(
      omitted_monitor
    )
]


# ============================================================
# B1. Overall full-sample Figure 3 reference
# ============================================================

OVERALL_WIND_PUB <- unique(
  OVERALL_LOO_PUB$wind_ms
)

OVERALL_WIND_PUB <- OVERALL_WIND_PUB[
  is.finite(
    OVERALL_WIND_PUB
  )
]

if (
  length(
    OVERALL_WIND_PUB
  ) !=
    1L
) {

  stop(
    "Overall LOO reference wind不是唯一值。",
    call. = FALSE
  )
}

FIG3_20_PUB <- FIG3A_PUB[
  radius_km ==
    20
]

FIG3_REP_PUB <- FIG3_20_PUB[
  which.min(
    abs(
      wind_speed_ms -
        OVERALL_WIND_PUB
    )
  )
]

if (
  nrow(
    FIG3_REP_PUB
  ) !=
    1L
) {

  stop(
    "无法取得Figure3 20-km representative-wind结果。",
    call. = FALSE
  )
}

if (
  abs(
    FIG3_REP_PUB$wind_speed_ms -
      OVERALL_WIND_PUB
  ) >
    0.02
) {

  stop(
    "Overall LOO wind与Figure3 representative wind不一致。",
    call. = FALSE
  )
}


OVERALL_FULL_EST_PUB <-
  FIG3_REP_PUB$estimate_mBq_m3

OVERALL_FULL_LOW_PUB <-
  FIG3_REP_PUB$ci95_low_mBq_m3

OVERALL_FULL_HIGH_PUB <-
  FIG3_REP_PUB$ci95_high_mBq_m3


# ============================================================
# B2. Overall index + influential monitor
# ============================================================

setorder(
  OVERALL_LOO_PUB,
  omitted_monitor
)

OVERALL_LOO_PUB[
  ,
  omitted_index :=
    seq_len(
      .N
    )
]

OVERALL_LOO_PUB[
  ,
  abs_shift :=
    abs(
      estimate -
        OVERALL_FULL_EST_PUB
    )
]

OVERALL_INF_PUB <- OVERALL_LOO_PUB[
  which.max(
    abs_shift
  )
]

OVERALL_INF_PUB[
  ,
  label_y :=
    pmin(
      135,
      pmax(
        6,
        omitted_index +
          fifelse(
            omitted_index >
              110,
            -8,
            8
          )
      )
    )
]

OVERALL_INF_PUB[
  ,
  label_text :=
    paste0(
      "Monitor ",
      omitted_monitor
    )
]


# ============================================================
# B3. Overall panel
# ============================================================

P_OVERALL_PUB <- ggplot(
  OVERALL_LOO_PUB,
  aes(
    x = estimate,
    y = omitted_index
  )
) +

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.40
  ) +

  geom_segment(
    aes(
      x = ci95_low,
      xend = ci95_high,
      y = omitted_index,
      yend = omitted_index
    ),
    linewidth = 0.30,
    color = "grey65"
  ) +

  geom_point(
    size = 0.95,
    color = "black"
  ) +

  geom_vline(
    xintercept =
      OVERALL_FULL_EST_PUB,
    color = "red",
    linewidth = 0.85
  ) +

  geom_point(
    data =
      OVERALL_INF_PUB,

    aes(
      x = estimate,
      y = omitted_index
    ),

    inherit.aes = FALSE,

    shape = 21,
    fill = "white",
    color = "black",
    size = 3,
    stroke = 0.9
  ) +

  geom_label(
    data =
      OVERALL_INF_PUB,

    aes(
      x = estimate,
      y = label_y,
      label = label_text
    ),

    inherit.aes = FALSE,

    size = 2.7,

    label.size = 0.15,

    fill = "white",

    hjust =
      if (
        OVERALL_INF_PUB$estimate >
          OVERALL_FULL_EST_PUB
      ) {
        1
      } else {
        0
      }
  ) +

  scale_y_continuous(
    breaks = c(
      1,
      35,
      70,
      105,
      140
    ),

    limits = c(
      0,
      141
    )
  ) +

  labs(
    title =
      "A. Overall upwind power-plant count",

    subtitle =
      paste0(
        "20 km, 90° upwind; 140-monitor overall analysis; effect evaluated at ",
        sprintf(
          "%.2f",
          OVERALL_WIND_PUB
        ),
        " m/s"
      ),

    x =
      "Change in gross beta per additional upwind power plant (mBq/m³)",

    y =
      "RadNet monitor omitted (ordered index)"
  ) +

  theme_bw(
    base_size = 9
  ) +

  theme(
    panel.grid.major.y =
      element_blank(),

    panel.grid.minor =
      element_blank(),

    plot.title =
      element_text(
        face = "bold",
        size = 11
      ),

    plot.subtitle =
      element_text(
        size = 8.4
      ),

    axis.text.y =
      element_text(
        size = 7
      )
  )


# ============================================================
# B4. Fuel-specific preparation
# ============================================================

FUEL_LABEL_PUB <- c(
  coal = "Coal",
  gas = "Natural gas",
  oil = "Oil",
  biomass = "Biomass"
)

METRIC_LABEL_PUB <- c(
  count = "Active plant count",
  heat = "Actual daily Heat Input"
)

FULL_KEY_PUB <- FUEL_FULL_PUB[
  ,
  .(
    fuel,
    metric,

    full_estimate =
      estimate,

    full_ci_low =
      ci95_low,

    full_ci_high =
      ci95_high,

    full_wind_ms =
      wind_ms
  )
]

FUEL_PLOT_PUB <- merge(
  FUEL_LOO_PUB,
  FULL_KEY_PUB,
  by = c(
    "fuel",
    "metric"
  ),
  all.x = TRUE,
  sort = FALSE
)


FUEL_PLOT_PUB[
  ,
  omitted_index :=
    frank(
      omitted_monitor,
      ties.method = "dense"
    ),
  by = .(
    fuel,
    metric
  )
]

FUEL_PLOT_PUB[
  ,
  abs_shift :=
    abs(
      estimate -
        full_estimate
    )
]

FUEL_INF_PUB <- FUEL_PLOT_PUB[
  ,
  .SD[
    which.max(
      abs_shift
    )
  ],
  by = .(
    fuel,
    metric
  )
]

FUEL_INF_PUB[
  ,
  label_y :=
    pmin(
      128,
      pmax(
        6,
        omitted_index +
          fifelse(
            omitted_index >
              105,
            -8,
            8
          )
      )
    )
]

FUEL_INF_PUB[
  ,
  label_text :=
    paste0(
      "ID ",
      omitted_monitor
    )
]


# ============================================================
# B5. B-I panel labels
# ============================================================

PANEL_MAP_PUB <- data.table(
  fuel = c(
    "coal",
    "gas",
    "oil",
    "biomass",
    "coal",
    "gas",
    "oil",
    "biomass"
  ),

  metric = c(
    rep(
      "count",
      4
    ),
    rep(
      "heat",
      4
    )
  ),

  panel = LETTERS[
    2:9
  ]
)


make_fuel_panel_pub <- function(
  fuel_i,
  metric_i
) {

  Z <- FUEL_PLOT_PUB[
    fuel ==
      fuel_i &
      metric ==
        metric_i
  ]

  R <- FULL_KEY_PUB[
    fuel ==
      fuel_i &
      metric ==
        metric_i
  ]

  I <- FUEL_INF_PUB[
    fuel ==
      fuel_i &
      metric ==
        metric_i
  ]

  PANEL_LETTER <- PANEL_MAP_PUB[
    fuel ==
      fuel_i &
      metric ==
        metric_i,
    panel
  ]

  if (
    nrow(Z) !=
      133L ||
    nrow(R) !=
      1L ||
    nrow(I) !=
      1L ||
    length(
      PANEL_LETTER
    ) !=
      1L
  ) {

    stop(
      paste0(
        "S4 fuel panel数据不完整：",
        fuel_i,
        " / ",
        metric_i
      ),
      call. = FALSE
    )
  }

  LABEL_HJUST <- if (
    I$estimate >
      R$full_estimate
  ) {
    1
  } else {
    0
  }

  ggplot(
    Z,
    aes(
      x = estimate,
      y = omitted_index
    )
  ) +

    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.32
    ) +

    geom_segment(
      aes(
        x = ci95_low,
        xend = ci95_high,
        y = omitted_index,
        yend = omitted_index
      ),
      linewidth = 0.22,
      color = "grey68"
    ) +

    geom_point(
      size = 0.72,
      color = "black"
    ) +

    geom_vline(
      xintercept =
        R$full_estimate,
      color = "red",
      linewidth = 0.68
    ) +

    geom_point(
      data = I,

      aes(
        x = estimate,
        y = omitted_index
      ),

      inherit.aes = FALSE,

      shape = 21,
      fill = "white",
      color = "black",
      stroke = 0.8,
      size = 2.55
    ) +

    geom_label(
      data = I,

      aes(
        x = estimate,
        y = label_y,
        label = label_text
      ),

      inherit.aes = FALSE,

      size = 2.15,

      label.size = 0.12,

      label.padding =
        unit(
          0.08,
          "lines"
        ),

      fill = "white",

      hjust =
        LABEL_HJUST
    ) +

    scale_y_continuous(
      breaks = c(
        1,
        34,
        67,
        100,
        133
      ),

      limits = c(
        0,
        134
      )
    ) +

    labs(
      title =
        paste0(
          PANEL_LETTER,
          ". ",
          FUEL_LABEL_PUB[
            fuel_i
          ],
          "\n",
          METRIC_LABEL_PUB[
            metric_i
          ]
        ),

      x = NULL,

      y = NULL
    ) +

    theme_bw(
      base_size = 8
    ) +

    theme(
      panel.grid.minor =
        element_blank(),

      panel.grid.major.y =
        element_blank(),

      plot.title =
        element_text(
          face = "bold",
          size = 8.4,
          hjust = 0,
          lineheight = 1.02
        ),

      axis.text =
        element_text(
          size = 6.4
        ),

      plot.margin =
        margin(
          4,
          5,
          4,
          5
        )
    )
}


FUEL_PANEL_PUB <- vector(
  "list",
  8
)

kk_pub <- 0L

for (
  metric_i in c(
    "count",
    "heat"
  )
) {

  for (
    fuel_i in c(
      "coal",
      "gas",
      "oil",
      "biomass"
    )
  ) {

    kk_pub <- kk_pub +
      1L

    FUEL_PANEL_PUB[[kk_pub]] <-
      make_fuel_panel_pub(
        fuel_i,
        metric_i
      )
  }
}


P_FUEL_PUB <- wrap_plots(
  FUEL_PANEL_PUB,
  ncol = 4
)


# ============================================================
# B6. Final S4
# ============================================================

S4_PUBLICATION <- (
  P_OVERALL_PUB /
    P_FUEL_PUB
) +

  plot_layout(
    heights = c(
      1.0,
      1.9
    )
  ) +

  plot_annotation(
    title =
      "Supplementary Figure S4. Influence of omitting one RadNet monitor on estimated power-plant associations",

    subtitle =
      paste0(
        "Each point and horizontal line represents the effect estimate and 95% confidence interval after omitting all observations from one monitor.\n",
        "Red vertical lines indicate full-sample estimates; dashed vertical lines indicate zero; open circles identify the most influential omission."
      ),

    caption =
      paste0(
        "Panel A evaluates overall 20-km 90° upwind power-plant count in the earlier 140-monitor random-intercept analysis. ",
        "Panels B–I evaluate 20-km fuel-specific actual-operation Count and Heat Input in the 133-monitor monitor × year fixed-effects analysis. ",
        "Fuel-specific panels use panel-specific x-axis scales."
      ),

    theme =
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 14
          ),

        plot.subtitle =
          element_text(
            size = 9,
            lineheight = 1.08
          ),

        plot.caption =
          element_text(
            size = 7.8,
            hjust = 0
          )
      )
  )


S4_PNG_PUB <- file.path(
  FIG_PUB,
  "Supplementary_Figure_S4_leave_one_monitor_out.png"
)

S4_PDF_PUB <- file.path(
  FIG_PUB,
  "Supplementary_Figure_S4_leave_one_monitor_out.pdf"
)

ggsave(
  S4_PNG_PUB,
  S4_PUBLICATION,
  width = 16,
  height = 15.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  S4_PDF_PUB,
  S4_PUBLICATION,
  width = 16,
  height = 15.2,
  units = "in",
  bg = "white"
)


# ============================================================
# C. Publication redraw audit
# ============================================================

S2_SIZE <- file.info(
  S2_PNG_PUB
)$size

S4_SIZE <- file.info(
  S4_PNG_PUB
)$size

AUDIT_REDRAW_PUB <- data.table(
  check = c(
    "S2_monitors_133",
    "S2_three_regions",
    "S2_png_exists",
    "S2_pdf_exists",
    "fuel_LOO_1064",
    "fuel_LOO_all_PASS",
    "overall_LOO_140",
    "overall_LOO_all_PASS",
    "S4_png_exists",
    "S4_pdf_exists",
    "S2_png_nonempty",
    "S4_png_nonempty"
  ),

  observed = c(
    as.character(
      uniqueN(
        M_PUB$loc_num
      )
    ),

    as.character(
      nrow(
        REGIONAL_COUNTS_PUB
      )
    ),

    as.character(
      file.exists(
        S2_PNG_PUB
      )
    ),

    as.character(
      file.exists(
        S2_PDF_PUB
      )
    ),

    as.character(
      nrow(
        FUEL_LOO_PUB
      )
    ),

    as.character(
      FUEL_LOO_PUB[
        status !=
          "PASS",
        .N
      ]
    ),

    as.character(
      nrow(
        OVERALL_LOO_PUB
      )
    ),

    as.character(
      OVERALL_LOO_PUB[
        status !=
          "PASS",
        .N
      ]
    ),

    as.character(
      file.exists(
        S4_PNG_PUB
      )
    ),

    as.character(
      file.exists(
        S4_PDF_PUB
      )
    ),

    as.character(
      is.finite(
        S2_SIZE
      ) &&
        S2_SIZE >
        10000
    ),

    as.character(
      is.finite(
        S4_SIZE
      ) &&
        S4_SIZE >
        10000
    )
  ),

  expected = c(
    "133",
    "3",
    "TRUE",
    "TRUE",
    "1064",
    "0",
    "140",
    "0",
    "TRUE",
    "TRUE",
    "TRUE",
    "TRUE"
  )
)

AUDIT_REDRAW_PUB[
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
  AUDIT_REDRAW_PUB,
  file.path(
    AUDIT_PUB,
    "06_publication_redraw_S2_S4_audit.csv"
  )
)


cat(
  "\n============================================================\n",
  "S2 REGIONAL COUNTS\n",
  "============================================================\n",
  sep = ""
)

print(
  REGIONAL_COUNTS_PUB,
  nrows = Inf
)


cat(
  "\n============================================================\n",
  "S4 MOST INFLUENTIAL MONITORS\n",
  "============================================================\n",
  sep = ""
)

cat(
  "\nOverall:\n"
)

print(
  OVERALL_INF_PUB[
    ,
    .(
      omitted_monitor,
      omitted_index,
      estimate,
      ci95_low,
      ci95_high,
      abs_shift
    )
  ]
)

cat(
  "\nFuel-specific:\n"
)

print(
  FUEL_INF_PUB[
    ,
    .(
      fuel,
      metric,
      omitted_monitor,
      omitted_index,
      full_estimate,
      estimate,
      ci95_low,
      ci95_high,
      abs_shift
    )
  ][
    order(
      match(
        metric,
        c(
          "count",
          "heat"
        )
      ),
      match(
        fuel,
        c(
          "coal",
          "gas",
          "oil",
          "biomass"
        )
      )
    )
  ],
  nrows = Inf
)


cat(
  "\n============================================================\n",
  "PUBLICATION REDRAW AUDIT\n",
  "============================================================\n",
  sep = ""
)

print(
  AUDIT_REDRAW_PUB,
  nrows = Inf
)

N_FAIL_REDRAW_PUB <- AUDIT_REDRAW_PUB[
  status ==
    "FAIL",
  .N
]

cat(
  "\n============================================================\n",
  ifelse(
    N_FAIL_REDRAW_PUB ==
      0L,
    "07H19 PUBLICATION S2/S4 REDRAW: PASS\n",
    "07H19 PUBLICATION S2/S4 REDRAW: FAIL\n"
  ),
  "Hard FAIL = ",
  N_FAIL_REDRAW_PUB,
  "\n",
  "S2 = ",
  S2_PNG_PUB,
  "\n",
  "S4 = ",
  S4_PNG_PUB,
  "\n",
  "============================================================\n",
  sep = ""
)

if (
  N_FAIL_REDRAW_PUB >
    0L
) {

  stop(
    "S2/S4 publication redraw存在FAIL。",
    call. = FALSE
  )
}

# >>> 07H19 PUBLICATION REDRAW S2 S4 END >>>
