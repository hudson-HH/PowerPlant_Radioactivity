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
  "ggplot2",
  "patchwork",
  "maps",
  "scales",
  "grid"
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
  library(ggplot2)
  library(patchwork)
  library(maps)
  library(scales)
  library(grid)
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

H10 <- file.path(
  HROOT,
  "07H10_燃料特异站点年份内综合正式分析"
)

H11 <- file.path(
  HROOT,
  "07H11_结果审计与论文结果整理"
)

H11_FREEZE <- file.path(
  H11,
  "07H11_FINAL_FREEZE.md"
)

H11_TABLE_DIR <- file.path(
  H11,
  "01_论文表"
)

# ------------------------------------------------------------
# 前期总体/大类分析的冻结结果
# 注意：
# 这是文章第一阶段分析，不与07H样本偷偷混为同一数据集。
# ------------------------------------------------------------

OLD09 <- file.path(
  ROOT,
  "model_build",
  "09_正式LME模型分析"
)

OLD10 <- file.path(
  ROOT,
  "model_build",
  "10_论文Table与Figure"
)

OLD11 <- file.path(
  ROOT,
  "model_build",
  "11_敏感性验证"
)

OLD09_EFFECTS <- file.path(
  OLD09,
  "02_所有暴露模型效应估计.csv"
)

OLD10_TABLE1 <- file.path(
  OLD10,
  "tables",
  "Table1_描述性统计.csv"
)

OLD11_CLASS <- file.path(
  OLD11,
  "08_combustion_vs_noncombustion_models.csv"
)

# ------------------------------------------------------------
# 地图输入
# ------------------------------------------------------------

STATION_FILE <- file.path(
  ROOT,
  "radnet命名标准化",
  "radnet_181站_CONUS主分析标记与PM覆盖.csv"
)

FUEL_MASTER <- file.path(
  ROOT,
  "powerplant_cleaned",
  "06C_powerplant_analysis_roles_frozen",
  "01_plant_fuel_group_year_analysis_role_master_06C_2001_2024.csv"
)

PLANT_COORD_FALLBACK <- file.path(
  ROOT,
  "powerplant_cleaned",
  "06C1_overall_plant_storage_exclusion_corrected",
  "02_all_plant_nonstorage_proximity_source_06C1_2001_2024.csv"
)

# ------------------------------------------------------------
# 07H10正式结果
# ------------------------------------------------------------

MAIN_EFFECTS_FILE <- file.path(
  H10,
  "01_主分析_监测站年份固定效应",
  "02_全部半径主模型代表风速效应.csv"
)

STRICT_EFFECTS_FILE <- file.path(
  H10,
  "02_严格时间固定效应",
  "02_全部半径严格时间固定效应结果.csv"
)

MAIN_CURVE_FILE <- file.path(
  H10,
  "03_连续风速分析",
  "01_主模型连续风速效应曲线.csv"
)

DIRECTION_EFFECTS_FILE <- file.path(
  H10,
  "04_上风向下风向直接比较",
  "02_上风向下风向及直接差值检验.csv"
)

WEIGHTED_EFFECTS_FILE <- file.path(
  H10,
  "05_距离加权敏感性",
  "02_距离加权全部结果.csv"
)

JOINT_EFFECTS_FILE <- file.path(
  H10,
  "06_数量热输入联合模型",
  "02_Count_Heat联合模型效应.csv"
)

JOINT_COLL_FILE <- file.path(
  H10,
  "06_数量热输入联合模型",
  "03_Count_Heat站点年份内相关性与VIF.csv"
)

PM_EFFECTS_FILE <- file.path(
  H10,
  "07_细颗粒物敏感性",
  "02_不调整PM2.5同样本结果.csv"
)

LAG_EFFECTS_FILE <- file.path(
  H10,
  "08_准逐日滞后分析",
  "02_准逐日滞后全部结果.csv"
)

DECOMP_COMPARE_FILE <- file.path(
  H10,
  "09_站点年份内外分解",
  "03_冻结原模型与站点年份内外方向解释.csv"
)

H10_HARD_FILE <- file.path(
  H10,
  "10_综合汇总",
  "04_07H10_hard_checks.csv"
)

H11_FUEL_SUMMARY_FILE <- file.path(
  H11_TABLE_DIR,
  "11_燃料层面全部分析综合审计.csv"
)

# ============================================================
# 2. Output
# ============================================================

OUT <- file.path(
  HROOT,
  "07H12_论文正式图表总整理"
)

# 同一任务同名覆盖
if (dir.exists(OUT)) {
  unlink(
    OUT,
    recursive = TRUE,
    force = TRUE
  )
}

MAIN_FIG_DIR <- file.path(
  OUT,
  "01_正文图"
)

MAIN_TABLE_DIR <- file.path(
  OUT,
  "02_正文表"
)

SUPP_FIG_DIR <- file.path(
  OUT,
  "03_补充图"
)

SUPP_TABLE_DIR <- file.path(
  OUT,
  "04_补充表"
)

SOURCE_DIR <- file.path(
  OUT,
  "05_绘图数据"
)

AUDIT_DIR <- file.path(
  OUT,
  "06_审计"
)

for (
  d in c(
    OUT,
    MAIN_FIG_DIR,
    MAIN_TABLE_DIR,
    SUPP_FIG_DIR,
    SUPP_TABLE_DIR,
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
# 3. Constants
# ============================================================

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

REP_RADII <- c(
  20,
  50,
  100,
  200
)

FUELS <- c(
  "coal",
  "gas",
  "oil",
  "other_fossil",
  "biomass",
  "solid_waste"
)

CORE_FUELS <- c(
  "coal",
  "gas",
  "oil",
  "biomass"
)

COMBUSTION <- c(
  "coal",
  "gas",
  "oil",
  "other_fossil",
  "biomass",
  "solid_waste"
)

NONCOMBUSTION_MAIN <- c(
  "nuclear",
  "geothermal"
)

NEGATIVE_CONTROL <- c(
  "wind",
  "solar",
  "hydro"
)

FUEL_LABEL <- c(
  coal = "Coal",
  gas = "Natural gas",
  oil = "Oil",
  other_fossil = "Other fossil",
  biomass = "Biomass",
  solid_waste = "Solid waste",
  nuclear = "Nuclear",
  geothermal = "Geothermal",
  wind = "Wind",
  solar = "Solar",
  hydro = "Hydro",
  other = "Other",
  mixed_multiple = "Mixed / multiple fuels"
)

METRIC_LABEL <- c(
  count = "Plant count",
  heat = "Heat Input",
  proximity = "Σ1/d",
  heat_proximity = "ΣHeat/d"
)

BRANCH_LABEL <- c(
  filter_period = "Filter-period",
  quasi_daily = "Quasi-daily"
)

# ============================================================
# 4. Helpers
# ============================================================

stopf <- function(...) {
  stop(
    paste0(...),
    call. = FALSE
  )
}

assert_file <- function(
  f,
  label = basename(f)
) {
  if (!file.exists(f)) {
    stopf(
      "缺少输入文件：",
      label,
      "\n",
      f
    )
  }
}

read_required <- function(
  f,
  cols = NULL
) {
  assert_file(f)

  D <- fread(
    f,
    showProgress = FALSE
  )

  if (!is.null(cols)) {
    miss <- setdiff(
      cols,
      names(D)
    )

    if (length(miss) > 0L) {
      stopf(
        basename(f),
        " 缺少字段：",
        paste(
          miss,
          collapse = ", "
        )
      )
    }
  }

  D
}

safe_num <- function(x) {
  suppressWarnings(
    as.numeric(x)
  )
}

as_flag <- function(x) {

  if (is.logical(x)) {
    return(x)
  }

  z <- tolower(
    trimws(
      as.character(x)
    )
  )

  z %in% c(
    "true",
    "t",
    "1",
    "yes",
    "y"
  )
}

pick_col <- function(
  nm,
  candidates,
  required = TRUE
) {

  z <- intersect(
    candidates,
    nm
  )

  if (length(z) > 0L) {
    return(z[1])
  }

  if (required) {
    stopf(
      "无法找到字段。候选：",
      paste(
        candidates,
        collapse = " | "
      )
    )
  }

  NA_character_
}

save_plot_pair <- function(
  p,
  filename_without_ext,
  width,
  height
) {

  png_file <- paste0(
    filename_without_ext,
    ".png"
  )

  pdf_file <- paste0(
    filename_without_ext,
    ".pdf"
  )

  ggsave(
    png_file,
    p,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  if (capabilities("cairo")) {

    ggsave(
      pdf_file,
      p,
      width = width,
      height = height,
      units = "in",
      device = cairo_pdf
    )

  } else {

    ggsave(
      pdf_file,
      p,
      width = width,
      height = height,
      units = "in",
      device = "pdf"
    )
  }

  invisible(
    c(
      png_file,
      pdf_file
    )
  )
}

base_theme <- theme_bw(
  base_size = 10
) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(
      fill = "white"
    ),
    strip.text = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold"
    ),
    legend.title = element_blank()
  )

add_labels <- function(D) {

  X <- copy(D)

  if ("fuel" %in% names(X)) {
    X[
      ,
      fuel_label :=
        unname(
          FUEL_LABEL[
            fuel
          ]
        )
    ]
  }

  if ("metric" %in% names(X)) {
    X[
      ,
      metric_label :=
        unname(
          METRIC_LABEL[
            metric
          ]
        )
    ]
  }

  if ("branch" %in% names(X)) {
    X[
      ,
      branch_label :=
        unname(
          BRANCH_LABEL[
            branch
          ]
        )
    ]
  }

  X
}

# ============================================================
# 5. Freeze/input audit
# ============================================================

assert_file(
  H11_FREEZE,
  "07H11冻结标记"
)

required_inputs <- c(
  MAIN_EFFECTS_FILE,
  STRICT_EFFECTS_FILE,
  MAIN_CURVE_FILE,
  DIRECTION_EFFECTS_FILE,
  WEIGHTED_EFFECTS_FILE,
  JOINT_EFFECTS_FILE,
  JOINT_COLL_FILE,
  PM_EFFECTS_FILE,
  LAG_EFFECTS_FILE,
  DECOMP_COMPARE_FILE,
  H10_HARD_FILE,
  H11_FUEL_SUMMARY_FILE,
  OLD09_EFFECTS,
  OLD10_TABLE1,
  OLD11_CLASS,
  STATION_FILE,
  FUEL_MASTER
)

for (f in required_inputs) {
  assert_file(f)
}

H10_HARD <- fread(
  H10_HARD_FILE
)

if (
  H10_HARD[
    status != "PASS",
    .N
  ] > 0L
) {
  stopf(
    "07H10硬性检查不是全部PASS。"
  )
}

cat(
  "\n============================================================\n",
  "07H12 INPUT AUDIT: PASS\n",
  "============================================================\n",
  sep = ""
)

# ============================================================
# 6. Read formal results
# ============================================================

MAIN <- add_labels(
  read_required(
    MAIN_EFFECTS_FILE
  )
)

STRICT <- add_labels(
  read_required(
    STRICT_EFFECTS_FILE
  )
)

CURVE <- add_labels(
  read_required(
    MAIN_CURVE_FILE
  )
)

DIRECTION <- add_labels(
  read_required(
    DIRECTION_EFFECTS_FILE
  )
)

WEIGHTED <- add_labels(
  read_required(
    WEIGHTED_EFFECTS_FILE
  )
)

JOINT <- add_labels(
  read_required(
    JOINT_EFFECTS_FILE
  )
)

JOINT_COLL <- read_required(
  JOINT_COLL_FILE
)

PM <- add_labels(
  read_required(
    PM_EFFECTS_FILE
  )
)

LAG <- add_labels(
  read_required(
    LAG_EFFECTS_FILE
  )
)

DECOMP <- add_labels(
  read_required(
    DECOMP_COMPARE_FILE
  )
)

FUEL_SUMMARY <- read_required(
  H11_FUEL_SUMMARY_FILE
)

# ============================================================
# 7. Build map data
# ============================================================

cat(
  "\n构建Figure 1地图数据...\n"
)

# ------------------------------------------------------------
# 7A. RadNet监测站
# ------------------------------------------------------------

ST <- fread(
  STATION_FILE,
  showProgress = FALSE
)

required_station_cols <- c(
  "monitor_latitude",
  "monitor_longitude"
)

if (
  !all(
    required_station_cols %in%
      names(ST)
  )
) {
  stopf(
    "RadNet站点表缺少经纬度字段。现有字段：",
    paste(
      names(ST),
      collapse = " | "
    )
  )
}

if (
  "main_analysis_conus" %in%
    names(ST)
) {
  ST <- ST[
    as_flag(
      main_analysis_conus
    )
  ]
}

ST[
  ,
  `:=`(
    monitor_latitude =
      safe_num(
        monitor_latitude
      ),
    monitor_longitude =
      safe_num(
        monitor_longitude
      )
  )
]

ST <- ST[
  is.finite(
    monitor_latitude
  ) &
  is.finite(
    monitor_longitude
  ) &
  monitor_latitude >=
    24 &
  monitor_latitude <=
    50 &
  monitor_longitude >=
    -125 &
  monitor_longitude <=
    -66
]

if (nrow(ST) == 0L) {
  stopf(
    "CONUS RadNet地图数据为空。"
  )
}

# ------------------------------------------------------------
# 7B. 读取plant-year燃料主表
# ------------------------------------------------------------

FM <- fread(
  FUEL_MASTER,
  showProgress = FALSE,
  colClasses = "character"
)

id_col <- pick_col(
  names(FM),
  c(
    "plant_id",
    "orispl",
    "orispl_code",
    "plant_id_eia"
  )
)

year_col <- pick_col(
  names(FM),
  c(
    "year",
    "analysis_year"
  )
)

fuel_col <- pick_col(
  names(FM),
  c(
    "fuel_group",
    "plant_fuel_group",
    "fuel_group_final",
    "analysis_fuel_group"
  )
)

FM[
  ,
  `:=`(
    plant_id_plot =
      trimws(
        as.character(
          get(id_col)
        )
      ),
    year_plot =
      suppressWarnings(
        as.integer(
          get(year_col)
        )
      ),
    fuel_group_plot_raw =
      trimws(
        as.character(
          get(fuel_col)
        )
      )
  )
]

if (
  FM[
    is.na(plant_id_plot) |
      plant_id_plot == "",
    .N
  ] > 0L
) {
  FM <- FM[
    !is.na(plant_id_plot) &
      plant_id_plot != ""
  ]
}

# ------------------------------------------------------------
# 7C. 先检查燃料主表本身是否已有坐标
#     与此前成功的论文地图脚本保持同样字段兼容范围
# ------------------------------------------------------------

fm_lat_col <- pick_col(
  names(FM),
  c(
    "analysis_latitude_final",
    "plant_latitude",
    "latitude_final",
    "latitude",
    "plant_lat",
    "lat"
  ),
  required = FALSE
)

fm_lon_col <- pick_col(
  names(FM),
  c(
    "analysis_longitude_final",
    "plant_longitude",
    "longitude_final",
    "longitude",
    "plant_lon",
    "lon"
  ),
  required = FALSE
)

if (
  !is.na(fm_lat_col) &&
  !is.na(fm_lon_col)
) {

  cat(
    "燃料主表内发现坐标字段：",
    fm_lat_col,
    " / ",
    fm_lon_col,
    "\n",
    sep = ""
  )

  FM[
    ,
    `:=`(
      plant_latitude =
        safe_num(
          get(fm_lat_col)
        ),
      plant_longitude =
        safe_num(
          get(fm_lon_col)
        )
    )
  ]

} else {

  # ----------------------------------------------------------
  # 7D. 燃料主表没有坐标：
  #     从冻结06C1 non-storage完整发电站源读取
  # ----------------------------------------------------------

  coord_candidates <- unique(
    c(
      PLANT_COORD_FALLBACK,
      file.path(
        ROOT,
        "发电站数据处理",
        "06C1_overall_plant_storage_exclusion_corrected",
        "02_all_plant_nonstorage_proximity_source_06C1_2001_2024.csv"
      )
    )
  )

  coord_candidates <- coord_candidates[
    file.exists(
      coord_candidates
    )
  ]

  if (length(coord_candidates) == 0L) {
    stopf(
      paste0(
        "燃料主表没有坐标，同时没有找到冻结06C1坐标源。\n",
        "已检查：\n",
        paste(
          unique(
            c(
              PLANT_COORD_FALLBACK,
              file.path(
                ROOT,
                "发电站数据处理",
                "06C1_overall_plant_storage_exclusion_corrected",
                "02_all_plant_nonstorage_proximity_source_06C1_2001_2024.csv"
              )
            )
          ),
          collapse = "\n"
        )
      )
    )
  }

  coord_file <- coord_candidates[1]

  cat(
    "使用冻结发电站坐标源：\n",
    coord_file,
    "\n",
    sep = ""
  )

  PC <- fread(
    coord_file,
    showProgress = FALSE,
    colClasses = "character"
  )

  pc_id <- pick_col(
    names(PC),
    c(
      "plant_id",
      "orispl",
      "orispl_code",
      "plant_id_eia"
    )
  )

  pc_lat <- pick_col(
    names(PC),
    c(
      "analysis_latitude_final",
      "plant_latitude",
      "latitude_final",
      "latitude",
      "plant_lat",
      "lat"
    )
  )

  pc_lon <- pick_col(
    names(PC),
    c(
      "analysis_longitude_final",
      "plant_longitude",
      "longitude_final",
      "longitude",
      "plant_lon",
      "lon"
    )
  )

  pc_state <- pick_col(
    names(PC),
    c(
      "state_final_06A3",
      "state_final",
      "state",
      "plant_state",
      "state_abbr"
    ),
    required = FALSE
  )

  cat(
    "识别到坐标字段：",
    pc_lat,
    " / ",
    pc_lon,
    "\n",
    sep = ""
  )

  PC[
    ,
    `:=`(
      plant_id_plot =
        trimws(
          as.character(
            get(pc_id)
          )
        ),
      plant_latitude =
        safe_num(
          get(pc_lat)
        ),
      plant_longitude =
        safe_num(
          get(pc_lon)
        )
    )
  ]

  if (!is.na(pc_state)) {

    PC[
      ,
      plant_state_plot :=
        toupper(
          trimws(
            as.character(
              get(pc_state)
            )
          )
        )
    ]

    lower48_states <- setdiff(
      state.abb,
      c(
        "AK",
        "HI"
      )
    )

    PC <- PC[
      plant_state_plot %in%
        c(
          lower48_states,
          "DC"
        )
    ]
  }

  PC <- PC[
    !is.na(plant_id_plot) &
      plant_id_plot != "" &
    is.finite(
      plant_latitude
    ) &
    is.finite(
      plant_longitude
    )
  ]

  # 同一物理电厂2001–2024可能出现多行。
  # 坐标取该电厂全部记录的中位数，不要求备用坐标文件必须含year。
  PC_COORD <- PC[
    ,
    .(
      plant_latitude =
        median(
          plant_latitude,
          na.rm = TRUE
        ),
      plant_longitude =
        median(
          plant_longitude,
          na.rm = TRUE
        )
    ),
    by =
      plant_id_plot
  ]

  if (
    PC_COORD[
      ,
      uniqueN(
        plant_id_plot
      )
    ] !=
      nrow(PC_COORD)
  ) {
    stopf(
      "发电站坐标聚合后plant_id仍不唯一。"
    )
  }

  FM <- merge(
    FM,
    PC_COORD,
    by =
      "plant_id_plot",
    all.x = TRUE,
    sort = FALSE
  )

  rm(
    PC,
    PC_COORD
  )
  invisible(
    gc()
  )
}

# ------------------------------------------------------------
# 7E. 地图燃料分类与CONUS筛选
# ------------------------------------------------------------

valid_fuels <- c(
  COMBUSTION,
  NONCOMBUSTION_MAIN,
  NEGATIVE_CONTROL,
  "other"
)

FM <- FM[
  fuel_group_plot_raw %in%
    valid_fuels &
  is.finite(
    plant_latitude
  ) &
  is.finite(
    plant_longitude
  ) &
  plant_latitude >=
    24 &
  plant_latitude <=
    50 &
  plant_longitude >=
    -125 &
  plant_longitude <=
    -66
]

if (nrow(FM) == 0L) {
  stopf(
    "合并坐标后没有可用于Figure 1的CONUS发电站记录。"
  )
}

coord_audit <- data.table(
  item = c(
    "fuel_master_rows_with_valid_map_coordinates",
    "unique_physical_plants_with_valid_map_coordinates",
    "unique_fuel_groups_in_map_source"
  ),
  value = c(
    nrow(FM),
    uniqueN(
      FM$plant_id_plot
    ),
    uniqueN(
      FM$fuel_group_plot_raw
    )
  )
)

fwrite(
  coord_audit,
  file.path(
    AUDIT_DIR,
    "00_Figure1_coordinate_input_audit.csv"
  )
)

print(
  coord_audit
)

PLANT_MAP <- FM[
  ,
  {

    u <- sort(
      unique(
        fuel_group_plot_raw
      )
    )

    exact_group <- if (
      length(u) ==
        1L
    ) {
      u
    } else {
      "mixed_multiple"
    }

    broad_group <- if (
      any(
        u %in%
          COMBUSTION
      )
    ) {

      "Combustion-capable"

    } else if (
      any(
        u %in%
          NONCOMBUSTION_MAIN
      )
    ) {

      "Non-combustion major energy"

    } else if (
      any(
        u %in%
          NEGATIVE_CONTROL
      )
    ) {

      "Wind / solar / hydro control"

    } else {

      "Other"
    }

    list(
      plant_latitude =
        median(
          plant_latitude,
          na.rm = TRUE
        ),
      plant_longitude =
        median(
          plant_longitude,
          na.rm = TRUE
        ),
      fuel_group_plot =
        exact_group,
      broad_group =
        broad_group,
      all_fuel_groups =
        paste(
          u,
          collapse = "|"
        )
    )
  },
  by =
    plant_id_plot
]

PLANT_MAP[
  ,
  fuel_label :=
    unname(
      FUEL_LABEL[
        fuel_group_plot
      ]
    )
]

if (
  PLANT_MAP[
    is.na(fuel_label),
    .N
  ] > 0L
) {
  PLANT_MAP[
    is.na(fuel_label),
    fuel_label :=
      "Other"
  ]
}

fwrite(
  PLANT_MAP,
  file.path(
    SOURCE_DIR,
    "Figure1_powerplant_locations_by_type.csv"
  )
)

fwrite(
  ST,
  file.path(
    SOURCE_DIR,
    "Figure1_radnet_monitor_locations.csv"
  )
)

cat(
  "Figure 1地图输入：",
  nrow(PLANT_MAP),
  "个物理发电站；",
  nrow(ST),
  "个RadNet监测站。\n",
  sep = ""
)

# 8. Figure 1
# ============================================================

US_MAP <- ggplot2::map_data(
  "state"
)

fuel_levels <- c(
  "Coal",
  "Natural gas",
  "Oil",
  "Other fossil",
  "Biomass",
  "Solid waste",
  "Nuclear",
  "Geothermal",
  "Wind",
  "Solar",
  "Hydro",
  "Other",
  "Mixed / multiple fuels"
)

PLANT_MAP[
  ,
  fuel_label :=
    factor(
      fuel_label,
      levels =
        fuel_levels
    )
]

p1a <- ggplot() +
  geom_polygon(
    data = US_MAP,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "white",
    color = "grey75",
    linewidth = 0.25
  ) +
  geom_point(
    data = PLANT_MAP,
    aes(
      x = plant_longitude,
      y = plant_latitude,
      color = fuel_label
    ),
    size = 0.75,
    alpha = 0.60
  ) +
  geom_point(
    data = ST,
    aes(
      x = monitor_longitude,
      y = monitor_latitude
    ),
    shape = 17,
    size = 1.7,
    color = "black"
  ) +
  coord_quickmap(
    xlim = c(
      -125,
      -66
    ),
    ylim = c(
      24,
      50
    ),
    expand = FALSE
  ) +
  scale_color_discrete(
    na.translate = FALSE
  ) +
  labs(
    title = "A. Power plants by energy / fuel type",
    x = NULL,
    y = NULL,
    color = NULL
  ) +
  theme_void(
    base_size = 10
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(
      face = "bold"
    )
  )

p1b <- ggplot() +
  geom_polygon(
    data = US_MAP,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "white",
    color = "grey75",
    linewidth = 0.25
  ) +
  geom_point(
    data = PLANT_MAP,
    aes(
      x = plant_longitude,
      y = plant_latitude,
      color = broad_group
    ),
    size = 0.8,
    alpha = 0.55
  ) +
  geom_point(
    data = ST,
    aes(
      x = monitor_longitude,
      y = monitor_latitude
    ),
    shape = 17,
    size = 1.7,
    color = "black"
  ) +
  coord_quickmap(
    xlim = c(
      -125,
      -66
    ),
    ylim = c(
      24,
      50
    ),
    expand = FALSE
  ) +
  labs(
    title = "B. Broad power-plant classes",
    x = NULL,
    y = NULL,
    color = NULL
  ) +
  theme_void(
    base_size = 10
  ) +
  theme(
    legend.position = "right",
    plot.title = element_text(
      face = "bold"
    )
  )

FIG1 <- p1a / p1b +
  plot_annotation(
    title = "Figure 1. Spatial distribution of U.S. power plants and RadNet monitoring stations",
    subtitle = "Triangles indicate RadNet monitoring stations; circles indicate power plants."
  )

save_plot_pair(
  FIG1,
  file.path(
    MAIN_FIG_DIR,
    "Figure1_US_powerplants_RadNet_map"
  ),
  13.5,
  11
)

# ============================================================
# 9. Figure 2 - study design
# ============================================================

flow_boxes <- data.table(
  xmin = c(
    0.3,
    3.0,
    5.7,
    8.4,
    2.0,
    5.1,
    3.6,
    3.6,
    3.6
  ),
  xmax = c(
    2.2,
    4.9,
    7.6,
    10.3,
    4.2,
    7.3,
    6.7,
    6.7,
    6.7
  ),
  ymin = c(
    8.0,
    8.0,
    8.0,
    8.0,
    5.8,
    5.8,
    3.7,
    1.9,
    0.1
  ),
  ymax = c(
    9.2,
    9.2,
    9.2,
    9.2,
    7.0,
    7.0,
    4.9,
    3.1,
    1.3
  ),
  label = c(
    "RadNet\ngross-beta monitoring",
    "Power-plant locations\nand fuel types",
    "EPA CAMPD\nactual daily operation",
    "Meteorology\nwind direction / speed",
    "Monitor–plant\nspatial linkage",
    "Daily directional exposure\n90° upwind / downwind",
    "Stage 1\nOverall and broad-class analysis",
    "Stage 2\nFuel-specific Count + Heat Input",
    "Stage 3\nWithin-monitor-year and\nrobustness / validity analyses"
  )
)

p2 <- ggplot() +
  geom_rect(
    data = flow_boxes,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    fill = "white",
    color = "black",
    linewidth = 0.5
  ) +
  geom_text(
    data = flow_boxes,
    aes(
      x = (
        xmin +
          xmax
      ) /
        2,
      y = (
        ymin +
          ymax
      ) /
        2,
      label = label
    ),
    size = 3.4,
    lineheight = 0.95
  ) +
  annotate(
    "segment",
    x = 1.25,
    xend = 3.1,
    y = 7.9,
    yend = 7.05,
    arrow = arrow(
      length = unit(
        0.15,
        "inches"
      )
    )
  ) +
  annotate(
    "segment",
    x = 3.95,
    xend = 3.1,
    y = 7.9,
    yend = 7.05,
    arrow = arrow(
      length = unit(
        0.15,
        "inches"
      )
    )
  ) +
  annotate(
    "segment",
    x = 6.65,
    xend = 6.2,
    y = 7.9,
    yend = 7.05,
    arrow = arrow(
      length = unit(
        0.15,
        "inches"
      )
    )
  ) +
  annotate(
    "segment",
    x = 9.35,
    xend = 6.2,
    y = 7.9,
    yend = 7.05,
    arrow = arrow(
      length = unit(
        0.15,
        "inches"
      )
    )
  ) +
  annotate(
    "segment",
    x = 4.2,
    xend = 5.1,
    y = 6.4,
    yend = 6.4,
    arrow = arrow(
      length = unit(
        0.15,
        "inches"
      )
    )
  ) +
  annotate(
    "segment",
    x = 6.2,
    xend = 5.15,
    y = 5.75,
    yend = 4.95,
    arrow = arrow(
      length = unit(
        0.15,
        "inches"
      )
    )
  ) +
  annotate(
    "segment",
    x = 5.15,
    xend = 5.15,
    y = 3.65,
    yend = 3.15,
    arrow = arrow(
      length = unit(
        0.15,
        "inches"
      )
    )
  ) +
  annotate(
    "segment",
    x = 5.15,
    xend = 5.15,
    y = 1.85,
    yend = 1.35,
    arrow = arrow(
      length = unit(
        0.15,
        "inches"
      )
    )
  ) +
  coord_cartesian(
    xlim = c(
      0,
      10.7
    ),
    ylim = c(
      0,
      9.6
    ),
    clip = "off"
  ) +
  labs(
    title = "Figure 2. Study design and analytical framework"
  ) +
  theme_void(
    base_size = 11
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

save_plot_pair(
  p2,
  file.path(
    MAIN_FIG_DIR,
    "Figure2_study_design_flowchart"
  ),
  11,
  9
)

# ============================================================
# 10. Table 1
# ============================================================

OLD_T1 <- fread(
  OLD10_TABLE1
)

old_desc <- melt(
  OLD_T1,
  id.vars = intersect(
    c(
      "variable",
      "Variable"
    ),
    names(OLD_T1)
  ),
  variable.name =
    "analysis_sample",
  value.name =
    "value"
)

old_desc[
  ,
  section :=
    "A. Earlier overall-analysis sample"
]

old_desc[
  ,
  note :=
    "Frozen first-stage overall/broad-class analysis"
]

get_h10_observed <- function(check_name) {
  z <- H10_HARD[
    check ==
      check_name,
    observed
  ]

  if (length(z) != 1L) {
    return(
      NA_character_
    )
  }

  as.character(
    z
  )
}

current_sample <- data.table(
  section =
    "B. Updated fuel-specific formal samples",
  variable = c(
    "quasi_daily_rows",
    "quasi_daily_monitors",
    "filter_period_rows",
    "filter_period_monitors"
  ),
  Variable = c(
    "Quasi-daily observations (n)",
    "Quasi-daily monitors (n)",
    "Filter-period observations (n)",
    "Filter-period monitors (n)"
  ),
  analysis_sample = c(
    "07H formal quasi-daily",
    "07H formal quasi-daily",
    "07H formal filter-period",
    "07H formal filter-period"
  ),
  value = c(
    get_h10_observed(
      "quasi_primary_rows_426657"
    ),
    get_h10_observed(
      "quasi_primary_monitors_133"
    ),
    get_h10_observed(
      "filter_primary_rows_107407"
    ),
    get_h10_observed(
      "filter_primary_monitors_133"
    )
  ),
  note =
    "Updated corrected 07F/07H formal sample"
)

plant_stats <- FM[
  fuel_group_plot_raw %in%
    valid_fuels,
  .(
    physical_plants =
      uniqueN(
        plant_id_plot
      ),
    plant_years =
      uniqueN(
        paste(
          plant_id_plot,
          year_plot,
          sep = "_"
        )
      )
  ),
  by =
    fuel_group_plot_raw
]

plant_stats[
  ,
  `:=`(
    section =
      "C. Power-plant characteristics",
    variable =
      paste0(
        "plants_",
        fuel_group_plot_raw
      ),
    Variable =
      unname(
        FUEL_LABEL[
          fuel_group_plot_raw
        ]
      ),
    analysis_sample =
      "2001–2024 plant master",
    value =
      paste0(
        physical_plants,
        " plants; ",
        plant_years,
        " plant-years"
      ),
    note =
      "Unique physical plants and plant-years"
  )
]

TABLE1 <- rbindlist(
  list(
    old_desc[
      ,
      .(
        section,
        variable,
        Variable,
        analysis_sample,
        value,
        note
      )
    ],
    current_sample,
    plant_stats[
      ,
      .(
        section,
        variable,
        Variable,
        analysis_sample,
        value,
        note
      )
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  TABLE1,
  file.path(
    MAIN_TABLE_DIR,
    "Table1_study_population_and_powerplant_characteristics.csv"
  )
)

# ============================================================
# 11. Figure 3 + Table 2
# ============================================================

OLD_EFFECTS <- fread(
  OLD09_EFFECTS
)

CLASS_MODELS <- fread(
  OLD11_CLASS
)

overall_row <- OLD_EFFECTS[
  model_id ==
    "overall_upwind_count_20km" &
  effect_type ==
    "per_raw_unit" &
  wind_speed_ms ==
    1
]

if (
  nrow(overall_row) !=
    1L
) {
  stopf(
    "无法唯一识别前期总体20 km上风向数量主效应。"
  )
}

FIG3_DATA <- data.table(
  category =
    "All power plants",
  analysis =
    "Overall model",
  estimate =
    overall_row$estimate_mBq_m3,
  ci_low =
    overall_row$ci95_low_mBq_m3,
  ci_high =
    overall_row$ci95_high_mBq_m3,
  unit =
    "Per +1 plant at wind = 1 m/s",
  source_stage =
    "Earlier frozen overall analysis"
)

class_rows <- CLASS_MODELS[
  status ==
    "PASS" &
  model_id %in%
    c(
      "combustion_capable_upwind_count_20km_90deg",
      "noncombustion_only_upwind_count_20km_90deg",
      "joint_combustion",
      "joint_noncombustion"
    )
]

if (
  nrow(class_rows) !=
    4L
) {
  stopf(
    "燃烧/非燃烧大类模型不是4个PASS结果。"
  )
}

class_rows[
  ,
  category :=
    fifelse(
      model_id %in%
        c(
          "combustion_capable_upwind_count_20km_90deg",
          "joint_combustion"
        ),
      "Combustion-capable",
      "Non-combustion-only"
    )
]

class_rows[
  ,
  analysis :=
    fifelse(
      grepl(
        "^joint_",
        model_id
      ),
      "Joint broad-class model",
      "Separate broad-class model"
    )
]

FIG3_DATA <- rbindlist(
  list(
    FIG3_DATA,
    class_rows[
      ,
      .(
        category,
        analysis,
        estimate =
          estimate_per_1,
        ci_low =
          ci95_low_per_1,
        ci_high =
          ci95_high_per_1,
        unit =
          "Per +1 plant at wind = 1 m/s",
        source_stage =
          "Earlier frozen broad-class sensitivity analysis"
      )
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)

detail_rows <- OLD_EFFECTS[
  exposure_group %in%
    c(
      "nuclear",
      "geothermal",
      "wind",
      "solar",
      "hydro"
    ) &
  metric ==
    "count" &
  radius_km ==
    20 &
  direction ==
    "upwind" &
  effect_type ==
    "per_raw_unit" &
  wind_speed_ms ==
    1
]

if (
  nrow(detail_rows) >
    0L
) {

  detail_rows[
    ,
    category :=
      unname(
        FUEL_LABEL[
          exposure_group
        ]
      )
  ]

  FIG3_DATA <- rbindlist(
    list(
      FIG3_DATA,
      detail_rows[
        ,
        .(
          category,
          analysis =
            "Individual non-combustion / control type",
          estimate =
            estimate_mBq_m3,
          ci_low =
            ci95_low_mBq_m3,
          ci_high =
            ci95_high_mBq_m3,
          unit =
            "Per +1 plant at wind = 1 m/s",
          source_stage =
            "Earlier frozen overall analysis"
        )
      ]
    ),
    use.names = TRUE,
    fill = TRUE
  )
}

FIG3_DATA[
  ,
  display_label :=
    paste0(
      category,
      " — ",
      analysis
    )
]

FIG3_DATA[
  ,
  display_label :=
    factor(
      display_label,
      levels =
        rev(
          unique(
            display_label
          )
        )
    )
]

fwrite(
  FIG3_DATA,
  file.path(
    SOURCE_DIR,
    "Figure3_Table2_overall_broadclass_results.csv"
  )
)

p3 <- ggplot(
  FIG3_DATA,
  aes(
    x = estimate,
    y = display_label
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = 2
  ) +
  geom_errorbar(
    aes(
      xmin = ci_low,
      xmax = ci_high
    ),
    width = 0.18,
    orientation = "y"
  ) +
  geom_point(
    size = 2.2
  ) +
  labs(
    title = "Figure 3. Overall and broad power-plant class associations with gross beta",
    subtitle = "20-km upwind plant-count effects at wind speed = 1 m/s; first-stage frozen analysis",
    x = "Change in gross beta (mBq/m³) per additional plant",
    y = NULL
  ) +
  base_theme

save_plot_pair(
  p3,
  file.path(
    MAIN_FIG_DIR,
    "Figure3_overall_broadclass_associations"
  ),
  10,
  6.5
)

TABLE2 <- copy(
  FIG3_DATA
)

TABLE2[
  ,
  display_label :=
    NULL
]

fwrite(
  TABLE2,
  file.path(
    MAIN_TABLE_DIR,
    "Table2_overall_and_broadclass_main_results.csv"
  )
)

# ============================================================
# 12. Figure 4
# ============================================================

make_radius_panel <- function(
  branch_value,
  metric_value,
  title_text
) {

  D <- MAIN[
    branch ==
      branch_value &
    metric ==
      metric_value &
    wind_reference ==
      "样本中位数"
  ]

  D[
    ,
    radius_factor :=
      factor(
        radius_km,
        levels =
          RADII
      )
  ]

  ggplot(
    D,
    aes(
      x = radius_factor,
      y = estimate,
      group = 1
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    geom_errorbar(
      aes(
        ymin = ci95_low,
        ymax = ci95_high
      ),
      width = 0.14
    ) +
    geom_line() +
    geom_point(
      size = 1.7
    ) +
    facet_wrap(
      ~fuel_label,
      ncol = 3,
      scales = "free_y"
    ) +
    labs(
      title = title_text,
      x = "Cumulative radius (km)",
      y = "Adjusted gross-beta change (mBq/m³)"
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
}

p4a <- make_radius_panel(
  "filter_period",
  "count",
  "A. Active plant count"
)

p4b <- make_radius_panel(
  "filter_period",
  "heat",
  "B. Actual daily Heat Input"
)

FIG4 <- p4a / p4b +
  plot_annotation(
    title = "Figure 4. Fuel-specific upwind associations across predefined cumulative radii",
    subtitle = "Filter-period analysis; monitor-by-year fixed effects; effects evaluated at the sample median wind speed."
  )

save_plot_pair(
  FIG4,
  file.path(
    MAIN_FIG_DIR,
    "Figure4_fuel_specific_radius_associations"
  ),
  13,
  12
)

# ============================================================
# 13. Figure 5
# ============================================================

make_direction_panel <- function(
  metric_value,
  title_text
) {

  D <- DIRECTION[
    branch ==
      "filter_period" &
    metric ==
      metric_value &
    component ==
      "上风向减下风向" &
    wind_reference ==
      "样本中位数"
  ]

  D[
    ,
    radius_factor :=
      factor(
        radius_km,
        levels =
          RADII
      )
  ]

  ggplot(
    D,
    aes(
      x = radius_factor,
      y = estimate,
      group = 1
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    geom_errorbar(
      aes(
        ymin = ci95_low,
        ymax = ci95_high
      ),
      width = 0.14
    ) +
    geom_line() +
    geom_point(
      size = 1.7
    ) +
    facet_wrap(
      ~fuel_label,
      ncol = 3,
      scales = "free_y"
    ) +
    labs(
      title = title_text,
      x = "Cumulative radius (km)",
      y = "Upwind effect − downwind effect (mBq/m³)"
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
}

p5a <- make_direction_panel(
  "count",
  "A. Active plant count"
)

p5b <- make_direction_panel(
  "heat",
  "B. Actual daily Heat Input"
)

FIG5 <- p5a / p5b +
  plot_annotation(
    title = "Figure 5. Direct directional contrast between upwind and downwind exposure",
    subtitle = "Positive values indicate a larger upwind than downwind adjusted association."
  )

save_plot_pair(
  FIG5,
  file.path(
    MAIN_FIG_DIR,
    "Figure5_upwind_minus_downwind"
  ),
  13,
  12
)

# ============================================================
# 14. Table 3
# ============================================================

T3_MAIN <- MAIN[
  branch ==
    "filter_period" &
  radius_km ==
    20 &
  wind_reference ==
    "样本中位数",
  .(
    fuel,
    fuel_label,
    metric,
    metric_label,
    effect_unit =
      unit,
    wind_ms,
    estimate,
    ci95_low,
    ci95_high,
    p_wald,
    direction,
    ci_status
  )
]

T3_DIR <- DIRECTION[
  branch ==
    "filter_period" &
  radius_km ==
    20 &
  component ==
    "上风向减下风向" &
  wind_reference ==
    "样本中位数",
  .(
    fuel,
    metric,
    upwind_minus_downwind =
      estimate,
    direction_ci_low =
      ci95_low,
    direction_ci_high =
      ci95_high,
    direction_p =
      p_wald
  )
]

TABLE3 <- merge(
  T3_MAIN,
  T3_DIR,
  by = c(
    "fuel",
    "metric"
  ),
  all.x = TRUE,
  sort = FALSE
)

TABLE3[
  ,
  fuel_order__ :=
    match(
      fuel,
      FUELS
    )
]

if (
  TABLE3[
    is.na(fuel_order__),
    .N
  ] > 0L
) {
  stopf(
    "Table 3出现未定义燃料类型：",
    paste(
      unique(
        TABLE3[
          is.na(fuel_order__),
          fuel
        ]
      ),
      collapse = " | "
    )
  )
}

setorder(
  TABLE3,
  fuel_order__,
  metric
)

TABLE3[
  ,
  fuel_order__ :=
    NULL
]

fwrite(
  TABLE3,
  file.path(
    MAIN_TABLE_DIR,
    "Table3_20km_fuel_specific_main_results.csv"
  )
)

# ============================================================
# 15. Table 4
# ============================================================

TABLE4 <- copy(
  FUEL_SUMMARY
)

fwrite(
  TABLE4,
  file.path(
    MAIN_TABLE_DIR,
    "Table4_robustness_and_internal_validity_summary.csv"
  )
)

# ============================================================
# 16. Figure 6A
#     main vs strict time FE
# ============================================================

M20 <- MAIN[
  branch ==
    "filter_period" &
  radius_km ==
    20 &
  fuel %in%
    CORE_FUELS &
  wind_reference ==
    "样本中位数"
]

S20 <- STRICT[
  branch ==
    "filter_period" &
  radius_km ==
    20 &
  fuel %in%
    CORE_FUELS &
  wind_reference ==
    "样本中位数"
]

D6A <- rbindlist(
  list(
    M20[
      ,
      .(
        fuel,
        fuel_label,
        metric,
        metric_label,
        model =
          "Monitor × year FE",
        estimate,
        ci_low =
          ci95_low,
        ci_high =
          ci95_high
      )
    ],
    S20[
      ,
      .(
        fuel,
        fuel_label,
        metric,
        metric_label,
        model =
          "Monitor × year + common-time FE",
        estimate,
        ci_low =
          ci95_low,
        ci_high =
          ci95_high
      )
    ]
  )
)

D6A[
  ,
  fuel_label :=
    factor(
      fuel_label,
      levels =
        unname(
          FUEL_LABEL[
            CORE_FUELS
          ]
        )
    )
]

p6a <- ggplot(
  D6A,
  aes(
    x = fuel_label,
    y = estimate,
    shape = model
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  geom_errorbar(
    aes(
      ymin = ci_low,
      ymax = ci_high
    ),
    width = 0.15,
    position =
      position_dodge(
        width = 0.35
      )
  ) +
  geom_point(
    size = 2,
    position =
      position_dodge(
        width = 0.35
      )
  ) +
  facet_wrap(
    ~metric_label,
    scales = "free_y"
  ) +
  labs(
    title = "A. Main vs stricter common-time control",
    x = NULL,
    y = "Effect (mBq/m³)",
    shape = NULL
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "bottom"
  )

# ============================================================
# 17. Figure 6B
#     old random-intercept vs within/between
# ============================================================

D6B0 <- DECOMP[
  branch ==
    "filter_period" &
  radius_km ==
    20 &
  fuel %in%
    CORE_FUELS
]

D6B <- rbindlist(
  list(
    D6B0[
      ,
      .(
        fuel,
        fuel_label,
        metric,
        metric_label,
        model =
          "Original random-intercept",
        estimate =
          old_estimate_wind1,
        ci_low =
          old_ci_low_wind1,
        ci_high =
          old_ci_high_wind1
      )
    ],
    D6B0[
      ,
      .(
        fuel,
        fuel_label,
        metric,
        metric_label,
        model =
          "Within monitor-year",
        estimate =
          within_estimate_wind1,
        ci_low =
          within_ci_low_wind1,
        ci_high =
          within_ci_high_wind1
      )
    ],
    D6B0[
      ,
      .(
        fuel,
        fuel_label,
        metric,
        metric_label,
        model =
          "Between monitor-years",
        estimate =
          between_estimate_wind1,
        ci_low =
          between_ci_low_wind1,
        ci_high =
          between_ci_high_wind1
      )
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)

D6B[
  ,
  fuel_label :=
    factor(
      fuel_label,
      levels =
        unname(
          FUEL_LABEL[
            CORE_FUELS
          ]
        )
    )
]

p6b <- ggplot(
  D6B,
  aes(
    x = fuel_label,
    y = estimate,
    shape = model
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  geom_errorbar(
    aes(
      ymin = ci_low,
      ymax = ci_high
    ),
    width = 0.12,
    position =
      position_dodge(
        width = 0.45
      )
  ) +
  geom_point(
    size = 1.8,
    position =
      position_dodge(
        width = 0.45
      )
  ) +
  facet_wrap(
    ~metric_label,
    scales = "free_y"
  ) +
  labs(
    title = "B. Original vs within- and between-monitor-year associations",
    subtitle = "Effects evaluated at wind speed = 1 m/s",
    x = NULL,
    y = "Effect (mBq/m³)",
    shape = NULL
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "bottom"
  )

# ============================================================
# 18. Figure 6C
#     raw vs weighted
# ============================================================

W20 <- WEIGHTED[
  branch ==
    "filter_period" &
  radius_km ==
    20 &
  direction ==
    "upwind" &
  fuel %in%
    CORE_FUELS &
  wind_reference ==
    "样本中位数"
]

W20[
  ,
  raw_metric :=
    fifelse(
      metric ==
        "proximity",
      "count",
      fifelse(
        metric ==
          "heat_proximity",
        "heat",
        NA_character_
      )
    )
]

D6C <- rbindlist(
  list(
    M20[
      ,
      .(
        fuel,
        fuel_label,
        metric,
        metric_label,
        model =
          "Cumulative exposure",
        estimate,
        ci_low =
          ci95_low,
        ci_high =
          ci95_high
      )
    ],
    W20[
      !is.na(
        raw_metric
      ),
      .(
        fuel,
        fuel_label,
        metric =
          raw_metric,
        metric_label =
          unname(
            METRIC_LABEL[
              raw_metric
            ]
          ),
        model =
          "Distance weighted",
        estimate,
        ci_low =
          ci95_low,
        ci_high =
          ci95_high
      )
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)

p6c <- ggplot(
  D6C,
  aes(
    x = fuel_label,
    y = estimate,
    shape = model
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  geom_errorbar(
    aes(
      ymin = ci_low,
      ymax = ci_high
    ),
    width = 0.12,
    position =
      position_dodge(
        width = 0.35
      )
  ) +
  geom_point(
    size = 1.8,
    position =
      position_dodge(
        width = 0.35
      )
  ) +
  facet_wrap(
    ~metric_label,
    scales = "free_y"
  ) +
  labs(
    title = "C. Cumulative vs distance-weighted exposure",
    x = NULL,
    y = "Effect (mBq/m³)",
    shape = NULL
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "bottom"
  )

# ============================================================
# 19. Figure 6D
#     PM-adjusted vs no-PM
# ============================================================

P20 <- PM[
  branch ==
    "filter_period" &
  radius_km ==
    20 &
  fuel %in%
    CORE_FUELS &
  wind_reference ==
    "样本中位数"
]

D6D <- rbindlist(
  list(
    M20[
      ,
      .(
        fuel,
        fuel_label,
        metric,
        metric_label,
        model =
          "PM2.5 adjusted",
        estimate,
        ci_low =
          ci95_low,
        ci_high =
          ci95_high
      )
    ],
    P20[
      ,
      .(
        fuel,
        fuel_label,
        metric,
        metric_label,
        model =
          "No PM2.5 adjustment",
        estimate,
        ci_low =
          ci95_low,
        ci_high =
          ci95_high
      )
    ]
  )
)

p6d <- ggplot(
  D6D,
  aes(
    x = fuel_label,
    y = estimate,
    shape = model
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2
  ) +
  geom_errorbar(
    aes(
      ymin = ci_low,
      ymax = ci_high
    ),
    width = 0.12,
    position =
      position_dodge(
        width = 0.35
      )
  ) +
  geom_point(
    size = 1.8,
    position =
      position_dodge(
        width = 0.35
      )
  ) +
  facet_wrap(
    ~metric_label,
    scales = "free_y"
  ) +
  labs(
    title = "D. PM2.5 adjustment sensitivity",
    x = NULL,
    y = "Effect (mBq/m³)",
    shape = NULL
  ) +
  base_theme +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "bottom"
  )

FIG6 <- (
  p6a |
    p6b
) /
  (
    p6c |
      p6d
  ) +
  plot_annotation(
    title = "Figure 6. Internal-validity and key robustness analyses",
    subtitle = "Panels use the predefined 20-km window; complete radius-specific results are provided in the Supplement."
  )

save_plot_pair(
  FIG6,
  file.path(
    MAIN_FIG_DIR,
    "Figure6_internal_validity_and_robustness"
  ),
  15,
  12
)

# ============================================================
# 20. Supplementary Figure S1
#     detailed fuel maps in facets
# ============================================================

PLANT_S1 <- PLANT_MAP[
  fuel_group_plot !=
    "mixed_multiple"
]

pS1 <- ggplot() +
  geom_polygon(
    data = US_MAP,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "white",
    color = "grey82",
    linewidth = 0.2
  ) +
  geom_point(
    data = PLANT_S1,
    aes(
      x = plant_longitude,
      y = plant_latitude
    ),
    size = 0.45,
    alpha = 0.55
  ) +
  facet_wrap(
    ~fuel_label,
    ncol = 3
  ) +
  coord_quickmap(
    xlim = c(
      -125,
      -66
    ),
    ylim = c(
      24,
      50
    ),
    expand = FALSE
  ) +
  labs(
    title = "Figure S1. Spatial distribution of individual power-plant energy / fuel types"
  ) +
  theme_void(
    base_size = 9
  ) +
  theme(
    strip.text = element_text(
      face = "bold"
    ),
    plot.title = element_text(
      face = "bold"
    )
  )

save_plot_pair(
  pS1,
  file.path(
    SUPP_FIG_DIR,
    "FigureS1_individual_powerplant_type_maps"
  ),
  13,
  13
)

# ============================================================
# 21. Supplementary Figure S2/S3
# ============================================================

pS2 <- make_radius_panel(
  "quasi_daily",
  "count",
  "Quasi-daily active plant count"
)

pS2 <- pS2 +
  labs(
    title = "Figure S2. Quasi-daily fuel-specific plant-count associations across all predefined radii"
  )

save_plot_pair(
  pS2,
  file.path(
    SUPP_FIG_DIR,
    "FigureS2_quasidaily_count_all_radii"
  ),
  13,
  7
)

pS3 <- make_radius_panel(
  "quasi_daily",
  "heat",
  "Quasi-daily actual Heat Input"
)

pS3 <- pS3 +
  labs(
    title = "Figure S3. Quasi-daily fuel-specific Heat-Input associations across all predefined radii"
  )

save_plot_pair(
  pS3,
  file.path(
    SUPP_FIG_DIR,
    "FigureS3_quasidaily_heat_all_radii"
  ),
  13,
  7
)

# ============================================================
# 22. Supplementary Figure S4
#     continuous wind curves
# ============================================================

make_wind_curve_panel <- function(
  branch_value,
  metric_value,
  panel_title
) {

  D <- CURVE[
    branch ==
      branch_value &
    metric ==
      metric_value &
    radius_km %in%
      REP_RADII
  ]

  D[
    ,
    radius_label :=
      factor(
        paste0(
          radius_km,
          " km"
        ),
        levels =
          paste0(
            REP_RADII,
            " km"
          )
      )
  ]

  ggplot(
    D,
    aes(
      x = wind_ms,
      y = estimate
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    geom_ribbon(
      aes(
        ymin = ci95_low,
        ymax = ci95_high
      ),
      alpha = 0.15
    ) +
    geom_line() +
    facet_grid(
      fuel_label ~ radius_label,
      scales = "free_y"
    ) +
    labs(
      title = panel_title,
      x = "Wind speed (m/s)",
      y = "Exposure effect (mBq/m³)"
    ) +
    base_theme +
    theme(
      strip.text = element_text(
        size = 7.5
      )
    )
}

s4a <- make_wind_curve_panel(
  "filter_period",
  "count",
  "A. Filter-period / Count"
)

s4b <- make_wind_curve_panel(
  "filter_period",
  "heat",
  "B. Filter-period / Heat Input"
)

s4c <- make_wind_curve_panel(
  "quasi_daily",
  "count",
  "C. Quasi-daily / Count"
)

s4d <- make_wind_curve_panel(
  "quasi_daily",
  "heat",
  "D. Quasi-daily / Heat Input"
)

FIGS4 <- (
  s4a /
    s4b /
    s4c /
    s4d
) +
  plot_annotation(
    title = "Figure S4. Continuous wind-speed effect modification",
    subtitle = "20, 50, 100, and 200 km are predefined representative radii for visualization only."
  )

save_plot_pair(
  FIGS4,
  file.path(
    SUPP_FIG_DIR,
    "FigureS4_continuous_wind_effect_curves"
  ),
  14,
  36
)

# ============================================================
# 23. General comparison plot helper
# ============================================================

make_two_model_radius_plot <- function(
  D1,
  D2,
  label1,
  label2,
  branch_value,
  metric_value,
  title_text
) {

  A <- D1[
    branch ==
      branch_value &
    metric ==
      metric_value &
    wind_reference ==
      "样本中位数",
    .(
      fuel,
      fuel_label,
      radius_km,
      model =
        label1,
      estimate,
      ci_low =
        ci95_low,
      ci_high =
        ci95_high
    )
  ]

  B <- D2[
    branch ==
      branch_value &
    metric ==
      metric_value &
    wind_reference ==
      "样本中位数",
    .(
      fuel,
      fuel_label,
      radius_km,
      model =
        label2,
      estimate,
      ci_low =
        ci95_low,
      ci_high =
        ci95_high
    )
  ]

  D <- rbindlist(
    list(
      A,
      B
    ),
    use.names = TRUE
  )

  D[
    ,
    radius_factor :=
      factor(
        radius_km,
        levels =
          RADII
      )
  ]

  ggplot(
    D,
    aes(
      x = radius_factor,
      y = estimate,
      group = model,
      linetype = model,
      shape = model
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    geom_line() +
    geom_point(
      size = 1.4
    ) +
    facet_wrap(
      ~fuel_label,
      ncol = 3,
      scales = "free_y"
    ) +
    labs(
      title = title_text,
      x = "Cumulative radius (km)",
      y = "Effect (mBq/m³)",
      shape = NULL,
      linetype = NULL
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "bottom"
    )
}

# ============================================================
# 24. Figure S5
# ============================================================

s5a <- make_two_model_radius_plot(
  MAIN,
  STRICT,
  "Monitor × year FE",
  "Monitor × year + time FE",
  "filter_period",
  "count",
  "A. Filter-period / Count"
)

s5b <- make_two_model_radius_plot(
  MAIN,
  STRICT,
  "Monitor × year FE",
  "Monitor × year + time FE",
  "filter_period",
  "heat",
  "B. Filter-period / Heat Input"
)

s5c <- make_two_model_radius_plot(
  MAIN,
  STRICT,
  "Monitor × year FE",
  "Monitor × year + time FE",
  "quasi_daily",
  "count",
  "C. Quasi-daily / Count"
)

s5d <- make_two_model_radius_plot(
  MAIN,
  STRICT,
  "Monitor × year FE",
  "Monitor × year + time FE",
  "quasi_daily",
  "heat",
  "D. Quasi-daily / Heat Input"
)

FIGS5 <- (
  s5a /
    s5b /
    s5c /
    s5d
) +
  plot_annotation(
    title = "Figure S5. Main monitor-by-year fixed-effects models versus stricter common-time control"
  )

save_plot_pair(
  FIGS5,
  file.path(
    SUPP_FIG_DIR,
    "FigureS5_strict_time_fixed_effects"
  ),
  13,
  28
)

# ============================================================
# 25. Figure S6
#     raw vs weighted
# ============================================================

make_weighted_panel <- function(
  branch_value,
  metric_value,
  title_text
) {

  rawD <- MAIN[
    branch ==
      branch_value &
    metric ==
      metric_value &
    wind_reference ==
      "样本中位数",
    .(
      fuel,
      fuel_label,
      radius_km,
      model =
        "Cumulative",
      estimate
    )
  ]

  weighted_metric <- if (
    metric_value ==
      "count"
  ) {
    "proximity"
  } else {
    "heat_proximity"
  }

  WD <- WEIGHTED[
    branch ==
      branch_value &
    metric ==
      weighted_metric &
    direction ==
      "upwind" &
    wind_reference ==
      "样本中位数",
    .(
      fuel,
      fuel_label,
      radius_km,
      model =
        "Distance weighted",
      estimate
    )
  ]

  D <- rbindlist(
    list(
      rawD,
      WD
    )
  )

  D[
    ,
    radius_factor :=
      factor(
        radius_km,
        levels =
          RADII
      )
  ]

  ggplot(
    D,
    aes(
      x = radius_factor,
      y = estimate,
      group = model,
      linetype = model,
      shape = model
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    geom_line() +
    geom_point(
      size = 1.4
    ) +
    facet_wrap(
      ~fuel_label,
      ncol = 3,
      scales = "free_y"
    ) +
    labs(
      title = title_text,
      x = "Cumulative radius (km)",
      y = "Effect estimate",
      linetype = NULL,
      shape = NULL
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "bottom"
    )
}

ws1 <- make_weighted_panel(
  "filter_period",
  "count",
  "A. Filter-period / Count vs Σ1/d"
)

ws2 <- make_weighted_panel(
  "filter_period",
  "heat",
  "B. Filter-period / Heat vs ΣHeat/d"
)

ws3 <- make_weighted_panel(
  "quasi_daily",
  "count",
  "C. Quasi-daily / Count vs Σ1/d"
)

ws4 <- make_weighted_panel(
  "quasi_daily",
  "heat",
  "D. Quasi-daily / Heat vs ΣHeat/d"
)

FIGS6 <- (
  ws1 /
    ws2 /
    ws3 /
    ws4
) +
  plot_annotation(
    title = "Figure S6. Cumulative versus distance-weighted exposure sensitivity"
  )

save_plot_pair(
  FIGS6,
  file.path(
    SUPP_FIG_DIR,
    "FigureS6_distance_weighted_sensitivity"
  ),
  13,
  28
)

# ============================================================
# 26. Figure S7
#     PM sensitivity
# ============================================================

ps1 <- make_two_model_radius_plot(
  MAIN,
  PM,
  "PM2.5 adjusted",
  "No PM2.5 adjustment",
  "filter_period",
  "count",
  "A. Filter-period / Count"
)

ps2 <- make_two_model_radius_plot(
  MAIN,
  PM,
  "PM2.5 adjusted",
  "No PM2.5 adjustment",
  "filter_period",
  "heat",
  "B. Filter-period / Heat Input"
)

ps3 <- make_two_model_radius_plot(
  MAIN,
  PM,
  "PM2.5 adjusted",
  "No PM2.5 adjustment",
  "quasi_daily",
  "count",
  "C. Quasi-daily / Count"
)

ps4 <- make_two_model_radius_plot(
  MAIN,
  PM,
  "PM2.5 adjusted",
  "No PM2.5 adjustment",
  "quasi_daily",
  "heat",
  "D. Quasi-daily / Heat Input"
)

FIGS7 <- (
  ps1 /
    ps2 /
    ps3 /
    ps4
) +
  plot_annotation(
    title = "Figure S7. PM2.5-adjustment sensitivity using identical analysis rows"
  )

save_plot_pair(
  FIGS7,
  file.path(
    SUPP_FIG_DIR,
    "FigureS7_PM25_sensitivity"
  ),
  13,
  28
)

# ============================================================
# 27. Figure S8
#     lag analysis
# ============================================================

LAG_ORDER <- c(
  "lag0",
  "lag1",
  "lag2",
  "lag3",
  "mean_0_1",
  "mean_0_2"
)

LAG_LABELS <- c(
  lag0 = "Lag 0",
  lag1 = "Lag 1",
  lag2 = "Lag 2",
  lag3 = "Lag 3",
  mean_0_1 = "Mean 0–1",
  mean_0_2 = "Mean 0–2"
)

make_lag_panel <- function(
  metric_value,
  title_text
) {

  D <- LAG[
    metric ==
      metric_value &
    radius_km %in%
      REP_RADII &
    wind_reference ==
      "样本中位数"
  ]

  D[
    ,
    lag_label :=
      factor(
        unname(
          LAG_LABELS[
            lag_spec
          ]
        ),
        levels =
          unname(
            LAG_LABELS[
              LAG_ORDER
            ]
          )
      )
  ]

  D[
    ,
    radius_label :=
      factor(
        paste0(
          radius_km,
          " km"
        ),
        levels =
          paste0(
            REP_RADII,
            " km"
          )
      )
  ]

  ggplot(
    D,
    aes(
      x = lag_label,
      y = estimate,
      group = radius_label,
      linetype = radius_label,
      shape = radius_label
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    geom_line() +
    geom_point(
      size = 1.5
    ) +
    facet_wrap(
      ~fuel_label,
      ncol = 3,
      scales = "free_y"
    ) +
    labs(
      title = title_text,
      x = "Exposure lag definition",
      y = "Effect (mBq/m³)",
      linetype = "Radius",
      shape = "Radius"
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),
      legend.position = "bottom"
    )
}

ls1 <- make_lag_panel(
  "count",
  "A. Plant count"
)

ls2 <- make_lag_panel(
  "heat",
  "B. Heat Input"
)

FIGS8 <- ls1 / ls2 +
  plot_annotation(
    title = "Figure S8. Prespecified quasi-daily lag analyses",
    subtitle = "No lag was selected on the basis of statistical significance."
  )

save_plot_pair(
  FIGS8,
  file.path(
    SUPP_FIG_DIR,
    "FigureS8_quasidaily_lag_analysis"
  ),
  13,
  14
)

# ============================================================
# 28. Supplementary Tables S1–S10
# ============================================================

for (ii in 1:10) {

  pattern <- sprintf(
    "^%02d_.*\\.csv$",
    ii
  )

  src <- list.files(
    H11_TABLE_DIR,
    pattern = pattern,
    full.names = TRUE
  )

  if (
    length(src) !=
      1L
  ) {
    stopf(
      "无法唯一定位H11论文表 ",
      ii,
      "；候选数=",
      length(src)
    )
  }

  dst <- file.path(
    SUPP_TABLE_DIR,
    paste0(
      "TableS",
      ii,
      "_",
      sub(
        paste0(
          "^",
          sprintf(
            "%02d_",
            ii
          )
        ),
        "",
        basename(src)
      )
    )
  )

  ok <- file.copy(
    src,
    dst,
    overwrite = TRUE
  )

  if (!ok) {
    stopf(
      "复制Table S",
      ii,
      "失败。"
    )
  }
}

# ============================================================
# 29. Additional source tables for Figure 6 / supplement
# ============================================================

fwrite(
  D6A,
  file.path(
    SOURCE_DIR,
    "Figure6A_main_vs_strict_20km.csv"
  )
)

fwrite(
  D6B,
  file.path(
    SOURCE_DIR,
    "Figure6B_original_within_between_20km.csv"
  )
)

fwrite(
  D6C,
  file.path(
    SOURCE_DIR,
    "Figure6C_raw_vs_weighted_20km.csv"
  )
)

fwrite(
  D6D,
  file.path(
    SOURCE_DIR,
    "Figure6D_PM_sensitivity_20km.csv"
  )
)

# ============================================================
# 30. Figure/Table captions and roles
# ============================================================

CAPTIONS <- c(
  "# 07H12 正文与补充材料图表说明",
  "",
  "## 正文图",
  "",
  "### Figure 1",
  "Spatial distribution of U.S. power plants and RadNet monitoring stations. RadNet stations are shown as triangles and power plants as circles. Panel A distinguishes individual energy/fuel types; Panel B groups power plants into broad analytical classes.",
  "",
  "### Figure 2",
  "Study design and analytical framework integrating RadNet gross-beta monitoring, power-plant locations and fuel types, actual daily CAMPD operation, meteorology, directional exposure assessment, broad-class analysis, and fuel-specific internal-validity analyses.",
  "",
  "### Figure 3",
  "Overall and broad power-plant class associations with gross beta. These estimates originate from the earlier frozen first-stage overall/broad-class analysis and are explicitly kept separate from the updated 07H fuel-specific analytical sample.",
  "",
  "### Figure 4",
  "Fuel-specific upwind associations across all ten predefined cumulative radii using filter-period exposure. Effects are evaluated at each model sample's median wind speed. Radii are not selected according to statistical significance.",
  "",
  "### Figure 5",
  "Direct model-based contrast between upwind and downwind exposure effects. Positive estimates indicate larger upwind than downwind associations.",
  "",
  "### Figure 6",
  "Internal-validity and key robustness analyses at the predefined 20-km anchor: stricter common-time control; original versus within/between monitor-year associations; cumulative versus distance-weighted exposure; and PM2.5-adjustment sensitivity.",
  "",
  "## 正文表",
  "",
  "### Table 1",
  "Study population and power-plant characteristics. Earlier overall-analysis and updated 07H samples are reported separately rather than being treated as identical.",
  "",
  "### Table 2",
  "Overall and broad-class 20-km upwind plant-count results from the frozen first-stage analysis.",
  "",
  "### Table 3",
  "Fuel-specific 20-km filter-period results from the monitor-by-year fixed-effects analysis, together with direct upwind-minus-downwind contrasts.",
  "",
  "### Table 4",
  "Fuel-level robustness and internal-validity summary across all predefined radii.",
  "",
  "## Supplementary Figures",
  "",
  "- Figure S1: Individual power-plant type maps.",
  "- Figure S2: Quasi-daily Count results across all radii.",
  "- Figure S3: Quasi-daily Heat Input results across all radii.",
  "- Figure S4: Continuous wind-speed effect curves at 20, 50, 100 and 200 km.",
  "- Figure S5: Main versus stricter common-time fixed-effects results.",
  "- Figure S6: Cumulative versus distance-weighted exposure results.",
  "- Figure S7: PM2.5-adjusted versus unadjusted identical-row sensitivity.",
  "- Figure S8: Prespecified quasi-daily lag results.",
  "",
  "## Supplementary Tables",
  "",
  "- Table S1–S10 are direct renamed copies of H11 formal Tables 1–10. No model is re-fitted in H12.",
  "",
  "## Interpretation rules",
  "",
  "- Positive/negative estimates are adjusted association directions, not proof of causality.",
  "- All ten cumulative radii are predefined and retained; no radius is selected according to statistical significance.",
  "- 20, 50, 100 and 200 km in wind/lag visualizations are predefined display radii only.",
  "- Σ1/d and ΣHeat/d are distance-weighted exposure sensitivities, not physical atmospheric concentration or radiation dose.",
  "- The first-stage overall/broad-class analysis and later corrected fuel-specific 07H analysis are displayed as separate analytical stages when their formal samples differ."
)

writeLines(
  CAPTIONS,
  file.path(
    OUT,
    "00_论文图表说明.md"
  ),
  useBytes = TRUE
)

# ============================================================
# 31. Manifest
# ============================================================

main_fig_files <- list.files(
  MAIN_FIG_DIR,
  pattern =
    "\\.(png|pdf)$",
  full.names =
    TRUE
)

supp_fig_files <- list.files(
  SUPP_FIG_DIR,
  pattern =
    "\\.(png|pdf)$",
  full.names =
    TRUE
)

main_table_files <- list.files(
  MAIN_TABLE_DIR,
  pattern =
    "\\.csv$",
  full.names =
    TRUE
)

supp_table_files <- list.files(
  SUPP_TABLE_DIR,
  pattern =
    "\\.csv$",
  full.names =
    TRUE
)

MANIFEST <- rbindlist(
  list(
    data.table(
      category =
        "Main figure",
      file =
        main_fig_files
    ),
    data.table(
      category =
        "Main table",
      file =
        main_table_files
    ),
    data.table(
      category =
        "Supplementary figure",
      file =
        supp_fig_files
    ),
    data.table(
      category =
        "Supplementary table",
      file =
        supp_table_files
    )
  ),
  use.names = TRUE
)

MANIFEST[
  ,
  filename :=
    basename(
      file
    )
]

fwrite(
  MANIFEST,
  file.path(
    AUDIT_DIR,
    "01_output_manifest.csv"
  )
)

# ============================================================
# 32. Hard checks
# ============================================================

N_MAIN_PNG <- length(
  list.files(
    MAIN_FIG_DIR,
    pattern =
      "\\.png$"
  )
)

N_MAIN_PDF <- length(
  list.files(
    MAIN_FIG_DIR,
    pattern =
      "\\.pdf$"
  )
)

N_SUPP_PNG <- length(
  list.files(
    SUPP_FIG_DIR,
    pattern =
      "\\.png$"
  )
)

N_SUPP_PDF <- length(
  list.files(
    SUPP_FIG_DIR,
    pattern =
      "\\.pdf$"
  )
)

N_MAIN_TABLE <- length(
  list.files(
    MAIN_TABLE_DIR,
    pattern =
      "\\.csv$"
  )
)

N_SUPP_TABLE <- length(
  list.files(
    SUPP_TABLE_DIR,
    pattern =
      "\\.csv$"
  )
)

HARD <- data.table(
  check = c(
    "H11_final_freeze_exists",
    "H10_hard_checks_all_PASS",
    "Figure1_main_created",
    "Figure2_main_created",
    "Figure3_main_created",
    "Figure4_main_created",
    "Figure5_main_created",
    "Figure6_main_created",
    "main_PNG_6",
    "main_PDF_6",
    "main_tables_4",
    "supplementary_PNG_8",
    "supplementary_PDF_8",
    "supplementary_tables_10",
    "Figure3_source_nonempty",
    "Table3_20km_nonempty",
    "map_powerplant_nonempty",
    "map_monitor_nonempty"
  ),
  observed = c(
    as.character(
      as.integer(
        file.exists(
          H11_FREEZE
        )
      )
    ),
    as.character(
      H10_HARD[
        status !=
          "PASS",
        .N
      ]
    ),
    as.character(
      as.integer(
        file.exists(
          file.path(
            MAIN_FIG_DIR,
            "Figure1_US_powerplants_RadNet_map.png"
          )
        )
      )
    ),
    as.character(
      as.integer(
        file.exists(
          file.path(
            MAIN_FIG_DIR,
            "Figure2_study_design_flowchart.png"
          )
        )
      )
    ),
    as.character(
      as.integer(
        file.exists(
          file.path(
            MAIN_FIG_DIR,
            "Figure3_overall_broadclass_associations.png"
          )
        )
      )
    ),
    as.character(
      as.integer(
        file.exists(
          file.path(
            MAIN_FIG_DIR,
            "Figure4_fuel_specific_radius_associations.png"
          )
        )
      )
    ),
    as.character(
      as.integer(
        file.exists(
          file.path(
            MAIN_FIG_DIR,
            "Figure5_upwind_minus_downwind.png"
          )
        )
      )
    ),
    as.character(
      as.integer(
        file.exists(
          file.path(
            MAIN_FIG_DIR,
            "Figure6_internal_validity_and_robustness.png"
          )
        )
      )
    ),
    as.character(
      N_MAIN_PNG
    ),
    as.character(
      N_MAIN_PDF
    ),
    as.character(
      N_MAIN_TABLE
    ),
    as.character(
      N_SUPP_PNG
    ),
    as.character(
      N_SUPP_PDF
    ),
    as.character(
      N_SUPP_TABLE
    ),
    as.character(
      as.integer(
        nrow(
          FIG3_DATA
        ) >
          0
      )
    ),
    as.character(
      as.integer(
        nrow(
          TABLE3
        ) >
          0
      )
    ),
    as.character(
      as.integer(
        nrow(
          PLANT_MAP
        ) >
          0
      )
    ),
    as.character(
      as.integer(
        nrow(
          ST
        ) >
          0
      )
    )
  ),
  expected = c(
    "1",
    "0",
    "1",
    "1",
    "1",
    "1",
    "1",
    "1",
    "6",
    "6",
    "4",
    "8",
    "8",
    "10",
    "1",
    "1",
    "1",
    "1"
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
    AUDIT_DIR,
    "02_07H12_hard_checks.csv"
  )
)

cat(
  "\n============================================================\n",
  "07H12 OUTPUT SUMMARY\n",
  "============================================================\n",
  sep = ""
)

cat(
  "Main figures: ",
  N_MAIN_PNG,
  " PNG + ",
  N_MAIN_PDF,
  " PDF\n",
  sep = ""
)

cat(
  "Main tables: ",
  N_MAIN_TABLE,
  "\n",
  sep = ""
)

cat(
  "Supplementary figures: ",
  N_SUPP_PNG,
  " PNG + ",
  N_SUPP_PDF,
  " PDF\n",
  sep = ""
)

cat(
  "Supplementary tables: ",
  N_SUPP_TABLE,
  "\n",
  sep = ""
)

cat(
  "\n============================================================\n",
  "07H12 HARD CHECKS\n",
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

REPORT <- c(
  "# 07H12 论文正式图表总整理",
  "",
  paste0(
    "- 生成时间：",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  "",
  "## 正文",
  "",
  "- Figure 1–6",
  "- Table 1–4",
  "",
  "## Supplementary",
  "",
  "- Figure S1–S8",
  "- Table S1–S10",
  "",
  "## 重要原则",
  "",
  "- 07H12不重新拟合任何模型。",
  "- 07H01–07H11保持只读。",
  "- 总体/燃烧-非燃烧分析与后期07H燃料特异分析按各自冻结样本明确区分。",
  "- 全部10个预设累计半径保留，不按照显著性挑选。",
  "- 20/50/100/200 km仅用于部分复杂曲线的预设展示。",
  "- 距离加权Σ1/d和ΣHeat/d不解释为实际浓度或剂量。",
  "",
  paste0(
    "- Hard FAIL = ",
    N_FAIL
  ),
  "",
  if (
    N_FAIL ==
      0L
  ) {
    "**07H12 PAPER FIGURE/TABLE PREPARATION: PASS**"
  } else {
    "**07H12 PAPER FIGURE/TABLE PREPARATION: FAIL**"
  }
)

writeLines(
  REPORT,
  file.path(
    OUT,
    "07H12_论文正式图表总整理报告.md"
  ),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    AUDIT_DIR,
    "03_sessionInfo.txt"
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
    "07H12 PAPER FIGURE/TABLE PREPARATION: PASS\n"
  )

} else {

  cat(
    "07H12 PAPER FIGURE/TABLE PREPARATION: FAIL\n"
  )
}

cat(
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
    "07H12存在Hard FAIL，请勿用于论文正式交付。"
  )
}

