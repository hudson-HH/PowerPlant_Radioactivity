#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999
)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(maps)
  library(patchwork)
  library(grid)
})

cat(
  "\n============================================================\n",
  "07H12 Figure 1：非互斥真实能源类型地图重绘\n",
  "============================================================\n",
  sep = ""
)

# ============================================================
# 1. 路径
# ============================================================

H12_DIR <- "/Users/hudson/Desktop/USA 辐射表征研究/发电站数据处理/06_07H_燃料特异CAMPD逐日暴露与模型/07H12_论文正式图表总整理"

SOURCE_DIR <- file.path(
  H12_DIR,
  "05_绘图数据"
)

FIG_DIR <- file.path(
  H12_DIR,
  "01_正文图"
)

AUDIT_DIR <- file.path(
  H12_DIR,
  "06_审计"
)

PLANT_FILE <- file.path(
  SOURCE_DIR,
  "Figure1_powerplant_locations_by_type.csv"
)

MONITOR_FILE <- file.path(
  SOURCE_DIR,
  "Figure1_radnet_monitor_locations.csv"
)

OUT_PNG <- file.path(
  FIG_DIR,
  "Figure1_US_powerplants_RadNet_map.png"
)

OUT_PDF <- file.path(
  FIG_DIR,
  "Figure1_US_powerplants_RadNet_map.pdf"
)

AUDIT_FILE <- file.path(
  AUDIT_DIR,
  "Figure1_energy_type_map_audit.csv"
)

EXPANDED_FILE <- file.path(
  SOURCE_DIR,
  "Figure1A_powerplant_locations_nonexclusive_energy_types.csv"
)

for (
  f in c(
    PLANT_FILE,
    MONITOR_FILE
  )
) {

  if (!file.exists(f)) {
    stop(
      paste0(
        "缺少Figure 1输入：\n",
        f
      ),
      call. = FALSE
    )
  }
}

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  AUDIT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 2. 读取07H12已审计的地图输入
# ============================================================

P0 <- fread(
  PLANT_FILE,
  showProgress = FALSE
)

M <- fread(
  MONITOR_FILE,
  showProgress = FALSE
)

required_P <- c(
  "plant_id_plot",
  "plant_longitude",
  "plant_latitude",
  "all_fuel_groups",
  "broad_group"
)

required_M <- c(
  "monitor_longitude",
  "monitor_latitude"
)

miss_P <- setdiff(
  required_P,
  names(P0)
)

miss_M <- setdiff(
  required_M,
  names(M)
)

if (length(miss_P) > 0L) {
  stop(
    paste0(
      "发电站地图源缺少字段：",
      paste(
        miss_P,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

if (length(miss_M) > 0L) {
  stop(
    paste0(
      "RadNet地图源缺少字段：",
      paste(
        miss_M,
        collapse = ", "
      )
    ),
    call. = FALSE
  )
}

# ============================================================
# 3. 基础清理
# ============================================================

clean_text <- function(x) {

  x <- as.character(x)

  x <- gsub(
    "\u00A0",
    " ",
    x,
    fixed = TRUE
  )

  x <- trimws(x)

  x <- gsub(
    "[[:space:]]+",
    " ",
    x
  )

  x
}

P0[
  ,
  `:=`(
    plant_id_plot =
      clean_text(
        plant_id_plot
      ),

    plant_longitude =
      suppressWarnings(
        as.numeric(
          plant_longitude
        )
      ),

    plant_latitude =
      suppressWarnings(
        as.numeric(
          plant_latitude
        )
      ),

    all_fuel_groups =
      clean_text(
        all_fuel_groups
      ),

    broad_group =
      clean_text(
        broad_group
      )
  )
]

M[
  ,
  `:=`(
    monitor_longitude =
      suppressWarnings(
        as.numeric(
          monitor_longitude
        )
      ),

    monitor_latitude =
      suppressWarnings(
        as.numeric(
          monitor_latitude
        )
      )
  )
]

P0 <- P0[
  !is.na(
    plant_id_plot
  ) &
  plant_id_plot != "" &
  is.finite(
    plant_longitude
  ) &
  is.finite(
    plant_latitude
  )
]

M <- M[
  is.finite(
    monitor_longitude
  ) &
  is.finite(
    monitor_latitude
  )
]

if (
  P0[
    ,
    uniqueN(
      plant_id_plot
    )
  ] !=
    nrow(P0)
) {

  stop(
    "Figure1_powerplant_locations_by_type.csv中的plant_id_plot不是每座物理电厂一行。",
    call. = FALSE
  )
}

cat(
  "物理发电站：",
  nrow(P0),
  "\n",
  "RadNet监测站：",
  nrow(M),
  "\n",
  sep = ""
)

# ============================================================
# 4. A面板：
#    从all_fuel_groups恢复真实、非互斥能源类型
#
# 例：
# all_fuel_groups = "coal|gas"
#
# 则同一物理电厂分别进入：
# Coal
# Natural gas
#
# 不再生成Mixed / multiple fuels。
# ============================================================

P_A <- P0[
  ,
  {

    fuels <- unlist(
      strsplit(
        all_fuel_groups,
        "\\|",
        perl = TRUE
      )
    )

    fuels <- clean_text(
      fuels
    )

    fuels <- tolower(
      fuels
    )

    fuels <- unique(
      fuels[
        !is.na(fuels) &
        fuels != ""
      ]
    )

    if (length(fuels) == 0L) {
      fuels <- "other"
    }

    data.table(
      fuel_raw =
        fuels
    )
  },
  by = .(
    plant_id_plot,
    plant_longitude,
    plant_latitude
  )
]

# ============================================================
# 5. 原始fuel group -> 正式英文显示名
# ============================================================

fuel_dictionary <- c(
  "coal" =
    "Coal",

  "gas" =
    "Natural gas",

  "natural_gas" =
    "Natural gas",

  "natural gas" =
    "Natural gas",

  "oil" =
    "Oil",

  "other_fossil" =
    "Other fossil",

  "other fossil" =
    "Other fossil",

  "biomass" =
    "Biomass",

  "solid_waste" =
    "Solid waste",

  "solid waste" =
    "Solid waste",

  "nuclear" =
    "Nuclear",

  "geothermal" =
    "Geothermal",

  "wind" =
    "Wind",

  "solar" =
    "Solar",

  "hydro" =
    "Hydro",

  "hydropower" =
    "Hydro",

  "other" =
    "Other"
)

P_A[
  ,
  fuel_label :=
    unname(
      fuel_dictionary[
        fuel_raw
      ]
    )
]

unknown_fuels <- P_A[
  is.na(
    fuel_label
  ),
  sort(
    unique(
      fuel_raw
    )
  )
]

if (length(unknown_fuels) > 0L) {

  cat(
    "\n出现未定义的原始能源类型：\n"
  )

  print(
    unknown_fuels
  )

  stop(
    "Figure 1A存在未识别能源类型，已停止以避免错误归类。",
    call. = FALSE
  )
}

P_A <- unique(
  P_A[
    ,
    .(
      plant_id_plot,
      plant_longitude,
      plant_latitude,
      fuel_label
    )
  ]
)

# ============================================================
# 6. A面板计数审计
#
# 注意：
# 这是非互斥计数。
# 同一物理电厂可进入多个能源类型。
# ============================================================

A_COUNTS <- P_A[
  ,
  .(
    n_physical_plants =
      uniqueN(
        plant_id_plot
      )
  ),
  by =
    fuel_label
]

DETAIL_LEVELS_MASTER <- c(
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
  "Other"
)

A_COUNTS[
  ,
  order__ :=
    match(
      fuel_label,
      DETAIL_LEVELS_MASTER
    )
]

setorder(
  A_COUNTS,
  order__
)

A_COUNTS[
  ,
  order__ :=
    NULL
]

cat(
  "\n============================================================\n",
  "A图：非互斥真实能源类型物理电厂数\n",
  "============================================================\n",
  sep = ""
)

print(
  A_COUNTS,
  nrows = Inf
)

# ============================================================
# 7. B面板继续保留07H12原来的互斥宽分类
# ============================================================

P_B <- copy(
  P0
)

broad_dictionary <- c(
  "combustion-capable" =
    "Combustion-capable",

  "non-combustion major energy" =
    "Non-combustion major energy",

  "wind / solar / hydro control" =
    "Wind / solar / hydro control",

  "other" =
    "Other"
)

P_B[
  ,
  broad_key :=
    tolower(
      broad_group
    )
]

P_B[
  ,
  broad_label :=
    unname(
      broad_dictionary[
        broad_key
      ]
    )
]

unknown_broad <- P_B[
  is.na(
    broad_label
  ),
  sort(
    unique(
      broad_group
    )
  )
]

if (length(unknown_broad) > 0L) {

  cat(
    "\n出现未定义的B面板宽分类：\n"
  )

  print(
    unknown_broad
  )

  stop(
    "Figure 1B存在未识别宽分类。",
    call. = FALSE
  )
}

B_COUNTS <- P_B[
  ,
  .(
    n_physical_plants =
      uniqueN(
        plant_id_plot
      )
  ),
  by =
    broad_label
]

cat(
  "\n============================================================\n",
  "B图：互斥宽分类物理电厂数\n",
  "============================================================\n",
  sep = ""
)

print(
  B_COUNTS,
  nrows = Inf
)

# ============================================================
# 8. 保存新的A面板正式绘图源
# ============================================================

fwrite(
  P_A,
  EXPANDED_FILE
)

# ============================================================
# 9. 高区分度颜色
# ============================================================

DETAIL_COLORS_MASTER <- c(
  "Coal" =
    "#4D4D4D",

  "Natural gas" =
    "#E69F00",

  "Oil" =
    "#D55E00",

  "Other fossil" =
    "#7F3C8D",

  "Biomass" =
    "#009E73",

  "Solid waste" =
    "#66A61E",

  "Nuclear" =
    "#CC0066",

  "Geothermal" =
    "#8B5CF6",

  "Wind" =
    "#00A6D6",

  "Solar" =
    "#C79A00",

  "Hydro" =
    "#0066CC",

  "Other" =
    "#9A9A9A",

  "RadNet monitoring station" =
    "#000000"
)

# 只显示数据里真实存在的类型
DETAIL_LEVELS <- DETAIL_LEVELS_MASTER[
  DETAIL_LEVELS_MASTER %in%
    A_COUNTS$fuel_label
]

DETAIL_LEVELS <- c(
  DETAIL_LEVELS,
  "RadNet monitoring station"
)

DETAIL_COLORS <- DETAIL_COLORS_MASTER[
  DETAIL_LEVELS
]

DETAIL_SHAPES <- setNames(
  rep(
    16,
    length(
      DETAIL_LEVELS
    )
  ),
  DETAIL_LEVELS
)

DETAIL_SHAPES[
  "RadNet monitoring station"
] <- 17

# ============================================================
# 10. B面板配色
# ============================================================

BROAD_LEVELS_MASTER <- c(
  "Combustion-capable",
  "Non-combustion major energy",
  "Wind / solar / hydro control",
  "Other"
)

BROAD_COLORS_MASTER <- c(
  "Combustion-capable" =
    "#D55E00",

  "Non-combustion major energy" =
    "#009E73",

  "Wind / solar / hydro control" =
    "#0072B2",

  "Other" =
    "#8A8A8A",

  "RadNet monitoring station" =
    "#000000"
)

BROAD_LEVELS <- BROAD_LEVELS_MASTER[
  BROAD_LEVELS_MASTER %in%
    B_COUNTS$broad_label
]

BROAD_LEVELS <- c(
  BROAD_LEVELS,
  "RadNet monitoring station"
)

BROAD_COLORS <- BROAD_COLORS_MASTER[
  BROAD_LEVELS
]

BROAD_SHAPES <- setNames(
  rep(
    16,
    length(
      BROAD_LEVELS
    )
  ),
  BROAD_LEVELS
)

BROAD_SHAPES[
  "RadNet monitoring station"
] <- 17

# ============================================================
# 11. 组合A面板绘图数据
# ============================================================

A_PLANT <- P_A[
  ,
  .(
    longitude =
      plant_longitude,

    latitude =
      plant_latitude,

    legend_label =
      fuel_label,

    point_group =
      "Power plant"
  )
]

A_MONITOR <- M[
  ,
  .(
    longitude =
      monitor_longitude,

    latitude =
      monitor_latitude,

    legend_label =
      "RadNet monitoring station",

    point_group =
      "RadNet"
  )
]

A_DATA <- rbindlist(
  list(
    A_PLANT,
    A_MONITOR
  ),
  use.names = TRUE
)

A_DATA[
  ,
  legend_label :=
    factor(
      legend_label,
      levels =
        DETAIL_LEVELS
    )
]

A_DATA[
  ,
  point_group :=
    factor(
      point_group,
      levels =
        c(
          "Power plant",
          "RadNet"
        )
    )
]

# ============================================================
# 12. 组合B面板绘图数据
# ============================================================

B_PLANT <- P_B[
  ,
  .(
    longitude =
      plant_longitude,

    latitude =
      plant_latitude,

    legend_label =
      broad_label,

    point_group =
      "Power plant"
  )
]

B_MONITOR <- M[
  ,
  .(
    longitude =
      monitor_longitude,

    latitude =
      monitor_latitude,

    legend_label =
      "RadNet monitoring station",

    point_group =
      "RadNet"
  )
]

B_DATA <- rbindlist(
  list(
    B_PLANT,
    B_MONITOR
  ),
  use.names = TRUE
)

B_DATA[
  ,
  legend_label :=
    factor(
      legend_label,
      levels =
        BROAD_LEVELS
    )
]

B_DATA[
  ,
  point_group :=
    factor(
      point_group,
      levels =
        c(
          "Power plant",
          "RadNet"
        )
    )
]

# ============================================================
# 13. 美国州界底图
# ============================================================

US_MAP <- map_data(
  "state"
)

MAP_THEME <- theme_void(
  base_size = 11
) +
  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 13,
        hjust = 0
      ),

    legend.position =
      "right",

    legend.title =
      element_blank(),

    legend.text =
      element_text(
        size = 9.5
      ),

    legend.key.height =
      unit(
        0.43,
        "cm"
      ),

    legend.key.width =
      unit(
        0.50,
        "cm"
      ),

    plot.margin =
      margin(
        5,
        10,
        5,
        5
      )
  )

# ============================================================
# 14. Figure 1A
#
# 注意：
# color guide手工指定每一个图例符号shape，
# shape guide关闭。
#
# 这样：
# - 发电站 = 彩色圆点
# - RadNet = 黑色三角形
# - 不再出现Duplicated override.aes警告
# ============================================================

A_LEGEND_SHAPES <- unname(
  DETAIL_SHAPES[
    DETAIL_LEVELS
  ]
)

pA <- ggplot() +

  geom_polygon(
    data = US_MAP,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "white",
    color = "#BEBEBE",
    linewidth = 0.32
  ) +

  geom_point(
    data = A_DATA,
    aes(
      x = longitude,
      y = latitude,
      color = legend_label,
      shape = legend_label,
      alpha = point_group,
      size = point_group
    )
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

  scale_color_manual(
    values =
      DETAIL_COLORS,
    breaks =
      DETAIL_LEVELS,
    drop =
      FALSE
  ) +

  scale_shape_manual(
    values =
      DETAIL_SHAPES,
    breaks =
      DETAIL_LEVELS,
    drop =
      FALSE
  ) +

  scale_alpha_manual(
    values = c(
      "Power plant" = 0.58,
      "RadNet" = 1
    ),
    guide = "none"
  ) +

  scale_size_manual(
    values = c(
      "Power plant" = 0.88,
      "RadNet" = 2.25
    ),
    guide = "none"
  ) +

  guides(
    color =
      guide_legend(
        order = 1,
        override.aes =
          list(
            shape =
              A_LEGEND_SHAPES,
            alpha =
              1,
            size =
              3
          )
      ),

    shape =
      "none"
  ) +

  labs(
    title =
      "A. Power plants by energy type"
  ) +

  MAP_THEME

# ============================================================
# 15. Figure 1B
# ============================================================

B_LEGEND_SHAPES <- unname(
  BROAD_SHAPES[
    BROAD_LEVELS
  ]
)

pB <- ggplot() +

  geom_polygon(
    data = US_MAP,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "white",
    color = "#BEBEBE",
    linewidth = 0.32
  ) +

  geom_point(
    data = B_DATA,
    aes(
      x = longitude,
      y = latitude,
      color = legend_label,
      shape = legend_label,
      alpha = point_group,
      size = point_group
    )
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

  scale_color_manual(
    values =
      BROAD_COLORS,
    breaks =
      BROAD_LEVELS,
    drop =
      FALSE
  ) +

  scale_shape_manual(
    values =
      BROAD_SHAPES,
    breaks =
      BROAD_LEVELS,
    drop =
      FALSE
  ) +

  scale_alpha_manual(
    values = c(
      "Power plant" = 0.62,
      "RadNet" = 1
    ),
    guide = "none"
  ) +

  scale_size_manual(
    values = c(
      "Power plant" = 0.92,
      "RadNet" = 2.25
    ),
    guide = "none"
  ) +

  guides(
    color =
      guide_legend(
        order = 1,
        override.aes =
          list(
            shape =
              B_LEGEND_SHAPES,
            alpha =
              1,
            size =
              3
          )
      ),

    shape =
      "none"
  ) +

  labs(
    title =
      "B. Broad power-plant classes"
  ) +

  MAP_THEME

# ============================================================
# 16. 整体Figure 1
# ============================================================

FIG1 <- pA / pB +

  plot_annotation(
    title =
      "Figure 1. Spatial distribution of U.S. power plants and RadNet monitoring stations",

    theme =
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 16,
            hjust = 0
          ),

        plot.margin =
          margin(
            8,
            8,
            8,
            8
          )
      )
  )

# ============================================================
# 17. 保存
# ============================================================

ggsave(
  filename =
    OUT_PNG,
  plot =
    FIG1,
  width =
    14,
  height =
    11,
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
    filename =
      OUT_PDF,
    plot =
      FIG1,
    width =
      14,
    height =
      11,
    units =
      "in",
    device =
      cairo_pdf
  )

} else {

  ggsave(
    filename =
      OUT_PDF,
    plot =
      FIG1,
    width =
      14,
    height =
      11,
    units =
      "in",
    device =
      "pdf"
  )
}

# ============================================================
# 18. 最终审计
# ============================================================

n_other_fossil <- A_COUNTS[
  fuel_label ==
    "Other fossil",
  n_physical_plants
]

if (
  length(
    n_other_fossil
  ) ==
    0L
) {
  n_other_fossil <- 0L
}

AUDIT <- rbindlist(
  list(

    data.table(
      section =
        "overall",

      item =
        c(
          "unique_physical_powerplants",
          "RadNet_monitors",
          "Figure1A_nonexclusive_energy_types",
          "Figure1B_mutually_exclusive_broad_classes",
          "Other_fossil_physical_plants",
          "Mixed_multiple_fuels_removed_from_A",
          "PNG_exists",
          "PDF_exists"
        ),

      value =
        c(
          nrow(P0),
          nrow(M),
          nrow(A_COUNTS),
          nrow(B_COUNTS),
          n_other_fossil,
          as.integer(
            !"Mixed / multiple fuels" %in%
              DETAIL_LEVELS
          ),
          as.integer(
            file.exists(
              OUT_PNG
            )
          ),
          as.integer(
            file.exists(
              OUT_PDF
            )
          )
        )
    ),

    A_COUNTS[
      ,
      .(
        section =
          "Figure1A_energy_type",

        item =
          fuel_label,

        value =
          n_physical_plants
      )
    ],

    B_COUNTS[
      ,
      .(
        section =
          "Figure1B_broad_class",

        item =
          broad_label,

        value =
          n_physical_plants
      )
    ]
  ),
  use.names = TRUE,
  fill = TRUE
)

fwrite(
  AUDIT,
  AUDIT_FILE
)

cat(
  "\n============================================================\n",
  "Figure 1 FINAL AUDIT\n",
  "============================================================\n",
  sep = ""
)

print(
  AUDIT,
  nrows = Inf
)

if (
  n_other_fossil ==
    0L
) {

  stop(
    paste0(
      "重新展开all_fuel_groups后Other fossil仍为0。\n",
      "这将说明源数据中确实没有Other fossil，而不是绘图问题。"
    ),
    call. = FALSE
  )
}

if (
  !file.exists(
    OUT_PNG
  ) ||
  !file.exists(
    OUT_PDF
  )
) {

  stop(
    "Figure 1输出失败。",
    call. = FALSE
  )
}

cat(
  "\n============================================================\n",
  "Figure 1 NON-EXCLUSIVE ENERGY TYPE REDRAW: PASS\n",
  "============================================================\n",
  "Other fossil physical plants = ",
  n_other_fossil,
  "\n",
  "PNG = ",
  OUT_PNG,
  "\n",
  "PDF = ",
  OUT_PDF,
  "\n",
  "Audit = ",
  AUDIT_FILE,
  "\n",
  sep = ""
)

