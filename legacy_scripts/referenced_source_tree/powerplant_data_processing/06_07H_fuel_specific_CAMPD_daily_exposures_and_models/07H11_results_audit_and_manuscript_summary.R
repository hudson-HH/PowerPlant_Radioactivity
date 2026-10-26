#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999,
  width = 220
)

required_packages <- c(
  "data.table",
  "ggplot2"
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
})

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

direction_text <- function(x) {
  fifelse(
    !is.finite(x),
    NA_character_,
    fifelse(
      x > 0,
      "正向",
      fifelse(
        x < 0,
        "负向",
        "零"
      )
    )
  )
}

ci_text <- function(
  low,
  high
) {
  fifelse(
    !is.finite(low) |
      !is.finite(high),
    NA_character_,
    fifelse(
      low > 0,
      "95%置信区间完全高于0",
      fifelse(
        high < 0,
        "95%置信区间完全低于0",
        "95%置信区间包含0"
      )
    )
  )
}

same_sign_text <- function(
  x,
  y
) {
  fifelse(
    is.finite(x) &
      is.finite(y),
    fifelse(
      sign(x) ==
        sign(y),
      "方向一致",
      "方向相反"
    ),
    "无法比较"
  )
}

safe_max <- function(x) {
  z <- x[
    is.finite(x)
  ]

  if (
    length(z) ==
      0L
  ) {
    return(
      NA_real_
    )
  }

  max(z)
}

safe_median <- function(x) {
  z <- x[
    is.finite(x)
  ]

  if (
    length(z) ==
      0L
  ) {
    return(
      NA_real_
    )
  }

  median(z)
}

safe_min <- function(x) {
  z <- x[
    is.finite(x)
  ]

  if (
    length(z) ==
      0L
  ) {
    return(
      NA_real_
    )
  }

  min(z)
}

assert_columns <- function(
  DT,
  cols,
  object_name
) {
  miss <- setdiff(
    cols,
    names(DT)
  )

  if (
    length(miss) >
      0L
  ) {
    stopf(
      object_name,
      " 缺少字段：",
      paste(
        miss,
        collapse = ", "
      )
    )
  }
}

read_required <- function(
  path,
  required_cols = NULL
) {
  if (!file.exists(path)) {
    stopf(
      "缺少07H10正式输出：\n",
      path
    )
  }

  DT <- fread(
    file = path,
    showProgress = FALSE
  )

  if (
    !is.null(
      required_cols
    )
  ) {
    assert_columns(
      DT,
      required_cols,
      basename(path)
    )
  }

  DT
}

branch_cn <- function(x) {
  fifelse(
    x ==
      "filter_period",
    "滤膜采样期",
    fifelse(
      x ==
        "quasi_daily",
      "准逐日",
      x
    )
  )
}

metric_cn <- function(x) {
  fifelse(
    x ==
      "count",
    "电厂数量",
    fifelse(
      x ==
        "heat",
      "实际热输入",
      fifelse(
        x ==
          "proximity",
        "Σ1/d",
        fifelse(
          x ==
            "heat_proximity",
          "ΣHeat/d",
          x
        )
      )
    )
  )
}

save_plot_pair <- function(
  p,
  base_file,
  width = 12,
  height = 8
) {
  png_file <- paste0(
    base_file,
    ".png"
  )

  pdf_file <- paste0(
    base_file,
    ".pdf"
  )

  ggsave(
    filename = png_file,
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  if (
    capabilities(
      "cairo"
    )
  ) {
    ggsave(
      filename = pdf_file,
      plot = p,
      width = width,
      height = height,
      units = "in",
      device = cairo_pdf
    )
  } else {
    ggsave(
      filename = pdf_file,
      plot = p,
      width = width,
      height = height,
      units = "in",
      device = "pdf"
    )
  }

  c(
    png_file,
    pdf_file
  )
}

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

D_MAIN <- file.path(
  H10,
  "01_主分析_监测站年份固定效应"
)

D_STRICT <- file.path(
  H10,
  "02_严格时间固定效应"
)

D_WIND <- file.path(
  H10,
  "03_连续风速分析"
)

D_DIRECTION <- file.path(
  H10,
  "04_上风向下风向直接比较"
)

D_WEIGHTED <- file.path(
  H10,
  "05_距离加权敏感性"
)

D_JOINT <- file.path(
  H10,
  "06_数量热输入联合模型"
)

D_PM <- file.path(
  H10,
  "07_细颗粒物敏感性"
)

D_LAG <- file.path(
  H10,
  "08_准逐日滞后分析"
)

D_DECOMP <- file.path(
  H10,
  "09_站点年份内外分解"
)

D_SUMMARY <- file.path(
  H10,
  "10_综合汇总"
)

H10_HARD_FILE <- file.path(
  D_SUMMARY,
  "04_07H10_hard_checks.csv"
)

OUT <- file.path(
  HROOT,
  "07H11_结果审计与论文结果整理"
)

TABLE_DIR <- file.path(
  OUT,
  "01_论文表"
)

PLOTDATA_DIR <- file.path(
  OUT,
  "02_绘图数据"
)

FIG_DIR <- file.path(
  OUT,
  "03_论文图"
)

SUPPFIG_DIR <- file.path(
  OUT,
  "04_补充图"
)

AUDIT_DIR <- file.path(
  OUT,
  "05_审计"
)

for (
  d in c(
    OUT,
    TABLE_DIR,
    PLOTDATA_DIR,
    FIG_DIR,
    SUPPFIG_DIR,
    AUDIT_DIR
  )
) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

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

FUELS <- c(
  "coal",
  "gas",
  "oil",
  "other_fossil",
  "biomass",
  "solid_waste"
)

FUEL_CN <- c(
  coal = "煤",
  gas = "天然气",
  oil = "燃油",
  other_fossil = "其他化石燃料",
  biomass = "生物质",
  solid_waste = "固体废物"
)

REPRESENTATIVE_RADII <- c(
  20,
  50,
  100,
  200
)

cat(
  "\n============================================================\n",
  "07H11：读取并验证07H10\n",
  "============================================================\n",
  sep = ""
)

H10_HARD <- read_required(
  H10_HARD_FILE,
  c(
    "check",
    "observed",
    "expected",
    "status"
  )
)

if (
  H10_HARD[
    status !=
      "PASS",
    .N
  ] >
    0L
) {
  print(
    H10_HARD[
      status !=
        "PASS"
    ]
  )

  stopf(
    "07H10存在未通过的硬性检查，禁止07H11继续。"
  )
}

MAIN_STATUS <- read_required(
  file.path(
    D_MAIN,
    "01_全部半径主模型状态.csv"
  )
)

MAIN_EFFECTS <- read_required(
  file.path(
    D_MAIN,
    "02_全部半径主模型代表风速效应.csv"
  ),
  c(
    "radius_km",
    "branch",
    "fuel",
    "metric",
    "wind_reference",
    "wind_ms",
    "estimate",
    "ci95_low",
    "ci95_high"
  )
)

MAIN_CURVE <- read_required(
  file.path(
    D_WIND,
    "01_主模型连续风速效应曲线.csv"
  ),
  c(
    "radius_km",
    "branch",
    "fuel",
    "metric",
    "wind_ms",
    "estimate",
    "ci95_low",
    "ci95_high"
  )
)

MAIN_WIND <- read_required(
  file.path(
    D_WIND,
    "02_主模型实际风速分布效应审计.csv"
  )
)

STRICT_STATUS <- read_required(
  file.path(
    D_STRICT,
    "01_全部半径严格时间固定效应模型状态.csv"
  )
)

STRICT_EFFECTS <- read_required(
  file.path(
    D_STRICT,
    "02_全部半径严格时间固定效应结果.csv"
  )
)

STRICT_CURVE <- read_required(
  file.path(
    D_WIND,
    "03_严格时间固定效应连续风速曲线.csv"
  )
)

STRICT_WIND <- read_required(
  file.path(
    D_WIND,
    "04_严格时间固定效应实际风速分布审计.csv"
  )
)

DIRECTION_STATUS <- read_required(
  file.path(
    D_DIRECTION,
    "01_上风向下风向联合模型状态.csv"
  )
)

DIRECTION_EFFECTS <- read_required(
  file.path(
    D_DIRECTION,
    "02_上风向下风向及直接差值检验.csv"
  )
)

WEIGHTED_STATUS <- read_required(
  file.path(
    D_WEIGHTED,
    "01_距离加权固定效应模型状态.csv"
  )
)

WEIGHTED_EFFECTS <- read_required(
  file.path(
    D_WEIGHTED,
    "02_距离加权全部结果.csv"
  )
)

JOINT_STATUS <- read_required(
  file.path(
    D_JOINT,
    "01_Count_Heat联合模型状态.csv"
  )
)

JOINT_EFFECTS <- read_required(
  file.path(
    D_JOINT,
    "02_Count_Heat联合模型效应.csv"
  )
)

JOINT_COLL <- read_required(
  file.path(
    D_JOINT,
    "03_Count_Heat站点年份内相关性与VIF.csv"
  )
)

PM_STATUS <- read_required(
  file.path(
    D_PM,
    "01_不调整PM2.5模型状态.csv"
  )
)

PM_EFFECTS <- read_required(
  file.path(
    D_PM,
    "02_不调整PM2.5同样本结果.csv"
  )
)

LAG_STATUS <- read_required(
  file.path(
    D_LAG,
    "01_准逐日滞后模型状态.csv"
  )
)

LAG_EFFECTS <- read_required(
  file.path(
    D_LAG,
    "02_准逐日滞后全部结果.csv"
  )
)

DECOMP_STATUS <- read_required(
  file.path(
    D_DECOMP,
    "01_站点年份内外分解模型状态.csv"
  )
)

DECOMP_EFFECTS <- read_required(
  file.path(
    D_DECOMP,
    "02_站点年份内外分解全部结果.csv"
  )
)

OLD_WITHIN_BETWEEN <- read_required(
  file.path(
    D_DECOMP,
    "03_冻结原模型与站点年份内外方向解释.csv"
  )
)

H10_INTEGRATED <- read_required(
  file.path(
    D_SUMMARY,
    "01_全部半径综合结果主表.csv"
  )
)

H10_MODEL_STATUS <- read_required(
  file.path(
    D_SUMMARY,
    "02_全部模型状态汇总.csv"
  )
)

H10_FAMILY_AUDIT <- read_required(
  file.path(
    D_SUMMARY,
    "03_分析家族完整性审计.csv"
  )
)

if (
  H10_FAMILY_AUDIT[
    status !=
      "PASS",
    .N
  ] >
    0L
) {
  stopf(
    "07H10分析家族完整性审计不是全PASS。"
  )
}

ALL_STATUS <- rbindlist(
  list(
    MAIN_STATUS[
      ,
      .(
        analysis_family =
          "主分析",
        status
      )
    ],
    STRICT_STATUS[
      ,
      .(
        analysis_family =
          "严格时间固定效应",
        status
      )
    ],
    DIRECTION_STATUS[
      ,
      .(
        analysis_family =
          "上风向下风向联合模型",
        status
      )
    ],
    WEIGHTED_STATUS[
      ,
      .(
        analysis_family =
          "距离加权",
        status
      )
    ],
    JOINT_STATUS[
      ,
      .(
        analysis_family =
          "Count+Heat联合模型",
        status
      )
    ],
    PM_STATUS[
      ,
      .(
        analysis_family =
          "PM2.5敏感性",
        status
      )
    ],
    DECOMP_STATUS[
      ,
      .(
        analysis_family =
          "站点年份内外分解",
        status
      )
    ],
    LAG_STATUS[
      ,
      .(
        analysis_family =
          "准逐日滞后",
        status
      )
    ]
  ),
  use.names = TRUE
)

if (
  ALL_STATUS[
    status ==
      "FAIL",
    .N
  ] >
    0L
) {
  stopf(
    "07H10存在FAIL模型，禁止进行论文结果整理。"
  )
}

cat(
  "07H10全部正式结果输入：PASS\n"
)

# ============================================================
# 一、标准化几个核心结果表
# ============================================================

MAIN_MEDIAN <- MAIN_EFFECTS[
  wind_reference ==
    "样本中位数"
]

MAIN_MEAN <- MAIN_EFFECTS[
  wind_reference ==
    "样本均值"
]

STRICT_MEDIAN <- STRICT_EFFECTS[
  wind_reference ==
    "样本中位数"
]

PM_MEDIAN <- PM_EFFECTS[
  wind_reference ==
    "样本中位数"
]

DIRECTION_MEDIAN <- DIRECTION_EFFECTS[
  wind_reference ==
    "样本中位数"
]

DIRECTION_DIFF <- DIRECTION_MEDIAN[
  component ==
    "上风向减下风向"
]

WEIGHTED_MEDIAN <- WEIGHTED_EFFECTS[
  wind_reference ==
    "样本中位数"
]

WEIGHTED_UP <- WEIGHTED_MEDIAN[
  direction ==
    "upwind"
]

WEIGHTED_UP[
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

JOINT_MEDIAN <- JOINT_EFFECTS[
  wind_reference ==
    "样本中位数"
]

JOINT_MEDIAN[
  ,
  metric :=
    fifelse(
      component ==
        "联合模型_Count",
      "count",
      fifelse(
        component ==
          "联合模型_Heat",
        "heat",
        NA_character_
      )
    )
]

LAG_MEDIAN <- LAG_EFFECTS[
  wind_reference ==
    "样本中位数"
]

MAIN_MEDIAN[
  ,
  `:=`(
    branch_label =
      branch_cn(branch),
    fuel_label =
      unname(
        FUEL_CN[fuel]
      ),
    metric_label =
      metric_cn(metric)
  )
]

MAIN_MEAN[
  ,
  `:=`(
    branch_label =
      branch_cn(branch),
    fuel_label =
      unname(
        FUEL_CN[fuel]
      ),
    metric_label =
      metric_cn(metric)
  )
]

STRICT_MEDIAN[
  ,
  `:=`(
    branch_label =
      branch_cn(branch),
    fuel_label =
      unname(
        FUEL_CN[fuel]
      ),
    metric_label =
      metric_cn(metric)
  )
]

# ============================================================
# 二、主模型 vs 严格时间固定效应
# ============================================================

MAIN_STRICT <- merge(
  MAIN_MEDIAN[
    ,
    .(
      branch,
      branch_label,
      fuel,
      fuel_label,
      radius_km,
      metric,
      metric_label,
      unit,
      main_wind_ms =
        wind_ms,
      main_estimate =
        estimate,
      main_se =
        se,
      main_ci_low =
        ci95_low,
      main_ci_high =
        ci95_high,
      main_p =
        p_wald,
      main_direction =
        direction,
      main_ci_status =
        ci_status
    )
  ],
  STRICT_MEDIAN[
    ,
    .(
      branch,
      fuel,
      radius_km,
      metric,
      strict_wind_ms =
        wind_ms,
      strict_estimate =
        estimate,
      strict_se =
        se,
      strict_ci_low =
        ci95_low,
      strict_ci_high =
        ci95_high,
      strict_p =
        p_wald,
      strict_direction =
        direction,
      strict_ci_status =
        ci_status
    )
  ],
  by = c(
    "branch",
    "fuel",
    "radius_km",
    "metric"
  ),
  all = TRUE,
  sort = FALSE
)

MAIN_STRICT[
  ,
  main_vs_strict :=
    same_sign_text(
      main_estimate,
      strict_estimate
    )
]

# ============================================================
# 三、PM2.5同样本敏感性
# ============================================================

PM_COMPARE <- merge(
  MAIN_MEDIAN[
    ,
    .(
      branch,
      branch_label,
      fuel,
      fuel_label,
      radius_km,
      metric,
      metric_label,
      unit,
      adjusted_estimate =
        estimate,
      adjusted_ci_low =
        ci95_low,
      adjusted_ci_high =
        ci95_high,
      adjusted_direction =
        direction
    )
  ],
  PM_MEDIAN[
    ,
    .(
      branch,
      fuel,
      radius_km,
      metric,
      no_pm_estimate =
        estimate,
      no_pm_ci_low =
        ci95_low,
      no_pm_ci_high =
        ci95_high,
      no_pm_direction =
        direction
    )
  ],
  by = c(
    "branch",
    "fuel",
    "radius_km",
    "metric"
  ),
  all = TRUE,
  sort = FALSE
)

PM_COMPARE[
  ,
  pm_adjustment_direction_compare :=
    same_sign_text(
      adjusted_estimate,
      no_pm_estimate
    )
]

# ============================================================
# 四、普通累计暴露 vs 距离加权
# ============================================================

WEIGHTED_COMPARE <- merge(
  MAIN_MEDIAN[
    ,
    .(
      branch,
      branch_label,
      fuel,
      fuel_label,
      radius_km,
      metric,
      metric_label,
      raw_estimate =
        estimate,
      raw_ci_low =
        ci95_low,
      raw_ci_high =
        ci95_high,
      raw_direction =
        direction
    )
  ],
  WEIGHTED_UP[
    !is.na(
      raw_metric
    ),
    .(
      branch,
      fuel,
      radius_km,
      metric =
        raw_metric,
      weighted_metric =
        metric,
      weighted_estimate =
        estimate,
      weighted_ci_low =
        ci95_low,
      weighted_ci_high =
        ci95_high,
      weighted_direction =
        direction
    )
  ],
  by = c(
    "branch",
    "fuel",
    "radius_km",
    "metric"
  ),
  all = TRUE,
  sort = FALSE
)

WEIGHTED_COMPARE[
  ,
  raw_vs_weighted :=
    same_sign_text(
      raw_estimate,
      weighted_estimate
    )
]

# ============================================================
# 五、连续风速曲线分类
# ============================================================

CURVE_PATTERN <- MAIN_CURVE[
  is.finite(
    estimate
  ),
  .(
    n_curve_points =
      .N,
    min_effect =
      safe_min(
        estimate
      ),
    max_effect =
      safe_max(
        estimate
      ),
    min_wind_ms =
      safe_min(
        wind_ms
      ),
    max_wind_ms =
      safe_max(
        wind_ms
      ),
    wind_pattern =
      if (
        all(
          estimate >
            0
        )
      ) {
        "实际主要风速范围内始终正向"
      } else if (
        all(
          estimate <
            0
        )
      ) {
        "实际主要风速范围内始终负向"
      } else {
        "实际主要风速范围内发生方向转换"
      }
  ),
  by = .(
    branch,
    fuel,
    radius_km,
    metric
  )
]

# ============================================================
# 六、上风向 vs 下风向直接检验
# ============================================================

DIRECTION_TABLE <- DIRECTION_DIFF[
  ,
  .(
    branch,
    branch_label =
      branch_cn(branch),
    fuel,
    fuel_label =
      unname(
        FUEL_CN[fuel]
      ),
    radius_km,
    metric,
    metric_label =
      metric_cn(metric),
    wind_ms,
    up_minus_down_estimate =
      estimate,
    up_minus_down_se =
      se,
    up_minus_down_ci_low =
      ci95_low,
    up_minus_down_ci_high =
      ci95_high,
    up_minus_down_p =
      p_wald,
    up_minus_down_direction =
      direction,
    up_minus_down_ci_status =
      ci_status,
    directionality_interpretation =
      fifelse(
        ci95_low >
          0,
        "上风向效应明确高于下风向",
        fifelse(
          ci95_high <
            0,
          "下风向效应明确高于上风向",
          "上风向与下风向差值置信区间包含0"
        )
      )
  )
]

# ============================================================
# 七、Count + Heat联合模型与共线性
# ============================================================

JOINT_TABLE <- merge(
  JOINT_MEDIAN[
    !is.na(
      metric
    ),
    .(
      branch,
      branch_label =
        branch_cn(branch),
      fuel,
      fuel_label =
        unname(
          FUEL_CN[fuel]
        ),
      radius_km,
      metric,
      metric_label =
        metric_cn(metric),
      joint_component =
        component,
      unit,
      wind_ms,
      estimate,
      se,
      ci95_low,
      ci95_high,
      p_wald,
      direction,
      ci_status
    )
  ],
  JOINT_COLL[
    ,
    .(
      branch,
      fuel,
      radius_km,
      within_count_heat_correlation,
      approximate_two_predictor_VIF
    )
  ],
  by = c(
    "branch",
    "fuel",
    "radius_km"
  ),
  all.x = TRUE,
  sort = FALSE
)

JOINT_TABLE[
  ,
  collinearity_category :=
    fifelse(
      is.infinite(
        approximate_two_predictor_VIF
      ),
      "近似VIF为无穷",
      fifelse(
        !is.finite(
          approximate_two_predictor_VIF
        ),
        "无法计算",
        fifelse(
          approximate_two_predictor_VIF >=
            10,
          "近似VIF≥10",
          fifelse(
            approximate_two_predictor_VIF >=
              5,
            "5≤近似VIF<10",
            "近似VIF<5"
          )
        )
      )
    )
]

# ============================================================
# 八、准逐日滞后分析
# ============================================================

LAG_ORDER <- c(
  "lag0",
  "lag1",
  "lag2",
  "lag3",
  "mean_0_1",
  "mean_0_2"
)

LAG_LABEL <- c(
  lag0 =
    "当天",
  lag1 =
    "前1天",
  lag2 =
    "前2天",
  lag3 =
    "前3天",
  mean_0_1 =
    "当天至前1天平均",
  mean_0_2 =
    "当天至前2天平均"
)

LAG_TABLE <- copy(
  LAG_MEDIAN
)

LAG_TABLE[
  ,
  `:=`(
    branch_label =
      "准逐日",
    fuel_label =
      unname(
        FUEL_CN[fuel]
      ),
    metric_label =
      metric_cn(metric),
    lag_label =
      unname(
        LAG_LABEL[lag_spec]
      ),
    lag_order =
      match(
        lag_spec,
        LAG_ORDER
      )
  )
]

setorder(
  LAG_TABLE,
  fuel,
  metric,
  radius_km,
  lag_order
)

LAG_PATTERN <- LAG_TABLE[
  ,
  .(
    n_lag_specs =
      .N,
    n_positive =
      sum(
        estimate >
          0,
        na.rm = TRUE
      ),
    n_negative =
      sum(
        estimate <
          0,
        na.rm = TRUE
      ),
    n_CI_above_zero =
      sum(
        ci95_low >
          0,
        na.rm = TRUE
      ),
    n_CI_below_zero =
      sum(
        ci95_high <
          0,
        na.rm = TRUE
      ),
    n_CI_include_zero =
      sum(
        ci95_low <=
          0 &
        ci95_high >=
          0,
        na.rm = TRUE
      ),
    lag_direction_sequence =
      paste(
        paste0(
          lag_label,
          ":",
          direction
        ),
        collapse = " | "
      )
  ),
  by = .(
    fuel,
    fuel_label,
    radius_km,
    metric,
    metric_label
  )
]

# ============================================================
# 九、原模型 vs 站点年份内/间
# ============================================================

DECOMP_TABLE <- copy(
  OLD_WITHIN_BETWEEN
)

DECOMP_TABLE[
  ,
  `:=`(
    branch_label =
      branch_cn(branch),
    fuel_label =
      unname(
        FUEL_CN[fuel]
      ),
    metric_label =
      metric_cn(metric),
    old_vs_within =
      same_sign_text(
        old_estimate_wind1,
        within_estimate_wind1
      ),
    within_vs_between =
      same_sign_text(
        within_estimate_wind1,
        between_estimate_wind1
      )
  )
]

DECOMP_TABLE[
  ,
  structural_pattern :=
    fifelse(
      is.finite(
        old_estimate_wind1
      ) &
      is.finite(
        within_estimate_wind1
      ) &
      is.finite(
        between_estimate_wind1
      ) &
      sign(
        within_estimate_wind1
      ) !=
        sign(
          between_estimate_wind1
        ),
      "站点年份内关联与站点年份间关联方向相反",
      fifelse(
        is.finite(
          within_estimate_wind1
        ) &
        is.finite(
          between_estimate_wind1
        ) &
        sign(
          within_estimate_wind1
        ) ==
          sign(
            between_estimate_wind1
          ),
        "站点年份内关联与站点年份间关联方向一致",
        "无法完整判断"
      )
    )
]

# ============================================================
# 十、滤膜采样期 vs 准逐日
# ============================================================

FILTER_MAIN <- MAIN_MEDIAN[
  branch ==
    "filter_period",
  .(
    fuel,
    fuel_label,
    radius_km,
    metric,
    metric_label,
    filter_estimate =
      estimate,
    filter_ci_low =
      ci95_low,
    filter_ci_high =
      ci95_high,
    filter_direction =
      direction,
    filter_ci_status =
      ci_status
  )
]

QUASI_MAIN <- MAIN_MEDIAN[
  branch ==
    "quasi_daily",
  .(
    fuel,
    radius_km,
    metric,
    quasi_estimate =
      estimate,
    quasi_ci_low =
      ci95_low,
    quasi_ci_high =
      ci95_high,
    quasi_direction =
      direction,
    quasi_ci_status =
      ci_status
  )
]

TIME_SCALE_COMPARE <- merge(
  FILTER_MAIN,
  QUASI_MAIN,
  by = c(
    "fuel",
    "radius_km",
    "metric"
  ),
  all = TRUE,
  sort = FALSE
)

TIME_SCALE_COMPARE[
  ,
  filter_vs_quasi :=
    same_sign_text(
      filter_estimate,
      quasi_estimate
    )
]

TIME_SCALE_SUMMARY <- TIME_SCALE_COMPARE[
  ,
  .(
    n_comparable_radii =
      sum(
        is.finite(
          filter_estimate
        ) &
        is.finite(
          quasi_estimate
        )
      ),
    n_same_direction =
      sum(
        sign(
          filter_estimate
        ) ==
          sign(
            quasi_estimate
          ),
        na.rm = TRUE
      ),
    n_opposite_direction =
      sum(
        sign(
          filter_estimate
        ) !=
          sign(
            quasi_estimate
          ),
        na.rm = TRUE
      ),
    filter_positive =
      sum(
        filter_estimate >
          0,
        na.rm = TRUE
      ),
    filter_negative =
      sum(
        filter_estimate <
          0,
        na.rm = TRUE
      ),
    quasi_positive =
      sum(
        quasi_estimate >
          0,
        na.rm = TRUE
      ),
    quasi_negative =
      sum(
        quasi_estimate <
          0,
        na.rm = TRUE
      )
  ),
  by = .(
    fuel,
    fuel_label,
    metric,
    metric_label
  )
]

# ============================================================
# 十一、综合燃料层面审计
# ============================================================

MAIN_GROUP <- MAIN_MEDIAN[
  ,
  .(
    n_main_radii =
      .N,
    n_main_positive =
      sum(
        estimate >
          0,
        na.rm = TRUE
      ),
    n_main_negative =
      sum(
        estimate <
          0,
        na.rm = TRUE
      ),
    n_main_CI_above_zero =
      sum(
        ci95_low >
          0,
        na.rm = TRUE
      ),
    n_main_CI_below_zero =
      sum(
        ci95_high <
          0,
        na.rm = TRUE
      ),
    n_main_CI_include_zero =
      sum(
        ci95_low <=
          0 &
        ci95_high >=
          0,
        na.rm = TRUE
      )
  ),
  by = .(
    branch,
    branch_label,
    fuel,
    fuel_label,
    metric,
    metric_label
  )
]

STRICT_GROUP <- MAIN_STRICT[
  ,
  .(
    n_main_strict_comparable =
      sum(
        is.finite(
          main_estimate
        ) &
        is.finite(
          strict_estimate
        )
      ),
    n_main_strict_same_direction =
      sum(
        sign(
          main_estimate
        ) ==
          sign(
            strict_estimate
        ),
        na.rm = TRUE
      ),
    n_main_strict_opposite_direction =
      sum(
        sign(
          main_estimate
        ) !=
          sign(
            strict_estimate
        ),
        na.rm = TRUE
      )
  ),
  by = .(
    branch,
    fuel,
    metric
  )
]

PM_GROUP <- PM_COMPARE[
  ,
  .(
    n_PM_comparable =
      sum(
        is.finite(
          adjusted_estimate
        ) &
        is.finite(
          no_pm_estimate
        )
      ),
    n_PM_same_direction =
      sum(
        sign(
          adjusted_estimate
        ) ==
          sign(
            no_pm_estimate
        ),
        na.rm = TRUE
      ),
    n_PM_opposite_direction =
      sum(
        sign(
          adjusted_estimate
        ) !=
          sign(
            no_pm_estimate
        ),
        na.rm = TRUE
      )
  ),
  by = .(
    branch,
    fuel,
    metric
  )
]

WEIGHT_GROUP <- WEIGHTED_COMPARE[
  ,
  .(
    n_weighted_comparable =
      sum(
        is.finite(
          raw_estimate
        ) &
        is.finite(
          weighted_estimate
        )
      ),
    n_weighted_same_direction =
      sum(
        sign(
          raw_estimate
        ) ==
          sign(
            weighted_estimate
        ),
        na.rm = TRUE
      ),
    n_weighted_opposite_direction =
      sum(
        sign(
          raw_estimate
        ) !=
          sign(
            weighted_estimate
        ),
        na.rm = TRUE
      )
  ),
  by = .(
    branch,
    fuel,
    metric
  )
]

DIR_GROUP <- DIRECTION_TABLE[
  ,
  .(
    n_direction_comparable =
      .N,
    n_up_minus_down_positive =
      sum(
        up_minus_down_estimate >
          0,
        na.rm = TRUE
      ),
    n_up_minus_down_negative =
      sum(
        up_minus_down_estimate <
          0,
        na.rm = TRUE
      ),
    n_up_minus_down_CI_above_zero =
      sum(
        up_minus_down_ci_low >
          0,
        na.rm = TRUE
      ),
    n_up_minus_down_CI_below_zero =
      sum(
        up_minus_down_ci_high <
          0,
        na.rm = TRUE
      )
  ),
  by = .(
    branch,
    fuel,
    metric
  )
]

WIND_GROUP <- CURVE_PATTERN[
  ,
  .(
    n_wind_curve_radii =
      .N,
    n_wind_always_positive =
      sum(
        wind_pattern ==
          "实际主要风速范围内始终正向"
      ),
    n_wind_always_negative =
      sum(
        wind_pattern ==
          "实际主要风速范围内始终负向"
      ),
    n_wind_direction_change =
      sum(
        wind_pattern ==
          "实际主要风速范围内发生方向转换"
      )
  ),
  by = .(
    branch,
    fuel,
    metric
  )
]

DECOMP_GROUP <- DECOMP_TABLE[
  ,
  .(
    n_decomp_comparable =
      sum(
        is.finite(
          within_estimate_wind1
        ) &
        is.finite(
          between_estimate_wind1
        )
      ),
    n_within_between_opposite =
      sum(
        sign(
          within_estimate_wind1
        ) !=
          sign(
            between_estimate_wind1
        ),
        na.rm = TRUE
      ),
    n_old_within_opposite =
      sum(
        sign(
          old_estimate_wind1
        ) !=
          sign(
            within_estimate_wind1
        ),
        na.rm = TRUE
      )
  ),
  by = .(
    branch,
    fuel,
    metric
  )
]

VIF_GROUP <- JOINT_COLL[
  ,
  .(
    joint_within_corr_median =
      safe_median(
        within_count_heat_correlation
      ),
    joint_within_corr_max_abs =
      safe_max(
        abs(
          within_count_heat_correlation
        )
      ),
    joint_VIF_median =
      safe_median(
        approximate_two_predictor_VIF
      ),
    joint_VIF_max =
      if (
        any(
          is.infinite(
            approximate_two_predictor_VIF
          )
        )
      ) {
        Inf
      } else {
        safe_max(
          approximate_two_predictor_VIF
        )
      },
    n_joint_VIF_ge5 =
      sum(
        approximate_two_predictor_VIF >=
          5,
        na.rm = TRUE
      ),
    n_joint_VIF_ge10 =
      sum(
        approximate_two_predictor_VIF >=
          10,
        na.rm = TRUE
      ),
    n_joint_VIF_infinite =
      sum(
        is.infinite(
          approximate_two_predictor_VIF
        )
      )
  ),
  by = .(
    branch,
    fuel
  )
]

LAG_GROUP <- LAG_TABLE[
  ,
  .(
    n_lag_effects =
      .N,
    n_lag_positive =
      sum(
        estimate >
          0,
        na.rm = TRUE
      ),
    n_lag_negative =
      sum(
        estimate <
          0,
        na.rm = TRUE
      ),
    n_lag_CI_above_zero =
      sum(
        ci95_low >
          0,
        na.rm = TRUE
      ),
    n_lag_CI_below_zero =
      sum(
        ci95_high <
          0,
        na.rm = TRUE
      )
  ),
  by = .(
    fuel,
    metric
  )
]

FUEL_SUMMARY <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = intersect(
        names(x),
        names(y)
      ),
      all.x = TRUE,
      sort = FALSE
    )
  },
  list(
    MAIN_GROUP,
    STRICT_GROUP,
    PM_GROUP,
    WEIGHT_GROUP,
    DIR_GROUP,
    WIND_GROUP,
    DECOMP_GROUP
  )
)

FUEL_SUMMARY <- merge(
  FUEL_SUMMARY,
  VIF_GROUP,
  by = c(
    "branch",
    "fuel"
  ),
  all.x = TRUE,
  sort = FALSE
)

FUEL_SUMMARY <- merge(
  FUEL_SUMMARY,
  LAG_GROUP,
  by = c(
    "fuel",
    "metric"
  ),
  all.x = TRUE,
  sort = FALSE
)

FUEL_SUMMARY[
  ,
  main_radius_pattern :=
    fifelse(
      n_main_positive /
        n_main_radii >=
        0.8,
      "多数预设半径为正向",
      fifelse(
        n_main_negative /
          n_main_radii >=
          0.8,
        "多数预设半径为负向",
        "不同预设半径间方向混合"
      )
    )
]

FUEL_SUMMARY[
  ,
  strict_time_robustness :=
    fifelse(
      n_main_strict_comparable >
        0,
      paste0(
        n_main_strict_same_direction,
        "/",
        n_main_strict_comparable,
        "个半径与严格时间控制方向一致"
      ),
      "无可比较结果"
    )
]

FUEL_SUMMARY[
  ,
  weighted_robustness :=
    fifelse(
      n_weighted_comparable >
        0,
      paste0(
        n_weighted_same_direction,
        "/",
        n_weighted_comparable,
        "个半径与距离加权方向一致"
      ),
      "无可比较结果"
    )
]

FUEL_SUMMARY[
  ,
  PM_sensitivity_text :=
    fifelse(
      n_PM_comparable >
        0,
      paste0(
        n_PM_same_direction,
        "/",
        n_PM_comparable,
        "个半径在调整/不调整PM2.5时方向一致"
      ),
      "无可比较结果"
    )
]

FUEL_SUMMARY[
  ,
  directionality_text :=
    fifelse(
      n_direction_comparable >
        0,
      paste0(
        n_up_minus_down_CI_above_zero,
        "/",
        n_direction_comparable,
        "个半径上风向减下风向的95%置信区间完全高于0"
      ),
      "无可比较结果"
    )
]

FUEL_SUMMARY[
  ,
  wind_text :=
    paste0(
      "始终正向",
      n_wind_always_positive,
      "；始终负向",
      n_wind_always_negative,
      "；实际风速范围内方向转换",
      n_wind_direction_change
    )
]

setorder(
  FUEL_SUMMARY,
  branch,
  fuel,
  metric
)

# ============================================================
# 十二、论文主表
# ============================================================

PAPER_MAIN <- MAIN_MEDIAN[
  ,
  .(
    分析分支 =
      branch_cn(branch),
    燃料 =
      unname(
        FUEL_CN[fuel]
      ),
    燃料代码 =
      fuel,
    累计半径_km =
      radius_km,
    暴露指标 =
      metric_cn(metric),
    效应单位 =
      unit,
    样本中位风速_m_s =
      wind_ms,
    效应估计_mBq_m3 =
      estimate,
    标准误 =
      se,
    `95%CI下限` =
      ci95_low,
    `95%CI上限` =
      ci95_high,
    P值 =
      p_wald,
    调整关联方向 =
      direction,
    置信区间状态 =
      ci_status
  )
]

PAPER_MEAN_WIND <- MAIN_MEAN[
  ,
  .(
    分析分支 =
      branch_cn(branch),
    燃料 =
      unname(
        FUEL_CN[fuel]
      ),
    燃料代码 =
      fuel,
    累计半径_km =
      radius_km,
    暴露指标 =
      metric_cn(metric),
    效应单位 =
      unit,
    样本平均风速_m_s =
      wind_ms,
    效应估计_mBq_m3 =
      estimate,
    `95%CI下限` =
      ci95_low,
    `95%CI上限` =
      ci95_high,
    调整关联方向 =
      direction,
    置信区间状态 =
      ci_status
  )
]

PAPER_STRICT <- MAIN_STRICT[
  ,
  .(
    分析分支 =
      branch_cn(branch),
    燃料 =
      unname(
        FUEL_CN[fuel]
      ),
    累计半径_km =
      radius_km,
    暴露指标 =
      metric_cn(metric),
    主分析效应 =
      main_estimate,
    主分析95CI下限 =
      main_ci_low,
    主分析95CI上限 =
      main_ci_high,
    严格时间控制效应 =
      strict_estimate,
    严格时间控制95CI下限 =
      strict_ci_low,
    严格时间控制95CI上限 =
      strict_ci_high,
    方向比较 =
      main_vs_strict
  )
]

PAPER_DIRECTION <- DIRECTION_TABLE[
  ,
  .(
    分析分支 =
      branch_label,
    燃料 =
      fuel_label,
    累计半径_km =
      radius_km,
    暴露指标 =
      metric_label,
    样本中位风速_m_s =
      wind_ms,
    上风向减下风向效应 =
      up_minus_down_estimate,
    `95%CI下限` =
      up_minus_down_ci_low,
    `95%CI上限` =
      up_minus_down_ci_high,
    P值 =
      up_minus_down_p,
    方向性解释 =
      directionality_interpretation
  )
]

PAPER_WEIGHTED <- WEIGHTED_COMPARE[
  ,
  .(
    分析分支 =
      branch_cn(branch),
    燃料 =
      unname(
        FUEL_CN[fuel]
      ),
    累计半径_km =
      radius_km,
    普通暴露指标 =
      metric_cn(metric),
    普通累计暴露效应 =
      raw_estimate,
    普通累计暴露方向 =
      raw_direction,
    距离加权暴露效应 =
      weighted_estimate,
    距离加权方向 =
      weighted_direction,
    方向比较 =
      raw_vs_weighted
  )
]

PAPER_PM <- PM_COMPARE[
  ,
  .(
    分析分支 =
      branch_cn(branch),
    燃料 =
      unname(
        FUEL_CN[fuel]
      ),
    累计半径_km =
      radius_km,
    暴露指标 =
      metric_cn(metric),
    调整PM2_5效应 =
      adjusted_estimate,
    调整PM2_5方向 =
      adjusted_direction,
    不调整PM2_5效应 =
      no_pm_estimate,
    不调整PM2_5方向 =
      no_pm_direction,
    方向比较 =
      pm_adjustment_direction_compare
  )
]

PAPER_JOINT <- JOINT_TABLE[
  ,
  .(
    分析分支 =
      branch_cn(branch),
    燃料 =
      fuel_label,
    累计半径_km =
      radius_km,
    联合模型分量 =
      joint_component,
    效应单位 =
      unit,
    效应估计 =
      estimate,
    `95%CI下限` =
      ci95_low,
    `95%CI上限` =
      ci95_high,
    P值 =
      p_wald,
    调整关联方向 =
      direction,
    Count_Heat站点年份内相关系数 =
      within_count_heat_correlation,
    近似VIF =
      approximate_two_predictor_VIF,
    共线性类别 =
      collinearity_category
  )
]

PAPER_LAG <- LAG_TABLE[
  ,
  .(
    燃料 =
      fuel_label,
    累计半径_km =
      radius_km,
    暴露指标 =
      metric_label,
    滞后定义 =
      lag_label,
    风速_m_s =
      wind_ms,
    效应估计 =
      estimate,
    `95%CI下限` =
      ci95_low,
    `95%CI上限` =
      ci95_high,
    P值 =
      p_wald,
    调整关联方向 =
      direction,
    置信区间状态 =
      ci_status
  )
]

PAPER_DECOMP <- DECOMP_TABLE[
  ,
  .(
    分析分支 =
      branch_label,
    燃料 =
      fuel_label,
    累计半径_km =
      radius_km,
    暴露指标 =
      metric_label,
    原随机截距模型_1m_s =
      old_estimate_wind1,
    站点年份内效应_1m_s =
      within_estimate_wind1,
    站点年份间效应_1m_s =
      between_estimate_wind1,
    原模型与站内方向比较 =
      old_vs_within,
    站内与站点年份间方向比较 =
      within_vs_between,
    结构解释 =
      structural_pattern
  )
]

PAPER_TIME <- TIME_SCALE_COMPARE[
  ,
  .(
    燃料 =
      fuel_label,
    累计半径_km =
      radius_km,
    暴露指标 =
      metric_label,
    滤膜采样期效应 =
      filter_estimate,
    滤膜采样期方向 =
      filter_direction,
    准逐日效应 =
      quasi_estimate,
    准逐日方向 =
      quasi_direction,
    时间尺度方向比较 =
      filter_vs_quasi
  )
]

fwrite(
  PAPER_MAIN,
  file.path(
    TABLE_DIR,
    "01_主分析_全部预设半径_样本中位风速.csv"
  )
)

fwrite(
  PAPER_MEAN_WIND,
  file.path(
    TABLE_DIR,
    "02_主分析_全部预设半径_样本平均风速.csv"
  )
)

fwrite(
  PAPER_STRICT,
  file.path(
    TABLE_DIR,
    "03_主分析与严格时间固定效应对照.csv"
  )
)

fwrite(
  PAPER_DIRECTION,
  file.path(
    TABLE_DIR,
    "04_上风向减下风向直接差值检验.csv"
  )
)

fwrite(
  PAPER_WEIGHTED,
  file.path(
    TABLE_DIR,
    "05_普通累计暴露与距离加权敏感性.csv"
  )
)

fwrite(
  PAPER_PM,
  file.path(
    TABLE_DIR,
    "06_PM2.5同样本敏感性.csv"
  )
)

fwrite(
  PAPER_JOINT,
  file.path(
    TABLE_DIR,
    "07_Count与Heat联合模型及共线性.csv"
  )
)

fwrite(
  PAPER_LAG,
  file.path(
    TABLE_DIR,
    "08_准逐日滞后分析_全部预设滞后.csv"
  )
)

fwrite(
  PAPER_DECOMP,
  file.path(
    TABLE_DIR,
    "09_原模型与站点年份内外分解.csv"
  )
)

fwrite(
  PAPER_TIME,
  file.path(
    TABLE_DIR,
    "10_滤膜采样期与准逐日结果对照.csv"
  )
)

fwrite(
  FUEL_SUMMARY,
  file.path(
    TABLE_DIR,
    "11_燃料层面全部分析综合审计.csv"
  )
)

fwrite(
  TIME_SCALE_SUMMARY,
  file.path(
    TABLE_DIR,
    "12_滤膜采样期与准逐日时间尺度汇总.csv"
  )
)

fwrite(
  LAG_PATTERN,
  file.path(
    TABLE_DIR,
    "13_滞后模式汇总_不选择最佳滞后.csv"
  )
)

# ============================================================
# 十三、绘图数据
# ============================================================

fwrite(
  MAIN_MEDIAN,
  file.path(
    PLOTDATA_DIR,
    "01_主分析半径效应绘图数据.csv"
  )
)

fwrite(
  MAIN_STRICT,
  file.path(
    PLOTDATA_DIR,
    "02_主分析与严格时间控制绘图数据.csv"
  )
)

fwrite(
  DIRECTION_TABLE,
  file.path(
    PLOTDATA_DIR,
    "03_上风向减下风向绘图数据.csv"
  )
)

fwrite(
  MAIN_CURVE[
    radius_km %in%
      REPRESENTATIVE_RADII
  ],
  file.path(
    PLOTDATA_DIR,
    "04_连续风速效应曲线_预设代表半径.csv"
  )
)

fwrite(
  LAG_TABLE,
  file.path(
    PLOTDATA_DIR,
    "05_准逐日滞后绘图数据.csv"
  )
)

# ============================================================
# 十四、论文图
# ============================================================

PLOT_MANIFEST <- list()
plot_index <- 0L

make_main_radius_plot <- function(
  branch_value,
  metric_value
) {
  D <- MAIN_MEDIAN[
    branch ==
      branch_value &
    metric ==
      metric_value
  ]

  D[
    ,
    radius_factor :=
      factor(
        radius_km,
        levels = RADII
      )
  ]

  p <- ggplot(
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
      width = 0.18
    ) +
    geom_line() +
    geom_point(
      size = 2
    ) +
    facet_wrap(
      ~fuel_label,
      scales = "free_y",
      ncol = 3
    ) +
    labs(
      x = "累计半径（km）",
      y = "总贝塔变化（mBq/m³）",
      title = paste0(
        branch_cn(
          branch_value
        ),
        "：",
        metric_cn(
          metric_value
        ),
        "的监测站×年份固定效应"
      ),
      subtitle = "点为样本中位风速下效应估计；误差线为95%置信区间"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor =
        element_blank(),
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1
        )
    )

  p
}

figure_counter <- 0L

for (
  branch_value in
  c(
    "filter_period",
    "quasi_daily"
  )
) {
  for (
    metric_value in
    c(
      "count",
      "heat"
    )
  ) {
    figure_counter <- figure_counter +
      1L

    p <- make_main_radius_plot(
      branch_value,
      metric_value
    )

    base <- file.path(
      FIG_DIR,
      sprintf(
        "%02d_主分析_%s_%s",
        figure_counter,
        branch_cn(
          branch_value
        ),
        metric_cn(
          metric_value
        )
      )
    )

    files <- save_plot_pair(
      p,
      base,
      12,
      8
    )

    plot_index <- plot_index +
      1L

    PLOT_MANIFEST[[plot_index]] <-
      data.table(
        figure_group =
          "主分析半径效应",
        branch =
          branch_value,
        metric =
          metric_value,
        base_file =
          base,
        png_file =
          files[1],
        pdf_file =
          files[2]
      )
  }
}

make_main_strict_plot <- function(
  branch_value,
  metric_value
) {
  D1 <- MAIN_STRICT[
    branch ==
      branch_value &
    metric ==
      metric_value
  ]

  A <- D1[
    ,
    .(
      fuel_label,
      radius_km,
      model =
        "监测站×年份固定效应",
      estimate =
        main_estimate,
      ci_low =
        main_ci_low,
      ci_high =
        main_ci_high
    )
  ]

  B <- D1[
    ,
    .(
      fuel_label,
      radius_km,
      model =
        "监测站×年份+共同时间固定效应",
      estimate =
        strict_estimate,
      ci_low =
        strict_ci_low,
      ci_high =
        strict_ci_high
    )
  ]

  D <- rbind(
    A,
    B
  )

  D[
    ,
    radius_factor :=
      factor(
        radius_km,
        levels = RADII
      )
  ]

  p <- ggplot(
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
    geom_line(
      position =
        position_dodge(
          width = 0.25
        )
    ) +
    geom_point(
      size = 1.8,
      position =
        position_dodge(
          width = 0.25
        )
    ) +
    facet_wrap(
      ~fuel_label,
      scales = "free_y",
      ncol = 3
    ) +
    labs(
      x = "累计半径（km）",
      y = "总贝塔变化（mBq/m³）",
      linetype = NULL,
      shape = NULL,
      title = paste0(
        branch_cn(
          branch_value
        ),
        "：主分析与严格时间控制比较（",
        metric_cn(
          metric_value
        ),
        "）"
      )
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor =
        element_blank(),
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1
        ),
      legend.position =
        "bottom"
    )

  p
}

for (
  branch_value in
  c(
    "filter_period",
    "quasi_daily"
  )
) {
  for (
    metric_value in
    c(
      "count",
      "heat"
    )
  ) {
    figure_counter <- figure_counter +
      1L

    p <- make_main_strict_plot(
      branch_value,
      metric_value
    )

    base <- file.path(
      FIG_DIR,
      sprintf(
        "%02d_主分析与严格时间控制_%s_%s",
        figure_counter,
        branch_cn(
          branch_value
        ),
        metric_cn(
          metric_value
        )
      )
    )

    files <- save_plot_pair(
      p,
      base,
      12,
      8
    )

    plot_index <- plot_index +
      1L

    PLOT_MANIFEST[[plot_index]] <-
      data.table(
        figure_group =
          "主分析与严格时间控制",
        branch =
          branch_value,
        metric =
          metric_value,
        base_file =
          base,
        png_file =
          files[1],
        pdf_file =
          files[2]
      )
  }
}

make_direction_plot <- function(
  branch_value,
  metric_value
) {
  D <- DIRECTION_TABLE[
    branch ==
      branch_value &
    metric ==
      metric_value
  ]

  D[
    ,
    radius_factor :=
      factor(
        radius_km,
        levels = RADII
      )
  ]

  p <- ggplot(
    D,
    aes(
      x = radius_factor,
      y = up_minus_down_estimate,
      group = 1
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = 2
    ) +
    geom_errorbar(
      aes(
        ymin =
          up_minus_down_ci_low,
        ymax =
          up_minus_down_ci_high
      ),
      width = 0.18
    ) +
    geom_line() +
    geom_point(
      size = 2
    ) +
    facet_wrap(
      ~fuel_label,
      scales = "free_y",
      ncol = 3
    ) +
    labs(
      x = "累计半径（km）",
      y = "上风向效应 − 下风向效应（mBq/m³）",
      title = paste0(
        branch_cn(
          branch_value
        ),
        "：方向性直接差值检验（",
        metric_cn(
          metric_value
        ),
        "）"
      ),
      subtitle = "大于0表示上风向效应高于下风向；误差线为95%置信区间"
    ) +
    theme_bw(
      base_size = 11
    ) +
    theme(
      panel.grid.minor =
        element_blank(),
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1
        )
    )

  p
}

for (
  branch_value in
  c(
    "filter_period",
    "quasi_daily"
  )
) {
  for (
    metric_value in
    c(
      "count",
      "heat"
    )
  ) {
    figure_counter <- figure_counter +
      1L

    p <- make_direction_plot(
      branch_value,
      metric_value
    )

    base <- file.path(
      FIG_DIR,
      sprintf(
        "%02d_方向性直接检验_%s_%s",
        figure_counter,
        branch_cn(
          branch_value
        ),
        metric_cn(
          metric_value
        )
      )
    )

    files <- save_plot_pair(
      p,
      base,
      12,
      8
    )

    plot_index <- plot_index +
      1L

    PLOT_MANIFEST[[plot_index]] <-
      data.table(
        figure_group =
          "方向性直接检验",
        branch =
          branch_value,
        metric =
          metric_value,
        base_file =
          base,
        png_file =
          files[1],
        pdf_file =
          files[2]
      )
  }
}

make_wind_plot <- function(
  branch_value,
  metric_value
) {
  D <- MAIN_CURVE[
    branch ==
      branch_value &
    metric ==
      metric_value &
    radius_km %in%
      REPRESENTATIVE_RADII
  ]

  D[
    ,
    `:=`(
      fuel_label =
        unname(
          FUEL_CN[fuel]
        ),
      radius_label =
        factor(
          paste0(
            radius_km,
            " km"
          ),
          levels =
            paste0(
              REPRESENTATIVE_RADII,
              " km"
            )
        )
    )
  ]

  p <- ggplot(
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
      alpha = 0.18
    ) +
    geom_line() +
    facet_grid(
      fuel_label ~ radius_label,
      scales = "free_y"
    ) +
    labs(
      x = "风速（m/s）",
      y = "暴露效应（mBq/m³）",
      title = paste0(
        branch_cn(
          branch_value
        ),
        "：连续风速效应曲线（",
        metric_cn(
          metric_value
        ),
        "）"
      ),
      subtitle = "20、50、100、200 km为预先指定的代表半径，仅用于图形展示，不用于选择最佳半径"
    ) +
    theme_bw(
      base_size = 9
    ) +
    theme(
      panel.grid.minor =
        element_blank(),
      strip.text =
        element_text(
          size = 8
        )
    )

  p
}

for (
  branch_value in
  c(
    "filter_period",
    "quasi_daily"
  )
) {
  for (
    metric_value in
    c(
      "count",
      "heat"
    )
  ) {
    figure_counter <- figure_counter +
      1L

    p <- make_wind_plot(
      branch_value,
      metric_value
    )

    base <- file.path(
      SUPPFIG_DIR,
      sprintf(
        "%02d_连续风速曲线_%s_%s",
        figure_counter,
        branch_cn(
          branch_value
        ),
        metric_cn(
          metric_value
        )
      )
    )

    files <- save_plot_pair(
      p,
      base,
      14,
      16
    )

    plot_index <- plot_index +
      1L

    PLOT_MANIFEST[[plot_index]] <-
      data.table(
        figure_group =
          "连续风速曲线",
        branch =
          branch_value,
        metric =
          metric_value,
        base_file =
          base,
        png_file =
          files[1],
        pdf_file =
          files[2]
      )
  }
}

PLOT_MANIFEST_DT <- rbindlist(
  PLOT_MANIFEST,
  use.names = TRUE
)

fwrite(
  PLOT_MANIFEST_DT,
  file.path(
    AUDIT_DIR,
    "01_图形输出清单.csv"
  )
)

# ============================================================
# 十五、自动生成可读结果摘要
# ============================================================

REPORT_LINES <- c(
  "# 07H11 结果审计与论文结果整理",
  "",
  paste0(
    "- 生成时间：",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  "",
  "## 1. 分析定位",
  "",
  "- 本步骤不重新拟合模型。",
  "- 所有统计结果均来自已经通过全部硬性检查的07H10。",
  "- 主要效应解释基于monitor×year固定效应，即使用2001–2024全部数据，但主要利用同一监测站同一年内部随时间发生的暴露变化。",
  "- 全部10个累计半径均保留，不根据P值选择最佳半径。",
  "- 风速作为连续效应修饰因素；1 m/s仅用于文献比较，正文主表同时保留样本中位风速和样本平均风速结果。",
  "- 距离加权Σ1/d和ΣHeat/d仅作为空间权重敏感性，不解释为实际浓度或剂量。",
  "",
  "## 2. 各燃料综合模式",
  ""
)

for (
  ii in
  seq_len(
    nrow(
      FUEL_SUMMARY
    )
  )
) {
  r <- FUEL_SUMMARY[ii]

  REPORT_LINES <- c(
    REPORT_LINES,
    paste0(
      "- ",
      branch_cn(
        r$branch
      ),
      " | ",
      r$fuel_label,
      " | ",
      r$metric_label,
      "：主模型",
      r$n_main_radii,
      "个可估计半径中",
      r$n_main_positive,
      "个正向、",
      r$n_main_negative,
      "个负向；",
      r$n_main_CI_above_zero,
      "个95%置信区间完全高于0、",
      r$n_main_CI_below_zero,
      "个完全低于0。严格时间控制：",
      r$strict_time_robustness,
      "；距离加权：",
      r$weighted_robustness,
      "；PM2.5敏感性：",
      r$PM_sensitivity_text,
      "；方向性：",
      r$directionality_text,
      "；连续风速：",
      r$wind_text,
      "。"
    )
  )
}

REPORT_LINES <- c(
  REPORT_LINES,
  "",
  "## 3. 时间尺度比较",
  ""
)

for (
  ii in
  seq_len(
    nrow(
      TIME_SCALE_SUMMARY
    )
  )
) {
  r <- TIME_SCALE_SUMMARY[ii]

  REPORT_LINES <- c(
    REPORT_LINES,
    paste0(
      "- ",
      r$fuel_label,
      " | ",
      r$metric_label,
      "：滤膜采样期与准逐日共有",
      r$n_comparable_radii,
      "个可比较半径；",
      r$n_same_direction,
      "个方向一致，",
      r$n_opposite_direction,
      "个方向相反。滤膜采样期为",
      r$filter_positive,
      "正/",
      r$filter_negative,
      "负；准逐日为",
      r$quasi_positive,
      "正/",
      r$quasi_negative,
      "负。"
    )
  )
}

REPORT_LINES <- c(
  REPORT_LINES,
  "",
  "## 4. 联合模型共线性",
  ""
)

for (
  ii in
  seq_len(
    nrow(
      VIF_GROUP
    )
  )
) {
  r <- VIF_GROUP[ii]

  REPORT_LINES <- c(
    REPORT_LINES,
    paste0(
      "- ",
      branch_cn(
        r$branch
      ),
      " | ",
      unname(
        FUEL_CN[
          r$fuel
        ]
      ),
      "：站点年份内Count–Heat相关系数绝对值最大值=",
      format(
        r$joint_within_corr_max_abs,
        digits = 4
      ),
      "；近似VIF中位数=",
      format(
        r$joint_VIF_median,
        digits = 4
      ),
      "；VIF≥5的半径数=",
      r$n_joint_VIF_ge5,
      "；VIF≥10的半径数=",
      r$n_joint_VIF_ge10,
      "；无穷VIF半径数=",
      r$n_joint_VIF_infinite,
      "。"
    )
  )
}

REPORT_LINES <- c(
  REPORT_LINES,
  "",
  "## 5. 结果解释原则",
  "",
  "- 正/负表示调整后的关联方向，不自动等同于因果效应。",
  "- 10个累计半径高度相关，不能把它们当成10次独立重复试验。",
  "- 上风向与下风向的正式比较使用同一模型中的直接差值及其95%置信区间。",
  "- 准逐日滞后分析保留全部预设滞后，不根据显著性挑选最佳滞后。",
  "- Count+Heat联合模型需与站点年份内相关系数和VIF共同解释。",
  "- 滤膜采样期与准逐日方向不一致时，应明确报告时间尺度差异，而不是只选择其中一套结果。",
  "- PM2.5调整前后方向明显变化时，应作为模型敏感性和潜在路径问题报告，而不是只选择更显著的一套模型。",
  ""
)

# ============================================================
# 十六、Hard checks
# ============================================================

EXPECTED_STATUS_ROWS <- data.table(
  analysis_family = c(
    "主分析",
    "严格时间固定效应",
    "上风向下风向联合模型",
    "距离加权",
    "Count+Heat联合模型",
    "PM2.5敏感性",
    "站点年份内外分解",
    "准逐日滞后"
  ),
  expected_rows = c(
    240L,
    240L,
    240L,
    480L,
    120L,
    240L,
    240L,
    720L
  )
)

OBSERVED_STATUS_ROWS <- data.table(
  analysis_family =
    EXPECTED_STATUS_ROWS$analysis_family,
  observed_rows = c(
    nrow(
      MAIN_STATUS
    ),
    nrow(
      STRICT_STATUS
    ),
    nrow(
      DIRECTION_STATUS
    ),
    nrow(
      WEIGHTED_STATUS
    ),
    nrow(
      JOINT_STATUS
    ),
    nrow(
      PM_STATUS
    ),
    nrow(
      DECOMP_STATUS
    ),
    nrow(
      LAG_STATUS
    )
  ),
  n_FAIL = c(
    MAIN_STATUS[
      status ==
        "FAIL",
      .N
    ],
    STRICT_STATUS[
      status ==
        "FAIL",
      .N
    ],
    DIRECTION_STATUS[
      status ==
        "FAIL",
      .N
    ],
    WEIGHTED_STATUS[
      status ==
        "FAIL",
      .N
    ],
    JOINT_STATUS[
      status ==
        "FAIL",
      .N
    ],
    PM_STATUS[
      status ==
        "FAIL",
      .N
    ],
    DECOMP_STATUS[
      status ==
        "FAIL",
      .N
    ],
    LAG_STATUS[
      status ==
        "FAIL",
      .N
    ]
  )
)

STATUS_AUDIT <- merge(
  EXPECTED_STATUS_ROWS,
  OBSERVED_STATUS_ROWS,
  by =
    "analysis_family"
)

STATUS_AUDIT[
  ,
  status :=
    fifelse(
      expected_rows ==
        observed_rows &
      n_FAIL ==
        0L,
      "PASS",
      "FAIL"
    )
]

fwrite(
  STATUS_AUDIT,
  file.path(
    AUDIT_DIR,
    "02_分析家族状态审计.csv"
  )
)

expected_main_effect_rows <-
  MAIN_STATUS[
    status ==
      "PASS",
    .N
  ] *
  5L

expected_strict_effect_rows <-
  STRICT_STATUS[
    status ==
      "PASS",
    .N
  ] *
  5L

expected_direction_effect_rows <-
  DIRECTION_STATUS[
    status ==
      "PASS",
    .N
  ] *
  15L

expected_weighted_effect_rows <-
  WEIGHTED_STATUS[
    status ==
      "PASS",
    .N
  ] *
  5L

expected_joint_effect_rows <-
  JOINT_STATUS[
    status ==
      "PASS",
    .N
  ] *
  10L

expected_PM_effect_rows <-
  PM_STATUS[
    status ==
      "PASS",
    .N
  ] *
  5L

expected_decomp_effect_rows <-
  DECOMP_STATUS[
    status ==
      "PASS",
    .N
  ] *
  10L

expected_lag_effect_rows <-
  LAG_STATUS[
    status ==
      "PASS",
    .N
  ] *
  5L

expected_main_curve_rows <-
  MAIN_STATUS[
    status ==
      "PASS",
    .N
  ] *
  201L

expected_strict_curve_rows <-
  STRICT_STATUS[
    status ==
      "PASS",
    .N
  ] *
  201L

N_PNG <- length(
  list.files(
    c(
      FIG_DIR,
      SUPPFIG_DIR
    ),
    pattern =
      "\\.png$",
    full.names =
      TRUE
  )
)

N_PDF <- length(
  list.files(
    c(
      FIG_DIR,
      SUPPFIG_DIR
    ),
    pattern =
      "\\.pdf$",
    full.names =
      TRUE
  )
)

HARD <- data.table(
  check = c(
    "H10_hard_checks_all_PASS",
    "H10_family_audit_all_PASS",
    "model_family_status_rows_complete",
    "all_model_FAIL_zero",
    "main_effect_rows_complete",
    "strict_effect_rows_complete",
    "direction_effect_rows_complete",
    "weighted_effect_rows_complete",
    "joint_effect_rows_complete",
    "PM_effect_rows_complete",
    "decomposition_effect_rows_complete",
    "lag_effect_rows_complete",
    "main_curve_201_points_per_PASS_model",
    "strict_curve_201_points_per_PASS_model",
    "main_median_table_nonempty",
    "fuel_summary_groups_24",
    "time_scale_summary_groups_12",
    "paper_tables_13",
    "plot_manifest_16_figures",
    "PNG_16",
    "PDF_16"
  ),
  observed = c(
    as.character(
      H10_HARD[
        status !=
          "PASS",
        .N
      ]
    ),
    as.character(
      H10_FAMILY_AUDIT[
        status !=
          "PASS",
        .N
      ]
    ),
    as.character(
      STATUS_AUDIT[
        status ==
          "PASS",
        .N
      ]
    ),
    as.character(
      ALL_STATUS[
        status ==
          "FAIL",
        .N
      ]
    ),
    as.character(
      nrow(
        MAIN_EFFECTS
      )
    ),
    as.character(
      nrow(
        STRICT_EFFECTS
      )
    ),
    as.character(
      nrow(
        DIRECTION_EFFECTS
      )
    ),
    as.character(
      nrow(
        WEIGHTED_EFFECTS
      )
    ),
    as.character(
      nrow(
        JOINT_EFFECTS
      )
    ),
    as.character(
      nrow(
        PM_EFFECTS
      )
    ),
    as.character(
      nrow(
        DECOMP_EFFECTS
      )
    ),
    as.character(
      nrow(
        LAG_EFFECTS
      )
    ),
    as.character(
      nrow(
        MAIN_CURVE
      )
    ),
    as.character(
      nrow(
        STRICT_CURVE
      )
    ),
    as.character(
      as.integer(
        nrow(
          MAIN_MEDIAN
        ) >
          0L
      )
    ),
    as.character(
      nrow(
        FUEL_SUMMARY
      )
    ),
    as.character(
      nrow(
        TIME_SCALE_SUMMARY
      )
    ),
    as.character(
      length(
        list.files(
          TABLE_DIR,
          pattern =
            "\\.csv$"
        )
      )
    ),
    as.character(
      nrow(
        PLOT_MANIFEST_DT
      )
    ),
    as.character(
      N_PNG
    ),
    as.character(
      N_PDF
    )
  ),
  expected = c(
    "0",
    "0",
    "8",
    "0",
    as.character(
      expected_main_effect_rows
    ),
    as.character(
      expected_strict_effect_rows
    ),
    as.character(
      expected_direction_effect_rows
    ),
    as.character(
      expected_weighted_effect_rows
    ),
    as.character(
      expected_joint_effect_rows
    ),
    as.character(
      expected_PM_effect_rows
    ),
    as.character(
      expected_decomp_effect_rows
    ),
    as.character(
      expected_lag_effect_rows
    ),
    as.character(
      expected_main_curve_rows
    ),
    as.character(
      expected_strict_curve_rows
    ),
    "1",
    "24",
    "12",
    "13",
    "16",
    "16",
    "16"
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
    "03_07H11_hard_checks.csv"
  )
)

HARD_FAIL <- HARD[
  status ==
    "FAIL",
  .N
]

REPORT_LINES <- c(
  REPORT_LINES,
  "",
  "## 6. 07H11运行状态",
  "",
  paste0(
    "- 论文表：",
    length(
      list.files(
        TABLE_DIR,
        pattern =
          "\\.csv$"
      )
    ),
    "个"
  ),
  paste0(
    "- 图形：",
    N_PNG,
    "个PNG + ",
    N_PDF,
    "个PDF"
  ),
  paste0(
    "- Hard FAIL：",
    HARD_FAIL
  ),
  "",
  if (
    HARD_FAIL ==
      0L
  ) {
    "**07H11 RESULT AUDIT AND PAPER PREPARATION: PASS**"
  } else {
    "**07H11 RESULT AUDIT AND PAPER PREPARATION: FAIL**"
  }
)

writeLines(
  REPORT_LINES,
  file.path(
    OUT,
    "07H11_结果审计与论文结果报告.md"
  ),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUT,
    "07H11_sessionInfo.txt"
  )
)

cat(
  "\n============================================================\n",
  "A. 07H11分析家族状态审计\n",
  "============================================================\n",
  sep = ""
)

print(
  STATUS_AUDIT,
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "B. 燃料层面综合结果\n",
  "============================================================\n",
  sep = ""
)

print(
  FUEL_SUMMARY[
    ,
    .(
      branch =
        branch_cn(branch),
      fuel_label,
      metric_label,
      n_main_radii,
      n_main_positive,
      n_main_negative,
      n_main_CI_above_zero,
      n_main_CI_below_zero,
      n_main_strict_same_direction,
      n_weighted_same_direction,
      n_PM_same_direction,
      n_up_minus_down_CI_above_zero,
      n_wind_always_positive,
      n_wind_always_negative,
      n_wind_direction_change
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "C. 滤膜采样期 vs 准逐日\n",
  "============================================================\n",
  sep = ""
)

print(
  TIME_SCALE_SUMMARY,
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "D. Count + Heat联合模型共线性\n",
  "============================================================\n",
  sep = ""
)

print(
  VIF_GROUP,
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "E. 站点年份内外方向相反情况\n",
  "============================================================\n",
  sep = ""
)

print(
  DECOMP_TABLE[
    structural_pattern ==
      "站点年份内关联与站点年份间关联方向相反",
    .(
      branch =
        branch_label,
      fuel =
        fuel_label,
      radius_km,
      metric =
        metric_label,
      old_estimate_wind1,
      within_estimate_wind1,
      between_estimate_wind1,
      old_vs_within,
      within_vs_between
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "F. 07H11 HARD CHECKS\n",
  "============================================================\n",
  sep = ""
)

print(
  HARD,
  nrows = Inf
)

cat(
  "\n============================================================\n"
)

if (
  HARD_FAIL ==
    0L
) {
  FREEZE <- c(
    "# 07H11 FINAL FREEZE",
    "",
    paste0(
      "- Freeze time: ",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    "- Upstream 07H10 hard checks: PASS",
    "- All 8 analysis families complete",
    "- Model FAIL count: 0",
    "- All 07H11 hard checks: PASS",
    "- No model was re-fitted in 07H11",
    "- All 10 predefined cumulative radii retained",
    "- No radius or lag was selected based on statistical significance",
    "",
    "**07H11 FINAL FREEZE: PASS**"
  )

  writeLines(
    FREEZE,
    file.path(
      OUT,
      "07H11_FINAL_FREEZE.md"
    ),
    useBytes = TRUE
  )

  cat(
    "07H11 RESULT AUDIT AND PAPER PREPARATION: PASS\n",
    "07H11 FINAL FREEZE: PASS\n"
  )

} else {

  cat(
    "07H11 RESULT AUDIT AND PAPER PREPARATION: FAIL\n"
  )
}

cat(
  "============================================================\n",
  "Hard FAIL = ",
  HARD_FAIL,
  "\n",
  sep = ""
)

if (
  HARD_FAIL >
    0L
) {
  stopf(
    "07H11存在硬性检查失败，禁止冻结。"
  )
}
