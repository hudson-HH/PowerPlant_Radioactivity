#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999,
  width = 220
)

required_packages <- c(
  "data.table",
  "ggplot2",
  "patchwork"
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
  library(ggplot2)
  library(patchwork)
})


# ============================================================
# 07H17
#
# 正文 Figure 6
# Internal validity and key robustness analyses
#
# A:
# Monitor × year FE
# vs
# Monitor × year FE + common-time FE
#
# B:
# Original random-intercept
# vs
# Within monitor-year
# vs
# Between monitor-years
#
# C:
# PM2.5 adjusted
# vs
# No PM2.5 adjustment
#
# 统一：
# - filter-period
# - 20 km anchor
# - 90° upwind actual operation
# - coal / gas / oil / biomass
# - Count / Heat
#
# Wind reference:
# A = sample median
# B = 1 m/s
#     原因：冻结的原随机截距模型是在1 m/s下保存，
#     为严格比较不能把新模型换成median而旧模型仍用1 m/s
# C = sample median
#
# 不重新拟合模型。
# 直接读取已经冻结并PASS的07H10结果。
# 并与当前07H16 Figure 4做交叉验证。
# ============================================================


# ============================================================
# 1. Paths
# ============================================================

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

HROOT <- file.path(
  ROOT,
  "发电站数据处理",
  "06_07H_燃料特异CAMPD逐日暴露与模型"
)

H10 <- file.path(
  HROOT,
  "07H10_燃料特异站点年份内综合正式分析"
)

H16 <- file.path(
  HROOT,
  "07H16_Figure4_Figure5正文与MonitorFE补充图"
)

OUT <- file.path(
  HROOT,
  "07H17_Figure6_内部有效性与关键稳健性"
)

FIG_DIR <- file.path(
  OUT,
  "figures"
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


MAIN_FILE <- file.path(
  H10,
  "01_主分析_监测站年份固定效应",
  "02_全部半径主模型代表风速效应.csv"
)

STRICT_FILE <- file.path(
  H10,
  "02_严格时间固定效应",
  "02_全部半径严格时间固定效应结果.csv"
)

PM_FILE <- file.path(
  H10,
  "07_细颗粒物敏感性",
  "02_不调整PM2.5同样本结果.csv"
)

DECOMP_FILE <- file.path(
  H10,
  "09_站点年份内外分解",
  "03_冻结原模型与站点年份内外方向解释.csv"
)

FIG4_SOURCE <- file.path(
  H16,
  "07_Figure4_source_monitorYearFE.csv"
)


# ============================================================
# 2. Parameters
# ============================================================

CORE_FUELS <- c(
  "coal",
  "gas",
  "oil",
  "biomass"
)

FUEL_LABEL <- c(
  coal = "Coal",
  gas = "Natural gas",
  oil = "Oil",
  biomass = "Biomass"
)

METRIC_LABEL <- c(
  count = "Active plant count\n(per additional active plant)",
  heat = "Actual daily Heat Input\n(per 10,000 MMBtu/day)"
)

ANCHOR_RADIUS <- 20L

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

check_columns <- function(
  D,
  required,
  object_name
) {

  miss <- setdiff(
    required,
    names(D)
  )

  if (
    length(miss) > 0L
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

prepare_labels <- function(D) {

  X <- copy(D)

  X[
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

  X[
    ,
    metric_label :=
      factor(
        METRIC_LABEL[
          metric
        ],
        levels =
          unname(
            METRIC_LABEL[
              c(
                "count",
                "heat"
              )
            ]
          )
      )
  ]

  X
}


# ============================================================
# 4. Input audit
# ============================================================

cat(
  "============================================================\n",
  "07H17 INPUT AUDIT\n",
  "============================================================\n",
  sep = ""
)

REQ_FILES <- c(
  MAIN_FILE,
  STRICT_FILE,
  PM_FILE,
  DECOMP_FILE,
  FIG4_SOURCE
)

MISSING_FILES <- REQ_FILES[
  !file.exists(
    REQ_FILES
  )
]

if (
  length(MISSING_FILES) > 0L
) {

  stopf(
    "缺少输入文件：\n",
    paste(
      MISSING_FILES,
      collapse = "\n"
    )
  )
}

cat(
  "07H10 main：PASS\n",
  "07H10 strict-time：PASS\n",
  "07H10 PM sensitivity：PASS\n",
  "07H10 within/between decomposition：PASS\n",
  "07H16 Figure 4 source：PASS\n",
  sep = ""
)


# ============================================================
# 5. Read frozen outputs
# ============================================================

MAIN <- fread(
  MAIN_FILE,
  showProgress = FALSE
)

STRICT <- fread(
  STRICT_FILE,
  showProgress = FALSE
)

PM <- fread(
  PM_FILE,
  showProgress = FALSE
)

DECOMP <- fread(
  DECOMP_FILE,
  showProgress = FALSE
)

FIG4 <- fread(
  FIG4_SOURCE,
  showProgress = FALSE
)


# ============================================================
# 6. Column audit
# ============================================================

check_columns(
  MAIN,
  c(
    "branch",
    "fuel",
    "metric",
    "radius_km",
    "wind_reference",
    "wind_ms",
    "estimate",
    "ci95_low",
    "ci95_high"
  ),
  "MAIN"
)

check_columns(
  STRICT,
  c(
    "branch",
    "fuel",
    "metric",
    "radius_km",
    "wind_reference",
    "wind_ms",
    "estimate",
    "ci95_low",
    "ci95_high"
  ),
  "STRICT"
)

check_columns(
  PM,
  c(
    "branch",
    "fuel",
    "metric",
    "radius_km",
    "wind_reference",
    "wind_ms",
    "estimate",
    "ci95_low",
    "ci95_high"
  ),
  "PM"
)

check_columns(
  DECOMP,
  c(
    "branch",
    "fuel",
    "metric",
    "radius_km",
    "old_estimate_wind1",
    "old_ci_low_wind1",
    "old_ci_high_wind1",
    "within_estimate_wind1",
    "within_ci_low_wind1",
    "within_ci_high_wind1",
    "between_estimate_wind1",
    "between_ci_low_wind1",
    "between_ci_high_wind1"
  ),
  "DECOMP"
)

check_columns(
  FIG4,
  c(
    "fe_structure",
    "fuel",
    "metric",
    "radius_km",
    "wind_reference",
    "wind_ms",
    "estimate",
    "ci95_low",
    "ci95_high"
  ),
  "FIG4"
)


# ============================================================
# 7. Current Figure 4 vs frozen H10 main cross-check
#
# 先确认Figure 6所调用的H10主模型
# 与现在正文Figure 4的20 km结果一致。
# ============================================================

H10_MAIN20 <- MAIN[
  branch ==
    "filter_period" &
  radius_km ==
    ANCHOR_RADIUS &
  fuel %in%
    CORE_FUELS &
  metric %in%
    c(
      "count",
      "heat"
    ) &
  wind_reference ==
    "样本中位数",
  .(
    fuel,
    metric,
    h10_wind_ms =
      wind_ms,
    h10_estimate =
      estimate,
    h10_ci_low =
      ci95_low,
    h10_ci_high =
      ci95_high
  )
]

H16_MAIN20 <- FIG4[
  fe_structure ==
    "monitor_year_FE" &
  radius_km ==
    ANCHOR_RADIUS &
  fuel %in%
    CORE_FUELS &
  metric %in%
    c(
      "count",
      "heat"
    ) &
  wind_reference ==
    "sample_median",
  .(
    fuel,
    metric,
    h16_wind_ms =
      wind_ms,
    h16_estimate =
      estimate,
    h16_ci_low =
      ci95_low,
    h16_ci_high =
      ci95_high
  )
]

CHECK_MAIN <- merge(
  H10_MAIN20,
  H16_MAIN20,
  by = c(
    "fuel",
    "metric"
  ),
  all = TRUE,
  sort = FALSE
)

CHECK_MAIN[
  ,
  `:=`(
    estimate_abs_diff =
      abs(
        h10_estimate -
        h16_estimate
      ),

    ci_low_abs_diff =
      abs(
        h10_ci_low -
        h16_ci_low
      ),

    ci_high_abs_diff =
      abs(
        h10_ci_high -
        h16_ci_high
      ),

    wind_abs_diff =
      abs(
        h10_wind_ms -
        h16_wind_ms
      )
  )
]

fwrite(
  CHECK_MAIN,
  file.path(
    OUT,
    "01_crosscheck_current_Figure4_vs_H10_main_20km.csv"
  )
)

if (
  nrow(CHECK_MAIN) !=
    8L
) {

  stopf(
    "Figure4 vs H10交叉验证不是8个组合；实际=",
    nrow(
      CHECK_MAIN
    )
  )
}

if (
  any(
    !is.finite(
      CHECK_MAIN$h10_estimate
    )
  ) ||
  any(
    !is.finite(
      CHECK_MAIN$h16_estimate
    )
  )
) {

  stopf(
    "Figure4 vs H10交叉验证存在缺失结果。"
  )
}

MAX_MAIN_DIFF <- max(
  c(
    CHECK_MAIN$estimate_abs_diff,
    CHECK_MAIN$ci_low_abs_diff,
    CHECK_MAIN$ci_high_abs_diff
  ),
  na.rm = TRUE
)

cat(
  "Figure4 vs H10 main max abs difference = ",
  format(
    MAX_MAIN_DIFF,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

if (
  !is.finite(
    MAX_MAIN_DIFF
  ) ||
  MAX_MAIN_DIFF >
    TOL
) {

  stopf(
    "当前Figure4与冻结H10主模型20 km不一致；禁止继续制作Figure6。"
  )
}

cat(
  "Current Figure4 ↔ H10 main 20 km：PASS\n"
)


# ============================================================
# 8. Figure 6A
# Main vs stricter common-time FE
#
# A/C使用样本中位风速
# ============================================================

A_MAIN <- MAIN[
  branch ==
    "filter_period" &
  radius_km ==
    ANCHOR_RADIUS &
  fuel %in%
    CORE_FUELS &
  metric %in%
    c(
      "count",
      "heat"
    ) &
  wind_reference ==
    "样本中位数",
  .(
    fuel,
    metric,
    model =
      "Monitor × year FE",
    wind_reference =
      "Sample median",
    wind_ms,
    estimate,
    ci_low =
      ci95_low,
    ci_high =
      ci95_high
  )
]

A_STRICT <- STRICT[
  branch ==
    "filter_period" &
  radius_km ==
    ANCHOR_RADIUS &
  fuel %in%
    CORE_FUELS &
  metric %in%
    c(
      "count",
      "heat"
    ) &
  wind_reference ==
    "样本中位数",
  .(
    fuel,
    metric,
    model =
      "Monitor × year + common-time FE",
    wind_reference =
      "Sample median",
    wind_ms,
    estimate,
    ci_low =
      ci95_low,
    ci_high =
      ci95_high
  )
]

D6A <- rbind(
  A_MAIN,
  A_STRICT,
  fill = TRUE
)

D6A <- prepare_labels(
  D6A
)

D6A[
  ,
  model :=
    factor(
      model,
      levels = c(
        "Monitor × year FE",
        "Monitor × year + common-time FE"
      )
    )
]


# ============================================================
# 9. Figure 6B
# Original random-intercept vs within/between
#
# 必须使用1 m/s：
# 原随机截距模型冻结结果就是1 m/s；
# 不能拿旧模型1 m/s和新模型median直接比较。
# ============================================================

B0 <- DECOMP[
  branch ==
    "filter_period" &
  radius_km ==
    ANCHOR_RADIUS &
  fuel %in%
    CORE_FUELS &
  metric %in%
    c(
      "count",
      "heat"
    )
]

D6B <- rbind(
  B0[
    ,
    .(
      fuel,
      metric,
      model =
        "Original random-intercept",
      wind_reference =
        "1 m/s",
      wind_ms =
        1,
      estimate =
        old_estimate_wind1,
      ci_low =
        old_ci_low_wind1,
      ci_high =
        old_ci_high_wind1
    )
  ],

  B0[
    ,
    .(
      fuel,
      metric,
      model =
        "Within monitor-year",
      wind_reference =
        "1 m/s",
      wind_ms =
        1,
      estimate =
        within_estimate_wind1,
      ci_low =
        within_ci_low_wind1,
      ci_high =
        within_ci_high_wind1
    )
  ],

  B0[
    ,
    .(
      fuel,
      metric,
      model =
        "Between monitor-years",
      wind_reference =
        "1 m/s",
      wind_ms =
        1,
      estimate =
        between_estimate_wind1,
      ci_low =
        between_ci_low_wind1,
      ci_high =
        between_ci_high_wind1
    )
  ],

  fill = TRUE
)

D6B <- prepare_labels(
  D6B
)

D6B[
  ,
  model :=
    factor(
      model,
      levels = c(
        "Original random-intercept",
        "Within monitor-year",
        "Between monitor-years"
      )
    )
]


# ============================================================
# 10. Figure 6C
# PM2.5 adjusted vs unadjusted
#
# H10的PM敏感性已经明确要求：
# no-PM模型与主模型保持完全相同行。
# 两套均使用sample median wind。
# ============================================================

C_MAIN <- MAIN[
  branch ==
    "filter_period" &
  radius_km ==
    ANCHOR_RADIUS &
  fuel %in%
    CORE_FUELS &
  metric %in%
    c(
      "count",
      "heat"
    ) &
  wind_reference ==
    "样本中位数",
  .(
    fuel,
    metric,
    model =
      "PM₂.₅ adjusted",
    wind_reference =
      "Sample median",
    wind_ms,
    estimate,
    ci_low =
      ci95_low,
    ci_high =
      ci95_high
  )
]

C_NOPM <- PM[
  branch ==
    "filter_period" &
  radius_km ==
    ANCHOR_RADIUS &
  fuel %in%
    CORE_FUELS &
  metric %in%
    c(
      "count",
      "heat"
    ) &
  wind_reference ==
    "样本中位数",
  .(
    fuel,
    metric,
    model =
      "No PM₂.₅ adjustment",
    wind_reference =
      "Sample median",
    wind_ms,
    estimate,
    ci_low =
      ci95_low,
    ci_high =
      ci95_high
  )
]

D6C <- rbind(
  C_MAIN,
  C_NOPM,
  fill = TRUE
)

D6C <- prepare_labels(
  D6C
)

D6C[
  ,
  model :=
    factor(
      model,
      levels = c(
        "PM₂.₅ adjusted",
        "No PM₂.₅ adjustment"
      )
    )
]


# ============================================================
# 11. Data integrity checks
# ============================================================

if (
  nrow(A_MAIN) !=
    8L ||
  nrow(A_STRICT) !=
    8L
) {

  stopf(
    "Figure6A预期主模型8行 + strict 8行；实际=",
    nrow(A_MAIN),
    " + ",
    nrow(A_STRICT)
  )
}

if (
  nrow(B0) !=
    8L
) {

  stopf(
    "Figure6B预期8个fuel×metric组合；实际=",
    nrow(B0)
  )
}

if (
  nrow(C_MAIN) !=
    8L ||
  nrow(C_NOPM) !=
    8L
) {

  stopf(
    "Figure6C预期PM adjusted 8行 + no-PM 8行；实际=",
    nrow(C_MAIN),
    " + ",
    nrow(C_NOPM)
  )
}

for (
  nm in c(
    "D6A",
    "D6B",
    "D6C"
  )
) {

  Z <- get(
    nm
  )

  if (
    any(
      !is.finite(
        Z$estimate
      )
    ) ||
    any(
      !is.finite(
        Z$ci_low
      )
    ) ||
    any(
      !is.finite(
        Z$ci_high
      )
    )
  ) {

    stopf(
      nm,
      "存在非有限estimate/CI。"
    )
  }
}


# ============================================================
# 12. Save plotting-source CSVs
# ============================================================

fwrite(
  D6A,
  file.path(
    OUT,
    "02_Figure6A_main_vs_strict_time_20km.csv"
  )
)

fwrite(
  D6B,
  file.path(
    OUT,
    "03_Figure6B_random_within_between_20km.csv"
  )
)

fwrite(
  D6C,
  file.path(
    OUT,
    "04_Figure6C_PM25_sensitivity_20km.csv"
  )
)


# ============================================================
# 13. Plot theme
# ============================================================

BASE_THEME <- theme_bw(
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
        face = "bold",
        size = 10
      ),

    legend.position =
      "bottom",

    legend.title =
      element_blank(),

    axis.text.x =
      element_text(
        angle = 25,
        hjust = 1
      ),

    plot.title =
      element_text(
        face = "bold",
        size = 11
      )
  )


# ============================================================
# 14. Generic panel maker
# ============================================================

make_panel <- function(
  D,
  title_text,
  y_title =
    "Adjusted gross-beta association (mBq/m³)",
  dodge_width =
    0.45
) {

  ggplot(
    D,
    aes(
      x = fuel_label,
      y = estimate,
      shape = model
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4
    ) +
    geom_errorbar(
      aes(
        ymin = ci_low,
        ymax = ci_high
      ),
      width = 0.12,
      linewidth = 0.45,
      position =
        position_dodge(
          width =
            dodge_width
        )
    ) +
    geom_point(
      size = 2.2,
      position =
        position_dodge(
          width =
            dodge_width
        )
    ) +
    facet_wrap(
      ~metric_label,
      nrow = 1,
      scales = "free_y"
    ) +
    labs(
      title =
        title_text,

      x = NULL,

      y =
        y_title,

      shape = NULL
    ) +
    BASE_THEME
}


# ============================================================
# 15. Panel A
# ============================================================

P6A <- make_panel(
  D6A,
  "A. Main model vs stricter common-time control"
)


# ============================================================
# 16. Panel B
# ============================================================

P6B <- make_panel(
  D6B,
  "B. Original random-intercept vs within- and between-monitor-year associations",
  dodge_width = 0.55
)


# ============================================================
# 17. Panel C
# ============================================================

P6C <- make_panel(
  D6C,
  "C. PM₂.₅-adjustment sensitivity"
)


# ============================================================
# 18. Combine Figure 6
# ============================================================

FIG6 <- (
  P6A /
  P6B /
  P6C
) +
  plot_annotation(
    title =
      "Figure 6. Internal validity and key robustness analyses",

    subtitle =
      paste0(
        "Filter-period analyses at the prespecified 20-km upwind anchor. ",
        "Panels A and C are evaluated at each model sample's median wind speed; ",
        "Panel B is evaluated at 1 m/s to preserve comparability with the frozen original random-intercept model."
      ),

    theme =
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 15
          ),

        plot.subtitle =
          element_text(
            size = 10
          )
      )
  )


# ============================================================
# 19. Save PNG + PDF
# ============================================================

PNG_FILE <- file.path(
  FIG_DIR,
  "Figure6_internal_validity_and_key_robustness.png"
)

PDF_FILE <- file.path(
  FIG_DIR,
  "Figure6_internal_validity_and_key_robustness.pdf"
)

ggsave(
  PNG_FILE,
  FIG6,
  width = 15,
  height = 15,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  PDF_FILE,
  FIG6,
  width = 15,
  height = 15,
  units = "in"
)


# ============================================================
# 20. Numerical summaries
# ============================================================

A_COMPARE <- dcast(
  D6A,
  fuel +
    metric ~
    model,
  value.var =
    "estimate"
)

C_COMPARE <- dcast(
  D6C,
  fuel +
    metric ~
    model,
  value.var =
    "estimate"
)

fwrite(
  A_COMPARE,
  file.path(
    OUT,
    "05_Figure6A_estimate_comparison.csv"
  )
)

fwrite(
  C_COMPARE,
  file.path(
    OUT,
    "06_Figure6C_estimate_comparison.csv"
  )
)


# ============================================================
# 21. Hard checks
# ============================================================

EXPECTED_FILES <- c(
  PNG_FILE,
  PDF_FILE,
  file.path(
    OUT,
    "02_Figure6A_main_vs_strict_time_20km.csv"
  ),
  file.path(
    OUT,
    "03_Figure6B_random_within_between_20km.csv"
  ),
  file.path(
    OUT,
    "04_Figure6C_PM25_sensitivity_20km.csv"
  )
)

HARD <- data.table(
  check = c(
    "current_Figure4_vs_H10_main_within_tolerance",
    "Figure6A_rows_16",
    "Figure6B_rows_24",
    "Figure6C_rows_16",
    "Figure6A_all_finite",
    "Figure6B_all_finite",
    "Figure6C_all_finite",
    "Figure6_output_files_5"
  ),

  observed = c(
    format(
      MAX_MAIN_DIFF,
      scientific = TRUE
    ),

    as.character(
      nrow(
        D6A
      )
    ),

    as.character(
      nrow(
        D6B
      )
    ),

    as.character(
      nrow(
        D6C
      )
    ),

    as.character(
      sum(
        is.finite(
          D6A$estimate
        ) &
        is.finite(
          D6A$ci_low
        ) &
        is.finite(
          D6A$ci_high
        )
      )
    ),

    as.character(
      sum(
        is.finite(
          D6B$estimate
        ) &
        is.finite(
          D6B$ci_low
        ) &
        is.finite(
          D6B$ci_high
        )
      )
    ),

    as.character(
      sum(
        is.finite(
          D6C$estimate
        ) &
        is.finite(
          D6C$ci_low
        ) &
        is.finite(
          D6C$ci_high
        )
      )
    ),

    as.character(
      sum(
        file.exists(
          EXPECTED_FILES
        )
      )
    )
  ),

  expected = c(
    paste0(
      "<=",
      TOL
    ),
    "16",
    "24",
    "16",
    "16",
    "24",
    "16",
    "5"
  )
)

HARD[
  ,
  status := c(
    ifelse(
      is.finite(
        MAX_MAIN_DIFF
      ) &&
      MAX_MAIN_DIFF <=
        TOL,
      "PASS",
      "FAIL"
    ),

    ifelse(
      nrow(
        D6A
      ) ==
        16L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      nrow(
        D6B
      ) ==
        24L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      nrow(
        D6C
      ) ==
        16L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      all(
        is.finite(
          D6A$estimate
        ) &
        is.finite(
          D6A$ci_low
        ) &
        is.finite(
          D6A$ci_high
        )
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      all(
        is.finite(
          D6B$estimate
        ) &
        is.finite(
          D6B$ci_low
        ) &
        is.finite(
          D6B$ci_high
        )
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      all(
        is.finite(
          D6C$estimate
        ) &
        is.finite(
          D6C$ci_low
        ) &
        is.finite(
          D6C$ci_high
        )
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      sum(
        file.exists(
          EXPECTED_FILES
        )
      ) ==
        5L,
      "PASS",
      "FAIL"
    )
  )
]

fwrite(
  HARD,
  file.path(
    OUT,
    "07_Figure6_hard_checks.csv"
  )
)


# ============================================================
# 22. Terminal summaries
# ============================================================

cat(
  "\n============================================================\n",
  "FIGURE 6A\n",
  "Main vs stricter common-time control\n",
  "20 km | sample median wind\n",
  "============================================================\n",
  sep = ""
)

print(
  D6A[
    ,
    .(
      fuel,
      metric,
      model,
      wind_ms,
      estimate,
      ci_low,
      ci_high
    )
  ][
    order(
      metric,
      fuel,
      model
    )
  ],
  nrows = Inf
)


cat(
  "\n============================================================\n",
  "FIGURE 6B\n",
  "Random-intercept vs within/between monitor-year\n",
  "20 km | 1 m/s\n",
  "============================================================\n",
  sep = ""
)

print(
  D6B[
    ,
    .(
      fuel,
      metric,
      model,
      estimate,
      ci_low,
      ci_high
    )
  ][
    order(
      metric,
      fuel,
      model
    )
  ],
  nrows = Inf
)


cat(
  "\n============================================================\n",
  "FIGURE 6C\n",
  "PM2.5 adjustment sensitivity\n",
  "20 km | sample median wind\n",
  "============================================================\n",
  sep = ""
)

print(
  D6C[
    ,
    .(
      fuel,
      metric,
      model,
      wind_ms,
      estimate,
      ci_low,
      ci_high
    )
  ][
    order(
      metric,
      fuel,
      model
    )
  ],
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
    "07H17 FIGURE 6 FINAL BUILD: PASS\n",
    "07H17 FIGURE 6 FINAL BUILD: FAIL\n"
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
    "07H17存在Hard FAIL，请检查07_Figure6_hard_checks.csv。"
  )
}

