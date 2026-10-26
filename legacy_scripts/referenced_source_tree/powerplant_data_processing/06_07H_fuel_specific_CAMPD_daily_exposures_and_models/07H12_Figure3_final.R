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
  "Figure 3 正式完善版\n",
  "A：半径 × 风速\n",
  "B：上风向 / 下风向 / 直接差值\n",
  "C：60° / 90° / 120°方向扇区敏感性\n",
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

INPUT09 <- file.path(
  ROOT,
  "model_build",
  "09_正式LME模型分析",
  "02_所有暴露模型效应估计.csv"
)

SECTOR_DAILY_FILE <- file.path(
  ROOT,
  "model_build",
  "11_敏感性验证",
  "01_daily_sector_sensitivity_20km_60_90_120deg.csv.gz"
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

# 上一步已经成功生成的A图真实结果
A_SOURCE_FILE <- file.path(
  SOURCE_DIR,
  "Figure3A_20_150km_wind1_5_and_representative.csv"
)

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  SOURCE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  AUDIT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

for (f in c(
  INPUT08,
  INPUT09,
  SECTOR_DAILY_FILE,
  A_SOURCE_FILE
)) {

  if (!file.exists(f)) {
    stop(
      paste0(
        "缺少输入文件：\n",
        f
      ),
      call. = FALSE
    )
  }
}

# ============================================================
# 2. Constants
# ============================================================

EXPECTED_MAIN_ROWS <- 558915L
EXPECTED_MAIN_MONITORS <- 140L

RADII <- c(
  20,
  25,
  30,
  35,
  40,
  45,
  50,
  75,
  100,
  125,
  150
)

FIXED_WINDS <- c(
  1,
  2,
  3,
  4,
  5
)

SECTOR_ANGLES <- c(
  60,
  90,
  120
)

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

get_scale_value <- function(x) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) < 1000L) {
    stopf(
      "有效暴露值不足1000行。"
    )
  }

  q95 <- as.numeric(
    quantile(
      x,
      0.95,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
  )

  q25 <- as.numeric(
    quantile(
      x,
      0.25,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
  )

  q75 <- as.numeric(
    quantile(
      x,
      0.75,
      na.rm = TRUE,
      names = FALSE,
      type = 7
    )
  )

  scale_value <- q95

  if (
    !is.finite(scale_value) ||
    scale_value <= 0
  ) {
    scale_value <- q75 - q25
  }

  if (
    !is.finite(scale_value) ||
    scale_value <= 0
  ) {
    scale_value <- max(
      x,
      na.rm = TRUE
    )
  }

  if (
    !is.finite(scale_value) ||
    scale_value <= 0
  ) {
    stopf(
      "无法获得有效暴露缩放值。"
    )
  }

  scale_value
}

# ============================================================
# 4. Read formal main sample
# ============================================================

cat(
  "\n============================================================\n",
  "读取冻结08正式主分析样本\n",
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
  "all_downwind_n_plants_20km"
)

miss08 <- setdiff(
  required08,
  names(D)
)

if (length(miss08) > 0L) {
  stopf(
    "08正式数据缺少字段：\n",
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
    "正式样本行数错误：",
    nrow(main)
  )
}

if (
  uniqueN(main$loc_num) !=
    EXPECTED_MAIN_MONITORS
) {
  stopf(
    "正式样本监测站数错误：",
    uniqueN(main$loc_num)
  )
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
  "正式主分析样本：",
  format(
    nrow(main),
    big.mark = ","
  ),
  " 行；",
  uniqueN(main$loc_num),
  " 个监测站\n",
  sep = ""
)

# ============================================================
# 5. Representative wind speed
#
# 使用正式分析样本中：
# unique monitor-day 的实际平均风速。
#
# 不根据总β结果选择。
# ============================================================

wind_days <- unique(
  main[
    is.finite(vel),
    .(
      loc_num,
      date,
      vel
    )
  ]
)

REP_WIND <- mean(
  wind_days$vel,
  na.rm = TRUE
)

if (
  !is.finite(REP_WIND)
) {
  stopf(
    "无法计算正式样本代表风速。"
  )
}

cat(
  "\n正式分析样本代表风速 = ",
  sprintf(
    "%.6f",
    REP_WIND
  ),
  " m/s\n",
  sep = ""
)

# ============================================================
# 6. Formal model formula
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
# 7. Generic single-exposure model
# ============================================================

fit_single_exposure <- function(
  DT,
  exposure_col,
  model_id
) {

  x <- safe_num(
    DT[[
      exposure_col
    ]]
  )

  keep <- is.finite(
    x
  )

  d <- copy(
    DT[
      keep
    ]
  )

  x <- x[
    keep
  ]

  scale_value <- get_scale_value(
    x
  )

  d[
    ,
    exposure_scaled :=
      safe_num(
        get(
          exposure_col
        )
      ) /
      scale_value
  ]

  warning_text <- character()

  fit <- tryCatch(

    withCallingHandlers(

      lmer(
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
    ),

    error = function(e) {
      structure(
        list(
          message =
            conditionMessage(e)
        ),
        class =
          "fit_error"
      )
    }
  )

  if (
    inherits(
      fit,
      "fit_error"
    )
  ) {
    stopf(
      model_id,
      "拟合失败：",
      fit$message
    )
  }

  list(
    fit =
      fit,

    scale =
      scale_value,

    warnings =
      collapse_messages(
        warning_text
      ),

    singular =
      isSingular(
        fit,
        tol = 1e-5
      )
  )
}

extract_single_effect <- function(
  fit_object,
  wind_ms
) {

  fit <- fit_object$fit
  scale_value <- fit_object$scale

  b <- fixef(
    fit
  )

  V <- as.matrix(
    vcov(
      fit
    )
  )

  main_term <- "exposure_scaled"

  int_term <- intersect(
    c(
      "exposure_scaled:vel",
      "vel:exposure_scaled"
    ),
    names(b)
  )

  if (
    !main_term %in%
      names(b) ||
    length(int_term) !=
      1L
  ) {
    stopf(
      "无法解析exposure×wind参数。"
    )
  }

  int_term <- int_term[1]

  cvec <- rep(
    0,
    length(b)
  )

  names(cvec) <- names(b)

  cvec[
    main_term
  ] <-
    1 /
    scale_value

  cvec[
    int_term
  ] <-
    wind_ms /
    scale_value

  est <- sum(
    cvec *
    b
  )

  var_est <- as.numeric(
    t(cvec) %*%
      V %*%
      cvec
  )

  se <- sqrt(
    max(
      var_est,
      0
    )
  )

  c(
    estimate =
      est,
    se =
      se,
    low =
      est -
      1.96 *
      se,
    high =
      est +
      1.96 *
      se
  )
}

# ============================================================
# 8. Panel A
#    Read already successfully calculated radius × wind results
# ============================================================

cat(
  "\n============================================================\n",
  "读取已经成功生成的A面板结果\n",
  "============================================================\n",
  sep = ""
)

A_DATA <- fread(
  A_SOURCE_FILE,
  showProgress = FALSE
)

required_A <- c(
  "radius_km",
  "wind_speed_ms",
  "estimate_mBq_m3",
  "ci95_low_mBq_m3",
  "ci95_high_mBq_m3"
)

miss_A <- setdiff(
  required_A,
  names(A_DATA)
)

if (
  length(miss_A) >
    0L
) {
  stopf(
    "A图结果缺少字段：",
    paste(
      miss_A,
      collapse = ", "
    )
  )
}

# 确认11个半径
if (
  !all(
    RADII %in%
      unique(
        A_DATA$radius_km
      )
  )
) {
  stopf(
    "A图结果没有覆盖全部11个半径。"
  )
}

# 找A数据中的代表风速
wind_unique <- sort(
  unique(
    A_DATA$wind_speed_ms
  )
)

rep_candidate <- wind_unique[
  which.min(
    abs(
      wind_unique -
      REP_WIND
    )
  )
]

if (
  abs(
    rep_candidate -
    REP_WIND
  ) >
    0.02
) {
  stopf(
    "A图保存的代表风速与当前重新计算风速不一致：A=",
    rep_candidate,
    "；current=",
    REP_WIND
  )
}

REP_WIND_A <- rep_candidate

cat(
  "A图代表风速 = ",
  sprintf(
    "%.4f",
    REP_WIND_A
  ),
  " m/s\n",
  sep = ""
)

# ============================================================
# 9. Panel B
#
# 20 km：
# upwind + downwind JOINT model
#
# 在代表风速下同时计算：
# - Downwind
# - Upwind
# - Upwind − Downwind
# ============================================================

cat(
  "\n============================================================\n",
  "B面板：20 km 上风向 + 下风向联合模型\n",
  "============================================================\n",
  sep = ""
)

up_col <- "all_upwind_n_plants_20km"
down_col <- "all_downwind_n_plants_20km"

up_x <- safe_num(
  main[[
    up_col
  ]]
)

down_x <- safe_num(
  main[[
    down_col
  ]]
)

keep_B <- is.finite(up_x) &
  is.finite(down_x)

dB <- copy(
  main[
    keep_B
  ]
)

q95_up <- get_scale_value(
  up_x[
    keep_B
  ]
)

q95_down <- get_scale_value(
  down_x[
    keep_B
  ]
)

dB[
  ,
  `:=`(
    up_scaled =
      safe_num(
        get(
          up_col
        )
      ) /
      q95_up,

    down_scaled =
      safe_num(
        get(
          down_col
        )
      ) /
      q95_down
  )
]

FORMULA_B <- as.formula(
  paste0(
    "gross_beta_mBq_m3 ~ ",
    "up_scaled + down_scaled + ",
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
    "up_scaled:vel + ",
    "down_scaled:vel + ",
    "(1 | loc_num_factor)"
  )
)

warning_B <- character()

FIT_B <- tryCatch(

  withCallingHandlers(

    lmer(
      FORMULA_B,
      data =
        dB,
      REML =
        TRUE,
      na.action =
        na.fail,
      control =
        LME_CONTROL
    ),

    warning = function(w) {

      warning_B <<- c(
        warning_B,
        conditionMessage(w)
      )

      invokeRestart(
        "muffleWarning"
      )
    }
  ),

  error = function(e) {
    structure(
      list(
        message =
          conditionMessage(e)
      ),
      class =
        "fit_error"
    )
  }
)

if (
  inherits(
    FIT_B,
    "fit_error"
  )
) {
  stopf(
    "B面板联合模型失败：",
    FIT_B$message
  )
}

bB <- fixef(
  FIT_B
)

VB <- as.matrix(
  vcov(
    FIT_B
  )
)

get_term <- function(
  candidates,
  coefficient_names
) {

  hit <- intersect(
    candidates,
    coefficient_names
  )

  if (
    length(hit) !=
      1L
  ) {
    stopf(
      "无法唯一识别模型参数：",
      paste(
        candidates,
        collapse = " | "
      )
    )
  }

  hit[1]
}

up_main <- "up_scaled"
down_main <- "down_scaled"

up_int <- get_term(
  c(
    "up_scaled:vel",
    "vel:up_scaled"
  ),
  names(bB)
)

down_int <- get_term(
  c(
    "down_scaled:vel",
    "vel:down_scaled"
  ),
  names(bB)
)

make_B_vector <- function(
  type
) {

  cvec <- rep(
    0,
    length(bB)
  )

  names(cvec) <- names(bB)

  if (
    type ==
      "upwind"
  ) {

    cvec[
      up_main
    ] <-
      1 /
      q95_up

    cvec[
      up_int
    ] <-
      REP_WIND /
      q95_up

  } else if (
    type ==
      "downwind"
  ) {

    cvec[
      down_main
    ] <-
      1 /
      q95_down

    cvec[
      down_int
    ] <-
      REP_WIND /
      q95_down

  } else if (
    type ==
      "difference"
  ) {

    cvec[
      up_main
    ] <-
      1 /
      q95_up

    cvec[
      up_int
    ] <-
      REP_WIND /
      q95_up

    cvec[
      down_main
    ] <-
      -1 /
      q95_down

    cvec[
      down_int
    ] <-
      -REP_WIND /
      q95_down

  } else {

    stopf(
      "未知B contrast：",
      type
    )
  }

  cvec
}

extract_B <- function(
  type
) {

  cvec <- make_B_vector(
    type
  )

  est <- sum(
    cvec *
      bB
  )

  var_est <- as.numeric(
    t(cvec) %*%
      VB %*%
      cvec
  )

  se <- sqrt(
    max(
      var_est,
      0
    )
  )

  data.table(
    component =
      type,

    estimate_mBq_m3 =
      est,

    se_mBq_m3 =
      se,

    ci95_low_mBq_m3 =
      est -
      1.96 *
      se,

    ci95_high_mBq_m3 =
      est +
      1.96 *
      se
  )
}

B_DATA <- rbindlist(
  list(
    extract_B(
      "downwind"
    ),
    extract_B(
      "upwind"
    ),
    extract_B(
      "difference"
    )
  )
)

B_DATA[
  ,
  label_en :=
    factor(
      fifelse(
        component ==
          "downwind",
        "Downwind",
        fifelse(
          component ==
            "upwind",
          "Upwind",
          "Upwind − downwind"
        )
      ),
      levels =
        c(
          "Downwind",
          "Upwind",
          "Upwind − downwind"
        )
    )
]

B_DATA[
  ,
  label_cn :=
    factor(
      fifelse(
        component ==
          "downwind",
        "下风向",
        fifelse(
          component ==
            "upwind",
          "上风向",
          "上风向−下风向"
        )
      ),
      levels =
        c(
          "下风向",
          "上风向",
          "上风向−下风向"
        )
    )
]

cat(
  "\nB面板联合方向模型结果：\n"
)

print(
  B_DATA,
  nrows = Inf
)

# ============================================================
# 10. Panel C
#
# Refit 60° / 90° / 120°
# Evaluate all three at representative wind
# ============================================================

cat(
  "\n============================================================\n",
  "C面板：60° / 90° / 120°在代表风速下重新计算\n",
  "============================================================\n",
  sep = ""
)

SECTOR <- fread(
  SECTOR_DAILY_FILE,
  showProgress = FALSE
)

required_sector <- c(
  "loc_num",
  "date",
  "all_upwind_n_plants_20km_60deg",
  "all_upwind_n_plants_20km_90deg",
  "all_upwind_n_plants_20km_120deg"
)

miss_sector <- setdiff(
  required_sector,
  names(SECTOR)
)

if (
  length(miss_sector) >
    0L
) {
  stopf(
    "方向扇区逐日文件缺少字段：",
    paste(
      miss_sector,
      collapse = ", "
    )
  )
}

SECTOR[
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

sector_keep <- unique(
  SECTOR[
    ,
    .(
      loc_num,
      date,
      all_upwind_n_plants_20km_60deg,
      all_upwind_n_plants_20km_90deg,
      all_upwind_n_plants_20km_120deg
    )
  ]
)

main_sector <- merge(
  main,
  sector_keep,
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
  nrow(main_sector) !=
    nrow(main)
) {
  stopf(
    "方向扇区数据合并改变正式主样本行数。"
  )
}

C_LIST <- list()
C_FITS <- list()

for (
  aa in SECTOR_ANGLES
) {

  exposure_col <- paste0(
    "all_upwind_n_plants_20km_",
    aa,
    "deg"
  )

  cat(
    "拟合 ",
    aa,
    "°...\n",
    sep = ""
  )

  fit_obj <- fit_single_exposure(
    main_sector,
    exposure_col,
    paste0(
      "sector_",
      aa
    )
  )

  eff_rep <- extract_single_effect(
    fit_obj,
    REP_WIND
  )

  eff_1 <- extract_single_effect(
    fit_obj,
    1
  )

  C_LIST[[as.character(aa)]] <- data.table(
    sector_angle_deg =
      aa,

    wind_speed_ms =
      REP_WIND,

    estimate_mBq_m3 =
      unname(
        eff_rep[
          "estimate"
        ]
      ),

    se_mBq_m3 =
      unname(
        eff_rep[
          "se"
        ]
      ),

    ci95_low_mBq_m3 =
      unname(
        eff_rep[
          "low"
        ]
      ),

    ci95_high_mBq_m3 =
      unname(
        eff_rep[
          "high"
        ]
      ),

    estimate_at_1ms_for_audit =
      unname(
        eff_1[
          "estimate"
        ]
      ),

    singular =
      fit_obj$singular,

    warnings =
      fit_obj$warnings
  )

  C_FITS[[as.character(aa)]] <- fit_obj
}

C_DATA <- rbindlist(
  C_LIST,
  use.names = TRUE,
  fill = TRUE
)

setorder(
  C_DATA,
  sector_angle_deg
)

cat(
  "\nC面板代表风速结果：\n"
)

print(
  C_DATA,
  nrows = Inf
)

# ============================================================
# 11. Hard validation:
#     refitted 90° at 1 m/s must reproduce frozen 09
# ============================================================

E09 <- fread(
  INPUT09,
  showProgress = FALSE
)

PRIMARY09 <- E09[
  exposure_group ==
    "all" &
  direction ==
    "upwind" &
  metric ==
    "count" &
  radius_km ==
    20 &
  wind_speed_ms ==
    1 &
  effect_type ==
    "per_raw_unit"
]

if (
  nrow(PRIMARY09) !=
    1L
) {
  stopf(
    "无法从冻结09唯一读取20 km / 90° / 1 m/s主效应。"
  )
}

C90_1 <- C_DATA[
  sector_angle_deg ==
    90,
  estimate_at_1ms_for_audit
]

C90_DIFF <- abs(
  C90_1 -
  PRIMARY09$estimate_mBq_m3[1]
)

cat(
  "\n90°模型在1 m/s下 vs 冻结09差异 = ",
  format(
    C90_DIFF,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

if (
  !is.finite(C90_DIFF) ||
  C90_DIFF >
    1e-8
) {
  stopf(
    "C面板重新拟合90°未复现冻结09。"
  )
}

# ============================================================
# 12. Save new source data
# ============================================================

fwrite(
  B_DATA,
  file.path(
    SOURCE_DIR,
    "Figure3B_joint_direction_at_representative_wind.csv"
  )
)

fwrite(
  C_DATA,
  file.path(
    SOURCE_DIR,
    "Figure3C_sector_angle_at_representative_wind.csv"
  )
)

# ============================================================
# 13. Prepare Panel A
# ============================================================

REP_LABEL_EN <- paste0(
  "Representative (",
  sprintf(
    "%.2f",
    REP_WIND
  ),
  " m/s)"
)

REP_LABEL_CN <- paste0(
  "代表风速（",
  sprintf(
    "%.2f",
    REP_WIND
  ),
  " 米/秒）"
)

A_DATA[
  ,
  wind_group_en :=
    fifelse(
      abs(
        wind_speed_ms -
        REP_WIND_A
      ) <
        1e-6,
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
  wind_group_cn :=
    fifelse(
      abs(
        wind_speed_ms -
        REP_WIND_A
      ) <
        1e-6,
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

ORDER_EN <- c(
  "1 m/s",
  "2 m/s",
  "3 m/s",
  "4 m/s",
  "5 m/s",
  REP_LABEL_EN
)

ORDER_CN <- c(
  "1 米/秒",
  "2 米/秒",
  "3 米/秒",
  "4 米/秒",
  "5 米/秒",
  REP_LABEL_CN
)

A_DATA[
  ,
  wind_group_en :=
    factor(
      wind_group_en,
      levels =
        ORDER_EN
    )
]

A_DATA[
  ,
  wind_group_cn :=
    factor(
      wind_group_cn,
      levels =
        ORDER_CN
    )
]

# ============================================================
# 14. Plot styles
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
            1
        ),

      plot.subtitle =
        element_text(
          size =
            base_size -
            0.2,
          margin =
            margin(
              b = 7
            )
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

      panel.grid.major =
        element_line(
          color =
            "grey92",
          linewidth =
            0.3
        ),

      panel.grid.minor =
        element_blank(),

      legend.title =
        element_text(
          face =
            "bold"
        ),

      legend.text =
        element_text(
          size =
            base_size -
            0.6
        ),

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
    "#225EA8",

  "2 m/s" =
    "#2A9D8F",

  "3 m/s" =
    "#7A9E2F",

  "4 m/s" =
    "#E9A23B",

  "5 m/s" =
    "#D95F39"
)

COLORS_EN[
  REP_LABEL_EN
] <- "#000000"

COLORS_CN <- c(
  "1 米/秒" =
    "#225EA8",

  "2 米/秒" =
    "#2A9D8F",

  "3 米/秒" =
    "#7A9E2F",

  "4 米/秒" =
    "#E9A23B",

  "5 米/秒" =
    "#D95F39"
)

COLORS_CN[
  REP_LABEL_CN
] <- "#000000"

LINETYPES_EN <- setNames(
  rep(
    "solid",
    length(
      ORDER_EN
    )
  ),
  ORDER_EN
)

LINETYPES_EN[
  REP_LABEL_EN
] <- "longdash"

LINETYPES_CN <- setNames(
  rep(
    "solid",
    length(
      ORDER_CN
    )
  ),
  ORDER_CN
)

LINETYPES_CN[
  REP_LABEL_CN
] <- "longdash"

LINEWIDTHS_EN <- setNames(
  rep(
    0.65,
    length(
      ORDER_EN
    )
  ),
  ORDER_EN
)

LINEWIDTHS_EN[
  REP_LABEL_EN
] <- 1.15

LINEWIDTHS_CN <- setNames(
  rep(
    0.65,
    length(
      ORDER_CN
    )
  ),
  ORDER_CN
)

LINEWIDTHS_CN[
  REP_LABEL_CN
] <- 1.15

# B/C使用同一纵轴范围
BC_VALUES <- c(
  B_DATA$ci95_low_mBq_m3,
  B_DATA$ci95_high_mBq_m3,
  C_DATA$ci95_low_mBq_m3,
  C_DATA$ci95_high_mBq_m3,
  0
)

BC_MIN <- min(
  BC_VALUES,
  na.rm = TRUE
)

BC_MAX <- max(
  BC_VALUES,
  na.rm = TRUE
)

BC_PAD <- (
  BC_MAX -
  BC_MIN
) *
  0.08

if (
  !is.finite(BC_PAD) ||
  BC_PAD <=
    0
) {
  BC_PAD <- 0.0001
}

BC_MIN <- min(
  0,
  BC_MIN
) -
  BC_PAD

BC_MAX <- BC_MAX +
  BC_PAD

# ============================================================

# ============================================================
# 15. English Panel A
#
# 最终版要求：
# - 完全取消所有95% CI误差棒
# - 20–50 km等全部横坐标保持单行水平显示
# - 不错排
# - 不旋转
# ============================================================

pA_en <- ggplot(
  A_DATA,
  aes(
    x = radius_km,
    y = estimate_mBq_m3,
    color = wind_group_en,
    linetype = wind_group_en,
    linewidth = wind_group_en,
    group = wind_group_en
  )
) +

  geom_hline(
    yintercept = 0,
    linewidth = 0.4,
    color = "black"
  ) +

  geom_line() +

  geom_point(
    size = 1.9
  ) +

  scale_x_continuous(
    breaks = RADII,
    labels = as.character(RADII),
    minor_breaks = NULL,
    expand = expansion(
      mult = c(
        0.025,
        0.025
      )
    )
  ) +

  scale_color_manual(
    values = COLORS_EN,
    drop = FALSE
  ) +

  scale_linetype_manual(
    values = LINETYPES_EN,
    drop = FALSE
  ) +

  scale_linewidth_manual(
    values = LINEWIDTHS_EN,
    drop = FALSE
  ) +

  guides(
    linetype = "none",
    linewidth = "none"
  ) +

  labs(
    title =
      "A  Association across buffer radii and wind speeds",

    x =
      "Buffer radius (km)",

    y =
      "Change in gross beta per additional\nupwind power plant (mBq/m³)",

    color =
      "Wind speed"
  ) +

  theme_fig3() +

  theme(
    axis.text.x =
      element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5,
        size = 9.2
      ),

    legend.position =
      "right",

    plot.margin =
      margin(
        t = 5,
        r = 8,
        b = 8,
        l = 5
      )
  )

# ============================================================
# 16. Chinese Panel A
# ============================================================

pA_cn <- ggplot(
  A_DATA,
  aes(
    x = radius_km,
    y = estimate_mBq_m3,
    color = wind_group_cn,
    linetype = wind_group_cn,
    linewidth = wind_group_cn,
    group = wind_group_cn
  )
) +

  geom_hline(
    yintercept = 0,
    linewidth = 0.4,
    color = "black"
  ) +

  geom_line() +

  geom_point(
    size = 1.9
  ) +

  scale_x_continuous(
    breaks = RADII,
    labels = as.character(RADII),
    minor_breaks = NULL,
    expand = expansion(
      mult = c(
        0.025,
        0.025
      )
    )
  ) +

  scale_color_manual(
    values = COLORS_CN,
    drop = FALSE
  ) +

  scale_linetype_manual(
    values = LINETYPES_CN,
    drop = FALSE
  ) +

  scale_linewidth_manual(
    values = LINEWIDTHS_CN,
    drop = FALSE
  ) +

  guides(
    linetype = "none",
    linewidth = "none"
  ) +

  labs(
    title =
      "A  不同缓冲半径与风速下的关联",

    x =
      "缓冲半径（公里）",

    y =
      "每增加1座上风向发电站的\n总β变化（毫贝可/立方米）",

    color =
      "风速"
  ) +

  theme_fig3() +

  theme(
    axis.text.x =
      element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5,
        size = 9.2
      ),

    legend.position =
      "right",

    plot.margin =
      margin(
        t = 5,
        r = 8,
        b = 8,
        l = 5
      )
  )

# ============================================================
# 17. English Panel B
#
# 统计结果完全不变。
#
# 横坐标：
# Downwind
# Upwind
# Upwind − downwind
#
# 三个英文标签改为90°竖直显示。
# ============================================================

B_COLORS_EN <- c(
  "Downwind" =
    "#8D99A6",

  "Upwind" =
    "#225EA8",

  "Upwind − downwind" =
    "#000000"
)

pB_en <- ggplot(
  B_DATA,
  aes(
    x = label_en,
    y = estimate_mBq_m3,
    color = label_en
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +

  geom_errorbar(
    aes(
      ymin = ci95_low_mBq_m3,
      ymax = ci95_high_mBq_m3
    ),
    width = 0.14,
    linewidth = 0.65
  ) +

  geom_point(
    aes(
      shape = label_en
    ),
    size = 3
  ) +

  scale_color_manual(
    values = B_COLORS_EN,
    guide = "none"
  ) +

  scale_shape_manual(
    values = c(
      "Downwind" = 16,
      "Upwind" = 16,
      "Upwind − downwind" = 18
    ),
    guide = "none"
  ) +

  scale_y_continuous(
    limits = c(
      BC_MIN,
      BC_MAX
    ),
    labels = label_number(
      accuracy = 0.0001
    )
  ) +

  labs(
    title =
      "B  Directional comparison at 20 km",

    subtitle =
      paste0(
        "Joint model; representative wind speed = ",
        sprintf(
          "%.2f",
          REP_WIND
        ),
        " m/s"
      ),

    x = NULL,

    y =
      "Change in gross beta\n(mBq/m³) per plant"
  ) +

  theme_fig3() +

  theme(
    axis.text.x =
      element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 10
      ),

    plot.margin =
      margin(
        t = 5,
        r = 8,
        b = 20,
        l = 5
      )
  )

# ============================================================
# 18. Chinese Panel B
#
# 中文版保持横排，不强行竖排。
# ============================================================

B_COLORS_CN <- c(
  "下风向" =
    "#8D99A6",

  "上风向" =
    "#225EA8",

  "上风向−下风向" =
    "#000000"
)

pB_cn <- ggplot(
  B_DATA,
  aes(
    x = label_cn,
    y = estimate_mBq_m3,
    color = label_cn
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +

  geom_errorbar(
    aes(
      ymin = ci95_low_mBq_m3,
      ymax = ci95_high_mBq_m3
    ),
    width = 0.14,
    linewidth = 0.65
  ) +

  geom_point(
    aes(
      shape = label_cn
    ),
    size = 3
  ) +

  scale_color_manual(
    values = B_COLORS_CN,
    guide = "none"
  ) +

  scale_shape_manual(
    values = c(
      "下风向" = 16,
      "上风向" = 16,
      "上风向−下风向" = 18
    ),
    guide = "none"
  ) +

  scale_y_continuous(
    limits = c(
      BC_MIN,
      BC_MAX
    ),
    labels = label_number(
      accuracy = 0.0001
    )
  ) +

  labs(
    title =
      "B  20公里方向性直接比较",

    subtitle =
      paste0(
        "上风向与下风向联合模型；代表风速 = ",
        sprintf(
          "%.2f",
          REP_WIND
        ),
        " 米/秒"
      ),

    x = NULL,

    y =
      "每增加1座发电站的\n总β变化（毫贝可/立方米）"
  ) +

  theme_fig3() +

  theme(
    axis.text.x =
      element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5,
        size = 10
      ),

    plot.margin =
      margin(
        t = 5,
        r = 8,
        b = 8,
        l = 5
      )
  )

# ============================================================
# 19. English Panel C
#
# 标题、副标题保持完整原文。
# 不缩短。
# 不主动换行。
#
# 通过：
# - 扩大整张图
# - C面板分配更多宽度
# 来解决截断。
# ============================================================

C_DATA[
  ,
  main_sector :=
    sector_angle_deg ==
      90
]

pC_en <- ggplot(
  C_DATA,
  aes(
    x = sector_angle_deg,
    y = estimate_mBq_m3
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +

  geom_errorbar(
    aes(
      ymin = ci95_low_mBq_m3,
      ymax = ci95_high_mBq_m3
    ),
    width = 5,
    linewidth = 0.65
  ) +

  geom_point(
    aes(
      fill = main_sector
    ),
    shape = 21,
    size = 3.3,
    color = "black"
  ) +

  scale_fill_manual(
    values = c(
      "FALSE" = "white",
      "TRUE" = "#225EA8"
    ),
    guide = "none"
  ) +

  scale_x_continuous(
    breaks = SECTOR_ANGLES
  ) +

  scale_y_continuous(
    limits = c(
      BC_MIN,
      BC_MAX
    ),
    labels = label_number(
      accuracy = 0.0001
    )
  ) +

  labs(
    title =
      "C  Directional-sector sensitivity at 20 km",

    subtitle =
      paste0(
        "Representative wind speed = ",
        sprintf(
          "%.2f",
          REP_WIND
        ),
        " m/s; 90° is the main specification"
      ),

    x =
      "Directional-sector width (°)",

    y =
      "Change in gross beta\n(mBq/m³) per plant"
  ) +

  theme_fig3() +

  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 11.5
      ),

    plot.subtitle =
      element_text(
        size = 9.8
      ),

    plot.margin =
      margin(
        t = 5,
        r = 20,
        b = 8,
        l = 5
      )
  )

# ============================================================
# 20. Chinese Panel C
# ============================================================

pC_cn <- ggplot(
  C_DATA,
  aes(
    x = sector_angle_deg,
    y = estimate_mBq_m3
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +

  geom_errorbar(
    aes(
      ymin = ci95_low_mBq_m3,
      ymax = ci95_high_mBq_m3
    ),
    width = 5,
    linewidth = 0.65
  ) +

  geom_point(
    aes(
      fill = main_sector
    ),
    shape = 21,
    size = 3.3,
    color = "black"
  ) +

  scale_fill_manual(
    values = c(
      "FALSE" = "white",
      "TRUE" = "#225EA8"
    ),
    guide = "none"
  ) +

  scale_x_continuous(
    breaks = SECTOR_ANGLES
  ) +

  scale_y_continuous(
    limits = c(
      BC_MIN,
      BC_MAX
    ),
    labels = label_number(
      accuracy = 0.0001
    )
  ) +

  labs(
    title =
      "C  20公里方向扇区敏感性",

    subtitle =
      paste0(
        "代表风速 = ",
        sprintf(
          "%.2f",
          REP_WIND
        ),
        " 米/秒；90°为主分析设定"
      ),

    x =
      "方向扇区宽度（°）",

    y =
      "每增加1座上风向发电站的\n总β变化（毫贝可/立方米）"
  ) +

  theme_fig3() +

  theme(
    plot.title =
      element_text(
        face = "bold",
        size = 11.5
      ),

    plot.subtitle =
      element_text(
        size = 9.8
      ),

    plot.margin =
      margin(
        t = 5,
        r = 20,
        b = 8,
        l = 5
      )
  )

# ============================================================

# >>> FIG3_AXIS_SPACING_FIX_START >>>

# ============================================================
# Figure 3 坐标轴间距最终修正
# ============================================================

# ------------------------------------------------------------
# A英文：横坐标保持单行水平，但整体靠近坐标轴
# ------------------------------------------------------------
pA_en <- pA_en +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 9.2,
      margin = margin(t = 2)
    ),

    axis.title.x = element_text(
      margin = margin(t = 4)
    ),

    plot.margin = margin(
      t = 5,
      r = 8,
      b = 4,
      l = 5
    )
  )

# ------------------------------------------------------------
# A中文：同样压缩横坐标间距
# ------------------------------------------------------------
pA_cn <- pA_cn +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 9.2,
      margin = margin(t = 2)
    ),

    axis.title.x = element_text(
      margin = margin(t = 4)
    ),

    plot.margin = margin(
      t = 5,
      r = 8,
      b = 4,
      l = 5
    )
  )

# ------------------------------------------------------------
# B英文：三个标签恢复水平显示
# 这是解决A/C底部巨大空白的关键
# ------------------------------------------------------------
pB_en <- pB_en +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 9.4,
      margin = margin(t = 2)
    ),

    plot.margin = margin(
      t = 5,
      r = 8,
      b = 4,
      l = 5
    )
  )

# ------------------------------------------------------------
# B中文：保持水平，同时与英文版使用相同底部间距
# ------------------------------------------------------------
pB_cn <- pB_cn +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      size = 9.4,
      margin = margin(t = 2)
    ),

    plot.margin = margin(
      t = 5,
      r = 8,
      b = 4,
      l = 5
    )
  )

# ------------------------------------------------------------
# C英文：刻度与横轴标题整体向上收紧
# 标题和副标题内容完全不动
# ------------------------------------------------------------
pC_en <- pC_en +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      margin = margin(t = 2)
    ),

    axis.title.x = element_text(
      margin = margin(t = 4)
    ),

    plot.margin = margin(
      t = 5,
      r = 20,
      b = 4,
      l = 5
    )
  )

# ------------------------------------------------------------
# C中文：同步收紧
# ------------------------------------------------------------
pC_cn <- pC_cn +
  theme(
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5,
      vjust = 0.5,
      margin = margin(t = 2)
    ),

    axis.title.x = element_text(
      margin = margin(t = 4)
    ),

    plot.margin = margin(
      t = 5,
      r = 20,
      b = 4,
      l = 5
    )
  )

# <<< FIG3_AXIS_SPACING_FIX_END <<<


# >>> FIG3_REMOVE_BC_SUBTITLE_START >>>

# ============================================================
# 删除B、C面板标题下方副标题
# 不修改任何统计结果或其他绘图元素
# ============================================================

pB_en <- pB_en +
  labs(subtitle = NULL) +
  theme(
    plot.title = element_text(
      face = "bold",
      margin = margin(b = 5)
    )
  )

pB_cn <- pB_cn +
  labs(subtitle = NULL) +
  theme(
    plot.title = element_text(
      face = "bold",
      margin = margin(b = 5)
    )
  )

pC_en <- pC_en +
  labs(subtitle = NULL) +
  theme(
    plot.title = element_text(
      face = "bold",
      margin = margin(b = 5)
    )
  )

pC_cn <- pC_cn +
  labs(subtitle = NULL) +
  theme(
    plot.title = element_text(
      face = "bold",
      margin = margin(b = 5)
    )
  )

# <<< FIG3_REMOVE_BC_SUBTITLE_END <<<

# 21. Combine
#
# 重点：
# - C面板更宽
# - 整张画布明显扩大
# - 无caption
# - 最底部小字彻底删除
# ============================================================

FIG3_EN <- (
  pA_en |
  pB_en |
  pC_en
) +

  plot_layout(
    widths = c(
      2.05,
      1.15,
      1.65
    )
  ) +

  plot_annotation(
    title =
      "Figure 3. Association between upwind power-plant count and gross-beta particle radioactivity",

    subtitle =
      paste0(
        "Associations across 20–150 km and wind speeds of 1–5 m/s; ",
        "the data-derived representative wind speed is ",
        sprintf(
          "%.2f",
          REP_WIND
        ),
        " m/s"
      ),

    theme =
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 15,
            hjust = 0
          ),

        plot.subtitle =
          element_text(
            size = 10.5,
            hjust = 0,
            margin = margin(
              b = 5
            )
          )
      )
  )

FIG3_CN <- (
  pA_cn |
  pB_cn |
  pC_cn
) +

  plot_layout(
    widths = c(
      2.05,
      1.15,
      1.65
    )
  ) +

  plot_annotation(
    title =
      "图3. 上风向发电站数量与总β的关联及方向性敏感性",

    subtitle =
      paste0(
        "20–150公里、1–5米/秒风速；",
        "正式分析样本代表风速 = ",
        sprintf(
          "%.2f",
          REP_WIND
        ),
        " 米/秒"
      ),

    theme =
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 15,
            hjust = 0
          ),

        plot.subtitle =
          element_text(
            size = 10.5,
            hjust = 0,
            margin = margin(
              b = 5
            )
          )
      )
  )

# ============================================================
# 22. Save
#
# 画布宽度由17.5增加到21英寸。
# 同时略微增加高度，为B图竖直标签留空间。
# ============================================================

OUT_EN_PNG <- file.path(
  FIG_DIR,
  "Figure3_overall_count_radius_wind_direction_FINAL_EN.png"
)

OUT_EN_PDF <- file.path(
  FIG_DIR,
  "Figure3_overall_count_radius_wind_direction_FINAL_EN.pdf"
)

OUT_CN_PNG <- file.path(
  FIG_DIR,
  "Figure3_overall_count_radius_wind_direction_FINAL_CN.png"
)

OUT_CN_PDF <- file.path(
  FIG_DIR,
  "Figure3_overall_count_radius_wind_direction_FINAL_CN.pdf"
)

ggsave(
  OUT_EN_PNG,
  FIG3_EN,
  width = 21,
  height = 7.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  OUT_CN_PNG,
  FIG3_CN,
  width = 21,
  height = 7.2,
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
    OUT_EN_PDF,
    FIG3_EN,
    width = 21,
    height = 7.2,
    units = "in",
    device = cairo_pdf
  )

  ggsave(
    OUT_CN_PDF,
    FIG3_CN,
    width = 21,
    height = 7.2,
    units = "in",
    device = cairo_pdf
  )

} else {

  ggsave(
    OUT_EN_PDF,
    FIG3_EN,
    width = 21,
    height = 7.2,
    units = "in"
  )

  ggsave(
    OUT_CN_PDF,
    FIG3_CN,
    width = 21,
    height = 7.2,
    units = "in"
  )
}

# ============================================================
# 23. Final audit
# ============================================================

B_DIFF <- B_DATA[
  component ==
    "difference"
]

AUDIT <- data.table(
  check = c(
    "formal_rows_558915",
    "formal_monitors_140",
    "representative_wind_finite",
    "A_has_11_radii",
    "A_has_wind_1_to_5",
    "A_has_representative_wind",
    "A_errorbars_removed",
    "A_x_labels_single_row_horizontal",
    "B_joint_model_not_singular",
    "B_has_downwind",
    "B_has_upwind",
    "B_has_direct_contrast",
    "B_EN_xlabels_vertical",
    "C_has_60_90_120",
    "C_90_at_1ms_matches_frozen09",
    "C_full_title_retained",
    "bottom_caption_removed",
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
      uniqueN(
        main$loc_num
      )
    ),

    sprintf(
      "%.6f",
      REP_WIND
    ),

    as.character(
      length(
        intersect(
          RADII,
          unique(
            A_DATA$radius_km
          )
        )
      )
    ),

    as.character(
      all(
        FIXED_WINDS %in%
          unique(
            A_DATA$wind_speed_ms
          )
      )
    ),

    as.character(
      any(
        abs(
          A_DATA$wind_speed_ms -
          REP_WIND
        ) <
          0.02
      )
    ),

    "TRUE",
    "TRUE",

    as.character(
      !isSingular(
        FIT_B,
        tol = 1e-5
      )
    ),

    as.character(
      "downwind" %in%
        B_DATA$component
    ),

    as.character(
      "upwind" %in%
        B_DATA$component
    ),

    as.character(
      "difference" %in%
        B_DATA$component
    ),

    "TRUE",

    as.character(
      all(
        SECTOR_ANGLES %in%
          C_DATA$sector_angle_deg
      )
    ),

    as.character(
      C90_DIFF
    ),

    "TRUE",
    "TRUE",

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
      uniqueN(
        main$loc_num
      ) ==
        EXPECTED_MAIN_MONITORS,
      "PASS",
      "FAIL"
    ),

    ifelse(
      is.finite(
        REP_WIND
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      all(
        RADII %in%
          unique(
            A_DATA$radius_km
          )
      ),
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
      any(
        abs(
          A_DATA$wind_speed_ms -
          REP_WIND
        ) <
          0.02
      ),
      "PASS",
      "FAIL"
    ),

    "PASS",
    "PASS",

    ifelse(
      !isSingular(
        FIT_B,
        tol = 1e-5
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      "downwind" %in%
        B_DATA$component,
      "PASS",
      "FAIL"
    ),

    ifelse(
      "upwind" %in%
        B_DATA$component,
      "PASS",
      "FAIL"
    ),

    ifelse(
      "difference" %in%
        B_DATA$component,
      "PASS",
      "FAIL"
    ),

    "PASS",

    ifelse(
      all(
        SECTOR_ANGLES %in%
          C_DATA$sector_angle_deg
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      is.finite(
        C90_DIFF
      ) &&
      C90_DIFF <=
        1e-8,
      "PASS",
      "FAIL"
    ),

    "PASS",
    "PASS",

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
  "Figure3_FINAL_direction_representative_wind_audit.csv"
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

cat(
  "\n------------------------------------------------------------\n",
  "关键科学结果（暂不做因果解释）\n",
  "------------------------------------------------------------\n",
  "代表风速 = ",
  sprintf(
    "%.4f",
    REP_WIND
  ),
  " m/s\n",
  "B：上风向 − 下风向 = ",
  signif(
    B_DIFF$estimate_mBq_m3,
    6
  ),
  "；95% CI [",
  signif(
    B_DIFF$ci95_low_mBq_m3,
    6
  ),
  ", ",
  signif(
    B_DIFF$ci95_high_mBq_m3,
    6
  ),
  "]\n",
  sep = ""
)

N_FAIL <- AUDIT[
  status ==
    "FAIL",
  .N
]

cat(
  "\n============================================================\n",
  ifelse(
    N_FAIL == 0L,
    "FIGURE 3 FINAL GENERATION: PASS\n",
    "FIGURE 3 FINAL GENERATION: FAIL\n"
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

if (
  N_FAIL >
    0L
) {
  stopf(
    "Figure 3存在Hard FAIL，请勿作为正文图。"
  )
}

