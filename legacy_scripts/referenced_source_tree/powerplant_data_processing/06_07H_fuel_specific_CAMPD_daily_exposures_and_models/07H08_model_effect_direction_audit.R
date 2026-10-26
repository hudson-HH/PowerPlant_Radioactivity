#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999
)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("缺少R包 data.table。", call. = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
})

stopf <- function(...) {
  stop(paste0(...), call. = FALSE)
}

assert_file <- function(file, label) {
  if (!file.exists(file)) {
    stopf("找不到", label, "：\n", file)
  }
}

assert_cols <- function(D, cols, label) {
  miss <- setdiff(cols, names(D))
  if (length(miss) > 0L) {
    stopf(
      label,
      "缺少字段：",
      paste(miss, collapse = ", ")
    )
  }
}

collapse_radius <- function(radius, keep) {
  z <- sort(
    unique(
      radius[
        keep %in% TRUE &
        !is.na(radius)
      ]
    )
  )

  if (length(z) == 0L) {
    return("")
  }

  paste(z, collapse = ";")
}

direction_pattern <- function(
  n_pass,
  n_positive,
  n_negative
) {

  if (n_pass == 0L) {
    return("NO_ESTIMABLE_MODEL")
  }

  if (n_positive == n_pass) {
    return("ALL_POSITIVE")
  }

  if (n_negative == n_pass) {
    return("ALL_NEGATIVE")
  }

  if (n_positive > n_negative) {
    return("POSITIVE_MAJORITY")
  }

  if (n_negative > n_positive) {
    return("NEGATIVE_MAJORITY")
  }

  "MIXED_BALANCED"
}

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

HROOT <- file.path(
  ROOT,
  "发电站数据处理",
  "06_07H_燃料特异CAMPD逐日暴露与模型"
)

H07 <- file.path(
  HROOT,
  "07H07_论文结果表与图形准备"
)

TABLE_DIR <- file.path(
  H07,
  "01_论文结果表"
)

FILE_SINGLE <- file.path(
  TABLE_DIR,
  "01_论文主结果_upwind单暴露_全部预设半径.csv"
)

FILE_JOINT <- file.path(
  TABLE_DIR,
  "02_论文补充结果_Joint_Count_Heat_upwind_全部半径.csv"
)

FILE_HARD <- file.path(
  H07,
  "07H07_hard_checks.csv"
)

FILE_FREEZE <- file.path(
  H07,
  "07H07_FINAL_FREEZE.md"
)

OUT <- file.path(
  HROOT,
  "07H08_模型效应方向与正负关联审计"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

assert_file(
  FILE_SINGLE,
  "07H07 upwind单暴露结果"
)

assert_file(
  FILE_JOINT,
  "07H07 joint结果"
)

assert_file(
  FILE_HARD,
  "07H07 hard checks"
)

assert_file(
  FILE_FREEZE,
  "07H07 FINAL FREEZE marker"
)

HARD07 <- fread(
  FILE_HARD,
  showProgress = FALSE
)

assert_cols(
  HARD07,
  c("check", "status"),
  "07H07 hard checks"
)

if (HARD07[status != "PASS", .N] > 0L) {
  stopf("07H07存在未PASS hard check。")
}

FREEZE_TEXT <- readLines(
  FILE_FREEZE,
  warn = FALSE,
  encoding = "UTF-8"
)

FREEZE_PASS <- any(
  grepl(
    "07H07 FINAL FREEZE: PASS",
    FREEZE_TEXT,
    fixed = TRUE
  )
)

if (!FREEZE_PASS) {
  stopf("未确认07H07 FINAL FREEZE: PASS。")
}

S <- fread(
  FILE_SINGLE,
  showProgress = FALSE
)

J <- fread(
  FILE_JOINT,
  showProgress = FALSE
)

REQ_SINGLE <- c(
  "branch",
  "branch_label",
  "fuel",
  "fuel_label",
  "radius_km",
  "component",
  "component_label",
  "status",
  "raw_increment",
  "raw_unit",
  "estimate_fixed",
  "ci_low_fixed",
  "ci_high_fixed",
  "estimate_1SD",
  "ci_low_1SD",
  "ci_high_1SD",
  "p_wald"
)

assert_cols(
  S,
  REQ_SINGLE,
  "07H07单暴露结果"
)

if (nrow(S) != 480L) {
  stopf(
    "07H07单暴露表应为480行，实际=",
    nrow(S)
  )
}

FUELS <- c(
  "coal",
  "gas",
  "oil",
  "other_fossil",
  "biomass",
  "solid_waste"
)

BRANCHES <- c(
  "filter_period",
  "quasi_daily"
)

RADII <- c(
  20L,
  25L,
  30L,
  35L,
  40L,
  45L,
  50L,
  100L,
  150L,
  200L
)

RAW_COMPONENTS <- c(
  "raw_count_only",
  "raw_heat_only"
)

WEIGHTED_COMPONENTS <- c(
  "weighted_proximity",
  "weighted_heat_proximity"
)

S[
  ,
  radius_km :=
    as.integer(radius_km)
]

NUM <- c(
  "estimate_fixed",
  "ci_low_fixed",
  "ci_high_fixed",
  "estimate_1SD",
  "ci_low_1SD",
  "ci_high_1SD",
  "p_wald"
)

for (cc in NUM) {
  S[
    ,
    (cc) :=
      suppressWarnings(
        as.numeric(get(cc))
      )
  ]
}

S[
  ,
  `:=`(
    effect_direction =
      fifelse(
        status != "PASS",
        "NOT_ESTIMABLE",
        fifelse(
          estimate_fixed > 0,
          "POSITIVE",
          fifelse(
            estimate_fixed < 0,
            "NEGATIVE",
            "ZERO"
          )
        )
      ),

    ci_direction =
      fifelse(
        status != "PASS",
        "NOT_ESTIMABLE",
        fifelse(
          is.finite(ci_low_fixed) &
          ci_low_fixed > 0,
          "CI_ABOVE_ZERO",
          fifelse(
            is.finite(ci_high_fixed) &
            ci_high_fixed < 0,
            "CI_BELOW_ZERO",
            "CI_INCLUDES_ZERO"
          )
        )
      )
  )
]

RAW <- S[
  component %in% RAW_COMPONENTS
]

WEIGHTED <- S[
  component %in% WEIGHTED_COMPONENTS
]

if (nrow(RAW) != 240L) {
  stopf(
    "raw Count/Heat upwind结果应为240行，实际=",
    nrow(RAW)
  )
}

if (nrow(WEIGHTED) != 240L) {
  stopf(
    "weighted upwind结果应为240行，实际=",
    nrow(WEIGHTED)
  )
}

RAW_SUM <- RAW[
  ,
  {

    np <- sum(status == "PASS")

    pos <- sum(
      status == "PASS" &
      estimate_fixed > 0,
      na.rm = TRUE
    )

    neg <- sum(
      status == "PASS" &
      estimate_fixed < 0,
      na.rm = TRUE
    )

    zero <- sum(
      status == "PASS" &
      estimate_fixed == 0,
      na.rm = TRUE
    )

    cip <- sum(
      status == "PASS" &
      is.finite(ci_low_fixed) &
      ci_low_fixed > 0,
      na.rm = TRUE
    )

    cin <- sum(
      status == "PASS" &
      is.finite(ci_high_fixed) &
      ci_high_fixed < 0,
      na.rm = TRUE
    )

    cic <- sum(
      status == "PASS" &
      is.finite(ci_low_fixed) &
      is.finite(ci_high_fixed) &
      ci_low_fixed <= 0 &
      ci_high_fixed >= 0,
      na.rm = TRUE
    )

    list(
      n_prespecified_radii = .N,
      n_pass = np,
      n_skip = sum(status == "SKIP"),
      n_positive = pos,
      n_negative = neg,
      n_zero = zero,
      positive_pct_among_PASS =
        if (np > 0L) 100 * pos / np else NA_real_,
      negative_pct_among_PASS =
        if (np > 0L) 100 * neg / np else NA_real_,
      n_CI_above_zero = cip,
      n_CI_below_zero = cin,
      n_CI_include_zero = cic,
      positive_radii_km =
        collapse_radius(
          radius_km,
          status == "PASS" &
          estimate_fixed > 0
        ),
      negative_radii_km =
        collapse_radius(
          radius_km,
          status == "PASS" &
          estimate_fixed < 0
        ),
      CI_above_zero_radii_km =
        collapse_radius(
          radius_km,
          status == "PASS" &
          is.finite(ci_low_fixed) &
          ci_low_fixed > 0
        ),
      CI_below_zero_radii_km =
        collapse_radius(
          radius_km,
          status == "PASS" &
          is.finite(ci_high_fixed) &
          ci_high_fixed < 0
        ),
      descriptive_direction_pattern =
        direction_pattern(
          np,
          pos,
          neg
        )
    )
  },
  by = .(
    branch,
    branch_label,
    fuel,
    fuel_label,
    component,
    component_label
  )
]

setorder(
  RAW_SUM,
  branch,
  fuel,
  component
)

fwrite(
  RAW_SUM,
  file.path(
    OUT,
    "01_07H08_raw主分析按燃料正负方向汇总.csv"
  )
)

RAW20 <- RAW[
  radius_km == 20L,
  .(
    branch,
    branch_label,
    fuel,
    fuel_label,
    component,
    component_label,
    radius_km,
    status,
    raw_increment,
    raw_unit,
    estimate_fixed,
    ci_low_fixed,
    ci_high_fixed,
    estimate_1SD,
    ci_low_1SD,
    ci_high_1SD,
    p_wald,
    effect_direction,
    ci_direction
  )
]

setorder(
  RAW20,
  branch,
  fuel,
  component
)

if (nrow(RAW20) != 24L) {
  stopf(
    "20 km raw主分析应为24行，实际=",
    nrow(RAW20)
  )
}

fwrite(
  RAW20,
  file.path(
    OUT,
    "02_07H08_20km_raw主分析逐燃料结果.csv"
  )
)

WEIGHTED_SUM <- WEIGHTED[
  ,
  {

    np <- sum(status == "PASS")

    pos <- sum(
      status == "PASS" &
      estimate_fixed > 0,
      na.rm = TRUE
    )

    neg <- sum(
      status == "PASS" &
      estimate_fixed < 0,
      na.rm = TRUE
    )

    cip <- sum(
      status == "PASS" &
      is.finite(ci_low_fixed) &
      ci_low_fixed > 0,
      na.rm = TRUE
    )

    cin <- sum(
      status == "PASS" &
      is.finite(ci_high_fixed) &
      ci_high_fixed < 0,
      na.rm = TRUE
    )

    list(
      n_prespecified_radii = .N,
      n_pass = np,
      n_skip = sum(status == "SKIP"),
      n_positive = pos,
      n_negative = neg,
      positive_pct_among_PASS =
        if (np > 0L) 100 * pos / np else NA_real_,
      negative_pct_among_PASS =
        if (np > 0L) 100 * neg / np else NA_real_,
      n_CI_above_zero = cip,
      n_CI_below_zero = cin,
      n_CI_include_zero =
        sum(
          status == "PASS" &
          is.finite(ci_low_fixed) &
          is.finite(ci_high_fixed) &
          ci_low_fixed <= 0 &
          ci_high_fixed >= 0,
          na.rm = TRUE
        ),
      positive_radii_km =
        collapse_radius(
          radius_km,
          status == "PASS" &
          estimate_fixed > 0
        ),
      negative_radii_km =
        collapse_radius(
          radius_km,
          status == "PASS" &
          estimate_fixed < 0
        ),
      CI_above_zero_radii_km =
        collapse_radius(
          radius_km,
          status == "PASS" &
          is.finite(ci_low_fixed) &
          ci_low_fixed > 0
        ),
      CI_below_zero_radii_km =
        collapse_radius(
          radius_km,
          status == "PASS" &
          is.finite(ci_high_fixed) &
          ci_high_fixed < 0
        ),
      descriptive_direction_pattern =
        direction_pattern(
          np,
          pos,
          neg
        )
    )
  },
  by = .(
    branch,
    branch_label,
    fuel,
    fuel_label,
    component,
    component_label
  )
]

setorder(
  WEIGHTED_SUM,
  branch,
  fuel,
  component
)

fwrite(
  WEIGHTED_SUM,
  file.path(
    OUT,
    "03_07H08_distance_weighted按燃料正负方向汇总.csv"
  )
)

NEGATIVE_RAW <- RAW[
  status == "PASS" &
  estimate_fixed < 0,
  .(
    branch,
    branch_label,
    fuel,
    fuel_label,
    radius_km,
    component,
    component_label,
    estimate_fixed,
    ci_low_fixed,
    ci_high_fixed,
    estimate_1SD,
    ci_low_1SD,
    ci_high_1SD,
    p_wald,
    ci_direction
  )
]

setorder(
  NEGATIVE_RAW,
  branch,
  fuel,
  component,
  radius_km
)

fwrite(
  NEGATIVE_RAW,
  file.path(
    OUT,
    "04_07H08_raw主分析全部负向结果.csv"
  )
)

NEGATIVE_CI_RAW <- RAW[
  status == "PASS" &
  is.finite(ci_high_fixed) &
  ci_high_fixed < 0,
  .(
    branch,
    branch_label,
    fuel,
    fuel_label,
    radius_km,
    component,
    component_label,
    estimate_fixed,
    ci_low_fixed,
    ci_high_fixed,
    estimate_1SD,
    ci_low_1SD,
    ci_high_1SD,
    p_wald
  )
]

setorder(
  NEGATIVE_CI_RAW,
  branch,
  fuel,
  component,
  radius_km
)

fwrite(
  NEGATIVE_CI_RAW,
  file.path(
    OUT,
    "05_07H08_raw主分析95CI完全低于0结果.csv"
  )
)

POSITIVE_CI_RAW <- RAW[
  status == "PASS" &
  is.finite(ci_low_fixed) &
  ci_low_fixed > 0,
  .(
    branch,
    branch_label,
    fuel,
    fuel_label,
    radius_km,
    component,
    component_label,
    estimate_fixed,
    ci_low_fixed,
    ci_high_fixed,
    estimate_1SD,
    ci_low_1SD,
    ci_high_1SD,
    p_wald
  )
]

setorder(
  POSITIVE_CI_RAW,
  branch,
  fuel,
  component,
  radius_km
)

fwrite(
  POSITIVE_CI_RAW,
  file.path(
    OUT,
    "06_07H08_raw主分析95CI完全高于0结果.csv"
  )
)

Q <- RAW[
  branch == "quasi_daily",
  .(
    fuel,
    fuel_label,
    radius_km,
    component,
    component_label,
    quasi_status = status,
    quasi_estimate = estimate_fixed,
    quasi_ci_low = ci_low_fixed,
    quasi_ci_high = ci_high_fixed
  )
]

F <- RAW[
  branch == "filter_period",
  .(
    fuel,
    fuel_label,
    radius_km,
    component,
    component_label,
    filter_status = status,
    filter_estimate = estimate_fixed,
    filter_ci_low = ci_low_fixed,
    filter_ci_high = ci_high_fixed
  )
]

BF <- merge(
  Q,
  F,
  by = c(
    "fuel",
    "fuel_label",
    "radius_km",
    "component",
    "component_label"
  ),
  all = TRUE,
  sort = FALSE
)

BF[
  ,
  `:=`(
    both_PASS =
      quasi_status == "PASS" &
      filter_status == "PASS",

    same_direction =
      fifelse(
        quasi_status == "PASS" &
        filter_status == "PASS" &
        is.finite(quasi_estimate) &
        is.finite(filter_estimate),
        sign(quasi_estimate) ==
          sign(filter_estimate),
        NA
      ),

    both_positive =
      quasi_status == "PASS" &
      filter_status == "PASS" &
      quasi_estimate > 0 &
      filter_estimate > 0,

    both_negative =
      quasi_status == "PASS" &
      filter_status == "PASS" &
      quasi_estimate < 0 &
      filter_estimate < 0
  )
]

if (nrow(BF) != 120L) {
  stopf(
    "quasi-filter raw主分析配对应为120行，实际=",
    nrow(BF)
  )
}

setorder(
  BF,
  fuel,
  component,
  radius_km
)

fwrite(
  BF,
  file.path(
    OUT,
    "07_07H08_quasi_filter正负方向一致性.csv"
  )
)

BF_SUM <- BF[
  ,
  .(
    n_radii = .N,
    n_both_PASS =
      sum(both_PASS, na.rm = TRUE),
    n_same_direction =
      sum(same_direction %in% TRUE, na.rm = TRUE),
    n_opposite_direction =
      sum(same_direction %in% FALSE, na.rm = TRUE),
    n_both_positive =
      sum(both_positive, na.rm = TRUE),
    n_both_negative =
      sum(both_negative, na.rm = TRUE)
  ),
  by = .(
    fuel,
    fuel_label,
    component,
    component_label
  )
]

setorder(
  BF_SUM,
  fuel,
  component
)

fwrite(
  BF_SUM,
  file.path(
    OUT,
    "08_07H08_quasi_filter按燃料方向一致性汇总.csv"
  )
)

OVERALL_RAW <- RAW[
  status == "PASS",
  .(
    n_PASS = .N,
    n_positive =
      sum(estimate_fixed > 0, na.rm = TRUE),
    n_negative =
      sum(estimate_fixed < 0, na.rm = TRUE),
    n_CI_above_zero =
      sum(
        is.finite(ci_low_fixed) &
        ci_low_fixed > 0,
        na.rm = TRUE
      ),
    n_CI_below_zero =
      sum(
        is.finite(ci_high_fixed) &
        ci_high_fixed < 0,
        na.rm = TRUE
      ),
    n_CI_include_zero =
      sum(
        is.finite(ci_low_fixed) &
        is.finite(ci_high_fixed) &
        ci_low_fixed <= 0 &
        ci_high_fixed >= 0,
        na.rm = TRUE
      )
  ),
  by = .(
    branch,
    branch_label,
    component,
    component_label
  )
]

fwrite(
  OVERALL_RAW,
  file.path(
    OUT,
    "09_07H08_raw主分析总体方向计数.csv"
  )
)

N_RAW_DUP <- RAW[
  ,
  .N,
  by = .(
    branch,
    fuel,
    radius_km,
    component
  )
][
  N > 1L,
  .N
]

N_WEIGHTED_DUP <- WEIGHTED[
  ,
  .N,
  by = .(
    branch,
    fuel,
    radius_km,
    component
  )
][
  N > 1L,
  .N
]

N_RAW_SIGN_ACCOUNT <- RAW[
  status == "PASS" &
  !(
    is.finite(estimate_fixed) &
    (
      estimate_fixed > 0 |
      estimate_fixed < 0 |
      estimate_fixed == 0
    )
  ),
  .N
]

N_RAW_CI_ACCOUNT <- RAW[
  status == "PASS" &
  !(
    (
      is.finite(ci_low_fixed) &
      ci_low_fixed > 0
    ) |
    (
      is.finite(ci_high_fixed) &
      ci_high_fixed < 0
    ) |
    (
      is.finite(ci_low_fixed) &
      is.finite(ci_high_fixed) &
      ci_low_fixed <= 0 &
      ci_high_fixed >= 0
    )
  ),
  .N
]

HARD08 <- data.table(
  check = c(
    "07H07_final_freeze_pass",
    "single_exposure_rows_480",
    "raw_main_rows_240",
    "weighted_rows_240",
    "raw_20km_rows_24",
    "raw_duplicate_keys_zero",
    "weighted_duplicate_keys_zero",
    "quasi_filter_pair_rows_120",
    "raw_PASS_effect_direction_accounted",
    "raw_PASS_CI_direction_accounted"
  ),

  observed = c(
    as.character(as.integer(FREEZE_PASS)),
    as.character(nrow(S)),
    as.character(nrow(RAW)),
    as.character(nrow(WEIGHTED)),
    as.character(nrow(RAW20)),
    as.character(N_RAW_DUP),
    as.character(N_WEIGHTED_DUP),
    as.character(nrow(BF)),
    as.character(N_RAW_SIGN_ACCOUNT),
    as.character(N_RAW_CI_ACCOUNT)
  ),

  expected = c(
    "1",
    "480",
    "240",
    "240",
    "24",
    "0",
    "0",
    "120",
    "0",
    "0"
  )
)

HARD08[
  ,
  status :=
    fifelse(
      observed == expected,
      "PASS",
      "FAIL"
    )
]

fwrite(
  HARD08,
  file.path(
    OUT,
    "10_07H08_hard_checks.csv"
  )
)

N_HARD_FAIL <- HARD08[
  status == "FAIL",
  .N
]

cat(
  "\n============================================================\n",
  "A. RAW MAIN ANALYSIS: POSITIVE vs NEGATIVE\n",
  "============================================================\n",
  sep = ""
)

print(
  RAW_SUM[
    ,
    .(
      branch,
      fuel,
      component,
      n_pass,
      n_positive,
      n_negative,
      n_CI_above_zero,
      n_CI_below_zero,
      n_CI_include_zero,
      descriptive_direction_pattern
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "B. NEGATIVE-DOMINANT RAW GROUPS\n",
  "============================================================\n",
  sep = ""
)

NEG_DOM <- RAW_SUM[
  n_negative > n_positive
]

if (nrow(NEG_DOM) == 0L) {

  cat(
    "没有发现负向结果占多数的 branch × fuel × exposure 组合。\n"
  )

} else {

  print(
    NEG_DOM[
      ,
      .(
        branch,
        fuel,
        component,
        n_pass,
        n_positive,
        n_negative,
        n_CI_above_zero,
        n_CI_below_zero,
        negative_radii_km,
        CI_below_zero_radii_km
      )
    ],
    nrows = Inf
  )
}

cat(
  "\n============================================================\n",
  "C. 20-km RAW MAIN RESULTS\n",
  "============================================================\n",
  sep = ""
)

print(
  RAW20[
    ,
    .(
      branch,
      fuel,
      component,
      status,
      estimate_fixed,
      ci_low_fixed,
      ci_high_fixed,
      effect_direction,
      ci_direction,
      p_wald
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "D. DISTANCE-WEIGHTED DIRECTION SUMMARY\n",
  "============================================================\n",
  sep = ""
)

print(
  WEIGHTED_SUM[
    ,
    .(
      branch,
      fuel,
      component,
      n_pass,
      n_positive,
      n_negative,
      n_CI_above_zero,
      n_CI_below_zero,
      n_CI_include_zero,
      descriptive_direction_pattern
    )
  ],
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "E. QUASI vs FILTER RAW DIRECTION CONSISTENCY\n",
  "============================================================\n",
  sep = ""
)

print(
  BF_SUM,
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "F. HARD CHECKS\n",
  "============================================================\n",
  sep = ""
)

print(
  HARD08,
  nrows = Inf
)

REPORT <- c(
  "# 07H08 模型效应方向与正负关联审计",
  "",
  paste0(
    "- Generated: ",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  ),
  "",
  "## Purpose",
  "",
  "- Determine whether fuel-specific power-plant exposure is associated with higher or lower RadNet gross beta.",
  "- Raw Count-only and Heat-only upwind models are treated as the main association results.",
  "- Distance-weighted Σ1/d and ΣHeat/d results are treated as sensitivity analyses.",
  "- All ten prespecified cumulative radii are retained.",
  "",
  "## Interpretation",
  "",
  "- POSITIVE means the fitted exposure coefficient is greater than zero.",
  "- NEGATIVE means the fitted exposure coefficient is less than zero.",
  "- CI_ABOVE_ZERO means the entire 95% CI is above zero.",
  "- CI_BELOW_ZERO means the entire 95% CI is below zero.",
  "- CI_INCLUDES_ZERO means the 95% CI crosses or touches zero.",
  "",
  "## Important caution",
  "",
  "- The ten cumulative radii are nested and highly correlated; counts such as 8/10 positive are descriptive radius-pattern summaries, not ten independent statistical replications.",
  "- No radius or fuel is selected according to P value.",
  "- Association direction does not by itself establish causality or atmospheric transport mechanism.",
  "",
  if (N_HARD_FAIL == 0L) {
    "**07H08 DIRECTION AUDIT: PASS**"
  } else {
    "**07H08 DIRECTION AUDIT: FAIL**"
  }
)

writeLines(
  REPORT,
  file.path(
    OUT,
    "11_07H08_正负关联审计说明.md"
  ),
  useBytes = TRUE
)

capture.output(
  sessionInfo(),
  file = file.path(
    OUT,
    "12_sessionInfo.txt"
  )
)

cat(
  "\n============================================================\n"
)

if (N_HARD_FAIL == 0L) {

  cat(
    "07H08 DIRECTION AUDIT: PASS\n"
  )

} else {

  cat(
    "07H08 DIRECTION AUDIT: FAIL\n"
  )
}

cat(
  "============================================================\n",
  "Raw main rows = ",
  nrow(RAW),
  " / 240\n",
  "20-km raw rows = ",
  nrow(RAW20),
  " / 24\n",
  "Weighted rows = ",
  nrow(WEIGHTED),
  " / 240\n",
  "Quasi-filter pairs = ",
  nrow(BF),
  " / 120\n",
  "Hard FAIL = ",
  N_HARD_FAIL,
  "\n",
  sep = ""
)

if (N_HARD_FAIL > 0L) {
  stopf("07H08存在hard failure。")
}

