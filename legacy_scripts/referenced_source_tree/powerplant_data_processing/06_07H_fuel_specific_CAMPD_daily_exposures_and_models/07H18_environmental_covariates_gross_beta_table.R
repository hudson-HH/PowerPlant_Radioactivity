#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999,
  width = 220
)

required_packages <- c(
  "data.table",
  "fixest",
  "lme4"
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
  library(lme4)
})

N_CORES <- suppressWarnings(
  parallel::detectCores(logical = TRUE)
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
# 07H18
#
# Environmental covariates and gross-beta radioactivity
#
# 目的：
# 生成类似参考论文 Table 2 的环境因素关联表。
#
# A. 正式推荐：
#    Monitor × year FE
#    与当前07H16 Figure4 / Figure5主模型一致
#
# B. 参考论文可比：
#    Monitor random-intercept model
#    可估计surface U、latitude、year等时间不变/长期项
#
# 注意：
# 这是环境协变量辅助分析。
# 不加入某一种特定发电站暴露，避免人为选择某个fuel/radius。
#
# 正式样本：
# 07F frozen filter-period sample
# 107,407 rows / 133 monitors
# ============================================================


# ============================================================
# 1. Paths
# ============================================================

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

PP <- file.path(
  ROOT,
  "发电站数据处理"
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

HROOT <- file.path(
  PP,
  "06_07H_燃料特异CAMPD逐日暴露与模型"
)

OUT <- file.path(
  HROOT,
  "07H18_环境协变量与gross_beta关联表"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Constants
# ============================================================

YEARS <- 2001:2024

EXPECTED_FP_ROWS <- 107407L
EXPECTED_MONITORS <- 133L


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

  ans[is.na(ans)] <- FALSE

  ans
}

normalize_term <- function(x) {

  y <- gsub(
    "`",
    "",
    x,
    fixed = TRUE
  )

  y[
    y ==
      "lat_diff_35:sunspot_monthly"
  ] <-
    "sunspot_monthly:lat_diff_35"

  y[
    y ==
      "soilm:surface_u_ppm_article"
  ] <-
    "surface_u_ppm_article:soilm"

  y
}

fmt_num <- function(x) {

  ans <- rep(
    NA_character_,
    length(x)
  )

  good <- is.finite(x)

  ans[good] <- formatC(
    x[good],
    digits = 4,
    format = "g"
  )

  ans
}

extract_fixest <- function(fit) {

  C <- as.matrix(
    coeftable(fit)
  )

  D <- data.table(
    term =
      rownames(C),

    estimate =
      as.numeric(
        C[, 1]
      ),

    se =
      as.numeric(
        C[, 2]
      )
  )

  if (
    ncol(C) >= 4L
  ) {
    D[
      ,
      p_wald :=
        as.numeric(
          C[, 4]
        )
    ]
  } else {
    D[
      ,
      p_wald :=
        2 *
        pnorm(
          abs(
            estimate /
            se
          ),
          lower.tail = FALSE
        )
    ]
  }

  D[
    ,
    `:=`(
      ci95_low =
        estimate -
        1.96 *
        se,

      ci95_high =
        estimate +
        1.96 *
        se,

      term_key =
        normalize_term(
          term
        )
    )
  ]

  D
}

extract_lmer <- function(fit) {

  C <- as.matrix(
    summary(fit)$coefficients
  )

  D <- data.table(
    term =
      rownames(C),

    estimate =
      as.numeric(
        C[, "Estimate"]
      ),

    se =
      as.numeric(
        C[, "Std. Error"]
      )
  )

  D <- D[
    term !=
      "(Intercept)"
  ]

  D[
    ,
    `:=`(
      ci95_low =
        estimate -
        1.96 *
        se,

      ci95_high =
        estimate +
        1.96 *
        se,

      p_wald =
        2 *
        pnorm(
          abs(
            estimate /
            se
          ),
          lower.tail = FALSE
        ),

      term_key =
        normalize_term(
          term
        )
    )
  ]

  D
}


# ============================================================
# 4. Exact 133-monitor universe
# ============================================================

cat(
  "============================================================\n",
  "A. 恢复正式133-monitor universe\n",
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

      if (!file.exists(f)) {
        stopf(
          "缺少07F quasi-daily文件：",
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
        "flag_07f_primary_quasi_sample"
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
          "年quasi-daily文件缺字段：",
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
    "正式monitor数不是133；实际=",
    length(
      ELIGIBLE_LOCS
    )
  )
}

cat(
  "Formal eligible monitors = 133：PASS\n"
)


# ============================================================
# 5. Read filter-period data
# ============================================================

cat(
  "\n============================================================\n",
  "B. 构建正式filter-period样本\n",
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
  "year_model",
  "temp",
  "flag_conus",
  "fpmean_campd_upwind_heat_observed_20km_90deg"
)

HEADER <- names(
  fread(
    FP_FILE,
    nrows = 0,
    showProgress = FALSE
  )
)

MISS <- setdiff(
  FP_COLS,
  HEADER
)

if (
  length(MISS) > 0L
) {
  stopf(
    "07F filter-period缺字段：",
    paste(
      MISS,
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
  "year_model",
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


# ============================================================
# 6. Exact frozen formal sample
# ============================================================

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
    "正式filter-period不是107407；实际=",
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
    "正式filter-period不是133 monitors。"
  )
}

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
# 7. Model variables
# ============================================================

FP[
  ,
  `:=`(
    pm_centered_10 =
      pm_filter_period_official -
      10,

    temp_centered =
      temp -
      273.15,

    temp_centered_sq =
      (temp - 273.15)^2,

    temp_centered_cu =
      (temp - 273.15)^3,

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

FP[
  ,
  loc_num_factor :=
    factor(
      loc_num
    )
]

MODEL_VARS <- c(
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
  "year_model",
  "temp_centered",
  "temp_centered_sq",
  "temp_centered_cu",
  "loc_num",
  "loc_num_factor",
  "calendar_year",
  "monitor_year_fe"
)

MODEL <- FP[
  complete.cases(
    FP[
      ,
      ..MODEL_VARS
    ]
  )
]

cat(
  "Environmental-model complete cases = ",
  format(
    nrow(MODEL),
    big.mark = ","
  ),
  "\n",
  "Monitors = ",
  uniqueN(
    MODEL$loc_num
  ),
  "\n",
  sep = ""
)

if (
  nrow(MODEL) ==
    0L
) {
  stopf(
    "环境模型没有complete cases。"
  )
}


# ============================================================
# 8. Covariate descriptive audit
# ============================================================

DESC_VARS <- c(
  "gross_beta_mBq_m3",
  "c_prop",
  "surface_u_ppm_article",
  "sunspot_monthly",
  "pm_filter_period_official",
  "inv_hpbl",
  "rhum",
  "soilm",
  "vel",
  "lat_diff_35",
  "year_model",
  "temp_centered"
)

DESC <- rbindlist(
  lapply(
    DESC_VARS,
    function(v) {

      x <- safe_num(
        MODEL[[v]]
      )

      data.table(
        variable = v,
        n = sum(
          is.finite(x)
        ),
        mean = mean(
          x,
          na.rm = TRUE
        ),
        sd = sd(
          x,
          na.rm = TRUE
        ),
        q25 = unname(
          quantile(
            x,
            0.25,
            na.rm = TRUE
          )
        ),
        median = median(
          x,
          na.rm = TRUE
        ),
        q75 = unname(
          quantile(
            x,
            0.75,
            na.rm = TRUE
          )
        ),
        min = min(
          x,
          na.rm = TRUE
        ),
        max = max(
          x,
          na.rm = TRUE
        )
      )
    }
  )
)

fwrite(
  DESC,
  file.path(
    OUT,
    "01_covariate_descriptive_audit.csv"
  )
)


# ============================================================
# 9. Term dictionary
# ============================================================

TERM_META <- data.table(
  term_key = c(
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
    "I(year_model^2)",
    "temp_centered",
    "temp_centered_sq",
    "temp_centered_cu",
    "sunspot_monthly:lat_diff_35",
    "surface_u_ppm_article:soilm"
  ),

  term_label = c(
    "Continental air-mass proportion",
    "Surface uranium",
    "Monthly sunspot number",
    "PM2.5",
    "Inverse planetary boundary-layer height",
    "Relative humidity",
    "Soil moisture",
    "Wind speed",
    "Latitude difference from 35 degrees",
    "Calendar time",
    "Calendar time squared",
    "Temperature",
    "Temperature squared",
    "Temperature cubed",
    "Sunspots × latitude",
    "Surface uranium × soil moisture"
  ),

  model_unit = c(
    "per 1.0 proportion",
    "per 1 ppm",
    "per 1 sunspot",
    "per 1 µg/m³",
    "per 1 native inverse-PBLH unit",
    "per 1 percentage point",
    "per 1 native soil-moisture unit",
    "per 1 m/s",
    "per 1 degree",
    "per 1 model-year unit",
    "quadratic calendar-time coefficient",
    "linear temperature coefficient",
    "quadratic temperature coefficient",
    "cubic temperature coefficient",
    "interaction coefficient",
    "interaction coefficient"
  ),

  description = c(
    "Continental air-mass contribution",
    "Natural radioactive background at the monitoring location",
    "Solar activity",
    "Filter-period fine particulate matter",
    "Atmospheric mixing and dilution",
    "Atmospheric moisture",
    "Surface soil-moisture condition",
    "Atmospheric transport",
    "Spatial modifier used in the sunspot interaction",
    "Long-term calendar trend",
    "Non-linear long-term calendar trend",
    "Meteorological temperature term",
    "Meteorological temperature polynomial term",
    "Meteorological temperature polynomial term",
    "Modification of the sunspot association by latitude",
    "Modification of the soil-moisture association by surface uranium"
  )
)


# ============================================================
# 10. Primary environmental model
#
# Monitor × year FE
#
# Important:
# surface U, latitude, year/year² main effects cannot be
# separately estimated because they are absorbed by the FE.
#
# Their interactions with time-varying variables remain
# estimable.
# ============================================================

cat(
  "\n============================================================\n",
  "C. Primary environmental model: Monitor × year FE\n",
  "============================================================\n",
  sep = ""
)

FORM_FE <- as.formula(
  paste0(
    "gross_beta_mBq_m3 ~ ",
    "c_prop + ",
    "sunspot_monthly + ",
    "pm_centered_10 + ",
    "inv_hpbl + ",
    "rhum + ",
    "soilm + ",
    "vel + ",
    "temp_centered + ",
    "temp_centered_sq + ",
    "temp_centered_cu + ",
    "sunspot_monthly:lat_diff_35 + ",
    "surface_u_ppm_article:soilm ",
    "| monitor_year_fe"
  )
)

FIT_FE <- feols(
  FORM_FE,
  data = MODEL,
  vcov = ~loc_num,
  notes = FALSE,
  warn = TRUE
)

FE_RAW <- extract_fixest(
  FIT_FE
)

FE_TABLE <- merge(
  FE_RAW,
  TERM_META,
  by = "term_key",
  all.x = TRUE,
  sort = FALSE
)

FE_TABLE[
  ,
  `:=`(
    model =
      "Primary: Monitor × year FE",

    n_rows =
      nobs(
        FIT_FE
      ),

    n_monitors =
      uniqueN(
        MODEL$loc_num
      ),

    estimate_text =
      fmt_num(
        estimate
      ),

    ci95_text =
      paste0(
        "(",
        fmt_num(
          ci95_low
        ),
        ", ",
        fmt_num(
          ci95_high
        ),
        ")"
      ),

    estimate_95CI =
      paste0(
        fmt_num(
          estimate
        ),
        " ",
        "(",
        fmt_num(
          ci95_low
        ),
        ", ",
        fmt_num(
          ci95_high
        ),
        ")"
      )
  )
]

FE_EXPECTED <- c(
  "c_prop",
  "sunspot_monthly",
  "pm_centered_10",
  "inv_hpbl",
  "rhum",
  "soilm",
  "vel",
  "temp_centered",
  "temp_centered_sq",
  "temp_centered_cu",
  "sunspot_monthly:lat_diff_35",
  "surface_u_ppm_article:soilm"
)

FE_MISSING <- setdiff(
  FE_EXPECTED,
  FE_TABLE$term_key
)

if (
  length(FE_MISSING) > 0L
) {

  stopf(
    "Primary FE模型缺少预期系数：",
    paste(
      FE_MISSING,
      collapse = ", "
    )
  )
}

FE_PUBLICATION <- FE_TABLE[
  match(
    FE_EXPECTED,
    term_key
  ),
  .(
    Environmental_factor =
      term_label,

    Model_unit =
      model_unit,

    Estimate =
      estimate,

    CI95_lower =
      ci95_low,

    CI95_upper =
      ci95_high,

    Estimate_95CI =
      estimate_95CI,

    Description =
      description
  )
]

fwrite(
  FE_RAW,
  file.path(
    OUT,
    "02_primary_monitorYearFE_raw_coefficients.csv"
  )
)

fwrite(
  FE_PUBLICATION,
  file.path(
    OUT,
    "03_TableS_environmental_covariates_primary_monitorYearFE.csv"
  )
)


# ============================================================
# 11. Reference-comparable random-intercept model
#
# This is NOT the current primary causal/associational model.
# It is generated only because the reference article's Table 2
# includes time-invariant surface U, latitude and calendar trend.
# ============================================================

cat(
  "\n============================================================\n",
  "D. Reference-comparable random-intercept model\n",
  "============================================================\n",
  sep = ""
)

FORM_RI <- as.formula(
  paste0(
    "gross_beta_mBq_m3 ~ ",
    "c_prop + ",
    "surface_u_ppm_article + ",
    "sunspot_monthly + ",
    "pm_centered_10 + ",
    "inv_hpbl + ",
    "rhum + ",
    "soilm + ",
    "vel + ",
    "lat_diff_35 + ",
    "year_model + ",
    "I(year_model^2) + ",
    "temp_centered + ",
    "temp_centered_sq + ",
    "temp_centered_cu + ",
    "sunspot_monthly:lat_diff_35 + ",
    "surface_u_ppm_article:soilm + ",
    "(1 | loc_num_factor)"
  )
)

RI_STATUS <- "PASS"
RI_NOTE <- ""

FIT_RI <- tryCatch(
  lmer(
    FORM_RI,
    data = MODEL,
    REML = FALSE,
    control =
      lmerControl(
        optimizer = "bobyqa",
        optCtrl =
          list(
            maxfun = 200000
          ),
        check.conv.singular =
          "ignore"
      )
  ),
  error = function(e) {
    e
  }
)

if (
  inherits(
    FIT_RI,
    "error"
  )
) {

  RI_STATUS <- "FAIL"

  RI_NOTE <- conditionMessage(
    FIT_RI
  )

  RI_RAW <- data.table()

  RI_PUBLICATION <- data.table()

} else {

  RI_RAW <- extract_lmer(
    FIT_RI
  )

  RI_TABLE <- merge(
    RI_RAW,
    TERM_META,
    by = "term_key",
    all.x = TRUE,
    sort = FALSE
  )

  RI_TABLE[
    ,
    `:=`(
      model =
        "Reference-comparable random-intercept",

      n_rows =
        nobs(
          FIT_RI
        ),

      n_monitors =
        uniqueN(
          MODEL$loc_num
        ),

      estimate_text =
        fmt_num(
          estimate
        ),

      ci95_text =
        paste0(
          "(",
          fmt_num(
            ci95_low
          ),
          ", ",
          fmt_num(
            ci95_high
          ),
          ")"
        ),

      estimate_95CI =
        paste0(
          fmt_num(
            estimate
          ),
          " ",
          "(",
          fmt_num(
            ci95_low
          ),
          ", ",
          fmt_num(
            ci95_high
          ),
          ")"
        )
    )
  ]

  RI_EXPECTED <- TERM_META$term_key

  RI_MISSING <- setdiff(
    RI_EXPECTED,
    RI_TABLE$term_key
  )

  if (
    length(RI_MISSING) > 0L
  ) {

    RI_STATUS <- "FAIL"

    RI_NOTE <- paste0(
      "Missing coefficients: ",
      paste(
        RI_MISSING,
        collapse = ", "
      )
    )
  }

  RI_PUBLICATION <- RI_TABLE[
    match(
      TERM_META$term_key,
      term_key
    ),
    .(
      Environmental_factor =
        term_label,

      Model_unit =
        model_unit,

      Estimate =
        estimate,

      CI95_lower =
        ci95_low,

      CI95_upper =
        ci95_high,

      Estimate_95CI =
        estimate_95CI,

      Description =
        description
    )
  ]

  fwrite(
    RI_RAW,
    file.path(
      OUT,
      "04_reference_randomIntercept_raw_coefficients.csv"
    )
  )

  fwrite(
    RI_PUBLICATION,
    file.path(
      OUT,
      "05_TableS_environmental_covariates_reference_randomIntercept.csv"
    )
  )
}


# ============================================================
# 12. Model audit
# ============================================================

AUDIT <- data.table(
  item = c(
    "formal_filter_rows",
    "formal_filter_monitors",
    "environment_model_complete_rows",
    "environment_model_monitors",
    "monitor_year_FE_groups",
    "primary_FE_terms",
    "reference_random_intercept_status"
  ),

  value = c(
    as.character(
      nrow(FP)
    ),

    as.character(
      uniqueN(
        FP$loc_num
      )
    ),

    as.character(
      nrow(MODEL)
    ),

    as.character(
      uniqueN(
        MODEL$loc_num
      )
    ),

    as.character(
      uniqueN(
        MODEL$monitor_year_fe
      )
    ),

    as.character(
      nrow(
        FE_RAW
      )
    ),

    RI_STATUS
  ),

  note = c(
    "Frozen formal filter-period sample",
    "Frozen eligible monitor universe",
    "Complete cases for all environmental covariates",
    "Monitors retained in environmental model",
    "Monitor × calendar-year fixed-effect groups",
    "Expected primary estimable environmental coefficients = 12",
    RI_NOTE
  )
)

fwrite(
  AUDIT,
  file.path(
    OUT,
    "06_model_audit.csv"
  )
)


# ============================================================
# 13. Hard checks
# ============================================================

HARD <- data.table(
  check = c(
    "formal_filter_rows_107407",
    "formal_filter_monitors_133",
    "environment_model_rows_positive",
    "environment_model_monitors_133",
    "primary_FE_expected_terms_12",
    "primary_FE_all_finite",
    "reference_random_intercept"
  ),

  severity = c(
    "HARD",
    "HARD",
    "HARD",
    "HARD",
    "HARD",
    "HARD",
    "INFO"
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
      nrow(MODEL)
    ),

    as.character(
      uniqueN(
        MODEL$loc_num
      )
    ),

    as.character(
      length(
        intersect(
          FE_EXPECTED,
          FE_TABLE$term_key
        )
      )
    ),

    as.character(
      sum(
        is.finite(
          FE_TABLE$estimate
        ) &
        is.finite(
          FE_TABLE$ci95_low
        ) &
        is.finite(
          FE_TABLE$ci95_high
        )
      )
    ),

    RI_STATUS
  ),

  expected = c(
    "107407",
    "133",
    ">0",
    "133",
    "12",
    "12",
    "informational"
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
      nrow(MODEL) >
        0L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      uniqueN(
        MODEL$loc_num
      ) ==
        EXPECTED_MONITORS,
      "PASS",
      "FAIL"
    ),

    ifelse(
      length(
        intersect(
          FE_EXPECTED,
          FE_TABLE$term_key
        )
      ) ==
        12L,
      "PASS",
      "FAIL"
    ),

    ifelse(
      all(
        is.finite(
          FE_TABLE$estimate
        ) &
        is.finite(
          FE_TABLE$ci95_low
        ) &
        is.finite(
          FE_TABLE$ci95_high
        )
      ),
      "PASS",
      "FAIL"
    ),

    ifelse(
      RI_STATUS ==
        "PASS",
      "PASS",
      "INFO"
    )
  )
]

fwrite(
  HARD,
  file.path(
    OUT,
    "07_07H18_hard_checks.csv"
  )
)


# ============================================================
# 14. Documentation
# ============================================================

NOTE <- c(
  "# 07H18 环境协变量与 gross-beta 关联表",
  "",
  "## 正式推荐结果",
  "",
  "- 03_TableS_environmental_covariates_primary_monitorYearFE.csv",
  "- 使用正式 filter-period 样本。",
  "- 使用 monitor × calendar-year fixed effects。",
  "- 标准误按 RadNet monitor 聚类。",
  "- 不加入某一种特定发电站暴露，避免人为选择 fuel / radius。",
  "- 这是环境协变量辅助关联模型，不替代发电站主模型。",
  "",
  "## 为什么没有 surface U、latitude、year 的主效应",
  "",
  "- Surface uranium 和 latitude 在站点内基本不随时间变化。",
  "- Calendar year 在 monitor × year fixed effects 中已经被完全吸收。",
  "- 因此这些主效应不能在 primary FE 模型中单独识别。",
  "- Sunspots × latitude 与 surface uranium × soil moisture 因包含时间变化成分，仍可估计。",
  "",
  "## 参考论文可比结果",
  "",
  "- 05_TableS_environmental_covariates_reference_randomIntercept.csv",
  "- 仅用于和参考论文类似的随机截距环境因素表进行比较。",
  "- 它不是当前正文 Figure 4–6 的主模型。",
  "",
  "## 解释原则",
  "",
  "- 表内结果为 adjusted associations，不解释为因果效应。",
  "- Temperature 多项式各项不能脱离整个多项式单独解释为温度边际效应。",
  "- 含交互项时，main coefficient 为交互变量取0时的条件系数。",
  "- 正式论文优先报告 primary Monitor × year FE 表。"
)

writeLines(
  NOTE,
  file.path(
    OUT,
    "00_07H18说明.md"
  ),
  useBytes = TRUE
)


# ============================================================
# 15. Terminal output
# ============================================================

cat(
  "\n============================================================\n",
  "PRIMARY TABLE: MONITOR × YEAR FE\n",
  "============================================================\n",
  sep = ""
)

print(
  FE_PUBLICATION,
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "REFERENCE-COMPARABLE RANDOM-INTERCEPT TABLE\n",
  "============================================================\n",
  sep = ""
)

if (
  RI_STATUS ==
    "PASS"
) {

  print(
    RI_PUBLICATION,
    nrows = Inf
  )

} else {

  cat(
    "Random-intercept comparator failed/skipped:\n",
    RI_NOTE,
    "\n",
    sep = ""
  )
}

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

N_HARD_FAIL <- HARD[
  severity ==
    "HARD" &
  status ==
    "FAIL",
  .N
]

cat(
  "\n============================================================\n",
  ifelse(
    N_HARD_FAIL ==
      0L,
    "07H18 ENVIRONMENTAL COVARIATE TABLE: PASS\n",
    "07H18 ENVIRONMENTAL COVARIATE TABLE: FAIL\n"
  ),
  "Hard FAIL = ",
  N_HARD_FAIL,
  "\n",
  "Reference random-intercept status = ",
  RI_STATUS,
  "\n",
  "Output = ",
  OUT,
  "\n",
  "============================================================\n",
  sep = ""
)

if (
  N_HARD_FAIL >
    0L
) {
  stopf(
    "07H18存在Hard FAIL，请检查07_07H18_hard_checks.csv。"
  )
}
