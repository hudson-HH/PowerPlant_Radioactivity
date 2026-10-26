#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999,
  width = 220
)

# ============================================================
# 07H12 Table 1 正式描述性统计
#
# 目的：
# 1. 仿照参考论文 Table 1 的结构；
# 2. 连续变量统一报告：
#       mean (Q25, Q75)
# 3. 正文比较：
#       - 全部 CONUS 正式数据
#       - 最终主分析样本
# 4. 不重新拟合模型；
# 5. 不修改任何冻结上游文件。
#
# 重要：
# - Table 1 是描述性统计，不做显著性检验。
# - 不进行“有电厂/无电厂”强行二分类。
# - 主分析样本必须严格等于：
#       558,915 rows
#       140 monitors
# ============================================================


# ============================================================
# 0. Package
# ============================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  install.packages(
    "data.table",
    repos = "https://cloud.r-project.org"
  )
}

suppressPackageStartupMessages({
  library(data.table)
})


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

H12 <- file.path(
  ROOT,
  "发电站数据处理",
  "06_07H_燃料特异CAMPD逐日暴露与模型",
  "07H12_论文正式图表总整理"
)

OUT_DIR <- file.path(
  H12,
  "02_正文表"
)

AUDIT_DIR <- file.path(
  H12,
  "06_审计"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  AUDIT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


OUT_EN_CSV <- file.path(
  OUT_DIR,
  "Table1_descriptive_characteristics_EN.csv"
)

OUT_CN_CSV <- file.path(
  OUT_DIR,
  "Table1_descriptive_characteristics_CN.csv"
)

OUT_SOURCE <- file.path(
  OUT_DIR,
  "Table1_source_statistics.csv"
)

OUT_EN_HTML <- file.path(
  OUT_DIR,
  "Table1_descriptive_characteristics_EN.html"
)

OUT_CN_HTML <- file.path(
  OUT_DIR,
  "Table1_descriptive_characteristics_CN.html"
)

OUT_NOTE <- file.path(
  OUT_DIR,
  "Table1_notes.txt"
)

OUT_AUDIT <- file.path(
  AUDIT_DIR,
  "Table1_descriptive_characteristics_audit.csv"
)


# ============================================================
# 2. Frozen expectations
# ============================================================

EXPECTED_ALL_ROWS <- 769770L
EXPECTED_ALL_MONITORS <- 181L

EXPECTED_CONUS_ROWS <- 741861L
EXPECTED_CONUS_MONITORS <- 168L

EXPECTED_MAIN_ROWS <- 558915L
EXPECTED_MAIN_MONITORS <- 140L


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


as_flag <- function(
  x,
  default = FALSE
) {

  if (is.logical(x)) {

    y <- x

  } else if (
    is.numeric(x) ||
    is.integer(x)
  ) {

    y <- x != 0

  } else {

    xx <- tolower(
      trimws(
        as.character(x)
      )
    )

    y <- xx %in%
      c(
        "true",
        "t",
        "1",
        "yes",
        "y"
      )
  }

  y[
    is.na(y)
  ] <- default

  y
}


fmt_int <- function(x) {

  if (
    length(x) == 0L ||
    is.na(x)
  ) {
    return("")
  }

  format(
    as.integer(x),
    big.mark = ",",
    scientific = FALSE,
    trim = TRUE
  )
}


fmt_num <- function(
  x,
  digits
) {

  if (
    length(x) == 0L ||
    !is.finite(x)
  ) {
    return("")
  }

  formatC(
    x,
    format = "f",
    digits = digits
  )
}


html_escape <- function(x) {

  x <- as.character(x)

  x <- gsub(
    "&",
    "&amp;",
    x,
    fixed = TRUE
  )

  x <- gsub(
    "<",
    "&lt;",
    x,
    fixed = TRUE
  )

  x <- gsub(
    ">",
    "&gt;",
    x,
    fixed = TRUE
  )

  x
}


# ============================================================
# 4. Read frozen 08 dataset
# ============================================================

cat(
  "\n============================================================\n",
  "Table 1：读取冻结08正式数据\n",
  "============================================================\n",
  sep = ""
)

if (!file.exists(INPUT08)) {
  stopf(
    "找不到冻结08正式数据：\n",
    INPUT08
  )
}

D <- fread(
  INPUT08,
  showProgress = TRUE
)


# ============================================================
# 5. Basic frozen checks
# ============================================================

required_base <- c(
  "loc_num",
  "flag_conus",
  "flag_reference_main_sample",
  "gross_beta_mBq_m3",
  "surface_u_ppm_article",
  "c_prop",
  "soilm",
  "rhum",
  "vel",
  "sunspot_monthly",
  "all_upwind_n_plants_20km"
)

missing_base <- setdiff(
  required_base,
  names(D)
)

if (length(missing_base) > 0L) {
  stopf(
    "08正式数据缺少Table 1字段：\n",
    paste(
      missing_base,
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
  flag_conus :=
    as_flag(
      flag_conus
    )
]

D[
  ,
  flag_reference_main_sample :=
    as_flag(
      flag_reference_main_sample
    )
]


# ============================================================
# 6. Reconstruct raw display variables if necessary
# ============================================================

# ------------------------------------------------------------
# PM2.5
# 正式08正常情况下已有 pm_c。
# 如果只有 pm_centered_10，则恢复原始 PM2.5：
# pm_c = pm_centered_10 + 10
# ------------------------------------------------------------

if (!"pm_c" %in% names(D)) {

  if (!"pm_centered_10" %in% names(D)) {
    stopf(
      "同时缺少pm_c与pm_centered_10，无法生成PM2.5描述统计。"
    )
  }

  D[
    ,
    pm_c :=
      safe_num(
        pm_centered_10
      ) +
      10
  ]
}


# ------------------------------------------------------------
# Temperature
# 正式08正常情况下已有 temp_C。
# 如果没有，则 temp(K) - 273.15。
# ------------------------------------------------------------

if (!"temp_C" %in% names(D)) {

  if (!"temp" %in% names(D)) {
    stopf(
      "同时缺少temp_C与temp，无法生成温度描述统计。"
    )
  }

  D[
    ,
    temp_C :=
      safe_num(
        temp
      ) -
      273.15
  ]
}


# ------------------------------------------------------------
# PBLH
# 正式08正常情况下已有 pblh_km。
# 如果没有，则 hpbl(m) / 1000。
# ------------------------------------------------------------

if (!"pblh_km" %in% names(D)) {

  if (!"hpbl" %in% names(D)) {
    stopf(
      "同时缺少pblh_km与hpbl，无法生成PBLH描述统计。"
    )
  }

  D[
    ,
    pblh_km :=
      safe_num(
        hpbl
      ) /
      1000
  ]
}


# ============================================================
# 7. Final variables
# ============================================================

TABLE_VARS <- data.table(

  variable = c(
    "gross_beta_mBq_m3",
    "surface_u_ppm_article",
    "c_prop",
    "pm_c",
    "soilm",
    "rhum",
    "temp_C",
    "pblh_km",
    "vel",
    "sunspot_monthly",
    "all_upwind_n_plants_20km"
  ),

  label_en = c(
    "Gross beta radioactivity (mBq/m³)",
    "Surface U-238 proxy (ppm)",
    "Continental air-mass proportion",
    "PM2.5 (µg/m³)",
    "Soil moisture (NARR)",
    "Relative humidity (%)",
    "Temperature (°C)",
    "Planetary boundary layer height (km)",
    "Wind speed (m/s)",
    "Monthly sunspots (n)",
    "Upwind power plants within 20 km (n)"
  ),

  label_cn = c(
    "总β放射性（mBq/m³）",
    "地表U-238代理值（ppm）",
    "大陆来源气团比例",
    "PM2.5（µg/m³）",
    "土壤湿度（NARR）",
    "相对湿度（%）",
    "温度（°C）",
    "行星边界层高度（km）",
    "风速（m/s）",
    "月太阳黑子数（n）",
    "20公里内上风向发电站数量（n）"
  ),

  digits = c(
    3L,
    2L,
    2L,
    2L,
    3L,
    1L,
    2L,
    2L,
    2L,
    2L,
    2L
  )
)


missing_vars <- setdiff(
  TABLE_VARS$variable,
  names(D)
)

if (length(missing_vars) > 0L) {
  stopf(
    "最终Table 1仍缺少变量：\n",
    paste(
      missing_vars,
      collapse = "\n"
    )
  )
}


# ============================================================
# 8. Define descriptive populations
# ============================================================

CONUS <- D[
  flag_conus %in%
    TRUE
]

MAIN <- D[
  flag_reference_main_sample %in%
    TRUE
]


# ============================================================
# 9. Frozen sample checks
# ============================================================

if (
  nrow(D) !=
    EXPECTED_ALL_ROWS
) {
  stopf(
    "08总行数不是冻结值：实际=",
    nrow(D),
    "；预期=",
    EXPECTED_ALL_ROWS
  )
}

if (
  uniqueN(
    D$loc_num
  ) !=
    EXPECTED_ALL_MONITORS
) {
  stopf(
    "08总监测站数不是冻结值。"
  )
}

if (
  nrow(CONUS) !=
    EXPECTED_CONUS_ROWS
) {
  stopf(
    "CONUS行数不是冻结值：实际=",
    nrow(CONUS),
    "；预期=",
    EXPECTED_CONUS_ROWS
  )
}

if (
  uniqueN(
    CONUS$loc_num
  ) !=
    EXPECTED_CONUS_MONITORS
) {
  stopf(
    "CONUS监测站数不是冻结值：实际=",
    uniqueN(CONUS$loc_num),
    "；预期=",
    EXPECTED_CONUS_MONITORS
  )
}

if (
  nrow(MAIN) !=
    EXPECTED_MAIN_ROWS
) {
  stopf(
    "主分析样本行数不是冻结值：实际=",
    nrow(MAIN),
    "；预期=",
    EXPECTED_MAIN_ROWS
  )
}

if (
  uniqueN(
    MAIN$loc_num
  ) !=
    EXPECTED_MAIN_MONITORS
) {
  stopf(
    "主分析监测站数不是冻结值：实际=",
    uniqueN(MAIN$loc_num),
    "；预期=",
    EXPECTED_MAIN_MONITORS
  )
}

cat(
  "全部正式数据：",
  fmt_int(
    nrow(D)
  ),
  " rows / ",
  uniqueN(
    D$loc_num
  ),
  " monitors\n",

  "CONUS：",
  fmt_int(
    nrow(CONUS)
  ),
  " rows / ",
  uniqueN(
    CONUS$loc_num
  ),
  " monitors\n",

  "最终主分析样本：",
  fmt_int(
    nrow(MAIN)
  ),
  " rows / ",
  uniqueN(
    MAIN$loc_num
  ),
  " monitors\n",
  sep = ""
)


# ============================================================
# 10. Main-sample missingness hard check
# ============================================================

MAIN_MISSING <- rbindlist(
  lapply(
    TABLE_VARS$variable,
    function(v) {

      x <- safe_num(
        MAIN[[
          v
        ]]
      )

      data.table(
        variable = v,

        n_total =
          length(x),

        n_finite =
          sum(
            is.finite(x)
          ),

        n_missing_or_nonfinite =
          sum(
            !is.finite(x)
          )
      )
    }
  )
)

cat(
  "\n主分析样本Table 1变量完整性：\n"
)

print(
  MAIN_MISSING,
  nrows = Inf
)

if (
  MAIN_MISSING[
    n_missing_or_nonfinite >
      0,
    .N
  ] >
    0L
) {

  stopf(
    "主分析样本中Table 1变量仍存在缺失/非有限值，禁止生成正式Table 1。"
  )
}


# ============================================================
# 11. Summary function
# ============================================================

summarize_variable <- function(
  DT,
  variable,
  digits,
  group_name
) {

  x <- safe_num(
    DT[[
      variable
    ]]
  )

  finite <- is.finite(
    x
  )

  xx <- x[
    finite
  ]

  if (length(xx) == 0L) {

    return(
      data.table(
        group =
          group_name,

        variable =
          variable,

        n_total =
          length(x),

        n_nonmissing =
          0L,

        mean =
          NA_real_,

        q25 =
          NA_real_,

        median =
          NA_real_,

        q75 =
          NA_real_,

        min =
          NA_real_,

        max =
          NA_real_,

        display =
          ""
      )
    )
  }

  m <- mean(
    xx
  )

  qs <- quantile(
    xx,
    probs = c(
      0.25,
      0.50,
      0.75
    ),
    na.rm = TRUE,
    names = FALSE,
    type = 7
  )

  lo <- min(
    xx
  )

  hi <- max(
    xx
  )

  display <- paste0(
    fmt_num(
      m,
      digits
    ),
    " (",
    fmt_num(
      qs[1],
      digits
    ),
    ", ",
    fmt_num(
      qs[3],
      digits
    ),
    ")"
  )

  data.table(
    group =
      group_name,

    variable =
      variable,

    n_total =
      length(x),

    n_nonmissing =
      length(xx),

    mean =
      m,

    q25 =
      qs[1],

    median =
      qs[2],

    q75 =
      qs[3],

    min =
      lo,

    max =
      hi,

    display =
      display
  )
}


# ============================================================
# 12. Generate source statistics
# ============================================================

SOURCE_LIST <- list()

for (
  i in seq_len(
    nrow(
      TABLE_VARS
    )
  )
) {

  v <- TABLE_VARS$variable[i]
  dg <- TABLE_VARS$digits[i]

  SOURCE_LIST[[length(SOURCE_LIST) + 1L]] <- summarize_variable(
    CONUS,
    v,
    dg,
    "CONUS_all"
  )

  SOURCE_LIST[[length(SOURCE_LIST) + 1L]] <- summarize_variable(
    MAIN,
    v,
    dg,
    "Primary_analysis"
  )
}

SOURCE_STATS <- rbindlist(
  SOURCE_LIST,
  use.names = TRUE,
  fill = TRUE
)

SOURCE_STATS <- merge(
  SOURCE_STATS,
  TABLE_VARS[
    ,
    .(
      variable,
      label_en,
      label_cn,
      digits
    )
  ],
  by =
    "variable",
  all.x =
    TRUE,
  sort =
    FALSE
)

setcolorder(
  SOURCE_STATS,
  c(
    "group",
    "variable",
    "label_en",
    "label_cn",
    "n_total",
    "n_nonmissing",
    "mean",
    "q25",
    "median",
    "q75",
    "min",
    "max",
    "display",
    "digits"
  )
)

fwrite(
  SOURCE_STATS,
  OUT_SOURCE
)


# ============================================================
# 13. Build publication Table 1
# ============================================================

get_display <- function(
  variable_name,
  group_name
) {

  idx <-
    SOURCE_STATS[["variable"]] == variable_name &
    SOURCE_STATS[["group"]] == group_name

  idx[is.na(idx)] <- FALSE

  z <- SOURCE_STATS[["display"]][idx]

  if (length(z) != 1L) {
    stopf(
      "Table 1无法唯一读取：",
      variable_name,
      " / ",
      group_name,
      "；匹配数=",
      length(z)
    )
  }

  z[[1L]]
}


TABLE1_EN <- data.table(

  Variable = c(
    "Monitors (n)",
    "Quasi-daily observations (n)",
    TABLE_VARS$label_en
  ),

  `All CONUS observations` = c(
    fmt_int(
      uniqueN(
        CONUS$loc_num
      )
    ),

    fmt_int(
      nrow(
        CONUS
      )
    ),

    vapply(
      TABLE_VARS$variable,
      function(v) {
        get_display(
          v,
          "CONUS_all"
        )
      },
      character(1)
    )
  ),

  `Primary analysis` = c(
    fmt_int(
      uniqueN(
        MAIN$loc_num
      )
    ),

    fmt_int(
      nrow(
        MAIN
      )
    ),

    vapply(
      TABLE_VARS$variable,
      function(v) {
        get_display(
          v,
          "Primary_analysis"
        )
      },
      character(1)
    )
  )
)


TABLE1_CN <- data.table(

  `变量` = c(
    "监测站数（n）",
    "准逐日观测数（n）",
    TABLE_VARS$label_cn
  ),

  `美国本土全部正式观测` = c(
    fmt_int(
      uniqueN(
        CONUS$loc_num
      )
    ),

    fmt_int(
      nrow(
        CONUS
      )
    ),

    vapply(
      TABLE_VARS$variable,
      function(v) {
        get_display(
          v,
          "CONUS_all"
        )
      },
      character(1)
    )
  ),

  `最终主分析` = c(
    fmt_int(
      uniqueN(
        MAIN$loc_num
      )
    ),

    fmt_int(
      nrow(
        MAIN
      )
    ),

    vapply(
      TABLE_VARS$variable,
      function(v) {
        get_display(
          v,
          "Primary_analysis"
        )
      },
      character(1)
    )
  )
)


# ============================================================
# 14. Save CSV
# ============================================================

fwrite(
  TABLE1_EN,
  OUT_EN_CSV,
  bom = TRUE
)

fwrite(
  TABLE1_CN,
  OUT_CN_CSV,
  bom = TRUE
)


# ============================================================
# 15. HTML generator
# ============================================================

write_table_html <- function(
  DT,
  path,
  title,
  subtitle,
  notes
) {

  headers <- names(
    DT
  )

  html <- c(
    "<!DOCTYPE html>",
    "<html>",
    "<head>",
    "<meta charset=\"UTF-8\">",
    "<style>",
    "body {",
    "  font-family: Arial, Helvetica, sans-serif;",
    "  margin: 36px;",
    "  color: #111111;",
    "}",
    ".title {",
    "  font-size: 19px;",
    "  font-weight: 700;",
    "  margin-bottom: 6px;",
    "}",
    ".subtitle {",
    "  font-size: 13px;",
    "  margin-bottom: 18px;",
    "}",
    "table {",
    "  border-collapse: collapse;",
    "  width: 100%;",
    "  max-width: 1100px;",
    "  font-size: 12.5px;",
    "}",
    "thead th {",
    "  background: #f1e3c4;",
    "  border-top: 2px solid #222222;",
    "  border-bottom: 1.5px solid #222222;",
    "  padding: 8px 10px;",
    "  text-align: center;",
    "}",
    "thead th:first-child {",
    "  text-align: left;",
    "}",
    "tbody td {",
    "  padding: 6px 10px;",
    "  border-bottom: 1px solid #dddddd;",
    "}",
    "tbody td:first-child {",
    "  text-align: left;",
    "}",
    "tbody td:not(:first-child) {",
    "  text-align: center;",
    "}",
    "tbody tr:last-child td {",
    "  border-bottom: 2px solid #222222;",
    "}",
    ".note {",
    "  margin-top: 10px;",
    "  max-width: 1100px;",
    "  font-size: 11px;",
    "  line-height: 1.45;",
    "}",
    "</style>",
    "</head>",
    "<body>",
    paste0(
      "<div class=\"title\">",
      html_escape(
        title
      ),
      "</div>"
    ),
    paste0(
      "<div class=\"subtitle\">",
      html_escape(
        subtitle
      ),
      "</div>"
    ),
    "<table>",
    "<thead>",
    "<tr>",
    paste0(
      "<th>",
      html_escape(
        headers
      ),
      "</th>",
      collapse = ""
    ),
    "</tr>",
    "</thead>",
    "<tbody>"
  )

  for (
    i in seq_len(
      nrow(DT)
    )
  ) {

    vals <- unlist(
      DT[
        i
      ],
      use.names = FALSE
    )

    html <- c(
      html,
      "<tr>",
      paste0(
        "<td>",
        html_escape(
          vals
        ),
        "</td>",
        collapse = ""
      ),
      "</tr>"
    )
  }

  html <- c(
    html,
    "</tbody>",
    "</table>",
    paste0(
      "<div class=\"note\">",
      paste(
        html_escape(
          notes
        ),
        collapse = "<br>"
      ),
      "</div>"
    ),
    "</body>",
    "</html>"
  )

  writeLines(
    html,
    path,
    useBytes = TRUE
  )
}


# ============================================================
# 16. English HTML
# ============================================================

EN_NOTES <- c(
  "Continuous variables are presented as mean (25th percentile, 75th percentile).",
  "CONUS formal data include all continental-U.S. quasi-daily observations in the frozen 08 dataset.",
  "The primary-analysis sample consists of observations with flag_reference_main_sample = TRUE.",
  "Gross-beta measurements from a multi-day RadNet collection period may be repeated across the corresponding quasi-daily rows by design.",
  "The 20-km upwind power-plant count is the daily count within the 90° directional sector.",
  "No hypothesis tests or p-values are presented because Table 1 is descriptive."
)

write_table_html(
  TABLE1_EN,
  OUT_EN_HTML,
  "Table 1. Descriptive characteristics of ambient gross-beta radioactivity, environmental covariates, and power-plant exposure",
  "Continental United States, 2001–2024",
  EN_NOTES
)


# ============================================================
# 17. Chinese HTML
# ============================================================

CN_NOTES <- c(
  "连续变量表示为：均值（第25百分位数，第75百分位数）。",
  "CONUS全部正式数据为冻结08数据集中美国本土范围内的全部准逐日观测。",
  "最终主分析样本定义为 flag_reference_main_sample = TRUE。",
  "RadNet一次多日滤膜采样得到的总β测量值按既定准逐日规则复制至对应日期，因此同一采样周期内可出现相同总β值。",
  "20公里上风向发电站数量为90°方向扇区内的逐日发电站数量。",
  "Table 1仅用于描述样本特征，不进行显著性检验或报告P值。"
)

write_table_html(
  TABLE1_CN,
  OUT_CN_HTML,
  "表1. 环境总β放射性、环境协变量及发电站暴露的描述性特征",
  "美国本土，2001–2024年",
  CN_NOTES
)


# ============================================================
# 18. Plain-text notes
# ============================================================

NOTE_LINES <- c(

  "TABLE 1 FORMAL NOTES",
  "============================================================",
  "",

  "English title:",
  paste0(
    "Table 1. Descriptive characteristics of ambient gross-beta ",
    "radioactivity, environmental covariates, and power-plant exposure"
  ),
  "",

  "Chinese title:",
  "表1. 环境总β放射性、环境协变量及发电站暴露的描述性特征",
  "",

  "Formal comparison columns:",
  paste0(
    "1. CONUS formal data: ",
    uniqueN(
      CONUS$loc_num
    ),
    " monitors; ",
    fmt_int(
      nrow(CONUS)
    ),
    " quasi-daily observations."
  ),
  paste0(
    "2. Primary-analysis sample: ",
    uniqueN(
      MAIN$loc_num
    ),
    " monitors; ",
    fmt_int(
      nrow(MAIN)
    ),
    " quasi-daily observations."
  ),
  "",

  "Continuous-variable format:",
  "Mean (25th percentile, 75th percentile)",
  "",

  "No hypothesis testing is performed in Table 1."
)

writeLines(
  NOTE_LINES,
  OUT_NOTE,
  useBytes = TRUE
)


# ============================================================
# 19. Coverage audit
# ============================================================

COVERAGE <- rbindlist(
  lapply(
    TABLE_VARS$variable,
    function(v) {

      x1 <- safe_num(
        CONUS[[
          v
        ]]
      )

      x2 <- safe_num(
        MAIN[[
          v
        ]]
      )

      data.table(
        variable =
          v,

        CONUS_n =
          length(x1),

        CONUS_nonmissing =
          sum(
            is.finite(x1)
          ),

        CONUS_complete_pct =
          100 *
          mean(
            is.finite(x1)
          ),

        PRIMARY_n =
          length(x2),

        PRIMARY_nonmissing =
          sum(
            is.finite(x2)
          ),

        PRIMARY_complete_pct =
          100 *
          mean(
            is.finite(x2)
          )
      )
    }
  )
)

fwrite(
  COVERAGE,
  file.path(
    AUDIT_DIR,
    "Table1_variable_coverage.csv"
  )
)


# ============================================================
# 20. Final hard audit
# ============================================================

AUDIT <- data.table(

  check = c(
    "08_input_exists",
    "all_rows_769770",
    "all_monitors_181",
    "CONUS_rows_741861",
    "CONUS_monitors_168",
    "primary_rows_558915",
    "primary_monitors_140",
    "primary_Table1_variables_complete",
    "EN_CSV_exists",
    "CN_CSV_exists",
    "source_statistics_exists",
    "EN_HTML_exists",
    "CN_HTML_exists",
    "notes_exists"
  ),

  observed = c(

    as.character(
      as.integer(
        file.exists(
          INPUT08
        )
      )
    ),

    as.character(
      nrow(D)
    ),

    as.character(
      uniqueN(
        D$loc_num
      )
    ),

    as.character(
      nrow(CONUS)
    ),

    as.character(
      uniqueN(
        CONUS$loc_num
      )
    ),

    as.character(
      nrow(MAIN)
    ),

    as.character(
      uniqueN(
        MAIN$loc_num
      )
    ),

    as.character(
      sum(
        MAIN_MISSING$n_missing_or_nonfinite
      )
    ),

    as.character(
      as.integer(
        file.exists(
          OUT_EN_CSV
        )
      )
    ),

    as.character(
      as.integer(
        file.exists(
          OUT_CN_CSV
        )
      )
    ),

    as.character(
      as.integer(
        file.exists(
          OUT_SOURCE
        )
      )
    ),

    as.character(
      as.integer(
        file.exists(
          OUT_EN_HTML
        )
      )
    ),

    as.character(
      as.integer(
        file.exists(
          OUT_CN_HTML
        )
      )
    ),

    as.character(
      as.integer(
        file.exists(
          OUT_NOTE
        )
      )
    )
  ),

  expected = c(
    "1",
    as.character(
      EXPECTED_ALL_ROWS
    ),
    as.character(
      EXPECTED_ALL_MONITORS
    ),
    as.character(
      EXPECTED_CONUS_ROWS
    ),
    as.character(
      EXPECTED_CONUS_MONITORS
    ),
    as.character(
      EXPECTED_MAIN_ROWS
    ),
    as.character(
      EXPECTED_MAIN_MONITORS
    ),
    "0",
    "1",
    "1",
    "1",
    "1",
    "1",
    "1"
  )
)

AUDIT[
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
  AUDIT,
  OUT_AUDIT
)


# ============================================================
# 21. Terminal preview
# ============================================================

cat(
  "\n============================================================\n",
  "TABLE 1 ENGLISH PREVIEW\n",
  "============================================================\n",
  sep = ""
)

print(
  TABLE1_EN,
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "TABLE 1 中文预览\n",
  "============================================================\n",
  sep = ""
)

print(
  TABLE1_CN,
  nrows = Inf
)


# ============================================================
# 22. Final audit
# ============================================================

cat(
  "\n============================================================\n",
  "TABLE 1 FINAL AUDIT\n",
  "============================================================\n",
  sep = ""
)

print(
  AUDIT,
  nrows = Inf
)

N_FAIL <- AUDIT[
  status ==
    "FAIL",
  .N
]

cat(
  "\n============================================================\n",
  ifelse(
    N_FAIL ==
      0L,
    "TABLE 1 GENERATION: PASS\n",
    "TABLE 1 GENERATION: FAIL\n"
  ),
  "Hard FAIL = ",
  N_FAIL,
  "\n\n",
  "English CSV = ",
  OUT_EN_CSV,
  "\n",
  "Chinese CSV = ",
  OUT_CN_CSV,
  "\n",
  "Source statistics = ",
  OUT_SOURCE,
  "\n",
  "English HTML = ",
  OUT_EN_HTML,
  "\n",
  "Chinese HTML = ",
  OUT_CN_HTML,
  "\n",
  "Audit = ",
  OUT_AUDIT,
  "\n",
  "============================================================\n",
  sep = ""
)

if (
  N_FAIL >
    0L
) {
  stopf(
    "Table 1存在Hard FAIL，请勿作为论文正文表。"
  )
}

