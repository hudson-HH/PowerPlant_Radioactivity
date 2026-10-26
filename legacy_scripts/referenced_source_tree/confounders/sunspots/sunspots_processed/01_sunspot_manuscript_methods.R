suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 0. 路径
# ============================================================

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

SUN_ROOT <- file.path(
  ROOT,
  "混杂因素",
  "太阳黑子"
)

OUT_DIR <- file.path(
  SUN_ROOT,
  "太阳黑子处理版"
)

QUASI_FILE <- file.path(
  ROOT,
  "radnet 监测数据",
  "04_正式输出",
  "radnet_gross_beta_quasi_daily_model_ready_2001_2024.csv.gz"
)

CW_FILE <- file.path(
  ROOT,
  "radnet命名标准化",
  "radnet_181站_loc_num_radnet_monitor_id对照表.csv"
)

OUT_MONTHLY <- file.path(
  OUT_DIR,
  "01_sunspot_monthly_文章方法_2001_2024.csv"
)

OUT_QUASI <- file.path(
  OUT_DIR,
  "02_radnet_sunspot_quasi_daily_文章方法_2001_2024.csv.gz"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 1. 自动查找SILSO文件
# ============================================================

all_files <- list.files(
  SUN_ROOT,
  recursive = TRUE,
  full.names = TRUE
)

monthly_hits <- all_files[
  basename(all_files) == "SN_m_tot_V2.0.csv"
]

daily_hits <- all_files[
  basename(all_files) == "SN_d_tot_V2.0.csv"
]

if (length(monthly_hits) > 1) {
  stop("发现多个SN_m_tot_V2.0.csv，请只保留一个正式原始文件。")
}

if (length(daily_hits) > 1) {
  stop("发现多个SN_d_tot_V2.0.csv，请只保留一个正式原始文件。")
}

cat(
  "\n============================================================\n",
  "太阳黑子 —— 文章月度方法\n",
  "============================================================\n",
  sep = ""
)

# ============================================================
# 2A. 优先读取SILSO官方月度数据
# ============================================================

if (length(monthly_hits) == 1) {

  source_file <- monthly_hits[1]

  cat(
    "使用SILSO官方月度文件：\n",
    source_file,
    "\n",
    sep = ""
  )

  x <- fread(
    source_file,
    sep = ";",
    header = FALSE,
    fill = TRUE,
    strip.white = TRUE
  )

  if (ncol(x) < 4) {
    stop("SN_m_tot_V2.0.csv字段数异常。")
  }

  monthly <- x[
    ,
    .(
      Year = as.integer(V1),
      Month = as.integer(V2),
      Num = as.numeric(V4)
    )
  ]

  monthly[
    Num < 0,
    Num := NA_real_
  ]

  source_method <- "SILSO_official_monthly"

# ============================================================
# 2B. 如果没有月度文件，用官方日数据计算月平均
# ============================================================

} else {

  if (length(daily_hits) != 1) {
    stop(
      "既没有SN_m_tot_V2.0.csv，也没有唯一的SN_d_tot_V2.0.csv。"
    )
  }

  source_file <- daily_hits[1]

  cat(
    "未找到SILSO官方月度文件。\n",
    "使用官方日数据按年月计算月平均：\n",
    source_file,
    "\n",
    sep = ""
  )

  x <- fread(
    source_file,
    sep = ";",
    header = FALSE,
    fill = TRUE,
    strip.white = TRUE
  )

  if (ncol(x) < 5) {
    stop("SN_d_tot_V2.0.csv字段数异常。")
  }

  daily <- x[
    ,
    .(
      Year = as.integer(V1),
      Month = as.integer(V2),
      Day = as.integer(V3),
      daily_sunspot = as.numeric(V5)
    )
  ]

  # SILSO中负值表示缺失，不当作真实太阳黑子数
  daily[
    daily_sunspot < 0,
    daily_sunspot := NA_real_
  ]

  monthly <- daily[
    !is.na(daily_sunspot),
    .(
      Num = mean(daily_sunspot),
      n_valid_days = .N
    ),
    by = .(
      Year,
      Month
    )
  ]

  source_method <- "SILSO_daily_to_monthly_mean"
}

# ============================================================
# 3. 限定研究期2001-2024
# ============================================================

monthly <- monthly[
  Year >= 2001 &
  Year <= 2024
]

setorder(
  monthly,
  Year,
  Month
)

# ============================================================
# 4. 结构检查
#
# 2001-2024 = 24年 × 12月 = 288个月
# ============================================================

expected_months <- CJ(
  Year = 2001:2024,
  Month = 1:12
)

check_months <- merge(
  expected_months,
  monthly[
    ,
    .(
      Year,
      Month,
      Num
    )
  ],
  by = c(
    "Year",
    "Month"
  ),
  all.x = TRUE
)

if (nrow(monthly) != 288) {

  cat(
    "实际月数 =",
    nrow(monthly),
    "\n"
  )

  stop(
    "2001-2024太阳黑子数据不是完整288个月。"
  )
}

if (
  monthly[
    ,
    .N,
    by = .(
      Year,
      Month
    )
  ][
    N > 1,
    .N
  ] > 0
) {

  stop(
    "太阳黑子存在重复Year+Month。"
  )
}

if (anyNA(check_months$Num)) {

  cat(
    "\n缺失月份：\n"
  )

  print(
    check_months[
      is.na(Num)
    ]
  )

  stop(
    "2001-2024存在缺失月度太阳黑子数。"
  )
}

monthly[
  ,
  sunspot_monthly :=
    Num
]

monthly[
  ,
  source :=
    source_method
]

# ============================================================
# 5. 输出正式月度太阳黑子表
# ============================================================

fwrite(
  monthly,
  OUT_MONTHLY
)

cat(
  "\n================ 月度太阳黑子数据 ================\n"
)

cat(
  "月份数 =",
  nrow(monthly),
  "\n"
)

cat(
  "日期范围 = 2001-01 至 2024-12\n"
)

cat(
  "缺失月份 =",
  monthly[
    is.na(Num),
    .N
  ],
  "\n"
)

cat(
  "\nNum分布：\n"
)

print(
  summary(
    monthly$Num
  )
)

# ============================================================
# 6. 读取正式RadNet quasi-daily
# ============================================================

if (!file.exists(QUASI_FILE)) {
  stop("找不到正式RadNet quasi-daily文件。")
}

quasi <- fread(
  QUASI_FILE
)

required_q <- c(
  "loc_num",
  "date",
  "collect_start",
  "collect_end"
)

missing_q <- setdiff(
  required_q,
  names(quasi)
)

if (length(missing_q) > 0) {

  stop(
    "quasi-daily缺少字段：",
    paste(
      missing_q,
      collapse = ", "
    )
  )
}

quasi[
  ,
  loc_num :=
    as.integer(loc_num)
]

quasi[
  ,
  date :=
    as.IDate(date)
]

quasi[
  ,
  collect_start :=
    as.IDate(collect_start)
]

quasi[
  ,
  collect_end :=
    as.IDate(collect_end)
]

quasi[
  ,
  Year :=
    as.integer(
      format(
        date,
        "%Y"
      )
    )
]

quasi[
  ,
  Month :=
    as.integer(
      format(
        date,
        "%m"
      )
    )
]

# ============================================================
# 7. 按作者方法：year + month匹配
#
# 作者代码：
#
# dataset$Month = month(dataset$Date)
# dataset <- left_join(
#   sun_spots,
#   by=c("year"="Year","Month"="Month")
# )
#
# Num即最终模型中的太阳黑子变量
# ============================================================

sun_for_join <- monthly[
  ,
  .(
    Year,
    Month,
    Num,
    sunspot_monthly
  )
]

quasi_sun <- merge(
  quasi[
    ,
    .(
      loc_num,
      date,
      collect_start,
      collect_end,
      Year,
      Month
    )
  ],
  sun_for_join,
  by = c(
    "Year",
    "Month"
  ),
  all.x = TRUE,
  sort = FALSE
)

if (nrow(quasi_sun) != nrow(quasi)) {

  stop(
    "太阳黑子join改变了quasi-daily行数：",
    nrow(quasi),
    " -> ",
    nrow(quasi_sun)
  )
}

# ============================================================
# 8. 加181站信息及CONUS标记
# ============================================================

cw <- fread(
  CW_FILE
)

cw[
  ,
  loc_num :=
    as.integer(loc_num)
]

CONUS_STATES <- c(
  setdiff(
    state.abb,
    c(
      "AK",
      "HI"
    )
  ),
  "DC"
)

cw[
  ,
  main_analysis_conus :=
    state_abbr %in% CONUS_STATES
]

if (
  cw[
    main_analysis_conus == TRUE,
    .N
  ] != 168
) {

  stop(
    "CONUS站数不是168。"
  )
}

station_info <- cw[
  ,
  .(
    loc_num,
    radnet_monitor_id,
    station,
    state_abbr,
    main_analysis_conus
  )
]

quasi_sun <- merge(
  quasi_sun,
  station_info,
  by = "loc_num",
  all.x = TRUE,
  sort = FALSE
)

if (nrow(quasi_sun) != nrow(quasi)) {

  stop(
    "站点信息join改变了quasi-daily行数。"
  )
}

setcolorder(
  quasi_sun,
  c(
    "loc_num",
    "radnet_monitor_id",
    "station",
    "state_abbr",
    "main_analysis_conus",
    "date",
    "collect_start",
    "collect_end",
    "Year",
    "Month",
    "Num",
    "sunspot_monthly"
  )
)

setorder(
  quasi_sun,
  loc_num,
  collect_start,
  collect_end,
  date
)

# ============================================================
# 9. 输出正式quasi-daily太阳黑子表
# ============================================================

fwrite(
  quasi_sun,
  OUT_QUASI
)

# ============================================================
# 10. 最终覆盖检查
# ============================================================

conus <- quasi_sun[
  main_analysis_conus == TRUE
]

cat(
  "\n============================================================\n"
)

cat(
  "太阳黑子文章方法结果\n"
)

cat(
  "============================================================\n"
)

cat(
  "181站总数 =",
  uniqueN(quasi_sun$loc_num),
  "\n"
)

cat(
  "CONUS站数 =",
  uniqueN(conus$loc_num),
  "\n"
)

cat(
  "全部quasi-daily行数 =",
  nrow(quasi_sun),
  "\n"
)

cat(
  "CONUS quasi-daily行数 =",
  nrow(conus),
  "\n"
)

cat(
  "CONUS有太阳黑子值行数 =",
  conus[
    !is.na(Num),
    .N
  ],
  "\n"
)

cat(
  "CONUS太阳黑子缺失行数 =",
  conus[
    is.na(Num),
    .N
  ],
  "\n"
)

cat(
  "CONUS太阳黑子完整率 =",
  sprintf(
    "%.3f%%",
    100 *
    mean(
      !is.na(
        conus$Num
      )
    )
  ),
  "\n"
)

if (
  conus[
    is.na(Num),
    .N
  ] != 0
) {

  stop(
    "CONUS存在太阳黑子缺失，暂不能冻结。"
  )
}

cat(
  "\n============================================================\n"
)

cat(
  "太阳黑子文章方法重建：PASS\n"
)

cat(
  "============================================================\n"
)

cat(
  "\n正式月度太阳黑子表：\n",
  OUT_MONTHLY,
  "\n",
  sep = ""
)

cat(
  "\n正式quasi-daily太阳黑子表：\n",
  OUT_QUASI,
  "\n",
  sep = ""
)

cat(
  "\n运行结束。\n"
)
