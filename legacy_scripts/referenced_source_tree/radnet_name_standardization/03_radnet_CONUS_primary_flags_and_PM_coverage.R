# ============================================================
# RadNet CONUS 主分析站点定义 + PM2.5覆盖与模型行损失评估
#
# 固定原则：
# 1. 不修改原181站crosswalk
# 2. 主分析地域 = contiguous 48 states + DC
# 3. AK / HI / PR / GU / CNMI / Canada 等不属于CONUS主分析
# 4. PM2.5要求RadNet周围50 km内存在PM站
# 5. 不进行最近站替代、扩大半径或PM均值填补
# 6. quasi-daily重复monitor-date若来自不同filter period，不去重
#    因为最终模型单位仍然是quasi-daily observation row
# 7. 只生成一个正式181站标记表；其他结果打印到终端
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})


# ============================================================
# 1. Paths
# ============================================================

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

CW_FILE <- file.path(
  ROOT,
  "radnet命名标准化",
  "radnet_181站_loc_num_radnet_monitor_id对照表.csv"
)

PM_LINK_FILE <- file.path(
  ROOT,
  "混杂因素",
  "PM2.5",
  "PM2.5处理版",
  "01_radnet_pm25_site_links_50km_文章方法.csv"
)

PM_DAILY_FILE <- file.path(
  ROOT,
  "混杂因素",
  "PM2.5",
  "PM2.5处理版",
  "02_radnet_pm25_daily_50km_文章方法_2001_2024.csv.gz"
)

RADNET_DAILY_FILE <- file.path(
  ROOT,
  "radnet 监测数据",
  "04_正式输出",
  "radnet_gross_beta_quasi_daily_model_ready_2001_2024.csv.gz"
)

OUT_FILE <- file.path(
  ROOT,
  "radnet命名标准化",
  "radnet_181站_CONUS主分析标记与PM覆盖.csv"
)


# ============================================================
# 2. Input checks
# ============================================================

input_files <- c(
  CW_FILE,
  PM_LINK_FILE,
  PM_DAILY_FILE,
  RADNET_DAILY_FILE
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {

  stop(
    "以下输入文件不存在：\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


# ============================================================
# 3. Read 181-station crosswalk
# ============================================================

cw <- fread(CW_FILE)

required_cw <- c(
  "loc_num",
  "radnet_monitor_id",
  "station",
  "city_name",
  "state_abbr",
  "monitor_latitude",
  "monitor_longitude"
)

missing_cw <- setdiff(
  required_cw,
  names(cw)
)

if (length(missing_cw) > 0) {

  stop(
    "crosswalk缺少字段：",
    paste(
      missing_cw,
      collapse = ", "
    )
  )
}

cw[, loc_num := as.integer(loc_num)]

if (nrow(cw) != 181) {
  stop(
    "crosswalk行数不是181：",
    nrow(cw)
  )
}

if (uniqueN(cw$loc_num) != 181) {
  stop(
    "crosswalk的loc_num不是181个唯一值。"
  )
}


# ============================================================
# 4. Define CONUS
#
# contiguous 48 states + DC
# ============================================================

CONUS_STATE_ABBR <- c(
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
    state_abbr %in% CONUS_STATE_ABBR
]

cw[
  ,
  main_analysis_scope :=
    fifelse(
      main_analysis_conus,
      "CONUS_main_analysis",
      "outside_CONUS_main_analysis"
    )
]


# ============================================================
# 5. Read PM 50-km link table
# ============================================================

pm_link <- fread(PM_LINK_FILE)

required_link <- c(
  "loc_num",
  "pm25_site_id"
)

missing_link <- setdiff(
  required_link,
  names(pm_link)
)

if (length(missing_link) > 0) {

  stop(
    "PM link表缺少字段：",
    paste(
      missing_link,
      collapse = ", "
    )
  )
}

pm_link[, loc_num := as.integer(loc_num)]

# 防止同一个RadNet-PM pair重复
pm_link_unique <- unique(
  pm_link[
    ,
    .(
      loc_num,
      pm25_site_id
    )
  ]
)

pm_station_summary <- pm_link_unique[
  ,
  .(
    pm50_n_sites =
      uniqueN(
        pm25_site_id
      )
  ),
  by = loc_num
]

cw <- merge(
  cw,
  pm_station_summary,
  by = "loc_num",
  all.x = TRUE,
  sort = FALSE
)

cw[
  is.na(pm50_n_sites),
  pm50_n_sites := 0L
]

cw[
  ,
  pm50_has_link :=
    pm50_n_sites > 0
]


# ============================================================
# 6. Read RadNet quasi-daily table
# ============================================================

rad <- fread(
  RADNET_DAILY_FILE
)

if (!"loc_num" %in% names(rad)) {
  stop(
    "RadNet quasi-daily表没有loc_num字段。"
  )
}

rad[, loc_num := as.integer(loc_num)]

date_candidates <- c(
  "date",
  "Date",
  "date_local"
)

date_col <- date_candidates[
  date_candidates %in% names(rad)
][1]

if (is.na(date_col)) {

  cat(
    "\nRadNet quasi-daily字段：\n"
  )

  print(
    names(rad)
  )

  stop(
    "无法找到quasi-daily日期字段。"
  )
}

if (date_col != "date") {

  setnames(
    rad,
    date_col,
    "date"
  )
}

rad[, date := as.IDate(date)]

# 再限制一次正式研究期
rad <- rad[
  date >= as.IDate("2001-01-01") &
  date <= as.IDate("2024-12-31")
]


# ============================================================
# 7. Validate quasi-daily station IDs
# ============================================================

unknown_rad_ids <- setdiff(
  unique(rad$loc_num),
  cw$loc_num
)

if (length(unknown_rad_ids) > 0) {

  stop(
    "quasi-daily中存在crosswalk没有的loc_num：",
    paste(
      unknown_rad_ids,
      collapse = ", "
    )
  )
}


# ============================================================
# 8. Per-station RadNet quasi-daily counts
#
# n_quasi_rows:
#   真正进入模型的数据行数口径
#
# n_unique_dates:
#   仅用于了解重复monitor-date情况
# ============================================================

rad_station_summary <- rad[
  ,
  .(
    quasi_rows_2001_2024 = .N,

    quasi_unique_dates_2001_2024 =
      uniqueN(
        date
      ),

    quasi_first_date =
      min(
        date,
        na.rm = TRUE
      ),

    quasi_last_date =
      max(
        date,
        na.rm = TRUE
      )
  ),
  by = loc_num
]

cw <- merge(
  cw,
  rad_station_summary,
  by = "loc_num",
  all.x = TRUE,
  sort = FALSE
)

cw[
  is.na(quasi_rows_2001_2024),
  quasi_rows_2001_2024 := 0L
]

cw[
  is.na(quasi_unique_dates_2001_2024),
  quasi_unique_dates_2001_2024 := 0L
]


# ============================================================
# 9. Read formal daily PM2.5
# ============================================================

pm_daily <- fread(
  PM_DAILY_FILE,
  select = c(
    "loc_num",
    "date",
    "pm_c"
  )
)

pm_daily[, loc_num := as.integer(loc_num)]
pm_daily[, date := as.IDate(date)]

if (
  pm_daily[
    ,
    anyDuplicated(
      paste(
        loc_num,
        date
      )
    )
  ] > 0
) {

  stop(
    "PM daily不是唯一loc_num + date。"
  )
}

if (anyNA(pm_daily$pm_c)) {

  stop(
    "正式PM daily中pm_c存在NA。"
  )
}


# ============================================================
# 10. Restrict quasi-daily to CONUS
# ============================================================

conus_ids <- cw[
  main_analysis_conus == TRUE,
  loc_num
]

rad_conus <- rad[
  loc_num %in% conus_ids
]


# ============================================================
# 11. Add PM spatial-link status to every CONUS quasi row
# ============================================================

station_pm_flag <- cw[
  main_analysis_conus == TRUE,
  .(
    loc_num,
    pm50_has_link,
    pm50_n_sites
  )
]

rad_conus_eval <- merge(
  rad_conus,
  station_pm_flag,
  by = "loc_num",
  all.x = TRUE,
  sort = FALSE
)


# ============================================================
# 12. Join daily PM
#
# PM table unique loc_num+date，
# 因此若RadNet同一天有不同filter-period observation，
# 每一行继承同一天PM，不进行RadNet去重。
# ============================================================

rad_conus_eval <- merge(
  rad_conus_eval,
  pm_daily,
  by = c(
    "loc_num",
    "date"
  ),
  all.x = TRUE,
  sort = FALSE
)


# ============================================================
# 13. Classify PM missingness
# ============================================================

rad_conus_eval[
  ,
  pm_status :=
    fifelse(
      !pm50_has_link,

      "01_no_PM_station_within_50km",

      fifelse(
        is.na(pm_c),

        "02_PM_station_exists_but_no_daily_measurement",

        "00_PM_available"
      )
    )
]


# ============================================================
# 14. Per-station actual PM availability on RadNet quasi-days
# ============================================================

station_model_summary <- rad_conus_eval[
  ,
  .(
    quasi_rows_pm_available =
      sum(
        pm_status ==
          "00_PM_available"
      ),

    quasi_rows_pm_missing_total =
      sum(
        pm_status !=
          "00_PM_available"
      ),

    quasi_rows_pm_missing_no50km =
      sum(
        pm_status ==
          "01_no_PM_station_within_50km"
      ),

    quasi_rows_pm_missing_daily =
      sum(
        pm_status ==
          "02_PM_station_exists_but_no_daily_measurement"
      )
  ),
  by = loc_num
]

cw <- merge(
  cw,
  station_model_summary,
  by = "loc_num",
  all.x = TRUE,
  sort = FALSE
)

count_fields <- c(
  "quasi_rows_pm_available",
  "quasi_rows_pm_missing_total",
  "quasi_rows_pm_missing_no50km",
  "quasi_rows_pm_missing_daily"
)

for (f in count_fields) {

  set(
    cw,
    which(
      is.na(
        cw[[f]]
      )
    ),
    f,
    0L
  )
}


# ============================================================
# 15. Restore canonical ordering
# ============================================================

setorder(
  cw,
  loc_num
)


# ============================================================
# 16. Core station counts
# ============================================================

n_all <- nrow(cw)

n_conus <- cw[
  main_analysis_conus == TRUE,
  uniqueN(loc_num)
]

n_nonconus <- cw[
  main_analysis_conus == FALSE,
  uniqueN(loc_num)
]

n_conus_pm <- cw[
  main_analysis_conus == TRUE &
  pm50_has_link == TRUE,
  uniqueN(loc_num)
]

n_conus_no_pm <- cw[
  main_analysis_conus == TRUE &
  pm50_has_link == FALSE,
  uniqueN(loc_num)
]


# ============================================================
# 17. Structural no-PM stations
# ============================================================

no_pm_conus <- cw[
  main_analysis_conus == TRUE &
  pm50_has_link == FALSE,
  .(
    loc_num,
    radnet_monitor_id,
    station,
    state_abbr,
    pm50_n_sites,
    quasi_rows_2001_2024,
    quasi_unique_dates_2001_2024,
    duplicated_quasi_date_rows =
      quasi_rows_2001_2024 -
      quasi_unique_dates_2001_2024,
    quasi_first_date,
    quasi_last_date
  )
]

setorder(
  no_pm_conus,
  state_abbr,
  station
)


# ============================================================
# 18. PM-related row losses
# ============================================================

pm_loss <- rad_conus_eval[
  ,
  .(
    rows = .N
  ),
  by = pm_status
]

total_conus_rows <- nrow(
  rad_conus_eval
)

pm_available_rows <- rad_conus_eval[
  pm_status ==
    "00_PM_available",
  .N
]

pm_structural_loss_rows <- rad_conus_eval[
  pm_status ==
    "01_no_PM_station_within_50km",
  .N
]

pm_daily_loss_rows <- rad_conus_eval[
  pm_status ==
    "02_PM_station_exists_but_no_daily_measurement",
  .N
]

pm_total_missing_rows <-
  pm_structural_loss_rows +
  pm_daily_loss_rows


# ============================================================
# 19. Station-level overlap
# ============================================================

stations_with_any_pm_on_radnet_days <- rad_conus_eval[
  pm_status ==
    "00_PM_available",
  uniqueN(loc_num)
]

linked_but_zero_overlap <- rad_conus_eval[
  pm50_has_link == TRUE,
  .(
    n_pm_available =
      sum(
        !is.na(pm_c)
      )
  ),
  by = loc_num
][
  n_pm_available == 0
]


# ============================================================
# 20. Validation against current known station universe
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "RadNet CONUS主分析站点定义\n"
)

cat(
  "============================================================\n"
)

cat(
  "181站总数 =",
  n_all,
  "\n"
)

cat(
  "CONUS loc_num =",
  n_conus,
  "\n"
)

cat(
  "非CONUS loc_num =",
  n_nonconus,
  "\n"
)

cat(
  "CONUS中有>=1个50 km PM站 =",
  n_conus_pm,
  "\n"
)

cat(
  "CONUS中没有50 km PM站 =",
  n_conus_no_pm,
  "\n"
)


# 已知crosswalk下应为168
if (n_conus != 168) {

  stop(
    "CONUS数量不是预期的168，请检查state_abbr。"
  )
}

# 当前PM链接结果应为163+5
if (
  n_conus_pm != 163 ||
  n_conus_no_pm != 5
) {

  stop(
    "CONUS PM覆盖不是预期的163有PM + 5无PM，请检查PM链接。"
  )
}


# ============================================================
# 21. Print non-CONUS
# ============================================================

cat(
  "\n================ 非CONUS的13站 ================\n"
)

print(
  cw[
    main_analysis_conus == FALSE,
    .(
      loc_num,
      radnet_monitor_id,
      station,
      state_abbr,
      quasi_rows_2001_2024
    )
  ][
    order(
      state_abbr,
      station
    )
  ]
)


# ============================================================
# 22. Print five structural no-PM CONUS stations
# ============================================================

cat(
  "\n================ CONUS中无50 km PM站的5站 ================\n"
)

print(
  no_pm_conus
)


# ============================================================
# 23. Print contribution of five stations
# ============================================================

five_rows <- sum(
  no_pm_conus$quasi_rows_2001_2024
)

cat(
  "\n================ 5个无PM站对quasi-daily的贡献 ================\n"
)

cat(
  "5站quasi-daily总行数 =",
  five_rows,
  "\n"
)

cat(
  "CONUS quasi-daily总行数 =",
  total_conus_rows,
  "\n"
)

cat(
  "仅因5站无50km PM造成的结构性行损失比例 =",
  sprintf(
    "%.3f%%",
    100 *
      five_rows /
      total_conus_rows
  ),
  "\n"
)


# ============================================================
# 24. Actual PM loss on RadNet quasi-days
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "PM2.5对CONUS quasi-daily模型行的实际影响\n"
)

cat(
  "============================================================\n"
)

print(
  pm_loss[
    order(
      pm_status
    )
  ]
)

cat(
  "\nCONUS quasi-daily总行数 =",
  total_conus_rows,
  "\n"
)

cat(
  "有PM2.5的行数 =",
  pm_available_rows,
  "\n"
)

cat(
  "因完全没有50km PM站而缺失 =",
  pm_structural_loss_rows,
  "\n"
)

cat(
  "有50km PM站但当天无PM测量而缺失 =",
  pm_daily_loss_rows,
  "\n"
)

cat(
  "PM2.5总缺失行数 =",
  pm_total_missing_rows,
  "\n"
)

cat(
  "仅考虑PM后保留比例 =",
  sprintf(
    "%.3f%%",
    100 *
      pm_available_rows /
      total_conus_rows
  ),
  "\n"
)

cat(
  "仅考虑PM后损失比例 =",
  sprintf(
    "%.3f%%",
    100 *
      pm_total_missing_rows /
      total_conus_rows
  ),
  "\n"
)

cat(
  "在实际RadNet quasi-days上至少有1天PM数据的CONUS loc_num =",
  stations_with_any_pm_on_radnet_days,
  "\n"
)


# ============================================================
# 25. Linked station but zero overlap check
# ============================================================

cat(
  "\n================ 有50km PM站但与RadNet日期完全无重叠 ================\n"
)

cat(
  "站数 =",
  nrow(
    linked_but_zero_overlap
  ),
  "\n"
)

if (nrow(linked_but_zero_overlap) > 0) {

  print(
    merge(
      linked_but_zero_overlap,
      cw[
        ,
        .(
          loc_num,
          radnet_monitor_id,
          station,
          state_abbr
        )
      ],
      by = "loc_num",
      all.x = TRUE
    )
  )
}


# ============================================================
# 26. PM missingness among linked stations only
# ============================================================

cat(
  "\n================ 已有50km PM站的163站：日期级PM缺失 ================\n"
)

linked_rows <- rad_conus_eval[
  pm50_has_link == TRUE,
  .N
]

linked_pm_missing_rows <- rad_conus_eval[
  pm50_has_link == TRUE &
  is.na(pm_c),
  .N
]

cat(
  "163个有PM链接站的quasi-daily行数 =",
  linked_rows,
  "\n"
)

cat(
  "其中当天PM缺失行数 =",
  linked_pm_missing_rows,
  "\n"
)

cat(
  "日期级PM缺失率 =",
  sprintf(
    "%.3f%%",
    100 *
      linked_pm_missing_rows /
      linked_rows
  ),
  "\n"
)


# ============================================================
# 27. Write ONE formal station table
# ============================================================

formal_cols <- c(
  "loc_num",
  "radnet_monitor_id",
  "station",
  "city_name",
  "state_abbr",
  "monitor_latitude",
  "monitor_longitude",

  "main_analysis_conus",
  "main_analysis_scope",

  "pm50_has_link",
  "pm50_n_sites",

  "quasi_rows_2001_2024",
  "quasi_unique_dates_2001_2024",
  "quasi_first_date",
  "quasi_last_date",

  "quasi_rows_pm_available",
  "quasi_rows_pm_missing_total",
  "quasi_rows_pm_missing_no50km",
  "quasi_rows_pm_missing_daily"
)

fwrite(
  cw[
    ,
    ..formal_cols
  ],
  OUT_FILE
)


# ============================================================
# 28. Final integrity checks
# ============================================================

out_check <- fread(
  OUT_FILE
)

if (nrow(out_check) != 181) {
  stop(
    "正式输出不是181行。"
  )
}

if (uniqueN(out_check$loc_num) != 181) {
  stop(
    "正式输出loc_num不是181个唯一值。"
  )
}

if (
  out_check[
    main_analysis_conus == TRUE,
    uniqueN(loc_num)
  ] != 168
) {

  stop(
    "正式输出CONUS不是168站。"
  )
}


# ============================================================
# 29. PASS
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "RadNet CONUS主分析定义与PM覆盖评估完成\n"
)

cat(
  "============================================================\n"
)

cat(
  "正式181站标记表：\n",
  OUT_FILE,
  "\n",
  sep = ""
)

cat(
  "\n总体状态：PASS\n"
)
