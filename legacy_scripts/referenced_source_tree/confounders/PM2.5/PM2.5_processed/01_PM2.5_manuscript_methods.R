# ============================================================
# PM2.5 —— 按参考文章最终方法重做
#
# 核心定义
# ------------------------------------------------------------
# 1. RadNet与EPA PM2.5监测站空间链接半径固定为 50 km
# 2. 50 km内使用全部PM2.5站，不选最近站
# 3. 不扩大半径
# 4. 不进行PM2.5空间填补或均值填补
# 5. RadNet时间单位首先是filter period：
#
#       [collect_start, collect_end)
#
#    即 collect_end 当天不属于该滤膜周期。
#
# 6. 对每一个RadNet filter period：
#      - 找50 km内全部PM站
#      - 提取该filter period内全部有效PM site-day
#      - 对全部匹配到的PM site-day均值求平均
#      - 得到该filter period唯一 pm_c
#
# 7. 再将该filter-period pm_c赋给该周期展开出的每一个
#    quasi-daily observation。
#
# 8. 同一个loc_num+date如果来自不同filter period，
#    不去重；正式观测键是：
#
#       loc_num + collect_start + collect_end + date
#
# 9. 181站底表全部保留；
#    主分析CONUS = contiguous 48 states + DC。
#
# 10. 不使用旧100 km版本。
# ============================================================


suppressPackageStartupMessages({
  library(data.table)
})


# ============================================================
# 0. Paths
# ============================================================

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

CW_FILE <- file.path(
  ROOT,
  "radnet命名标准化",
  "radnet_181站_loc_num_radnet_monitor_id对照表.csv"
)

FILTER_FILE <- file.path(
  ROOT,
  "radnet 监测数据",
  "02_滤膜周期数据",
  "radnet_gross_beta_filter_period_2001_2024.csv.gz"
)

QUASI_FILE <- file.path(
  ROOT,
  "radnet 监测数据",
  "04_正式输出",
  "radnet_gross_beta_quasi_daily_model_ready_2001_2024.csv.gz"
)

PM_RAW_FILE <- file.path(
  ROOT,
  "混杂因素",
  "PM2.5",
  "outputs",
  "pm25_site_day_2001_2024.csv"
)

OUT_DIR <- file.path(
  ROOT,
  "混杂因素",
  "PM2.5",
  "PM2.5处理版"
)

OUT_LINK <- file.path(
  OUT_DIR,
  "01_radnet_pm25_site_links_50km_文章方法.csv"
)

OUT_PERIOD <- file.path(
  OUT_DIR,
  "02_radnet_pm25_filter_period_50km_文章方法_2001_2024.csv.gz"
)

OUT_QUASI <- file.path(
  OUT_DIR,
  "03_radnet_pm25_quasi_daily_50km_文章方法_2001_2024.csv.gz"
)

OUT_STATION <- file.path(
  ROOT,
  "radnet命名标准化",
  "radnet_181站_CONUS主分析标记与PM覆盖.csv"
)

OLD_DAILY_FILE <- file.path(
  OUT_DIR,
  "02_radnet_pm25_daily_50km_文章方法_2001_2024.csv.gz"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 1. Input checks
# ============================================================

required_files <- c(
  CW_FILE,
  FILTER_FILE,
  QUASI_FILE,
  PM_RAW_FILE
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  stop(
    "缺少输入文件：\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}


# ============================================================
# 2. Read 181-station RadNet master
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
    "crosswalk不是181行：",
    nrow(cw)
  )
}

if (uniqueN(cw$loc_num) != 181) {
  stop(
    "crosswalk loc_num不是181个唯一值。"
  )
}


# ============================================================
# 3. CONUS main-analysis flag
#
# contiguous 48 states + DC
# ============================================================

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
# 4. Read PM2.5 site-day source
# ============================================================

cat(
  "\n================ 读取PM2.5 site-day ================\n"
)

pm <- fread(PM_RAW_FILE)

required_pm <- c(
  "pm25_site_id",
  "date_local",
  "pm25_daily_mean",
  "n_pm25_valid",
  "latitude",
  "longitude"
)

missing_pm <- setdiff(
  required_pm,
  names(pm)
)

if (length(missing_pm) > 0) {

  stop(
    "PM2.5 site-day缺少字段：",
    paste(
      missing_pm,
      collapse = ", "
    )
  )
}

pm[
  ,
  pm25_site_id :=
    as.character(pm25_site_id)
]

pm[
  ,
  date_local :=
    as.IDate(date_local)
]

pm <- pm[
  date_local >= as.IDate("2001-01-01") &
  date_local <= as.IDate("2024-12-31")
]

cat(
  "PM site-day行数 =",
  nrow(pm),
  "\n"
)

cat(
  "PM站数 =",
  uniqueN(pm$pm25_site_id),
  "\n"
)

cat(
  "日期范围 =",
  as.character(min(pm$date_local)),
  "至",
  as.character(max(pm$date_local)),
  "\n"
)


# ============================================================
# 5. Validate one PM site-day row per site-date
# ============================================================

dup_pm <- pm[
  ,
  .N,
  by = .(
    pm25_site_id,
    date_local
  )
][
  N > 1
]

if (nrow(dup_pm) > 0) {

  stop(
    "PM输入不是唯一pm25_site_id + date_local；",
    "重复组合数=",
    nrow(dup_pm)
  )
}


# ============================================================
# 6. Representative PM station coordinates
#
# 原输入每天都带站点坐标。
# 对同一个AQS站用全部site-day坐标中位数作为站点坐标。
# ============================================================

pm_sites <- pm[
  is.finite(latitude) &
  is.finite(longitude),
  .(
    pm_site_latitude =
      median(
        latitude,
        na.rm = TRUE
      ),

    pm_site_longitude =
      median(
        longitude,
        na.rm = TRUE
      )
  ),
  by = pm25_site_id
]

cat(
  "有有效坐标的PM站数 =",
  nrow(pm_sites),
  "\n"
)


# ============================================================
# 7. Haversine distance
# ============================================================

haversine_km <- function(
  lat1,
  lon1,
  lat2,
  lon2
) {

  R <- 6371.0088
  rad <- pi / 180

  lat1 <- lat1 * rad
  lon1 <- lon1 * rad
  lat2 <- lat2 * rad
  lon2 <- lon2 * rad

  dlat <- lat2 - lat1
  dlon <- lon2 - lon1

  a <-
    sin(dlat / 2)^2 +
    cos(lat1) *
    cos(lat2) *
    sin(dlon / 2)^2

  a <- pmin(
    1,
    pmax(
      0,
      a
    )
  )

  2 * R * asin(
    sqrt(a)
  )
}


# ============================================================
# 8. Rebuild RadNet–PM 50 km links
#
# 与作者方法一致：
# 所有 <=50 km PM站全部进入link table。
# ============================================================

cat(
  "\n================ 建立RadNet-PM 50 km链接 ================\n"
)

link_list <- vector(
  "list",
  nrow(cw)
)

for (i in seq_len(nrow(cw))) {

  rr <- cw[i]

  d <- haversine_km(
    rr$monitor_latitude,
    rr$monitor_longitude,
    pm_sites$pm_site_latitude,
    pm_sites$pm_site_longitude
  )

  keep <- which(
    is.finite(d) &
    d <= 50
  )

  if (length(keep) == 0) {

    link_list[[i]] <- NULL
    next
  }

  x <- pm_sites[keep]

  x[
    ,
    `:=`(
      loc_num =
        rr$loc_num,

      radnet_monitor_id =
        rr$radnet_monitor_id,

      station =
        rr$station,

      state_abbr =
        rr$state_abbr,

      monitor_latitude =
        rr$monitor_latitude,

      monitor_longitude =
        rr$monitor_longitude,

      distance_km =
        d[keep]
    )
  ]

  link_list[[i]] <- x[
    ,
    .(
      loc_num,
      radnet_monitor_id,
      station,
      state_abbr,
      monitor_latitude,
      monitor_longitude,

      pm25_site_id,
      pm_site_latitude,
      pm_site_longitude,

      distance_km
    )
  ]
}

links <- rbindlist(
  link_list,
  use.names = TRUE,
  fill = TRUE
)

setorder(
  links,
  loc_num,
  distance_km,
  pm25_site_id
)

if (
  links[
    ,
    anyDuplicated(
      paste(
        loc_num,
        pm25_site_id
      )
    )
  ] > 0
) {

  stop(
    "RadNet-PM link表存在重复pair。"
  )
}

fwrite(
  links,
  OUT_LINK
)

cat(
  "50 km RadNet-PM pair数 =",
  nrow(links),
  "\n"
)

cat(
  "至少有1个50 km PM站的RadNet loc_num =",
  uniqueN(links$loc_num),
  "\n"
)


# ============================================================
# 9. Per-station spatial PM status
# ============================================================

link_station <- links[
  ,
  .(
    pm50_n_sites =
      uniqueN(pm25_site_id)
  ),
  by = loc_num
]

station_info <- merge(
  cw,
  link_station,
  by = "loc_num",
  all.x = TRUE,
  sort = FALSE
)

station_info[
  is.na(pm50_n_sites),
  pm50_n_sites := 0L
]

station_info[
  ,
  pm50_has_link :=
    pm50_n_sites > 0
]


# ============================================================
# 10. Valid PM observations
#
# 保留0值。
# 不做PM>0截断。
# 仅要求：
#   pm25_daily_mean finite
#   n_pm25_valid > 0
# ============================================================

pm_valid <- pm[
  is.finite(pm25_daily_mean) &
  !is.na(n_pm25_valid) &
  n_pm25_valid > 0,
  .(
    pm25_site_id,
    date_local,
    pm25_daily_mean
  )
]

setkey(
  pm_valid,
  pm25_site_id
)

cat(
  "\n有效PM site-day数 =",
  nrow(pm_valid),
  "\n"
)


# ============================================================
# 11. 从正式RadNet quasi-daily outcome表提取有效filter periods
#
# PM最终是给正式quasi-daily outcome服务的。
#
# 因此这里不再从更上游filter-period文件决定正式分析周期，
# 而是直接从正式quasi-daily表提取：
#
#   loc_num + collect_start + collect_end
#
# 只有真正形成至少1条quasi-daily observation的filter period
# 才进入PM匹配。
# ============================================================

cat(
  "\n================ 从正式quasi-daily提取filter periods ================\n"
)

quasi_period_source <- fread(
  QUASI_FILE
)

required_quasi_period <- c(
  "loc_num",
  "date",
  "collect_start",
  "collect_end"
)

missing_quasi_period <- setdiff(
  required_quasi_period,
  names(quasi_period_source)
)

if (length(missing_quasi_period) > 0) {

  cat(
    "\nquasi-daily字段：\n"
  )

  print(
    names(quasi_period_source)
  )

  stop(
    "正式quasi-daily文件缺少字段：",
    paste(
      missing_quasi_period,
      collapse = ", "
    )
  )
}

quasi_period_source[, loc_num := as.integer(loc_num)]
quasi_period_source[, date := as.IDate(date)]
quasi_period_source[, collect_start := as.IDate(collect_start)]
quasi_period_source[, collect_end := as.IDate(collect_end)]


# ------------------------------------------------------------
# 正式quasi-daily必须满足：
#
# collect_start <= date < collect_end
# ------------------------------------------------------------

bad_quasi_interval <- quasi_period_source[
  is.na(collect_start) |
  is.na(collect_end) |
  is.na(date) |
  collect_end <= collect_start |
  date < collect_start |
  date >= collect_end
]

cat(
  "正式quasi-daily总行数 =",
  nrow(quasi_period_source),
  "\n"
)

cat(
  "违反[collect_start, collect_end)定义的quasi-daily行数 =",
  nrow(bad_quasi_interval),
  "\n"
)

if (nrow(bad_quasi_interval) > 0) {

  print(
    bad_quasi_interval[
      1:min(30, .N)
    ]
  )

  stop(
    "正式quasi-daily表存在无效周期/日期关系，不能继续PM匹配。"
  )
}


# ------------------------------------------------------------
# 提取真正进入最终outcome universe的唯一filter periods
# ------------------------------------------------------------

periods <- unique(
  quasi_period_source[
    ,
    .(
      loc_num,
      collect_start,
      collect_end
    )
  ]
)

cat(
  "有效唯一filter periods =",
  nrow(periods),
  "\n"
)

cat(
  "有效filter-period loc_num数 =",
  uniqueN(periods$loc_num),
  "\n"
)


# ------------------------------------------------------------
# 上游filter-period文件仅做终端审计。
#
# 不修改原始文件，也不让无效周期进入PM匹配。
# ------------------------------------------------------------

if (file.exists(FILTER_FILE)) {

  filter_audit <- fread(
    FILTER_FILE,
    select = c(
      "loc_num",
      "collect_start",
      "collect_end"
    )
  )

  filter_audit[, loc_num := as.integer(loc_num)]
  filter_audit[, collect_start := as.IDate(collect_start)]
  filter_audit[, collect_end := as.IDate(collect_end)]

  bad_upstream <- filter_audit[
    is.na(collect_start) |
    is.na(collect_end) |
    collect_end <= collect_start
  ]

  cat(
    "上游filter-period表无效/非正长度记录 =",
    nrow(bad_upstream),
    "\n"
  )

  if (nrow(bad_upstream) > 0) {

    cat(
      "这些记录不会用于PM匹配，因为没有形成正式quasi-daily observation。\n"
    )

    cat(
      "其中 collect_end == collect_start =",
      bad_upstream[
        !is.na(collect_start) &
        !is.na(collect_end) &
        collect_end == collect_start,
        .N
      ],
      "\n"
    )

    cat(
      "其中 collect_end < collect_start =",
      bad_upstream[
        !is.na(collect_start) &
        !is.na(collect_end) &
        collect_end < collect_start,
        .N
      ],
      "\n"
    )

    cat(
      "其中日期缺失 =",
      bad_upstream[
        is.na(collect_start) |
        is.na(collect_end),
        .N
      ],
      "\n"
    )
  }
}


# ============================================================
# 12. Filter-period PM matching
#
# 核心：
#
# 对每个RadNet loc_num：
#   1) 找其50km PM stations
#   2) PM day作为长度为1天的时间点
#   3) 与RadNet [collect_start, collect_end) overlap
#   4) 同一filter period内所有有效PM site-day求mean
#
# 这里不先要求same-day quasi-daily匹配。
# ============================================================

cat(
  "\n================ 按filter period匹配PM2.5 ================\n"
)

loc_ids <- sort(
  unique(
    periods$loc_num
  )
)

period_pm_list <- vector(
  "list",
  length(loc_ids)
)

for (k in seq_along(loc_ids)) {

  id <- loc_ids[k]

  if (
    k == 1L ||
    k %% 10L == 0L ||
    k == length(loc_ids)
  ) {

    cat(
      sprintf(
        "处理RadNet %d/%d，loc_num=%s\n",
        k,
        length(loc_ids),
        id
      )
    )
  }

  p_i <- periods[
    loc_num == id
  ]

  site_ids <- links[
    loc_num == id,
    unique(pm25_site_id)
  ]

  # ----------------------------------------------------------
  # 12A. 50km内根本没有PM站
  # ----------------------------------------------------------

  if (length(site_ids) == 0) {

    out_i <- copy(p_i)

    out_i[
      ,
      `:=`(
        pm_c =
          NA_real_,

        n_pm_site_days_period =
          0L,

        n_pm_unique_dates_period =
          0L,

        n_pm_sites_reporting_period =
          0L
      )
    ]

    period_pm_list[[k]] <- out_i

    next
  }


  # ----------------------------------------------------------
  # 12B. 取该RadNet关联的所有有效PM site-day
  # ----------------------------------------------------------

  pm_i <- pm_valid[
    pm25_site_id %in% site_ids
  ]

  if (nrow(pm_i) == 0) {

    out_i <- copy(p_i)

    out_i[
      ,
      `:=`(
        pm_c =
          NA_real_,

        n_pm_site_days_period =
          0L,

        n_pm_unique_dates_period =
          0L,

        n_pm_sites_reporting_period =
          0L
      )
    ]

    period_pm_list[[k]] <- out_i

    next
  }


  # ----------------------------------------------------------
  # 12C. foverlaps
  #
  # RadNet period：
  #   period_start = collect_start
  #   period_end   = collect_end - 1
  #
  # PM：
  #   pm_start = pm_end = date_local
  # ----------------------------------------------------------

  p_int <- copy(p_i)

  p_int[
    ,
    `:=`(
      period_start =
        collect_start,

      period_end =
        collect_end - 1L
    )
  ]

  setkey(
    p_int,
    period_start,
    period_end
  )

  pm_int <- copy(pm_i)

  pm_int[
    ,
    `:=`(
      pm_start =
        date_local,

      pm_end =
        date_local
    )
  ]

  matched <- foverlaps(
    x = pm_int,
    y = p_int,

    by.x = c(
      "pm_start",
      "pm_end"
    ),

    by.y = c(
      "period_start",
      "period_end"
    ),

    type = "within",
    nomatch = 0L
  )


  # ----------------------------------------------------------
  # 12D. 按filter period求PM平均
  #
  # pm_c = filter period内全部匹配PM site-day的均值
  #
  # 这对应作者代码中：
  #
  # group_by(collect_start, collect_end)
  # mean(mass)
  #
  # 不按距离加权。
  # 不选最近站。
  # ----------------------------------------------------------

  if (nrow(matched) == 0) {

    agg <- data.table(
      loc_num = integer(),
      collect_start = as.IDate(character()),
      collect_end = as.IDate(character()),
      pm_c = numeric(),
      n_pm_site_days_period = integer(),
      n_pm_unique_dates_period = integer(),
      n_pm_sites_reporting_period = integer()
    )

  } else {

    agg <- matched[
      ,
      .(
        pm_c =
          mean(
            pm25_daily_mean
          ),

        n_pm_site_days_period =
          .N,

        n_pm_unique_dates_period =
          uniqueN(
            date_local
          ),

        n_pm_sites_reporting_period =
          uniqueN(
            pm25_site_id
          )
      ),
      by = .(
        loc_num,
        collect_start,
        collect_end
      )
    ]
  }


  # ----------------------------------------------------------
  # 12E. Keep every filter period, including missing PM
  # ----------------------------------------------------------

  out_i <- merge(
    p_i,
    agg,
    by = c(
      "loc_num",
      "collect_start",
      "collect_end"
    ),
    all.x = TRUE,
    sort = FALSE
  )

  out_i[
    is.na(n_pm_site_days_period),
    n_pm_site_days_period := 0L
  ]

  out_i[
    is.na(n_pm_unique_dates_period),
    n_pm_unique_dates_period := 0L
  ]

  out_i[
    is.na(n_pm_sites_reporting_period),
    n_pm_sites_reporting_period := 0L
  ]

  period_pm_list[[k]] <- out_i
}


# ============================================================
# 13. Combine filter-period PM
# ============================================================

period_pm <- rbindlist(
  period_pm_list,
  use.names = TRUE,
  fill = TRUE
)

if (nrow(period_pm) != nrow(periods)) {

  stop(
    "filter-period PM输出行数发生变化：",
    nrow(periods),
    " -> ",
    nrow(period_pm)
  )
}

if (
  period_pm[
    ,
    anyDuplicated(
      paste(
        loc_num,
        collect_start,
        collect_end
      )
    )
  ] > 0
) {

  stop(
    "filter-period PM输出主键不唯一。"
  )
}


# ============================================================
# 14. Add station metadata + spatial status
# ============================================================

station_small <- station_info[
  ,
  .(
    loc_num,
    radnet_monitor_id,
    station,
    city_name,
    state_abbr,
    main_analysis_conus,
    pm50_has_link,
    pm50_n_sites
  )
]

period_pm <- merge(
  period_pm,
  station_small,
  by = "loc_num",
  all.x = TRUE,
  sort = FALSE
)

period_pm[
  ,
  pm_status :=
    fifelse(
      !pm50_has_link,

      "01_no_PM_station_within_50km",

      fifelse(
        is.na(pm_c),

        "02_no_PM_measurement_within_filter_period",

        "00_PM_available"
      )
    )
]


# ============================================================
# 15. Formal filter-period output
# ============================================================

setcolorder(
  period_pm,
  c(
    "loc_num",
    "radnet_monitor_id",
    "station",
    "city_name",
    "state_abbr",

    "collect_start",
    "collect_end",

    "main_analysis_conus",

    "pm50_has_link",
    "pm50_n_sites",

    "pm_c",

    "n_pm_site_days_period",
    "n_pm_unique_dates_period",
    "n_pm_sites_reporting_period",

    "pm_status"
  )
)

setorder(
  period_pm,
  loc_num,
  collect_start,
  collect_end
)

fwrite(
  period_pm,
  OUT_PERIOD
)


# ============================================================
# 16. Read RadNet quasi-daily formal table
# ============================================================

cat(
  "\n================ 构建quasi-daily PM表 ================\n"
)

quasi <- fread(
  QUASI_FILE
)

required_quasi <- c(
  "loc_num",
  "date",
  "collect_start",
  "collect_end"
)

missing_quasi <- setdiff(
  required_quasi,
  names(quasi)
)

if (length(missing_quasi) > 0) {

  cat(
    "\nquasi-daily字段：\n"
  )

  print(
    names(quasi)
  )

  stop(
    "quasi-daily缺少字段：",
    paste(
      missing_quasi,
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


# ============================================================
# 17. Quasi-daily observation key
# ============================================================

quasi_key <- quasi[
  ,
  .(
    loc_num,
    date,
    collect_start,
    collect_end
  )
]

dup_quasi_key <- quasi_key[
  ,
  .N,
  by = .(
    loc_num,
    date,
    collect_start,
    collect_end
  )
][
  N > 1
]

if (nrow(dup_quasi_key) > 0) {

  stop(
    "quasi-daily observation key不唯一：",
    nrow(dup_quasi_key)
  )
}


# ============================================================
# 18. Propagate filter-period PM to every quasi-daily row
#
# PM不是按quasi-day重新匹配；
# 而是同一filter period所有quasi-day共享同一个pm_c。
# ============================================================

pm_period_for_join <- period_pm[
  ,
  .(
    loc_num,
    collect_start,
    collect_end,

    pm_c,

    pm50_has_link,
    pm50_n_sites,

    n_pm_site_days_period,
    n_pm_unique_dates_period,
    n_pm_sites_reporting_period,

    pm_status
  )
]

quasi_pm <- merge(
  quasi_key,
  pm_period_for_join,

  by = c(
    "loc_num",
    "collect_start",
    "collect_end"
  ),

  all.x = TRUE,
  sort = FALSE
)

if (nrow(quasi_pm) != nrow(quasi_key)) {

  stop(
    "PM merge改变了quasi-daily行数：",
    nrow(quasi_key),
    " -> ",
    nrow(quasi_pm)
  )
}


# ============================================================
# 19. Add station metadata
# ============================================================

quasi_pm <- merge(
  quasi_pm,
  station_small,
  by = c(
    "loc_num",
    "pm50_has_link",
    "pm50_n_sites"
  ),
  all.x = TRUE,
  sort = FALSE
)

if (nrow(quasi_pm) != nrow(quasi_key)) {

  stop(
    "station metadata merge改变了quasi-daily行数。"
  )
}

setcolorder(
  quasi_pm,
  c(
    "loc_num",
    "radnet_monitor_id",
    "station",
    "city_name",
    "state_abbr",

    "date",
    "collect_start",
    "collect_end",

    "main_analysis_conus",

    "pm50_has_link",
    "pm50_n_sites",

    "pm_c",

    "n_pm_site_days_period",
    "n_pm_unique_dates_period",
    "n_pm_sites_reporting_period",

    "pm_status"
  )
)

setorder(
  quasi_pm,
  loc_num,
  collect_start,
  collect_end,
  date
)

fwrite(
  quasi_pm,
  OUT_QUASI
)


# ============================================================
# 20. Update formal 181-station CONUS/PM coverage table
#
# 覆盖旧same-day统计，使该表同步到filter-period方法。
# ============================================================

quasi_station <- quasi_pm[
  ,
  .(
    quasi_rows_2001_2024 =
      .N,

    quasi_unique_dates_2001_2024 =
      uniqueN(date),

    quasi_first_date =
      min(date),

    quasi_last_date =
      max(date),

    quasi_rows_pm_available =
      sum(
        !is.na(pm_c)
      ),

    quasi_rows_pm_missing_total =
      sum(
        is.na(pm_c)
      ),

    quasi_rows_pm_missing_no50km =
      sum(
        pm_status ==
          "01_no_PM_station_within_50km"
      ),

    quasi_rows_pm_missing_filter_period =
      sum(
        pm_status ==
          "02_no_PM_measurement_within_filter_period"
      )
  ),
  by = loc_num
]

station_out <- merge(
  station_info,
  quasi_station,
  by = "loc_num",
  all.x = TRUE,
  sort = FALSE
)

count_cols <- c(
  "quasi_rows_2001_2024",
  "quasi_unique_dates_2001_2024",
  "quasi_rows_pm_available",
  "quasi_rows_pm_missing_total",
  "quasi_rows_pm_missing_no50km",
  "quasi_rows_pm_missing_filter_period"
)

for (f in count_cols) {

  set(
    station_out,
    which(
      is.na(
        station_out[[f]]
      )
    ),
    f,
    0L
  )
}

formal_station_cols <- c(
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
  "quasi_rows_pm_missing_filter_period"
)

setorder(
  station_out,
  loc_num
)

fwrite(
  station_out[
    ,
    ..formal_station_cols
  ],
  OUT_STATION
)


# ============================================================
# 21. CONUS audit
# ============================================================

conus_period <- period_pm[
  main_analysis_conus == TRUE
]

conus_quasi <- quasi_pm[
  main_analysis_conus == TRUE
]

cat(
  "\n\n============================================================\n"
)

cat(
  "PM2.5 —— 文章50 km + filter-period方法结果\n"
)

cat(
  "============================================================\n"
)

cat(
  "RadNet总数 =",
  nrow(cw),
  "\n"
)

cat(
  "CONUS loc_num =",
  cw[
    main_analysis_conus == TRUE,
    uniqueN(loc_num)
  ],
  "\n"
)

cat(
  "CONUS中有>=1个50 km PM站 =",
  station_info[
    main_analysis_conus == TRUE &
    pm50_has_link == TRUE,
    uniqueN(loc_num)
  ],
  "\n"
)

cat(
  "CONUS中无50 km PM站 =",
  station_info[
    main_analysis_conus == TRUE &
    pm50_has_link == FALSE,
    uniqueN(loc_num)
  ],
  "\n"
)


# ============================================================
# 22. Five structural no-PM stations
# ============================================================

cat(
  "\n================ CONUS无50 km PM站 ================\n"
)

print(
  station_out[
    main_analysis_conus == TRUE &
    pm50_has_link == FALSE,
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
# 23. Filter-period retention
# ============================================================

period_status <- conus_period[
  ,
  .(
    filter_periods = .N
  ),
  by = pm_status
]

cat(
  "\n================ CONUS filter-period PM覆盖 ================\n"
)

print(
  period_status[
    order(pm_status)
  ]
)

n_conus_period <-
  nrow(conus_period)

n_conus_period_pm <-
  conus_period[
    !is.na(pm_c),
    .N
  ]

cat(
  "\nCONUS filter periods总数 =",
  n_conus_period,
  "\n"
)

cat(
  "有filter-period PM的period数 =",
  n_conus_period_pm,
  "\n"
)

cat(
  "filter-period PM保留率 =",
  sprintf(
    "%.3f%%",
    100 *
    n_conus_period_pm /
    n_conus_period
  ),
  "\n"
)


# ============================================================
# 24. Quasi-daily retention
# ============================================================

quasi_status <- conus_quasi[
  ,
  .(
    quasi_rows = .N
  ),
  by = pm_status
]

cat(
  "\n================ CONUS quasi-daily PM覆盖 ================\n"
)

print(
  quasi_status[
    order(pm_status)
  ]
)

n_conus_quasi <-
  nrow(conus_quasi)

n_conus_quasi_pm <-
  conus_quasi[
    !is.na(pm_c),
    .N
  ]

n_conus_quasi_missing <-
  conus_quasi[
    is.na(pm_c),
    .N
  ]

cat(
  "\nCONUS quasi-daily总行数 =",
  n_conus_quasi,
  "\n"
)

cat(
  "有PM的quasi-daily行数 =",
  n_conus_quasi_pm,
  "\n"
)

cat(
  "PM缺失quasi-daily行数 =",
  n_conus_quasi_missing,
  "\n"
)

cat(
  "按作者filter-period方法的PM保留率 =",
  sprintf(
    "%.3f%%",
    100 *
    n_conus_quasi_pm /
    n_conus_quasi
  ),
  "\n"
)

cat(
  "按作者filter-period方法的PM损失率 =",
  sprintf(
    "%.3f%%",
    100 *
    n_conus_quasi_missing /
    n_conus_quasi
  ),
  "\n"
)


# ============================================================
# 25. Split structural vs temporal missing
# ============================================================

structural_missing <- conus_quasi[
  pm_status ==
    "01_no_PM_station_within_50km",
  .N
]

temporal_missing <- conus_quasi[
  pm_status ==
    "02_no_PM_measurement_within_filter_period",
  .N
]

cat(
  "\n其中：\n"
)

cat(
  "50 km内根本无PM站导致缺失 =",
  structural_missing,
  "\n"
)

cat(
  "50 km有PM站但整个filter period无PM测量 =",
  temporal_missing,
  "\n"
)


# ============================================================
# 26. Stations with any usable PM
# ============================================================

cat(
  "\nCONUS中至少有一个filter period获得PM的loc_num =",
  conus_period[
    !is.na(pm_c),
    uniqueN(loc_num)
  ],
  "\n"
)


# ============================================================
# 27. Spatially linked but never any usable filter-period PM
# ============================================================

station_temporal <- conus_period[
  pm50_has_link == TRUE,
  .(
    n_filter_periods =
      .N,

    n_filter_periods_with_pm =
      sum(
        !is.na(pm_c)
      )
  ),
  by = .(
    loc_num,
    radnet_monitor_id,
    station,
    state_abbr
  )
][
  n_filter_periods_with_pm == 0
]

cat(
  "\n================ 有50 km PM站但所有filter period均无PM ================\n"
)

cat(
  "站数 =",
  nrow(station_temporal),
  "\n"
)

if (nrow(station_temporal) > 0) {

  print(
    station_temporal[
      order(
        state_abbr,
        station
      )
    ]
  )
}


# ============================================================
# 28. PM distribution
# ============================================================

cat(
  "\n================ CONUS filter-period pm_c分布 ================\n"
)

print(
  summary(
    conus_period$pm_c
  )
)

cat(
  "\n================ CONUS quasi-daily pm_c分布 ================\n"
)

print(
  summary(
    conus_quasi$pm_c
  )
)


# ============================================================
# 29. Core integrity checks
# ============================================================

if (
  cw[
    main_analysis_conus == TRUE,
    uniqueN(loc_num)
  ] != 168
) {

  stop(
    "CONUS loc_num不是168。"
  )
}

if (
  station_info[
    main_analysis_conus == TRUE &
    pm50_has_link == TRUE,
    uniqueN(loc_num)
  ] != 163
) {

  stop(
    "CONUS有50km PM链接站不是163。"
  )
}

if (
  station_info[
    main_analysis_conus == TRUE &
    pm50_has_link == FALSE,
    uniqueN(loc_num)
  ] != 5
) {

  stop(
    "CONUS无50km PM链接站不是5。"
  )
}

if (nrow(quasi_pm) != nrow(quasi)) {

  stop(
    "quasi-daily正式输出行数与RadNet底表不同。"
  )
}

if (
  period_pm[
    !is.na(pm_c) &
    n_pm_site_days_period <= 0,
    .N
  ] > 0
) {

  stop(
    "存在pm_c非NA但n_pm_site_days_period<=0。"
  )
}

if (
  period_pm[
    is.na(pm_c) &
    n_pm_site_days_period > 0,
    .N
  ] > 0
) {

  stop(
    "存在有PM site-day但pm_c仍为NA。"
  )
}


# ============================================================
# 30. Old same-day file warning
# ============================================================

if (file.exists(OLD_DAILY_FILE)) {

  cat(
    "\n注意：旧same-day PM文件仍存在：\n",
    OLD_DAILY_FILE,
    "\n",
    sep = ""
  )

  cat(
    "它不再作为正式PM输入；",
    "正式模型请使用新的03_quasi_daily文件。\n"
  )
}


# ============================================================
# 31. PASS
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "PM2.5文章方法重做完成：PASS\n"
)

cat(
  "============================================================\n"
)

cat(
  "\n正式50 km链接表：\n",
  OUT_LINK,
  "\n",
  sep = ""
)

cat(
  "\n正式filter-period PM表：\n",
  OUT_PERIOD,
  "\n",
  sep = ""
)

cat(
  "\n正式quasi-daily PM表：\n",
  OUT_QUASI,
  "\n",
  sep = ""
)

cat(
  "\n更新后的181站CONUS/PM覆盖表：\n",
  OUT_STATION,
  "\n",
  sep = ""
)

cat(
  "\n运行结束。\n"
)
