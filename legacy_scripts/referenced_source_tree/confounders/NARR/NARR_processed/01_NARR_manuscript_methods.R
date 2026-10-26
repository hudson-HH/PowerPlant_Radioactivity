suppressPackageStartupMessages({
  library(data.table)
  library(raster)
  library(sp)
})

# ============================================================
# 0. 路径
# ============================================================

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

NARR_ROOT <- file.path(
  ROOT,
  "混杂因素",
  "NARR"
)

OUT_DIR <- file.path(
  NARR_ROOT,
  "NARR处理版"
)

CW_FILE <- file.path(
  ROOT,
  "radnet命名标准化",
  "radnet_181站_loc_num_radnet_monitor_id对照表.csv"
)

QUASI_FILE <- file.path(
  ROOT,
  "radnet 监测数据",
  "04_正式输出",
  "radnet_gross_beta_quasi_daily_model_ready_2001_2024.csv.gz"
)

OUT_DAILY <- file.path(
  OUT_DIR,
  "01_radnet_narr_daily_文章方法_2001_2024.csv.gz"
)

OUT_QUASI <- file.path(
  OUT_DIR,
  "02_radnet_narr_quasi_daily_文章方法_2001_2024.csv.gz"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

YEARS <- 2001:2024

VARS <- c(
  "uwnd.10m",
  "vwnd.10m",
  "rhum.2m",
  "air.2m",
  "hpbl",
  "soilm"
)


# ============================================================
# 1. 检查输入
# ============================================================

if (!file.exists(CW_FILE)) {
  stop("缺少181站crosswalk。")
}

if (!file.exists(QUASI_FILE)) {
  stop("缺少正式RadNet quasi-daily文件。")
}

# ============================================================
# 自动递归查找NARR NetCDF
# ============================================================

all_nc <- list.files(
  NARR_ROOT,
  pattern = "\\.nc$",
  recursive = TRUE,
  full.names = TRUE
)

cat(
  "NARR目录下找到NetCDF文件总数 =",
  length(all_nc),
  "\n"
)

find_narr_file <- function(var, year) {

  target <- paste0(
    var,
    ".",
    year,
    ".nc"
  )

  hit <- all_nc[
    basename(all_nc) == target
  ]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  if (length(hit) > 1) {

    cat(
      "\n发现同名NARR文件：",
      target,
      "\n",
      sep = ""
    )

    print(hit)

    stop(
      "同一个NARR变量年份存在多个同名文件。"
    )
  }

  hit
}

need <- CJ(
  var = VARS,
  year = YEARS
)

need[
  ,
  file :=
    mapply(
      find_narr_file,
      var,
      year,
      USE.NAMES = FALSE
    )
]

missing_files <- need[
  is.na(file)
]

if (nrow(missing_files) > 0) {

  cat(
    "\n真正缺失的NARR文件：\n"
  )

  print(missing_files)

  stop(
    "真正缺失NARR NetCDF文件，数量=",
    nrow(missing_files)
  )
}

cat(
  "\nNARR NetCDF检查：PASS\n"
)

cat(
  "文件数 =",
  nrow(need),
  "\n"
)


# ============================================================
# 2. 181站crosswalk
# ============================================================

cw <- fread(
  CW_FILE
)

required_cw <- c(
  "loc_num",
  "radnet_monitor_id",
  "station",
  "state_abbr",
  "monitor_latitude",
  "monitor_longitude"
)

miss_cw <- setdiff(
  required_cw,
  names(cw)
)

if (length(miss_cw) > 0) {

  stop(
    "crosswalk缺少字段：",
    paste(
      miss_cw,
      collapse = ", "
    )
  )
}

cw[
  ,
  loc_num :=
    as.integer(loc_num)
]

if (nrow(cw) != 181) {
  stop("crosswalk不是181行。")
}

if (uniqueN(cw$loc_num) != 181) {
  stop("loc_num不是181个唯一站。")
}

setorder(
  cw,
  loc_num
)


# ============================================================
# 3. CONUS标记
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


# ============================================================
# 4. 建立RadNet空间点
#
# 与作者代码一致：
# 原始点为WGS84，
# 每次提取前再spTransform到NARR raster CRS。
# ============================================================

radnet_points <- SpatialPointsDataFrame(
  coords = cbind(
    cw$monitor_longitude,
    cw$monitor_latitude
  ),

  data = data.frame(
    point_index =
      seq_len(
        nrow(cw)
      ),

    loc_num =
      cw$loc_num
  ),

  proj4string =
    CRS(
      "+proj=longlat +datum=WGS84 +no_defs"
    )
)

N_POINTS <- nrow(cw)


# ============================================================
# 5. 作者式单变量、单年份NARR提取函数
#
# 对应作者：
#
# variable <- stack(...)
# points <- spTransform(points, proj4string(variable))
# ext_value <- extract(variable, points)
# ext_value <- c(ext_value)
#
# 然后按layer顺序生成DOY。
# ============================================================

extract_narr_one <- function(
  var,
  year
) {

  f <- find_narr_file(
    var,
    year
  )

  if (is.na(f)) {
    stop(
      "找不到NARR文件：",
      var,
      ".",
      year,
      ".nc"
    )
  }

  s <- raster::stack(
    f
  )

  pts <- sp::spTransform(
    radnet_points,
    raster::crs(s)
  )

  m <- raster::extract(
    s,
    pts
  )

  if (is.null(dim(m))) {

    m <- matrix(
      m,
      nrow = N_POINTS,
      ncol = raster::nlayers(s)
    )
  }

  n_days <- raster::nlayers(
    s
  )

  if (nrow(m) != N_POINTS) {

    stop(
      var,
      " ",
      year,
      "：提取站点数不等于181。"
    )
  }

  if (ncol(m) != n_days) {

    stop(
      var,
      " ",
      year,
      "：提取layer数量错误。"
    )
  }

  # raster::extract返回：
  # 行 = RadNet站
  # 列 = 日期layer
  #
  # 作者使用c(ext_value)，R按列展开：
  # day1全部181站
  # day2全部181站
  # ...

  value <- c(
    m
  )

  out <- data.table(
    loc_num =
      rep(
        cw$loc_num,
        times = n_days
      ),

    doy =
      rep(
        seq_len(n_days),
        each = N_POINTS
      ),

    value =
      value
  )

  setnames(
    out,
    "value",
    var
  )

  out[
    ,
    year :=
      year
  ]

  out[
    ,
    date :=
      as.IDate(
        sprintf(
          "%04d-01-01",
          year
        )
      ) +
      doy -
      1L
  ]

  out
}


# ============================================================
# 6. 逐年重建六个NARR变量
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "开始按作者方法提取NARR\n"
)

cat(
  "============================================================\n"
)

year_list <- vector(
  "list",
  length(YEARS)
)

for (
  yy_i in seq_along(YEARS)
) {

  yy <- YEARS[yy_i]

  cat(
    "\n---------------- ",
    yy,
    " ----------------\n",
    sep = ""
  )

  annual <- NULL

  for (
    vv in VARS
  ) {

    cat(
      "提取 ",
      vv,
      " ... ",
      sep = ""
    )

    t <- extract_narr_one(
      vv,
      yy
    )

    if (is.null(annual)) {

      annual <- t

    } else {

      if (
        !identical(
          annual$loc_num,
          t$loc_num
        ) ||
        !identical(
          annual$doy,
          t$doy
        )
      ) {

        stop(
          yy,
          "年不同NARR变量行顺序不一致。"
        )
      }

      annual[
        ,
        (vv) :=
          t[[vv]]
      ]
    }

    cat(
      "完成\n"
    )
  }


  # ==========================================================
  # 作者原始计算
  # ==========================================================

  annual[
    ,
    wdir :=
      (
        270 -
        atan2(
          vwnd.10m,
          uwnd.10m
        ) *
        180 /
        pi
      ) %%
      360
  ]

  annual[
    ,
    vel :=
      sqrt(
        vwnd.10m^2 +
        uwnd.10m^2
      )
  ]


  # ==========================================================
  # 作者正式模型中的别名
  #
  # 注意：
  # temp不转摄氏度。
  # ==========================================================

  annual[
    ,
    temp :=
      air.2m
  ]

  annual[
    ,
    inv_hpbl :=
      1 /
      hpbl
  ]

  annual[
    ,
    rhum :=
      rhum.2m
  ]

  year_list[[yy_i]] <- annual

  cat(
    yy,
    "完成，行数=",
    nrow(annual),
    "\n"
  )
}


# ============================================================
# 7. 合并2001-2024
# ============================================================

narr <- rbindlist(
  year_list,
  use.names = TRUE,
  fill = TRUE
)

setkey(
  narr,
  loc_num,
  date
)

dup <- narr[
  ,
  .N,
  by = .(
    loc_num,
    date
  )
][
  N > 1
]

if (nrow(dup) > 0) {

  stop(
    "NARR存在重复loc_num+date：",
    nrow(dup)
  )
}


# ============================================================
# 8. 加站点信息
# ============================================================

station_info <- cw[
  ,
  .(
    loc_num,
    radnet_monitor_id,
    station,
    state_abbr,
    monitor_latitude,
    monitor_longitude,
    main_analysis_conus
  )
]

narr <- merge(
  narr,
  station_info,
  by = "loc_num",
  all.x = TRUE,
  sort = FALSE
)

setcolorder(
  narr,
  c(
    "loc_num",
    "radnet_monitor_id",
    "station",
    "state_abbr",
    "monitor_latitude",
    "monitor_longitude",
    "main_analysis_conus",

    "date",
    "year",
    "doy",

    "air.2m",
    "temp",

    "hpbl",
    "inv_hpbl",

    "rhum.2m",
    "rhum",

    "soilm",

    "uwnd.10m",
    "vwnd.10m",

    "vel",
    "wdir"
  )
)

setorder(
  narr,
  loc_num,
  date
)


# ============================================================
# 9. 写正式daily NARR
# ============================================================

fwrite(
  narr,
  OUT_DAILY
)


# ============================================================
# 10. 与正式RadNet quasi-daily匹配
#
# NARR按：
# loc_num + date
#
# 逐日匹配。
#
# 不做filter-period平均。
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "匹配正式RadNet quasi-daily\n"
)

cat(
  "============================================================\n"
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

miss_quasi <- setdiff(
  required_quasi,
  names(quasi)
)

if (length(miss_quasi) > 0) {

  stop(
    "quasi-daily缺少字段：",
    paste(
      miss_quasi,
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
# 11. 保留正式quasi-daily观测键
#
# overlapping filter periods不去重。
# ============================================================

qkey <- quasi[
  ,
  .(
    loc_num,
    date,
    collect_start,
    collect_end
  )
]

if (
  qkey[
    ,
    .N,
    by = .(
      loc_num,
      date,
      collect_start,
      collect_end
    )
  ][
    N > 1,
    .N
  ] > 0
) {

  stop(
    "正式quasi-daily observation key不唯一。"
  )
}


# ============================================================
# 12. NARR daily join
# ============================================================

narr_join <- narr[
  ,
  .(
    loc_num,
    date,

    air.2m,
    temp,

    hpbl,
    inv_hpbl,

    rhum.2m,
    rhum,

    soilm,

    uwnd.10m,
    vwnd.10m,

    vel,
    wdir
  )
]

q_narr <- merge(
  qkey,
  narr_join,

  by = c(
    "loc_num",
    "date"
  ),

  all.x = TRUE,
  sort = FALSE
)

if (nrow(q_narr) != nrow(qkey)) {

  stop(
    "NARR join改变了quasi-daily行数：",
    nrow(qkey),
    " -> ",
    nrow(q_narr)
  )
}


# ============================================================
# 13. 加站点信息与CONUS
# ============================================================

q_narr <- merge(
  q_narr,
  station_info,

  by = "loc_num",

  all.x = TRUE,
  sort = FALSE
)

if (nrow(q_narr) != nrow(qkey)) {

  stop(
    "station metadata join改变了行数。"
  )
}

setcolorder(
  q_narr,
  c(
    "loc_num",
    "radnet_monitor_id",
    "station",
    "state_abbr",
    "main_analysis_conus",

    "date",
    "collect_start",
    "collect_end",

    "air.2m",
    "temp",

    "hpbl",
    "inv_hpbl",

    "rhum.2m",
    "rhum",

    "soilm",

    "uwnd.10m",
    "vwnd.10m",

    "vel",
    "wdir"
  )
)

setorder(
  q_narr,
  loc_num,
  collect_start,
  collect_end,
  date
)

fwrite(
  q_narr,
  OUT_QUASI
)


# ============================================================
# 14. CONUS覆盖检查
# ============================================================

core_vars <- c(
  "air.2m",
  "hpbl",
  "rhum.2m",
  "soilm",
  "uwnd.10m",
  "vwnd.10m",
  "vel"
)

q_narr[
  ,
  narr_complete :=
    complete.cases(
      .SD
    ),
  .SDcols = core_vars
]

conus <- q_narr[
  main_analysis_conus == TRUE
]

cat(
  "\n============================================================\n"
)

cat(
  "NARR文章方法结果\n"
)

cat(
  "============================================================\n"
)

cat(
  "181站总数 =",
  uniqueN(q_narr$loc_num),
  "\n"
)

cat(
  "CONUS站数 =",
  uniqueN(conus$loc_num),
  "\n"
)

cat(
  "CONUS quasi-daily总行数 =",
  nrow(conus),
  "\n"
)

cat(
  "CONUS NARR完整行数 =",
  conus[
    narr_complete == TRUE,
    .N
  ],
  "\n"
)

cat(
  "CONUS NARR缺失行数 =",
  conus[
    narr_complete == FALSE,
    .N
  ],
  "\n"
)

cat(
  "CONUS NARR完整率 =",
  sprintf(
    "%.3f%%",
    100 *
    mean(
      conus$narr_complete
    )
  ),
  "\n"
)


# ============================================================
# 15. 若有缺失，打印站点
# ============================================================

missing_station <- conus[
  narr_complete == FALSE,
  .(
    missing_rows = .N,
    first_missing = min(date),
    last_missing = max(date)
  ),
  by = .(
    loc_num,
    radnet_monitor_id,
    station,
    state_abbr
  )
]

cat(
  "\n================ CONUS存在NARR缺失的站 ================\n"
)

cat(
  "站数 =",
  nrow(missing_station),
  "\n"
)

if (nrow(missing_station) > 0) {

  print(
    missing_station[
      order(
        -missing_rows
      )
    ]
  )
}


# ============================================================
# 16. 各变量分布
# ============================================================

cat(
  "\n================ CONUS NARR变量分布 ================\n"
)

for (
  v in c(
    "air.2m",
    "hpbl",
    "rhum.2m",
    "soilm",
    "uwnd.10m",
    "vwnd.10m",
    "vel"
  )
) {

  cat(
    "\n",
    v,
    "\n",
    sep = ""
  )

  print(
    summary(
      conus[[v]]
    )
  )
}


# ============================================================
# 17. 作者后续模型筛选条件，仅统计，不删除
# ============================================================

cat(
  "\n================ 作者后续模型筛选条件，仅审计 ================\n"
)

cat(
  "vel >= 6 的CONUS quasi-daily行数 =",
  conus[
    !is.na(vel) &
    vel >= 6,
    .N
  ],
  "\n"
)

cat(
  "hpbl <= 150 的CONUS quasi-daily行数 =",
  conus[
    !is.na(hpbl) &
    hpbl <= 150,
    .N
  ],
  "\n"
)

cat(
  "hpbl >= 3000 的CONUS quasi-daily行数 =",
  conus[
    !is.na(hpbl) &
    hpbl >= 3000,
    .N
  ],
  "\n"
)

cat(
  "注意：以上行没有在NARR重建阶段删除。\n"
)


# ============================================================
# 18. 最终状态
# ============================================================

cat(
  "\n============================================================\n"
)

if (
  conus[
    narr_complete == FALSE,
    .N
  ] == 0
) {

  cat(
    "NARR文章方法重建：PASS\n"
  )

} else {

  cat(
    "NARR文章方法重建完成，但CONUS存在缺失，需要检查。\n"
  )
}

cat(
  "============================================================\n"
)

cat(
  "\n正式daily NARR：\n",
  OUT_DAILY,
  "\n",
  sep = ""
)

cat(
  "\n正式quasi-daily NARR：\n",
  OUT_QUASI,
  "\n",
  sep = ""
)

cat(
  "\n运行结束。\n"
)
