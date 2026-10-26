suppressPackageStartupMessages({
  library(data.table)
})

# ============================================================
# 0. 项目路径
# ============================================================

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

OUT_DIR <- file.path(
  ROOT,
  "混杂因素",
  "混杂因素统一合并"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

OUT_FILE <- file.path(
  OUT_DIR,
  "01_radnet_181站_环境混杂因素统一表_2001_2024.csv.gz"
)


# ============================================================
# 1. 正式输入文件
# ============================================================

RADNET_FILE <- file.path(
  ROOT,
  "radnet 监测数据",
  "04_正式输出",
  "radnet_gross_beta_quasi_daily_model_ready_2001_2024.csv.gz"
)

CROSSWALK_FILE <- file.path(
  ROOT,
  "radnet命名标准化",
  "radnet_181站_loc_num_radnet_monitor_id对照表.csv"
)

PM_FILE <- file.path(
  ROOT,
  "混杂因素",
  "PM2.5",
  "PM2.5处理版",
  "03_radnet_pm25_quasi_daily_50km_文章方法_2001_2024.csv.gz"
)

NARR_FILE <- file.path(
  ROOT,
  "混杂因素",
  "NARR",
  "NARR处理版",
  "02_radnet_narr_quasi_daily_文章方法_2001_2024.csv.gz"
)

SUN_FILE <- file.path(
  ROOT,
  "混杂因素",
  "太阳黑子",
  "太阳黑子处理版",
  "02_radnet_sunspot_quasi_daily_文章方法_2001_2024.csv.gz"
)

HYSPLIT_FILE <- file.path(
  ROOT,
  "混杂因素",
  "08A4_HYSPLIT_air_mass_origin",
  "HYSPLIT处理版",
  "radnet_181站_HYSPLIT_daily_文章方法_2001_2024.csv.gz"
)


# ============================================================
# 2. 自动寻找“地表 U ”正式文件
#
# 因为该文件夹名称末尾可能包含空格。
# ============================================================

conf_root <- file.path(
  ROOT,
  "混杂因素"
)

u_dirs <- list.dirs(
  conf_root,
  recursive = FALSE,
  full.names = TRUE
)

u_dirs <- u_dirs[
  trimws(
    basename(u_dirs)
  ) == "地表 U"
]

if (length(u_dirs) != 1) {

  cat(
    "\n找到的地表U目录：\n"
  )

  print(
    u_dirs
  )

  stop(
    "无法唯一确定地表 U 文件夹。"
  )
}

U_DIR <- u_dirs[1]

U_FILE <- file.path(
  U_DIR,
  "地表 U处理版",
  "radnet_181站_地表U_文章方法.csv"
)

# 如果实际目录名是“地表U处理版”，自动兼容
if (!file.exists(U_FILE)) {

  U_FILE2 <- file.path(
    U_DIR,
    "地表U处理版",
    "radnet_181站_地表U_文章方法.csv"
  )

  if (file.exists(U_FILE2)) {
    U_FILE <- U_FILE2
  }
}


# ============================================================
# 3. 输入文件完整性检查
# ============================================================

FILES <- c(
  RadNet = RADNET_FILE,
  crosswalk = CROSSWALK_FILE,
  PM25 = PM_FILE,
  NARR = NARR_FILE,
  sunspot = SUN_FILE,
  HYSPLIT = HYSPLIT_FILE,
  surface_U = U_FILE
)

missing_files <- FILES[
  !file.exists(
    FILES
  )
]

if (length(missing_files) > 0) {

  cat(
    "\n缺失正式输入文件：\n"
  )

  print(
    missing_files
  )

  stop(
    "存在缺失输入文件，停止合并。"
  )
}

cat(
  "\n============================================================\n"
)

cat(
  "全部正式输入文件：PASS\n"
)

cat(
  "============================================================\n"
)


# ============================================================
# 4. 工具函数
# ============================================================

as_date_cols <- function(
  dt
) {

  for (
    v in intersect(
      c(
        "date",
        "Date",
        "collect_start",
        "collect_end"
      ),
      names(dt)
    )
  ) {

    set(
      dt,
      j = v,
      value = as.IDate(
        dt[[v]]
      )
    )
  }

  dt
}


assert_unique <- function(
  dt,
  keys,
  label
) {

  bad <- dt[
    ,
    .N,
    by = keys
  ][
    N > 1
  ]

  if (nrow(bad) > 0) {

    cat(
      "\n",
      label,
      " 存在重复键：\n",
      sep = ""
    )

    print(
      head(
        bad,
        20
      )
    )

    stop(
      label,
      " 连接键不唯一。"
    )
  }
}


# ============================================================
# 5. 读取正式RadNet outcome
# ============================================================

cat(
  "\n读取RadNet正式quasi-daily...\n"
)

rad <- fread(
  RADNET_FILE
)

rad <- as_date_cols(
  rad
)

if ("loc_num" %in% names(rad)) {
  rad[
    ,
    loc_num :=
      as.integer(loc_num)
  ]
}

RAD_KEY <- c(
  "loc_num",
  "collect_start",
  "collect_end",
  "date"
)

missing_key <- setdiff(
  RAD_KEY,
  names(rad)
)

if (length(missing_key) > 0) {

  stop(
    "RadNet正式表缺少观察键：",
    paste(
      missing_key,
      collapse = ", "
    )
  )
}

assert_unique(
  rad,
  RAD_KEY,
  "RadNet正式quasi-daily"
)

N0 <- nrow(rad)

cat(
  "RadNet行数 =",
  format(
    N0,
    big.mark = ","
  ),
  "\n"
)

cat(
  "RadNet站数 =",
  uniqueN(
    rad$loc_num
  ),
  "\n"
)

if (N0 != 769770) {

  warning(
    "RadNet正式quasi-daily行数不是预期769,770，请确认。"
  )
}


# ============================================================
# 6. 加正式181站metadata和CONUS标记
# ============================================================

cw <- fread(
  CROSSWALK_FILE
)

cw[
  ,
  loc_num :=
    as.integer(loc_num)
]

if (nrow(cw) != 181) {
  stop("crosswalk行数不是181。")
}

if (uniqueN(cw$loc_num) != 181) {
  stop("crosswalk中的loc_num不是181个唯一值。")
}

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
    state_abbr %in%
    CONUS_STATES
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

cw_keep <- intersect(
  c(
    "loc_num",
    "radnet_monitor_id",
    "station",
    "city_name",
    "state_abbr",
    "monitor_latitude",
    "monitor_longitude",
    "main_analysis_conus"
  ),
  names(cw)
)

# 避免RadNet本身已有同名metadata导致.x/.y
new_cw_cols <- setdiff(
  cw_keep,
  names(rad)
)

if (length(new_cw_cols) > 0) {

  rad <- merge(
    rad,
    cw[
      ,
      c(
        "loc_num",
        new_cw_cols
      ),
      with = FALSE
    ],
    by = "loc_num",
    all.x = TRUE,
    sort = FALSE
  )
}

if (nrow(rad) != N0) {
  stop("crosswalk连接改变了RadNet行数。")
}


# ============================================================
# 7. 地表U
#
# 连接键：loc_num
# ============================================================

cat(
  "\n---------------- 地表U ----------------\n"
)

u <- fread(
  U_FILE
)

u[
  ,
  loc_num :=
    as.integer(loc_num)
]

assert_unique(
  u,
  "loc_num",
  "地表U"
)

# 正式下游字段优先
if (
  !"surface_u_ppm_article" %in%
  names(u)
) {

  # 自动寻找可能的U字段
  candidates <- c(
    "U",
    "surface_u_ppm",
    "u_ppm",
    "uranium_ppm"
  )

  hit <- intersect(
    candidates,
    names(u)
  )

  if (length(hit) != 1) {

    cat(
      "地表U字段：\n"
    )

    print(
      names(u)
    )

    stop(
      "无法唯一识别地表U变量。"
    )
  }

  setnames(
    u,
    hit,
    "surface_u_ppm_article"
  )
}

u_join <- u[
  ,
  .(
    loc_num,
    surface_u_ppm_article
  )
]

rad <- merge(
  rad,
  u_join,
  by = "loc_num",
  all.x = TRUE,
  sort = FALSE
)

if (nrow(rad) != N0) {
  stop("地表U连接改变了主表行数。")
}

cat(
  "地表U缺失行数 =",
  rad[
    is.na(
      surface_u_ppm_article
    ),
    .N
  ],
  "\n"
)


# ============================================================
# 8. HYSPLIT
#
# 连接键：loc_num + date
# ============================================================

cat(
  "\n---------------- HYSPLIT ----------------\n"
)

hys <- fread(
  HYSPLIT_FILE
)

hys <- as_date_cols(
  hys
)

hys[
  ,
  loc_num :=
    as.integer(loc_num)
]

if (
  !"date" %in% names(hys) &&
  "Date" %in% names(hys)
) {

  setnames(
    hys,
    "Date",
    "date"
  )
}

if (
  !"c_prop" %in% names(hys)
) {

  candidates <- c(
    "continental_trajectory_fraction",
    "continental_fraction",
    "continental_origin_fraction"
  )

  hit <- intersect(
    candidates,
    names(hys)
  )

  if (length(hit) != 1) {

    cat(
      "HYSPLIT字段：\n"
    )

    print(
      names(hys)
    )

    stop(
      "无法唯一识别HYSPLIT c_prop。"
    )
  }

  setnames(
    hys,
    hit,
    "c_prop"
  )
}

assert_unique(
  hys,
  c(
    "loc_num",
    "date"
  ),
  "HYSPLIT daily"
)

hys_join <- hys[
  ,
  .(
    loc_num,
    date,
    c_prop
  )
]

rad <- merge(
  rad,
  hys_join,
  by = c(
    "loc_num",
    "date"
  ),
  all.x = TRUE,
  sort = FALSE
)

if (nrow(rad) != N0) {
  stop("HYSPLIT连接改变了主表行数。")
}

cat(
  "HYSPLIT c_prop缺失行数 =",
  rad[
    is.na(c_prop),
    .N
  ],
  "\n"
)


# ============================================================
# 9. NARR
#
# NARR物理上按 loc_num + date，
# 使用formal quasi表前先压缩成daily唯一值。
#
# 这样可以明确保证不同filter period重叠日
# 不会造成many-to-many。
# ============================================================

cat(
  "\n---------------- NARR ----------------\n"
)

narr <- fread(
  NARR_FILE
)

narr <- as_date_cols(
  narr
)

narr[
  ,
  loc_num :=
    as.integer(loc_num)
]

if (
  !"date" %in% names(narr) &&
  "Date" %in% names(narr)
) {

  setnames(
    narr,
    "Date",
    "date"
  )
}

NARR_VARS <- c(
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

missing_narr_vars <- setdiff(
  NARR_VARS,
  names(narr)
)

if (length(missing_narr_vars) > 0) {

  stop(
    "NARR正式表缺少字段：",
    paste(
      missing_narr_vars,
      collapse = ", "
    )
  )
}

# 检查同一loc_num+date下，
# 不同filter-period复制出来的NARR值是否完全一致
narr_check <- narr[
  ,
  lapply(
    .SD,
    uniqueN,
    na.rm = FALSE
  ),
  by = .(
    loc_num,
    date
  ),
  .SDcols = NARR_VARS
]

bad_narr <- narr_check[
  ,
  Reduce(
    `|`,
    lapply(
      .SD,
      function(x) x > 1
    )
  ),
  .SDcols = NARR_VARS
]

if (any(bad_narr)) {

  stop(
    "同一个loc_num+date在不同filter period中存在不一致NARR值。"
  )
}

narr_daily <- unique(
  narr[
    ,
    c(
      "loc_num",
      "date",
      NARR_VARS
    ),
    with = FALSE
  ],
  by = c(
    "loc_num",
    "date"
  )
)

assert_unique(
  narr_daily,
  c(
    "loc_num",
    "date"
  ),
  "NARR daily"
)

rad <- merge(
  rad,
  narr_daily,
  by = c(
    "loc_num",
    "date"
  ),
  all.x = TRUE,
  sort = FALSE
)

if (nrow(rad) != N0) {
  stop("NARR连接改变了主表行数。")
}

cat(
  "NARR soilm缺失行数 =",
  rad[
    is.na(soilm),
    .N
  ],
  "\n"
)


# ============================================================
# 10. PM2.5
#
# 关键：
# 必须按完整filter-period observation key连接。
# ============================================================

cat(
  "\n---------------- PM2.5 ----------------\n"
)

pm <- fread(
  PM_FILE
)

pm <- as_date_cols(
  pm
)

pm[
  ,
  loc_num :=
    as.integer(loc_num)
]

if (
  !"pm_c" %in% names(pm)
) {

  candidates <- c(
    "pm25",
    "pm",
    "pm25_filter_period_mean"
  )

  hit <- intersect(
    candidates,
    names(pm)
  )

  if (length(hit) != 1) {

    cat(
      "PM2.5字段：\n"
    )

    print(
      names(pm)
    )

    stop(
      "无法唯一识别PM2.5变量pm_c。"
    )
  }

  setnames(
    pm,
    hit,
    "pm_c"
  )
}

assert_unique(
  pm,
  RAD_KEY,
  "PM2.5 quasi-daily"
)

pm_join <- pm[
  ,
  .(
    loc_num,
    collect_start,
    collect_end,
    date,
    pm_c
  )
]

rad <- merge(
  rad,
  pm_join,
  by = RAD_KEY,
  all.x = TRUE,
  sort = FALSE
)

if (nrow(rad) != N0) {
  stop("PM2.5连接改变了主表行数。")
}

cat(
  "PM2.5 pm_c缺失行数 =",
  rad[
    is.na(pm_c),
    .N
  ],
  "\n"
)


# ============================================================
# 11. 太阳黑子
#
# 使用已经冻结的formal quasi-daily表。
#
# 仍然按完整观察键连接，
# 确保与769,770条RadNet observation完全一一对应。
# ============================================================

cat(
  "\n---------------- 太阳黑子 ----------------\n"
)

sun <- fread(
  SUN_FILE
)

sun <- as_date_cols(
  sun
)

sun[
  ,
  loc_num :=
    as.integer(loc_num)
]

if (
  !"Num" %in% names(sun)
) {

  if (
    "sunspot_monthly" %in%
    names(sun)
  ) {

    sun[
      ,
      Num :=
        sunspot_monthly
    ]

  } else {

    stop(
      "太阳黑子正式表没有Num或sunspot_monthly。"
    )
  }
}

assert_unique(
  sun,
  RAD_KEY,
  "太阳黑子 quasi-daily"
)

sun_join <- sun[
  ,
  .(
    loc_num,
    collect_start,
    collect_end,
    date,
    sunspot_monthly = Num
  )
]

rad <- merge(
  rad,
  sun_join,
  by = RAD_KEY,
  all.x = TRUE,
  sort = FALSE
)

if (nrow(rad) != N0) {
  stop("太阳黑子连接改变了主表行数。")
}

cat(
  "太阳黑子缺失行数 =",
  rad[
    is.na(
      sunspot_monthly
    ),
    .N
  ],
  "\n"
)


# ============================================================
# 12. 最终变量统一命名
# ============================================================

# 纬度统一提供一个模型字段latitude
if (
  !"latitude" %in% names(rad) &&
  "monitor_latitude" %in%
  names(rad)
) {

  rad[
    ,
    latitude :=
      monitor_latitude
  ]
}

if (
  !"longitude" %in% names(rad) &&
  "monitor_longitude" %in%
  names(rad)
) {

  rad[
    ,
    longitude :=
      monitor_longitude
  ]
}

# 模型常用别名
rad[
  ,
  surface_u_ppm :=
    surface_u_ppm_article
]

rad[
  ,
  continental_trajectory_fraction :=
    c_prop
]

rad[
  ,
  sun :=
    sunspot_monthly
]


# ============================================================
# 13. 再次确认观察键
# ============================================================

assert_unique(
  rad,
  RAD_KEY,
  "最终统一表"
)

if (nrow(rad) != N0) {
  stop("最终统一表行数发生变化。")
}

if (
  uniqueN(
    rad$loc_num
  ) != 181
) {

  stop(
    "最终统一表不是181个RadNet站。"
  )
}


# ============================================================
# 14. CONUS检查
# ============================================================

if (
  !"main_analysis_conus" %in%
  names(rad)
) {

  stop(
    "最终统一表缺少main_analysis_conus。"
  )
}

conus <- rad[
  main_analysis_conus == TRUE
]

if (
  uniqueN(
    conus$loc_num
  ) != 168
) {

  stop(
    "最终统一表CONUS站数不是168。"
  )
}

if (
  nrow(conus) != 741861
) {

  warning(
    "CONUS quasi-daily行数不是预期741,861。"
  )
}


# ============================================================
# 15. 各混杂因素覆盖统计
#
# 这里只统计，不删除任何行。
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "CONUS环境混杂因素覆盖情况\n"
)

cat(
  "============================================================\n"
)

cover_vars <- c(
  "surface_u_ppm_article",
  "c_prop",
  "temp",
  "hpbl",
  "inv_hpbl",
  "rhum",
  "soilm",
  "vel",
  "pm_c",
  "sunspot_monthly"
)

for (
  v in cover_vars
) {

  n_ok <- conus[
    !is.na(
      get(v)
    ),
    .N
  ]

  n_miss <- nrow(conus) - n_ok

  cat(
    sprintf(
      "%-28s 有值=%7d  缺失=%7d  完整率=%7.3f%%\n",
      v,
      n_ok,
      n_miss,
      100 *
      n_ok /
      nrow(conus)
    )
  )
}


# ============================================================
# 16. 完整环境混杂因素样本，仅统计
#
# 注意：
# 不在这里删除。
# ============================================================

core_complete_vars <- c(
  "surface_u_ppm_article",
  "c_prop",
  "temp",
  "inv_hpbl",
  "rhum",
  "soilm",
  "vel",
  "pm_c",
  "sunspot_monthly"
)

conus[
  ,
  environmental_complete :=
    complete.cases(
      .SD
    ),
  .SDcols = core_complete_vars
]

n_complete <- conus[
  environmental_complete == TRUE,
  .N
]

cat(
  "\nCONUS全部核心环境混杂因素同时完整：",
  format(
    n_complete,
    big.mark = ","
  ),
  " / ",
  format(
    nrow(conus),
    big.mark = ","
  ),
  " (",
  sprintf(
    "%.3f%%",
    100 *
    n_complete /
    nrow(conus)
  ),
  ")\n",
  sep = ""
)


# ============================================================
# 17. 作者后续模型筛选条件
#
# 仅做审计，不在统一合并阶段删除。
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "参考文章后续模型筛选条件，仅统计，不执行\n"
)

cat(
  "============================================================\n"
)

cat(
  "pm_c >= 35：",
  conus[
    !is.na(pm_c) &
    pm_c >= 35,
    .N
  ],
  "\n"
)

cat(
  "vel >= 6：",
  conus[
    !is.na(vel) &
    vel >= 6,
    .N
  ],
  "\n"
)

cat(
  "hpbl <= 150：",
  conus[
    !is.na(hpbl) &
    hpbl <= 150,
    .N
  ],
  "\n"
)

cat(
  "hpbl >= 3000：",
  conus[
    !is.na(hpbl) &
    hpbl >= 3000,
    .N
  ],
  "\n"
)

if (
  "gross_beta" %in%
  names(conus)
) {

  cat(
    "gross_beta >= 0.04：",
    conus[
      !is.na(gross_beta) &
      gross_beta >= 0.04,
      .N
    ],
    "\n"
  )
}


# ============================================================
# 18. 最终排序
# ============================================================

setorder(
  rad,
  loc_num,
  collect_start,
  collect_end,
  date
)


# ============================================================
# 19. 输出正式统一表
# ============================================================

fwrite(
  rad,
  OUT_FILE
)


# ============================================================
# 20. 回读检查
# ============================================================

check <- fread(
  OUT_FILE,
  select = c(
    "loc_num",
    "collect_start",
    "collect_end",
    "date"
  )
)

if (nrow(check) != N0) {
  stop("输出文件回读后行数不一致。")
}

if (uniqueN(check$loc_num) != 181) {
  stop("输出文件回读后站数不一致。")
}


# ============================================================
# 21. 最终结果
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "环境混杂因素统一合并：PASS\n"
)

cat(
  "============================================================\n"
)

cat(
  "正式总行数 =",
  format(
    nrow(rad),
    big.mark = ","
  ),
  "\n"
)

cat(
  "正式RadNet站数 =",
  uniqueN(
    rad$loc_num
  ),
  "\n"
)

cat(
  "CONUS站数 =",
  uniqueN(
    conus$loc_num
  ),
  "\n"
)

cat(
  "CONUS quasi-daily行数 =",
  format(
    nrow(conus),
    big.mark = ","
  ),
  "\n"
)

cat(
  "\n正式输出文件：\n",
  OUT_FILE,
  "\n",
  sep = ""
)

cat(
  "\n统一合并阶段没有删除任何RadNet observation。\n"
)

cat(
  "后续模型筛选将在模型构建阶段单独执行。\n"
)

cat(
  "\n运行结束。\n"
)
