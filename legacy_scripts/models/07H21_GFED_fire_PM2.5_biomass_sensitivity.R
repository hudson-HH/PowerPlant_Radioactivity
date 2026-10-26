# ============================================================
# 07H21
# GFED fire-sourced PM2.5 与 Biomass–gross beta 敏感性分析
#
# 目的：
# 1. 从2001–2023 GFED fire-PM2.5 NetCDF年度压缩包中
#    提取RadNet监测站逐日fire-sourced PM2.5
# 2. 按RadNet filter period [collect_start, collect_end)
#    计算火源PM2.5平均浓度
# 3. 严格复现07H16 Figure 4 Biomass正式模型
# 4. 比较：
#       M0 = 原07H16模型
#       M1 = 原07H16模型 + GFED fire-sourced PM2.5
# 5. M0/M1强制使用完全相同的complete-case样本
#
# 正式范围：
# - branch: filter-period
# - fuel: biomass
# - direction: 90° upwind
# - metric: Count / Heat
# - radius: 20,25,30,35,40,45,50 km
# - FE: monitor × calendar year
# - clustered SE: loc_num
# - effect evaluated at each model sample median wind speed
# - study period for this sensitivity analysis: 2001–2023
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

required_packages <- c(
  "data.table",
  "ncdf4",
  "fixest"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  install.packages(missing_packages)
}

suppressPackageStartupMessages({
  library(data.table)
  library(ncdf4)
  library(fixest)
})


# ============================================================
# 1. Paths
# ============================================================

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

PP_ROOT <- file.path(
  ROOT,
  "发电站数据处理"
)

HROOT <- file.path(
  PP_ROOT,
  "06_07H_燃料特异CAMPD逐日暴露与模型"
)

OUT <- file.path(
  HROOT,
  "07H21_GFED火源PM2.5与Biomass敏感性分析"
)

dir.create(
  OUT,
  recursive = TRUE,
  showWarnings = FALSE
)

DAILY_FIRE_OUT <- file.path(
  OUT,
  "00_GFED逐日RadNet匹配"
)

dir.create(
  DAILY_FIRE_OUT,
  recursive = TRUE,
  showWarnings = FALSE
)


# ---------- GFED ----------
FIRE_DIR <- "/Volumes/NO NAME/野火"

YEARS_FIRE <- 2001:2023

FIRE_ARCHIVES <- file.path(
  FIRE_DIR,
  paste0(
    YEARS_FIRE,
    "gfedfirepm25.nc.tar.gz"
  )
)


# ---------- 07F filter-period ----------
FP_FILE <- file.path(
  PP_ROOT,
  "04_07F_正式统计分析数据构建",
  "03_07F_filter_period_analysis_base_2001_2024.csv.gz"
)


# ---------- 07F quasi-daily yearly ----------
# 用于严格恢复133个正式eligible monitors
QDIR <- file.path(
  PP_ROOT,
  "04_07F_正式统计分析数据构建",
  "quasi_daily_yearly"
)


# ---------- 07H fuel-specific daily exposure ----------
FUEL_DAILY_DIR <- file.path(
  HROOT,
  "fuel_daily_yearly"
)


# ---------- 07A正式冻结monitor坐标 ----------
# 只读取已经冻结、后续07H正式使用的07A输入
GEOM_FILE <- file.path(
  PP_ROOT,
  "00_既往发电站空间暴露与模型流程_冻结",
  "02_07A_RadNet发电站空间链接",
  "07A_radnet_powerplant_same_year_250km",
  "03_monitor_year_all_plant_link_250km_07A_2001_2024.csv"
)

if (!file.exists(GEOM_FILE)) {
  stop(
    paste0(
      "找不到正式冻结07A geometry文件：\n",
      GEOM_FILE
    )
  )
}

cat(
  "07A正式冻结geometry：PASS\n",
  GEOM_FILE,
  "\n",
  sep = ""
)


# ============================================================
# 2. Basic input audit
# ============================================================

cat(
  "\n============================================================\n",
  "07H21 输入文件检查\n",
  "============================================================\n",
  sep = ""
)

if (!file.exists(FP_FILE)) {
  stop(
    "找不到正式07F filter-period文件：\n",
    FP_FILE
  )
}

if (!dir.exists(QDIR)) {
  stop(
    "找不到07F quasi-daily yearly目录：\n",
    QDIR
  )
}

if (!dir.exists(FUEL_DAILY_DIR)) {
  stop(
    "找不到07H fuel_daily_yearly目录：\n",
    FUEL_DAILY_DIR
  )
}

missing_fire <- FIRE_ARCHIVES[
  !file.exists(FIRE_ARCHIVES)
]

if (length(missing_fire) > 0L) {
  stop(
    paste0(
      "以下GFED压缩包不存在：\n",
      paste(
        missing_fire,
        collapse = "\n"
      )
    )
  )
}

cat(
  "GFED archives = ",
  length(FIRE_ARCHIVES),
  "/23\n",
  sep = ""
)

cat(
  "07F filter-period：PASS\n",
  "07F quasi-daily yearly：PASS\n",
  "07H fuel_daily_yearly：PASS\n",
  "07A geometry：PASS\n",
  sep = ""
)


# ============================================================
# 3. Helper functions
# ============================================================

as_bool <- function(x) {
  
  if (is.logical(x)) {
    return(x)
  }
  
  z <- toupper(
    trimws(
      as.character(x)
    )
  )
  
  z %in% c(
    "TRUE",
    "T",
    "1",
    "YES",
    "Y"
  )
}


strict_mean <- function(x) {
  
  if (
    length(x) == 0L ||
    any(!is.finite(x))
  ) {
    return(NA_real_)
  }
  
  mean(x)
}


# ------------------------------------------------------------
# NetCDF time decoder
# ------------------------------------------------------------

decode_nc_time <- function(
    vals,
    units,
    year
) {
  
  units <- as.character(units)
  
  if (
    length(units) == 1L &&
    !is.na(units) &&
    nzchar(units)
  ) {
    
    u <- trimws(units)
    
    unit_type <- NA_character_
    
    if (grepl(
      "^days since",
      u,
      ignore.case = TRUE
    )) {
      unit_type <- "days"
    }
    
    if (grepl(
      "^hours since",
      u,
      ignore.case = TRUE
    )) {
      unit_type <- "hours"
    }
    
    if (grepl(
      "^minutes since",
      u,
      ignore.case = TRUE
    )) {
      unit_type <- "minutes"
    }
    
    if (grepl(
      "^seconds since",
      u,
      ignore.case = TRUE
    )) {
      unit_type <- "seconds"
    }
    
    if (!is.na(unit_type)) {
      
      origin_text <- sub(
        "^[A-Za-z]+ since[[:space:]]+",
        "",
        u
      )
      
      origin_text <- sub(
        "[[:space:]]+UTC$",
        "",
        origin_text,
        ignore.case = TRUE
      )
      
      origin <- as.POSIXct(
        origin_text,
        tz = "UTC",
        tryFormats = c(
          "%Y-%m-%d %H:%M:%S",
          "%Y-%m-%d %H:%M",
          "%Y-%m-%d"
        )
      )
      
      if (!is.na(origin)) {
        
        multiplier <- switch(
          unit_type,
          days = 86400,
          hours = 3600,
          minutes = 60,
          seconds = 1
        )
        
        d <- as.Date(
          origin +
            as.numeric(vals) *
            multiplier,
          tz = "UTC"
        )
        
        if (
          length(d) ==
          length(vals) &&
          all(!is.na(d))
        ) {
          return(d)
        }
      }
    }
  }
  
  
  # ---------- fallback ----------
  first_day <- as.Date(
    sprintf(
      "%d-01-01",
      year
    )
  )
  
  last_day <- as.Date(
    sprintf(
      "%d-12-31",
      year
    )
  )
  
  expected <- seq(
    first_day,
    last_day,
    by = "day"
  )
  
  if (
    length(vals) ==
    length(expected)
  ) {
    return(expected)
  }
  
  stop(
    paste0(
      "无法解析NetCDF time：year=",
      year,
      ", n_time=",
      length(vals),
      ", units=",
      units
    )
  )
}


# ------------------------------------------------------------
# PM2.5 unit normalization to µg/m3
# ------------------------------------------------------------

normalize_pm_units <- function(
    x,
    units
) {
  
  if (
    length(units) == 0L ||
    is.na(units) ||
    !nzchar(units)
  ) {
    return(x)
  }
  
  u <- tolower(
    as.character(units)
  )
  
  u <- gsub(
    "µ|μ",
    "u",
    u
  )
  
  # kg/m3 -> µg/m3
  if (grepl("kg", u)) {
    return(
      x * 1e9
    )
  }
  
  # mg/m3 -> µg/m3
  if (grepl("mg", u)) {
    return(
      x * 1e3
    )
  }
  
  # µg/m3 / ug/m3 / microgram
  if (
    grepl("ug", u) ||
    grepl("micro", u)
  ) {
    return(x)
  }
  
  x
}


# ------------------------------------------------------------
# Automatically identify PM2.5 variable and dimensions
# ------------------------------------------------------------

inspect_nc_variable <- function(nc) {
  
  vars <- names(
    nc$var
  )
  
  if (length(vars) == 0L) {
    stop(
      "NetCDF中没有可读取变量。"
    )
  }
  
  INFO <- rbindlist(
    lapply(
      vars,
      function(v) {
        
        dims <- vapply(
          nc$var[[v]]$dim,
          function(x) x$name,
          character(1)
        )
        
        data.table(
          variable = v,
          dimensions = paste(
            dims,
            collapse = ","
          ),
          ndims = length(dims)
        )
      }
    )
  )
  
  cat(
    "\nNetCDF variables:\n"
  )
  
  print(INFO)
  
  
  score_variable <- function(v) {
    
    obj <- nc$var[[v]]
    
    dn <- tolower(
      vapply(
        obj$dim,
        function(x) x$name,
        character(1)
      )
    )
    
    has_lon <- any(
      grepl(
        "lon|longitude",
        dn
      )
    )
    
    has_lat <- any(
      grepl(
        "lat|latitude",
        dn
      )
    )
    
    has_time <- any(
      grepl(
        "time|date|day",
        dn
      )
    )
    
    if (
      !has_lon ||
      !has_lat ||
      !has_time
    ) {
      return(-Inf)
    }
    
    nm <- tolower(v)
    
    s <- 1
    
    if (grepl(
      "pm2.?5|pm25",
      nm
    )) {
      s <- s + 10
    }
    
    if (grepl(
      "fire",
      nm
    )) {
      s <- s + 5
    }
    
    if (length(dn) == 3L) {
      s <- s + 2
    }
    
    s
  }
  
  
  scores <- vapply(
    vars,
    score_variable,
    numeric(1)
  )
  
  if (
    all(!is.finite(scores))
  ) {
    stop(
      paste0(
        "无法自动识别同时包含lon/lat/time维度的PM2.5变量。\n",
        "上方已打印NetCDF变量结构。"
      )
    )
  }
  
  vname <- vars[
    which.max(scores)
  ]
  
  obj <- nc$var[[vname]]
  
  dn <- tolower(
    vapply(
      obj$dim,
      function(x) x$name,
      character(1)
    )
  )
  
  lon_pos <- grep(
    "lon|longitude",
    dn
  )[1]
  
  lat_pos <- grep(
    "lat|latitude",
    dn
  )[1]
  
  time_pos <- grep(
    "time|date|day",
    dn
  )[1]
  
  if (
    any(
      is.na(
        c(
          lon_pos,
          lat_pos,
          time_pos
        )
      )
    )
  ) {
    stop(
      "已找到变量，但无法识别lon/lat/time位置。"
    )
  }
  
  extra_pos <- setdiff(
    seq_along(dn),
    c(
      lon_pos,
      lat_pos,
      time_pos
    )
  )
  
  if (length(extra_pos) > 0L) {
    
    extra_len <- vapply(
      obj$dim[extra_pos],
      function(z) length(z$vals),
      integer(1)
    )
    
    if (any(extra_len != 1L)) {
      stop(
        paste0(
          "目标变量存在非单例额外维度：",
          paste(
            dn[extra_pos],
            collapse = ", "
          )
        )
      )
    }
  }
  
  list(
    variable = vname,
    object = obj,
    dim_names = dn,
    lon_pos = lon_pos,
    lat_pos = lat_pos,
    time_pos = time_pos
  )
}


# ============================================================
# 4. Recover 133 formal eligible monitors
# ============================================================

cat(
  "\n============================================================\n",
  "恢复07F正式133个eligible monitors\n",
  "============================================================\n",
  sep = ""
)

Q_FILES <- file.path(
  QDIR,
  paste0(
    "07F_quasi_daily_analysis_base_",
    2001:2024,
    ".csv.gz"
  )
)

if (any(!file.exists(Q_FILES))) {
  stop(
    "07F quasi-daily yearly文件不完整。"
  )
}

ELIGIBLE_LOCS <- sort(
  unique(
    unlist(
      lapply(
        Q_FILES,
        function(f) {
          
          d <- fread(
            f,
            select = c(
              "loc_num",
              "flag_07f_primary_quasi_sample"
            )
          )
          
          d[
            as_bool(
              flag_07f_primary_quasi_sample
            ),
            as.integer(
              loc_num
            )
          ]
        }
      )
    )
  )
)

ELIGIBLE_LOCS <- ELIGIBLE_LOCS[
  !is.na(
    ELIGIBLE_LOCS
  )
]

cat(
  "eligible monitors = ",
  length(ELIGIBLE_LOCS),
  "\n",
  sep = ""
)

if (
  length(ELIGIBLE_LOCS) !=
  133L
) {
  stop(
    paste0(
      "正式eligible monitor应为133，实际=",
      length(ELIGIBLE_LOCS)
    )
  )
}


# ============================================================
# 5. Reproduce formal 107,407-row filter-period base
# ============================================================

cat(
  "\n============================================================\n",
  "恢复07G/07H正式filter-period分析样本\n",
  "============================================================\n",
  sep = ""
)

FP <- fread(
  FP_FILE
)

FP[
  ,
  collect_start :=
    as.IDate(
      collect_start
    )
]

FP[
  ,
  collect_end :=
    as.IDate(
      collect_end
    )
]

PRIMARY_FP <-
  "fpmean_campd_upwind_heat_observed_20km_90deg"

required_fp <- c(
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
  PRIMARY_FP
)

missing_fp <- setdiff(
  required_fp,
  names(FP)
)

if (length(missing_fp) > 0L) {
  stop(
    paste0(
      "07F filter-period缺少字段：\n",
      paste(
        missing_fp,
        collapse = "\n"
      )
    )
  )
}

FP[
  ,
  flag_conus_bool :=
    as_bool(
      flag_conus
    )
]

FORMAL_FP <- FP[
  flag_conus_bool &
    loc_num %in% ELIGIBLE_LOCS &
    is.finite(gross_beta_bq_m3) &
    gross_beta_bq_m3 < 0.04 &
    is.finite(pm_filter_period_official) &
    pm_filter_period_official < 35 &
    is.finite(vel) &
    vel < 6 &
    is.finite(hpbl) &
    hpbl > 150 &
    hpbl < 3000 &
    is.finite(c_prop) &
    is.finite(soilm) &
    is.finite(
      get(PRIMARY_FP)
    )
]

cat(
  "formal filter-period rows = ",
  nrow(FORMAL_FP),
  "\n",
  sep = ""
)

cat(
  "formal monitors = ",
  uniqueN(
    FORMAL_FP$loc_num
  ),
  "\n",
  sep = ""
)

if (
  nrow(FORMAL_FP) !=
  107407L
) {
  stop(
    paste0(
      "正式filter-period应为107407行，实际=",
      nrow(FORMAL_FP)
    )
  )
}

if (
  uniqueN(
    FORMAL_FP$loc_num
  ) !=
  133L
) {
  stop(
    "正式filter-period monitor数不是133。"
  )
}


# ============================================================
# 6. Restrict fire sensitivity analysis to 2001–2023
# ============================================================

FORMAL_FP[
  ,
  calendar_year :=
    as.integer(
      format(
        collect_start,
        "%Y"
      )
    )
]

FP23 <- copy(
  FORMAL_FP[
    calendar_year >= 2001 &
      calendar_year <= 2023
  ]
)

FP23[
  ,
  period_id :=
    .I
]

cat(
  "\n2001–2023 filter-period rows = ",
  nrow(FP23),
  "\n",
  sep = ""
)

cat(
  "2001–2023 monitors = ",
  uniqueN(
    FP23$loc_num
  ),
  "\n",
  sep = ""
)

# ============================================================
# 7. Read formal RadNet coordinates
#
# 正式逻辑与既有07H一致：
# 1. loc_num <-> radnet_monitor_id 从冻结07F逐年文件恢复
# 2. monitor坐标从冻结07A monitor-year基础表读取
# 3. 最终只保留正式133个eligible monitors
# ============================================================

cat(
  "\n============================================================\n",
  "读取RadNet正式经纬度\n",
  "============================================================\n",
  sep = ""
)


# ------------------------------------------------------------
# 7A. 从冻结07F建立 loc_num <-> radnet_monitor_id crosswalk
# ------------------------------------------------------------

CROSS_PARTS <- lapply(
  Q_FILES,
  function(f) {
    
    H <- names(
      fread(
        f,
        nrows = 0,
        showProgress = FALSE
      )
    )
    
    miss <- setdiff(
      c(
        "loc_num",
        "radnet_monitor_id"
      ),
      H
    )
    
    if (length(miss) > 0L) {
      stop(
        paste0(
          basename(f),
          " 缺少字段：",
          paste(
            miss,
            collapse = ", "
          )
        )
      )
    }
    
    d <- fread(
      f,
      select = c(
        "loc_num",
        "radnet_monitor_id"
      ),
      showProgress = FALSE
    )
    
    d[
      ,
      `:=`(
        loc_num =
          as.integer(
            loc_num
          ),
        
        radnet_monitor_id =
          trimws(
            as.character(
              radnet_monitor_id
            )
          )
      )
    ]
    
    unique(
      d[
        !is.na(loc_num) &
          !is.na(radnet_monitor_id) &
          nzchar(radnet_monitor_id),
        .(
          loc_num,
          radnet_monitor_id
        )
      ]
    )
  }
)


CROSS <- unique(
  rbindlist(
    CROSS_PARTS,
    use.names = TRUE,
    fill = TRUE
  )
)


DUP_LOC <- CROSS[
  ,
  .N,
  by = loc_num
][
  N > 1L
]

DUP_ID <- CROSS[
  ,
  .N,
  by = radnet_monitor_id
][
  N > 1L
]


if (
  nrow(DUP_LOC) > 0L ||
  nrow(DUP_ID) > 0L
) {
  stop(
    "07F loc_num <-> radnet_monitor_id不是一一对应。"
  )
}


cat(
  "07F monitor crosswalk = ",
  nrow(CROSS),
  "\n",
  sep = ""
)


if (
  nrow(CROSS) != 181L
) {
  stop(
    paste0(
      "07F monitor crosswalk应为181，实际=",
      nrow(CROSS)
    )
  )
}


# ------------------------------------------------------------
# 7B. 从正式冻结07A monitor-year基础表读取坐标
# ------------------------------------------------------------

COORD_FILE <- file.path(
  dirname(GEOM_FILE),
  "01_radnet_monitor_year_analysis_base_07A_2001_2024.csv"
)


if (!file.exists(COORD_FILE)) {
  stop(
    paste0(
      "找不到正式07A monitor-year基础表：\n",
      COORD_FILE
    )
  )
}


COORD_HEADER <- names(
  fread(
    COORD_FILE,
    nrows = 0,
    showProgress = FALSE
  )
)


required_coord <- c(
  "radnet_monitor_id",
  "monitor_latitude",
  "monitor_longitude"
)


missing_coord <- setdiff(
  required_coord,
  COORD_HEADER
)


if (length(missing_coord) > 0L) {
  stop(
    paste0(
      "07A monitor-year基础表缺少字段：\n",
      paste(
        missing_coord,
        collapse = "\n"
      )
    )
  )
}


COORD0 <- fread(
  COORD_FILE,
  select = required_coord,
  showProgress = FALSE
)


COORD0[
  ,
  `:=`(
    radnet_monitor_id =
      trimws(
        as.character(
          radnet_monitor_id
        )
      ),
    
    monitor_latitude =
      as.numeric(
        monitor_latitude
      ),
    
    monitor_longitude =
      as.numeric(
        monitor_longitude
      )
  )
]


COORD0 <- COORD0[
  !is.na(radnet_monitor_id) &
    nzchar(radnet_monitor_id) &
    is.finite(monitor_latitude) &
    is.finite(monitor_longitude)
]


COORD_CHECK <- COORD0[
  ,
  .(
    n_lat =
      uniqueN(
        round(
          monitor_latitude,
          6
        )
      ),
    
    n_lon =
      uniqueN(
        round(
          monitor_longitude,
          6
        )
      )
  ),
  by = radnet_monitor_id
]


BAD_COORD <- COORD_CHECK[
  n_lat > 1L |
    n_lon > 1L
]


if (nrow(BAD_COORD) > 0L) {
  stop(
    paste0(
      "发现同一个radnet_monitor_id对应多个不同坐标，站数=",
      nrow(BAD_COORD)
    )
  )
}


COORD <- COORD0[
  ,
  .(
    monitor_latitude =
      first(
        monitor_latitude
      ),
    
    monitor_longitude =
      first(
        monitor_longitude
      )
  ),
  by = radnet_monitor_id
]


cat(
  "07A coordinate monitors = ",
  nrow(COORD),
  "\n",
  sep = ""
)


# ------------------------------------------------------------
# 7C. crosswalk + coordinate
# ------------------------------------------------------------

MON <- merge(
  CROSS,
  COORD,
  by = "radnet_monitor_id",
  all.x = TRUE,
  sort = FALSE
)


MON <- MON[
  loc_num %in%
    ELIGIBLE_LOCS
]


if (
  any(
    !is.finite(
      MON$monitor_latitude
    )
  ) ||
  any(
    !is.finite(
      MON$monitor_longitude
    )
  )
) {
  
  BAD <- MON[
    !is.finite(
      monitor_latitude
    ) |
      !is.finite(
        monitor_longitude
      )
  ]
  
  print(BAD)
  
  stop(
    "正式eligible monitor中存在无法获得07A坐标的站点。"
  )
}


MON <- unique(
  MON[
    ,
    .(
      loc_num,
      radnet_monitor_id,
      monitor_latitude,
      monitor_longitude
    )
  ],
  by = "loc_num"
)


cat(
  "RadNet coordinate monitors = ",
  nrow(MON),
  "\n",
  sep = ""
)


if (
  nrow(MON) != 133L
) {
  stop(
    paste0(
      "正式RadNet coordinate monitors应为133，实际=",
      nrow(MON)
    )
  )
}


if (
  !setequal(
    MON$loc_num,
    ELIGIBLE_LOCS
  )
) {
  stop(
    "坐标表与正式133个eligible monitors不完全一致。"
  )
}


cat(
  "正式133个RadNet坐标恢复：PASS\n"
)


# ============================================================
# 8. Extract GFED daily fire-PM2.5 at each RadNet monitor
# ============================================================

process_gfed_year <- function(
    year,
    archive,
    MON
) {
  
  outfile <- file.path(
    DAILY_FIRE_OUT,
    paste0(
      "GFED_fire_PM25_RadNet_daily_",
      year,
      ".csv.gz"
    )
  )
  
  # ---------- checkpoint ----------
  if (file.exists(outfile)) {
    
    cat(
      year,
      "：已有逐日结果，跳过NetCDF重新读取。\n"
    )
    
    return(outfile)
  }
  
  
  cat(
    "\n------------------------------------------------------------\n",
    "GFED year ",
    year,
    "\n",
    "------------------------------------------------------------\n",
    sep = ""
  )
  
  
  tmpdir_year <- file.path(
    tempdir(),
    paste0(
      "GFED_",
      year
    )
  )
  
  if (dir.exists(tmpdir_year)) {
    unlink(
      tmpdir_year,
      recursive = TRUE,
      force = TRUE
    )
  }
  
  dir.create(
    tmpdir_year,
    recursive = TRUE
  )
  
  
  # ---------- extract tar.gz ----------
  utils::untar(
    archive,
    exdir = tmpdir_year
  )
  
  nc_files <- list.files(
    tmpdir_year,
    pattern = "\\.nc$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(nc_files) == 0L) {
    stop(
      paste0(
        year,
        "压缩包解压后没有.nc文件。"
      )
    )
  }
  
  # 如果archive里有多个nc，优先最大文件
  sizes <- file.info(
    nc_files
  )$size
  
  nc_file <- nc_files[
    which.max(sizes)
  ]
  
  cat(
    "NetCDF = ",
    nc_file,
    "\n",
    sep = ""
  )
  
  
  nc <- nc_open(
    nc_file
  )
  
  on.exit(
    {
      try(
        nc_close(nc),
        silent = TRUE
      )
      
      unlink(
        tmpdir_year,
        recursive = TRUE,
        force = TRUE
      )
    },
    add = TRUE
  )
  
  
  meta <- inspect_nc_variable(
    nc
  )
  
  vname <- meta$variable
  obj <- meta$object
  
  lon_pos <- meta$lon_pos
  lat_pos <- meta$lat_pos
  time_pos <- meta$time_pos
  
  
  lon_vals <- as.numeric(
    obj$dim[[lon_pos]]$vals
  )
  
  lat_vals <- as.numeric(
    obj$dim[[lat_pos]]$vals
  )
  
  time_vals <- as.numeric(
    obj$dim[[time_pos]]$vals
  )
  
  time_units <- obj$dim[[time_pos]]$units
  
  dates <- decode_nc_time(
    time_vals,
    time_units,
    year
  )
  
  
  cat(
    "Selected variable = ",
    vname,
    "\n",
    sep = ""
  )
  
  cat(
    "Variable units = ",
    obj$units,
    "\n",
    sep = ""
  )
  
  cat(
    "lon n = ",
    length(lon_vals),
    ", lat n = ",
    length(lat_vals),
    ", time n = ",
    length(time_vals),
    "\n",
    sep = ""
  )
  
  
  # ---------- map monitor longitude ----------
  M <- copy(MON)
  
  M[
    ,
    nc_lon :=
      monitor_longitude
  ]
  
  if (
    min(
      lon_vals,
      na.rm = TRUE
    ) >= 0 &&
    max(
      lon_vals,
      na.rm = TRUE
    ) > 180
  ) {
    
    M[
      monitor_longitude < 0,
      nc_lon :=
        monitor_longitude +
        360
    ]
  }
  
  
  M[
    ,
    lon_index :=
      vapply(
        nc_lon,
        function(z) {
          which.min(
            abs(
              lon_vals -
                z
            )
          )
        },
        integer(1)
      )
  ]
  
  M[
    ,
    lat_index :=
      vapply(
        monitor_latitude,
        function(z) {
          which.min(
            abs(
              lat_vals -
                z
            )
          )
        },
        integer(1)
      )
  ]
  
  M[
    ,
    gfed_grid_lon :=
      lon_vals[
        lon_index
      ]
  ]
  
  M[
    ,
    gfed_grid_lat :=
      lat_vals[
        lat_index
      ]
  ]
  
  
  # ---------- unique grid cells ----------
  CELLS <- unique(
    M[
      ,
      .(
        lon_index,
        lat_index
      )
    ]
  )
  
  
  read_one_cell <- function(
    ix,
    iy
  ) {
    
    nd <- length(
      obj$dim
    )
    
    start <- rep(
      1L,
      nd
    )
    
    count <- rep(
      1L,
      nd
    )
    
    start[
      lon_pos
    ] <- ix
    
    start[
      lat_pos
    ] <- iy
    
    start[
      time_pos
    ] <- 1L
    
    count[
      time_pos
    ] <- -1L
    
    
    x <- ncvar_get(
      nc,
      vname,
      start = start,
      count = count,
      collapse_degen = TRUE
    )
    
    x <- as.numeric(x)
    
    x[
      abs(x) >
        1e30
    ] <- NA_real_
    
    x <- normalize_pm_units(
      x,
      obj$units
    )
    
    if (
      length(x) !=
      length(dates)
    ) {
      stop(
        paste0(
          "GFED ",
          year,
          " grid cell time length不一致。"
        )
      )
    }
    
    data.table(
      lon_index = ix,
      lat_index = iy,
      date = as.IDate(dates),
      fire_pm25_gfed =
        x
    )
  }
  
  
  CELL_TS <- rbindlist(
    lapply(
      seq_len(
        nrow(CELLS)
      ),
      function(i) {
        
        read_one_cell(
          CELLS$lon_index[i],
          CELLS$lat_index[i]
        )
      }
    ),
    use.names = TRUE
  )
  
  
  OUT_Y <- merge(
    M,
    CELL_TS,
    by = c(
      "lon_index",
      "lat_index"
    ),
    allow.cartesian = TRUE,
    sort = FALSE
  )
  
  
  OUT_Y <- OUT_Y[
    as.integer(
      format(
        as.Date(date),
        "%Y"
      )
    ) == year
  ]
  
  
  keep_cols <- c(
    "loc_num",
    intersect(
      "radnet_monitor_id",
      names(OUT_Y)
    ),
    "date",
    "monitor_latitude",
    "monitor_longitude",
    "gfed_grid_lat",
    "gfed_grid_lon",
    "fire_pm25_gfed"
  )
  
  OUT_Y <- OUT_Y[
    ,
    ..keep_cols
  ]
  
  
  # ---------- hard checks ----------
  if (
    OUT_Y[
      ,
      anyDuplicated(
        paste(
          loc_num,
          date,
          sep = "_"
        )
      )
    ] > 0L
  ) {
    stop(
      paste0(
        year,
        " GFED RadNet daily存在重复loc_num+date。"
      )
    )
  }
  
  
  fwrite(
    OUT_Y,
    outfile,
    compress = "gzip"
  )
  
  cat(
    year,
    "：rows = ",
    nrow(OUT_Y),
    ", monitors = ",
    uniqueN(
      OUT_Y$loc_num
    ),
    "\n",
    sep = ""
  )
  
  outfile
}


GFED_DAILY_FILES <- vapply(
  seq_along(
    YEARS_FIRE
  ),
  function(i) {
    
    process_gfed_year(
      YEARS_FIRE[i],
      FIRE_ARCHIVES[i],
      MON
    )
  },
  character(1)
)


# ============================================================
# 9. Combine 2001–2023 GFED daily RadNet data
# ============================================================

cat(
  "\n============================================================\n",
  "合并GFED 2001–2023逐日RadNet数据\n",
  "============================================================\n",
  sep = ""
)

GFED_DAILY <- rbindlist(
  lapply(
    GFED_DAILY_FILES,
    fread
  ),
  use.names = TRUE,
  fill = TRUE
)

GFED_DAILY[
  ,
  date :=
    as.IDate(
      date
    )
]

setkey(
  GFED_DAILY,
  loc_num,
  date
)

if (
  GFED_DAILY[
    ,
    anyDuplicated(
      paste(
        loc_num,
        date,
        sep = "_"
      )
    )
  ] > 0L
) {
  stop(
    "GFED combined daily存在重复loc_num+date。"
  )
}

cat(
  "GFED daily rows = ",
  nrow(GFED_DAILY),
  "\n",
  sep = ""
)

cat(
  "GFED daily monitors = ",
  uniqueN(
    GFED_DAILY$loc_num
  ),
  "\n",
  sep = ""
)

cat(
  "GFED date range = ",
  as.character(
    min(
      GFED_DAILY$date
    )
  ),
  " ~ ",
  as.character(
    max(
      GFED_DAILY$date
    )
  ),
  "\n",
  sep = ""
)


# ============================================================
# 10. Reconstruct Biomass filter-period exposure
#     exactly from frozen fuel_daily_yearly
# ============================================================

cat(
  "\n============================================================\n",
  "重构Biomass 20–50 km filter-period暴露\n",
  "============================================================\n",
  sep = ""
)

RADII <- seq(
  20,
  50,
  by = 5
)

bio_count_raw <- paste0(
  "campd_biomass_upwind_n_plants_",
  RADII,
  "km_90deg"
)

bio_heat_raw <- paste0(
  "campd_biomass_upwind_heat_observed_",
  RADII,
  "km_90deg"
)

BIO_RAW_COLS <- c(
  bio_count_raw,
  bio_heat_raw
)


fuel_daily_files <- file.path(
  FUEL_DAILY_DIR,
  paste0(
    "07H_fuel_daily_exposure_",
    2001:2024,
    ".csv.gz"
  )
)

if (
  any(
    !file.exists(
      fuel_daily_files
    )
  )
) {
  stop(
    paste0(
      "fuel_daily_yearly缺少年度文件：\n",
      paste(
        fuel_daily_files[
          !file.exists(
            fuel_daily_files
          )
        ],
        collapse = "\n"
      )
    )
  )
}


read_bio_daily <- function(f) {
  
  header <- names(
    fread(
      f,
      nrows = 0
    )
  )
  
  need <- c(
    "loc_num",
    "date",
    BIO_RAW_COLS
  )
  
  miss <- setdiff(
    need,
    header
  )
  
  if (
    length(miss) > 0L
  ) {
    stop(
      paste0(
        basename(f),
        " 缺少Biomass字段：\n",
        paste(
          miss,
          collapse = "\n"
        )
      )
    )
  }
  
  d <- fread(
    f,
    select = need
  )
  
  d[
    ,
    date :=
      as.IDate(
        date
      )
  ]
  
  d
}


BIO_DAILY <- rbindlist(
  lapply(
    fuel_daily_files,
    read_bio_daily
  ),
  use.names = TRUE
)

setkey(
  BIO_DAILY,
  loc_num,
  date
)

if (
  BIO_DAILY[
    ,
    anyDuplicated(
      paste(
        loc_num,
        date,
        sep = "_"
      )
    )
  ] > 0L
) {
  stop(
    "Biomass daily exposure存在重复loc_num+date。"
  )
}


# ============================================================
# 11. Build filter-period day map
#     IMPORTANT: [collect_start, collect_end)
# ============================================================

cat(
  "\n============================================================\n",
  "构建filter-period逐日映射\n",
  "============================================================\n",
  sep = ""
)

if (
  any(
    FP23$collect_end <=
    FP23$collect_start
  )
) {
  stop(
    "存在collect_end <= collect_start。"
  )
}


DAY_MAP <- FP23[
  ,
  .(
    date =
      seq(
        collect_start,
        collect_end - 1L,
        by = "day"
      )
  ),
  by = .(
    period_id,
    loc_num
  )
]

cat(
  "period-day rows = ",
  nrow(DAY_MAP),
  "\n",
  sep = ""
)


# ============================================================
# 12. Aggregate Biomass daily exposure to filter period
# ============================================================

MAP_BIO <- merge(
  DAY_MAP,
  BIO_DAILY,
  by = c(
    "loc_num",
    "date"
  ),
  all.x = TRUE,
  sort = FALSE
)

FP_BIO <- MAP_BIO[
  ,
  lapply(
    .SD,
    strict_mean
  ),
  by = period_id,
  .SDcols = BIO_RAW_COLS
]

setnames(
  FP_BIO,
  BIO_RAW_COLS,
  paste0(
    "fpmean_",
    BIO_RAW_COLS
  )
)

FP23 <- merge(
  FP23,
  FP_BIO,
  by = "period_id",
  all.x = TRUE,
  sort = FALSE
)


# ============================================================
# 13. Aggregate GFED daily fire PM2.5 to filter period
# ============================================================

cat(
  "\n============================================================\n",
  "计算GFED filter-period fire-sourced PM2.5\n",
  "============================================================\n",
  sep = ""
)

MAP_FIRE <- merge(
  DAY_MAP,
  GFED_DAILY[
    ,
    .(
      loc_num,
      date,
      fire_pm25_gfed
    )
  ],
  by = c(
    "loc_num",
    "date"
  ),
  all.x = TRUE,
  sort = FALSE
)


FP_FIRE <- MAP_FIRE[
  ,
  .(
    fire_expected_days =
      .N,
    
    fire_observed_days =
      sum(
        is.finite(
          fire_pm25_gfed
        )
      ),
    
    fire_pm25_gfed_mean =
      if (
        all(
          is.finite(
            fire_pm25_gfed
          )
        )
      ) {
        mean(
          fire_pm25_gfed
        )
      } else {
        NA_real_
      },
    
    fire_pm25_gfed_max =
      if (
        any(
          is.finite(
            fire_pm25_gfed
          )
        )
      ) {
        max(
          fire_pm25_gfed,
          na.rm = TRUE
        )
      } else {
        NA_real_
      }
  ),
  by = period_id
]


FP23 <- merge(
  FP23,
  FP_FIRE,
  by = "period_id",
  all.x = TRUE,
  sort = FALSE
)


# ============================================================
# 14. Fire-PM2.5 QC
# ============================================================

cat(
  "\n============================================================\n",
  "GFED fire-PM2.5覆盖QC\n",
  "============================================================\n",
  sep = ""
)

FP23[
  ,
  fire_complete :=
    fire_expected_days ==
    fire_observed_days &
    is.finite(
      fire_pm25_gfed_mean
    )
]

cat(
  "2001–2023 filter periods = ",
  nrow(FP23),
  "\n",
  sep = ""
)

cat(
  "GFED complete periods = ",
  sum(
    FP23$fire_complete
  ),
  "\n",
  sep = ""
)

cat(
  "GFED incomplete periods = ",
  sum(
    !FP23$fire_complete
  ),
  "\n",
  sep = ""
)

cat(
  "GFED complete rate = ",
  sprintf(
    "%.3f%%",
    100 *
      mean(
        FP23$fire_complete
      )
  ),
  "\n",
  sep = ""
)


FIRE_QC_BY_YEAR <- FP23[
  ,
  .(
    n_periods = .N,
    n_complete =
      sum(
        fire_complete
      ),
    complete_pct =
      100 *
      mean(
        fire_complete
      ),
    mean_fire_pm25 =
      mean(
        fire_pm25_gfed_mean,
        na.rm = TRUE
      ),
    median_fire_pm25 =
      median(
        fire_pm25_gfed_mean,
        na.rm = TRUE
      ),
    p95_fire_pm25 =
      quantile(
        fire_pm25_gfed_mean,
        0.95,
        na.rm = TRUE
      ),
    max_fire_pm25 =
      max(
        fire_pm25_gfed_max,
        na.rm = TRUE
      )
  ),
  by = calendar_year
]

fwrite(
  FIRE_QC_BY_YEAR,
  file.path(
    OUT,
    "01_GFED_fire_PM25_filter_period_年度覆盖QC.csv"
  )
)


FIRE_SUMMARY <- FP23[
  is.finite(
    fire_pm25_gfed_mean
  ),
  .(
    n = .N,
    mean =
      mean(
        fire_pm25_gfed_mean
      ),
    sd =
      sd(
        fire_pm25_gfed_mean
      ),
    min =
      min(
        fire_pm25_gfed_mean
      ),
    p25 =
      quantile(
        fire_pm25_gfed_mean,
        0.25
      ),
    median =
      median(
        fire_pm25_gfed_mean
      ),
    p75 =
      quantile(
        fire_pm25_gfed_mean,
        0.75
      ),
    p95 =
      quantile(
        fire_pm25_gfed_mean,
        0.95
      ),
    p99 =
      quantile(
        fire_pm25_gfed_mean,
        0.99
      ),
    max =
      max(
        fire_pm25_gfed_mean
      ),
    pct_zero =
      100 *
      mean(
        fire_pm25_gfed_mean == 0
      )
  )
]

fwrite(
  FIRE_SUMMARY,
  file.path(
    OUT,
    "02_GFED_fire_PM25_filter_period_分布.csv"
  )
)


# total PM2.5 versus fire PM2.5
COR_DATA <- FP23[
  fire_complete &
    is.finite(
      pm_filter_period_official
    )
]

FIRE_PM_COR <- data.table(
  n = nrow(COR_DATA),
  
  pearson =
    cor(
      COR_DATA$pm_filter_period_official,
      COR_DATA$fire_pm25_gfed_mean,
      method = "pearson"
    ),
  
  spearman =
    cor(
      COR_DATA$pm_filter_period_official,
      COR_DATA$fire_pm25_gfed_mean,
      method = "spearman"
    )
)

fwrite(
  FIRE_PM_COR,
  file.path(
    OUT,
    "03_GFED_fire_PM25与官方总PM25相关性.csv"
  )
)


# ============================================================
# 15. Save merged analysis data
# ============================================================

fwrite(
  FP23,
  file.path(
    OUT,
    "04_filter_period_2001_2023_加入GFED_fire_PM25.csv.gz"
  ),
  compress = "gzip"
)


# ============================================================
# 16. Prepare exact 07H16 covariates
# ============================================================

FP23[
  ,
  pm_centered_10 :=
    pm_filter_period_official -
    10
]

FP23[
  ,
  temp_centered :=
    temp -
    273.15
]

FP23[
  ,
  monitor_year_fe :=
    interaction(
      loc_num,
      calendar_year,
      drop = TRUE,
      lex.order = TRUE
    )
]


BASE_MODEL_VARS <- c(
  "gross_beta_mBq_m3",
  "c_prop",
  "sunspot_monthly",
  "pm_centered_10",
  "inv_hpbl",
  "rhum",
  "soilm",
  "vel",
  "temp_centered",
  "lat_diff_35",
  "surface_u_ppm_article",
  "loc_num",
  "monitor_year_fe",
  "fire_pm25_gfed_mean"
)

missing_model_vars <- setdiff(
  BASE_MODEL_VARS,
  names(FP23)
)

if (
  length(
    missing_model_vars
  ) > 0L
) {
  stop(
    paste0(
      "正式模型缺少字段：\n",
      paste(
        missing_model_vars,
        collapse = "\n"
      )
    )
  )
}


# ============================================================
# 17. Formula builders
# ============================================================

BASE_RHS <- c(
  "c_prop",
  "sunspot_monthly",
  "pm_centered_10",
  "inv_hpbl",
  "rhum",
  "soilm",
  "vel",
  "temp_centered",
  "I(temp_centered^2)",
  "I(temp_centered^3)",
  "sunspot_monthly:lat_diff_35",
  "surface_u_ppm_article:soilm"
)


make_formula <- function(
    include_fire = FALSE
) {
  
  rhs <- c(
    "exposure_model",
    "exposure_model:vel",
    BASE_RHS
  )
  
  if (include_fire) {
    rhs <- c(
      rhs,
      "fire_pm25_gfed_centered"
    )
  }
  
  as.formula(
    paste0(
      "gross_beta_mBq_m3 ~ ",
      paste(
        rhs,
        collapse = " + "
      ),
      " | monitor_year_fe"
    )
  )
}


# ============================================================
# 18. Extract conditional exposure effect at median wind
# ============================================================

extract_effect <- function(
    fit,
    wind_value
) {
  
  b <- coef(fit)
  V <- vcov(fit)
  
  main_name <-
    "exposure_model"
  
  interaction_candidates <- c(
    "exposure_model:vel",
    "vel:exposure_model"
  )
  
  interaction_name <-
    interaction_candidates[
      interaction_candidates %in%
        names(b)
    ][1]
  
  if (
    !main_name %in%
    names(b) ||
    length(
      interaction_name
    ) == 0L ||
    is.na(
      interaction_name
    )
  ) {
    stop(
      "无法找到exposure或exposure×vel系数。"
    )
  }
  
  est <-
    b[
      main_name
    ] +
    wind_value *
    b[
      interaction_name
    ]
  
  var_est <-
    V[
      main_name,
      main_name
    ] +
    wind_value^2 *
    V[
      interaction_name,
      interaction_name
    ] +
    2 *
    wind_value *
    V[
      main_name,
      interaction_name
    ]
  
  se <-
    sqrt(
      max(
        var_est,
        0
      )
    )
  
  z <-
    est /
    se
  
  p <-
    2 *
    pnorm(
      -abs(z)
    )
  
  data.table(
    estimate =
      as.numeric(est),
    se =
      as.numeric(se),
    ci95_low =
      as.numeric(
        est -
          1.96 *
          se
      ),
    ci95_high =
      as.numeric(
        est +
          1.96 *
          se
      ),
    p_value =
      as.numeric(p)
  )
}


# ============================================================
# 19. Biomass Figure 4 sensitivity models
# ============================================================

cat(
  "\n============================================================\n",
  "运行Biomass Figure 4 GFED fire-PM2.5敏感性分析\n",
  "============================================================\n",
  sep = ""
)

FORMULA_M0 <- make_formula(
  FALSE
)

FORMULA_M1 <- make_formula(
  TRUE
)


RESULTS <- list()

k <- 1L


for (metric in c(
  "count",
  "heat"
)) {
  
  for (radius in RADII) {
    
    if (
      metric ==
      "count"
    ) {
      
      exposure_col <- paste0(
        "fpmean_campd_biomass_upwind_n_plants_",
        radius,
        "km_90deg"
      )
      
      exposure_unit <-
        "per 1 operating biomass plant"
      
      scale_factor <- 1
      
    } else {
      
      exposure_col <- paste0(
        "fpmean_campd_biomass_upwind_heat_observed_",
        radius,
        "km_90deg"
      )
      
      exposure_unit <-
        "per 10,000 MMBtu/day biomass heat input"
      
      scale_factor <- 10000
    }
    
    
    if (
      !exposure_col %in%
      names(FP23)
    ) {
      stop(
        paste0(
          "找不到Biomass exposure：",
          exposure_col
        )
      )
    }
    
    
    # --------------------------------------------------------
    # IMPORTANT:
    # M0和M1完全使用同一complete-case样本
    # --------------------------------------------------------
    
    MODEL_VARS <- c(
      BASE_MODEL_VARS,
      exposure_col
    )
    
    DS <- FP23[
      fire_complete &
        complete.cases(
          FP23[
            ,
            ..MODEL_VARS
          ]
        )
    ]
    
    
    DS[
      ,
      exposure_model :=
        get(
          exposure_col
        ) /
        scale_factor
    ]
    
    
    DS <- DS[
      is.finite(
        exposure_model
      )
    ]
    
    
    # fire PM中心化只用于数值稳定
    fire_center <-
      mean(
        DS$fire_pm25_gfed_mean
      )
    
    DS[
      ,
      fire_pm25_gfed_centered :=
        fire_pm25_gfed_mean -
        fire_center
    ]
    
    
    wind_median <-
      median(
        DS$vel
      )
    
    
    cat(
      metric,
      " | ",
      radius,
      " km | rows=",
      nrow(DS),
      " | monitors=",
      uniqueN(
        DS$loc_num
      ),
      " | wind median=",
      sprintf(
        "%.4f",
        wind_median
      ),
      "\n",
      sep = ""
    )
    
    
    # ---------- M0 original ----------
    FIT_M0 <- feols(
      FORMULA_M0,
      data = DS,
      vcov = ~loc_num,
      notes = FALSE,
      warn = FALSE
    )
    
    
    # ---------- M1 + fire-PM2.5 ----------
    FIT_M1 <- feols(
      FORMULA_M1,
      data = DS,
      vcov = ~loc_num,
      notes = FALSE,
      warn = FALSE
    )
    
    
    E0 <- extract_effect(
      FIT_M0,
      wind_median
    )
    
    E1 <- extract_effect(
      FIT_M1,
      wind_median
    )
    
    
    # fire-PM2.5 coefficient
    b1 <- coef(
      FIT_M1
    )
    
    v1 <- vcov(
      FIT_M1
    )
    
    fire_name <-
      "fire_pm25_gfed_centered"
    
    if (
      fire_name %in%
      names(b1)
    ) {
      
      fire_est <-
        b1[
          fire_name
        ]
      
      fire_se <-
        sqrt(
          v1[
            fire_name,
            fire_name
          ]
        )
      
      fire_low <-
        fire_est -
        1.96 *
        fire_se
      
      fire_high <-
        fire_est +
        1.96 *
        fire_se
      
      fire_p <-
        2 *
        pnorm(
          -abs(
            fire_est /
              fire_se
          )
        )
      
    } else {
      
      fire_est <-
        NA_real_
      
      fire_se <-
        NA_real_
      
      fire_low <-
        NA_real_
      
      fire_high <-
        NA_real_
      
      fire_p <-
        NA_real_
    }
    
    
    RESULTS[[k]] <- data.table(
      fuel = "biomass",
      metric = metric,
      radius_km = radius,
      exposure_column =
        exposure_col,
      exposure_unit =
        exposure_unit,
      
      n_complete_case =
        nrow(DS),
      
      n_model_M0 =
        nobs(
          FIT_M0
        ),
      
      n_model_M1 =
        nobs(
          FIT_M1
        ),
      
      n_monitors =
        uniqueN(
          DS$loc_num
        ),
      
      wind_reference =
        "sample median wind",
      
      wind_ms =
        wind_median,
      
      fire_pm25_center =
        fire_center,
      
      original_estimate =
        E0$estimate,
      
      original_ci95_low =
        E0$ci95_low,
      
      original_ci95_high =
        E0$ci95_high,
      
      original_p =
        E0$p_value,
      
      fire_adjusted_estimate =
        E1$estimate,
      
      fire_adjusted_ci95_low =
        E1$ci95_low,
      
      fire_adjusted_ci95_high =
        E1$ci95_high,
      
      fire_adjusted_p =
        E1$p_value,
      
      estimate_change =
        E1$estimate -
        E0$estimate,
      
      relative_change_pct =
        ifelse(
          abs(
            E0$estimate
          ) >
            .Machine$double.eps,
          100 *
            (
              E1$estimate -
                E0$estimate
            ) /
            E0$estimate,
          NA_real_
        ),
      
      attenuation_pct =
        ifelse(
          abs(
            E0$estimate
          ) >
            .Machine$double.eps,
          100 *
            (
              E0$estimate -
                E1$estimate
            ) /
            abs(
              E0$estimate
            ),
          NA_real_
        ),
      
      fire_pm25_coef =
        as.numeric(
          fire_est
        ),
      
      fire_pm25_ci95_low =
        as.numeric(
          fire_low
        ),
      
      fire_pm25_ci95_high =
        as.numeric(
          fire_high
        ),
      
      fire_pm25_p =
        as.numeric(
          fire_p
        )
    )
    
    k <- k + 1L
  }
}


RESULT <- rbindlist(
  RESULTS,
  use.names = TRUE
)


# ============================================================
# 20. Hard checks
# ============================================================

cat(
  "\n============================================================\n",
  "07H21模型审计\n",
  "============================================================\n",
  sep = ""
)

cat(
  "model pairs = ",
  nrow(RESULT),
  "\n",
  sep = ""
)

if (
  nrow(RESULT) !=
  14L
) {
  stop(
    "应有14组Biomass模型：2 metrics × 7 radii。"
  )
}


if (
  any(
    RESULT$n_model_M0 !=
    RESULT$n_model_M1
  )
) {
  stop(
    "M0与M1模型实际N不一致。"
  )
}


if (
  any(
    !is.finite(
      RESULT$original_estimate
    )
  ) ||
  any(
    !is.finite(
      RESULT$fire_adjusted_estimate
    )
  )
) {
  stop(
    "存在模型效应估计失败。"
  )
}


cat(
  "14/14 Biomass model pairs：PASS\n"
)

cat(
  "M0/M1 same-sample audit：PASS\n"
)


# ============================================================
# 21. Save model results
# ============================================================

fwrite(
  RESULT,
  file.path(
    OUT,
    "05_Biomass_Figure4_GFED_fire_PM25敏感性结果.csv"
  )
)


# concise comparison table
COMPARISON <- RESULT[
  ,
  .(
    metric,
    radius_km,
    n_model =
      n_model_M1,
    n_monitors,
    wind_ms,
    
    original_estimate,
    original_ci95_low,
    original_ci95_high,
    
    fire_adjusted_estimate,
    fire_adjusted_ci95_low,
    fire_adjusted_ci95_high,
    
    estimate_change,
    relative_change_pct,
    attenuation_pct,
    
    fire_pm25_coef,
    fire_pm25_ci95_low,
    fire_pm25_ci95_high,
    fire_pm25_p
  )
]

fwrite(
  COMPARISON,
  file.path(
    OUT,
    "06_Biomass加入GFED火源PM25前后比较.csv"
  )
)


# ============================================================
# 22. Optional figure
# ============================================================

if (
  requireNamespace(
    "ggplot2",
    quietly = TRUE
  )
) {
  
  library(
    ggplot2
  )
  
  PLOT_DATA <- rbind(
    RESULT[
      ,
      .(
        metric,
        radius_km,
        model =
          "Original Figure 4",
        estimate =
          original_estimate,
        ci95_low =
          original_ci95_low,
        ci95_high =
          original_ci95_high
      )
    ],
    RESULT[
      ,
      .(
        metric,
        radius_km,
        model =
          "+ GFED fire PM2.5",
        estimate =
          fire_adjusted_estimate,
        ci95_low =
          fire_adjusted_ci95_low,
        ci95_high =
          fire_adjusted_ci95_high
      )
    ]
  )
  
  
  p <- ggplot(
    PLOT_DATA,
    aes(
      x = radius_km,
      y = estimate,
      group = model,
      linetype = model,
      shape = model
    )
  ) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.4
    ) +
    geom_errorbar(
      aes(
        ymin = ci95_low,
        ymax = ci95_high
      ),
      width = 0.8,
      position =
        position_dodge(
          width = 1
        )
    ) +
    geom_line(
      position =
        position_dodge(
          width = 1
        )
    ) +
    geom_point(
      size = 2.2,
      position =
        position_dodge(
          width = 1
        )
    ) +
    facet_wrap(
      ~metric,
      scales = "free_y",
      ncol = 1
    ) +
    scale_x_continuous(
      breaks = RADII
    ) +
    labs(
      x = "Radius (km)",
      y = "Biomass association with gross beta",
      linetype = NULL,
      shape = NULL,
      title =
        "Biomass association before and after adjustment for GFED fire-sourced PM2.5"
    ) +
    theme_bw(
      base_size = 12
    ) +
    theme(
      legend.position = "top"
    )
  
  
  ggsave(
    filename =
      file.path(
        OUT,
        "07_Biomass加入GFED火源PM25前后比较.png"
      ),
    plot = p,
    width = 8,
    height = 8,
    dpi = 300
  )
}


# ============================================================
# 23. Final console summary
# ============================================================

cat(
  "\n============================================================\n",
  "07H21 GFED fire-sourced PM2.5 sensitivity：完成\n",
  "============================================================\n",
  sep = ""
)

cat(
  "\n正式分析范围：2001–2023\n"
)

cat(
  "Fuel：Biomass\n"
)

cat(
  "Direction：90° upwind\n"
)

cat(
  "Metrics：Count + Heat\n"
)

cat(
  "Radii：20–50 km, step 5 km\n"
)

cat(
  "FE：Monitor × year\n"
)

cat(
  "SE：clustered by loc_num\n"
)

cat(
  "Wind reference：each model complete-case sample median\n"
)

cat(
  "M0/M1：严格相同complete-case样本\n\n"
)


print(
  COMPARISON
)


cat(
  "\n输出目录：\n",
  OUT,
  "\n",
  sep = ""
)

cat(
  "\n核心查看文件：\n",
  "05_Biomass_Figure4_GFED_fire_PM25敏感性结果.csv\n",
  "06_Biomass加入GFED火源PM25前后比较.csv\n",
  "07_Biomass加入GFED火源PM25前后比较.png\n",
  sep = ""
)

cat(
  "\n============================================================\n",
  "运行结束\n",
  "============================================================\n",
  sep = ""
)
