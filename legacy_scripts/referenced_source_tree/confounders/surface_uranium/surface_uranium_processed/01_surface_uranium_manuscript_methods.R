# ============================================================
# RadNet 地表 U：按照参考文章公开代码重新构建
#
# 文章核心方法：
#
# radnet_basic$U <- raster::extract(
#     rad_u,
#     radnet,
#     buffer = 30000,
#     fun = median,
#     na.rm = TRUE
# )
#
# radnet_basic[is.na(radnet_basic$U),"U"] <-
#     mean(radnet_basic[!is.na(radnet_basic$U),"U"])
#
# 即：
# 1. 每个 RadNet 站点周围 30 km
# 2. U 栅格有效像元取 median
# 3. 仍缺失的站点使用所有非缺失 RadNet 站点 U 的 mean 填补
#
# 正式主键：EPA loc_num
# ============================================================


# =========================
# 0. 环境
# =========================

required_packages <- c(
  "data.table",
  "raster",
  "sp"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  install.packages(
    missing_packages,
    repos = "https://cloud.r-project.org"
  )
}

library(data.table)
library(raster)
library(sp)


# =========================
# 1. 路径
# =========================

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

U_ROOT <- file.path(
  ROOT,
  "混杂因素",
  "地表 U "
)

OUT_DIR <- file.path(
  U_ROOT,
  "地表 U处理版"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

CROSSWALK_FILE <- file.path(
  ROOT,
  "radnet命名标准化",
  "radnet_181站_loc_num_radnet_monitor_id对照表.csv"
)

U_FILE <- file.path(
  U_ROOT,
  "NAMrad_U",
  "NAMrad_U.flt"
)

OLD_U_FILE <- file.path(
  ROOT,
  "confounder_build",
  "08A1B_surface_radiometric_fill",
  "01_radnet_monitor_surface_radiometric_08A1B.csv"
)

OUT_FILE <- file.path(
  OUT_DIR,
  "radnet_181站_地表U_文章方法.csv"
)


# =========================
# 2. 检查输入
# =========================

if (!file.exists(CROSSWALK_FILE)) {
  stop(
    "找不到 RadNet 181站正式对照表：\n",
    CROSSWALK_FILE
  )
}

if (!file.exists(U_FILE)) {
  stop(
    "找不到 USGS U 栅格：\n",
    U_FILE
  )
}


# =========================
# 3. 读取181站正式对照表
# =========================

radnet <- fread(
  CROSSWALK_FILE
)

required_cols <- c(
  "loc_num",
  "radnet_monitor_id",
  "station",
  "city_name",
  "state_abbr",
  "monitor_latitude",
  "monitor_longitude"
)

missing_cols <- setdiff(
  required_cols,
  names(radnet)
)

if (length(missing_cols) > 0) {
  stop(
    "crosswalk缺少字段：",
    paste(missing_cols, collapse = ", ")
  )
}

if (nrow(radnet) != 181) {
  stop(
    "正式RadNet crosswalk不是181行。"
  )
}

if (uniqueN(radnet$loc_num) != 181) {
  stop(
    "loc_num不是181个唯一值。"
  )
}

if (
  radnet[
    is.na(monitor_latitude) |
    is.na(monitor_longitude),
    .N
  ] > 0
) {
  stop(
    "存在缺失RadNet坐标。"
  )
}


# =========================
# 4. 读取原始 USGS U 栅格
# =========================

cat(
  "\n读取 USGS U 栅格：\n",
  U_FILE,
  "\n",
  sep = ""
)

u_raster <- raster(
  U_FILE
)

cat(
  "\n================ U raster ================\n"
)

cat(
  "rows =",
  nrow(u_raster),
  "\n"
)

cat(
  "cols =",
  ncol(u_raster),
  "\n"
)

cat(
  "resolution =",
  paste(
    res(u_raster),
    collapse = " × "
  ),
  "\n"
)

cat(
  "CRS =\n",
  projection(u_raster),
  "\n"
)

cat(
  "extent =\n"
)

print(
  extent(u_raster)
)


# =========================
# 5. 定义 USGS DNAG CRS
#
# 原始 FLT/HDR 没有内嵌 CRS，
# 但 USGS OFR 2005-1413 明确说明该数据采用：
#
# Spherical Transverse Mercator
# central meridian = 100 W
# latitude of origin = 0
# scale factor = 0.926
# false easting/northing = 0
# spherical Earth radius = 6,371,204 m
# units = meters
#
# 因此这里不是“猜测 CRS”，
# 而是按照 USGS 数据说明给原始栅格补回其正式 DNAG 投影。
# =========================

DNAG_CRS <- "+proj=tmerc +lat_0=0 +lon_0=-100 +k=0.926 +x_0=0 +y_0=0 +R=6371204 +units=m +no_defs"

if (
  is.na(projection(u_raster)) ||
  projection(u_raster) == ""
) {

  projection(u_raster) <- DNAG_CRS

  cat(
    "\n原始 FLT 未内嵌 CRS；",
    "已按照 USGS 官方说明定义为 DNAG Spherical Transverse Mercator。\n",
    sep = ""
  )
}

if (
  isTRUE(
    raster::isLonLat(u_raster)
  )
) {
  stop(
    "U raster仍被识别为经纬度坐标系，停止执行。"
  )
}

cat(
  "正式使用的 U raster CRS：\n",
  projection(u_raster),
  "\n"
)


# =========================
# 6. 建立 WGS84 RadNet 空间点
# =========================

radnet_sp <- SpatialPointsDataFrame(
  coords = radnet[
    ,
    .(
      monitor_longitude,
      monitor_latitude
    )
  ],
  data = as.data.frame(radnet),
  proj4string = CRS(
    "+proj=longlat +datum=WGS84 +no_defs"
  )
)


# =========================
# 7. 转换到 USGS U raster CRS
# =========================

radnet_u_crs <- spTransform(
  radnet_sp,
  CRS(
    projection(
      u_raster
    )
  )
)


# =========================
# 8. 完全按照文章方法：
#
# 30 km buffer
# median
# na.rm = TRUE
# =========================

cat(
  "\n开始按文章方法提取：",
  "30 km buffer + median...\n",
  sep = ""
)

surface_u_30km_median <- raster::extract(
  u_raster,
  radnet_u_crs,
  buffer = 30000,
  fun = median,
  na.rm = TRUE
)


if (
  length(surface_u_30km_median) !=
  nrow(radnet)
) {
  stop(
    "U提取结果长度与181站不一致。"
  )
}


# =========================
# 9. 文章式缺失值处理
#
# 所有非缺失站点U的mean
# 用于替换仍为NA的站点
# =========================

n_available_before <- sum(
  !is.na(
    surface_u_30km_median
  )
)

n_missing_before <- sum(
  is.na(
    surface_u_30km_median
  )
)

if (n_available_before == 0) {
  stop(
    "所有RadNet站点U均为NA，说明空间提取失败。"
  )
}

surface_u_mean_nonmissing <- mean(
  surface_u_30km_median[
    !is.na(
      surface_u_30km_median
    )
  ]
)

surface_u_article <- surface_u_30km_median

surface_u_imputed_article <- is.na(
  surface_u_article
)

surface_u_article[
  surface_u_imputed_article
] <- surface_u_mean_nonmissing


# =========================
# 10. 构建正式181站结果
# =========================

result <- copy(
  radnet[
    ,
    .(
      loc_num,
      radnet_monitor_id,
      station,
      city_name,
      state_abbr,
      monitor_latitude,
      monitor_longitude
    )
  ]
)

result[
  ,
  surface_u_ppm_30km_median :=
    surface_u_30km_median
]

result[
  ,
  surface_u_ppm_article :=
    surface_u_article
]

result[
  ,
  surface_u_imputed_article :=
    surface_u_imputed_article
]

result[
  ,
  surface_u_method :=
    fifelse(
      surface_u_imputed_article,
      "mean_of_nonmissing_radnet_U_after_30km_median",
      "30km_buffer_median"
    )
]


# =========================
# 11. 硬性检查
# =========================

if (nrow(result) != 181) {
  stop(
    "正式U结果不是181行。"
  )
}

if (uniqueN(result$loc_num) != 181) {
  stop(
    "正式U结果loc_num不是181个唯一值。"
  )
}

if (
  result[
    is.na(surface_u_ppm_article),
    .N
  ] != 0
) {
  stop(
    "文章式均值填补后仍存在U缺失。"
  )
}

if (
  any(
    !is.finite(
      result$surface_u_ppm_article
    )
  )
) {
  stop(
    "正式U中存在非有限数值。"
  )
}


# =========================
# 12. 与旧08A1B结果做终端比较
#
# 仅作为验证参考；
# 旧文件不参与新U的计算。
# =========================

if (
  file.exists(
    OLD_U_FILE
  )
) {

  old_u <- fread(
    OLD_U_FILE
  )

  if (
    all(
      c(
        "radnet_monitor_id",
        "surface_u_ppm_buffer30km_median"
      ) %in%
        names(old_u)
    )
  ) {

    cmp <- merge(
      result[
        ,
        .(
          radnet_monitor_id,
          new_u =
            surface_u_ppm_30km_median
        )
      ],
      old_u[
        ,
        .(
          radnet_monitor_id,
          old_u =
            surface_u_ppm_buffer30km_median
        )
      ],
      by = "radnet_monitor_id",
      all = FALSE
    )

    both <- cmp[
      !is.na(new_u) &
      !is.na(old_u)
    ]

    cat(
      "\n================ 与旧08A1B比较，仅供验证 ================\n"
    )

    cat(
      "共同有值站数 =",
      nrow(both),
      "\n"
    )

    if (nrow(both) > 0) {

      cat(
        "mean absolute difference =",
        mean(
          abs(
            both$new_u -
            both$old_u
          )
        ),
        "ppm\n"
      )

      cat(
        "max absolute difference =",
        max(
          abs(
            both$new_u -
            both$old_u
          )
        ),
        "ppm\n"
      )
    }

    cat(
      "新结果缺失站数 =",
      cmp[
        is.na(new_u),
        .N
      ],
      "\n"
    )

    cat(
      "旧结果缺失站数 =",
      cmp[
        is.na(old_u),
        .N
      ],
      "\n"
    )
  }
}


# =========================
# 13. 写正式结果
# =========================

setorder(
  result,
  loc_num
)

fwrite(
  result,
  OUT_FILE,
  na = ""
)


# =========================
# 14. 最终终端摘要
# =========================

cat(
  "\n\n",
  "============================================================\n",
  "RadNet 地表 U：文章方法重做完成\n",
  "============================================================\n",
  sep = ""
)

cat(
  "RadNet站点数：",
  nrow(result),
  "\n",
  sep = ""
)

cat(
  "30 km median直接有值：",
  n_available_before,
  "\n",
  sep = ""
)

cat(
  "30 km median缺失：",
  n_missing_before,
  "\n",
  sep = ""
)

cat(
  "文章式均值填补值：",
  surface_u_mean_nonmissing,
  " ppm\n",
  sep = ""
)

cat(
  "最终U非缺失：",
  result[
    !is.na(
      surface_u_ppm_article
    ),
    .N
  ],
  "\n",
  sep = ""
)

cat(
  "最终U缺失：",
  result[
    is.na(
      surface_u_ppm_article
    ),
    .N
  ],
  "\n",
  sep = ""
)

cat(
  "\n被文章式均值填补的站点：\n"
)

print(
  result[
    surface_u_imputed_article == TRUE,
    .(
      loc_num,
      radnet_monitor_id,
      station,
      state_abbr,
      surface_u_ppm_30km_median,
      surface_u_ppm_article
    )
  ]
)

cat(
  "\n正式U描述：\n"
)

print(
  summary(
    result$surface_u_ppm_article
  )
)

cat(
  "\n正式输出：\n",
  OUT_FILE,
  "\n",
  sep = ""
)

cat(
  "\n总体状态：PASS\n"
)
