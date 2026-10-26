library(data.table)

# ============================================================
# HYSPLIT daily air-mass origin
#
# 不重新计算轨迹。
# 现有 HYSPLIT 已按参考文章方法构建：
# - RadNet 逐站
# - daily trajectory
# - continental trajectory fraction
#
# 本脚本只做：
# 1. 统一到新的 EPA loc_num
# 2. 保留原 continental_trajectory_fraction
# 3. 增加文章代码变量名 c_prop
# 4. 做必要硬检查
# ============================================================


# =========================
# 1. 路径
# =========================

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

HYSPLIT_ROOT <- file.path(
  ROOT,
  "混杂因素",
  "08A4_HYSPLIT_air_mass_origin"
)

OUT_DIR <- file.path(
  HYSPLIT_ROOT,
  "HYSPLIT处理版"
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

HYSPLIT_FILE <- file.path(
  HYSPLIT_ROOT,
  "05_final",
  "01_daily_continental_trajectory_fraction_2001_2024.csv.gz"
)

OUT_FILE <- file.path(
  OUT_DIR,
  "radnet_181站_HYSPLIT_daily_文章方法_2001_2024.csv.gz"
)


# =========================
# 2. 输入检查
# =========================

if (!file.exists(CROSSWALK_FILE)) {
  stop(
    "找不到 RadNet 命名标准化对照表：\n",
    CROSSWALK_FILE
  )
}

if (!file.exists(HYSPLIT_FILE)) {
  stop(
    "找不到 HYSPLIT 正式 daily 文件：\n",
    HYSPLIT_FILE
  )
}


# =========================
# 3. 读取 crosswalk
# =========================

cw <- fread(
  CROSSWALK_FILE
)

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
    paste(missing_cw, collapse = ", ")
  )
}

if (nrow(cw) != 181) {
  stop("crosswalk不是181行。")
}

if (uniqueN(cw$loc_num) != 181) {
  stop("loc_num不是181个唯一值。")
}

if (uniqueN(cw$radnet_monitor_id) != 181) {
  stop("radnet_monitor_id不是181个唯一值。")
}


# =========================
# 4. 读取 HYSPLIT daily
# =========================

h <- fread(
  HYSPLIT_FILE
)

required_h <- c(
  "radnet_monitor_id",
  "date",
  "year",
  "continental_trajectory_fraction",
  "n_arrivals_valid",
  "n_total_endpoints",
  "n_land_endpoints",
  "n_full_72h_trajectories"
)

missing_h <- setdiff(
  required_h,
  names(h)
)

if (length(missing_h) > 0) {
  stop(
    "HYSPLIT文件缺少字段：",
    paste(missing_h, collapse = ", ")
  )
}

h[, date := as.IDate(date)]


# =========================
# 5. 输入结构检查
# =========================

cat("\n================ HYSPLIT 输入 ================\n")
cat("rows =", nrow(h), "\n")
cat("monitors =", uniqueN(h$radnet_monitor_id), "\n")
cat("date =", as.character(min(h$date)), "to", as.character(max(h$date)), "\n")

if (uniqueN(h$radnet_monitor_id) != 181) {
  stop("HYSPLIT输入不是181个RadNet监测站。")
}

if (min(h$date, na.rm = TRUE) != as.IDate("2001-01-01")) {
  stop("HYSPLIT起始日期不是2001-01-01。")
}

if (max(h$date, na.rm = TRUE) != as.IDate("2024-12-31")) {
  stop("HYSPLIT结束日期不是2024-12-31。")
}

dup_old <- h[
  ,
  .N,
  by = .(
    radnet_monitor_id,
    date
  )
][N > 1]

if (nrow(dup_old) > 0) {
  print(head(dup_old, 20))
  stop("HYSPLIT输入存在重复 monitor-date。")
}


# =========================
# 6. 检查所有旧ID都能映射到 loc_num
# =========================

unmatched_old <- setdiff(
  unique(h$radnet_monitor_id),
  cw$radnet_monitor_id
)

if (length(unmatched_old) > 0) {
  cat("\n未匹配的HYSPLIT旧ID：\n")
  print(unmatched_old)
  stop("存在HYSPLIT旧ID无法映射到loc_num。")
}


# =========================
# 7. 合并新的正式 loc_num
# =========================

result <- merge(
  h,
  cw[
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
  ],
  by = "radnet_monitor_id",
  all.x = TRUE,
  sort = FALSE
)

if (result[is.na(loc_num), .N] > 0) {
  stop("合并后仍存在缺失loc_num。")
}


# =========================
# 8. 增加文章变量名 c_prop
#
# c_prop = continental trajectory fraction
# 不修改原值。
# =========================

result[
  ,
  c_prop :=
    continental_trajectory_fraction
]


# =========================
# 9. HYSPLIT值域检查
# =========================

bad_fraction <- result[
  !is.na(c_prop) &
  (
    c_prop < 0 |
    c_prop > 1
  )
]

if (nrow(bad_fraction) > 0) {
  print(head(bad_fraction, 20))
  stop("存在 c_prop 超出0–1范围。")
}

bad_arrivals <- result[
  !is.na(n_arrivals_valid) &
  (
    n_arrivals_valid < 0 |
    n_arrivals_valid > 4
  )
]

if (nrow(bad_arrivals) > 0) {
  print(head(bad_arrivals, 20))
  stop("存在 n_arrivals_valid 超出0–4范围。")
}

bad_full72 <- result[
  !is.na(n_full_72h_trajectories) &
  (
    n_full_72h_trajectories < 0 |
    n_full_72h_trajectories > 4
  )
]

if (nrow(bad_full72) > 0) {
  print(head(bad_full72, 20))
  stop("存在 n_full_72h_trajectories 超出0–4范围。")
}


# =========================
# 10. 新正式主键重复检查
# =========================

dup_new <- result[
  ,
  .N,
  by = .(
    loc_num,
    date
  )
][N > 1]

if (nrow(dup_new) > 0) {
  print(head(dup_new, 20))
  stop("转换为loc_num后出现重复loc_num-date。")
}


# =========================
# 11. 字段整理
# =========================

result <- result[
  ,
  .(
    loc_num,
    radnet_monitor_id,
    station,
    city_name,
    state_abbr,
    monitor_latitude,
    monitor_longitude,
    date,
    year,
    c_prop,
    continental_trajectory_fraction,
    n_arrivals_valid,
    n_total_endpoints,
    n_land_endpoints,
    n_full_72h_trajectories
  )
]

setorder(
  result,
  loc_num,
  date
)


# =========================
# 12. 最终硬检查
# =========================

if (uniqueN(result$loc_num) != 181) {
  stop("最终HYSPLIT不是181个loc_num。")
}

if (nrow(result) != nrow(h)) {
  stop("处理前后HYSPLIT行数发生变化。")
}

if (
  !isTRUE(
    all.equal(
      result$c_prop,
      result$continental_trajectory_fraction,
      check.attributes = FALSE
    )
  )
) {
  stop("c_prop与原continental_trajectory_fraction不一致。")
}


# =========================
# 13. 输出
# =========================

fwrite(
  result,
  OUT_FILE
)


# =========================
# 14. 终端摘要
# =========================

cat("\n\n============================================================\n")
cat("HYSPLIT 文章方法标准化完成\n")
cat("============================================================\n")

cat("rows =", nrow(result), "\n")
cat("loc_num =", uniqueN(result$loc_num), "\n")
cat("radnet_monitor_id =", uniqueN(result$radnet_monitor_id), "\n")
cat("date =", as.character(min(result$date)), "to", as.character(max(result$date)), "\n")

cat("\nc_prop非缺失 =", result[!is.na(c_prop), .N], "\n")
cat("c_prop缺失 =", result[is.na(c_prop), .N], "\n")

cat("\nn_arrivals_valid分布：\n")
print(
  result[
    ,
    .N,
    by = n_arrivals_valid
  ][
    order(n_arrivals_valid)
  ]
)

cat("\nn_full_72h_trajectories分布：\n")
print(
  result[
    ,
    .N,
    by = n_full_72h_trajectories
  ][
    order(n_full_72h_trajectories)
  ]
)

cat("\nc_prop描述：\n")
print(
  summary(result$c_prop)
)

cat("\n前5行：\n")
print(
  head(result, 5)
)

cat("\n正式输出：\n")
cat(OUT_FILE, "\n")

cat("\n总体状态：PASS\n")
