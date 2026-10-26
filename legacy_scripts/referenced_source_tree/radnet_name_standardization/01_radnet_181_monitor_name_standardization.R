library(data.table)

root <- "/Users/hudson/Desktop/USA 辐射表征研究"

out_dir <- file.path(root, "radnet命名标准化")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

epa_file <- file.path(
  root,
  "radnet 监测数据",
  "02_滤膜周期数据",
  "radnet_gross_beta_filter_period_2001_2024.csv.gz"
)

old_file <- file.path(
  root,
  "confounder_build",
  "08A0B_deep_metadata",
  "05_radnet_unique_monitor_coordinates.csv"
)

out_file <- file.path(
  out_dir,
  "radnet_181站_loc_num_radnet_monitor_id对照表.csv"
)

epa <- fread(epa_file)
old <- fread(old_file)

epa_station <- unique(
  epa[, .(
    loc_num,
    station,
    city_name,
    state_abbr
  )]
)

epa_station <- epa_station[
  ,
  .(
    station = station[1],
    city_name = city_name[1],
    state_abbr = state_abbr[1]
  ),
  by = loc_num
]

old[, old_station := sub("^.* - ", "", radnet_monitor_id)]

normalize_station <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "")
  x <- gsub("[^A-Z0-9]", "", x)
  x
}

normalize_state <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x <- gsub("[^A-Z]", "", x)
  x
}

epa_station[, match_station := normalize_station(station)]
epa_station[, match_state := normalize_state(state_abbr)]
old[, match_station := normalize_station(old_station)]
old[, match_state := normalize_state(radnet_state_abbr)]

epa_station[, match_key := paste(match_state, match_station, sep = "__")]
old[, match_key := paste(match_state, match_station, sep = "__")]

epa_dup <- epa_station[, .N, by = match_key][N > 1]
old_dup <- old[, .N, by = match_key][N > 1]

if (nrow(epa_dup) > 0) {
  print(epa_dup)
  stop("EPA标准化后存在重复匹配键。")
}

if (nrow(old_dup) > 0) {
  print(old_dup)
  stop("旧RadNet标准化后存在重复匹配键。")
}

crosswalk <- merge(
  epa_station,
  old[
    ,
    .(
      match_key,
      radnet_monitor_id,
      old_station,
      monitor_latitude,
      monitor_longitude,
      radnet_monitor_label,
      radnet_state_abbr
    )
  ],
  by = "match_key",
  all = TRUE,
  sort = FALSE
)

epa_unmatched <- crosswalk[
  !is.na(loc_num) & is.na(radnet_monitor_id)
]

old_unmatched <- crosswalk[
  is.na(loc_num) & !is.na(radnet_monitor_id)
]

cat("\n============================================================\n")
cat("RadNet 命名标准化检查\n")
cat("============================================================\n")

cat("EPA loc_num 数：", uniqueN(epa_station$loc_num), "\n")
cat("旧 radnet_monitor_id 数：", uniqueN(old$radnet_monitor_id), "\n")
cat("EPA未匹配：", nrow(epa_unmatched), "\n")
cat("旧ID未匹配：", nrow(old_unmatched), "\n")

if (uniqueN(epa_station$loc_num) != 181) stop("EPA loc_num不是181个。")
if (uniqueN(old$radnet_monitor_id) != 181) stop("旧radnet_monitor_id不是181个。")
if (nrow(epa_unmatched) != 0) stop("存在EPA loc_num未匹配。")
if (nrow(old_unmatched) != 0) stop("存在旧radnet_monitor_id未匹配。")

crosswalk <- crosswalk[
  !is.na(loc_num) & !is.na(radnet_monitor_id),
  .(
    loc_num,
    radnet_monitor_id,
    station,
    city_name,
    state_abbr,
    monitor_latitude,
    monitor_longitude,
    radnet_monitor_label,
    old_station
  )
]

setorder(crosswalk, loc_num)

if (nrow(crosswalk) != 181) stop("最终crosswalk不是181行。")
if (uniqueN(crosswalk$loc_num) != 181) stop("loc_num不是181个唯一值。")
if (uniqueN(crosswalk$radnet_monitor_id) != 181) stop("radnet_monitor_id不是181个唯一值。")
if (crosswalk[is.na(monitor_latitude) | is.na(monitor_longitude), .N] != 0) stop("存在缺失坐标。")

fwrite(
  crosswalk,
  out_file
)

cat("\n============================================================\n")
cat("RadNet 181站命名标准化完成\n")
cat("============================================================\n")
cat("正式行数：", nrow(crosswalk), "\n")
cat("唯一 loc_num：", uniqueN(crosswalk$loc_num), "\n")
cat("唯一 radnet_monitor_id：", uniqueN(crosswalk$radnet_monitor_id), "\n")
cat("缺失纬度：", crosswalk[is.na(monitor_latitude), .N], "\n")
cat("缺失经度：", crosswalk[is.na(monitor_longitude), .N], "\n")
cat("\n正式输出：\n", out_file, "\n", sep = "")
cat("\n总体状态：PASS\n")
