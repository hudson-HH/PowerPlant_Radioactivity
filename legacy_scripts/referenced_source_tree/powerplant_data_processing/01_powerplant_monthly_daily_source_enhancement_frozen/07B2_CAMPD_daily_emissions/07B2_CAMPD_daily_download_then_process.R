#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999,
  timeout = 3600
)

# ============================================================
# 07B2_CAMPD_daily_download_then_process.R
#
# 两阶段工作流
# ------------------------------------------------------------
# PHASE 1：先把 2001–2024 共 96 个季度的 CAMPD daily
#          apportioned emissions 全部下载并原样 gzip 保存。
#
#          只有 96/96 全部下载完成后，才允许进入 PHASE 2。
#
# PHASE 2：再逐季度读取 raw gzip，完成：
#          unit-day selected
#          -> plant-day
#          -> plant-day × EPA primary fuel
#          -> QC
#
# 重要原则
# ------------------------------------------------------------
# 1. 07B2 不直接 merge 到 07B1。
# 2. CAMPD facilityId 保存为 plant_id，后续 07B3 再与
#    EIA/07A/07C 统一。
# 3. NA 不自动改成 0。
# 4. EPA primaryFuelInfo 先原样保存，不在 07B2 强制映射。
# 5. 下载与处理都支持断点续跑。
# 6. raw .json.gz 默认保留；如将来磁盘吃紧，可在全部处理
#    PASS 后再手动删除，当前代码不会自动删 raw。
#
# API
# ------------------------------------------------------------
# 使用 EPA CAMPD Streaming Services daily apportioned emissions。
# Streaming Services 与普通 emissions API 提供同类数据，
# 但不要求分页，适合大批量连续下载。
#
# 运行前：
#   export CAM_API_KEY='你的API_KEY'
#
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

required_packages <- c(
  "data.table",
  "httr",
  "jsonlite"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(httr)
  library(jsonlite)
})


# ============================================================
# 1. Paths
# ============================================================

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

CATALOG_FILE <- file.path(
  ROOT,
  "发电站数据处理",
  "07B2_CAMPD_hourly_emissions",
  "02_daily_emissions_bulk_candidates_2001_2024.csv"
)

OUT_DIR <- file.path(
  ROOT,
  "发电站数据处理",
  "07B2_CAMPD_daily_emissions"
)

RAW_DIR <- file.path(
  OUT_DIR,
  "00_raw_download_first"
)

DOWNLOAD_QC_DIR <- file.path(
  OUT_DIR,
  "01_download_qc"
)

DOWNLOAD_DONE_DIR <- file.path(
  OUT_DIR,
  "02_download_done"
)

UNIT_DIR <- file.path(
  OUT_DIR,
  "03_unit_day_selected"
)

PLANT_DIR <- file.path(
  OUT_DIR,
  "04_plant_day"
)

PLANT_FUEL_DIR <- file.path(
  OUT_DIR,
  "05_plant_day_primary_fuel"
)

PROCESS_QC_DIR <- file.path(
  OUT_DIR,
  "06_process_qc"
)

PROCESS_DONE_DIR <- file.path(
  OUT_DIR,
  "07_process_done"
)

TMP_DIR <- file.path(
  OUT_DIR,
  "99_temp"
)

for (d in c(
  OUT_DIR,
  RAW_DIR,
  DOWNLOAD_QC_DIR,
  DOWNLOAD_DONE_DIR,
  UNIT_DIR,
  PLANT_DIR,
  PLANT_FUEL_DIR,
  PROCESS_QC_DIR,
  PROCESS_DONE_DIR,
  TMP_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}


# ============================================================
# 2. Constants
# ============================================================

YEARS <- 2001:2024
QUARTERS <- 1:4

STREAM_URL <-
  "https://api.epa.gov/easey/streaming-services/emissions/apportioned/daily"

API_KEY <- Sys.getenv(
  "CAM_API_KEY"
)

if (!nzchar(API_KEY)) {
  stop(
    paste0(
      "CAM_API_KEY 未设置。\n",
      "请先在 Terminal 运行：\n",
      "export CAM_API_KEY='你的API_KEY'\n"
    ),
    call. = FALSE
  )
}

# FALSE = 直接执行完整两阶段流程：
# 先下载96季度；96/96全部完成后才进入处理。
# 如果第一个季度接口参数仍有问题，代码会立即停止，不会继续后续季度。
TEST_DOWNLOAD_ONLY_2001_Q1 <- FALSE

# 只有全部 96 个季度下载完成后才会执行处理。
PROCESS_AFTER_ALL_DOWNLOADS <- TRUE

# raw gzip 默认永久保留，保证可重复审计。
DELETE_RAW_AFTER_PROCESS <- FALSE

# 96 次 streaming 请求远低于 1000 requests/hour，
# 每个季度之间仍留一点间隔。
REQUEST_PAUSE_SECONDS <- 1


# ============================================================
# 3. Small helpers
# ============================================================

stopf <- function(...) {
  stop(
    paste0(...),
    call. = FALSE
  )
}


clean_text <- function(x) {
  y <- trimws(
    as.character(x)
  )

  y[
    y %in% c(
      "",
      "NA",
      "N/A",
      "NULL",
      "null"
    )
  ] <- NA_character_

  y
}


clean_plant_id <- function(x) {
  y <- clean_text(x)

  y <- sub(
    "\\.0+$",
    "",
    y
  )

  y
}


safe_num <- function(x) {
  suppressWarnings(
    as.numeric(x)
  )
}


sum_keep_all_na <- function(x) {
  x <- safe_num(x)

  if (
    length(x) == 0L ||
      all(is.na(x))
  ) {
    return(NA_real_)
  }

  sum(
    x,
    na.rm = TRUE
  )
}


first_nonmissing <- function(x) {
  y <- clean_text(x)
  y <- y[
    !is.na(y)
  ]

  if (length(y) == 0L) {
    return(NA_character_)
  }

  y[1L]
}


paste_unique <- function(x) {
  y <- sort(
    unique(
      clean_text(x)
    )
  )

  y <- y[
    !is.na(y)
  ]

  if (length(y) == 0L) {
    return(NA_character_)
  }

  paste(
    y,
    collapse = ";"
  )
}


file_ok <- function(path) {
  file.exists(path) &&
    is.finite(
      file.info(path)$size
    ) &&
    file.info(path)$size > 0
}


human_mb <- function(bytes) {
  round(
    bytes / 1024^2,
    2
  )
}


quarter_start_end <- function(
  year,
  quarter
) {

  start_month <- c(
    1L,
    4L,
    7L,
    10L
  )[quarter]

  start <- as.IDate(
    sprintf(
      "%04d-%02d-01",
      year,
      start_month
    )
  )

  if (quarter < 4L) {

    end <- as.IDate(
      sprintf(
        "%04d-%02d-01",
        year,
        start_month + 3L
      )
    ) - 1L

  } else {

    end <- as.IDate(
      sprintf(
        "%04d-12-31",
        year
      )
    )
  }

  list(
    start = start,
    end = end
  )
}


optional_col <- function(
  DT,
  name,
  default = NA_character_
) {

  if (name %in% names(DT)) {
    return(
      DT[[name]]
    )
  }

  rep(
    default,
    nrow(DT)
  )
}


optional_num <- function(
  DT,
  name
) {

  if (name %in% names(DT)) {
    return(
      safe_num(
        DT[[name]]
      )
    )
  }

  rep(
    NA_real_,
    nrow(DT)
  )
}


# ============================================================
# 4. Binary gzip copy helpers
# ============================================================

gzip_file_keep_source <- function(
  source,
  destination
) {

  if (file.exists(destination)) {
    unlink(
      destination,
      force = TRUE
    )
  }

  con_in <- file(
    source,
    open = "rb"
  )

  con_out <- gzfile(
    destination,
    open = "wb",
    compression = 6
  )

  on.exit({
    try(
      close(con_in),
      silent = TRUE
    )

    try(
      close(con_out),
      silent = TRUE
    )
  }, add = TRUE)

  repeat {

    buf <- readBin(
      con_in,
      what = "raw",
      n = 1024L * 1024L
    )

    if (length(buf) == 0L) {
      break
    }

    writeBin(
      buf,
      con_out
    )
  }

  close(con_in)
  close(con_out)

  if (!file_ok(destination)) {
    stopf(
      "gzip输出失败：",
      destination
    )
  }

  invisible(
    destination
  )
}


gunzip_file_keep_source <- function(
  source_gz,
  destination_json
) {

  if (file.exists(destination_json)) {
    unlink(
      destination_json,
      force = TRUE
    )
  }

  con_in <- gzfile(
    source_gz,
    open = "rb"
  )

  con_out <- file(
    destination_json,
    open = "wb"
  )

  on.exit({
    try(
      close(con_in),
      silent = TRUE
    )

    try(
      close(con_out),
      silent = TRUE
    )
  }, add = TRUE)

  repeat {

    buf <- readBin(
      con_in,
      what = "raw",
      n = 1024L * 1024L
    )

    if (length(buf) == 0L) {
      break
    }

    writeBin(
      buf,
      con_out
    )
  }

  close(con_in)
  close(con_out)

  if (!file_ok(destination_json)) {
    stopf(
      "gunzip输出失败：",
      destination_json
    )
  }

  invisible(
    destination_json
  )
}


# ============================================================
# 5. Verify the already audited daily catalog
# ============================================================

if (!file.exists(CATALOG_FILE)) {
  stopf(
    "找不到已通过PASS的daily catalog：\n",
    CATALOG_FILE
  )
}

catalog <- fread(
  CATALOG_FILE,
  showProgress = FALSE
)

if (!all(
  c(
    "filename",
    "bytes"
  ) %in% names(catalog)
)) {
  stopf(
    "daily catalog缺少 filename / bytes 字段。"
  )
}

parts <- tstrsplit(
  catalog$filename,
  "-",
  fixed = TRUE,
  fill = NA_character_
)

catalog[
  ,
  year := suppressWarnings(
    as.integer(
      parts[[3]]
    )
  )
]

catalog[
  ,
  quarter := suppressWarnings(
    as.integer(
      sub(
        ".csv",
        "",
        sub(
          "q",
          "",
          parts[[4]],
          fixed = TRUE
        ),
        fixed = TRUE
      )
    )
  )
]

catalog <- catalog[
  year %in% YEARS &
    quarter %in% QUARTERS
]

dup_catalog <- catalog[
  ,
  .N,
  by = .(
    year,
    quarter
  )
][
  N != 1L
]

if (
  nrow(catalog) != 96L ||
    uniqueN(catalog$year) != 24L ||
    nrow(dup_catalog) != 0L
) {
  stopf(
    "daily catalog结构错误：不是2001–2024 × 4季度的96个唯一文件。"
  )
}

setorder(
  catalog,
  year,
  quarter
)

catalog[
  ,
  expected_start :=
    as.character(
      vapply(
        seq_len(.N),
        function(i) {
          as.character(
            quarter_start_end(
              year[i],
              quarter[i]
            )$start
          )
        },
        character(1)
      )
    )
]

catalog[
  ,
  expected_end :=
    as.character(
      vapply(
        seq_len(.N),
        function(i) {
          as.character(
            quarter_start_end(
              year[i],
              quarter[i]
            )$end
          )
        },
        character(1)
      )
    )
]

fwrite(
  catalog,
  file.path(
    OUT_DIR,
    "00_daily_catalog_verified_2001_2024.csv"
  )
)

cat(
  "\n============================================================\n",
  "07B2 CAMPD DAILY：catalog verified\n",
  "============================================================\n",
  "文件数：",
  nrow(catalog),
  "\n年份：",
  min(catalog$year),
  "–",
  max(catalog$year),
  "\n官方 bulk CSV 总大小：",
  round(
    sum(
      catalog$bytes,
      na.rm = TRUE
    ) /
      1024^3,
    3
  ),
  " GB\n",
  sep = ""
)


# ============================================================
# 6. PHASE 1 download helper
# ============================================================

download_one_quarter_raw <- function(
  year,
  quarter,
  temp_json
) {

  dr <- quarter_start_end(
    year,
    quarter
  )

  begin_date <- as.character(
    dr$start
  )

  end_date <- as.character(
    dr$end
  )

  if (file.exists(temp_json)) {
    unlink(
      temp_json,
      force = TRUE
    )
  }

  query_base <- list(
    beginDate = begin_date,
    endDate = end_date
  )

  cat(
    "[DOWNLOAD] ",
    year,
    " Q",
    quarter,
    " -> ",
    begin_date,
    " ... ",
    end_date,
    "\n",
    sep = ""
  )

  # First try x-api-key header.
  r <- tryCatch(
    GET(
      STREAM_URL,
      query = query_base,
      add_headers(
        `x-api-key` = API_KEY,
        Accept = "application/json"
      ),
      write_disk(
        temp_json,
        overwrite = TRUE
      ),
      timeout(
        3600
      )
    ),
    error = function(e) e
  )

  if (
    inherits(
      r,
      "error"
    )
  ) {
    return(
      list(
        ok = FALSE,
        status = NA_integer_,
        auth_mode = "x-api-key_header",
        content_type = NA_character_,
        message = conditionMessage(r)
      )
    )
  }

  status1 <- status_code(
    r
  )

  auth_mode <- "x-api-key_header"

  # Authentication fallback.
  if (status1 %in% c(
    401L,
    403L
  )) {

    if (file.exists(temp_json)) {
      unlink(
        temp_json,
        force = TRUE
      )
    }

    r <- tryCatch(
      GET(
        STREAM_URL,
        query = c(
          query_base,
          list(
            api_key = API_KEY
          )
        ),
        add_headers(
          Accept = "application/json"
        ),
        write_disk(
          temp_json,
          overwrite = TRUE
        ),
        timeout(
          3600
        )
      ),
      error = function(e) e
    )

    if (
      inherits(
        r,
        "error"
      )
    ) {
      return(
        list(
          ok = FALSE,
          status = NA_integer_,
          auth_mode = "api_key_query",
          content_type = NA_character_,
          message = conditionMessage(r)
        )
      )
    }

    status1 <- status_code(
      r
    )

    auth_mode <- "api_key_query"
  }

  ct <- headers(r)[["content-type"]]

  if (is.null(ct)) {
    ct <- NA_character_
  }

  if (
    status1 < 200L ||
      status1 >= 300L ||
      !file_ok(temp_json)
  ) {

    response_preview <- ""

    if (file.exists(temp_json)) {
      response_preview <- tryCatch(
        paste(
          readLines(
            temp_json,
            n = 40,
            warn = FALSE,
            encoding = "UTF-8"
          ),
          collapse = "\n"
        ),
        error = function(e) ""
      )
    }

    return(
      list(
        ok = FALSE,
        status = status1,
        auth_mode = auth_mode,
        content_type = as.character(ct),
        message = substr(
          response_preview,
          1,
          4000
        )
      )
    )
  }

  list(
    ok = TRUE,
    status = status1,
    auth_mode = auth_mode,
    content_type = as.character(ct),
    message = "",
    begin_date = begin_date,
    end_date = end_date
  )
}


# ============================================================
# 7. PHASE 1: DOWNLOAD ALL FIRST
# ============================================================

cat(
  "\n============================================================\n",
  "PHASE 1：先下载，再处理\n",
  "============================================================\n",
  sep = ""
)

download_catalog <- if (
  TEST_DOWNLOAD_ONLY_2001_Q1
) {
  catalog[
    year == 2001L &
      quarter == 1L
  ]
} else {
  copy(
    catalog
  )
}

if (TEST_DOWNLOAD_ONLY_2001_Q1) {
  cat(
    "当前为下载接口安全测试：只下载 2001 Q1。\n",
    "2001 Q1 PASS 后，将 TEST_DOWNLOAD_ONLY_2001_Q1 <- FALSE，",
    "重新运行即可下载完整96季度。\n\n",
    sep = ""
  )
}

for (i in seq_len(
  nrow(download_catalog)
)) {

  yy <- download_catalog$year[i]
  qq <- download_catalog$quarter[i]

  tag <- sprintf(
    "%04d_q%d",
    yy,
    qq
  )

  raw_gz <- file.path(
    RAW_DIR,
    paste0(
      "CAMPD_daily_raw_",
      tag,
      ".json.gz"
    )
  )

  temp_json <- file.path(
    TMP_DIR,
    paste0(
      "TEMP_download_",
      tag,
      ".json"
    )
  )

  qc_file <- file.path(
    DOWNLOAD_QC_DIR,
    paste0(
      "download_qc_",
      tag,
      ".csv"
    )
  )

  done_file <- file.path(
    DOWNLOAD_DONE_DIR,
    paste0(
      tag,
      ".download.done"
    )
  )

  cat(
    "\n------------------------------------------------------------\n",
    "[DOWNLOAD ",
    i,
    "/",
    nrow(download_catalog),
    "] ",
    yy,
    " Q",
    qq,
    "\n",
    sep = ""
  )

  if (
    file.exists(done_file) &&
      file_ok(raw_gz) &&
      file_ok(qc_file)
  ) {
    cat(
      "已有download PASS checkpoint，跳过。\n"
    )

    next
  }

  if (file.exists(done_file)) {
    unlink(
      done_file,
      force = TRUE
    )
  }

  dl <- download_one_quarter_raw(
    year = yy,
    quarter = qq,
    temp_json = temp_json
  )

  if (!dl$ok) {

    stopf(
      "CAMPD streaming下载失败：",
      yy,
      " Q",
      qq,
      "\nHTTP status=",
      dl$status,
      "\nauth=",
      dl$auth_mode,
      "\ncontent-type=",
      dl$content_type,
      "\nEPA response / error:\n",
      dl$message,
      "\n\n程序停止；已有成功季度不会丢失。"
    )
  }

  raw_uncompressed_bytes <- file.info(
    temp_json
  )$size

  # Minimal format gate only.
  # This is not scientific processing; it only prevents saving an HTML/error page.
  first_char <- tryCatch(
    {
      con <- file(
        temp_json,
        open = "rb"
      )

      on.exit(
        close(con),
        add = TRUE
      )

      b <- readBin(
        con,
        "raw",
        n = 256L
      )

      txt <- rawToChar(
        b
      )

      sub(
        "^\\s*",
        "",
        txt
      )
    },
    error = function(e) ""
  )

  if (
    !startsWith(first_char, "[") &&
      !startsWith(first_char, "{")
  ) {
    stopf(
      "下载得到的内容不像JSON：",
      yy,
      " Q",
      qq,
      "\n临时文件保留：",
      temp_json
    )
  }

  gzip_file_keep_source(
    source = temp_json,
    destination = raw_gz
  )

  raw_compressed_bytes <- file.info(
    raw_gz
  )$size

  # Verify gzip can actually be opened.
  gz_test <- tryCatch(
    {
      cc <- gzfile(
        raw_gz,
        open = "rb"
      )

      on.exit(
        close(cc),
        add = TRUE
      )

      length(
        readBin(
          cc,
          "raw",
          n = 32L
        )
      ) > 0L
    },
    error = function(e) FALSE
  )

  if (!gz_test) {
    stopf(
      "raw gzip可读性检查失败：",
      raw_gz
    )
  }

  download_qc <- data.table(
    year = yy,
    quarter = qq,
    begin_date =
      dl$begin_date,
    end_date =
      dl$end_date,
    http_status =
      dl$status,
    auth_mode =
      dl$auth_mode,
    content_type =
      dl$content_type,
    uncompressed_json_bytes =
      raw_uncompressed_bytes,
    compressed_json_gz_bytes =
      raw_compressed_bytes,
    compression_ratio =
      raw_compressed_bytes /
      raw_uncompressed_bytes,
    raw_gz_path =
      raw_gz,
    status =
      "PASS"
  )

  fwrite(
    download_qc,
    qc_file
  )

  # Only after gzip + QC PASS do we remove uncompressed temp.
  unlink(
    temp_json,
    force = TRUE
  )

  if (file.exists(temp_json)) {
    stopf(
      "临时未压缩JSON无法删除：",
      temp_json
    )
  }

  writeLines(
    c(
      paste0(
        "year=",
        yy
      ),
      paste0(
        "quarter=",
        qq
      ),
      "download_status=PASS",
      paste0(
        "raw_gz=",
        raw_gz
      ),
      paste0(
        "completed_at=",
        format(
          Sys.time(),
          "%Y-%m-%d %H:%M:%S"
        )
      )
    ),
    done_file,
    useBytes = TRUE
  )

  cat(
    "DOWNLOAD PASS：raw gzip = ",
    human_mb(
      raw_compressed_bytes
    ),
    " MB\n",
    sep = ""
  )

  if (
    REQUEST_PAUSE_SECONDS > 0
  ) {
    Sys.sleep(
      REQUEST_PAUSE_SECONDS
    )
  }
}


# ============================================================
# 8. Download completion gate
# ============================================================

expected <- CJ(
  year = YEARS,
  quarter = QUARTERS
)

expected[
  ,
  tag := sprintf(
    "%04d_q%d",
    year,
    quarter
  )
]

expected[
  ,
  raw_gz := file.path(
    RAW_DIR,
    paste0(
      "CAMPD_daily_raw_",
      tag,
      ".json.gz"
    )
  )
]

expected[
  ,
  done_file := file.path(
    DOWNLOAD_DONE_DIR,
    paste0(
      tag,
      ".download.done"
    )
  )
]

expected[
  ,
  download_complete :=
    file.exists(
      raw_gz
    ) &
    file.exists(
      done_file
    )
]

fwrite(
  expected,
  file.path(
    OUT_DIR,
    "08_download_completion_status.csv"
  )
)

n_downloaded <- expected[
  download_complete == TRUE,
  .N
]

cat(
  "\n============================================================\n",
  "PHASE 1 下载状态\n",
  "============================================================\n",
  "完成：",
  n_downloaded,
  " / 96\n",
  sep = ""
)

# In test mode, intentionally stop here.
if (TEST_DOWNLOAD_ONLY_2001_Q1) {

  cat(
    "\n2001 Q1下载测试完成后，请先检查终端是否显示 DOWNLOAD PASS。\n",
    "确认后将代码顶部：\n",
    "TEST_DOWNLOAD_ONLY_2001_Q1 <- TRUE\n",
    "改为：\n",
    "TEST_DOWNLOAD_ONLY_2001_Q1 <- FALSE\n",
    "再运行同一个脚本即可下载全部96季度。\n",
    "\n当前不会进入数据处理阶段。\n",
    sep = ""
  )

  quit(
    save = "no",
    status = 0
  )
}

if (
  n_downloaded != 96L
) {
  stopf(
    "只有 ",
    n_downloaded,
    "/96 个季度下载完成。",
    "\n按“先全部下载、后处理”原则，PHASE 2不会启动。",
    "\n重新运行同一脚本即可断点续下。"
  )
}

cat(
  "\n96/96 DOWNLOAD PASS。\n",
  "现在才进入 PHASE 2 数据处理。\n",
  sep = ""
)


# ============================================================
# 9. Normalize daily raw data
# ============================================================

normalize_daily <- function(
  raw_obj,
  year,
  quarter
) {

  if (
    is.list(raw_obj) &&
      !is.data.frame(raw_obj) &&
      "items" %in% names(raw_obj)
  ) {
    raw_obj <- raw_obj$items
  }

  DT <- as.data.table(
    raw_obj
  )

  if (nrow(DT) == 0L) {
    stopf(
      year,
      " Q",
      quarter,
      " raw数据0行。"
    )
  }

  required <- c(
    "stateCode",
    "facilityName",
    "facilityId",
    "unitId",
    "date",
    "primaryFuelInfo",
    "heatInput",
    "so2Mass",
    "noxMass",
    "co2Mass",
    "grossLoad"
  )

  missing_required <- setdiff(
    required,
    names(DT)
  )

  if (length(missing_required) > 0L) {
    stopf(
      year,
      " Q",
      quarter,
      " 缺少必要字段：",
      paste(
        missing_required,
        collapse = ", "
      ),
      "\n实际字段：",
      paste(
        names(DT),
        collapse = ", "
      )
    )
  }

  out <- data.table(
    state_code =
      clean_text(
        DT$stateCode
      ),

    facility_name =
      clean_text(
        DT$facilityName
      ),

    plant_id =
      clean_plant_id(
        DT$facilityId
      ),

    unit_id =
      clean_text(
        DT$unitId
      ),

    database_unit_id =
      clean_plant_id(
        optional_col(
          DT,
          "unit_id"
        )
      ),

    unit_type =
      clean_text(
        optional_col(
          DT,
          "unitType"
        )
      ),

    associated_stacks =
      clean_text(
        optional_col(
          DT,
          "associatedStacks"
        )
      ),

    primary_fuel_info =
      clean_text(
        DT$primaryFuelInfo
      ),

    secondary_fuel_info =
      clean_text(
        optional_col(
          DT,
          "secondaryFuelInfo"
        )
      ),

    so2_control_info =
      clean_text(
        optional_col(
          DT,
          "so2ControlInfo"
        )
      ),

    pm_control_info =
      clean_text(
        optional_col(
          DT,
          "pmControlInfo"
        )
      ),

    nox_control_info =
      clean_text(
        optional_col(
          DT,
          "noxControlInfo"
        )
      ),

    hg_control_info =
      clean_text(
        optional_col(
          DT,
          "hgControlInfo"
        )
      ),

    program_code_info =
      clean_text(
        optional_col(
          DT,
          "programCodeInfo"
        )
      ),

    date =
      as.IDate(
        substr(
          clean_text(
            DT$date
          ),
          1L,
          10L
        )
      ),

    gross_load =
      optional_num(
        DT,
        "grossLoad"
      ),

    steam_load =
      optional_num(
        DT,
        "steamLoad"
      ),

    so2_mass =
      optional_num(
        DT,
        "so2Mass"
      ),

    so2_rate =
      optional_num(
        DT,
        "so2Rate"
      ),

    nox_mass =
      optional_num(
        DT,
        "noxMass"
      ),

    nox_rate =
      optional_num(
        DT,
        "noxRate"
      ),

    co2_mass =
      optional_num(
        DT,
        "co2Mass"
      ),

    co2_rate =
      optional_num(
        DT,
        "co2Rate"
      ),

    heat_input_mmbtu =
      optional_num(
        DT,
        "heatInput"
      ),

    heat_input_measure_flag =
      clean_text(
        optional_col(
          DT,
          "heatInputMeasureFlg"
        )
      ),

    sum_op_time =
      optional_num(
        DT,
        "sumOpTime"
      ),

    count_op_time =
      optional_num(
        DT,
        "countOpTime"
      )
  )

  out[
    ,
    `:=`(
      year =
        as.integer(
          format(
            date,
            "%Y"
          )
        ),

      month =
        as.integer(
          format(
            date,
            "%m"
          )
        )
    )
  ]

  out[
    ,
    quarter :=
      as.integer(
        (month - 1L) %/% 3L +
          1L
      )
  ]

  setcolorder(
    out,
    c(
      "year",
      "quarter",
      "month",
      "date",
      setdiff(
        names(out),
        c(
          "year",
          "quarter",
          "month",
          "date"
        )
      )
    )
  )

  out
}


# ============================================================
# 10. Aggregation
# ============================================================

aggregate_plant_day <- function(DT) {

  DT[
    ,
    .(
      state_code =
        first_nonmissing(
          state_code
        ),

      facility_name =
        first_nonmissing(
          facility_name
        ),

      primary_fuels =
        paste_unique(
          primary_fuel_info
        ),

      secondary_fuels =
        paste_unique(
          secondary_fuel_info
        ),

      unit_types =
        paste_unique(
          unit_type
        ),

      n_units =
        uniqueN(
          unit_id[
            !is.na(unit_id)
          ]
        ),

      n_rows =
        .N,

      n_heat_input_nonmissing =
        sum(
          is.finite(
            heat_input_mmbtu
          )
        ),

      n_so2_nonmissing =
        sum(
          is.finite(
            so2_mass
          )
        ),

      n_nox_nonmissing =
        sum(
          is.finite(
            nox_mass
          )
        ),

      n_co2_nonmissing =
        sum(
          is.finite(
            co2_mass
          )
        ),

      n_gross_load_nonmissing =
        sum(
          is.finite(
            gross_load
          )
        ),

      daily_heat_input_mmbtu_observed =
        sum_keep_all_na(
          heat_input_mmbtu
        ),

      daily_so2_mass_observed =
        sum_keep_all_na(
          so2_mass
        ),

      daily_nox_mass_observed =
        sum_keep_all_na(
          nox_mass
        ),

      daily_co2_mass_observed =
        sum_keep_all_na(
          co2_mass
        ),

      daily_gross_load_observed =
        sum_keep_all_na(
          gross_load
        ),

      daily_steam_load_observed =
        sum_keep_all_na(
          steam_load
        ),

      daily_sum_op_time =
        sum_keep_all_na(
          sum_op_time
        ),

      daily_count_op_time =
        sum_keep_all_na(
          count_op_time
        )
    ),
    by = .(
      year,
      quarter,
      month,
      date,
      plant_id
    )
  ]
}


aggregate_plant_day_primary_fuel <- function(DT) {

  DT[
    ,
    .(
      state_code =
        first_nonmissing(
          state_code
        ),

      facility_name =
        first_nonmissing(
          facility_name
        ),

      n_units =
        uniqueN(
          unit_id[
            !is.na(unit_id)
          ]
        ),

      n_rows =
        .N,

      n_heat_input_nonmissing =
        sum(
          is.finite(
            heat_input_mmbtu
          )
        ),

      n_so2_nonmissing =
        sum(
          is.finite(
            so2_mass
          )
        ),

      n_nox_nonmissing =
        sum(
          is.finite(
            nox_mass
          )
        ),

      n_co2_nonmissing =
        sum(
          is.finite(
            co2_mass
          )
        ),

      daily_heat_input_mmbtu_observed =
        sum_keep_all_na(
          heat_input_mmbtu
        ),

      daily_so2_mass_observed =
        sum_keep_all_na(
          so2_mass
        ),

      daily_nox_mass_observed =
        sum_keep_all_na(
          nox_mass
        ),

      daily_co2_mass_observed =
        sum_keep_all_na(
          co2_mass
        ),

      daily_gross_load_observed =
        sum_keep_all_na(
          gross_load
        ),

      daily_sum_op_time =
        sum_keep_all_na(
          sum_op_time
        ),

      daily_count_op_time =
        sum_keep_all_na(
          count_op_time
        )
    ),
    by = .(
      year,
      quarter,
      month,
      date,
      plant_id,
      primary_fuel_info
    )
  ]
}


# ============================================================
# 11. PHASE 2: PROCESS ONLY AFTER ALL 96 DOWNLOADS PASS
# ============================================================

if (!PROCESS_AFTER_ALL_DOWNLOADS) {

  cat(
    "\nPROCESS_AFTER_ALL_DOWNLOADS = FALSE。\n",
    "96季度raw已经下载完成，但本次不处理。\n",
    sep = ""
  )

  quit(
    save = "no",
    status = 0
  )
}

cat(
  "\n============================================================\n",
  "PHASE 2：96季度全部下载完成，现在开始处理\n",
  "============================================================\n",
  sep = ""
)

for (i in seq_len(
  nrow(catalog)
)) {

  yy <- catalog$year[i]
  qq <- catalog$quarter[i]

  tag <- sprintf(
    "%04d_q%d",
    yy,
    qq
  )

  raw_gz <- file.path(
    RAW_DIR,
    paste0(
      "CAMPD_daily_raw_",
      tag,
      ".json.gz"
    )
  )

  temp_json <- file.path(
    TMP_DIR,
    paste0(
      "TEMP_process_",
      tag,
      ".json"
    )
  )

  unit_out <- file.path(
    UNIT_DIR,
    paste0(
      "CAMPD_unit_day_selected_",
      tag,
      ".csv.gz"
    )
  )

  plant_out <- file.path(
    PLANT_DIR,
    paste0(
      "CAMPD_plant_day_",
      tag,
      ".csv.gz"
    )
  )

  plant_fuel_out <- file.path(
    PLANT_FUEL_DIR,
    paste0(
      "CAMPD_plant_day_primary_fuel_",
      tag,
      ".csv.gz"
    )
  )

  qc_out <- file.path(
    PROCESS_QC_DIR,
    paste0(
      "process_qc_",
      tag,
      ".csv"
    )
  )

  done_file <- file.path(
    PROCESS_DONE_DIR,
    paste0(
      tag,
      ".process.done"
    )
  )

  cat(
    "\n------------------------------------------------------------\n",
    "[PROCESS ",
    i,
    "/96] ",
    yy,
    " Q",
    qq,
    "\n",
    sep = ""
  )

  if (
    file.exists(done_file) &&
      file_ok(unit_out) &&
      file_ok(plant_out) &&
      file_ok(plant_fuel_out) &&
      file_ok(qc_out)
  ) {
    cat(
      "已有process PASS checkpoint，跳过。\n"
    )

    next
  }

  if (!file_ok(raw_gz)) {
    stopf(
      "缺少raw gzip：",
      raw_gz
    )
  }

  if (file.exists(done_file)) {
    unlink(
      done_file,
      force = TRUE
    )
  }

  # ----------------------------
  # 11A. Decompress one quarter
  # ----------------------------

  gunzip_file_keep_source(
    source_gz = raw_gz,
    destination_json = temp_json
  )

  # ----------------------------
  # 11B. Parse
  # ----------------------------

  raw_obj <- tryCatch(
    fromJSON(
      temp_json,
      simplifyDataFrame = TRUE
    ),
    error = function(e) e
  )

  if (
    inherits(
      raw_obj,
      "error"
    )
  ) {
    stopf(
      "JSON解析失败：",
      yy,
      " Q",
      qq,
      "\n",
      conditionMessage(raw_obj),
      "\n临时解压文件保留：",
      temp_json
    )
  }

  D <- normalize_daily(
    raw_obj = raw_obj,
    year = yy,
    quarter = qq
  )

  rm(
    raw_obj
  )

  invisible(
    gc()
  )

  # 解析成功后即可删临时解压版本；
  # 原始 gzip 仍永久保留。
  unlink(
    temp_json,
    force = TRUE
  )

  # ----------------------------
  # 11C. Hard QC
  # ----------------------------

  expected_dates <- quarter_start_end(
    yy,
    qq
  )

  n_missing_plant_id <-
    D[
      is.na(plant_id) |
        plant_id == "",
      .N
    ]

  n_missing_date <-
    D[
      is.na(date),
      .N
    ]

  n_wrong_year <-
    D[
      !is.na(year) &
        year != yy,
      .N
    ]

  n_wrong_quarter <-
    D[
      !is.na(quarter) &
        quarter != qq,
      .N
    ]

  n_outside_date_range <-
    D[
      !is.na(date) &
        (
          date < expected_dates$start |
            date > expected_dates$end
        ),
      .N
    ]

  dup_unit_day <- D[
    ,
    .N,
    by = .(
      plant_id,
      unit_id,
      date
    )
  ][
    N > 1L
  ]

  hard_fail <-
    n_missing_plant_id +
    n_missing_date +
    n_wrong_year +
    n_wrong_quarter +
    n_outside_date_range +
    nrow(
      dup_unit_day
    )

  if (hard_fail > 0L) {

    if (nrow(
      dup_unit_day
    ) > 0L) {
      fwrite(
        dup_unit_day,
        file.path(
          PROCESS_QC_DIR,
          paste0(
            "DUPLICATE_unit_day_",
            tag,
            ".csv"
          )
        )
      )
    }

    stopf(
      "硬性QC失败：",
      yy,
      " Q",
      qq,
      "\nmissing plant_id=",
      n_missing_plant_id,
      "\nmissing date=",
      n_missing_date,
      "\nwrong year=",
      n_wrong_year,
      "\nwrong quarter=",
      n_wrong_quarter,
      "\noutside date range=",
      n_outside_date_range,
      "\nduplicate unit-day keys=",
      nrow(
        dup_unit_day
      )
    )
  }

  # ----------------------------
  # 11D. Aggregate
  # ----------------------------

  P <- aggregate_plant_day(
    D
  )

  PF <- aggregate_plant_day_primary_fuel(
    D
  )

  # ----------------------------
  # 11E. Save compressed outputs
  # ----------------------------

  fwrite(
    D,
    unit_out,
    compress = "gzip"
  )

  fwrite(
    P,
    plant_out,
    compress = "gzip"
  )

  fwrite(
    PF,
    plant_fuel_out,
    compress = "gzip"
  )

  # ----------------------------
  # 11F. Process QC
  # ----------------------------

  qc <- data.table(
    year = yy,
    quarter = qq,

    n_unit_day_rows =
      nrow(D),

    n_unique_plants =
      uniqueN(
        D$plant_id
      ),

    n_unique_units =
      uniqueN(
        paste(
          D$plant_id,
          D$unit_id,
          sep = "::"
        )
      ),

    n_unique_dates =
      uniqueN(
        D$date
      ),

    date_min =
      as.character(
        min(
          D$date,
          na.rm = TRUE
        )
      ),

    date_max =
      as.character(
        max(
          D$date,
          na.rm = TRUE
        )
      ),

    n_plant_day_rows =
      nrow(P),

    n_plant_day_primary_fuel_rows =
      nrow(PF),

    n_primary_fuel_categories =
      uniqueN(
        D$primary_fuel_info,
        na.rm = TRUE
      ),

    n_heat_input_missing =
      sum(
        !is.finite(
          D$heat_input_mmbtu
        )
      ),

    pct_heat_input_missing =
      mean(
        !is.finite(
          D$heat_input_mmbtu
        )
      ),

    n_so2_missing =
      sum(
        !is.finite(
          D$so2_mass
        )
      ),

    pct_so2_missing =
      mean(
        !is.finite(
          D$so2_mass
        )
      ),

    n_nox_missing =
      sum(
        !is.finite(
          D$nox_mass
        )
      ),

    pct_nox_missing =
      mean(
        !is.finite(
          D$nox_mass
        )
      ),

    n_co2_missing =
      sum(
        !is.finite(
          D$co2_mass
        )
      ),

    pct_co2_missing =
      mean(
        !is.finite(
          D$co2_mass
        )
      ),

    n_gross_load_missing =
      sum(
        !is.finite(
          D$gross_load
        )
      ),

    pct_gross_load_missing =
      mean(
        !is.finite(
          D$gross_load
        )
      ),

    n_heat_input_negative =
      sum(
        D$heat_input_mmbtu < 0,
        na.rm = TRUE
      ),

    n_so2_negative =
      sum(
        D$so2_mass < 0,
        na.rm = TRUE
      ),

    n_nox_negative =
      sum(
        D$nox_mass < 0,
        na.rm = TRUE
      ),

    n_co2_negative =
      sum(
        D$co2_mass < 0,
        na.rm = TRUE
      ),

    n_duplicate_unit_day_keys =
      nrow(
        dup_unit_day
      ),

    unit_day_output_bytes =
      file.info(
        unit_out
      )$size,

    plant_day_output_bytes =
      file.info(
        plant_out
      )$size,

    plant_day_primary_fuel_output_bytes =
      file.info(
        plant_fuel_out
      )$size,

    hard_status =
      "PASS"
  )

  fwrite(
    qc,
    qc_out
  )

  # ----------------------------
  # 11G. Verify outputs
  # ----------------------------

  if (
    !file_ok(unit_out) ||
      !file_ok(plant_out) ||
      !file_ok(plant_fuel_out) ||
      !file_ok(qc_out)
  ) {
    stopf(
      "输出文件验证失败：",
      yy,
      " Q",
      qq
    )
  }

  test_unit <- tryCatch(
    fread(
      unit_out,
      nrows = 2L,
      showProgress = FALSE
    ),
    error = function(e) e
  )

  test_plant <- tryCatch(
    fread(
      plant_out,
      nrows = 2L,
      showProgress = FALSE
    ),
    error = function(e) e
  )

  test_pf <- tryCatch(
    fread(
      plant_fuel_out,
      nrows = 2L,
      showProgress = FALSE
    ),
    error = function(e) e
  )

  if (
    inherits(
      test_unit,
      "error"
    ) ||
      inherits(
        test_plant,
        "error"
      ) ||
      inherits(
        test_pf,
        "error"
      )
  ) {
    stopf(
      "压缩输出可读性验证失败：",
      yy,
      " Q",
      qq
    )
  }

  writeLines(
    c(
      paste0(
        "year=",
        yy
      ),
      paste0(
        "quarter=",
        qq
      ),
      "process_status=PASS",
      paste0(
        "completed_at=",
        format(
          Sys.time(),
          "%Y-%m-%d %H:%M:%S"
        )
      )
    ),
    done_file,
    useBytes = TRUE
  )

  if (DELETE_RAW_AFTER_PROCESS) {
    unlink(
      raw_gz,
      force = TRUE
    )
  }

  cat(
    "PROCESS PASS\n"
  )

  rm(
    D,
    P,
    PF,
    qc,
    test_unit,
    test_plant,
    test_pf
  )

  invisible(
    gc()
  )
}


# ============================================================
# 12. Consolidate download QC
# ============================================================

download_qc_files <- list.files(
  DOWNLOAD_QC_DIR,
  pattern = "^download_qc_[0-9]{4}_q[1-4]\\.csv$",
  full.names = TRUE
)

download_qc_all <- rbindlist(
  lapply(
    download_qc_files,
    function(f) {
      fread(
        f,
        showProgress = FALSE
      )
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

setorder(
  download_qc_all,
  year,
  quarter
)

fwrite(
  download_qc_all,
  file.path(
    OUT_DIR,
    "09_download_qc_all_quarters_2001_2024.csv"
  )
)


# ============================================================
# 13. Consolidate processing QC
# ============================================================

process_qc_files <- list.files(
  PROCESS_QC_DIR,
  pattern = "^process_qc_[0-9]{4}_q[1-4]\\.csv$",
  full.names = TRUE
)

process_qc_all <- rbindlist(
  lapply(
    process_qc_files,
    function(f) {
      fread(
        f,
        showProgress = FALSE
      )
    }
  ),
  use.names = TRUE,
  fill = TRUE
)

setorder(
  process_qc_all,
  year,
  quarter
)

fwrite(
  process_qc_all,
  file.path(
    OUT_DIR,
    "10_process_qc_all_quarters_2001_2024.csv"
  )
)


# ============================================================
# 14. Primary fuel audit
# ============================================================

plant_fuel_files <- list.files(
  PLANT_FUEL_DIR,
  pattern = "^CAMPD_plant_day_primary_fuel_[0-9]{4}_q[1-4]\\.csv\\.gz$",
  full.names = TRUE
)

fuel_dict_list <- lapply(
  plant_fuel_files,
  function(f) {

    x <- fread(
      f,
      select = c(
        "year",
        "primary_fuel_info",
        "plant_id",
        "daily_heat_input_mmbtu_observed"
      ),
      showProgress = FALSE
    )

    x[
      ,
      .(
        n_plant_days =
          .N,

        n_unique_plants =
          uniqueN(
            plant_id
          ),

        total_heat_input_mmbtu_observed =
          sum(
            daily_heat_input_mmbtu_observed,
            na.rm = TRUE
          )
      ),
      by = .(
        year,
        primary_fuel_info
      )
    ]
  }
)

fuel_dict <- rbindlist(
  fuel_dict_list,
  use.names = TRUE,
  fill = TRUE
)

setorder(
  fuel_dict,
  year,
  primary_fuel_info
)

fwrite(
  fuel_dict,
  file.path(
    OUT_DIR,
    "11_primary_fuel_info_audit_by_year.csv"
  )
)


# ============================================================
# 15. Final completion
# ============================================================

completion <- CJ(
  year = YEARS,
  quarter = QUARTERS
)

completion[
  ,
  tag := sprintf(
    "%04d_q%d",
    year,
    quarter
  )
]

completion[
  ,
  download_done :=
    file.exists(
      file.path(
        DOWNLOAD_DONE_DIR,
        paste0(
          tag,
          ".download.done"
        )
      )
    )
]

completion[
  ,
  process_done :=
    file.exists(
      file.path(
        PROCESS_DONE_DIR,
        paste0(
          tag,
          ".process.done"
        )
      )
    )
]

fwrite(
  completion,
  file.path(
    OUT_DIR,
    "12_final_quarter_completion_status.csv"
  )
)

n_download_done <- completion[
  download_done == TRUE,
  .N
]

n_process_done <- completion[
  process_done == TRUE,
  .N
]

hard_final_status <- if (
  n_download_done == 96L &&
    n_process_done == 96L
) {
  "PASS"
} else {
  "INCOMPLETE"
}


# ============================================================
# 16. Report
# ============================================================

report_file <- file.path(
  OUT_DIR,
  "13_07B2_CAMPD_daily_download_then_process_report.md"
)

report_lines <- c(
  "# 07B2 CAMPD daily emissions — download first, process second",
  "",
  paste0(
    "- Generated: ",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  "- Study period: 2001–2024",
  "- Expected quarters: 96",
  paste0(
    "- Download complete: ",
    n_download_done,
    "/96"
  ),
  paste0(
    "- Processing complete: ",
    n_process_done,
    "/96"
  ),
  paste0(
    "- Final status: ",
    hard_final_status
  ),
  "",
  "## Architecture",
  "",
  "1. All quarter-level raw CAMPD daily responses are downloaded first.",
  "2. Raw responses are gzip-compressed without scientific transformation.",
  "3. Processing starts only when all 96 quarters have download checkpoints.",
  "4. Unit-day selected, plant-day and plant-day × primary fuel outputs are then built.",
  "",
  "## Scientific handling",
  "",
  "- facilityId is preserved as plant_id.",
  "- NA is never automatically converted to zero.",
  "- EPA primaryFuelInfo is preserved in raw form in 07B2.",
  "- Plant-day sums carry nonmissing-unit counts for key metrics.",
  "- 07B2 is not directly merged into 07B1; harmonization is deferred to 07B3.",
  "",
  "## Raw-data retention",
  "",
  paste0(
    "- DELETE_RAW_AFTER_PROCESS = ",
    DELETE_RAW_AFTER_PROCESS
  ),
  "- Default is FALSE, so raw .json.gz files are retained for reproducibility."
)

writeLines(
  report_lines,
  report_file,
  useBytes = TRUE
)

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    OUT_DIR,
    "14_sessionInfo.txt"
  ),
  useBytes = TRUE
)


# ============================================================
# 17. Summary
# ============================================================

cat(
  "\n============================================================\n",
  "07B2 CAMPD DAILY：完成\n",
  "============================================================\n",
  "Download：",
  n_download_done,
  "/96\n",
  "Process：",
  n_process_done,
  "/96\n",
  "Final status：",
  hard_final_status,
  "\n",
  sep = ""
)

cat(
  "\n核心输出目录：\n",
  OUT_DIR,
  "\n\n",
  "优先检查：\n",
  file.path(
    OUT_DIR,
    "09_download_qc_all_quarters_2001_2024.csv"
  ),
  "\n",
  file.path(
    OUT_DIR,
    "10_process_qc_all_quarters_2001_2024.csv"
  ),
  "\n",
  file.path(
    OUT_DIR,
    "11_primary_fuel_info_audit_by_year.csv"
  ),
  "\n",
  file.path(
    OUT_DIR,
    "12_final_quarter_completion_status.csv"
  ),
  "\n",
  report_file,
  "\n",
  sep = ""
)
