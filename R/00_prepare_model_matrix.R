#!/usr/bin/env Rscript

# Optional curator route: rebuild the portable analysis matrix from the local,
# read-only frozen 07F/07H products. This is not invoked by run_all.R because
# the prepared matrix is distributed with this code package.

suppressPackageStartupMessages(library(data.table))
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value=TRUE)[1]))), "common.R"))

research_root <- Sys.getenv("RADNET_RESEARCH_ROOT", unset = "")
if (!nzchar(research_root)) stop("Set RADNET_RESEARCH_ROOT to the 'USA 辐射表征研究' directory.", call. = FALSE)
research_root <- normalizePath(research_root, mustWork = TRUE)

pp <- file.path(research_root, "发电站数据处理")
hroot <- file.path(pp, "06_07H_燃料特异CAMPD逐日暴露与模型")
fp_file <- file.path(pp, "04_07F_正式统计分析数据构建", "03_07F_filter_period_analysis_base_2001_2024.csv.gz")
qdir <- file.path(pp, "04_07F_正式统计分析数据构建", "quasi_daily_yearly")
daily_dir <- file.path(hroot, "fuel_daily_yearly")

years <- 2001:2024
radii <- seq(20L, 50L, 5L)
fuels <- c("coal", "gas", "oil", "biomass")
directions <- c("upwind", "downwind")

clean_chr <- function(x) { y <- trimws(as.character(x)); y[y %in% c("", "NA", "N/A", "NULL", ".")] <- NA_character_; y }
as_bool <- function(x) tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
safe_num <- function(x) suppressWarnings(as.numeric(x))
strict_mean <- function(x) if (all(is.finite(x))) mean(x) else NA_real_

qfiles <- file.path(qdir, paste0("07F_quasi_daily_analysis_base_", years, ".csv.gz"))
assert_true(all(file.exists(qfiles)), "One or more formal 07F quasi-daily files is missing.")
eligible <- unlist(lapply(qfiles, function(f) {
  z <- fread(f, select = c("loc_num", "flag_07f_primary_quasi_sample"), showProgress = FALSE)
  z[, `:=`(loc_num = clean_chr(loc_num), flag_07f_primary_quasi_sample = as_bool(flag_07f_primary_quasi_sample))]
  z[flag_07f_primary_quasi_sample == TRUE, unique(loc_num)]
}), use.names = FALSE)
eligible_locs <- sort(unique(eligible))
assert_true(length(eligible_locs) == 133L, "Formal quasi-daily universe must contain 133 monitors.")

base_required <- c(
  "loc_num", "radnet_monitor_id", "collect_start", "collect_end",
  "gross_beta_bq_m3", "gross_beta_mBq_m3", "c_prop", "surface_u_ppm_article",
  "sunspot_monthly", "pm_filter_period_official", "inv_hpbl", "hpbl", "rhum",
  "soilm", "vel", "lat_diff_35", "year_model", "temp", "flag_conus",
  "fpmean_campd_upwind_heat_observed_20km_90deg"
)

fp <- fread(fp_file, select = base_required, showProgress = TRUE)
fp[, `:=`(loc_num = clean_chr(loc_num), collect_start = as.IDate(collect_start),
          collect_end = as.IDate(collect_end), flag_conus = as_bool(flag_conus))]
numeric_cols <- setdiff(base_required, c("loc_num", "radnet_monitor_id", "collect_start", "collect_end", "flag_conus"))
for (cc in numeric_cols) set(fp, j = cc, value = safe_num(fp[[cc]]))

primary_fp <- "fpmean_campd_upwind_heat_observed_20km_90deg"
fp <- fp[flag_conus & loc_num %in% eligible_locs & is.finite(gross_beta_bq_m3) & gross_beta_bq_m3 < 0.04 &
           is.finite(pm_filter_period_official) & pm_filter_period_official < 35 &
           is.finite(vel) & vel < 6 & is.finite(hpbl) & hpbl > 150 & hpbl < 3000 &
           is.finite(c_prop) & is.finite(soilm) & is.finite(get(primary_fp))]
assert_true(nrow(fp) == 107407L, paste0("Expected 107407 formal periods; observed ", nrow(fp)))
assert_true(uniqueN(fp$loc_num) == 133L, "Expected 133 monitors in the formal sample.")

fp[, `:=`(period_id = .I, period_days = as.integer(collect_end - collect_start),
          pm_model_centered = pm_filter_period_official - 10,
          temp_centered = temp - 273.15, temp2 = (temp - 273.15)^2,
          temp3 = (temp - 273.15)^3,
          calendar_year = as.integer(format(collect_start, "%Y")))]
fp[, monitor_year_fe := interaction(loc_num, calendar_year, drop = TRUE, lex.order = TRUE)]

exposure_cols <- unlist(lapply(directions, function(direction) {
  unlist(lapply(fuels, function(fuel) c(
    paste0("campd_", fuel, "_", direction, "_n_plants_", radii, "km_90deg"),
    paste0("campd_", fuel, "_", direction, "_heat_observed_", radii, "km_90deg")
  )), use.names = FALSE)
}), use.names = FALSE)

daily_files <- file.path(daily_dir, paste0("07H_fuel_daily_exposure_", years, ".csv.gz"))
assert_true(all(file.exists(daily_files)), "One or more of the 24 frozen daily exposure files is missing.")
daily <- rbindlist(lapply(daily_files, function(f) fread(f, select = c("loc_num", "date", exposure_cols), showProgress = FALSE)))
daily[, `:=`(loc_num = clean_chr(loc_num), date = as.IDate(date), daily_match = TRUE)]
for (cc in exposure_cols) set(daily, j = cc, value = safe_num(daily[[cc]]))
assert_true(daily[, anyDuplicated(paste(loc_num, date))] == 0L, "Duplicate monitor-date rows in daily exposure files.")

fp_days <- fp[, .(date = seq(collect_start, collect_end - 1L, by = "day")), by = .(period_id, loc_num)]
fp_days[, date := as.IDate(date)]
fp_days <- merge(fp_days, daily, by = c("loc_num", "date"), all.x = TRUE, sort = FALSE)
assert_true(fp_days[is.na(daily_match), .N] == 0L, "Some formal period-days did not match daily exposure data.")
fp_exp <- fp_days[, lapply(.SD, strict_mean), by = period_id, .SDcols = exposure_cols]
setnames(fp_exp, exposure_cols, paste0("fpmean_", exposure_cols))
fp <- merge(fp, fp_exp, by = "period_id", all.x = TRUE, sort = FALSE)

keep <- c("period_id", "loc_num", "radnet_monitor_id", "collect_start", "collect_end", "period_days",
          "calendar_year", "monitor_year_fe", "gross_beta_mBq_m3", "c_prop",
          "surface_u_ppm_article", "sunspot_monthly", "pm_filter_period_official",
          "pm_model_centered", "inv_hpbl", "rhum", "soilm", "vel", "lat_diff_35",
          "year_model", "temp", "temp_centered", "temp2", "temp3", paste0("fpmean_", exposure_cols))
dir.create(MODEL_DATA_DIR, recursive = TRUE, showWarnings = FALSE)
out <- file.path(MODEL_DATA_DIR, "filter_period_primary_2001_2024.csv.gz")
fwrite(fp[, ..keep], out, compress = "gzip")
cat("Wrote ", out, " (", nrow(fp), " rows; ", length(keep), " columns)\n", sep = "")
