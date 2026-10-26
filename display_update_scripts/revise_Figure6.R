#!/usr/bin/env Rscript
# Figure 6: identical prepared data, models and fuel-by-metric baseline samples.
# Only panel B's evaluation wind changes; all A/C published values are retained.
# Usage: Rscript revise_Figure6.R [fit|plot|all] [project_root]
options(stringsAsFactors = FALSE, warn = 1, scipen = 999)
args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args)) args[1] else "all"
stopifnot(mode %in% c("fit", "plot", "all"))
sf <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
source(file.path(dirname(normalizePath(sf, mustWork = TRUE)), "..", "R", "common.R"))
project <- if (length(args) > 1L) normalizePath(args[2], mustWork = TRUE) else ROOT
version <- "figure6"
old <- file.path(project, "provenance", "revision_v3", "legacy_1ms")
work <- file.path(project, "outputs", "display_updates", version, "working")
logs <- file.path(project, "validation", "display_updates", version)
delivery <- file.path(project, "outputs", "display_updates", version)
source_dir <- file.path(delivery, "source_data")
figure_dir <- file.path(delivery, "figures")
for (p in c(work, logs, source_dir, figure_dir)) dir.create(p, recursive = TRUE, showWarnings = FALSE)
needed <- c("data.table", "fixest", "lme4", "jsonlite")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(data.table); library(fixest)})
matrix_file <- file.path(project, "data", "prepared_model_data", "filter_period_primary_2001_2024.csv.gz")
frozen_file <- file.path(old, "Figure_6_panel_sample_sizes_and_estimates.csv")
matrix_md5 <- unname(tools::md5sum(matrix_file))
frozen <- fread(frozen_file)
stopifnot(nrow(frozen) == 56L, uniqueN(frozen, by = c("panel", "fuel", "metric", "model")) == 56L)
expected_winds <- c(coal_count = 3.15442053456934, gas_count = 3.15442053456934,
  oil_count = 3.15442053456934, biomass_count = 3.15442053456934,
  coal_heat = 3.15275005725255, gas_heat = 3.13870352480861,
  oil_heat = 3.13587557596978, biomass_heat = 3.15442053456934)
env <- c("c_prop", "sunspot_monthly", "pm_model_centered", "inv_hpbl", "rhum", "soilm", "vel", "temp_centered", "temp2", "temp3", "sunspot_monthly:lat_diff_35", "surface_u_ppm_article:soilm")
mixed_env <- c(env, "surface_u_ppm_article", "lat_diff_35", "year_model", "I(year_model^2)")
opts <- list(optimizer = "bobyqa", optCtrl = list(maxfun = 200000),
  check.conv.singular = "ignore", check.rankX = "message+drop.cols")
if ("autoscale" %in% names(formals(lme4::lmerControl))) opts$autoscale <- TRUE
control <- do.call(lme4::lmerControl, opts)
model_record <- function(fit, mixed = FALSE) {
  list(coefficients = if (mixed) lme4::fixef(fit) else coef(fit),
    covariance = as.matrix(vcov(fit)), N_model = nobs(fit), formula = formula(fit),
    covariance_definition = if (mixed) "lme4 fixed-effect covariance; REML random-intercept model" else "fixest loc_num-clustered covariance",
    convergence = if (mixed) fit@optinfo else NULL,
    random_intercept_singular = if (mixed) lme4::isSingular(fit) else NULL)
}
slope <- function(record, term, wind) {
  b <- record$coefficients; V <- record$covariance; ti <- paste0(term, ":vel")
  stopifnot(all(c(term, ti) %in% names(b)), all(is.finite(V[c(term, ti), c(term, ti)])))
  estimate <- b[term] + wind * b[ti]
  variance <- V[term, term] + wind^2 * V[ti, ti] + 2 * wind * V[term, ti]
  stopifnot(is.finite(variance), variance >= 0)
  se <- sqrt(variance)
  data.table(estimate = unname(estimate), ci_low = unname(estimate - 1.96 * se),
    ci_high = unname(estimate + 1.96 * se), N_model = record$N_model)
}
num_file <- file.path(work, "Figure6_recomputed_results.rds")
if (mode %in% c("fit", "all")) {
  cat("Reading prepared matrix: ", matrix_file, "\n", sep = "")
  d <- fread(matrix_file)
  d[, time_fe := interaction(collect_start, collect_end, drop = TRUE, lex.order = TRUE)]
  stopifnot(nrow(d) == 107407L, uniqueN(d$loc_num) == 133L)
  grid <- CJ(fuel = c("coal", "gas", "oil", "biomass"), metric = c("count", "heat"))
  legacy_out <- new_out <- list(); records <- list(); samples <- list()
  for (i in seq_len(nrow(grid))) {
    fuel <- grid$fuel[i]; metric <- grid$metric[i]; key <- paste(fuel, metric, sep = "_")
    cache <- file.path(work, paste0("model_covariance_", key, ".rds"))
    ex <- paste0("fpmean_campd_", fuel, "_upwind_", if (metric == "count") "n_plants" else "heat_observed", "_20km_90deg")
    x <- copy(d[is.finite(get(ex))])
    x[, exposure := get(ex) / if (metric == "heat") 10000 else 1]
    wind <- median(x$vel)
    stopifnot(abs(wind - expected_winds[key]) < 1e-12)
    sample <- list(fuel = fuel, metric = metric, exposure_column = ex,
      N_pre = nrow(x), N_monitors = uniqueN(x$loc_num), wind_ms = wind,
      wind_definition = "Median filter-period vel in the fuel-by-metric exposure-finite baseline sample, before FE singleton removal; shared across A/B/C models",
      period_ids = x$period_id, matrix_md5 = matrix_md5)
    cat(format(Sys.time()), " START ", key, " N_pre=", nrow(x), " wind=", format(wind, digits = 16), "\n", sep = "")
    if (file.exists(cache)) {
      z <- readRDS(cache)
      stopifnot(identical(z$sample$matrix_md5, matrix_md5), identical(z$sample$period_ids, x$period_id), abs(z$sample$wind_ms - wind) < 1e-12)
      rec <- z$models
      cat("  Reusing verified model coefficient/covariance checkpoint\n")
    } else {
      main <- feols(as.formula(paste0("gross_beta_mBq_m3 ~ exposure + exposure:vel + ", paste(env, collapse = " + "), " | monitor_year_fe")), data = x, vcov = ~loc_num, notes = FALSE)
      strict <- feols(as.formula(paste0("gross_beta_mBq_m3 ~ exposure + exposure:vel + ", paste(env, collapse = " + "), " | monitor_year_fe + time_fe")), data = x, vcov = ~loc_num, notes = FALSE)
      nopm <- feols(as.formula(paste0("gross_beta_mBq_m3 ~ exposure + exposure:vel + ", paste(setdiff(env, "pm_model_centered"), collapse = " + "), " | monitor_year_fe")), data = x, vcov = ~loc_num, notes = FALSE)
      cat("  FE models completed\n")
      x[, loc_num_factor := factor(loc_num)]
      original <- lme4::lmer(as.formula(paste0("gross_beta_mBq_m3 ~ exposure + exposure:vel + ", paste(mixed_env, collapse = " + "), " + (1 | loc_num_factor)")), data = x, REML = TRUE, na.action = na.fail, control = control)
      cat("  Original random-intercept completed\n")
      x[, exposure_between := mean(exposure), by = monitor_year_fe]
      x[, exposure_within := exposure - exposure_between]
      mundlak <- lme4::lmer(as.formula(paste0("gross_beta_mBq_m3 ~ exposure_within + exposure_between + exposure_within:vel + exposure_between:vel + ", paste(mixed_env, collapse = " + "), " + (1 | loc_num_factor)")), data = x, REML = TRUE, na.action = na.fail, control = control)
      cat("  Within/between random-intercept completed\n")
      rec <- list(main = model_record(main), strict = model_record(strict), nopm = model_record(nopm),
        original = model_record(original, TRUE), mundlak = model_record(mundlak, TRUE))
      saveRDS(list(sample = sample, models = rec, R_session = sessionInfo()), cache)
      rm(main, strict, nopm, original, mundlak); gc(FALSE)
    }
    a <- rbind(slope(rec$main, "exposure", wind)[, model := "Monitor × year FE"],
      slope(rec$strict, "exposure", wind)[, model := "Monitor × year + common-time FE"])[, panel := "A"]
    cc <- rbind(slope(rec$main, "exposure", wind)[, model := "PM₂.₅ adjusted"],
      slope(rec$nopm, "exposure", wind)[, model := "No PM₂.₅ adjustment"])[, panel := "C"]
    b_at <- function(w) rbind(slope(rec$original, "exposure", w)[, model := "Original random-intercept"],
      slope(rec$mundlak, "exposure_within", w)[, model := "Within monitor-year"],
      slope(rec$mundlak, "exposure_between", w)[, model := "Between monitor-years"])[, panel := "B"]
    legacy <- rbind(a, b_at(1), cc)
    updated <- rbind(a, b_at(wind), cc)
    for (r in list(legacy, updated)) r[, `:=`(fuel = fuel, metric = metric, N_pre = nrow(x), N_monitors = uniqueN(x$loc_num))]
    legacy[, wind_ms := ifelse(panel == "B", 1, wind)]
    updated[, wind_ms := wind]
    legacy_out[[i]] <- legacy; new_out[[i]] <- updated; records[[key]] <- rec; samples[[key]] <- sample
    fwrite(rbindlist(new_out), file.path(logs, "Figure6_new_median_progress.csv"))
    cat(format(Sys.time()), " FINISHED ", key, "\n", sep = "")
    print(updated[panel == "B", .(fuel, metric, model, estimate, ci_low, ci_high, N_model, wind_ms)])
  }
  legacy <- rbindlist(legacy_out); recomputed <- rbindlist(new_out)
  check <- merge(legacy, frozen[, .(panel, fuel, metric, model, frozen_estimate = estimate,
    frozen_low = ci_low, frozen_high = ci_high, frozen_N_model = N_model,
    frozen_N_pre = N_pre, frozen_N_monitors = N_monitors, frozen_wind = wind_ms)],
    by = c("panel", "fuel", "metric", "model"), all = TRUE)
  check[, max_abs_diff := pmax(abs(estimate - frozen_estimate), abs(ci_low - frozen_low), abs(ci_high - frozen_high))]
  check[, tolerance := ifelse(panel == "B", 1e-5, 1e-7)]
  check[, status := ifelse(is.finite(max_abs_diff) & max_abs_diff <= tolerance &
    N_model == frozen_N_model & N_pre == frozen_N_pre & N_monitors == frozen_N_monitors & abs(wind_ms - frozen_wind) < 1e-12, "PASS", "FAIL")]
  fwrite(check, file.path(logs, "Figure6_legacy_1ms_refit_vs_frozen_56rows.csv"))
  stopifnot(nrow(check) == 56L, all(check$status == "PASS"))
  # A/C are copied only after their independent numerical and sample checks pass.
  results <- copy(frozen)
  b <- recomputed[panel == "B"]
  results[b, on = .(panel, fuel, metric, model), `:=`(estimate = i.estimate,
    ci_low = i.ci_low, ci_high = i.ci_high, wind_ms = i.wind_ms)]
  stopifnot(identical(results[panel != "B"], frozen[panel != "B"]))
  results[, wind_reference := "Shared fuel-by-metric baseline sample median"]
  results[, wind_definition := "Median filter-period vel in exposure-finite baseline sample before FE singleton removal; not each fitted model's independent median"]
  fwrite(results, file.path(source_dir, "Figure_6_panel_sample_sizes_and_estimates.csv"))
  fwrite(results, file.path(source_dir, "Figure_6_internal_validity_robustness.csv"))
  for (p in c("A", "B", "C")) {
    suffix <- c(A = "main_vs_strict_time", B = "within_between", C = "PM25_sensitivity")[p]
    fwrite(results[panel == p], file.path(source_dir, paste0("Figure_6", p, "_", suffix, ".csv")))
  }
  oldb <- frozen[panel == "B"]
  comparison <- merge(results[panel == "B"], oldb[, .(fuel, metric, model,
    old_estimate = estimate, old_ci_low = ci_low, old_ci_high = ci_high, old_wind_ms = wind_ms)], by = c("fuel", "metric", "model"))
  direction <- function(lo, hi) fifelse(lo > 0, "positive_CI_excludes_zero", fifelse(hi < 0, "negative_CI_excludes_zero", "CI_includes_zero"))
  comparison[, `:=`(old_CI_direction = direction(old_ci_low, old_ci_high), new_CI_direction = direction(ci_low, ci_high))]
  comparison[, CI_interpretation_changed := old_CI_direction != new_CI_direction]
  fwrite(comparison, file.path(source_dir, "Figure_6B_1ms_vs_shared_median_comparison.csv"))
  sample_table <- rbindlist(lapply(samples, function(s) as.data.table(s[setdiff(names(s), "period_ids")])))
  fwrite(sample_table, file.path(source_dir, "Figure_6_reference_winds_and_samples.csv"))
  saveRDS(list(results = results, recomputed = recomputed, legacy = legacy, covariance_records = records,
    baseline_samples = samples, source_matrix = matrix_file, matrix_md5 = matrix_md5,
    frozen_file = frozen_file, session = sessionInfo()), num_file)
  writeLines(capture.output(sessionInfo()), file.path(logs, "Figure6_sessionInfo.txt"))
  jsonlite::write_json(list(matrix_file = matrix_file, matrix_md5 = matrix_md5,
    frozen_file = frozen_file, frozen_md5 = unname(tools::md5sum(frozen_file)),
    validation_rows = nrow(check), validation_pass = sum(check$status == "PASS"),
    max_abs_diff_FE = max(check[panel != "B", max_abs_diff]), max_abs_diff_mixed = max(check[panel == "B", max_abs_diff]),
    A_C_preserved_exactly = TRUE, B_evaluation_only_changed = TRUE,
    within_between_definition = "Exposure is split by monitor-year mean and deviation; each component interacts with the same observed filter-period vel, not a separately decomposed wind variable",
    confidence_interval = "beta_X + w*beta_X:vel +/- 1.96*sqrt(V_XX + w^2*V_interaction_interaction + 2*w*V_X_interaction)"),
    file.path(logs, "Figure6_numerical_validation.json"), pretty = TRUE, auto_unbox = TRUE)
  cat("NUMERICAL VALIDATION: all 56 legacy rows PASS; A/C retained exactly; B re-evaluated at shared median.\n")
  print(comparison[, .(fuel, metric, model, estimate, ci_low, ci_high, N_model, wind_ms, CI_interpretation_changed)])
}

if (mode %in% c("plot", "all")) {
  numerical <- readRDS(num_file)
  stopifnot(identical(numerical$matrix_md5, matrix_md5))
  source(file.path(dirname(normalizePath(sf)), "Figure6_display_v3.R"))
  render_figure6_v3(numerical$results,
    output_stem = file.path(figure_dir, "Figure_6_internal_validity_robustness"),
    qa_dir = logs, work_dir = work)
}
