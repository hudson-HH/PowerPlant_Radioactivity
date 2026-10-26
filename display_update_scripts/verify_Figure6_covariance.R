#!/usr/bin/env Rscript
# Independent contrast-vector verification; does not refit any model.
# CLI: Rscript verify_Figure6_covariance.R results.rds report.json [56row_reference.csv]
verify_figure6_covariance <- function(rds_file, report_file, reference_csv = NULL) {
  stopifnot(requireNamespace("data.table", quietly = TRUE), requireNamespace("jsonlite", quietly = TRUE))
  z <- readRDS(rds_file)
  errors <- numeric(); convergence <- list(); rows <- list()
  for (key in names(z$covariance_records)) {
    s <- z$baseline_samples[[key]]; records <- z$covariance_records[[key]]
    maps <- data.table::data.table(panel = c("A", "A", "B", "B", "B", "C", "C"),
      fit = c("main", "strict", "original", "mundlak", "mundlak", "main", "nopm"),
      term = c("exposure", "exposure", "exposure", "exposure_within", "exposure_between", "exposure", "exposure"),
      model = c("Monitor × year FE", "Monitor × year + common-time FE", "Original random-intercept",
        "Within monitor-year", "Between monitor-years", "PM₂.₅ adjusted", "No PM₂.₅ adjustment"))
    for (i in seq_len(nrow(maps))) {
      record <- records[[maps$fit[i]]]; term <- maps$term[i]
      a <- setNames(rep(0, length(record$coefficients)), names(record$coefficients))
      a[term] <- 1; a[paste0(term, ":vel")] <- s$wind_ms
      est <- sum(a * record$coefficients)
      se <- sqrt(drop(t(a) %*% record$covariance %*% a))
      row <- data.table::data.table(panel = maps$panel[i], fuel = s$fuel, metric = s$metric,
        model = maps$model[i], estimate = est, ci_low = est - 1.96 * se, ci_high = est + 1.96 * se,
        N_model = record$N_model, N_pre = s$N_pre, N_monitors = s$N_monitors, wind_ms = s$wind_ms)
      rows[[length(rows) + 1L]] <- row
    }
    for (m in c("original", "mundlak")) {
      convergence[[paste(key, m, sep = "_")]] <- list(optimizer_code = records[[m]]$convergence$conv$opt,
        lme4_messages = records[[m]]$convergence$conv$lme4,
        singular = records[[m]]$random_intercept_singular)
    }
  }
  reconstructed <- data.table::rbindlist(rows)
  ref <- if (is.null(reference_csv)) z$results else data.table::fread(reference_csv)
  compare <- merge(reconstructed, ref, by = c("panel", "fuel", "metric", "model"), suffixes = c("_vector", "_reference"), all = TRUE)
  compare[, max_abs_error := pmax(abs(estimate_vector - estimate_reference),
    abs(ci_low_vector - ci_low_reference), abs(ci_high_vector - ci_high_reference))]
  max_error <- max(compare$max_abs_error)
  stopifnot(nrow(compare) == 56L, is.finite(max_error), max_error < 1e-12,
    all(compare$N_model_vector == compare$N_model_reference),
    all(compare$N_pre_vector == compare$N_pre_reference),
    all(abs(compare$wind_ms_vector - compare$wind_ms_reference) < 1e-12))
  report <- list(method = "Independent full contrast vector a: estimate=a^T beta; variance=a^T V a",
    compared_rows = 56L, max_abs_error = max_error, baseline_matrix_md5 = z$matrix_md5,
    covariance_terms_retained = TRUE, mixed_models = convergence, verdict = "PASS")
  jsonlite::write_json(report, report_file, auto_unbox = TRUE, pretty = TRUE, digits = NA, null = "null")
  cat("Independent covariance-vector verification: PASS; max error ", format(max_error), "\n", sep = "")
  invisible(reconstructed)
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2L) stop("Usage: Rscript verify_Figure6_covariance.R results.rds report.json [reference.csv]")
  verify_figure6_covariance(args[1], args[2], if (length(args) > 2L) args[3] else NULL)
}
