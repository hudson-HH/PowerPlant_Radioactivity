#!/usr/bin/env Rscript

required_packages(c("data.table", "digest"))

spec <- list(
  list("main_figures/Figure4.csv", 56L, c("fuel","metric","radius_km","wind_reference","wind_ms","estimate","ci95_low","ci95_high","n_rows")),
  list("main_figures/Figure5.csv", 168L, c("fuel","metric","radius_km","component","wind_reference","wind_ms","estimate","ci95_low","ci95_high","n_rows")),
  list("main_figures/Figure6A.csv", 16L, c("fuel","metric","model","estimate","ci_low","ci_high")),
  list("main_figures/Figure6B.csv", 24L, c("fuel","metric","model","estimate","ci_low","ci_high")),
  list("main_figures/Figure6C.csv", 16L, c("fuel","metric","model","estimate","ci_low","ci_high")),
  list("sensitivity/03_current_and_natural_spline_effects.csv", 112L, c("fuel","metric","radius_km","model","N_model","wind_reference","wind_ms","estimate")),
  list("sensitivity/H21_biomass_GFED_full.csv", 14L, c("metric","radius_km","n_model_M0","n_model_M1","wind_reference","wind_ms","original_estimate","fire_adjusted_estimate"))
)

rows <- vector("list", length(spec))
for (i in seq_along(spec)) {
  rel <- spec[[i]][[1]]; expected_n <- spec[[i]][[2]]; cols <- spec[[i]][[3]]
  path <- file.path(SOURCE_DIR, rel); d <- read_csv(path); assert_columns(d, cols, rel)
  ok <- nrow(d) == expected_n && all(vapply(d[, ..cols], function(z) !all(is.na(z)), logical(1)))
  rows[[i]] <- data.table::data.table(file = rel, rows = nrow(d), expected_rows = expected_n,
                                      sha256 = digest::digest(file = path, algo = "sha256"), status = if (ok) "PASS" else "FAIL")
}

f4 <- read_csv(file.path(SOURCE_DIR, "main_figures", "Figure4.csv"))
assert_true(all(f4$wind_reference == "sample_median"), "Figure 4 source must use sample-median wind.")
f6b <- read_csv(file.path(SOURCE_DIR, "main_figures", "Figure6B.csv"))
f6wind <- read_csv(file.path(SOURCE_DIR, "display_v3", "Figure_6_reference_winds_and_samples.csv"))
f6check <- merge(f6b[, .(fuel, metric, wind_ms)], f6wind[, .(fuel, metric, expected_wind=wind_ms)],
                 by=c("fuel", "metric"), all.x=TRUE)
assert_true(nrow(f6check)==24L && all(is.finite(f6check$expected_wind)) &&
              all(abs(f6check$wind_ms-f6check$expected_wind)<1e-12),
            "Figure 6B must use the shared fuel-by-metric baseline sample median.")
h21 <- read_csv(file.path(SOURCE_DIR, "sensitivity", "H21_biomass_GFED_full.csv"))
assert_true(all(h21$wind_reference == "sample median wind"), "H21 must be labelled sample-median wind.")
assert_true(all(h21$n_model_M0 == h21$n_model_M1), "H21 M0/M1 samples must match.")

rows[[length(rows) + 1L]] <- data.table::data.table(file = "semantic_checks", rows = NA_integer_, expected_rows = NA_integer_,
  sha256 = NA_character_, status = "PASS",
  note = "F4/H21 use sample-median wind; F6B uses shared baseline sample median; H21 M0/M1 n are equal")
write_validation(rows, "source_data_validation.csv")
assert_true(all(vapply(rows, function(x) x$status[[1]] == "PASS", logical(1))), "Source-data validation failed.")
cat("Source-data validation: PASS\n")
