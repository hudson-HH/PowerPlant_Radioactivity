options(stringsAsFactors = FALSE, warn = 1, scipen = 999)

script_file <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep("^--file=", args, value = TRUE)
  if (length(hit)) normalizePath(sub("^--file=", "", hit[[1]]), mustWork = TRUE) else NA_character_
}

package_root <- function() {
  override <- Sys.getenv("RADNET_CODE_PACKAGE", unset = "")
  if (nzchar(override)) return(normalizePath(override, mustWork = TRUE))
  sf <- script_file()
  if (!is.na(sf)) return(normalizePath(file.path(dirname(sf), ".."), mustWork = TRUE))
  normalizePath(getwd(), mustWork = TRUE)
}

ROOT <- package_root()
SOURCE_DIR <- file.path(ROOT, "data", "source_data")
MODEL_DATA_DIR <- file.path(ROOT, "data", "prepared_model_data")
OUTPUT_DIR <- file.path(ROOT, "outputs", "v1_intermediates")
VALIDATION_DIR <- file.path(ROOT, "validation")

dir.create(file.path(OUTPUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(VALIDATION_DIR, recursive = TRUE, showWarnings = FALSE)

required_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing R packages: ", paste(missing, collapse = ", "),
         ". Install them before running; this package never installs software automatically.",
         call. = FALSE)
  }
}

read_csv <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
  data.table::fread(path, showProgress = FALSE, check.names = TRUE)
}

assert_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
  invisible(TRUE)
}

write_validation <- function(rows, filename) {
  data.table::fwrite(data.table::rbindlist(rows, fill = TRUE), file.path(VALIDATION_DIR, filename))
}

fuel_labels <- c(coal = "Coal", gas = "Natural gas", oil = "Oil", biomass = "Biomass")
metric_labels <- c(count = "CAMPD plant count", heat = "Heat input")
