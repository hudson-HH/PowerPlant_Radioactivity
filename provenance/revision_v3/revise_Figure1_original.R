#!/usr/bin/env Rscript
# Descriptive map revision only; does not change any analytical exposure or model.
options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(maps)
  library(patchwork)
  library(grid)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
project <- Sys.getenv("RADNET_PROJECT", "<ORIGINAL_PROJECT_ROOT>")
input_dir <- if (length(args) >= 1L) args[[1]] else file.path(project, "05_Deliverables/2026-09-29_Nature_manuscript/source_data")
output_dir <- if (length(args) >= 2L) args[[2]] else file.path(project, "05_Deliverables/2026-10-01_Figure1与样本风速核对")
qa_dir <- if (length(args) >= 3L) args[[3]] else file.path(project, "07_Logs/2026-10-01_Figure1与样本风速核对")
skill_root <- Sys.getenv("NATURE_FIGURE_SKILL", "<ORIGINAL_FIGURE_SKILL>")
qa_python <- Sys.getenv("NATURE_FIGURE_QA_PYTHON", "<ORIGINAL_QA_PYTHON>")
for (pkg in c("ragg", "svglite")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required R export package is missing: ", pkg)
}
source(file.path(skill_root, "scripts/panel_alignment.R"))
for (d in c(output_dir, qa_dir, file.path(output_dir, "source_data"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)

plant_file <- file.path(input_dir, "Figure_1A_powerplant_locations.csv")
monitor_file <- file.path(input_dir, "Figure_1_radnet_monitor_locations.csv")
stopifnot(file.exists(plant_file), file.exists(monitor_file))
P <- fread(plant_file)
M <- fread(monitor_file)
stopifnot(all(c("plant_id_plot", "plant_longitude", "plant_latitude", "fuel_label") %in% names(P)))
stopifnot(all(c("loc_num", "monitor_longitude", "monitor_latitude") %in% names(M)))
stopifnot(all(is.finite(P$plant_longitude)), all(is.finite(P$plant_latitude)))
stopifnot(all(is.finite(M$monitor_longitude)), all(is.finite(M$monitor_latitude)))
stopifnot(!anyDuplicated(P[, .(plant_id_plot, fuel_label)]), !anyDuplicated(M$loc_num))
stopifnot(P[, uniqueN(paste(plant_longitude, plant_latitude)), by = plant_id_plot][V1 > 1L, .N] == 0L)

fuel_levels <- c("Coal", "Natural gas", "Oil", "Other fossil", "Biomass", "Solid waste", "Nuclear", "Geothermal", "Wind", "Solar", "Hydro", "Other")
combustion <- c("Coal", "Natural gas", "Oil", "Other fossil", "Biomass", "Solid waste")
noncombustion <- c("Nuclear", "Geothermal", "Wind", "Solar", "Hydro")
stopifnot(all(P$fuel_label %in% fuel_levels))

# One row per physical plant; all roles remain in panel a.
# Any recorded combustion role wins in panel b, including mixed-energy plants.
# An unknown role without a combustion role is not silently called non-combustion.
B <- P[, .(
  plant_longitude = first(plant_longitude),
  plant_latitude = first(plant_latitude),
  all_energy_roles = paste(sort(unique(fuel_label)), collapse = "|"),
  has_combustion_role = any(fuel_label %in% combustion),
  has_noncombustion_role = any(fuel_label %in% noncombustion),
  has_other_role = any(!fuel_label %in% c(combustion, noncombustion))
), by = plant_id_plot]
B[, process_class := fifelse(has_combustion_role, "Combustion-capable",
  fifelse(has_noncombustion_role & !has_other_role, "Non-combustion only", "Other / unclassified"))]
stopifnot(nrow(B) == uniqueN(P$plant_id_plot))
stopifnot(B[, .N] == sum(B[, .N, by = process_class]$N))

copy_dest <- file.path(output_dir, "source_data")
stopifnot(file.copy(plant_file, file.path(copy_dest, basename(plant_file)), overwrite = TRUE))
stopifnot(file.copy(monitor_file, file.path(copy_dest, basename(monitor_file)), overwrite = TRUE))
fwrite(B, file.path(copy_dest, "Figure_1B_combustion_based_classification.csv"))
fwrite(P[, .(physical_plants = uniqueN(plant_id_plot)), by = fuel_label], file.path(copy_dest, "Figure_1A_energy_role_counts.csv"))
fwrite(B[, .(physical_plants = .N), by = process_class], file.path(copy_dest, "Figure_1B_process_class_counts.csv"))

monitor_label <- "RadNet location"
palette_a <- c("Coal" = "#4D4D4D", "Natural gas" = "#D39D2B", "Oil" = "#D77B45",
  "Other fossil" = "#7E668D", "Biomass" = "#D81B60", "Solid waste" = "#478F5F",
  "Nuclear" = "#784599", "Geothermal" = "#B287BD", "Wind" = "#74B8CA",
  "Solar" = "#E5D4A6", "Hydro" = "#508FCC", "Other" = "#A6A6A6", "RadNet location" = "#111111")
levels_a <- c(fuel_levels, monitor_label)
shape_a <- setNames(rep(16, length(levels_a)), levels_a)
shape_a[c("Biomass", monitor_label)] <- c(18, 17)

draw_order_a <- c("Solar", "Other", "Wind", "Hydro", "Geothermal", "Natural gas", "Oil", "Other fossil", "Coal", "Nuclear", "Solid waste", "Biomass")
AP <- P[, .(longitude = plant_longitude, latitude = plant_latitude, category = fuel_label)]
AP[, `:=`(draw_order = match(category, draw_order_a), alpha_draw = 0.55, size_draw = 0.43)]
AP[category == "Solar", `:=`(alpha_draw = 0.68, size_draw = 0.30)]
AP[category == "Biomass", `:=`(alpha_draw = 0.96, size_draw = 0.77)]
setorder(AP, draw_order)
MP <- M[, .(longitude = monitor_longitude, latitude = monitor_latitude,
  category = monitor_label, draw_order = 100L, alpha_draw = 1, size_draw = 1.20)]
AD <- rbindlist(list(AP, MP))
AD[, category := factor(category, levels = levels_a)]
stopifnot(nrow(AD) == nrow(P) + nrow(M))

levels_b <- c("Combustion-capable", "Non-combustion only", "Other / unclassified", monitor_label)
palette_b <- c("Combustion-capable" = "#D55E00", "Non-combustion only" = "#A6CADB",
  "Other / unclassified" = "#888888", "RadNet location" = "#111111")
shape_b <- setNames(c(16, 16, 16, 17), levels_b)
BP <- B[, .(longitude = plant_longitude, latitude = plant_latitude, category = process_class)]
BP[, `:=`(draw_order = match(category, c("Non-combustion only", "Other / unclassified", "Combustion-capable")), alpha_draw = 0.60, size_draw = 0.40)]
BP[category == "Combustion-capable", `:=`(alpha_draw = 0.75, size_draw = 0.46)]
setorder(BP, draw_order)
BD <- rbindlist(list(BP, MP))
BD[, category := factor(category, levels = levels_b)]
stopifnot(nrow(BD) == nrow(B) + nrow(M))

US_MAP <- map_data("state")
map_theme <- theme_void(base_size = 6.5, base_family = "Helvetica") + theme(
  plot.title = element_text(size = 8, face = "bold", hjust = 0, margin = margin(b = 3)),
  legend.position = "right", legend.title = element_blank(),
  legend.text = element_text(size = 6.2),
  legend.key.height = unit(3.1, "mm"), legend.key.width = unit(3.6, "mm"),
  legend.margin = margin(0, 0, 0, 0),
  legend.box.spacing = unit(2.3, "mm"),
  plot.margin = margin(2, 2, 2, 2)
)

draw_map <- function(data, palette, shapes, categories, title) {
  legend_sizes <- ifelse(categories == monitor_label, 1.6, ifelse(categories == "Biomass", 1.8, 1.5))
  ggplot() +
    geom_polygon(data = US_MAP, aes(x = long, y = lat, group = group),
      fill = "white", colour = "#CCCCCC", linewidth = 0.18) +
    geom_point(data = data, aes(x = longitude, y = latitude, colour = category,
      shape = category, alpha = alpha_draw, size = size_draw)) +
    coord_quickmap(xlim = c(-125, -66), ylim = c(24, 50), expand = FALSE) +
    scale_colour_manual(values = palette, breaks = categories, limits = categories, drop = FALSE) +
    scale_shape_manual(values = shapes, breaks = categories, limits = categories, drop = FALSE) +
    scale_alpha_identity(guide = "none") + scale_size_identity(guide = "none") +
    guides(colour = guide_legend(order = 1, override.aes = list(shape = unname(shapes[categories]),
      alpha = 1, size = legend_sizes)), shape = "none") +
    labs(title = title) + map_theme
}
pA <- draw_map(AD, palette_a, shape_a, levels_a, "a  Energy types (non-exclusive)")
pB <- draw_map(BD, palette_b, shape_b, levels_b, "b  Combustion-based plant classes")
figure <- pA / pB + plot_layout(heights = c(1, 1)) + plot_annotation(
  title = "Power-plant inventory and RadNet locations",
  theme = theme(plot.title = element_text(size = 9, face = "bold", family = "Helvetica", hjust = 0,
    margin = margin(b = 3)), plot.margin = margin(3, 3, 3, 3)))
width_mm <- 183
height_mm <- 146
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

# Fixed-aspect maps have null gtable units; unresolved unit sums are not actual
# panel rectangles. Measure the fully drawn PDF viewports, as in this project's
# validated rendering workflow, then apply the same strict alignment auditor.
probe_path <- file.path(qa_dir, "Figure_1_geometry_probe.pdf")
grDevices::pdf(probe_path, width = width_in, height = height_in, family = "Helvetica", useDingbats = FALSE)
grob <- patchworkGrob(figure)
grid.newpage(); grid.draw(grob); grid.force()
vps <- grid.ls(viewports = TRUE, grobs = FALSE, print = FALSE)
hits <- which(vps$type == "vpListing" & grepl("^panel(?:-[0-9]+)*[.]", vps$name, perl = TRUE))
stopifnot(length(hits) == 2L)
boxes <- lapply(hits, function(j) {
  path <- c(strsplit(vps$vpPath[j], "::", fixed = TRUE)[[1]], vps$name[j])
  path <- path[path != "ROOT" & nzchar(path)]
  upViewport(0); downViewport(do.call(vpPath, as.list(path)))
  a <- deviceLoc(unit(0, "npc"), unit(0, "npc"), valueOnly = TRUE)
  b <- deviceLoc(unit(1, "npc"), unit(1, "npc"), valueOnly = TRUE)
  c(a$x, a$y, b$x, b$y) * 72
})
upViewport(0); dev.off()
rectangles <- do.call(rbind, boxes)
stopifnot(is.matrix(rectangles), nrow(rectangles) == 2L, all(rectangles[, 3] > rectangles[, 1]), all(rectangles[, 4] > rectangles[, 2]))
rectangles <- rectangles[order(-rectangles[, 4]), , drop = FALSE]
layout_manifest <- list(schema_version = 1L, backend = "r-forced-PDF-viewport",
  figure = list(width_pt = width_in * 72, height_pt = height_in * 72),
  panels = lapply(1:2, function(i) list(id = letters[i], bbox_pt = rectangles[i, ])),
  column_groups = list(c("a", "b")), exemptions = list())
layout_path <- file.path(qa_dir, "Figure_1.alignment-layout.json")
write_json(layout_manifest, layout_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
audit_path <- file.path(qa_dir, "Figure_1.alignment.json")
alignment_status <- system2(qa_python, c(shQuote(file.path(skill_root, "scripts/audit_panel_alignment.py")),
  shQuote(layout_path), "--json-out", shQuote(audit_path), "--overlay-svg",
  shQuote(file.path(qa_dir, "Figure_1.alignment.svg")), "--strict"))
if (alignment_status != 0L) stop("Panel alignment failed; see ", audit_path)

stem <- file.path(output_dir, "Figure_1_powerplants_RadNet_map_revised")
grDevices::pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Helvetica", useDingbats = FALSE)
print(figure); dev.off()
svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
print(figure); dev.off()
ragg::agg_png(paste0(stem, ".png"), width = width_mm, height = height_mm, units = "mm", res = 600, background = "white")
print(figure); dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = width_mm, height = height_mm, units = "mm", res = 600, compression = "lzw", background = "white")
print(figure); dev.off()
# R-native QA previews, not replacement scientific panels.
for (i in seq_along(list(pA, pB))) {
  ragg::agg_png(file.path(qa_dir, paste0("Figure_1_panel_", letters[i], "_preview.png")),
    width = width_mm, height = 74, units = "mm", res = 300, background = "white")
  print(list(pA, pB)[[i]]); dev.off()
}

data_audit <- list(
  input_energy_role_rows = nrow(P), input_physical_plants = uniqueN(P$plant_id_plot),
  input_monitor_rows = nrow(M), input_unique_monitors = uniqueN(M$loc_num),
  excluded_source_rows = 0L, panel_a_plant_rows = nrow(AP), panel_b_plant_rows = nrow(BP),
  mixed_combustion_noncombustion_plants = B[has_combustion_role & has_noncombustion_role, .N],
  unknown_and_noncombustion_without_combustion = B[has_other_role & has_noncombustion_role & !has_combustion_role, .N],
  process_class_counts = B[, .(physical_plants = .N), by = process_class],
  width_mm = width_mm, height_mm = height_mm, raster_dpi = 600,
  drawing_backend = "R", R_version = as.character(getRversion()),
  source_md5 = as.list(tools::md5sum(c(plant_file, monitor_file))),
  colour_a = as.list(palette_a), colour_b = as.list(palette_b))
write_json(data_audit, file.path(qa_dir, "Figure_1_data_audit.json"), auto_unbox = TRUE, pretty = TRUE)
print(data_audit[c("input_energy_role_rows", "input_physical_plants", "input_unique_monitors", "process_class_counts")])
cat("Figure1 revision rendered. Final PDF text/collision audits and visual review are still required.\n")
