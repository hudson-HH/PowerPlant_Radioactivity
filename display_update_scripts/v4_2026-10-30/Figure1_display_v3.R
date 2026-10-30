# Adapted from the supplied Figure1_display_v3.R; display-only combustion subset.
P_all <- fread('source_data_manuscript_v5/Figure_1A_powerplant_locations.csv')
M <- fread('source_data_manuscript_v5/Figure_1_radnet_monitor_locations.csv')
combustion <- c('Coal','Natural gas','Oil','Other fossil','Biomass','Solid waste')
P <- P_all[fuel_label %in% combustion]
stopifnot(!anyDuplicated(P[,.(plant_id_plot,fuel_label)]))
fwrite(P,'data/display_assets/v4_2026-10-30/source_data/Figure_1_combustion_locations.csv')
counts <- P[,.(unique_plants=uniqueN(plant_id_plot),plant_role_points=.N),by=fuel_label]
counts <- counts[match(combustion,fuel_label)]
fwrite(counts,'data/display_assets/v4_2026-10-30/source_data/Figure_1_combustion_counts.csv')
jsonlite::write_json(list(unique_plants=uniqueN(P$plant_id_plot),plant_role_points=nrow(P),
 counts=setNames(as.list(counts$unique_plants),counts$fuel_label),radnet_locations=nrow(M),
 original_role_rows=nrow(P_all),excluded_noncombustion_role_rows=nrow(P_all)-nrow(P)),
 'data/display_assets/v4_2026-10-30/source_data/Figure_1_counts.json',auto_unbox=TRUE,pretty=TRUE)
monitor_label <- "RadNet location"
palette_a <- c("Coal" = "#4D4D4D", "Natural gas" = "#D39D2B", "Oil" = "#D77B45",
  "Other fossil" = "#7E668D", "Biomass" = "#D81B60", "Solid waste" = "#478F5F",
  "Nuclear" = "#784599", "Geothermal" = "#B287BD", "Wind" = "#74B8CA",
  "Solar" = "#E5D4A6", "Hydro" = "#508FCC", "Other" = "#A6A6A6", "RadNet location" = "#111111")
levels_a <- c(combustion, monitor_label)
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

pA <- draw_map(AD,palette_a,shape_a,levels_a,NULL)
stopifnot(nrow(AD)==nrow(P)+nrow(M),all(AP$category %in% combustion))
# Keep the original width; remove the former second map and title allocation.
export_figure(pA,'Figure_1',183/25.4,74/25.4,baseline=NULL)
fwrite(AD,'data/display_assets/v4_2026-10-30/qa/audits/Figure_1_plotted_points.csv')
