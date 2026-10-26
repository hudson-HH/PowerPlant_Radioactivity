#!/usr/bin/env Rscript

required_packages(c("data.table", "ggplot2", "patchwork", "svglite"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

theme_paper <- function() {
  theme_bw(base_size = 10,base_family="Helvetica") +
    theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
          strip.text = element_text(face = "bold"), legend.position = "top",
          plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 9))
}

save_both <- function(plot, stem, width, height) {
  p_png <- file.path(OUTPUT_DIR, "figures", paste0(stem, ".png"))
  p_pdf <- file.path(OUTPUT_DIR, "figures", paste0(stem, ".pdf"))
  ggsave(p_png, plot, width = width, height = height, dpi = 300, bg = "white")
  ggsave(p_pdf, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  svglite::svglite(file.path(OUTPUT_DIR,"figures",paste0(stem,".svg")),width=width,height=height,bg="white")
  if(inherits(plot,"grob")){grid::grid.newpage();grid::grid.draw(plot)}else print(plot)
  grDevices::dev.off()
  c(p_png, p_pdf)
}

harmonize_facet_axes <- function(plot) {
  g<-ggplotGrob(plot)
  ax<-unique(g$layout$l[grepl("^axis-l",g$layout$name)])
  if(length(ax)>1){
    widest<-do.call(grid::unit.pmax,lapply(ax,function(j)g$widths[j]))
    for(j in ax)g$widths[j]<-widest
  }
  g
}

source(file.path(ROOT,"R","02a_redraw_overall_inventory.R"),local=FALSE)

f4 <- read_csv(file.path(SOURCE_DIR, "main_figures", "Figure4.csv"))
f4[, `:=`(fuel_label = factor(fuel_labels[fuel], levels = unname(fuel_labels)),
          metric_label = factor(metric_labels[metric], levels = unname(metric_labels)))]
p4 <- ggplot(f4, aes(radius_km, estimate)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35) +
  geom_ribbon(aes(ymin = ci95_low, ymax = ci95_high), fill = "#4472C4", alpha = 0.16) +
  geom_line(color = "#234A84", linewidth = 0.6) + geom_point(color = "#234A84", size = 1.5) +
  facet_wrap(vars(metric_label,fuel_label),ncol=4, scales = "free_y") + scale_x_continuous(breaks = seq(20, 50, 10)) +
  labs(title = "Fuel-specific upwind actual-operation associations",
       subtitle = "Monitor × year fixed effects; effects evaluated at each model's sample-median wind speed",
       x = "Cumulative radius (km)", y = "Change in gross-beta radioactivity (mBq/m³)") + theme_paper()
p4<-harmonize_facet_axes(p4)
save_both(p4, "Figure4_upwind_actual_operation_monitorYearFE", 11, 6.7)

f5 <- read_csv(file.path(SOURCE_DIR, "main_figures", "Figure5.csv"))[component == "upwind_minus_downwind"]
f5[, `:=`(fuel_label = factor(fuel_labels[fuel], levels = unname(fuel_labels)),
          metric_label = factor(metric_labels[metric], levels = unname(metric_labels)))]
p5 <- ggplot(f5, aes(radius_km, estimate)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35) +
  geom_ribbon(aes(ymin = ci95_low, ymax = ci95_high), fill = "#C55A11", alpha = 0.16) +
  geom_line(color = "#8A3C08", linewidth = 0.6) + geom_point(color = "#8A3C08", size = 1.5) +
  facet_wrap(vars(metric_label,fuel_label),ncol=4, scales = "free_y") + scale_x_continuous(breaks = seq(20, 50, 10)) +
  labs(title = "Directional contrast: upwind minus downwind",
       subtitle = "Joint directional models with monitor × year fixed effects; sample-median wind speed",
       x = "Cumulative radius (km)", y = "Upwind − downwind estimate (mBq/m³)") + theme_paper()
p5<-harmonize_facet_axes(p5)
save_both(p5, "Figure5_upwind_minus_downwind_monitorYearFE", 11, 6.7)

plot_f6_panel <- function(file, title, model_levels) {
  d <- read_csv(file)
  d[, `:=`(fuel_label = factor(fuel_labels[fuel], levels = unname(fuel_labels)),
            metric_label = factor(metric_labels[metric], levels = unname(metric_labels)),
            model = factor(model, levels = model_levels))]
  ggplot(d, aes(fuel_label, estimate, ymin = ci_low, ymax = ci_high, shape=model, group=model)) +
    geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
    geom_errorbar(width = 0.15, linewidth = 0.45,position=position_dodge(width=0.5)) + geom_point(size = 1.6,position=position_dodge(width=0.5)) +
    facet_wrap(~metric_label, scales = "free_y",ncol=2) +
    labs(title = title, x = NULL, y = "Estimate (mBq/m³)") +
    theme_paper() + theme(axis.text.x = element_text(angle = 25, hjust = 1),legend.title=element_blank(),legend.position="bottom")
}

p6a <- plot_f6_panel(file.path(SOURCE_DIR, "main_figures", "Figure6A.csv"),
  "A  Common-time fixed effects",
  c("Monitor × year FE", "Monitor × year + common-time FE"))
p6b <- plot_f6_panel(file.path(SOURCE_DIR, "main_figures", "Figure6B.csv"),
  "B  Within- and between-monitor-year decomposition (baseline median wind)",
  c("Original random-intercept", "Within monitor-year", "Between monitor-years"))
p6c <- plot_f6_panel(file.path(SOURCE_DIR, "main_figures", "Figure6C.csv"),
  "C  PM₂.₅ adjustment",
  c("PM₂.₅ adjusted", "No PM₂.₅ adjustment"))
p6 <- p6a / p6b / p6c + plot_annotation(title = "Internal-validity and key robustness analyses")
save_both(p6, "Figure6_internal_validity_and_key_robustness", 11, 11)

h21 <- read_csv(file.path(SOURCE_DIR, "sensitivity", "H21_biomass_GFED_full.csv"))
assert_true(all(h21$wind_reference == "sample median wind"), "H21 source is not labelled sample-median wind.")
h21_long <- rbind(
  h21[, .(metric, radius_km, model = "Reference model", estimate = original_estimate,
           low = original_ci95_low, high = original_ci95_high)],
  h21[, .(metric, radius_km, model = "Adjusted for GFED fire-sourced PM2.5", estimate = fire_adjusted_estimate,
           low = fire_adjusted_ci95_low, high = fire_adjusted_ci95_high)]
)
h21_long[, metric_label := factor(metric_labels[metric], levels = unname(metric_labels))]
p21 <- ggplot(h21_long, aes(radius_km, estimate, group = model, linetype = model, shape = model)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.35) +
  geom_errorbar(aes(ymin = low, ymax = high), width = 0.7, position = position_dodge(width = 0.8)) +
  geom_line(position = position_dodge(width = 0.8)) +
  geom_point(size = 2, position = position_dodge(width = 0.8)) +
  facet_wrap(~metric_label, ncol = 1, scales = "free_y") + scale_x_continuous(breaks = seq(20, 50, 5)) +
  labs(title = "Biomass associations before and after adjustment for GFED fire-sourced PM2.5",
       subtitle = "2001–2023 matched samples; monitor × year fixed effects; effects evaluated at sample-median wind speed",
       x = "Cumulative radius (km)", y = "Adjusted gross-beta change (mBq/m³)", linetype = NULL, shape = NULL) + theme_paper()
p21<-harmonize_facet_axes(p21)
save_both(p21, "Figure_GFED_biomass_sensitivity", 10, 8.5)

source(file.path(ROOT,"R","02b_redraw_supplementary_diagnostics.R"),local=FALSE)

# Tables are copied from immutable source data into a clean output namespace.
table_map <- c(
  Table1 = file.path(SOURCE_DIR, "main_tables", "Table1.csv"),
  Table_environmental_covariates = file.path(SOURCE_DIR, "supplementary", "Table_environmental_covariates.csv"),
  Supplementary_Table_S2 = file.path(SOURCE_DIR, "supplementary", "Supplementary_Table_S2_Figure4_source_data.csv"),
  Supplementary_Table_S3 = file.path(SOURCE_DIR, "supplementary", "Supplementary_Table_S3_temporal_stability_20km.csv"),
  Supplementary_Table_S4 = file.path(SOURCE_DIR, "supplementary", "Supplementary_Table_S4_direction_sector_and_downwind_control.csv"),
  H20B_natural_spline_effect_change = file.path(SOURCE_DIR, "sensitivity", "04_current_vs_natural_spline_effect_change.csv"),
  H21_biomass_GFED_comparison = file.path(SOURCE_DIR, "sensitivity", "H21_biomass_GFED_comparison.csv")
)
for (nm in names(table_map)) {
  # Publication headings contain spaces; do not normalize them to R names.
  tab<-data.table::fread(table_map[[nm]],showProgress=FALSE,check.names=FALSE)
  if(nm=="Table1"){
    if("Primary analysis" %in% names(tab))
      setnames(tab,"Primary analysis","Overall inventory analysis")
    assert_columns(tab,c("Variable","Overall inventory analysis"),"Table 1")
    tab[grepl("^Soil moisture",Variable),Variable:="Soil moisture (kg/m²)"]
  }
  fwrite(tab,file.path(OUTPUT_DIR,"tables",paste0(nm,".csv")))
}

cat("Figure/table redraw: PASS\n")
