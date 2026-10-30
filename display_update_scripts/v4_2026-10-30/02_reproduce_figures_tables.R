#!/usr/bin/env Rscript

required_packages(c("data.table", "ggplot2", "patchwork", "svglite"))
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })

theme_paper <- function() {
  theme_bw(base_size = 10,base_family="Helvetica") +
    theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey92"),
          strip.text = element_text(face = "bold"), legend.position = "top",
          plot.title = element_text(face = "bold"), plot.subtitle = element_text(size = 9))
}

# The local export helper preserves the v3 ggplot geometry and records data invariants.
harmonize_facet_axes <- function(plot) {
  old <- plot
  plot <- plot + labs(title=NULL,subtitle=NULL,caption=NULL)
  g <- ggplotGrob(plot)
  ax <- unique(g$layout$l[grepl('^axis-l',g$layout$name)])
  if(length(ax)>1) {
    widest <- do.call(grid::unit.pmax,lapply(ax,function(j)g$widths[j]))
    for(j in ax)g$widths[j] <- widest
  }
  attr(g,'baseline_plot') <- old
  attr(g,'new_plot') <- plot
  g
}
source(file.path(ROOT,"display_update_scripts/v4_2026-10-30","02a_redraw_overall_inventory.R"),local=FALSE)

f4 <- read_csv(file.path(SOURCE_DIR,"Figure_4_estimates.csv"))
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

f5 <- read_csv(file.path(SOURCE_DIR,"Figure_5_estimates.csv"))[component == "upwind_minus_downwind"]
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

h21 <- read_csv(file.path(SOURCE_DIR,"Supplementary_Figure_S8_GFED_estimates.csv"))
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

source(file.path(ROOT,"display_update_scripts/v4_2026-10-30","02b_redraw_supplementary_diagnostics.R"),local=FALSE)

