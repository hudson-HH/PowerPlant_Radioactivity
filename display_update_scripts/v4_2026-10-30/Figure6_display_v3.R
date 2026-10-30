#!/usr/bin/env Rscript
# Portable final-display module. All drawing/export/previews are R-only.
# Loaded by display_update_scripts/v4_2026-10-30/run_all.R with local source-data and relative export paths.
.figure6_module_file <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) {
  f <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(f)) normalizePath(sub("^--file=", "", f[1]), mustWork = TRUE) else NA_character_
})
.figure6_font_audit <- file.path(dirname(.figure6_module_file), "audit_effective_pdf_fonts.py")

render_figure6_v3 <- function(data, output_stem, qa_dir, work_dir,
  skill_root = Sys.getenv("NATURE_FIGURE_SKILL", file.path(dirname(.figure6_module_file), "..", "tools", "nature-figure")),
  qa_python = Sys.getenv("NATURE_FIGURE_QA_PYTHON", Sys.which("python3")),
  font_audit_script = .figure6_font_audit) {
  needed <- c("data.table", "ggplot2", "patchwork", "grid", "jsonlite", "ragg", "svglite")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R plotting packages: ", paste(missing, collapse = ", "))
  suppressPackageStartupMessages({library(data.table); library(ggplot2); library(patchwork); library(grid)})
  results <- data.table::copy(data.table::as.data.table(data))
  stopifnot(all(c("panel", "fuel", "metric", "model", "estimate", "ci_low", "ci_high", "wind_ms") %in% names(results)))
  stopifnot(nrow(results) == 56L, all(is.finite(results$estimate)),
    all(results$ci_low <= results$estimate), all(results$ci_high >= results$estimate),
    data.table::uniqueN(results, by = c("panel", "fuel", "metric", "model")) == 56L)
  reference_winds <- c(coal_count = 3.15442053456934, gas_count = 3.15442053456934,
    oil_count = 3.15442053456934, biomass_count = 3.15442053456934,
    coal_heat = 3.15275005725255, gas_heat = 3.13870352480861,
    oil_heat = 3.13587557596978, biomass_heat = 3.15442053456934)
  stopifnot(all(abs(results$wind_ms - reference_winds[paste(results$fuel, results$metric, sep = "_")]) < 1e-12))
  logs <- qa_dir; work <- work_dir
  for (p in c(dirname(output_stem), logs, work)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

  fuel_labels <- c(coal = "Coal", gas = "Natural gas", oil = "Oil", biomass = "Biomass")
  metric_labels <- c(count = "CAMPD plant count\n(per additional plant)",
    heat = "Reported heat input\n(per 10,000 MMBtu/day)")
  models <- list(A = c("Monitor × year FE", "Monitor × year + common-time FE"),
    B = c("Original random-intercept", "Within monitor-year", "Between monitor-years"),
    C = c("PM₂.₅ adjusted", "No PM₂.₅ adjustment"))
  display <- c("Monitor × year FE" = "Monitor × year FE", "Monitor × year + common-time FE" = "Monitor × year + common-time FE",
    "Original random-intercept" = "Original random-intercept", "Within monitor-year" = "Within monitor-year", "Between monitor-years" = "Between monitor-years",
    "PM₂.₅ adjusted" = "PM2.5 adjusted", "No PM₂.₅ adjustment" = "No PM2.5 adjustment")
  palettes <- list(A = c("#2C628E", "#B87529"), B = c("#595959", "#2C628E", "#B87529"), C = c("#2C628E", "#B87529"))
  shapes <- list(A = c(16, 17), B = c(16, 17, 15), C = c(16, 17))
  theme_fig <- theme_classic(base_size = 6.5, base_family = "Helvetica") + theme(
    axis.line = element_line(linewidth = 0.28, colour = "#333333"),
    axis.ticks = element_line(linewidth = 0.25, colour = "#333333"),
    axis.text = element_text(size = 6.2, colour = "#222222"),
    axis.title = element_text(size = 6.5),
    axis.title.y = element_text(margin = margin(r = 5)),
    strip.text = element_text(size = 6.5, face = "bold", margin = margin(t = 2, b = 4)),
    strip.background = element_rect(fill = "#F3F3F3", colour = NA),
    panel.grid = element_blank(), panel.spacing = unit(8, "mm"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0, margin = margin(b = 4)),
    legend.position = "bottom", legend.title = element_blank(),
    legend.text = element_text(size = 6.2), legend.key.width = unit(3.5, "mm"),
    legend.key.height = unit(3, "mm"), legend.spacing.x = unit(1.5, "mm"),
    legend.margin = margin(t = 2, b = 2), legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(3, 4, 3, 4))
  make_panel <- function(p, title) {
    z <- copy(results[panel == p])
    z[, `:=`(fuel_label = factor(fuel_labels[fuel], levels = unname(fuel_labels)),
      metric_label = factor(metric_labels[metric], levels = unname(metric_labels)),
      model = factor(model, levels = models[[p]]))]
    pd <- position_dodge(width = 0.58)
    ggplot(z, aes(fuel_label, estimate, ymin = ci_low, ymax = ci_high,
      colour = model, shape = model, group = model)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "#777777", linewidth = 0.25) +
      geom_errorbar(width = 0.12, linewidth = 0.38, position = pd) +
      geom_point(size = 1.3, position = pd) +
      scale_colour_manual(values = setNames(palettes[[p]], models[[p]]), labels = display[models[[p]]], drop = FALSE) +
      scale_shape_manual(values = setNames(shapes[[p]], models[[p]]), labels = display[models[[p]]], drop = FALSE) +
      scale_y_continuous(n.breaks = 4, expand = expansion(mult = c(0.07, 0.09))) +
      facet_wrap(vars(metric_label), ncol = 2, scales = "free_y") +
      labs(title = title, x = NULL, y = "Gross-beta estimate (mBq/m³)") + theme_fig +
      guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1))
  }
  pA <- make_panel("A", "a  Common-time fixed effects")
  pB <- make_panel("B", "b  Within- and between-monitor-year decomposition")
  pC <- make_panel("C", "c  PM2.5 adjustment")
  baseline <- pA / pB / pC + plot_layout(heights = c(1, 1, 1)) + plot_annotation(
    title = "Internal validity and robustness",
    subtitle = "20 km upwind exposure; shared baseline sample-median wind; 95% confidence intervals",
    caption = paste0("Count: 3.154 m/s (all fuels). Heat: coal 3.153, gas 3.139, oil 3.136, biomass 3.154 m/s.\n",
      "Reference medians are shared across a-c; model-specific sample sizes are provided in Source Data."),
    theme = theme(plot.title = element_text(size = 9, face = "bold", family = "Helvetica", hjust = 0),
      plot.subtitle = element_text(size = 6.2, family = "Helvetica", hjust = 0, margin = margin(b = 5)),
      plot.caption = element_text(size = 6, family = "Helvetica", hjust = 0, margin = margin(t = 4)),
      plot.margin = margin(5, 5, 5, 5)))
  fig <- pA/pB/pC+plot_layout(heights=c(1,1,1))
  export_figure(fig,'Figure_6',183/25.4,200/25.4,baseline=baseline)
  invisible(list(figure=fig,panels=list(pA,pB,pC),source_data=results))
}
