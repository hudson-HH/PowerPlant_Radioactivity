#!/usr/bin/env Rscript
# Portable final-display module. All drawing/export/previews are R-only.
# CLI: Rscript Figure6_display_v3.R input.csv output_stem qa_dir work_dir
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
  source(file.path(skill_root, "scripts", "panel_alignment.R"))
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
  fig <- pA / pB / pC + plot_layout(heights = c(1, 1, 1)) + plot_annotation(
    title = "Internal validity and robustness",
    subtitle = "20 km upwind exposure; shared baseline sample-median wind; 95% confidence intervals",
    caption = paste0("Count: 3.154 m/s (all fuels). Heat: coal 3.153, gas 3.139, oil 3.136, biomass 3.154 m/s.\n",
      "Reference medians are shared across a-c; model-specific sample sizes are provided in Source Data."),
    theme = theme(plot.title = element_text(size = 9, face = "bold", family = "Helvetica", hjust = 0),
      plot.subtitle = element_text(size = 6.2, family = "Helvetica", hjust = 0, margin = margin(b = 5)),
      plot.caption = element_text(size = 6, family = "Helvetica", hjust = 0, margin = margin(t = 4)),
      plot.margin = margin(5, 5, 5, 5)))
  width_mm = 183; height_mm = 200
  width_in <- width_mm / 25.4; height_in <- height_mm / 25.4
  stem <- output_stem
  # The generic helper is loaded above. Facet/null units are instead measured
  # from fully drawn grid viewports on the same device and physical dimensions
  # as the final PDF; unit-sum placeholders are never accepted as panel areas.
  probe <- file.path(work, "Figure6_geometry_probe.pdf")
  grDevices::pdf(probe, width = width_in, height = height_in, family = "Helvetica", useDingbats = FALSE)
  grob <- patchworkGrob(fig)
  grid.newpage(); grid.draw(grob); grid.force()
  vps <- grid.ls(viewports = TRUE, grobs = FALSE, print = FALSE)
  hits <- which(vps$type == "vpListing" & grepl("^panel(?:-[0-9]+)*[.]", vps$name, perl = TRUE))
  # One viewport per actual data facet, not each outer patchwork cell.
  boxes <- lapply(hits, function(j) {
    path <- c(strsplit(vps$vpPath[j], "::", fixed = TRUE)[[1]], vps$name[j])
    path <- path[path != "ROOT" & nzchar(path)]
    upViewport(0); downViewport(do.call(vpPath, as.list(path)))
    a <- deviceLoc(unit(0, "npc"), unit(0, "npc"), valueOnly = TRUE)
    b <- deviceLoc(unit(1, "npc"), unit(1, "npc"), valueOnly = TRUE)
    c(a$x, a$y, b$x, b$y) * 72
  })
  upViewport(0); dev.off()
  rectangles <- unique(do.call(rbind, boxes))
  stopifnot(nrow(rectangles) == 6L, all(rectangles[, 3] > rectangles[, 1]), all(rectangles[, 4] > rectangles[, 2]))
  rectangles <- rectangles[order(-round(rectangles[, 4], 4), rectangles[, 1]), , drop = FALSE]
  ids <- c("a_count", "a_heat", "b_count", "b_heat", "c_count", "c_heat")
  manifest <- list(schema_version = 1L, backend = "r-forced-PDF-viewport",
    figure = list(width_pt = width_in * 72, height_pt = height_in * 72),
    panels = lapply(1:6, function(i) list(id = ids[i], bbox_pt = rectangles[i, ])),
    row_groups = list(ids[1:2], ids[3:4], ids[5:6]),
    column_groups = list(ids[c(1, 3, 5)], ids[c(2, 4, 6)]), exemptions = list())
  layout_file <- file.path(logs, "Figure6.alignment-layout.json")
  jsonlite::write_json(manifest, layout_file, auto_unbox = TRUE, pretty = TRUE, digits = NA)
  audit_args <- c(shQuote(file.path(skill_root, "scripts", "audit_panel_alignment.py")), shQuote(layout_file),
    "--json-out", shQuote(file.path(logs, "Figure6.alignment.json")),
    "--overlay-svg", shQuote(file.path(logs, "Figure6.alignment.svg")),
    "--tolerance-pt", "1.5", "--gutter-tolerance-pt", "1.5", "--strict")
  alignment_status <- system2(qa_python, audit_args)
  if (alignment_status != 0L) stop("Figure6 physical panel alignment failed")
  # Native PDF preserves editable Helvetica text and simple independent traces.
  grDevices::pdf(paste0(stem, ".pdf"), width = width_in, height = height_in, family = "Helvetica", useDingbats = FALSE)
  print(fig); dev.off()
  svglite::svglite(paste0(stem, ".svg"), width = width_in, height = height_in)
  print(fig); dev.off()
  ragg::agg_png(paste0(stem, ".png"), width = width_mm, height = height_mm, units = "mm", res = 600, background = "white")
  print(fig); dev.off()
  ragg::agg_tiff(paste0(stem, ".tiff"), width = width_mm, height = height_mm, units = "mm", res = 600, compression = "lzw", background = "white")
  print(fig); dev.off()
  for (i in seq_along(list(pA, pB, pC))) {
    ragg::agg_png(file.path(logs, paste0("Figure6_panel_", letters[i], "_preview.png")),
      width = width_mm, height = 64, units = "mm", res = 300, background = "white")
    print(list(pA, pB, pC)[[i]]); dev.off()
  }
  saveRDS(list(figure = fig, panels = list(pA, pB, pC), dimensions_mm = c(width_mm, height_mm)), file.path(work, "Figure6_plot_objects.rds"))
  collision_status <- system2(qa_python, c(shQuote(file.path(skill_root, "scripts", "audit_figure_collisions.py")),
    shQuote(paste0(stem, ".pdf")), "--json-out", shQuote(file.path(logs, "Figure6.collision-audit.json")), "--strict"))
  if (collision_status != 0L) stop("Figure6 final-PDF collision gate failed; inspect report")
  font_status <- system2(qa_python, c(shQuote(font_audit_script), shQuote(paste0(stem, ".pdf")),
    shQuote(paste0(stem, ".svg")), "--json-out", shQuote(file.path(logs, "Figure6.effective-font-audit.json")), "--min-pt", "5"))
  if (font_status != 0L) stop("Figure6 effective-font/editable-text gate failed")
  cat("Figure6 R export and rendered alignment/collision/effective-font gates completed.\n")
  invisible(list(figure = fig, panels = list(pA, pB, pC), manifest = manifest,
    source_data = results, width_mm = width_mm, height_mm = height_mm))
}
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 4L) stop("Usage: Rscript Figure6_display_v3.R input.csv output_stem qa_dir work_dir")
  render_figure6_v3(data.table::fread(args[1]), args[2], args[3], args[4])
}
