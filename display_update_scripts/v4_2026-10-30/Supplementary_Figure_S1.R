# Descriptive geographic grid: every plant appears once in each combustion role
# it has, using the exact non-exclusive Figure 1 display table and coordinates.
s1_roles <- c('Coal','Natural gas','Oil','Other fossil','Biomass','Solid waste')
PLANT_S1 <- fread('data/display_assets/v4_2026-10-30/source_data/Figure_1_combustion_locations.csv')
stopifnot(setequal(unique(PLANT_S1$fuel_label),s1_roles),
          !anyDuplicated(PLANT_S1[,.(plant_id_plot,fuel_label)]),
          all(is.finite(PLANT_S1$plant_longitude)),
          all(is.finite(PLANT_S1$plant_latitude)),
          all(PLANT_S1$plant_longitude >= -125 & PLANT_S1$plant_longitude <= -66),
          all(PLANT_S1$plant_latitude >= 24 & PLANT_S1$plant_latitude <= 50))
PLANT_S1[,fuel_label:=factor(fuel_label,levels=s1_roles)]
US_MAP <- ggplot2::map_data('state')

pS1 <- ggplot() +
  geom_polygon(
    data = US_MAP,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "white",
    color = "grey82",
    linewidth = 0.2
  ) +
  geom_point(
    data = PLANT_S1,
    aes(
      x = plant_longitude,
      y = plant_latitude
    ),
    size = 0.45,
    alpha = 0.55,
    colour = 'black'
  ) +
  facet_wrap(
    ~fuel_label,
    ncol = 3,
    drop = FALSE
  ) +
  coord_quickmap(
    xlim = c(
      -125,
      -66
    ),
    ylim = c(
      24,
      50
    ),
    expand = FALSE
  ) +
  labs(
    title = "Figure S1. Spatial distribution of individual power-plant energy / fuel types"
  ) +
  theme_void(
    base_size = 9
  ) +
  theme(
    strip.text = element_text(
      face = "bold"
    ),
    plot.margin = margin(5.76,5.76,5.76,5.76,unit='pt'),
    plot.title = element_text(
      face = "bold"
    )
  )


baseline <- pS1
pS1 <- pS1+labs(title=NULL,subtitle=NULL,caption=NULL)
fwrite(PLANT_S1,'data/display_assets/v4_2026-10-30/source_data/Supplementary_Figure_S1_combustion_locations.csv')
s1counts <- PLANT_S1[,.(unique_plants=uniqueN(plant_id_plot),plant_role_points=.N),by=fuel_label]
s1counts <- s1counts[order(fuel_label)]
stopifnot(identical(as.integer(s1counts$unique_plants),c(798L,3500L,2675L,108L,1057L,188L)))
fwrite(s1counts,'data/display_assets/v4_2026-10-30/source_data/Supplementary_Figure_S1_counts.csv')
jsonlite::write_json(setNames(as.list(s1counts$unique_plants),as.character(s1counts$fuel_label)),
 'data/display_assets/v4_2026-10-30/source_data/Supplementary_Figure_S1_counts.json',auto_unbox=TRUE,pretty=TRUE)
# Fixed-aspect map panels need a content-derived canvas height. Measure the
# gtable on the final PDF device so strip typography is included, then add only
# the explicit uniform plot margins; no post-export crop or coordinate change.
s1_width <- 13
grDevices::cairo_pdf('data/display_assets/v4_2026-10-30/.tmp/Supplementary_Figure_S1_sizing.pdf',width=s1_width,height=7,family='Helvetica')
s1_grob <- ggplotGrob(pS1)
s1_build <- ggplot_build(pS1)
s1_aspect <- s1_build$layout$coord$aspect(s1_build$layout$panel_params[[1]])
s1_fixed_width <- convertWidth(sum(s1_grob$widths),'in',valueOnly=TRUE)
s1_fixed_height <- convertHeight(sum(s1_grob$heights),'in',valueOnly=TRUE)
s1_height <- s1_fixed_height + 2 * (s1_width-s1_fixed_width)/3 * s1_aspect
dev.off()
# Cairo uses whole PDF points; round upward so vector and raster canvases agree.
s1_height <- ceiling(s1_height*72)/72
export_figure(pS1,'Supplementary_Figure_S1',s1_width,s1_height,baseline=baseline)
