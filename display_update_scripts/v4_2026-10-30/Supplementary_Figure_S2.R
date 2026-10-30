# Final publication section of the supplied regional map script.
suppressPackageStartupMessages(library(sf))
P_PUB <- unique(fread('source_data_manuscript_v5/Figure_1A_powerplant_locations.csv')[,.(plant_id_plot,plant_longitude,plant_latitude)])
M_PUB <- fread('source_data_manuscript_v5/Figure_1_radnet_monitor_locations.csv')
eligible <- unique(fread('source_data_manuscript_v5/Supplementary_Figure_S5_fuel_specific_LOO.csv')$omitted_monitor)
stopifnot(length(eligible)==133)
M_PUB <- M_PUB[loc_num %in% eligible]
stopifnot(nrow(M_PUB)==133)
assign_region_pub <- function(lon) {

  fifelse(
    lon < -105,
    "A. Western CONUS",
    fifelse(
      lon < -90,
      "B. Central CONUS",
      "C. Eastern CONUS"
    )
  )
}

M_PUB[
  ,
  region :=
    assign_region_pub(
      monitor_longitude
    )
]


P_SF_PUB <- st_as_sf(
  P_PUB,
  coords = c(
    "plant_longitude",
    "plant_latitude"
  ),
  crs = 4326,
  remove = FALSE
)

M_SF_PUB <- st_as_sf(
  M_PUB,
  coords = c(
    "monitor_longitude",
    "monitor_latitude"
  ),
  crs = 4326,
  remove = FALSE
)

P_5070_PUB <- st_transform(
  P_SF_PUB,
  5070
)

M_5070_PUB <- st_transform(
  M_SF_PUB,
  5070
)

B20_PUB <- st_transform(
  st_buffer(
    M_5070_PUB,
    dist = 20000
  ),
  4326
)

B50_PUB <- st_transform(
  st_buffer(
    M_5070_PUB,
    dist = 50000
  ),
  4326
)

US_RAW_PUB <- maps::map(
  "state",
  plot = FALSE,
  fill = TRUE
)

US_SF_PUB <- st_as_sf(
  US_RAW_PUB
)

if (
  is.na(
    st_crs(
      US_SF_PUB
    )
  )
) {

  US_SF_PUB <- st_set_crs(
    US_SF_PUB,
    4326
  )
}


REGIONS_PUB <- c(
  "A. Western CONUS",
  "B. Central CONUS",
  "C. Eastern CONUS"
)

REGIONAL_OBJECTS_PUB <- list()

REGIONAL_COUNTS_PUB <- data.table()

for (
  reg in REGIONS_PUB
) {

  MM_5070 <- M_5070_PUB[
    M_5070_PUB$region ==
      reg,
  ]

  MM <- M_SF_PUB[
    M_SF_PUB$region ==
      reg,
  ]

  BB20 <- B20_PUB[
    B20_PUB$region ==
      reg,
  ]

  BB50 <- B50_PUB[
    B50_PUB$region ==
      reg,
  ]

  if (
    nrow(
      MM
    ) ==
      0L
  ) {

    stop(
      paste0(
        "区域无monitor：",
        reg
      ),
      call. = FALSE
    )
  }

  PLANT_MATCH <- lengths(
    st_is_within_distance(
      P_5070_PUB,
      MM_5070,
      dist = 50000
    )
  ) > 0L

  PP <- P_SF_PUB[
    PLANT_MATCH,
  ]

  REGIONAL_COUNTS_PUB <- rbind(
    REGIONAL_COUNTS_PUB,

    data.table(
      region =
        reg,

      n_monitors =
        nrow(
          MM
        ),

      n_unique_plants_within_50km =
        nrow(
          PP
        )
    )
  )

  if (
    nrow(PP) > 0L
  ) {

    GEOM_FOR_BBOX <- c(
      st_geometry(MM),
      st_geometry(PP)
    )

  } else {

    GEOM_FOR_BBOX <- st_geometry(MM)
  }

  BB <- st_bbox(
    GEOM_FOR_BBOX
  )

  XSPAN <- as.numeric(
    BB["xmax"] -
      BB["xmin"]
  )

  YSPAN <- as.numeric(
    BB["ymax"] -
      BB["ymin"]
  )

  XPAD <- max(
    0.9,
    0.07 *
      XSPAN
  )

  YPAD <- max(
    0.7,
    0.07 *
      YSPAN
  )

  XLIM <- c(
    max(
      -125,
      as.numeric(
        BB["xmin"]
      ) -
        XPAD
    ),

    min(
      -66,
      as.numeric(
        BB["xmax"]
      ) +
        XPAD
    )
  )

  YLIM <- c(
    max(
      24,
      as.numeric(
        BB["ymin"]
      ) -
        YPAD
    ),

    min(
      50,
      as.numeric(
        BB["ymax"]
      ) +
        YPAD
    )
  )

  REGIONAL_OBJECTS_PUB[[reg]] <- list(
    monitors =
      MM,

    plants =
      PP,

    buffer20 =
      BB20,

    buffer50 =
      BB50,

    xlim =
      XLIM,

    ylim =
      YLIM
  )
}


make_S2_panel_pub <- function(reg) {

  OBJ <- REGIONAL_OBJECTS_PUB[[reg]]

  CC <- REGIONAL_COUNTS_PUB[
    region ==
      reg
  ]

  XLIM <- OBJ$xlim
  YLIM <- OBJ$ylim

  DX <- diff(
    XLIM
  )

  DY <- diff(
    YLIM
  )

  LABEL_DATA <- data.frame(
    x =
      XLIM[1] +
      0.025 *
      DX,

    y =
      YLIM[2] -
      0.025 *
      DY,

    label =
      paste0(
        CC$n_monitors,
        " RadNet monitors\n",
        format(
          CC$n_unique_plants_within_50km,
          big.mark = ","
        ),
        " unique power plants within 50 km"
      )
  )

  ggplot() +

    geom_sf(
      data = US_SF_PUB,
      fill = "white",
      color = "grey82",
      linewidth = 0.25
    ) +

    geom_sf(
      data = OBJ$buffer50,
      fill = NA,
      color = "grey58",
      linewidth = 0.38,
      linetype = "dashed"
    ) +

    geom_sf(
      data = OBJ$buffer20,
      fill = NA,
      color = "grey20",
      linewidth = 0.48
    ) +

    geom_sf(
      data = OBJ$plants,
      color = "grey48",
      size = 0.48,
      alpha = 0.34
    ) +

    geom_sf(
      data = OBJ$monitors,
      shape = 17,
      color = "black",
      size = 1.7
    ) +

    geom_label(
      data =
        LABEL_DATA,

      aes(
        x = x,
        y = y,
        label = label
      ),

      inherit.aes = FALSE,

      hjust = 0,
      vjust = 1,

      size = 2.65,

      lineheight = 1.03,

      fill = "white",

      color = "black",

      label.size = 0.15,

      label.padding =
        unit(
          0.11,
          "lines"
        )
    ) +

    coord_sf(
      xlim = XLIM,
      ylim = YLIM,
      expand = FALSE
    ) +

    labs(
      title =
        reg
    ) +

    theme_void(
      base_size = 10
    ) +

    theme(
      plot.title =
        element_text(
          face = "bold",
          size = 11,
          hjust = 0
        ),

      panel.border =
        element_rect(
          color = "grey75",
          fill = NA,
          linewidth = 0.45
        ),

      plot.margin =
        margin(
          5,
          6,
          5,
          6
        )
    )
}


S2_A_PUB <- make_S2_panel_pub(
  "A. Western CONUS"
)

S2_B_PUB <- make_S2_panel_pub(
  "B. Central CONUS"
)

S2_C_PUB <- make_S2_panel_pub(
  "C. Eastern CONUS"
)

S2_PUBLICATION <- (
  S2_A_PUB |
    S2_B_PUB |
    S2_C_PUB
) +

  plot_layout(
    widths = c(
      1,
      1,
      1
    )
  ) +

  plot_annotation(
    title =
      "Supplementary Figure S2. Regional distribution of RadNet monitors and nearby power plants",

    subtitle =
      paste0(
        "Solid circles indicate 20-km monitor-centered buffers; dashed circles indicate 50-km buffers.\n",
        "Black triangles represent RadNet monitors; grey points represent unique nearby power plants. ",
        "Regional divisions are descriptive and are not statistical analysis strata."
      ),

    caption =
      paste0(
        "Regional plant counts represent unique power plants located within 50 km of at least one RadNet monitor shown in that panel. ",
        "A plant may contribute to more than one descriptive panel if it lies within 50 km of monitors on opposite sides of a regional boundary."
      ),

    theme =
      theme(
        plot.title =
          element_text(
            face = "bold",
            size = 14
          ),

        plot.subtitle =
          element_text(
            size = 9.5,
            lineheight = 1.08
          ),

        plot.caption =
          element_text(
            size = 7.5,
            hjust = 0
          )
      )
  )



reference_counts <- fread('source_data_manuscript_v5/Supplementary_Figure_S2_regional_counts.csv')
stopifnot(isTRUE(all.equal(as.data.frame(REGIONAL_COUNTS_PUB),as.data.frame(reference_counts),check.attributes=FALSE)))
fwrite(REGIONAL_COUNTS_PUB,'data/display_assets/v4_2026-10-30/qa/audits/S2_regional_count_recheck.csv')
baseline <- S2_PUBLICATION
remove_sample_label <- function(p) {
 p$layers <- Filter(function(l)!inherits(l$geom,'GeomLabel'),p$layers)
 p
}
S2_A_clean <- remove_sample_label(S2_A_PUB)
S2_B_clean <- remove_sample_label(S2_B_PUB)
S2_C_clean <- remove_sample_label(S2_C_PUB)
clean <- (S2_A_clean|S2_B_clean|S2_C_clean)+plot_layout(widths=c(1,1,1))
# Compare every scientific layer with the original after excluding the expressly removed labels.
scientific_baseline <- (remove_sample_label(S2_A_PUB)|remove_sample_label(S2_B_PUB)|remove_sample_label(S2_C_PUB))+plot_layout(widths=c(1,1,1))
export_figure(clean,'Supplementary_Figure_S2',18,7.2,baseline=scientific_baseline)
saveRDS(list(before=baseline,after=clean),'data/display_assets/v4_2026-10-30/qa/audits/S2_plot_objects.rds')
