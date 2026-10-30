# Display-only extraction of the supplied S7 drawing code.
FUELS <- c('coal','gas','oil','biomass'); METRICS <- c('count','heat'); RADII <- seq(20,50,5)
FUEL_LABEL <- c(coal='Coal',gas='Natural gas',oil='Oil',biomass='Biomass')
METRIC_LABEL <- c(count='Active plant count',heat='Actual daily Heat Input')
CMP <- fread('source_data_manuscript_v5/Supplementary_Figure_S7_natural_spline_effect_change.csv')
RESULTS <- rbind(CMP[,.(fuel,metric,radius_km,model='current',estimate=current_estimate,ci95_low=current_ci95_low,ci95_high=current_ci95_high)],
 CMP[,.(fuel,metric,radius_km,model='natural_spline',estimate=ns_estimate,ci95_low=ns_ci95_low,ci95_high=ns_ci95_high)])
PLOT <- copy(
  RESULTS
)

PLOT[
  ,
  fuel_label :=
    factor(
      FUEL_LABEL[
        fuel
      ],
      levels =
        unname(
          FUEL_LABEL[
            FUELS
          ]
        )
    )
]

PLOT[
  ,
  metric_label :=
    factor(
      METRIC_LABEL[
        metric
      ],
      levels =
        unname(
          METRIC_LABEL[
            METRICS
          ]
        )
    )
]

PLOT[
  ,
  model_label :=
    factor(
      fifelse(
        model ==
          "current",
        "Current main model",
        "Natural spline sensitivity"
      ),
      levels = c(
        "Current main model",
        "Natural spline sensitivity"
      )
    )
]

pd <- position_dodge(
  width = 1.4
)

p_all <- ggplot(
  PLOT,
  aes(
    x =
      radius_km,
    y =
      estimate,
    group =
      model_label,
    linetype =
      model_label,
    shape =
      model_label
  )
) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35
  ) +

  geom_errorbar(
    aes(
      ymin =
        ci95_low,
      ymax =
        ci95_high
    ),
    position =
      pd,
    width =
      1.2,
    linewidth =
      0.4
  ) +

  geom_line(
    position =
      pd,
    linewidth =
      0.65
  ) +

  geom_point(
    position =
      pd,
    size =
      1.8
  ) +

  facet_grid(
    metric_label ~
      fuel_label,
    scales =
      "free_y"
  ) +

  scale_x_continuous(
    breaks =
      RADII
  ) +

  labs(
    title =
      "Fuel-specific associations: current vs natural-spline covariate models",

    subtitle =
      paste0(
        "20–50 km; Monitor × year FE; ",
        "natural cubic splines with df=3 for continuous background covariates"
      ),

    x =
      "Cumulative radius (km)",

    y =
      "Adjusted gross-beta change (mBq/m³)",

    linetype =
      NULL,

    shape =
      NULL
  ) +

  theme_bw(
    base_size = 10
  ) +

  theme(
    panel.grid.minor =
      element_blank(),

    strip.text =
      element_text(
        face = "bold"
      ),

    legend.position =
      "bottom"
  )


baseline <- p_all
p_all <- p_all+labs(title=NULL,subtitle=NULL,caption=NULL)
export_figure(p_all,'Supplementary_Figure_S7',13,7.5,baseline=baseline)
