# Display-only extraction of the supplied S6 drawing code.
COMB_FUELS <- c('coal','gas','oil','biomass')
COMB_LABEL <- c(coal='Coal',gas='Natural gas',oil='Oil',biomass='Biomass')
METRICS <- c('count','heat')
METRIC_LABEL <- c(count='Active plant Count',heat='Actual daily Heat Input')
FE_LABEL <- c(monitor_FE='Monitor FE',monitor_year_FE='Monitor × year FE')
MODEL_LABEL <- c(M0_noncomb_matched='Unadjusted · matched sample',M1_noncomb_adjusted='Adjusted for non-combustion energy')
RADII <- seq(20,50,5)
MED <- fread('source_data_manuscript_v5/Supplementary_Figure_S6_adjusted_combustion_effects.csv')[wind_reference=='sample_median']
decorate <- function(D) {

  D[
    ,
    fuel_label :=
      factor(
        COMB_LABEL[
          fuel
        ],
        levels =
          unname(
            COMB_LABEL[
              COMB_FUELS
            ]
          )
      )
  ]


  D[
    ,
    fe_label :=
      factor(
        FE_LABEL[
          fe_structure
        ],
        levels =
          unname(
            FE_LABEL[
              c(
                "monitor_FE",
                "monitor_year_FE"
              )
            ]
          )
      )
  ]


  D[
    ,
    model_label :=
      factor(
        MODEL_LABEL[
          model_id
        ],
        levels =
          unname(
            MODEL_LABEL[
              c(
                "M0_noncomb_matched",
                "M1_noncomb_adjusted"
              )
            ]
          )
      )
  ]


  D
}


theme_my <- theme_bw(
  base_size = 10
) +
  theme(
    panel.grid.minor =
      element_blank(),

    strip.background =
      element_rect(
        fill = "white"
      ),

    strip.text =
      element_text(
        face = "bold"
      ),

    legend.position =
      "bottom"
  )


FINAL <- decorate(
  copy(
    MED[
      model_id ==
        "M1_noncomb_adjusted"
    ]
  )
)


FINAL[
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


p4 <- ggplot(
  FINAL,
  aes(
    x = radius_km,
    y = estimate,
    group = fuel_label,
    color = fuel_label,
    linetype = fuel_label,
    shape = fuel_label
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  geom_errorbar(
    aes(
      ymin = ci95_low,
      ymax = ci95_high
    ),
    width = 0.7,
    linewidth = 0.35
  ) +
  geom_line(
    linewidth = 0.7
  ) +
  geom_point(
    size = 1.8
  ) +
  facet_grid(
    fe_label ~ metric_label,
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = RADII
  ) +
  labs(
    title =
      "Combustion associations after adjustment for non-combustion energy",

    subtitle =
      "Wind, solar, hydro, nuclear and geothermal upwind Count simultaneously adjusted.",

    x =
      "Cumulative radius (km)",

    y =
      "Adjusted gross-beta effect",

    color =
      "Fuel",

    linetype =
      "Fuel",

    shape =
      "Fuel"
  ) +
  theme_my



baseline <- p4
p4 <- p4+labs(title=NULL,subtitle=NULL,caption=NULL)
export_figure(p4,'Supplementary_Figure_S6',13,8.5,baseline=baseline)
