#!/usr/bin/env Rscript

options(
  stringsAsFactors = FALSE,
  warn = 1,
  scipen = 999,
  width = 220
)

required_packages <- c(
  "data.table",
  "ggplot2",
  "gridExtra",
  "grid"
)

for (p in required_packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(
      p,
      repos = "https://cloud.r-project.org"
    )
  }
}

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})


# ============================================================
# 07H18B
#
# Supplementary Table:
# Environmental covariates and gross-beta particle radioactivity
#
# 正式来源：
# 03_TableS_environmental_covariates_primary_monitorYearFE.csv
#
# 输出：
# 1. publication CSV
# 2. PNG
# 3. PDF
# ============================================================


# ============================================================
# 1. Paths
# ============================================================

ROOT <- "/Users/hudson/Desktop/USA 辐射表征研究"

HROOT <- file.path(
  ROOT,
  "发电站数据处理",
  "06_07H_燃料特异CAMPD逐日暴露与模型"
)

IN_DIR <- file.path(
  HROOT,
  "07H18_环境协变量与gross_beta关联表"
)

IN_FILE <- file.path(
  IN_DIR,
  "03_TableS_environmental_covariates_primary_monitorYearFE.csv"
)

AUDIT_FILE <- file.path(
  IN_DIR,
  "06_model_audit.csv"
)

HARD_FILE <- file.path(
  IN_DIR,
  "07_07H18_hard_checks.csv"
)

OUT_CSV <- file.path(
  IN_DIR,
  "08_TableS_environmental_covariates_publication.csv"
)

OUT_PNG <- file.path(
  IN_DIR,
  "08_TableS_environmental_covariates_publication.png"
)

OUT_PDF <- file.path(
  IN_DIR,
  "08_TableS_environmental_covariates_publication.pdf"
)


# ============================================================
# 2. Input audit
# ============================================================

REQ <- c(
  IN_FILE,
  AUDIT_FILE,
  HARD_FILE
)

MISS <- REQ[
  !file.exists(
    REQ
  )
]

if (
  length(MISS) > 0L
) {
  stop(
    paste0(
      "缺少输入文件：\n",
      paste(
        MISS,
        collapse = "\n"
      )
    ),
    call. = FALSE
  )
}

HARD <- fread(
  HARD_FILE,
  showProgress = FALSE
)

N_HARD_FAIL <- HARD[
  severity ==
    "HARD" &
  status ==
    "FAIL",
  .N
]

if (
  N_HARD_FAIL > 0L
) {
  stop(
    "07H18仍存在Hard FAIL，禁止制作论文表。",
    call. = FALSE
  )
}

AUDIT <- fread(
  AUDIT_FILE,
  showProgress = FALSE
)

D <- fread(
  IN_FILE,
  showProgress = FALSE
)

if (
  nrow(D) != 12L
) {
  stop(
    paste0(
      "正式环境协变量表应有12行，实际=",
      nrow(D)
    ),
    call. = FALSE
  )
}


# ============================================================
# 3. Helper
# ============================================================

fmt <- function(x) {

  ans <- rep(
    NA_character_,
    length(x)
  )

  for (
    i in seq_along(x)
  ) {

    z <- x[i]

    if (
      !is.finite(z)
    ) {

      ans[i] <- NA_character_

    } else if (
      z != 0 &&
      abs(z) < 0.001
    ) {

      ans[i] <- formatC(
        z,
        format = "e",
        digits = 2
      )

    } else if (
      abs(z) >= 100
    ) {

      ans[i] <- formatC(
        z,
        format = "f",
        digits = 2
      )

    } else {

      ans[i] <- formatC(
        z,
        format = "f",
        digits = 4
      )
    }
  }

  ans
}


# ============================================================
# 4. Publication labels
# ============================================================

LABELS <- data.table(
  Environmental_factor = c(
    "Continental air-mass proportion",
    "Monthly sunspot number",
    "PM2.5",
    "Inverse planetary boundary-layer height",
    "Relative humidity",
    "Soil moisture",
    "Wind speed",
    "Temperature",
    "Temperature squared",
    "Temperature cubed",
    "Sunspots × latitude",
    "Surface uranium × soil moisture"
  ),

  Display_factor = c(
    "Continental air-mass proportion",
    "Monthly sunspot number",
    "PM2.5",
    "Inverse PBLH",
    "Relative humidity",
    "Soil moisture",
    "Wind speed",
    "Temperature",
    "Temperature²",
    "Temperature³",
    "Sunspots × latitude",
    "Surface uranium × soil moisture"
  ),

  Interpretation = c(
    "Continental air-mass contribution",
    "Solar activity",
    "Fine particulate matter",
    "Atmospheric mixing and dilution",
    "Atmospheric moisture",
    "Surface soil-moisture condition",
    "Atmospheric transport",
    "Temperature polynomial term",
    "Temperature polynomial term",
    "Temperature polynomial term",
    "Modification of solar-activity association by latitude",
    "Modification of soil-moisture association by surface uranium"
  )
)


# ============================================================
# 5. Merge
# ============================================================

D <- merge(
  D,
  LABELS,
  by = "Environmental_factor",
  all.x = TRUE,
  sort = FALSE
)

if (
  any(
    is.na(
      D$Display_factor
    )
  )
) {
  stop(
    "存在未匹配的变量标签。",
    call. = FALSE
  )
}


# ============================================================
# 6. Desired order
# ============================================================

ORDER <- c(
  "Continental air-mass proportion",
  "Monthly sunspot number",
  "PM2.5",
  "Inverse planetary boundary-layer height",
  "Relative humidity",
  "Soil moisture",
  "Wind speed",
  "Temperature",
  "Temperature squared",
  "Temperature cubed",
  "Sunspots × latitude",
  "Surface uranium × soil moisture"
)

D[
  ,
  order_id :=
    match(
      Environmental_factor,
      ORDER
    )
]

setorder(
  D,
  order_id
)


# ============================================================
# 7. Format estimate and CI
# ============================================================

D[
  ,
  `:=`(
    Estimate_display =
      fmt(
        Estimate
      ),

    CI_display =
      paste0(
        "(",
        fmt(
          CI95_lower
        ),
        ", ",
        fmt(
          CI95_upper
        ),
        ")"
      )
  )
]


# ============================================================
# 8. Final publication table
# ============================================================

TAB <- D[
  ,
  .(
    `Environmental factor` =
      Display_factor,

    `Estimate` =
      Estimate_display,

    `95% CI` =
      CI_display,

    `Interpretation` =
      Interpretation
  )
]

fwrite(
  TAB,
  OUT_CSV
)


# ============================================================
# 9. Retrieve sample information
# ============================================================

N_ROWS <- AUDIT[
  item ==
    "environment_model_complete_rows",
  value
]

N_MON <- AUDIT[
  item ==
    "environment_model_monitors",
  value
]

N_FE <- AUDIT[
  item ==
    "monitor_year_FE_groups",
  value
]


# ============================================================
# 10. Draw publication-style table
# ============================================================

TITLE <- paste0(
  "Supplementary Table. Associations of environmental covariates ",
  "with ambient gross-beta particle radioactivity"
)

SUBTITLE <- paste0(
  "Filter-period analysis; Monitor × year fixed effects; ",
  "N = ",
  format(
    as.numeric(N_ROWS),
    big.mark = ","
  ),
  " observations from ",
  N_MON,
  " RadNet monitoring stations."
)

FOOTNOTE <- paste0(
  "Estimates are adjusted regression coefficients from the environmental-covariate model. ",
  "Standard errors are clustered by monitoring station. ",
  "Temperature, temperature², and temperature³ should be interpreted jointly as a polynomial. ",
  "Because Monitor × year fixed effects absorb time-invariant monitor characteristics and calendar-year main effects, ",
  "surface uranium, latitude, and calendar year are not separately estimated as main effects. ",
  "Interaction terms involving surface uranium or latitude remain estimable."
)


# ------------------------------------------------------------
# Table theme
# ------------------------------------------------------------

TTHEME <- ttheme_minimal(
  base_size = 10,

  core = list(
    fg_params = list(
      hjust = 0,
      x = 0.02
    ),

    bg_params = list(
      fill = "white",
      col = "grey85"
    )
  ),

  colhead = list(
    fg_params = list(
      fontface = "bold",
      hjust = 0,
      x = 0.02
    ),

    bg_params = list(
      fill = "grey95",
      col = "grey70"
    )
  )
)

TG <- tableGrob(
  TAB,
  rows = NULL,
  theme = TTHEME
)

# ------------------------------------------------------------
# Adjust column widths
# ------------------------------------------------------------

TG$widths[1] <- unit(
  5.2,
  "cm"
)

TG$widths[2] <- unit(
  2.5,
  "cm"
)

TG$widths[3] <- unit(
  4.2,
  "cm"
)

TG$widths[4] <- unit(
  8.0,
  "cm"
)


# ============================================================
# 11. Compose page
# ============================================================

make_page <- function() {

  grid.newpage()

  pushViewport(
    viewport(
      layout =
        grid.layout(
          nrow = 5,
          ncol = 1,
          heights = unit(
            c(
              0.8,
              0.6,
              8.4,
              1.5,
              0.3
            ),
            c(
              "cm",
              "cm",
              "null",
              "cm",
              "cm"
            )
          )
        )
    )
  )

  grid.text(
    TITLE,
    vp =
      viewport(
        layout.pos.row = 1
      ),
    x = unit(
      0.02,
      "npc"
    ),
    just = "left",
    gp =
      gpar(
        fontsize = 15,
        fontface = "bold"
      )
  )

  grid.text(
    SUBTITLE,
    vp =
      viewport(
        layout.pos.row = 2
      ),
    x = unit(
      0.02,
      "npc"
    ),
    just = "left",
    gp =
      gpar(
        fontsize = 10
      )
  )

  pushViewport(
    viewport(
      layout.pos.row = 3
    )
  )

  grid.draw(
    TG
  )

  upViewport()

  grid.text(
    FOOTNOTE,
    vp =
      viewport(
        layout.pos.row = 4
      ),
    x = unit(
      0.02,
      "npc"
    ),
    y = unit(
      0.95,
      "npc"
    ),
    just = c(
      "left",
      "top"
    ),
    gp =
      gpar(
        fontsize = 8.5
      )
  )

  upViewport()
}


# ============================================================
# 12. Save PNG
# ============================================================

png(
  filename = OUT_PNG,
  width = 3600,
  height = 2300,
  res = 300,
  bg = "white"
)

make_page()

dev.off()


# ============================================================
# 13. Save PDF
# ============================================================

pdf(
  file = OUT_PDF,
  width = 12,
  height = 7.7,
  useDingbats = FALSE
)

make_page()

dev.off()


# ============================================================
# 14. Audit
# ============================================================

CHECK <- data.table(
  item = c(
    "input_rows",
    "output_rows",
    "all_estimates_finite",
    "all_CI_finite",
    "publication_csv_exists",
    "publication_png_exists",
    "publication_pdf_exists"
  ),

  observed = c(
    nrow(D),
    nrow(TAB),
    sum(
      is.finite(
        D$Estimate
      )
    ),
    sum(
      is.finite(
        D$CI95_lower
      ) &
      is.finite(
        D$CI95_upper
      )
    ),
    file.exists(
      OUT_CSV
    ),
    file.exists(
      OUT_PNG
    ),
    file.exists(
      OUT_PDF
    )
  ),

  expected = c(
    12,
    12,
    12,
    12,
    TRUE,
    TRUE,
    TRUE
  )
)

CHECK[
  ,
  status :=
    ifelse(
      as.character(
        observed
      ) ==
        as.character(
          expected
        ),
      "PASS",
      "FAIL"
    )
]

fwrite(
  CHECK,
  file.path(
    IN_DIR,
    "09_TableS_publication_audit.csv"
  )
)

print(
  TAB,
  nrows = Inf
)

cat(
  "\n============================================================\n",
  "07H18B PUBLICATION TABLE\n",
  "============================================================\n",
  sep = ""
)

print(
  CHECK,
  nrows = Inf
)

N_FAIL <- CHECK[
  status ==
    "FAIL",
  .N
]

cat(
  "\n",
  ifelse(
    N_FAIL == 0L,
    "07H18B TABLE BUILD: PASS\n",
    "07H18B TABLE BUILD: FAIL\n"
  ),
  "Output CSV = ",
  OUT_CSV,
  "\nOutput PNG = ",
  OUT_PNG,
  "\nOutput PDF = ",
  OUT_PDF,
  "\n",
  sep = ""
)

if (
  N_FAIL > 0L
) {
  stop(
    "07H18B存在FAIL。",
    call. = FALSE
  )
}

