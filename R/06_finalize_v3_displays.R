# Entry-point wiring only: invoke the existing final Figure 1/6 renderers.
# No data construction, model fitting, or statistical transformations here.
source(file.path(ROOT, "R", "Figure1_display_v3.R"))
source(file.path(ROOT, "R", "Figure6_display_v3.R"))

render_figure1_v3(
  input_dir = file.path(ROOT, "source_data_manuscript_v5"),
  output_stem = file.path(ROOT, "outputs", "figures", "Figure_1_powerplants_RadNet_map"),
  qa_dir = file.path(VALIDATION_DIR, "figure1"))
render_figure6_v3(
  data = data.table::fread(file.path(ROOT, "source_data_manuscript_v5",
    "Figure_6_panel_sample_sizes_and_estimates.csv")),
  output_stem = file.path(ROOT, "outputs", "figures", "Figure_6_internal_validity_robustness"),
  qa_dir = file.path(VALIDATION_DIR, "figure6"),
  work_dir = file.path(ROOT, "outputs", "display_work", "figure6"))
