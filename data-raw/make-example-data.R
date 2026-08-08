# Rebuild both package datasets and both CSV examples deterministically.
#
# Design parameters:
# - seed: 20260004;
# - sessions: 5;
# - increments per session: 390;
# - covariates: 25 in five correlated clusters;
# - active coefficients: 0.50, -0.52, 0.40, -0.35, 0.30;
# - correlation of the cancellation pair: 0.92;
# - response-noise scale: 0.15;
# - finite-activity jump rows: 83, 318, and 370;
# - target_index: 131;
# - candidate windows: 32, 64, and 128 increments.
source(file.path("R", "data.R"), local = TRUE)
example_seed <- 20260004L
set.seed(example_seed)
examples <- .make_hfsis_examples(seed = example_seed)
hf_sis_example <- examples$hf_sis_example
hf_sis_panel_example <- examples$hf_sis_panel_example
save(hf_sis_example,
  file = file.path("data", "hf_sis_example.rda"), compress = "xz")
save(hf_sis_panel_example,
  file = file.path("data", "hf_sis_panel_example.rda"), compress = "xz")
utils::write.csv(
  hf_sis_example$data,
  file.path("inst", "extdata", "hf_sis_example.csv"),
  row.names = FALSE
)
utils::write.csv(
  hf_sis_example$metadata,
  file.path("inst", "extdata", "factor_metadata_example.csv"),
  row.names = FALSE
)
