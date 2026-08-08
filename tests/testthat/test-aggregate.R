test_that("repeated-fit frequencies retain method denominators and metadata order", {
  data(hf_sis_panel_example, package = "hfsis")
  fits <- lapply(hf_sis_panel_example$sessions[1:3], function(dat) {
    m_sis(as.matrix(dat[hf_sis_panel_example$x_names]), dat$y,
      hf_sis_panel_example$target_index, list_size = 5,
      control = hfsis_control(hf_sis_panel_example$windows))
  })
  frequency <- selection_frequency(fits, metadata = hf_sis_panel_example$metadata)
  expect_s3_class(frequency, "hfsis_selection_frequency")
  expect_equal(frequency$variable, hf_sis_panel_example$x_names)
  expect_true(all(frequency$fit_count == 3L))
  expect_equal(sum(frequency$selected_count), 15L)
  expect_equal(frequency$frequency, frequency$selected_count / 3)
  expect_true(all(is.finite(frequency$mean_rank)))
  expect_equal(nrow(selection_frequency(fits, top = 4L)), 4L)
  share <- cluster_share(fits, hf_sis_panel_example$metadata)
  expect_s3_class(share, "hfsis_cluster_share")
  expect_equal(sum(share$selection_share), 1)
})

test_that("aggregation rejects ambiguous inputs and metadata", {
  fits <- example_fits()
  expect_error(selection_frequency(list()), "nonempty")
  bad <- fits$m
  bad$inputs_summary$variable_names <- rev(bad$inputs_summary$variable_names)
  expect_error(selection_frequency(list(fits$m, bad)), "same variable")
  expect_error(selection_frequency(list(fits$m),
    metadata = data.frame(variable = c("x01", "x01"))), "unique")
  expect_error(cluster_share(list(fits$m),
    data.frame(variable = fits$m$inputs_summary$variable_names)),
    "cluster_column")
})
