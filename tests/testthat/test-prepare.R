test_that("prepared inputs expose matrix-free window summaries", {
  inputs <- example_inputs()
  expect_s3_class(inputs, "hfsis_inputs")
  expect_named(inputs, c(
    "call", "n", "p", "target_index", "variable_names", "control",
    "threshold_scales", "thresholds", "candidate_windows", "truncated_x",
    "truncated_y", "x_variance_by_window", "y_variance_by_window",
    "r_by_window", "variance_floor_active_by_window", "input_id"
  ))
  expect_equal(dim(inputs$truncated_x), c(max(inputs$candidate_windows), inputs$p))
  expect_length(inputs$truncated_y, max(inputs$candidate_windows))
  expect_equal(
    dim(inputs$x_variance_by_window),
    c(length(inputs$candidate_windows), inputs$p)
  )
  expect_equal(dim(inputs$r_by_window), dim(inputs$x_variance_by_window))
  expect_equal(
    dim(inputs$variance_floor_active_by_window),
    c(length(inputs$candidate_windows), inputs$p + 1L)
  )
  expect_output(print(inputs), "Prepared HF-SIS inputs")
})

test_that("fixed thresholds and no truncation are explicit", {
  data(hf_sis_example, package = "hfsis")
  dat <- hf_sis_example$data
  x <- as.matrix(dat[hf_sis_example$x_names])
  fixed <- hfsis_control(hf_sis_example$windows, truncation = "fixed",
    scales = rep(2, ncol(x) + 1L))
  fixed_inputs <- prepare_hfsis(x, dat$y, hf_sis_example$target_index, fixed)
  expect_equal(unname(fixed_inputs$threshold_scales), rep(2, ncol(x) + 1L))
  expect_error(prepare_hfsis(x, dat$y, hf_sis_example$target_index,
    hfsis_control(hf_sis_example$windows, truncation = "fixed", scales = 1)),
    "length")

  none <- prepare_hfsis(x, dat$y, hf_sis_example$target_index,
    hfsis_control(hf_sis_example$windows, truncation = "none"))
  expect_true(all(is.infinite(none$thresholds)))
})
