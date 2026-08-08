test_that("prepared input and raw input interfaces retain the same convention", {
  data(hf_sis_example, package = "hfsis")
  dat <- hf_sis_example$data
  x <- as.matrix(dat[hf_sis_example$x_names])
  control <- hfsis_control(hf_sis_example$windows)
  prepared <- prepare_hfsis(
    x, dat[[hf_sis_example$response]],
    hf_sis_example$target_index, control
  )
  expect_identical(prepared$target_index, hf_sis_example$target_index)
  expect_error(prepare_hfsis(x, dat$y,
    nrow(x) - max(hf_sis_example$windows) + 2L, control), "extends beyond")
})
