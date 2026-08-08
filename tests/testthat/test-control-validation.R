test_that("control validates every public option", {
  expect_s3_class(hfsis_control(), "hfsis_control")
  expect_error(hfsis_control(c(32, 16)), "strictly increasing")
  expect_error(hfsis_control(c(16, 16)), "strictly increasing")
  expect_error(hfsis_control(delta = 0), "positive")
  expect_error(hfsis_control(confidence_index = -1), "nonnegative")
  expect_error(hfsis_control(lepski_constant = 0), "positive")
  expect_error(hfsis_control(variance_floor = 0), "positive")
  expect_error(hfsis_control(threshold_multiplier = 0), "positive")
  expect_error(hfsis_control(threshold_exponent = -0.1), "nonnegative")
  expect_error(hfsis_control(truncation = "fixed"), "requires")
  expect_output(print(hfsis_control(truncation = "none")), "HF-SIS control")
})

test_that("raw data validation is strict and does not impute", {
  x <- matrix(rnorm(80), 20, 4)
  y <- rnorm(20)
  control <- hfsis_control(c(4L, 8L), truncation = "none")
  expect_error(prepare_hfsis(x, y[-1], 2, control), "same number")
  x[1, 1] <- NA_real_
  expect_error(prepare_hfsis(x, y, 2, control), "finite")
  x[1, 1] <- 0
  expect_error(prepare_hfsis(x, y, 14, control), "extends beyond")
  expect_error(prepare_hfsis(x, y, 2, control,
    variable_names = rep("x", 4)), "duplicates")
  bad <- as.data.frame(x)
  bad[[1]] <- letters[seq_len(nrow(bad))]
  expect_error(prepare_hfsis(bad, y, 2, control), "Every column")
  expect_error(prepare_hfsis(x, y, 2, list()), "hfsis_control")
})
