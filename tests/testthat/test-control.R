test_that("control defaults match the implementation brief", {
  control <- hfsis_control()
  expect_identical(control$windows, c(32L, 64L, 128L))
  expect_null(control$delta)
  expect_identical(control$truncation, "bipower")
  expect_equal(control$threshold_multiplier, 4)
  expect_equal(control$threshold_exponent, 0.47)
})
