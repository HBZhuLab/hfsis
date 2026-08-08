test_that("frequency output exposes counts and rank summaries", {
  fits <- example_fits()
  result <- selection_frequency(list(fits$m, fits$pr))
  required <- c(
    "variable", "method", "selected_count", "fit_count",
    "frequency", "mean_rank", "median_rank"
  )
  expect_identical(names(result), required)
})
