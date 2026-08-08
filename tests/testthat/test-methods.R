test_that("summary keeps controls, windows, paths, and diagnostics", {
  fits <- example_fits()
  result <- summary(fits$pr)
  expect_true(all(c(
    "control", "candidate_windows", "first_window", "final_window",
    "iteration_windows", "window_path", "active_columns_computed",
    "top_scores", "retained_variables", "iteration_sizes", "path",
    "diagnostics"
  ) %in% names(result)))
})
