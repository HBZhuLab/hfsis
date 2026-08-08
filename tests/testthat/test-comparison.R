test_that("comparison exposes the contracted public fields", {
  fits <- example_fits()
  comparison <- compare_screens(fits$m, fits$pr)
  required <- c(
    "selected_m", "selected_pr", "overlap", "overlap_size",
    "only_m", "only_pr", "rank_table", "same_initial_window",
    "selected_k_m", "pr_first_k", "pr_final_k"
  )
  expect_true(all(required %in% names(comparison)))
})
