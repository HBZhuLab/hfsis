test_that("screen comparisons preserve rank semantics", {
  fits <- example_fits()
  comparison <- compare_screens(fits$m, fits$pr)
  expect_s3_class(comparison, "hfsis_comparison")
  expect_equal(nrow(comparison$rank_table), 25L)
  expect_equal(comparison$overlap_size,
    length(intersect(fits$m$selected, fits$pr$selected)))
  expect_true(all(is.finite(comparison$rank_table$rank_pr)))
  expect_true(comparison$same_initial_window)
  expect_output(print(comparison), "overlap")
  expect_s3_class(summary(comparison), "summary_hfsis_comparison")
  expect_output(print(summary(comparison)), "Summary")
})

test_that("comparison rejects mismatched inputs", {
  fits <- example_fits()
  altered <- fits$pr
  altered$inputs_summary$target_index <- altered$inputs_summary$target_index + 1L
  expect_error(compare_screens(fits$m, altered), "same input")
  expect_error(compare_screens(fits$pr, fits$m), "m_sis_fit")
})
