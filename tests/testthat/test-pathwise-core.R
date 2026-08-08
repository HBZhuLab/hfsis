test_that("marginal summaries equal the full reference", {
  inputs <- .make_pathwise_test_inputs()
  reference <- .reference_correlations(inputs)
  for (index in seq_along(inputs$candidate_windows)) {
    expected <- reference[[index]][seq_len(inputs$p), inputs$p + 1L]
    expect_equal(
      unname(inputs$r_by_window[index, ]),
      unname(expected),
      tolerance = 1e-12
    )
  }
})

test_that("requested columns equal the full reference", {
  inputs <- .make_pathwise_test_inputs()
  reference <- .reference_correlations(inputs)
  for (active_set in list(integer(0L), 1L, c(1L, 2L), c(1L, 2L, 3L, 4L))) {
    cache <- hfsis:::.new_column_cache()
    columns <- hfsis:::.get_active_columns(inputs, active_set, cache)
    for (index in seq_along(inputs$candidate_windows)) {
      expected <- reference[[index]][
        seq_len(inputs$p), active_set, drop = FALSE
      ]
      expect_equal(unname(columns[[index]]), unname(expected), tolerance = 1e-12)
    }
  }
})

test_that("active-set selector equals the restricted full reference", {
  inputs <- .make_pathwise_test_inputs()
  for (active_set in list(integer(0L), 1L, c(1L, 2L), c(1L, 2L, 3L))) {
    matrix_free <- hfsis:::.select_pathwise_window(inputs, active_set)
    reference <- .reference_pathwise_selector(inputs, active_set)
    expect_identical(matrix_free$selected_k, reference$selected_k)
    expect_equal(
      matrix_free$diagnostics$max_discrepancy_to_smaller,
      reference$discrepancies,
      tolerance = 1e-12
    )
  }
})

test_that("partial-residual scores equal the restricted full reference", {
  inputs <- .make_pathwise_test_inputs()
  active_set <- c(1L, 3L)
  selector <- hfsis:::.select_pathwise_window(inputs, active_set)
  cache <- hfsis:::.new_column_cache()
  columns <- hfsis:::.active_columns_at_window(
    inputs, active_set, selector$index, cache
  )
  r <- inputs$r_by_window[selector$index, ]
  matrix_free <- hfsis:::.partial_residual_scores_internal(
    r, columns, active_set
  )$scores
  correlation <- .reference_correlations(inputs)[[selector$index]]
  reference <- .reference_partial_residual_scores(r, correlation, active_set)
  expect_equal(matrix_free, reference, tolerance = 1e-12)
})
