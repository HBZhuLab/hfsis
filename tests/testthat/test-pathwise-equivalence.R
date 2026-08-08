test_that("matrix-free PR-SISIS equals the restricted full path", {
  inputs <- .make_pathwise_test_inputs(p = 10L)
  matrix_free <- pr_sisis(inputs, model_cap = 6L, block_size = 2L)
  reference <- .reference_pathwise_runner(inputs, model_cap = 6L, block_size = 2L)

  expect_identical(matrix_free$selected_index, reference$selected)
  expect_identical(matrix_free$iteration_windows, reference$iteration_windows)
  expect_equal(
    unname(matrix_free$entry_scores),
    reference$entry_scores,
    tolerance = 1e-12
  )
  expect_identical(matrix_free$final_k, reference$final_k)
  expect_equal(
    unname(matrix_free$scores),
    unname(reference$final_scores),
    tolerance = 1e-12
  )
  expect_identical(matrix_free$ranking_index, reference$ranking)
  rownames(reference$path) <- NULL
  expect_equal(matrix_free$path, reference$path, tolerance = 1e-12)
})
