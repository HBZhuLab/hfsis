test_that("installed examples have documented dimensions and finite values", {
  data(hf_sis_example, package = "hfsis")
  data(hf_sis_panel_example, package = "hfsis")
  expect_equal(dim(hf_sis_example$data), c(390L, 27L))
  expect_length(hf_sis_example$x_names, 25L)
  expect_identical(hf_sis_example$response, "y")
  expect_true(nzchar(hf_sis_example$description))
  expect_equal(hf_sis_example$true_support, sprintf("x%02d", 1:5))
  expect_true(all(hf_sis_example$true_support %in% names(hf_sis_example$data)))
  expect_lte(hf_sis_example$target_index + max(hf_sis_example$windows) - 1L,
    nrow(hf_sis_example$data))
  expect_true(all(vapply(hf_sis_example$data, function(x) all(is.finite(x)), logical(1L))))
  expect_length(hf_sis_panel_example$sessions, 5L)
  expect_true(all(vapply(hf_sis_panel_example$sessions,
    function(x) identical(dim(x), c(390L, 27L)), logical(1L))))

  regenerated <- hfsis:::.make_hfsis_examples(seed = 20260004L)
  expect_identical(regenerated$hf_sis_example, hf_sis_example)
  expect_identical(regenerated$hf_sis_panel_example, hf_sis_panel_example)

  prepared <- prepare_hfsis(
    as.matrix(hf_sis_example$data[hf_sis_example$x_names]),
    hf_sis_example$data[[hf_sis_example$response]],
    hf_sis_example$target_index,
    hfsis_control(hf_sis_example$windows)
  )
  selector <- hfsis:::.select_pathwise_window(prepared, integer(0L))
  marginal_scores <- prepared$r_by_window[selector$index, ]
  active_column <- hfsis:::.compute_active_column(prepared, 2L)[selector$index, ]
  marginal <- abs(marginal_scores[1L])
  residual <- abs(marginal_scores[1L] - active_column[1L] * marginal_scores[2L])
  expect_lt(marginal, 0.05)
  expect_gt(residual, 2 * marginal)
})
