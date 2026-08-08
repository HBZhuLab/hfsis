test_that("PR-SISIS follows block_size and model_cap", {
  inputs <- .make_pathwise_test_inputs(p = 24L)
  cases <- list(
    list(model_cap = 10L, block_size = 2L, updates = 5L, sizes = rep(2L, 5L)),
    list(model_cap = 20L, block_size = 2L, updates = 10L, sizes = rep(2L, 10L)),
    list(model_cap = 9L, block_size = 2L, updates = 5L, sizes = c(2L, 2L, 2L, 2L, 1L))
  )
  for (case in cases) {
    fit <- pr_sisis(
      inputs,
      model_cap = case$model_cap,
      block_size = case$block_size
    )
    expect_equal(fit$iterations_completed, case$updates)
    expect_equal(as.integer(table(fit$path$iteration)), case$sizes)
    expect_length(fit$selected, case$model_cap)
    expect_length(fit$ranking, inputs$p)
    expect_true(fit$diagnostics$complete)
  }
})

test_that("partial-residual helper uses requested columns", {
  active_columns <- matrix(c(
    1, 0.3, -0.2,
    0.3, 1, 0.4
  ), nrow = 3L, ncol = 2L)
  r <- c(0.5, -0.1, 0.2)
  selected <- c(1L, 2L)
  result <- hfsis:::.partial_residual_scores_internal(
    r,
    active_columns,
    selected
  )
  expected <- r - drop(active_columns %*%
    solve(active_columns[selected, , drop = FALSE], r[selected]))
  expect_false(result$singularity)
  expect_equal(result$scores, expected)
})

test_that("active columns are cached once", {
  inputs <- .make_pathwise_test_inputs()
  cache <- hfsis:::.new_column_cache()
  hfsis:::.get_active_column(inputs, 2L, cache)
  hfsis:::.get_active_column(inputs, 2L, cache)
  hfsis:::.get_active_columns(inputs, c(2L, 3L), cache)
  hfsis:::.get_active_columns(inputs, c(2L, 3L), cache)
  expect_identical(hfsis:::.cached_active_indices(cache), c(2L, 3L))

  fit <- pr_sisis(inputs, model_cap = 6L, block_size = 2L)
  expect_lte(fit$active_columns_computed, fit$model_cap)
  expect_identical(fit$active_columns_computed, length(fit$selected_index))
})

test_that("singular active blocks stop without a complete ranking", {
  set.seed(42)
  n <- 100L
  first <- rnorm(n)
  x <- cbind(first, first, rnorm(n), rnorm(n))
  y <- first + rnorm(n, sd = 0.05)
  fit <- pr_sisis(
    x,
    y,
    target_index = 1L,
    model_cap = 4L,
    block_size = 2L,
    control = hfsis_control(c(24L, 48L), truncation = "none"),
    variable_names = paste0("x", seq_len(ncol(x)))
  )
  expect_length(fit$selected, 2L)
  expect_true(fit$singularity)
  expect_true(fit$early_termination)
  expect_false(fit$diagnostics$complete)
  expect_length(fit$ranking, 2L)
})

test_that("prepared and fitted objects contain no dense p by p object", {
  set.seed(11)
  n <- 80L
  p <- 1000L
  x <- matrix(rnorm(n * p), nrow = n)
  y <- rnorm(n)
  inputs <- prepare_hfsis(
    x,
    y,
    target_index = 1L,
    control = hfsis_control(c(16L, 32L), truncation = "none")
  )
  fit <- pr_sisis(inputs, model_cap = 4L, block_size = 2L)
  expect_false(.contains_dense_square(inputs, p))
  expect_false(.contains_dense_square(fit, p))
})

test_that("raw and prepared PR-SISIS interfaces agree", {
  inputs <- .make_pathwise_test_inputs()
  raw <- pr_sisis(
    inputs$truncated_x,
    inputs$truncated_y,
    target_index = 1L,
    model_cap = 6L,
    block_size = 2L,
    control = hfsis_control(
      inputs$candidate_windows,
      delta = inputs$control$delta,
      truncation = "none",
      lepski_constant = inputs$control$lepski_constant
    ),
    variable_names = inputs$variable_names
  )
  prepared <- prepare_hfsis(
    inputs$truncated_x,
    inputs$truncated_y,
    target_index = 1L,
    control = hfsis_control(
      inputs$candidate_windows,
      delta = inputs$control$delta,
      truncation = "none",
      lepski_constant = inputs$control$lepski_constant
    ),
    variable_names = inputs$variable_names
  )
  fit <- pr_sisis(prepared, model_cap = 6L, block_size = 2L)
  expect_identical(raw$selected_index, fit$selected_index)
  expect_identical(raw$iteration_windows, fit$iteration_windows)
  expect_equal(raw$path, fit$path, tolerance = 1e-12)
  expect_identical(raw$ranking_index, fit$ranking_index)
})

test_that("pr_sis is an alias of PR-SISIS", {
  inputs <- .make_pathwise_test_inputs()
  first <- pr_sisis(inputs, model_cap = 4L, block_size = 2L)
  second <- pr_sis(inputs, model_cap = 4L, block_size = 2L)
  first$call <- NULL
  second$call <- NULL
  expect_identical(first, second)
})
