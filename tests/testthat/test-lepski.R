test_that("stochastic radius decreases and the smallest window is admissible", {
  windows <- c(16L, 32L, 64L, 128L)
  radius <- hfsis:::.stochastic_radius(windows, p = 25L)
  expect_true(all(diff(radius) < 0))

  fit <- m_sis(example_inputs(), list_size = 5L)
  expect_true(fit$window_diagnostics$admissible[1L])
  expect_true(fit$window_diagnostics$selected[
    match(fit$selected_k, fit$window_diagnostics$window)
  ])
})
