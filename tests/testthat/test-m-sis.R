test_that("M-SIS produces a full deterministic ranking", {
  inputs <- example_inputs()
  fit <- m_sis(inputs, list_size = 7L)
  expect_s3_class(fit, "m_sis_fit")
  expect_length(fit$selected, 7L)
  expect_length(fit$ranking, inputs$p)
  expect_identical(fit$selected, fit$ranking[seq_len(7L)])
  selector <- hfsis:::.select_pathwise_window(inputs, integer(0L))
  r <- inputs$r_by_window[selector$index, ]
  expect_equal(fit$ranking_index, order(-abs(r), seq_along(r)))
  expect_true(fit$diagnostics$complete)
  expect_error(m_sis(inputs, list_size = 0), "between 1")
  expect_error(m_sis(inputs, list_size = inputs$p + 1), "between 1")
})

test_that("M-SIS is exact on a hand-computable prepared input", {
  inputs <- example_inputs()
  inputs$p <- 3L
  inputs$variable_names <- c("a", "b", "c")
  inputs$r_by_window <- matrix(
    rep(c(-0.2, 0.8, 0.5), each = length(inputs$candidate_windows)),
    nrow = length(inputs$candidate_windows),
    dimnames = list(as.character(inputs$candidate_windows), inputs$variable_names)
  )
  inputs$input_id <- "tiny"
  one <- m_sis(inputs, list_size = 1L)
  middle <- m_sis(inputs, list_size = 2L)
  full <- m_sis(inputs, list_size = 3L)
  expect_identical(one$selected, "b")
  expect_identical(middle$selected, c("b", "c"))
  expect_identical(full$selected, c("b", "c", "a"))
})

test_that("M-SIS accepts raw data and preserves user names", {
  data(hf_sis_example, package = "hfsis")
  dat <- hf_sis_example$data
  names <- paste0("factor_", seq_along(hf_sis_example$x_names))
  fit <- m_sis(as.matrix(dat[hf_sis_example$x_names]), dat$y,
    hf_sis_example$target_index, list_size = 3,
    control = hfsis_control(hf_sis_example$windows), variable_names = names)
  expect_true(all(fit$selected %in% names))
  prepared <- prepare_hfsis(
    as.matrix(dat[hf_sis_example$x_names]), dat$y,
    hf_sis_example$target_index, hfsis_control(hf_sis_example$windows),
    variable_names = names
  )
  expect_identical(fit$selected, m_sis(prepared, list_size = 3)$selected)
})

test_that("exact score ties are broken by original position", {
  set.seed(10)
  base <- rnorm(80)
  x <- cbind(first = base, second = base, third = rnorm(80))
  y <- base + rnorm(80, sd = 0.1)
  fit <- m_sis(x, y, 20, list_size = 2,
    control = hfsis_control(c(16L, 32L), truncation = "none"))
  expect_identical(fit$ranking[1:2], c("first", "second"))
})
