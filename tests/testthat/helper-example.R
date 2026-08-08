example_inputs <- function(truncation = "bipower") {
  data(hf_sis_example, package = "hfsis", envir = environment())
  dat <- hf_sis_example$data
  control <- hfsis_control(
    hf_sis_example$windows,
    truncation = truncation
  )
  prepare_hfsis(
    as.matrix(dat[hf_sis_example$x_names]), dat$y,
    hf_sis_example$target_index, control
  )
}

example_fits <- function() {
  inputs <- example_inputs()
  list(
    m = m_sis(inputs, list_size = 5L),
    pr = pr_sisis(inputs, model_cap = 5L, block_size = 2L)
  )
}

fixture_correlation <- function(covariance, variance_floor = 1e-12) {
  diagonal <- pmax(diag(covariance), variance_floor)
  covariance / sqrt(outer(diagonal, diagonal))
}
