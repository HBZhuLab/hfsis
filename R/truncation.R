.estimate_bipower_scales <- function(z, span) {
  products <- abs(z[-1L, , drop = FALSE] * z[-nrow(z), , drop = FALSE])
  sqrt(pmax((pi / 2) * colSums(products) / span, 1e-12))
}

.resolve_thresholds <- function(z, control, delta) {
  if (control$truncation == "none") {
    return(list(scales = rep(Inf, ncol(z)), thresholds = rep(Inf, ncol(z))))
  }
  if (control$truncation == "fixed") {
    if (length(control$scales) != ncol(z)) {
      .hfsis_stop("Fixed `scales` must have length ncol(x) + 1 for the response.")
    }
    scales <- control$scales
  } else {
    span <- nrow(z) * delta
    scales <- .estimate_bipower_scales(z, span)
  }
  thresholds <- control$threshold_multiplier * scales *
    delta^control$threshold_exponent
  list(scales = scales, thresholds = thresholds)
}

.truncate_coordinates <- function(block, thresholds) {
  keep <- abs(block) <= matrix(
    thresholds, nrow = nrow(block), ncol = ncol(block), byrow = TRUE
  )
  block * keep
}
