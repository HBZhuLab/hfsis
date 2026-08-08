#' Control adaptive HF-SIS estimation
#'
#' Creates a validated control object shared by [prepare_hfsis()], [m_sis()],
#' and [pr_sisis()]. Inputs are equally spaced increments or returns. No control
#' value is changed silently.
#'
#' @param windows Strictly increasing positive integer candidate window sizes.
#' @param delta Observation spacing in normalized time. `NULL` uses `1 / n`,
#'   where `n` is the number of input rows.
#' @param confidence_index Nonnegative confidence index used in the simultaneous
#'   stochastic radius.
#' @param lepski_constant Positive multiplier in the Lepski comparison.
#' @param variance_floor Positive floor used only in correlation normalization.
#' @param truncation One of `"bipower"`, `"fixed"`, or `"none"`. The last
#'   option explicitly removes jump protection.
#' @param threshold_multiplier Positive multiplier for coordinate thresholds.
#' @param threshold_exponent Nonnegative power of `delta` in the thresholds.
#' @param scales For fixed truncation, one positive scale per augmented
#'   coordinate `(X, Y)`. Ignored otherwise.
#'
#' @return An object of class `hfsis_control`.
#' @references Zhu, H. *Ultra-high-dimensional spot support screening in
#'   continuous-time regression*. Manuscript.
#' @seealso [prepare_hfsis()]
#' @examples
#' hfsis_control(windows = c(16L, 32L), truncation = "none")
#' @export
hfsis_control <- function(
    windows = c(32L, 64L, 128L),
    delta = NULL,
    confidence_index = 0,
    lepski_constant = 1,
    variance_floor = 1e-12,
    truncation = c("bipower", "fixed", "none"),
    threshold_multiplier = 4,
    threshold_exponent = 0.47,
    scales = NULL) {
  windows <- .validate_windows(windows)
  truncation <- match.arg(truncation)
  scalar_nonnegative <- function(value, name) {
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) || value < 0) {
      .hfsis_stop(sprintf("`%s` must be one nonnegative finite number.", name))
    }
    as.numeric(value)
  }
  confidence_index <- scalar_nonnegative(confidence_index, "confidence_index")
  threshold_exponent <- scalar_nonnegative(threshold_exponent, "threshold_exponent")
  if (length(lepski_constant) != 1L || !is.numeric(lepski_constant) ||
      !is.finite(lepski_constant) || lepski_constant <= 0) {
    .hfsis_stop("`lepski_constant` must be one positive finite number.")
  }
  if (length(variance_floor) != 1L || !is.numeric(variance_floor) ||
      !is.finite(variance_floor) || variance_floor <= 0) {
    .hfsis_stop("`variance_floor` must be one positive finite number.")
  }
  if (length(threshold_multiplier) != 1L || !is.numeric(threshold_multiplier) ||
      !is.finite(threshold_multiplier) || threshold_multiplier <= 0) {
    .hfsis_stop("`threshold_multiplier` must be one positive finite number.")
  }
  if (!is.null(delta)) .resolve_delta(delta, 1L)
  if (truncation == "fixed") {
    if (is.null(scales) || !is.numeric(scales) || any(!is.finite(scales)) ||
        any(scales <= 0)) {
      .hfsis_stop("Fixed truncation requires positive finite user-supplied `scales`.")
    }
    scales <- as.numeric(scales)
  }
  structure(
    list(
      windows = windows,
      delta = delta,
      confidence_index = confidence_index,
      lepski_constant = as.numeric(lepski_constant),
      variance_floor = as.numeric(variance_floor),
      truncation = truncation,
      threshold_multiplier = as.numeric(threshold_multiplier),
      threshold_exponent = threshold_exponent,
      scales = scales
    ),
    class = "hfsis_control"
  )
}

#' @export
print.hfsis_control <- function(x, ...) {
  cat("HF-SIS control\n")
  cat("  windows:", paste(x$windows, collapse = ", "), "\n")
  cat("  delta:", if (is.null(x$delta)) "1 / n (resolved at fit)" else x$delta, "\n")
  cat("  truncation:", x$truncation, "\n")
  if (identical(x$truncation, "fixed")) {
    cat("  fixed scales:", paste(x$scales, collapse = ", "), "\n")
  }
  cat("  threshold multiplier/exponent:", x$threshold_multiplier, "/",
      x$threshold_exponent, "\n")
  cat("  confidence index:", x$confidence_index, "\n")
  cat("  Lepski constant:", x$lepski_constant, "\n")
  cat("  variance floor:", format(x$variance_floor, scientific = TRUE), "\n")
  invisible(x)
}
