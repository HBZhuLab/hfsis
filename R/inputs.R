#' Prepare shared adaptive spot-correlation inputs
#'
#' Validates equally spaced increment or return data, applies coordinatewise
#' truncation, and computes the marginal and diagonal summaries used by the
#' candidate forward windows. `target_index` is the row of the first increment
#' included in every forward window.
#'
#' @param x Numeric matrix or data frame of covariate increments; rows are
#'   consecutive equally spaced observations.
#' @param y Numeric response-increment vector.
#' @param target_index Row index of the first increment in the forward window.
#' @param control An [hfsis_control()] object.
#' @param variable_names Unique covariate names in deterministic tie-breaking
#'   order.
#'
#' @return An object with classes `c("hfsis_inputs", "hfsis")` containing
#'   truncated local increments, thresholds, and marginal summaries for each
#'   candidate window.
#' @references Zhu, H. *Ultra-high-dimensional spot support screening in
#'   continuous-time regression*. Manuscript.
#' @seealso [hfsis_control()], [m_sis()], [pr_sisis()]
#' @examples
#' data(hf_sis_example)
#' dat <- hf_sis_example$data
#' x <- as.matrix(dat[hf_sis_example$x_names])
#' prepare_hfsis(x, dat$y, hf_sis_example$target_index,
#'   hfsis_control(hf_sis_example$windows))
#' @export
prepare_hfsis <- function(
    x,
    y,
    target_index,
    control = hfsis_control(),
    variable_names = colnames(x)) {
  call <- match.call()
  if (!inherits(control, "hfsis_control")) {
    .hfsis_stop("`control` must be created by hfsis_control().")
  }
  validated <- .validate_hfsis_data(
    x, y, target_index, variable_names, control$windows
  )
  x <- validated$x
  y <- validated$y
  target_index <- validated$target_index
  variable_names <- validated$variable_names
  n <- nrow(x)
  p <- ncol(x)
  delta <- .resolve_delta(control$delta, n)
  resolved_control <- control
  resolved_control$delta <- delta
  z <- cbind(x, y)
  threshold_info <- .resolve_thresholds(z, resolved_control, delta)
  k_max <- max(control$windows)
  rows <- target_index + seq_len(k_max) - 1L
  local_z <- z[rows, , drop = FALSE]
  truncated_z <- .truncate_coordinates(
    local_z,
    threshold_info$thresholds
  )
  truncated_x <- truncated_z[, seq_len(p), drop = FALSE]
  truncated_y <- truncated_z[, p + 1L]
  summaries <- .compute_base_window_summaries(
    truncated_x,
    truncated_y,
    control$windows,
    delta,
    control$variance_floor
  )
  dimnames(summaries$x_variance_by_window) <- list(
    as.character(control$windows), variable_names
  )
  dimnames(summaries$r_by_window) <- list(
    as.character(control$windows), variable_names
  )
  dimnames(summaries$variance_floor_active_by_window) <- list(
    as.character(control$windows), c(variable_names, "response")
  )
  output <- list(
    call = call,
    n = n,
    p = p,
    target_index = target_index,
    variable_names = variable_names,
    control = resolved_control,
    threshold_scales = stats::setNames(threshold_info$scales, c(variable_names, "response")),
    thresholds = stats::setNames(threshold_info$thresholds, c(variable_names, "response")),
    candidate_windows = control$windows,
    truncated_x = truncated_x,
    truncated_y = truncated_y,
    x_variance_by_window = summaries$x_variance_by_window,
    y_variance_by_window = summaries$y_variance_by_window,
    r_by_window = summaries$r_by_window,
    variance_floor_active_by_window =
      summaries$variance_floor_active_by_window,
    input_id = .make_input_id(z, target_index, variable_names)
  )
  class(output) <- c("hfsis_inputs", "hfsis")
  output
}
