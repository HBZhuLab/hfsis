.stochastic_radius <- function(windows, p, confidence_index = 0) {
  logarithm <- 2 * log(exp(1) * (p + 1L)) +
    confidence_index + log(2 * length(windows))
  sqrt(logarithm / windows) + logarithm / windows
}

.active_set_discrepancy <- function(
    first_r,
    first_columns,
    second_r,
    second_columns) {
  r_distance <- max(abs(first_r - second_r))
  if (is.null(first_columns) || ncol(first_columns) == 0L) {
    return(r_distance)
  }
  max(r_distance, max(abs(first_columns - second_columns)))
}

.select_pathwise_window <- function(
    inputs,
    active_set,
    column_cache = NULL) {
  if (is.null(column_cache)) column_cache <- .new_column_cache()
  windows <- inputs$candidate_windows
  radius <- .stochastic_radius(
    windows,
    inputs$p,
    inputs$control$confidence_index
  )
  active_columns <- .get_active_columns(inputs, active_set, column_cache)
  accepted <- rep(TRUE, length(windows))
  max_discrepancy <- numeric(length(windows))
  for (k_index in seq_along(windows)) {
    distances <- vapply(seq_len(k_index), function(l_index) {
      .active_set_discrepancy(
        inputs$r_by_window[k_index, ], active_columns[[k_index]],
        inputs$r_by_window[l_index, ], active_columns[[l_index]]
      )
    }, numeric(1L))
    max_discrepancy[k_index] <- max(distances)
    if (k_index > 1L) {
      bounds <- inputs$control$lepski_constant *
        (radius[k_index] + radius[seq_len(k_index)])
      accepted[k_index] <- all(distances <= bounds)
    }
  }
  selected_index <- max(which(accepted))
  list(
    index = selected_index,
    selected_k = windows[selected_index],
    diagnostics = data.frame(
      window = windows,
      stochastic_radius = radius,
      max_discrepancy_to_smaller = max_discrepancy,
      admissible = accepted,
      selected = seq_along(windows) == selected_index,
      stringsAsFactors = FALSE
    )
  )
}
