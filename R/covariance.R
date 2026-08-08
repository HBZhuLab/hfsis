.local_truncated_covariance <- function(z, target_index, k, delta, thresholds) {
  rows <- target_index + seq_len(k) - 1L
  block <- z[rows, , drop = FALSE]
  truncated <- .truncate_coordinates(block, thresholds)
  covariance <- crossprod(truncated) / (k * delta)
  (covariance + t(covariance)) / 2
}

.cumulative_columns <- function(x) {
  result <- apply(x, 2L, cumsum)
  if (is.null(dim(result))) {
    result <- matrix(result, ncol = ncol(x))
  }
  result
}

.compute_base_window_summaries <- function(
    truncated_x,
    truncated_y,
    windows,
    delta,
    variance_floor) {
  p <- ncol(truncated_x)
  denominators <- windows * delta
  cumulative_x_squared <- .cumulative_columns(truncated_x^2)
  cumulative_xy <- .cumulative_columns(truncated_x * truncated_y)
  cumulative_y_squared <- cumsum(truncated_y^2)

  x_variance <- cumulative_x_squared[windows, , drop = FALSE] / denominators
  y_variance <- cumulative_y_squared[windows] / denominators
  xy_covariance <- cumulative_xy[windows, , drop = FALSE] / denominators
  x_variance_floored <- pmax(x_variance, variance_floor)
  y_variance_floored <- pmax(y_variance, variance_floor)
  r_by_window <- xy_covariance /
    sqrt(x_variance_floored * y_variance_floored)
  floor_active <- cbind(
    x_variance <= variance_floor,
    response = y_variance <= variance_floor
  )

  dim(x_variance) <- c(length(windows), p)
  dim(r_by_window) <- c(length(windows), p)
  list(
    x_variance_by_window = x_variance,
    y_variance_by_window = as.numeric(y_variance),
    r_by_window = r_by_window,
    variance_floor_active_by_window = floor_active
  )
}

.compute_active_column <- function(inputs, active_index) {
  if (!.is_count(active_index) || active_index > inputs$p) {
    .hfsis_stop("`active_index` must identify one input variable.")
  }
  active_index <- as.integer(active_index)
  products <- inputs$truncated_x * inputs$truncated_x[, active_index]
  cumulative_products <- .cumulative_columns(products)
  covariance_by_window <-
    cumulative_products[inputs$candidate_windows, , drop = FALSE] /
    (inputs$candidate_windows * inputs$control$delta)
  denominator <- sqrt(
    pmax(inputs$x_variance_by_window, inputs$control$variance_floor) *
      pmax(
        inputs$x_variance_by_window[, active_index],
        inputs$control$variance_floor
      )
  )
  result <- covariance_by_window / denominator
  dim(result) <- c(length(inputs$candidate_windows), inputs$p)
  colnames(result) <- inputs$variable_names
  result
}

.new_column_cache <- function() {
  cache <- new.env(parent = emptyenv())
  assign(".computed_indices", integer(0L), envir = cache)
  cache
}

.get_active_column <- function(inputs, active_index, cache) {
  key <- as.character(active_index)
  if (!exists(key, envir = cache, inherits = FALSE)) {
    assign(key, .compute_active_column(inputs, active_index), envir = cache)
    computed <- get(".computed_indices", envir = cache, inherits = FALSE)
    assign(
      ".computed_indices",
      c(computed, as.integer(active_index)),
      envir = cache
    )
  }
  get(key, envir = cache, inherits = FALSE)
}

.get_active_columns <- function(inputs, active_set, cache) {
  if (!length(active_set)) {
    return(lapply(seq_along(inputs$candidate_windows), function(index) {
      matrix(
        numeric(0L),
        nrow = inputs$p,
        ncol = 0L,
        dimnames = list(inputs$variable_names, character(0L))
      )
    }))
  }
  columns <- lapply(active_set, function(active_index) {
    .get_active_column(inputs, active_index, cache)
  })
  lapply(seq_along(inputs$candidate_windows), function(window_index) {
    result <- vapply(columns, function(column) {
      column[window_index, ]
    }, numeric(inputs$p))
    if (length(active_set) == 1L) {
      result <- matrix(result, ncol = 1L)
    }
    dimnames(result) <- list(
      inputs$variable_names,
      inputs$variable_names[active_set]
    )
    result
  })
}

.active_columns_at_window <- function(
    inputs,
    active_set,
    window_index,
    cache) {
  if (!length(active_set)) {
    return(matrix(
      numeric(0L),
      nrow = inputs$p,
      ncol = 0L,
      dimnames = list(inputs$variable_names, character(0L))
    ))
  }
  columns <- lapply(active_set, function(active_index) {
    .get_active_column(inputs, active_index, cache)[window_index, ]
  })
  result <- do.call(cbind, columns)
  dimnames(result) <- list(
    inputs$variable_names,
    inputs$variable_names[active_set]
  )
  result
}

.cached_active_indices <- function(cache) {
  get(".computed_indices", envir = cache, inherits = FALSE)
}

.normalize_augmented_correlation <- function(covariance, variance_floor) {
  raw_diagonal <- diag(covariance)
  diagonal <- pmax(raw_diagonal, variance_floor)
  correlation <- covariance / sqrt(outer(diagonal, diagonal))
  correlation <- (correlation + t(correlation)) / 2
  list(
    correlation = correlation,
    raw_diagonal = raw_diagonal,
    floored_diagonal = diagonal,
    floor_active = raw_diagonal <= variance_floor
  )
}
