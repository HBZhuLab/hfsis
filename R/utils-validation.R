.hfsis_stop <- function(..., call. = FALSE) {
  stop(..., call. = call.)
}

.is_count <- function(x) {
  length(x) == 1L && is.numeric(x) && is.finite(x) &&
    x >= 1 && x == floor(x)
}

.validate_windows <- function(windows) {
  if (!is.numeric(windows) || length(windows) < 1L || any(!is.finite(windows)) ||
      any(windows <= 0) || any(windows != floor(windows)) ||
      (length(windows) > 1L && any(diff(windows) <= 0))) {
    .hfsis_stop("`windows` must be strictly increasing positive integers.")
  }
  as.integer(windows)
}

.validate_hfsis_data <- function(x, y, target_index, variable_names, windows) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1L)))) {
      .hfsis_stop("Every column of `x` must be numeric.")
    }
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    .hfsis_stop("`x` must be a numeric matrix or numeric data frame.")
  }
  storage.mode(x) <- "double"
  if (!is.numeric(y) || is.matrix(y) || is.data.frame(y)) {
    .hfsis_stop("`y` must be a numeric vector.")
  }
  y <- as.numeric(y)
  if (nrow(x) != length(y)) {
    .hfsis_stop("`x` and `y` must have the same number of rows.")
  }
  if (nrow(x) < 2L || ncol(x) < 1L) {
    .hfsis_stop("`x` must contain at least two rows and one column.")
  }
  if (any(!is.finite(x)) || any(!is.finite(y))) {
    .hfsis_stop("`x` and `y` must contain only finite values; rows are not dropped or imputed.")
  }
  if (!.is_count(target_index)) {
    .hfsis_stop("`target_index` must be one positive integer.")
  }
  target_index <- as.integer(target_index)
  if (any(windows < 2L)) {
    .hfsis_stop("Every candidate window must contain at least two increments.")
  }
  if (target_index + max(windows) - 1L > nrow(x)) {
    .hfsis_stop("A candidate window extends beyond the available rows.")
  }
  if (is.null(variable_names)) {
    variable_names <- paste0("x", seq_len(ncol(x)))
  }
  if (!is.character(variable_names) || length(variable_names) != ncol(x) ||
      anyNA(variable_names) || any(!nzchar(variable_names))) {
    .hfsis_stop("`variable_names` must provide one nonempty name per column of `x`.")
  }
  if (anyDuplicated(variable_names)) {
    .hfsis_stop("`variable_names` must not contain duplicates.")
  }
  list(
    x = x,
    y = y,
    target_index = target_index,
    variable_names = variable_names
  )
}

.resolve_delta <- function(delta, n) {
  if (is.null(delta)) return(1 / n)
  if (length(delta) != 1L || !is.numeric(delta) || !is.finite(delta) || delta <= 0) {
    .hfsis_stop("`delta` must be NULL or one positive finite number.")
  }
  as.numeric(delta)
}

.validate_list_size <- function(list_size, p) {
  if (!.is_count(list_size) || list_size > p) {
    .hfsis_stop("`list_size` must be an integer between 1 and the number of variables.")
  }
  as.integer(list_size)
}

.make_input_id <- function(z, target_index, variable_names) {
  weights <- ((seq_along(z) - 1L) %% 104729L) + 1L
  c1 <- sum(z)
  c2 <- sum(as.numeric(z) * weights)
  paste(
    nrow(z), ncol(z), target_index, paste(variable_names, collapse = "\r"),
    format(c1, digits = 17L, scientific = TRUE),
    format(c2, digits = 17L, scientific = TRUE),
    sep = "|"
  )
}

.fit_inputs_summary <- function(inputs) {
  list(
    n = inputs$n,
    p = inputs$p,
    target_index = inputs$target_index,
    variable_names = inputs$variable_names,
    input_id = inputs$input_id,
    candidate_windows = inputs$candidate_windows,
    control = inputs$control
  )
}

.as_prepared_inputs <- function(x, y, target_index, control, variable_names, call) {
  if (inherits(x, "hfsis_inputs")) return(x)
  if (is.null(y) || is.null(target_index)) {
    .hfsis_stop("Raw input requires both `y` and `target_index`.")
  }
  prepare_hfsis(
    x = x,
    y = y,
    target_index = target_index,
    control = control,
    variable_names = variable_names
  )
}
