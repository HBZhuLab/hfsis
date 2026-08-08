.partial_residual_scores_internal <- function(
    r,
    active_columns,
    selected) {
  if (length(selected) == 0L) {
    return(list(
      scores = r,
      singularity = FALSE,
      coefficients = numeric(0L)
    ))
  }
  R_AA <- active_columns[selected, , drop = FALSE]
  coefficients <- tryCatch(
    solve(R_AA, r[selected]),
    error = function(e) NULL
  )
  if (is.null(coefficients)) {
    return(list(
      scores = NULL,
      singularity = TRUE,
      coefficients = NULL
    ))
  }
  list(
    scores = r - drop(active_columns %*% coefficients),
    singularity = FALSE,
    coefficients = coefficients
  )
}

.empty_pr_path <- function() {
  data.frame(
    iteration = integer(),
    variable = character(),
    variable_index = integer(),
    score = double(),
    abs_score = double(),
    selected_k = integer(),
    active_size_before = integer(),
    model_size = integer(),
    stringsAsFactors = FALSE
  )
}

#' Partial-residual spot iterative screening (PR-SISIS)
#'
#' Performs iterative screening with active-set adaptive window selection.
#' At each update, the local window is selected using the marginal correlation
#' vector and the correlation columns associated with the current active set.
#' New correlation columns are computed on demand and cached.
#'
#' @param x Numeric covariate input or prepared `hfsis_inputs` object.
#' @param y Response increments for raw input; ignored for prepared input.
#' @param target_index First forward-window row for raw input; ignored for
#'   prepared input.
#' @param model_cap Maximum number of variables added to the active set.
#' @param block_size Positive number of variables added at each update.
#' @param control [hfsis_control()] for raw input.
#' @param variable_names Variable names in deterministic order.
#'
#' @return An object with classes `c("pr_sis_fit", "hfsis_fit", "hfsis")`.
#'   A complete fit contains a full variable ranking. If an active correlation
#'   block is singular, the ranking contains the variables selected before
#'   termination and `diagnostics$complete` is `FALSE`.
#' @references Zhu, H. *Ultra-high-dimensional spot support screening in
#'   continuous-time regression*. Manuscript.
#' @seealso [m_sis()], [screen_path()], [selected_window()]
#' @examples
#' data(hf_sis_example)
#' dat <- hf_sis_example$data
#' inputs <- prepare_hfsis(as.matrix(dat[hf_sis_example$x_names]), dat$y,
#'   hf_sis_example$target_index, hfsis_control(hf_sis_example$windows))
#' pr_sisis(inputs, model_cap = 6L, block_size = 2L)
#' @export
pr_sisis <- function(
    x,
    y = NULL,
    target_index = NULL,
    model_cap = 10L,
    block_size = 2L,
    control = hfsis_control(),
    variable_names = colnames(x)) {
  call <- match.call()
  inputs <- .as_prepared_inputs(
    x, y, target_index, control, variable_names, call
  )
  if (!.is_count(model_cap)) {
    .hfsis_stop("`model_cap` must be one positive integer.")
  }
  if (!.is_count(block_size)) {
    .hfsis_stop("`block_size` must be one positive integer.")
  }
  model_cap <- as.integer(min(model_cap, inputs$p))
  block_size <- as.integer(block_size)
  max_iterations <- ceiling(model_cap / block_size)

  p <- inputs$p
  selected <- integer(0L)
  column_cache <- .new_column_cache()
  entry_scores <- rep(NA_real_, p)
  selected_blocks <- list()
  score_history <- list()
  ranking_history <- list()
  path_rows <- list()
  window_rows <- list()
  iteration_windows <- integer(0L)
  singularity <- FALSE
  early_termination <- FALSE
  iterations_completed <- 0L
  first_k <- NA_integer_
  final_k <- NA_integer_
  final_scores <- rep(NA_real_, p)

  for (iteration in seq_len(max_iterations)) {
    if (length(selected) >= model_cap) break

    previous_size <- length(selected)
    selector <- .select_pathwise_window(
      inputs,
      active_set = selected,
      column_cache = column_cache
    )
    iteration_windows <- c(iteration_windows, selector$selected_k)
    if (is.na(first_k)) first_k <- selector$selected_k
    window_rows[[length(window_rows) + 1L]] <- data.frame(
      stage = "update",
      iteration = iteration,
      active_size = previous_size,
      selected_k = selector$selected_k,
      stringsAsFactors = FALSE
    )

    window_index <- selector$index
    r <- inputs$r_by_window[window_index, , drop = TRUE]
    active_columns <- .active_columns_at_window(
      inputs,
      active_set = selected,
      window_index = window_index,
      cache = column_cache
    )
    residual <- .partial_residual_scores_internal(
      r,
      active_columns,
      selected
    )
    if (residual$singularity) {
      singularity <- TRUE
      early_termination <- TRUE
      final_k <- selector$selected_k
      break
    }

    candidates <- setdiff(seq_len(p), selected)
    candidate_ranking <- .rank_abs_scores(residual$scores, candidates)
    number_to_add <- min(
      block_size,
      model_cap - length(selected),
      length(candidates)
    )
    if (number_to_add <= 0L) {
      early_termination <- length(selected) < model_cap
      break
    }
    added <- candidate_ranking[seq_len(number_to_add)]
    score_history[[iteration]] <- stats::setNames(
      residual$scores,
      inputs$variable_names
    )
    ranking_history[[iteration]] <- inputs$variable_names[candidate_ranking]
    entry_scores[added] <- residual$scores[added]
    selected_blocks[[iteration]] <- added
    path_rows[[length(path_rows) + 1L]] <- data.frame(
      iteration = iteration,
      variable = inputs$variable_names[added],
      variable_index = added,
      score = residual$scores[added],
      abs_score = abs(residual$scores[added]),
      selected_k = selector$selected_k,
      active_size_before = previous_size,
      model_size = previous_size + seq_along(added),
      stringsAsFactors = FALSE
    )
    selected <- c(selected, added)
    iterations_completed <- iteration
  }

  if (!singularity) {
    final_selector <- .select_pathwise_window(
      inputs,
      active_set = selected,
      column_cache = column_cache
    )
    final_k <- final_selector$selected_k
    window_rows[[length(window_rows) + 1L]] <- data.frame(
      stage = "final",
      iteration = NA_integer_,
      active_size = length(selected),
      selected_k = final_k,
      stringsAsFactors = FALSE
    )
    final_r <- inputs$r_by_window[final_selector$index, , drop = TRUE]
    final_columns <- .active_columns_at_window(
      inputs,
      active_set = selected,
      window_index = final_selector$index,
      cache = column_cache
    )
    final_residual <- .partial_residual_scores_internal(
      final_r,
      final_columns,
      selected
    )
    if (final_residual$singularity) {
      singularity <- TRUE
      early_termination <- TRUE
    } else {
      final_scores <- final_residual$scores
    }
  } else if (length(window_rows)) {
    window_rows[[length(window_rows) + 1L]] <- data.frame(
      stage = "final",
      iteration = NA_integer_,
      active_size = length(selected),
      selected_k = final_k,
      stringsAsFactors = FALSE
    )
  }

  if (!singularity) {
    remaining <- setdiff(seq_len(p), selected)
    ranking_index <- c(selected, .rank_abs_scores(final_scores, remaining))
    ranking <- inputs$variable_names[ranking_index]
  } else {
    ranking_index <- selected
    ranking <- inputs$variable_names[selected]
  }
  path <- if (length(path_rows)) do.call(rbind, path_rows) else .empty_pr_path()
  rownames(path) <- NULL
  window_path <- if (length(window_rows)) {
    do.call(rbind, window_rows)
  } else {
    data.frame(
      stage = character(),
      iteration = integer(),
      active_size = integer(),
      selected_k = integer(),
      stringsAsFactors = FALSE
    )
  }
  rownames(window_path) <- NULL
  complete <- !singularity && length(ranking_index) == p
  if (!complete) early_termination <- TRUE
  names(final_scores) <- inputs$variable_names
  names(entry_scores) <- inputs$variable_names
  active_columns_computed <- length(.cached_active_indices(column_cache))

  output <- list(
    method = "PR-SISIS",
    call = call,
    selected = inputs$variable_names[selected],
    selected_index = selected,
    ranking = ranking,
    ranking_index = ranking_index,
    scores = final_scores,
    entry_scores = entry_scores,
    path = path,
    window_path = window_path,
    iteration_windows = iteration_windows,
    first_k = first_k,
    final_k = final_k,
    score_history = score_history,
    ranking_history = ranking_history,
    selected_blocks = selected_blocks,
    iterations_completed = iterations_completed,
    singularity = singularity,
    early_termination = early_termination,
    model_cap = model_cap,
    block_size = block_size,
    active_columns_computed = active_columns_computed,
    inputs_summary = .fit_inputs_summary(inputs),
    diagnostics = list(
      complete = complete,
      output_size = length(selected),
      singularity = singularity,
      early_termination = early_termination,
      max_iterations = max_iterations
    )
  )
  class(output) <- c("pr_sis_fit", "hfsis_fit", "hfsis")
  output
}

#' @rdname pr_sisis
#' @export
pr_sis <- pr_sisis
