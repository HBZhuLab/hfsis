.validate_fit <- function(object) {
  if (!inherits(object, "hfsis_fit")) {
    .hfsis_stop("`object` must inherit from `hfsis_fit`.")
  }
  invisible(object)
}

.validate_iteration <- function(iteration, maximum) {
  if (!.is_count(iteration) || iteration > maximum) {
    .hfsis_stop(sprintf("`iteration` must be an integer between 1 and %d.", maximum))
  }
  as.integer(iteration)
}

#' Extract selected variables and screening results
#'
#' These accessors return ordinary base-R objects while preserving the input
#' variable order used for deterministic tie breaking.
#'
#' @param object A fitted `hfsis_fit` object.
#' @param iteration For PR-SISIS, an optional completed iteration. For M-SIS,
#'   only `NULL` or `1` is admissible.
#'
#' @return `selected_variables()` returns a character vector.
#'   `screen_scores()` returns a data frame with one row per variable.
#'   `screen_ranking()` returns a ranked data frame. Complete PR-SISIS fits
#'   provide a full ordering; after singular termination, undefined ranks are
#'   omitted. `selected_window()` returns one integer. For PR-SISIS,
#'   `iteration = NULL` returns the final-ranking window. `screen_path()`
#'   returns a data frame describing selection over iterations.
#' @references Zhu, H. *Ultra-high-dimensional spot support screening in
#'   continuous-time regression*. Manuscript.
#' @name hfsis_extractors
NULL

#' @rdname hfsis_extractors
#' @export
selected_variables <- function(object) {
  .validate_fit(object)
  unname(object$selected)
}

#' @rdname hfsis_extractors
#' @export
screen_scores <- function(object, iteration = NULL) {
  .validate_fit(object)
  variables <- object$inputs_summary$variable_names
  if (inherits(object, "m_sis_fit")) {
    if (!is.null(iteration)) .validate_iteration(iteration, 1L)
    scores <- object$scores
  } else if (is.null(iteration)) {
    scores <- object$scores
  } else {
    iteration <- .validate_iteration(iteration, object$iterations_completed)
    scores <- object$score_history[[iteration]]
  }
  data.frame(
    variable = variables,
    score = unname(scores[variables]),
    abs_score = abs(unname(scores[variables])),
    selected = variables %in% object$selected,
    stringsAsFactors = FALSE
  )
}

#' @rdname hfsis_extractors
#' @export
screen_ranking <- function(object, iteration = NULL) {
  .validate_fit(object)
  if (inherits(object, "m_sis_fit")) {
    if (!is.null(iteration)) .validate_iteration(iteration, 1L)
    ranking <- object$ranking
    scores <- object$scores[ranking]
  } else if (is.null(iteration)) {
    ranking <- object$ranking
    scores <- object$scores[ranking]
  } else {
    iteration <- .validate_iteration(iteration, object$iterations_completed)
    ranking <- object$ranking_history[[iteration]]
    scores <- object$score_history[[iteration]][ranking]
  }
  data.frame(
    rank = seq_along(ranking),
    variable = unname(ranking),
    score = unname(scores),
    abs_score = abs(unname(scores)),
    selected = unname(ranking) %in% object$selected,
    stringsAsFactors = FALSE
  )
}

#' @rdname hfsis_extractors
#' @export
selected_window <- function(object, iteration = NULL) {
  .validate_fit(object)
  if (inherits(object, "m_sis_fit")) {
    if (!is.null(iteration)) .validate_iteration(iteration, 1L)
    return(as.integer(object$selected_k))
  }
  if (is.null(iteration)) return(as.integer(object$final_k))
  iteration <- .validate_iteration(iteration, length(object$iteration_windows))
  as.integer(object$iteration_windows[iteration])
}

#' @rdname hfsis_extractors
#' @export
screen_path <- function(object) {
  .validate_fit(object)
  if (inherits(object, "pr_sis_fit")) return(object$path)
  data.frame(
    iteration = 1L,
    variables_added = paste(object$selected, collapse = ", "),
    model_size = length(object$selected),
    selected_k = object$selected_k,
    stringsAsFactors = FALSE
  )
}
