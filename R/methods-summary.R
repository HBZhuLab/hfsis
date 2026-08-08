#' Summarize an HF-SIS fit
#'
#' @param object A fitted hfsis_fit object.
#' @param ... Unused.
#'
#' @return A programmatically accessible object of class summary_hfsis_fit
#'   containing controls, window diagnostics, leading scores, retained
#'   variables, PR-SISIS path information, and numerical diagnostics.
#' @references Zhu, H. *Ultra-high-dimensional spot support screening in
#'   continuous-time regression*. Manuscript.
#' @seealso [m_sis()], [pr_sisis()], [screen_scores()]
#' @export
summary.hfsis_fit <- function(object, ...) {
  path <- screen_path(object)
  iteration_sizes <- if (inherits(object, "pr_sis_fit") && nrow(path)) {
    as.integer(table(path$iteration))
  } else {
    length(object$selected)
  }
  output <- list(
    method = object$method,
    dimensions = c(n = object$inputs_summary$n, p = object$inputs_summary$p),
    target_index = object$inputs_summary$target_index,
    control = object$inputs_summary$control,
    candidate_windows = object$inputs_summary$candidate_windows,
    top_scores = utils::head(screen_ranking(object), 10L),
    retained_variables = object$selected,
    retained_size = length(object$selected),
    iteration_sizes = iteration_sizes,
    path = path,
    diagnostics = object$diagnostics
  )
  if (inherits(object, "pr_sis_fit")) {
    output$first_window <- object$first_k
    output$final_window <- object$final_k
    output$iteration_windows <- object$iteration_windows
    output$window_path <- object$window_path
    output$active_columns_computed <- object$active_columns_computed
  } else {
    output$selected_window <- object$selected_k
    output$window_diagnostics <- object$window_diagnostics
  }
  class(output) <- "summary_hfsis_fit"
  output
}

#' @export
print.summary_hfsis_fit <- function(x, ...) {
  cat("Summary of", x$method, "fit\n")
  cat("  n / p:", x$dimensions["n"], "/", x$dimensions["p"], "\n")
  cat("  target index:", x$target_index, "\n")
  cat("  candidate windows:", paste(x$candidate_windows, collapse = ", "), "\n")
  cat("  retained variables:", paste(x$retained_variables, collapse = ", "), "\n")
  if (identical(x$method, "PR-SISIS")) {
    cat("  first / final window:", x$first_window, "/", x$final_window, "\n")
    cat("  additions by iteration:", paste(x$iteration_sizes, collapse = ", "), "\n")
    cat("  active columns computed:", x$active_columns_computed, "\n")
    cat("  singularity / early termination:", x$diagnostics$singularity, "/",
        x$diagnostics$early_termination, "\n")
  } else {
    cat("  selected window:", x$selected_window, "\n")
  }
  invisible(x)
}
