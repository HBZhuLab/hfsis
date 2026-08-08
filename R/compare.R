.comparison_rank <- function(fit, variables) {
  result <- rep(NA_integer_, length(variables))
  names(result) <- variables
  result[fit$ranking] <- seq_along(fit$ranking)
  result
}

#' Compare M-SIS and PR-SISIS screening fits
#'
#' The two fits must use exactly the same prepared input, variable order, and
#' target index. Complete PR-SISIS fits provide a full ordering. If PR-SISIS
#' terminates because of a singular active block, ranks not defined after
#' termination are returned as missing.
#'
#' @param m_sis_fit A fitted m_sis_fit object.
#' @param pr_sis_fit A fitted pr_sis_fit object.
#'
#' @return An object of class hfsis_comparison containing the two selected
#'   sets, their overlap and differences, paired ranks and scores, and window
#'   agreement diagnostics.
#' @references Zhu, H. *Ultra-high-dimensional spot support screening in
#'   continuous-time regression*. Manuscript.
#' @seealso [m_sis()], [pr_sisis()], [screen_ranking()]
#' @examples
#' data(hf_sis_example)
#' dat <- hf_sis_example$data
#' inputs <- prepare_hfsis(as.matrix(dat[hf_sis_example$x_names]), dat$y,
#'   hf_sis_example$target_index, hfsis_control(hf_sis_example$windows))
#' compare_screens(
#'   m_sis(inputs, list_size = 5),
#'   pr_sisis(inputs, model_cap = 5L, block_size = 2L)
#' )
#' @export
compare_screens <- function(m_sis_fit, pr_sis_fit) {
  if (!inherits(m_sis_fit, "m_sis_fit")) {
    .hfsis_stop("'m_sis_fit' must be returned by m_sis().")
  }
  if (!inherits(pr_sis_fit, "pr_sis_fit")) {
    .hfsis_stop("'pr_sis_fit' must be returned by pr_sisis().")
  }
  m_info <- m_sis_fit$inputs_summary
  pr_info <- pr_sis_fit$inputs_summary
  if (!identical(m_info$input_id, pr_info$input_id) ||
      !identical(m_info$variable_names, pr_info$variable_names) ||
      !identical(m_info$target_index, pr_info$target_index)) {
    .hfsis_stop("Fits must use the same input, variable order, and target index.")
  }
  variables <- m_info$variable_names
  rank_m <- .comparison_rank(m_sis_fit, variables)
  rank_pr <- .comparison_rank(pr_sis_fit, variables)
  selected_m <- unname(m_sis_fit$selected)
  selected_pr <- unname(pr_sis_fit$selected)
  overlap <- intersect(selected_m, selected_pr)
  rank_table <- data.frame(
    variable = variables,
    rank_m = unname(rank_m),
    rank_pr = unname(rank_pr),
    score_m = unname(m_sis_fit$scores[variables]),
    score_pr = unname(pr_sis_fit$scores[variables]),
    selected_m = variables %in% selected_m,
    selected_pr = variables %in% selected_pr,
    stringsAsFactors = FALSE
  )
  output <- list(
    call = match.call(),
    selected_m = selected_m,
    selected_pr = selected_pr,
    overlap = overlap,
    overlap_size = length(overlap),
    only_m = setdiff(selected_m, selected_pr),
    only_pr = setdiff(selected_pr, selected_m),
    rank_table = rank_table,
    selected_k_m = m_sis_fit$selected_k,
    pr_first_k = pr_sis_fit$first_k,
    pr_final_k = pr_sis_fit$final_k,
    same_initial_window = identical(m_sis_fit$selected_k, pr_sis_fit$first_k),
    input_id = m_info$input_id
  )
  class(output) <- c("hfsis_comparison", "hfsis")
  output
}

#' @export
print.hfsis_comparison <- function(x, ...) {
  cat("HF-SIS screening comparison\n")
  cat("  selected variables: M-SIS", length(x$selected_m),
      "/ PR-SISIS", length(x$selected_pr), "\n")
  cat("  overlap:", x$overlap_size, "\n")
  cat("  M-SIS / PR-SISIS first windows:", x$selected_k_m, "vs",
      x$pr_first_k,
      if (x$same_initial_window) "(same)" else "(different)", "\n")
  cat("  PR-SISIS final window:", x$pr_final_k, "\n")
  invisible(x)
}

#' @export
summary.hfsis_comparison <- function(object, ...) {
  union_size <- length(union(object$selected_m, object$selected_pr))
  output <- list(
    selected_sizes = c(M_SIS = length(object$selected_m),
      PR_SISIS = length(object$selected_pr)),
    overlap = object$overlap,
    overlap_size = object$overlap_size,
    only_m = object$only_m,
    only_pr = object$only_pr,
    jaccard = if (union_size) object$overlap_size / union_size else 1,
    same_initial_window = object$same_initial_window,
    windows = c(
      M_SIS = object$selected_k_m,
      PR_SISIS_first = object$pr_first_k,
      PR_SISIS_final = object$pr_final_k
    ),
    rank_table = object$rank_table
  )
  class(output) <- "summary_hfsis_comparison"
  output
}

#' @export
print.summary_hfsis_comparison <- function(x, ...) {
  cat("Summary of HF-SIS screening comparison\n")
  cat("  selected sizes:", paste(x$selected_sizes, collapse = " vs "), "\n")
  cat("  overlap:", x$overlap_size, "\n")
  cat("  Jaccard:", format(round(x$jaccard, 3), nsmall = 3), "\n")
  cat("  same initial window:", x$same_initial_window, "\n")
  invisible(x)
}

#' @export
plot.hfsis_comparison <- function(x, ...) {
  table <- x$rank_table
  keep <- is.finite(table$rank_m) & is.finite(table$rank_pr)
  if (!any(keep)) .hfsis_stop("There are no jointly ranked variables to plot.")
  graphics::plot(
    table$rank_m[keep], table$rank_pr[keep],
    xlab = "M-SIS rank", ylab = "PR-SISIS rank",
    main = "HF-SIS paired ranks", ...
  )
  graphics::abline(a = 0, b = 1, lty = 2, col = "grey50")
  invisible(x)
}
