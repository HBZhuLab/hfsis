#' @export
print.hfsis_inputs <- function(x, ...) {
  cat("Prepared HF-SIS inputs\n")
  cat("  observations / variables:", x$n, "/", x$p, "\n")
  cat("  target index:", x$target_index, "\n")
  cat("  candidate windows:", paste(x$candidate_windows, collapse = ", "), "\n")
  cat("  truncation:", x$control$truncation, "\n")
  invisible(x)
}

#' @export
print.hfsis_fit <- function(x, ...) {
  cat(x$method, "fit\n")
  cat("  n / p:", x$inputs_summary$n, "/", x$inputs_summary$p, "\n")
  cat("  target index:", x$inputs_summary$target_index, "\n")
  cat("  retained list size:", length(x$selected), "\n")
  cat("  selected variables:", paste(x$selected, collapse = ", "), "\n")
  if (inherits(x, "pr_sis_fit")) {
    cat("  first / final window:", x$first_k, "/", x$final_k, "\n")
    cat("  updates completed:", x$iterations_completed, "\n")
    cat("  model cap / block size:", x$model_cap, "/", x$block_size, "\n")
    cat("  singularity:", x$singularity, "\n")
  } else {
    cat("  selected window:", x$selected_k, "\n")
  }
  invisible(x)
}
