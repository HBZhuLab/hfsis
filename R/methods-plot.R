#' Plot an HF-SIS fit
#'
#' @param x A fitted hfsis_fit object.
#' @param type Plot type: leading absolute scores, adaptive-window
#'   diagnostics, or the selection path.
#' @param top Positive number of score labels displayed.
#' @param ... Additional graphical parameters passed to [graphics::barplot()]
#'   for score plots and [graphics::plot()] for window and path plots.
#'
#' @return x, invisibly. No graphics device is opened automatically.
#' @references Zhu, H. *Ultra-high-dimensional spot support screening in
#'   continuous-time regression*. Manuscript.
#' @seealso [screen_scores()], [screen_path()], [selected_window()]
#' @export
plot.hfsis_fit <- function(
    x, type = c("scores", "window", "path"), top = 20L, ...) {
  .validate_fit(x)
  type <- match.arg(type)
  if (!.is_count(top)) .hfsis_stop("'top' must be one positive integer.")
  if (type == "scores") {
    scores <- screen_scores(x)
    rows <- order(-scores$abs_score, seq_len(nrow(scores)), na.last = TRUE)
    rows <- utils::head(rows, top)
    colors <- ifelse(scores$selected[rows], "black", "grey75")
    graphics::barplot(
      rev(scores$abs_score[rows]), names.arg = rev(scores$variable[rows]),
      horiz = TRUE, las = 1, col = rev(colors),
      xlab = "Absolute spot-correlation score",
      main = paste(x$method, "screening scores"), ...
    )
  } else if (type == "window") {
    if (inherits(x, "pr_sis_fit")) {
      windows <- c(x$iteration_windows, x$final_k)
      stages <- seq_along(windows)
      graphics::plot(
        stages, windows,
        type = "b", pch = c(rep(19, length(x$iteration_windows)), 21),
        xlab = "PR-SISIS update",
        ylab = "Selected window",
        main = "PR-SISIS selected windows", ...
      )
    } else {
      diagnostics <- x$window_diagnostics
      colors <- ifelse(diagnostics$admissible, "black", "grey70")
      graphics::plot(
        diagnostics$window, diagnostics$max_discrepancy_to_smaller,
        type = "b", pch = 19, col = colors,
        xlab = "Candidate window (larger uses more increments)",
        ylab = "Maximum marginal discrepancy to smaller windows",
        main = "M-SIS Lepski window diagnostics", ...
      )
      selected <- diagnostics$selected
      graphics::points(diagnostics$window[selected],
        diagnostics$max_discrepancy_to_smaller[selected],
        pch = 21, bg = "white", cex = 1.5)
      graphics::legend("topright",
        legend = c("admissible", "not admissible", "selected"),
        pch = c(19, 19, 21), col = c("black", "grey70", "black"),
        pt.bg = c(NA, NA, "white"), bty = "n")
    }
  } else if (inherits(x, "pr_sis_fit")) {
    path <- screen_path(x)
    if (!nrow(path)) .hfsis_stop("The PR-SISIS path is empty.")
    graphics::plot(
      path$model_size, path$abs_score, type = "b",
      xlab = "Selected model size", ylab = "Absolute residual score at entry",
      main = "PR-SISIS selection path", ...
    )
    graphics::text(path$model_size, path$abs_score,
      labels = path$variable, pos = 3, cex = 0.7)
  } else {
    ranking <- screen_ranking(x)
    graphics::plot(
      ranking$rank, ranking$abs_score, type = "h",
      xlab = "One-stage M-SIS rank", ylab = "Absolute score",
      main = "M-SIS one-stage path", ...
    )
  }
  invisible(x)
}
