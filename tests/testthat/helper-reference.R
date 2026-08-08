.reference_correlations <- function(inputs) {
  z <- cbind(inputs$truncated_x, inputs$truncated_y)
  lapply(inputs$candidate_windows, function(window) {
    covariance <- crossprod(z[seq_len(window), , drop = FALSE]) /
      (window * inputs$control$delta)
    diagonal <- pmax(diag(covariance), inputs$control$variance_floor)
    correlation <- covariance / sqrt(outer(diagonal, diagonal))
    (correlation + t(correlation)) / 2
  })
}

.reference_pathwise_selector <- function(inputs, active_set) {
  correlations <- .reference_correlations(inputs)
  radius <- hfsis:::.stochastic_radius(
    inputs$candidate_windows,
    inputs$p,
    inputs$control$confidence_index
  )
  accepted <- rep(TRUE, length(inputs$candidate_windows))
  discrepancies <- numeric(length(inputs$candidate_windows))
  for (k_index in seq_along(inputs$candidate_windows)) {
    distances <- vapply(seq_len(k_index), function(l_index) {
      first <- correlations[[k_index]]
      second <- correlations[[l_index]]
      marginal_distance <- max(abs(
        first[seq_len(inputs$p), inputs$p + 1L] -
          second[seq_len(inputs$p), inputs$p + 1L]
      ))
      if (!length(active_set)) return(marginal_distance)
      max(
        marginal_distance,
        max(abs(
          first[seq_len(inputs$p), active_set, drop = FALSE] -
            second[seq_len(inputs$p), active_set, drop = FALSE]
        ))
      )
    }, numeric(1L))
    discrepancies[k_index] <- max(distances)
    if (k_index > 1L) {
      bounds <- inputs$control$lepski_constant *
        (radius[k_index] + radius[seq_len(k_index)])
      accepted[k_index] <- all(distances <= bounds)
    }
  }
  index <- max(which(accepted))
  list(
    index = index,
    selected_k = inputs$candidate_windows[index],
    correlations = correlations,
    discrepancies = discrepancies
  )
}

.reference_partial_residual_scores <- function(r, correlation, active_set) {
  if (!length(active_set)) return(r)
  coefficients <- solve(
    correlation[active_set, active_set, drop = FALSE],
    r[active_set]
  )
  r - drop(
    correlation[seq_along(r), active_set, drop = FALSE] %*% coefficients
  )
}

.reference_pathwise_runner <- function(inputs, model_cap, block_size) {
  selected <- integer(0L)
  entry_scores <- rep(NA_real_, inputs$p)
  iteration_windows <- integer(0L)
  path_rows <- list()
  max_iterations <- ceiling(model_cap / block_size)
  for (iteration in seq_len(max_iterations)) {
    selector <- .reference_pathwise_selector(inputs, selected)
    iteration_windows <- c(iteration_windows, selector$selected_k)
    correlation <- selector$correlations[[selector$index]]
    r <- correlation[seq_len(inputs$p), inputs$p + 1L]
    scores <- .reference_partial_residual_scores(r, correlation, selected)
    candidates <- setdiff(seq_len(inputs$p), selected)
    ranking <- candidates[order(-abs(scores[candidates]), candidates)]
    number_to_add <- min(
      block_size,
      model_cap - length(selected),
      length(candidates)
    )
    added <- ranking[seq_len(number_to_add)]
    previous_size <- length(selected)
    entry_scores[added] <- scores[added]
    path_rows[[iteration]] <- data.frame(
      iteration = iteration,
      variable = inputs$variable_names[added],
      variable_index = added,
      score = scores[added],
      abs_score = abs(scores[added]),
      selected_k = selector$selected_k,
      active_size_before = previous_size,
      model_size = previous_size + seq_along(added),
      stringsAsFactors = FALSE
    )
    selected <- c(selected, added)
    if (length(selected) >= model_cap) break
  }
  final_selector <- .reference_pathwise_selector(inputs, selected)
  final_correlation <- final_selector$correlations[[final_selector$index]]
  final_r <- final_correlation[seq_len(inputs$p), inputs$p + 1L]
  final_scores <- .reference_partial_residual_scores(
    final_r,
    final_correlation,
    selected
  )
  remaining <- setdiff(seq_len(inputs$p), selected)
  ranking <- c(
    selected,
    remaining[order(-abs(final_scores[remaining]), remaining)]
  )
  list(
    selected = selected,
    iteration_windows = iteration_windows,
    entry_scores = entry_scores,
    final_k = final_selector$selected_k,
    final_scores = final_scores,
    ranking = ranking,
    path = do.call(rbind, path_rows)
  )
}

.contains_dense_square <- function(object, p) {
  if ((is.matrix(object) || is.array(object)) && identical(dim(object), c(p, p))) {
    return(TRUE)
  }
  if (!is.list(object)) return(FALSE)
  any(vapply(object, .contains_dense_square, logical(1L), p = p))
}

.make_pathwise_test_inputs <- function(p = 8L) {
  set.seed(20260808)
  n <- 160L
  x <- matrix(stats::rnorm(n * p), nrow = n)
  x[, 2L] <- 0.55 * x[, 1L] + sqrt(1 - 0.55^2) * x[, 2L]
  y <- 0.7 * x[, 1L] - 0.4 * x[, 3L] + stats::rnorm(n, sd = 0.7)
  prepare_hfsis(
    x,
    y,
    target_index = 21L,
    control = hfsis_control(
      windows = c(24L, 48L, 96L),
      delta = 1 / n,
      truncation = "none",
      lepski_constant = 0.2
    ),
    variable_names = paste0("x", seq_len(p))
  )
}
