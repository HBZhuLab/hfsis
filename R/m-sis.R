.rank_abs_scores <- function(scores, candidates = seq_along(scores)) {
  candidates[order(-abs(scores[candidates]), candidates)]
}

#' Marginal spot independence screening (M-SIS)
#'
#' Selects a local window using the marginal correlation vector and ranks
#' variables by the absolute marginal scores at that window. Equal scores are
#' resolved by original variable order. Supply raw increments or an object
#' returned by [prepare_hfsis()].
#'
#' @param x Numeric covariate input or prepared `hfsis_inputs` object.
#' @param y Response increments for raw input; ignored for prepared input.
#' @param target_index First forward-window row for raw input; ignored for
#'   prepared input.
#' @param list_size Number of variables retained.
#' @param control [hfsis_control()] for raw input.
#' @param variable_names Variable names in deterministic order.
#'
#' @return An object with classes `c("m_sis_fit", "hfsis_fit", "hfsis")`.
#' @references Zhu, H. *Ultra-high-dimensional spot support screening in
#'   continuous-time regression*. Manuscript.
#' @seealso [prepare_hfsis()], [pr_sisis()], [selected_variables()]
#' @examples
#' data(hf_sis_example)
#' dat <- hf_sis_example$data
#' fit <- m_sis(as.matrix(dat[hf_sis_example$x_names]), dat$y,
#'   hf_sis_example$target_index, list_size = 5,
#'   control = hfsis_control(hf_sis_example$windows))
#' selected_variables(fit)
#' @export
m_sis <- function(
    x,
    y = NULL,
    target_index = NULL,
    list_size = 10L,
    control = hfsis_control(),
    variable_names = colnames(x)) {
  call <- match.call()
  inputs <- .as_prepared_inputs(
    x, y, target_index, control, variable_names, call
  )
  list_size <- .validate_list_size(list_size, inputs$p)
  selector <- .select_pathwise_window(
    inputs,
    active_set = integer(0L)
  )
  r_hat <- inputs$r_by_window[selector$index, , drop = TRUE]
  ranking_index <- .rank_abs_scores(r_hat)
  selected_index <- ranking_index[seq_len(list_size)]
  output <- list(
    method = "M-SIS",
    call = call,
    selected = inputs$variable_names[selected_index],
    selected_index = selected_index,
    ranking = inputs$variable_names[ranking_index],
    ranking_index = ranking_index,
    scores = stats::setNames(as.numeric(r_hat), inputs$variable_names),
    selected_k = selector$selected_k,
    window_diagnostics = selector$diagnostics,
    inputs_summary = .fit_inputs_summary(inputs),
    diagnostics = list(
      complete = TRUE,
      list_size = list_size,
      maximum_absolute_score = max(abs(r_hat))
    )
  )
  class(output) <- c("m_sis_fit", "hfsis_fit", "hfsis")
  output
}
