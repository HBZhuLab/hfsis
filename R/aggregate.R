.validate_fit_list <- function(fits) {
  if (!is.list(fits) || length(fits) < 1L ||
      !all(vapply(fits, inherits, logical(1L), what = "hfsis_fit"))) {
    .hfsis_stop("'fits' must be a nonempty list of hfsis_fit objects.")
  }
  variables <- fits[[1L]]$inputs_summary$variable_names
  if (!all(vapply(fits, function(fit) {
    identical(fit$inputs_summary$variable_names, variables)
  }, logical(1L)))) {
    .hfsis_stop("Every fit must use the same variable names and order.")
  }
  variables
}

.validate_metadata <- function(metadata, variables) {
  if (!is.data.frame(metadata) || !"variable" %in% names(metadata) ||
      anyNA(metadata$variable) || anyDuplicated(metadata$variable)) {
    .hfsis_stop("'metadata' must have a unique, nonmissing 'variable' column.")
  }
  locations <- match(variables, metadata$variable)
  if (anyNA(locations)) {
    .hfsis_stop("Every fitted variable must occur in 'metadata$variable'.")
  }
  metadata[locations, , drop = FALSE]
}

#' Selection frequencies across repeated HF-SIS fits
#'
#' Frequencies describe empirical recurrence only. They are neither support
#' probabilities nor significance measures.
#'
#' @param fits A nonempty list of fitted hfsis_fit objects with a common
#'   variable universe and order.
#' @param top Optional positive number of variables retained per method after
#'   ordering by decreasing frequency, mean rank, and original position.
#' @param metadata Optional data frame whose unique variable column matches
#'   variables to annotations. It is joined without changing variable order.
#'
#' @return A data frame of class hfsis_selection_frequency with columns
#'   variable, method, selected_count, fit_count, frequency, mean_rank, and
#'   median_rank, followed by optional metadata.
#' @references Zhu, H. *Ultra-high-dimensional spot support screening in
#'   continuous-time regression*. Manuscript.
#' @seealso [cluster_share()], [plot_selection_frequency()]
#' @examples
#' data(hf_sis_panel_example)
#' fits <- lapply(hf_sis_panel_example$sessions[1:2], function(dat) {
#'   m_sis(as.matrix(dat[hf_sis_panel_example$x_names]), dat$y,
#'     hf_sis_panel_example$target_index, list_size = 5,
#'     control = hfsis_control(hf_sis_panel_example$windows))
#' })
#' selection_frequency(fits)
#' @export
selection_frequency <- function(fits, top = NULL, metadata = NULL) {
  variables <- .validate_fit_list(fits)
  methods <- unique(vapply(fits, function(fit) fit$method, character(1L)))
  result <- do.call(rbind, lapply(methods, function(method) {
    method_fits <- fits[vapply(fits, function(fit) {
      identical(fit$method, method)
    }, logical(1L))]
    selected_count <- vapply(variables, function(variable) {
      sum(vapply(method_fits, function(fit) {
        variable %in% fit$selected
      }, logical(1L)))
    }, integer(1L))
    ranks <- vapply(method_fits, function(fit) {
      match(variables, fit$ranking)
    }, integer(length(variables)))
    if (is.null(dim(ranks))) ranks <- matrix(ranks, ncol = 1L)
    mean_rank <- apply(ranks, 1L, function(value) {
      if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE)
    })
    median_rank <- apply(ranks, 1L, function(value) {
      if (all(is.na(value))) NA_real_ else stats::median(value, na.rm = TRUE)
    })
    data.frame(
      variable = variables,
      method = method,
      selected_count = selected_count,
      fit_count = length(method_fits),
      frequency = selected_count / length(method_fits),
      mean_rank = mean_rank,
      median_rank = median_rank,
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL
  if (!is.null(metadata)) {
    ordered_metadata <- .validate_metadata(metadata, variables)
    locations <- match(result$variable, ordered_metadata$variable)
    additions <- ordered_metadata[locations,
      setdiff(names(ordered_metadata), "variable"), drop = FALSE]
    rownames(additions) <- NULL
    result <- cbind(result, additions)
  }
  if (!is.null(top)) {
    if (!.is_count(top)) .hfsis_stop("'top' must be NULL or one positive integer.")
    keep <- unlist(lapply(methods, function(method) {
      rows <- which(result$method == method)
      ordered <- order(-result$frequency[rows], result$mean_rank[rows], rows,
        na.last = TRUE)
      rows[utils::head(ordered, top)]
    }), use.names = FALSE)
    result <- result[keep, , drop = FALSE]
    rownames(result) <- NULL
  }
  class(result) <- c("hfsis_selection_frequency", "data.frame")
  result
}

#' Plot repeated-sample selection frequencies
#'
#' @param x A result from [selection_frequency()].
#' @param top Positive number of variables shown within each method.
#' @param group Optional metadata column used in variable labels.
#' @param methods Optional character subset of methods to display.
#' @param ... Additional arguments passed to [graphics::barplot()].
#'
#' @return x, invisibly.
#' @export
plot_selection_frequency <- function(
    x, top = 15L, group = NULL, methods = NULL, ...) {
  if (!inherits(x, "hfsis_selection_frequency")) {
    .hfsis_stop("'x' must be returned by selection_frequency().")
  }
  if (!.is_count(top)) .hfsis_stop("'top' must be one positive integer.")
  available_methods <- unique(x$method)
  if (is.null(methods)) methods <- available_methods
  if (!is.character(methods) || !length(methods) ||
      any(!methods %in% available_methods)) {
    .hfsis_stop("'methods' must select methods present in 'x'.")
  }
  if (!is.null(group) && (!is.character(group) || length(group) != 1L ||
      !group %in% names(x))) {
    .hfsis_stop("'group' must be NULL or name one metadata column in 'x'.")
  }
  old <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old), add = TRUE)
  graphics::par(mfrow = c(length(methods), 1L), mar = c(8, 4, 3, 1))
  for (method in methods) {
    rows <- which(x$method == method)
    rows <- rows[order(-x$frequency[rows], x$mean_rank[rows], rows,
      na.last = TRUE)]
    rows <- utils::head(rows, top)
    labels <- x$variable[rows]
    if (!is.null(group)) labels <- paste0(labels, " [", x[[group]][rows], "]")
    graphics::barplot(
      x$frequency[rows], names.arg = labels, las = 2,
      ylim = c(0, 1), ylab = "Descriptive selection frequency",
      main = method, ...
    )
  }
  invisible(x)
}

#' Aggregate repeated screening selections by cluster
#'
#' @param fits A nonempty list of fitted hfsis_fit objects.
#' @param metadata Data frame with unique variable names and a cluster column.
#' @param cluster_column Name of the metadata column defining clusters.
#'
#' @return A data frame of class hfsis_cluster_share. Shares describe the
#'   composition of selections and are not inferential probabilities.
#' @references Zhu, H. *Ultra-high-dimensional spot support screening in
#'   continuous-time regression*. Manuscript.
#' @seealso [selection_frequency()], [plot_cluster_share()]
#' @export
cluster_share <- function(fits, metadata, cluster_column = "cluster") {
  variables <- .validate_fit_list(fits)
  ordered_metadata <- .validate_metadata(metadata, variables)
  if (!is.character(cluster_column) || length(cluster_column) != 1L ||
      !cluster_column %in% names(ordered_metadata)) {
    .hfsis_stop("'cluster_column' must name one column in 'metadata'.")
  }
  if (anyNA(ordered_metadata[[cluster_column]])) {
    .hfsis_stop("Cluster labels must not be missing.")
  }
  frequency <- selection_frequency(fits, metadata = ordered_metadata)
  key <- interaction(
    frequency$method, frequency[[cluster_column]], drop = TRUE, lex.order = TRUE
  )
  pieces <- split(seq_len(nrow(frequency)), key)
  result <- do.call(rbind, lapply(pieces, function(rows) {
    method <- frequency$method[rows[1L]]
    method_total <- sum(frequency$selected_count[frequency$method == method])
    selection_count <- sum(frequency$selected_count[rows])
    data.frame(
      method = method,
      cluster = as.character(frequency[[cluster_column]][rows[1L]]),
      variable_count = length(rows),
      selection_count = selection_count,
      fit_count = unique(frequency$fit_count[rows]),
      mean_frequency = mean(frequency$frequency[rows]),
      selection_share = if (method_total > 0) selection_count / method_total else 0,
      stringsAsFactors = FALSE
    )
  }))
  rownames(result) <- NULL
  class(result) <- c("hfsis_cluster_share", "data.frame")
  result
}

#' Plot descriptive cluster selection shares
#'
#' @param x A result from [cluster_share()].
#' @param ... Additional arguments passed to [graphics::barplot()].
#'
#' @return x, invisibly.
#' @export
plot_cluster_share <- function(x, ...) {
  if (!inherits(x, "hfsis_cluster_share")) {
    .hfsis_stop("'x' must be returned by cluster_share().")
  }
  methods <- unique(x$method)
  clusters <- unique(x$cluster)
  values <- matrix(0, nrow = length(clusters), ncol = length(methods),
    dimnames = list(clusters, methods))
  for (row in seq_len(nrow(x))) {
    values[x$cluster[row], x$method[row]] <- x$selection_share[row]
  }
  graphics::barplot(
    values, beside = TRUE, legend.text = rownames(values),
    ylab = "Descriptive share of selections",
    main = "Selection composition by cluster", ...
  )
  invisible(x)
}
