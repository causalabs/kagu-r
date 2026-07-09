#' Per-node model fitting
#'
#' @name inference
NULL

#' Fit a single node's conditional distribution
#'
#' Thin dispatcher that delegates to the node's [Mechanism]. Kept as a
#' convenience/entry point; `KaguModel$fit()` calls it once per node.
#'
#' @param node Character scalar - name of the node to fit.
#' @param parents Character vector of parent node names (empty for root nodes).
#' @param data A `data.frame` with columns for all nodes.
#' @param mechanism A `Mechanism` instance (defaults to [GPMechanism] in
#'   `KaguModel`).
#' @param ... Additional arguments forwarded to the mechanism's `$fit()`.
#' @return A mechanism-specific fit object.
#' @export
fit_node <- function(node, parents, data, mechanism, ...) {
  mechanism$fit(node, parents, data, ...)
}

#' Fit a node and return both its fit and its log marginal likelihood
#'
#' Used by [kagu_discover()]: the fit is reused for effect propagation, the log
#' marginal likelihood is the per-node factor in a DAG's marginal likelihood. A
#' failed fit or marginal likelihood yields `NA` (and `NULL` fit) with a warning,
#' so a single problematic node does not abort an entire search.
#'
#' @inheritParams fit_node
#' @return A list with elements `fit` and `logml`.
#' @noRd
.fit_and_marglik <- function(node, parents, data, mechanism, ...) {
  label <- sprintf("%s ~ %s", node,
                   if (length(parents)) paste(parents, collapse = " + ") else "1")

  fit <- tryCatch(mechanism$fit(node, parents, data, ...),
                  error = function(e) e)
  if (inherits(fit, "error")) {
    warning(sprintf("Fit failed for %s: %s", label, conditionMessage(fit)))
    return(list(fit = NULL, logml = NA_real_))
  }

  # The fit already carries its log marginal likelihood; fall back to computing
  # it if a mechanism does not.
  logml <- if (!is.null(fit$logml)) fit$logml else {
    tryCatch(mechanism$log_marglik(node, parents, data),
             error = function(e) {
               warning(sprintf("Marginal likelihood failed for %s: %s",
                               label, conditionMessage(e)))
               NA_real_
             })
  }
  list(fit = fit, logml = logml)
}
