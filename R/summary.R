#' Summary table for a fitted KaguModel
#'
#' @name summary
NULL

#' Build a per-node summary table across all nodes
#'
#' For each node, the node's mechanism contributes summary rows - for a Gaussian
#' process, the direct local effect of each parent (the fitted function's
#' gradient at the parents' means, comparable to a regression coefficient) and
#' the residual noise sd.
#'
#' @param model A fitted `KaguModel`.
#' @return A `tibble` with columns `node`, `term`, `mean`, `sd`, `hdi_lower`,
#'   `hdi_upper`.
#' @export
build_summary_table <- function(model) {
  order   <- topological_sort(model$dag)
  results <- lapply(order, function(node) {
    model$mechanisms[[node]]$node_terms(
      node, model$dag[[node]], model$data, model$traces[[node]]
    )
  })
  tibble::as_tibble(do.call(rbind, results))
}
