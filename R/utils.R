# Internal utilities — not exported

#' Check that every DAG node has a matching column in the data
#'
#' Raises an informative error if any node (and hence any node referenced as a
#' parent) is absent from `data`. The message distinguishes a simply forgotten
#' column from a genuinely latent variable, since the latter is a substantive
#' identification problem rather than a typo.
#'
#' @param dag Named list — the DAG specification.
#' @param data A `data.frame`.
#' @noRd
.check_data_nodes <- function(dag, data) {
  missing <- setdiff(names(dag), names(data))
  if (length(missing) == 0) return(invisible(NULL))

  cli::cli_abort(c(
    "{cli::qty(missing)}DAG node{?s} {.val {missing}} {?is/are} not present in {.arg data}.",
    "i" = paste(
      "{cli::qty(missing)}If {?this variable was/these variables were} simply",
      "left out, add {?it/them} as {?a column/columns} and re-fit."
    ),
    "!" = paste(
      "{cli::qty(missing)}If {?it is/they are} {.emph latent} — an unmeasured",
      "common cause, for instance — Kagu cannot estimate {?its/their}",
      "mechanism, because there are no observations to condition on."
    ),
    "i" = paste(
      "Unobserved confounding biases every effect that flows through the",
      "missing node, and no amount of modelling recovers it from the data",
      "alone. Where a confounder is truly latent, establish that your target",
      "effect is identifiable (e.g. via an instrument, a valid adjustment set,",
      "or a front-door path) before trusting the estimates."
    )
  ))
}

#' Convert a [n_chains, n_draws] matrix to a posterior draws_array
#'
#' Produces an array of class `draws_array` with dimensions
#' [iteration, chain, variable] as expected by the `posterior` package.
#'
#' @param samples Numeric matrix with dim `[n_chains, n_draws]`.
#' @param varname Character scalar — name for the variable dimension.
#' @return A `posterior::draws_array`.
#' @noRd
.to_draws_array <- function(samples, varname = "effect") {
  n_chains <- nrow(samples)
  n_draws  <- ncol(samples)

  # posterior draws_array convention: [iteration, chain, variable]
  arr <- array(
    aperm(samples),               # [n_draws, n_chains]
    dim      = c(n_draws, n_chains, 1L),
    dimnames = list(NULL, NULL, varname)
  )
  class(arr) <- c("draws_array", "draws", "array")
  arr
}
