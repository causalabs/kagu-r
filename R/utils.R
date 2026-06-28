# Internal utilities — not exported

#' Extract (n_chains, n_draws) from a brmsfit
#'
#' @param fit A `brmsfit` object.
#' @return Named integer vector with elements `n_chains` and `n_draws`.
#' @noRd
.chain_draw_shape <- function(fit) {
  draws <- posterior::as_draws_array(fit)
  c(n_chains = dim(draws)[[2]], n_draws = dim(draws)[[1]])
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

#' Suppress brms/Stan console output during sampling
#'
#' @param expr Expression to evaluate silently.
#' @noRd
.quietly <- function(expr) {
  suppressMessages(suppressWarnings(
    utils::capture.output(result <- expr, type = "message")
  ))
  result
}
