#' Per-node model fitting via brms
#'
#' @name inference
NULL

#' Fit a single node's conditional distribution
#'
#' Builds the brms formula from the mechanism, fits with `brms::brm()`, and
#' returns the resulting `brmsfit`. Stan/brms console output is suppressed;
#' on error the captured output is printed before re-raising.
#'
#' @param node Character scalar — name of the node to fit.
#' @param parents Character vector of parent node names (empty for root nodes).
#' @param data A `data.frame` with columns for all nodes.
#' @param mechanism A `Mechanism` instance.
#' @param draws Integer — number of post-warmup draws per chain (default 1000).
#' @param tune Integer — number of warmup draws per chain (default 1000).
#' @param chains Integer — number of MCMC chains (default 4).
#' @param backend Character — Stan backend passed to `brms::brm()`. Defaults
#'   to `"cmdstanr"` for faster compilation and sampling. Set to `"rstan"` if
#'   cmdstanr is not installed.
#' @param silent Integer passed to `brms::brm()` — 2 suppresses all output
#'   (default), 0 shows everything.
#' @param ... Additional arguments forwarded to `brms::brm()`.
#' @return A `brmsfit` object.
#' @export
fit_node <- function(node, parents, data, mechanism,
                     draws = 1000L, tune = 1000L, chains = 4L,
                     backend = "cmdstanr", silent = 2L, ...) {
  formula <- mechanism$build_formula(node, parents)
  priors  <- .build_priors(mechanism, parents)

  .quietly(brms::brm(
    formula  = formula,
    data     = data,
    family   = mechanism$family(),
    prior    = priors,
    chains   = chains,
    iter     = draws + tune,
    warmup   = tune,
    backend  = backend,
    silent   = silent,
    refresh  = 0,
    ...
  ))
}

# --- Internal helpers --------------------------------------------------------

#' Build brms priors from a LinearMechanism (or fallback defaults)
#' @noRd
.build_priors <- function(mechanism, parents) {
  if (!inherits(mechanism, "LinearMechanism")) {
    return(brms::empty_prior())
  }

  p <- c(
    brms::prior_string(
      sprintf("normal(0, %g)", mechanism$prior_alpha),
      class = "Intercept"
    ),
    brms::prior_string(
      sprintf("normal(0, %g)", mechanism$prior_sigma),
      class = "sigma"
    )
  )

  for (parent in parents) {
    p <- c(p, brms::prior_string(
      sprintf("normal(0, %g)", mechanism$prior_beta),
      class = "b", coef = parent
    ))
  }

  p
}
