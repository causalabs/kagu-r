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
#' @param save_all Logical — if `TRUE`, fit with `save_pars = save_pars(all =
#'   TRUE)` so the fit is suitable for bridge sampling (marginal likelihood).
#'   Required by [.node_log_marglik()]; off by default.
#' @param ... Additional arguments forwarded to `brms::brm()`.
#' @return A `brmsfit` object.
#' @export
fit_node <- function(node, parents, data, mechanism,
                     draws = 1000L, tune = 1000L, chains = 4L,
                     backend = "cmdstanr", silent = 2L, save_all = FALSE, ...) {
  args <- c(
    list(
      formula = mechanism$build_formula(node, parents),
      data    = data,
      family  = mechanism$family(),
      prior   = .build_priors(mechanism, parents),
      chains  = chains,
      iter    = draws + tune,
      warmup  = tune,
      backend = backend,
      silent  = silent,
      refresh = 0
    ),
    list(...)
  )
  if (save_all) args$save_pars <- brms::save_pars(all = TRUE)

  .quietly(do.call(brms::brm, args))
}

#' Fit a node's local model and bridge-sample its log marginal likelihood
#'
#' Fits `node ~ parents` (with all parameters saved) and bridge-samples the log
#' marginal likelihood `log P(node | parents)`. Returns *both* the fit and the
#' log marginal likelihood: the fit is reused for effect propagation (so
#' [kagu_discover()] never has to refit), while the log marginal likelihood is
#' the per-node factor in the Markov factorisation of a DAG's marginal
#' likelihood. A failed fit or bridge sample yields `NA` (and `NULL` fit) with a
#' warning, so a single problematic node does not abort an entire search.
#'
#' @inheritParams fit_node
#' @return A list with elements `fit` (a `brmsfit` or `NULL`) and `logml`
#'   (numeric, possibly `NA`).
#' @noRd
.fit_and_marglik <- function(node, parents, data, mechanism,
                             draws = 2000L, tune = 1000L, chains = 4L,
                             backend = "cmdstanr", ...) {
  label <- sprintf("%s ~ %s", node,
                   if (length(parents)) paste(parents, collapse = " + ") else "1")

  fit <- tryCatch(
    fit_node(node, parents, data, mechanism, draws = draws, tune = tune,
             chains = chains, backend = backend, save_all = TRUE, ...),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    warning(sprintf("Fit failed for %s: %s", label, conditionMessage(fit)))
    return(list(fit = NULL, logml = NA_real_))
  }

  bs <- tryCatch(
    .quietly(bridgesampling::bridge_sampler(fit)),
    error = function(e) e
  )
  if (inherits(bs, "error")) {
    warning(sprintf("Bridge sampling failed for %s: %s", label,
                    conditionMessage(bs)))
    return(list(fit = fit, logml = NA_real_))
  }
  list(fit = fit, logml = bs$logml)
}

#' Log marginal likelihood of a single node's local model
#'
#' Thin wrapper around [.fit_and_marglik()] that returns only the log marginal
#' likelihood (or `NA`).
#'
#' @inheritParams fit_node
#' @return A single numeric.
#' @noRd
.node_log_marglik <- function(node, parents, data, mechanism, ...) {
  .fit_and_marglik(node, parents, data, mechanism, ...)$logml
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
