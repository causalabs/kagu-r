#' Mechanism base class and the Gaussian-process implementation
#'
#' A `Mechanism` encapsulates the probabilistic model for a single node's
#' conditional distribution `P(X | parents)`. Kagu's sole built-in mechanism is
#' the Gaussian process ([GPMechanism]), a flexible nonparametric default. To
#' add a new backend, subclass `Mechanism` and implement the methods below.
#'
#' @name mechanisms
NULL

# =============================================================================
# Mechanism ABC
# =============================================================================

#' Mechanism abstract base class
#'
#' @description
#' A mechanism owns the full lifecycle of a node's conditional model, decoupled
#' from any particular inference backend. Subclasses implement:
#' - `$fit(node, parents, data, ...)` - fit the local model, returning an opaque
#'   fit object.
#' - `$predict_mean(node, parents, parent_values, fit)` - the conditional mean
#'   `E[node | parents]` for each posterior draw, as a `[n_chains, n_draws]`
#'   matrix (Gaussian processes use a single chain, so `n_chains = 1`).
#' - `$log_marglik(node, parents, data)` - the log marginal likelihood
#'   `log P(node | parents)`, used by structure discovery.
#' - `$posterior_shape(fit)` - `c(n_chains, n_draws)` for the fit.
#' - `$node_terms(node, parents, data, fit)` - summary rows for `model$summary()`.
#' - `$simulate_prior(node, parents, parent_values, data, ndraws)` - replicated
#'   node values drawn from the prior predictive distribution.
#' - `$simulate_posterior(node, parents, parent_values, fit, draws, n)` -
#'   replicated node values drawn from the posterior predictive distribution.
#' - `$function_draws(node, parents, grid, data, fit, ndraws)` - prior or
#'   posterior draws of the node's conditional-mean function on a grid.
#'
#' @export
Mechanism <- R6::R6Class("Mechanism",
  public = list(
    #' @description Fit the node's local model.
    #' @param node Character scalar - the node name (response).
    #' @param parents Character vector of parent node names.
    #' @param data A `data.frame` with columns for the node and its parents.
    #' @param ... Backend-specific arguments.
    #' @return An opaque fit object.
    fit = function(node, parents, data, ...) {
      stop("Mechanism$fit() is abstract - implement in a subclass.")
    },

    #' @description Conditional mean for each posterior draw.
    #' @param node Character scalar - the node name.
    #' @param parents Character vector of parent node names.
    #' @param parent_values Named list of `[n_chains, n_draws]` matrices, one per
    #'   parent, giving the parent values for each posterior draw.
    #' @param fit A fit object from `$fit()`.
    #' @return A `[n_chains, n_draws]` numeric matrix.
    predict_mean = function(node, parents, parent_values, fit) {
      stop("Mechanism$predict_mean() is abstract - implement in a subclass.")
    },

    #' @description Log marginal likelihood `log P(node | parents)`.
    #' @param node Character scalar - the node name.
    #' @param parents Character vector of parent node names.
    #' @param data A `data.frame`.
    #' @return A single numeric.
    log_marglik = function(node, parents, data) {
      stop("Mechanism$log_marglik() is abstract - implement in a subclass.")
    },

    #' @description Posterior shape of a fit.
    #' @param fit A fit object from `$fit()`.
    #' @return Named integer vector `c(n_chains, n_draws)`.
    posterior_shape = function(fit) {
      stop("Mechanism$posterior_shape() is abstract - implement in a subclass.")
    },

    #' @description Summary rows for this node (used by `model$summary()`).
    #' @param node Character scalar - the node name.
    #' @param parents Character vector of parent node names.
    #' @param data A `data.frame`.
    #' @param fit A fit object from `$fit()`.
    #' @return A `tibble` with columns `node`, `term`, `mean`, `sd`,
    #'   `hdi_lower`, `hdi_upper`.
    node_terms = function(node, parents, data, fit) {
      stop("Mechanism$node_terms() is abstract - implement in a subclass.")
    },

    #' @description Simulate node values from the prior predictive distribution.
    #' @param node Character scalar - the node name.
    #' @param parents Character vector of parent node names.
    #' @param parent_values Named list of `[ndraws, n]` matrices of simulated
    #'   parent values, one per parent (empty for a root node).
    #' @param data A `data.frame` - sets the node's location/scale and `n`.
    #' @param ndraws Integer - number of prior predictive draws.
    #' @return A `[ndraws, n]` numeric matrix.
    simulate_prior = function(node, parents, parent_values, data, ndraws) {
      stop("Mechanism$simulate_prior() is abstract - implement in a subclass.")
    },

    #' @description Simulate node values from the posterior predictive
    #'   distribution (conditional mean plus observation noise).
    #' @param node Character scalar - the node name.
    #' @param parents Character vector of parent node names.
    #' @param parent_values Named list, one entry per parent: either a length-`n`
    #'   vector shared by every draw (e.g. the observed parents) or a
    #'   `[length(draws), n]` matrix with one row per draw.
    #' @param fit A fit object from `$fit()`.
    #' @param draws Integer vector - which posterior draws to use.
    #' @param n Integer - number of observations to simulate per draw.
    #' @return A `[length(draws), n]` numeric matrix.
    simulate_posterior = function(node, parents, parent_values, fit, draws, n) {
      stop("Mechanism$simulate_posterior() is abstract - implement in a subclass.")
    },

    #' @description Prior (`fit = NULL`) or posterior draws of the node's
    #'   conditional-mean function, evaluated on a grid of parent values.
    #' @param node Character scalar - the node name.
    #' @param parents Character vector of parent node names.
    #' @param grid A `data.frame` with one column per parent.
    #' @param data A `data.frame` - sets the node's location/scale.
    #' @param fit A fit object from `$fit()`, or `NULL` for prior draws.
    #' @param ndraws Integer - number of function draws.
    #' @return A `[ndraws, nrow(grid)]` numeric matrix.
    function_draws = function(node, parents, grid, data, fit = NULL, ndraws = 50L) {
      stop("Mechanism$function_draws() is abstract - implement in a subclass.")
    }
  )
)

# =============================================================================
# GPMechanism
# =============================================================================

#' Gaussian-process mechanism (default)
#'
#' @description
#' Models each node's conditional mean as an **exact Gaussian process** of its
#' parents, with a squared-exponential (ARD) kernel - one lengthscale per parent,
#' which captures non-linearity and interactions automatically. Kernel
#' hyperparameters (lengthscales, signal and noise variance) are point-estimated
#' by type-II maximum likelihood (`hyper = "ml"`) or maximum a posteriori
#' (`hyper = "map"`) under the priors in [gp_prior()]; there is no MCMC. A node
#' with no parents is modelled by its marginal (a Normal). This is Kagu's
#' default and only mechanism.
#'
#' The default [gp_prior()] is log-uniform over the box the likelihood search is
#' bounded to, so the default ML fit *is* the MAP fit under that prior. Supplying
#' any other prior switches the default to `hyper = "map"`, so the prior you
#' check with `$prior_predictive()` is the prior the fit actually uses.
#'
#' Structure discovery uses the **exact** GP log marginal likelihood, which - in
#' contrast to fast basis/eigenfunction approximations - is well calibrated: for
#' a genuinely unidentifiable (e.g. linear-Gaussian) edge it does not manufacture
#' spurious confidence about direction.
#'
#' Posterior function samples for effect propagation are drawn by **pathwise /
#' decoupled sampling** (a random-feature prior plus the exact data update), so
#' each draw is a coherent function evaluable at any point - giving correctly
#' correlated uncertainty for do-calculus contrasts.
#'
#' @export
#'
#' @examples
#' mech <- GPMechanism$new()
#'
#' # Weakly informative priors, fitted by MAP
#' mech <- GPMechanism$new(prior = gp_prior(lengthscale = prior_lognormal(0, 0.5)))
#' mech$hyper
GPMechanism <- R6::R6Class("GPMechanism",
  inherit = Mechanism,
  public = list(
    #' @field num_results Integer - number of posterior draws to keep (default 1000).
    num_results = 1000L,
    #' @field n_features Integer - random Fourier features for pathwise sampling
    #'   (default 300).
    n_features = 300L,
    #' @field jitter Numeric - diagonal jitter for numerical stability (default 1e-6).
    jitter = 1e-6,
    #' @field prior A [gp_prior()] - priors on the (standardised) hyperparameters.
    prior = NULL,
    #' @field hyper Character - `"ml"` (type-II maximum likelihood) or `"map"`.
    hyper = "ml",

    #' @description Create a new GPMechanism.
    #' @param num_results Integer - number of posterior draws to keep.
    #' @param n_features Integer - number of random Fourier features.
    #' @param jitter Numeric - diagonal jitter added to the kernel.
    #' @param prior A [gp_prior()] object (default: the implicit log-uniform box).
    #' @param hyper `"ml"` or `"map"`. `NULL` (default) picks `"ml"` when every
    #'   kernel prior is log-uniform (where ML and MAP coincide) and `"map"`
    #'   otherwise.
    initialize = function(num_results = 1000L, n_features = 300L, jitter = 1e-6,
                          prior = gp_prior(), hyper = NULL) {
      if (!inherits(prior, "kagu_gp_prior")) {
        stop("`prior` must be created with gp_prior().")
      }
      self$num_results <- as.integer(num_results)
      self$n_features  <- as.integer(n_features)
      self$jitter      <- jitter
      self$prior       <- prior
      self$hyper <- if (is.null(hyper)) {
        all_loguniform <- all(vapply(prior[c("lengthscale", "signal_sd", "noise_sd")],
                                     function(p) p$family == "loguniform", logical(1)))
        if (all_loguniform) "ml" else "map"
      } else {
        match.arg(hyper, c("ml", "map"))
      }
    },

    #' @description Fit the node (exact GP if it has parents, marginal Normal if not).
    fit = function(node, parents, data, ...) {
      y <- data[[node]]

      if (length(parents) == 0L) {
        # Root node: marginal Normal with a flat prior on the mean.
        n   <- length(y)
        mu  <- stats::rnorm(self$num_results, mean(y), stats::sd(y) / sqrt(n))
        sig <- sqrt((n - 1) * stats::var(y) / stats::rchisq(self$num_results, n - 1))
        return(list(kind = "root", n_draws = self$num_results,
                    mu_draws = mu, sigma_draws = sig, logml = .gaussian_evidence(y)))
      }

      fit <- .gp_core(y, as.matrix(data[parents]), self$jitter,
                      prior = self$prior, hyper = self$hyper)
      fit <- .gp_add_sampler(fit, self$num_results, self$n_features)
      fit$parents <- parents
      fit
    },

    #' @description Conditional mean draws (see [Mechanism]).
    predict_mean = function(node, parents, parent_values, fit) {
      if (identical(fit$kind, "root")) {
        return(matrix(fit$mu_draws, nrow = 1L))
      }
      # One query point per posterior draw; evaluate the matched pathwise sample.
      Q <- do.call(cbind, lapply(parents, function(p) as.vector(parent_values[[p]])))
      matrix(.gp_eval(fit, Q, matched = TRUE), nrow = 1L)
    },

    #' @description Exact GP log marginal likelihood (see [Mechanism]).
    log_marglik = function(node, parents, data) {
      y <- data[[node]]
      if (length(parents) == 0L) return(.gaussian_evidence(y))
      .gp_core(y, as.matrix(data[parents]), self$jitter,
               prior = self$prior, hyper = self$hyper)$logml
    },

    #' @description Posterior shape (single chain).
    posterior_shape = function(fit) {
      c(n_chains = 1L, n_draws = fit$n_draws)
    },

    #' @description Per-node summary rows: each parent's direct local effect
    #'   (function gradient at the parent means) plus the residual noise sd.
    node_terms = function(node, parents, data, fit) {
      sigma_draws <- if (identical(fit$kind, "root")) {
        fit$sigma_draws
      } else {
        s2 <- fit$sn2 * fit$sy^2                       # noise variance, original scale
        sqrt(fit$n * s2 / stats::rchisq(fit$n_draws, fit$n))
      }
      rows <- list(.term_row(node, "sigma (noise)", sigma_draws))

      if (length(parents) > 0L) {
        means <- vapply(parents, function(p) mean(data[[p]]), numeric(1))
        eps   <- vapply(parents, function(p) stats::sd(data[[p]]) * 1e-3, numeric(1))
        base  <- matrix(means, nrow = 1L)
        f0    <- as.vector(.gp_eval(fit, base, matched = FALSE))
        for (j in seq_along(parents)) {
          pert <- base; pert[, j] <- pert[, j] + eps[[j]]
          f1   <- as.vector(.gp_eval(fit, pert, matched = FALSE))
          rows[[length(rows) + 1L]] <- .term_row(node, parents[[j]], (f1 - f0) / eps[[j]])
        }
      }
      do.call(rbind, rows)
    },

    #' @description Prior predictive draws (see [Mechanism]). Each draw samples
    #'   hyperparameters from `prior`, an exact GP function at the parents' rows,
    #'   and observation noise - on the standardised scale - then maps back to the
    #'   data's location and scale.
    simulate_prior = function(node, parents, parent_values, data, ndraws) {
      n  <- nrow(data)
      sc <- .y_scale(data[[node]])
      p  <- self$prior

      if (length(parents) == 0L) {
        mu <- .prior_sample(p$root_mean, ndraws)
        sd <- .prior_sample(p$root_sd, ndraws)
        z  <- matrix(stats::rnorm(ndraws * n), ndraws, n)
        return(sc$my + sc$sy * (mu + sd * z))   # mu, sd recycle down the rows
      }

      xs  <- .x_scale(as.matrix(data[parents]))
      out <- matrix(NA_real_, ndraws, n)
      for (s in seq_len(ndraws)) {
        X  <- do.call(cbind, lapply(parents, function(q) parent_values[[q]][s, ]))
        f  <- .rgp_prior(.standardise(X, xs), p, self$jitter)
        sn <- .prior_sample(p$noise_sd, 1L)
        out[s, ] <- sc$my + sc$sy * (f + sn * stats::rnorm(n))
      }
      out
    },

    #' @description Posterior predictive draws (see [Mechanism]): the matched
    #'   pathwise function draw at each row's parents plus noise, with the noise
    #'   sd drawn as in `$node_terms()`.
    simulate_posterior = function(node, parents, parent_values, fit, draws, n) {
      k <- length(draws)
      if (identical(fit$kind, "root")) {
        z <- matrix(stats::rnorm(k * n), k, n)
        return(fit$mu_draws[draws] + fit$sigma_draws[draws] * z)
      }

      sig <- fit$sy * sqrt(fit$n * fit$sn2 / stats::rchisq(k, fit$n))
      shared <- !is.matrix(parent_values[[parents[[1]]]])
      if (shared) {
        # Same parents for every draw (e.g. the observed data): evaluate once.
        X <- do.call(cbind, lapply(parents, function(q) parent_values[[q]]))
        f <- t(.gp_eval(fit, X, matched = FALSE)[, draws, drop = FALSE])   # [k x n]
      } else {
        f <- t(vapply(seq_len(k), function(i) {
          X <- do.call(cbind, lapply(parents, function(q) parent_values[[q]][i, ]))
          .gp_eval_draw(fit, X, draws[[i]])
        }, numeric(n)))
      }
      f + sig * matrix(stats::rnorm(k * n), k, n)   # sig recycles down the rows
    },

    #' @description Prior or posterior function draws on a grid (see [Mechanism]).
    function_draws = function(node, parents, grid, data, fit = NULL, ndraws = 50L) {
      if (length(parents) == 0L) {
        stop(sprintf("Node '%s' has no parents, so it has no function to draw.", node))
      }
      G <- as.matrix(grid[parents])
      if (is.null(fit)) {
        sc <- .y_scale(data[[node]])
        Gs <- .standardise(G, .x_scale(as.matrix(data[parents])))
        return(t(vapply(seq_len(ndraws), function(s) {
          sc$my + sc$sy * .rgp_prior(Gs, self$prior, self$jitter)
        }, numeric(nrow(G)))))
      }
      draws <- sort(sample.int(fit$n_draws, ndraws))
      t(.gp_eval(fit, G, matched = FALSE)[, draws, drop = FALSE])
    }
  )
)

# =============================================================================
# Internal exact-GP helpers
# =============================================================================

#' Squared-exponential (ARD) kernel on standardised inputs
#' @noRd
.rbf <- function(X1, X2, ls, sf2) {
  Z1 <- sweep(X1, 2, ls, "/"); Z2 <- sweep(X2, 2, ls, "/")
  d2 <- outer(rowSums(Z1^2), rowSums(Z2^2), "+") - 2 * tcrossprod(Z1, Z2)
  sf2 * exp(-0.5 * pmax(d2, 0))
}

#' Fit an exact GP by type-II ML or MAP; return hypers, Cholesky, log evidence
#'
#' Hyperparameters are optimised on `p = (log ls_1..d, log sf2, log sn2)`. The
#' search box is each log-uniform prior's support (the default [gp_prior()]
#' reproduces the historical box exactly) or a fixed default box otherwise. With
#' `hyper = "map"` the log prior density of `p` is added to the objective. The
#' returned `logml` is always the exact log marginal likelihood at the optimum.
#' @noRd
.gp_core <- function(y, X, jitter = 1e-6, prior = gp_prior(), hyper = "ml") {
  X <- as.matrix(X); n <- nrow(X); d <- ncol(X)
  xs <- .x_scale(X); Xs <- .standardise(X, xs)
  sc <- .y_scale(y); ys <- (y - sc$my) / sc$sy

  nll <- function(p) {
    ls <- exp(p[seq_len(d)]); sf2 <- exp(p[d + 1]); sn2 <- exp(p[d + 2])
    K  <- .rbf(Xs, Xs, ls, sf2) + diag(sn2 + jitter, n)
    ch <- tryCatch(chol(K), error = function(e) NULL); if (is.null(ch)) return(1e10)
    al <- backsolve(ch, forwardsolve(t(ch), ys))
    0.5 * sum(ys * al) + sum(log(diag(ch))) + 0.5 * n * log(2 * pi)
  }

  # sd priors map onto variances via log(sf2) = 2 * log(sf)
  b_ls  <- .or_default(.prior_log_bounds(prior$lengthscale), c(-3, 4))
  b_sf2 <- 2 * .or_default(.prior_log_bounds(prior$signal_sd), c(-2, 2))
  b_sn2 <- 2 * .or_default(.prior_log_bounds(prior$noise_sd), c(-3, 1))
  lower <- c(rep(b_ls[1], d), b_sf2[1], b_sn2[1])
  upper <- c(rep(b_ls[2], d), b_sf2[2], b_sn2[2])
  start <- pmin(pmax(c(rep(0, d), 0, -1), lower), upper)

  objective <- if (identical(hyper, "map")) {
    function(p) {
      lp <- sum(.prior_logdens_log(prior$lengthscale, p[seq_len(d)])) +
        .prior_logdens_log(prior$signal_sd, p[d + 1] / 2) +
        .prior_logdens_log(prior$noise_sd,  p[d + 2] / 2)
      if (!is.finite(lp)) return(1e10)
      nll(p) - lp
    }
  } else {
    nll
  }

  opt <- stats::optim(start, objective, method = "L-BFGS-B",
                      lower = lower, upper = upper)
  ls <- exp(opt$par[seq_len(d)]); sf2 <- exp(opt$par[d + 1]); sn2 <- exp(opt$par[d + 2])
  K  <- .rbf(Xs, Xs, ls, sf2) + diag(sn2 + jitter, n)
  ch <- chol(K); alpha <- backsolve(ch, forwardsolve(t(ch), ys))

  list(kind = "gp", Xs = Xs, ys = ys, ch = ch, alpha = alpha,
       ls = ls, sf2 = sf2, sn2 = sn2, mx = xs$mx, sx = xs$sx, my = sc$my, sy = sc$sy,
       n = n, d = d, hyper = hyper, logml = -nll(opt$par))
}

#' Column location/scale used to standardise GP inputs
#' @noRd
.x_scale <- function(X) {
  mx <- colMeans(X); sx <- apply(X, 2, stats::sd)
  sx[!is.finite(sx) | sx == 0] <- 1
  list(mx = mx, sx = sx)
}

#' Location/scale used to standardise a GP response
#' @noRd
.y_scale <- function(y) {
  my <- mean(y); sy <- stats::sd(y)
  if (!is.finite(sy) || sy == 0) sy <- 1
  list(my = my, sy = sy)
}

#' @noRd
.standardise <- function(X, xs) {
  sweep(sweep(as.matrix(X), 2, xs$mx, "-"), 2, xs$sx, "/")
}

#' @noRd
.or_default <- function(x, default) if (is.null(x)) default else x

#' One exact GP prior function draw at standardised inputs (standardised scale)
#'
#' Samples hyperparameters from `prior`, then `f ~ N(0, K)`. Duplicate input
#' rows (e.g. a discrete parent) share one function value, which is both exact
#' and avoids a singular kernel.
#' @noRd
.rgp_prior <- function(Xs, prior, jitter = 1e-6) {
  Xs  <- as.matrix(Xs)
  key <- do.call(paste, c(as.data.frame(Xs), sep = "\r"))
  u   <- !duplicated(key)
  U   <- Xs[u, , drop = FALSE]
  ls  <- .prior_sample(prior$lengthscale, ncol(Xs))
  sf2 <- .prior_sample(prior$signal_sd, 1L)^2
  R   <- .safe_chol(.rbf(U, U, ls, sf2), jitter * sf2)
  f_u <- as.vector(crossprod(R, stats::rnorm(nrow(U))))
  f_u[match(key, key[u])]
}

#' Cholesky factor with escalating diagonal jitter
#' @noRd
.safe_chol <- function(K, eps) {
  eps <- max(eps, 1e-10)
  for (i in seq_len(8)) {
    R <- tryCatch(chol(K + diag(eps, nrow(K))), error = function(e) NULL)
    if (!is.null(R)) return(R)
    eps <- eps * 10
  }
  stop("Kernel matrix is not positive definite, even with added jitter.")
}

#' Attach pathwise (decoupled) sampling state to an exact-GP fit
#'
#' Random-feature prior + exact data update (Matheron's rule), giving `S`
#' coherent posterior function samples evaluable at any query point.
#' @noRd
.gp_add_sampler <- function(fit, S, R) {
  n <- fit$n; d <- fit$d
  Om <- matrix(stats::rnorm(R * d), R, d) / rep(fit$ls, each = R)   # spectral freqs
  bb <- stats::runif(R, 0, 2 * pi)
  W  <- matrix(stats::rnorm(R * S), R, S)                          # prior weights
  phi_tr <- sqrt(2 * fit$sf2 / R) * cos(fit$Xs %*% t(Om) + rep(bb, each = n))
  eps    <- matrix(stats::rnorm(n * S, sd = sqrt(fit$sn2)), n, S)
  resid  <- fit$ys - phi_tr %*% W - eps
  fit$V  <- backsolve(fit$ch, forwardsolve(t(fit$ch), resid))       # K^{-1}(y - prior)
  fit$Om <- Om; fit$bb <- bb; fit$W <- W; fit$R <- R; fit$n_draws <- S
  fit
}

#' Evaluate the pathwise samples at query points (original scale)
#'
#' `matched = TRUE`: `Xnew` has one row per draw; returns the length-S vector
#' `f^(s)(Xnew[s, ])` (coherent per-draw propagation). `matched = FALSE`: returns
#' the full `[nrow(Xnew) x S]` matrix (every draw at every point).
#' @noRd
.gp_eval <- function(fit, Xnew, matched = FALSE) {
  Xn  <- sweep(sweep(as.matrix(Xnew), 2, fit$mx, "-"), 2, fit$sx, "/")
  Phi <- sqrt(2 * fit$sf2 / fit$R) * cos(Xn %*% t(fit$Om) + rep(fit$bb, each = nrow(Xn)))
  Ks  <- .rbf(Xn, fit$Xs, fit$ls, fit$sf2)
  if (matched) {
    val <- rowSums(Phi * t(fit$W)) + rowSums(Ks * t(fit$V))          # [S]
  } else {
    val <- Phi %*% fit$W + Ks %*% fit$V                              # [m x S]
  }
  fit$my + fit$sy * val
}

#' Evaluate a single pathwise draw `s` at many query points (original scale)
#'
#' Used for joint (whole-graph) posterior predictive simulation, where each draw
#' has its own simulated parent values.
#' @noRd
.gp_eval_draw <- function(fit, Xnew, s) {
  Xn  <- .standardise(Xnew, list(mx = fit$mx, sx = fit$sx))
  Phi <- sqrt(2 * fit$sf2 / fit$R) * cos(Xn %*% t(fit$Om) + rep(fit$bb, each = nrow(Xn)))
  Ks  <- .rbf(Xn, fit$Xs, fit$ls, fit$sf2)
  as.vector(fit$my + fit$sy * (Phi %*% fit$W[, s] + Ks %*% fit$V[, s]))
}

#' Gaussian marginal evidence of a root node
#'
#' Computed on the variable standardised to unit variance (as are the conditional
#' evidences), so discovery scores are scale-free and comparable across nodes.
#' @noRd
.gaussian_evidence <- function(y) {
  n <- length(y)
  -0.5 * n * (log(2 * pi) + 1)
}

#' Summarise a vector of posterior draws into one summary-table row
#' @noRd
.term_row <- function(node, term, draws) {
  h <- bayestestR::hdi(as.numeric(draws), ci = 0.90)
  tibble::tibble(
    node = node, term = term,
    mean = mean(draws), sd = stats::sd(draws),
    hdi_lower = h$CI_low, hdi_upper = h$CI_high
  )
}
