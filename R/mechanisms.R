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
#' - `$fit(node, parents, data, ...)` — fit the local model, returning an opaque
#'   fit object.
#' - `$predict_mean(node, parents, parent_values, fit)` — the conditional mean
#'   `E[node | parents]` for each posterior draw, as a `[n_chains, n_draws]`
#'   matrix (Gaussian processes use a single chain, so `n_chains = 1`).
#' - `$log_marglik(node, parents, data)` — the log marginal likelihood
#'   `log P(node | parents)`, used by structure discovery.
#' - `$posterior_shape(fit)` — `c(n_chains, n_draws)` for the fit.
#' - `$node_terms(node, parents, data, fit)` — summary rows for `model$summary()`.
#'
#' @export
Mechanism <- R6::R6Class("Mechanism",
  public = list(
    #' @description Fit the node's local model.
    #' @param node Character scalar — the node name (response).
    #' @param parents Character vector of parent node names.
    #' @param data A `data.frame` with columns for the node and its parents.
    #' @param ... Backend-specific arguments.
    #' @return An opaque fit object.
    fit = function(node, parents, data, ...) {
      stop("Mechanism$fit() is abstract — implement in a subclass.")
    },

    #' @description Conditional mean for each posterior draw.
    #' @param node Character scalar — the node name.
    #' @param parents Character vector of parent node names.
    #' @param parent_values Named list of `[n_chains, n_draws]` matrices, one per
    #'   parent, giving the parent values for each posterior draw.
    #' @param fit A fit object from `$fit()`.
    #' @return A `[n_chains, n_draws]` numeric matrix.
    predict_mean = function(node, parents, parent_values, fit) {
      stop("Mechanism$predict_mean() is abstract — implement in a subclass.")
    },

    #' @description Log marginal likelihood `log P(node | parents)`.
    #' @param node Character scalar — the node name.
    #' @param parents Character vector of parent node names.
    #' @param data A `data.frame`.
    #' @return A single numeric.
    log_marglik = function(node, parents, data) {
      stop("Mechanism$log_marglik() is abstract — implement in a subclass.")
    },

    #' @description Posterior shape of a fit.
    #' @param fit A fit object from `$fit()`.
    #' @return Named integer vector `c(n_chains, n_draws)`.
    posterior_shape = function(fit) {
      stop("Mechanism$posterior_shape() is abstract — implement in a subclass.")
    },

    #' @description Summary rows for this node (used by `model$summary()`).
    #' @param node Character scalar — the node name.
    #' @param parents Character vector of parent node names.
    #' @param data A `data.frame`.
    #' @param fit A fit object from `$fit()`.
    #' @return A `tibble` with columns `node`, `term`, `mean`, `sd`,
    #'   `hdi_lower`, `hdi_upper`.
    node_terms = function(node, parents, data, fit) {
      stop("Mechanism$node_terms() is abstract — implement in a subclass.")
    }
  )
)

# =============================================================================
# GPMechanism
# =============================================================================

#' Gaussian-process mechanism (default)
#'
#' @description
#' Models each node's conditional mean as a Gaussian process of its parents. The
#' GP is a Bayesian linear model over the smooth eigen-basis of the \pkg{BayesGPfit}
#' package, fitted in closed form (conjugate Normal posterior, with the noise and
#' amplitude hyperparameters set by type-II maximum likelihood). A node with no
#' parents is modelled by its marginal (a Normal). This is Kagu's default and
#' only mechanism.
#'
#' Because the posterior over the fitted function is exact and Gaussian, full
#' posterior uncertainty — not merely a credible interval — flows through effect
#' propagation, and the log marginal likelihood used by structure discovery is
#' available in closed form (no MCMC, no bridge sampling).
#'
#' @export
#'
#' @examples
#' mech <- GPMechanism$new()
GPMechanism <- R6::R6Class("GPMechanism",
  inherit = Mechanism,
  public = list(
    #' @field poly_degree Integer — GP basis polynomial degree (default 10).
    poly_degree = 10L,
    #' @field num_results Integer — number of posterior draws to keep (default 1000).
    num_results = 1000L,
    #' @field a,b Numeric — GP kernel hyperparameters (defaults 0.01, 1).
    a = 0.01,
    #' @field b b kernel hyperparameter.
    b = 1,

    #' @description Create a new GPMechanism.
    #' @param poly_degree Integer — GP basis polynomial degree.
    #' @param num_results Integer — number of posterior draws to keep.
    #' @param a,b Numeric — GP kernel hyperparameters.
    initialize = function(poly_degree = 10L, num_results = 1000L,
                          a = 0.01, b = 1) {
      self$poly_degree <- as.integer(poly_degree)
      self$num_results <- as.integer(num_results)
      self$a <- a
      self$b <- b
    },

    #' @description Fit the node (GP if it has parents, marginal Normal if not).
    fit = function(node, parents, data, ...) {
      y <- data[[node]]

      if (length(parents) == 0L) {
        # Root node: model the marginal as a Normal with a flat prior on the
        # mean, giving posterior draws of the mean and the residual sd.
        n   <- length(y)
        mu  <- stats::rnorm(self$num_results, mean(y), stats::sd(y) / sqrt(n))
        sig <- sqrt((n - 1) * stats::var(y) / stats::rchisq(self$num_results, n - 1))
        return(list(kind = "root", n_draws = self$num_results,
                    mu_draws = mu, sigma_draws = sig,
                    logml = .gaussian_evidence(y)))
      }

      pd  <- .cap_poly_degree(self$poly_degree, length(parents))
      fit <- .gp_conj_fit(y, as.matrix(data[parents]), poly_degree = pd,
                          a = self$a, b = self$b, num_results = self$num_results)
      fit$kind    <- "gp"
      fit$parents <- parents
      fit
    },

    #' @description Conditional mean draws (see [Mechanism]).
    predict_mean = function(node, parents, parent_values, fit) {
      if (identical(fit$kind, "root")) {
        return(matrix(fit$mu_draws, nrow = 1L))
      }
      n_draws <- ncol(parent_values[[parents[[1]]]])
      # One query point per posterior draw: row d = parent values at draw d.
      Q <- do.call(cbind, lapply(parents, function(p) as.vector(parent_values[[p]])))
      matrix(.gp_function_draws(fit, Q, matched = TRUE, n_draws = n_draws), nrow = 1L)
    },

    #' @description Type-II log marginal likelihood (see [Mechanism]).
    log_marglik = function(node, parents, data) {
      y <- data[[node]]
      if (length(parents) == 0L) return(.gaussian_evidence(y))
      .gp_type2_evidence(y, as.matrix(data[parents]),
                         poly_degree = .cap_poly_degree(self$poly_degree,
                                                        length(parents)),
                         a = self$a, b = self$b)
    },

    #' @description Posterior shape (single chain).
    posterior_shape = function(fit) {
      c(n_chains = 1L, n_draws = fit$n_draws)
    },

    #' @description Per-node summary rows: each parent's direct local effect
    #'   (function gradient at the parent means) plus the residual noise sd.
    node_terms = function(node, parents, data, fit) {
      rows <- list(.term_row(node, "sigma (noise)", fit$sigma_draws))

      if (length(parents) > 0L) {
        means <- vapply(parents, function(p) mean(data[[p]]), numeric(1))
        eps   <- vapply(parents, function(p) stats::sd(data[[p]]) * 1e-3, numeric(1))
        base  <- matrix(means, nrow = 1L)                       # [1 x d]
        f0    <- .gp_function_draws(fit, base, matched = FALSE)  # [1 x n_draws]
        for (j in seq_along(parents)) {
          pert      <- base; pert[, j] <- pert[, j] + eps[[j]]
          f1        <- .gp_function_draws(fit, pert, matched = FALSE)
          grad      <- as.vector(f1 - f0) / eps[[j]]
          rows[[length(rows) + 1L]] <- .term_row(node, parents[[j]], grad)
        }
      }
      do.call(rbind, rows)
    }
  )
)

# =============================================================================
# Internal GP helpers
# =============================================================================

#' Build the standardised GP eigen-basis for inputs `x`
#' @noRd
.gp_basis <- function(x, poly_degree, a, b, center = NULL, scale = NULL) {
  x      <- cbind(x)
  if (is.null(center)) center <- colMeans(x)
  if (is.null(scale))  scale  <- apply(x, 2, stats::sd)
  wx  <- BayesGPfit::GP.std.grids(x, center = center, scale = scale, max_range = 6)
  Psi <- BayesGPfit::GP.eigen.funcs.fast(wx, poly_degree, a, b)
  list(Psi = Psi, center = center, scale = scale,
       lambda = BayesGPfit::GP.eigen.value(poly_degree, a, b, d = ncol(x)))
}

#' Conjugate Bayesian linear regression on the GP eigen-basis
#'
#' Fits `y = Xmat theta + noise`, `theta_k ~ N(0, tau2 lambda_k)`, with
#' `(sigma2, tau2)` set by type-II maximum likelihood. Returns the exact Gaussian
#' posterior over `theta`, pre-drawn coefficient samples (for coherent per-draw
#' propagation), the residual-sd draws, and the log marginal likelihood.
#' @noRd
.gp_conj_fit <- function(y, x, poly_degree = 10L, a = 0.01, b = 1,
                         num_results = 1000L) {
  bs   <- .gp_basis(x, poly_degree, a, b)
  Xmat <- bs$Psi; lam <- bs$lambda
  ybar <- mean(y); yc <- y - ybar
  n <- length(yc); L <- ncol(Xmat)
  XtX <- crossprod(Xmat); Xty <- as.vector(crossprod(Xmat, yc)); yty <- sum(yc^2)

  neg_ll <- function(par) {
    s2 <- exp(par[[1]]); t2 <- exp(par[[2]])
    M  <- XtX / s2 + diag(1 / (t2 * lam), L)
    ch <- tryCatch(chol(M), error = function(e) NULL)
    if (is.null(ch)) return(1e10)
    sol  <- backsolve(ch, forwardsolve(t(ch), Xty))
    quad <- (yty - sum(Xty * sol) / s2) / s2
    ldC  <- n * log(s2) + sum(log(t2 * lam)) + 2 * sum(log(diag(ch)))
    0.5 * (n * log(2 * pi) + ldC + quad)
  }
  opt <- stats::optim(c(log(stats::var(yc)), 0), neg_ll, method = "Nelder-Mead")
  s2  <- exp(opt$par[[1]]); t2 <- exp(opt$par[[2]])

  # Exact Gaussian posterior over theta, and pre-drawn coefficient samples.
  Prec  <- XtX / s2 + diag(1 / (t2 * lam), L)
  chP   <- chol(Prec)
  Sig   <- chol2inv(chP)
  mu    <- as.vector(Sig %*% Xty) / s2
  theta <- t(mu + backsolve(chP, matrix(stats::rnorm(L * num_results), L, num_results)))

  # Residual-sd draws (scaled inverse chi-square around the type-II noise).
  sigma_draws <- sqrt(n * s2 / stats::rchisq(num_results, n))

  list(theta = theta, ybar = ybar, sigma_draws = sigma_draws,
       center = bs$center, scale = bs$scale, poly_degree = poly_degree,
       a = a, b = b, n_draws = num_results, logml = -opt$value)
}

#' Posterior draws of the GP function at query points
#'
#' `matched = TRUE` pairs query-point row `d` with coefficient draw `d` (coherent
#' per-draw propagation); otherwise every draw is evaluated at every point.
#' @noRd
.gp_function_draws <- function(fit, newx, matched = FALSE, n_draws = NULL) {
  bs  <- .gp_basis(newx, fit$poly_degree, fit$a, fit$b,
                   center = fit$center, scale = fit$scale)
  Psi <- bs$Psi                                    # [m x L]

  if (matched) {
    theta <- fit$theta[seq_len(n_draws), , drop = FALSE]
    return(fit$ybar + rowSums(Psi * theta))        # f^(d)(point_d), length n_draws
  }
  fit$ybar + Psi %*% t(fit$theta)                  # [m x n_draws]
}

#' Cap the polynomial degree so the basis size stays manageable in high dims
#'
#' The eigen-basis has `choose(poly_degree + d, d)` functions (monomials of total
#' degree <= poly_degree in `d` inputs), which grows fast with the number of
#' parents; shrink `poly_degree` until it is under `max_basis`.
#' @noRd
.cap_poly_degree <- function(poly_degree, d, max_basis = 400L) {
  pd <- as.integer(poly_degree)
  while (pd > 2L && choose(pd + d, d) > max_basis) pd <- pd - 1L
  pd
}

#' Gaussian marginal evidence of a root node (mean profiled out by centring)
#' @noRd
.gaussian_evidence <- function(y) {
  yc <- y - mean(y)
  n  <- length(yc)
  s2 <- mean(yc^2)
  if (s2 <= 0) s2 <- .Machine$double.eps
  -0.5 * (n * log(2 * pi) + n * log(s2) + n)
}

#' Type-II (empirical-Bayes) log marginal likelihood of the basis-linear GP model
#'
#' Evidence of `y ~ N(0, sigma2 I + tau2 * Xmat diag(lambda) Xmat^T)` with the
#' node mean removed by centring and `(sigma2, tau2)` maximised. Uses the
#' matrix-determinant lemma / Woodbury so the cost is O(L^3) rather than O(n^3).
#' @noRd
.gp_type2_evidence <- function(y, x, poly_degree = 10L, a = 0.01, b = 1) {
  bs   <- .gp_basis(x, poly_degree, a, b)
  Xmat <- bs$Psi; lam <- bs$lambda
  yc  <- y - mean(y)
  n   <- length(yc); L <- ncol(Xmat)
  XtX <- crossprod(Xmat); Xty <- as.vector(crossprod(Xmat, yc)); yty <- sum(yc^2)

  neg_ll <- function(par) {
    s2 <- exp(par[[1]]); t2 <- exp(par[[2]])
    M  <- diag(1 / (t2 * lam), L) + XtX / s2
    ch <- tryCatch(chol(M), error = function(e) NULL)
    if (is.null(ch)) return(1e10)
    sol  <- backsolve(ch, forwardsolve(t(ch), Xty))
    quad <- (yty - sum(Xty * sol) / s2) / s2
    ldC  <- n * log(s2) + sum(log(t2 * lam)) + 2 * sum(log(diag(ch)))
    0.5 * (n * log(2 * pi) + ldC + quad)
  }
  opt <- stats::optim(c(log(stats::var(yc)), 0), neg_ll, method = "Nelder-Mead")
  -opt$value
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
