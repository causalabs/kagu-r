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
#' hyperparameters (lengthscales, signal and noise variance) are set by type-II
#' maximum likelihood; there is no MCMC. A node with no parents is modelled by
#' its marginal (a Normal). This is Kagu's default and only mechanism.
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

    #' @description Create a new GPMechanism.
    #' @param num_results Integer - number of posterior draws to keep.
    #' @param n_features Integer - number of random Fourier features.
    #' @param jitter Numeric - diagonal jitter added to the kernel.
    initialize = function(num_results = 1000L, n_features = 300L, jitter = 1e-6) {
      self$num_results <- as.integer(num_results)
      self$n_features  <- as.integer(n_features)
      self$jitter      <- jitter
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

      fit <- .gp_core(y, as.matrix(data[parents]), self$jitter)
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
      .gp_core(y, as.matrix(data[parents]), self$jitter)$logml
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

#' Fit an exact GP by type-II ML; return hypers, Cholesky, and log evidence
#' @noRd
.gp_core <- function(y, X, jitter = 1e-6) {
  X <- as.matrix(X); n <- nrow(X); d <- ncol(X)
  mx <- colMeans(X); sx <- apply(X, 2, stats::sd); sx[sx == 0] <- 1
  Xs <- sweep(sweep(X, 2, mx, "-"), 2, sx, "/")
  my <- mean(y); sy <- stats::sd(y); if (sy == 0) sy <- 1
  ys <- (y - my) / sy

  nll <- function(p) {
    ls <- exp(p[seq_len(d)]); sf2 <- exp(p[d + 1]); sn2 <- exp(p[d + 2])
    K  <- .rbf(Xs, Xs, ls, sf2) + diag(sn2 + jitter, n)
    ch <- tryCatch(chol(K), error = function(e) NULL); if (is.null(ch)) return(1e10)
    al <- backsolve(ch, forwardsolve(t(ch), ys))
    0.5 * sum(ys * al) + sum(log(diag(ch))) + 0.5 * n * log(2 * pi)
  }
  opt <- stats::optim(c(rep(0, d), 0, -1), nll, method = "L-BFGS-B",
                      lower = c(rep(-3, d), -4, -6), upper = c(rep(4, d), 4, 2))
  ls <- exp(opt$par[seq_len(d)]); sf2 <- exp(opt$par[d + 1]); sn2 <- exp(opt$par[d + 2])
  K  <- .rbf(Xs, Xs, ls, sf2) + diag(sn2 + jitter, n)
  ch <- chol(K); alpha <- backsolve(ch, forwardsolve(t(ch), ys))

  list(kind = "gp", Xs = Xs, ys = ys, ch = ch, alpha = alpha,
       ls = ls, sf2 = sf2, sn2 = sn2, mx = mx, sx = sx, my = my, sy = sy,
       n = n, d = d, logml = -opt$value)
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
