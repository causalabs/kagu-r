#' Prior distributions for Gaussian-process hyperparameters
#'
#' @description
#' Constructors for the prior families accepted by [gp_prior()]. Each returns a
#' small `kagu_prior` object (a family name plus its parameters) that Kagu can
#' both sample from (prior predictive simulation) and evaluate (MAP fitting).
#'
#' - `prior_loguniform(lower, upper)` - uniform on `log(x)` over `[lower, upper]`.
#' - `prior_lognormal(meanlog, sdlog)` - `log(x) ~ Normal(meanlog, sdlog)`.
#' - `prior_invgamma(shape, scale)` - inverse-gamma; a common lengthscale prior
#'   because it puts little mass near zero (very wiggly functions).
#' - `prior_halfnormal(scale)` - `|Normal(0, scale)|`.
#' - `prior_normal(mean, sd)` - unbounded; only valid for `root_mean`.
#'
#' @param lower,upper Positive bounds of the log-uniform prior.
#' @param meanlog,sdlog Mean and sd of `log(x)` for the log-normal prior.
#' @param shape Shape of the inverse-gamma prior.
#' @param scale Scale of the inverse-gamma or half-normal prior.
#' @param mean,sd Mean and sd of the normal prior.
#' @return A `kagu_prior` object.
#' @seealso [gp_prior()]
#' @examples
#' prior_lognormal(0, 0.5)
#' prior_invgamma(5, 5)
#' @name prior_families
NULL

#' @rdname prior_families
#' @export
prior_loguniform <- function(lower, upper) {
  stopifnot(is.numeric(lower), is.numeric(upper), lower > 0, upper > lower)
  .new_prior("loguniform", list(lower = lower, upper = upper))
}

#' @rdname prior_families
#' @export
prior_lognormal <- function(meanlog = 0, sdlog = 1) {
  stopifnot(is.numeric(meanlog), is.numeric(sdlog), sdlog > 0)
  .new_prior("lognormal", list(meanlog = meanlog, sdlog = sdlog))
}

#' @rdname prior_families
#' @export
prior_invgamma <- function(shape, scale) {
  stopifnot(is.numeric(shape), is.numeric(scale), shape > 0, scale > 0)
  .new_prior("invgamma", list(shape = shape, scale = scale))
}

#' @rdname prior_families
#' @export
prior_halfnormal <- function(scale = 1) {
  stopifnot(is.numeric(scale), scale > 0)
  .new_prior("halfnormal", list(scale = scale))
}

#' @rdname prior_families
#' @export
prior_normal <- function(mean = 0, sd = 1) {
  stopifnot(is.numeric(mean), is.numeric(sd), sd > 0)
  .new_prior("normal", list(mean = mean, sd = sd))
}

#' @export
format.kagu_prior <- function(x, ...) {
  args <- paste(names(x$params), signif(unlist(x$params), 4), sep = " = ",
                collapse = ", ")
  sprintf("%s(%s)", x$family, args)
}

#' @export
print.kagu_prior <- function(x, ...) {
  cat("<kagu_prior>", format(x), "\n")
  invisible(x)
}

#' Priors for a GPMechanism's hyperparameters
#'
#' @description
#' Bundles the priors a [GPMechanism] places on its hyperparameters. All are on
#' the **standardised** scale: before fitting, each parent and the node itself
#' are centred and scaled to unit sd, so `signal_sd = 1` means "functions vary
#' about as much as the data", and a lengthscale of 1 means "the function
#' changes over roughly one sd of the parent".
#'
#' The defaults are log-uniform over exactly the box the type-II maximum
#' likelihood search is bounded to. Maximising the likelihood inside that box
#' is the same as finding the posterior mode under this prior, so the defaults
#' make Kagu's implicit prior explicit without changing any fit. They are
#' deliberately vague; see `vignette("predictive_checks")` for what they imply
#' about the data and how to choose something more informative.
#'
#' Root nodes (no parents) are modelled by a Normal; `root_mean` and `root_sd`
#' are used for prior predictive simulation only - a root node's posterior uses
#' a non-informative reference prior, which these weak priors barely differ
#' from once there is any data.
#'
#' @param lengthscale Prior on each parent's kernel lengthscale.
#' @param signal_sd Prior on the kernel signal sd (the function's scale).
#' @param noise_sd Prior on the observation noise sd.
#' @param root_mean Prior on a root node's (standardised) mean.
#' @param root_sd Prior on a root node's (standardised) sd.
#' @return A `kagu_gp_prior` object.
#' @seealso [prior_families], [GPMechanism]
#' @examples
#' gp_prior()   # the implicit default
#'
#' # Weakly informative: smooth functions, signal and noise on the data's scale
#' gp_prior(
#'   lengthscale = prior_lognormal(0, 0.5),
#'   signal_sd   = prior_halfnormal(1),
#'   noise_sd    = prior_halfnormal(1)
#' )
#' @export
gp_prior <- function(lengthscale = prior_loguniform(exp(-3), exp(4)),
                     signal_sd   = prior_loguniform(exp(-2), exp(2)),
                     noise_sd    = prior_loguniform(exp(-3), exp(1)),
                     root_mean   = prior_normal(0, 1),
                     root_sd     = prior_lognormal(0, 0.5)) {
  priors <- list(lengthscale = lengthscale, signal_sd = signal_sd,
                 noise_sd = noise_sd, root_mean = root_mean, root_sd = root_sd)
  for (nm in names(priors)) {
    if (!inherits(priors[[nm]], "kagu_prior")) {
      stop(sprintf("`%s` must be a prior, e.g. prior_lognormal(0, 1).", nm))
    }
    if (nm != "root_mean" && !.prior_positive(priors[[nm]])) {
      stop(sprintf("`%s` needs a prior on positive values (not %s).",
                   nm, priors[[nm]]$family))
    }
  }
  structure(priors, class = "kagu_gp_prior")
}

#' @export
print.kagu_gp_prior <- function(x, ...) {
  cat("<kagu_gp_prior> (standardised scale)\n")
  for (nm in names(x)) cat(sprintf("  %-12s %s\n", nm, format(x[[nm]])))
  invisible(x)
}

# -----------------------------------------------------------------------------
# Internal helpers
# -----------------------------------------------------------------------------

#' @noRd
.new_prior <- function(family, params) {
  structure(list(family = family, params = params), class = "kagu_prior")
}

#' Is the prior supported on positive values only?
#' @noRd
.prior_positive <- function(p) p$family != "normal"

#' Draw `n` values from a prior
#' @noRd
.prior_sample <- function(p, n) {
  a <- p$params
  switch(p$family,
    loguniform = exp(stats::runif(n, log(a$lower), log(a$upper))),
    lognormal  = stats::rlnorm(n, a$meanlog, a$sdlog),
    invgamma   = 1 / stats::rgamma(n, shape = a$shape, rate = a$scale),
    halfnormal = abs(stats::rnorm(n, 0, a$scale)),
    normal     = stats::rnorm(n, a$mean, a$sd)
  )
}

#' Log density of a prior at `x`
#' @noRd
.prior_logdens <- function(p, x) {
  a <- p$params
  switch(p$family,
    loguniform = ifelse(x >= a$lower & x <= a$upper,
                        -log(x) - log(log(a$upper / a$lower)), -Inf),
    lognormal  = stats::dlnorm(x, a$meanlog, a$sdlog, log = TRUE),
    invgamma   = ifelse(x > 0, a$shape * log(a$scale) - lgamma(a$shape) -
                          (a$shape + 1) * log(x) - a$scale / x, -Inf),
    halfnormal = ifelse(x >= 0, log(2) + stats::dnorm(x, 0, a$scale, log = TRUE),
                        -Inf),
    normal     = stats::dnorm(x, a$mean, a$sd, log = TRUE)
  )
}

#' Log density of a positive prior on `u = log(x)` (includes the Jacobian)
#' @noRd
.prior_logdens_log <- function(p, u) .prior_logdens(p, exp(u)) + u

#' `c(lower, upper)` on the log scale if the prior is log-uniform, else `NULL`
#' @noRd
.prior_log_bounds <- function(p) {
  if (p$family == "loguniform") log(c(p$params$lower, p$params$upper)) else NULL
}
