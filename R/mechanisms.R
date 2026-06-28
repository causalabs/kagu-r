#' Mechanism base class and built-in implementations
#'
#' A `Mechanism` encapsulates the probabilistic model for a single node's
#' conditional distribution `P(X | parents)`. Subclass and implement
#' `build_formula()` and `predict_mean()` to add new mechanism types.
#'
#' @name mechanisms
NULL

# =============================================================================
# Mechanism ABC
# =============================================================================

#' Mechanism abstract base class
#'
#' @description
#' All mechanisms must implement two methods:
#' - `$build_formula(node, parents)` — returns a `brms::brmsformula`.
#' - `$predict_mean(node, parents, parent_values, fit)` — given posterior draws
#'   and parent values as `[n_chains, n_draws]` matrices, returns the predicted
#'   conditional mean as a `[n_chains, n_draws]` matrix.
#'
#' @export
Mechanism <- R6::R6Class("Mechanism",
  public = list(
    #' @description Build the brms formula for this node.
    #' @param node Character scalar — the node name (response variable).
    #' @param parents Character vector of parent node names.
    #' @return A `brms::brmsformula`.
    build_formula = function(node, parents) {
      stop("Mechanism$build_formula() is abstract — implement in a subclass.")
    },

    #' @description Predict the conditional mean for each posterior draw.
    #' @param node Character scalar — the node name.
    #' @param parents Character vector of parent node names.
    #' @param parent_values Named list of `[n_chains, n_draws]` matrices,
    #'   one per parent, giving the parent values for each posterior draw.
    #' @param fit A `brmsfit` object from fitting this node.
    #' @return A `[n_chains, n_draws]` numeric matrix.
    predict_mean = function(node, parents, parent_values, fit) {
      stop("Mechanism$predict_mean() is abstract — implement in a subclass.")
    },

    #' @description brms family for this mechanism.
    #' @return A `brmsfamily` object (e.g. `brms::gaussian()`).
    family = function() {
      stop("Mechanism$family() is abstract — implement in a subclass.")
    }
  )
)

# =============================================================================
# LinearMechanism
# =============================================================================

#' Linear Gaussian mechanism
#'
#' @description
#' Models the conditional distribution as:
#' ```
#' alpha    ~ Normal(0, 10)
#' beta_j   ~ Normal(0, 2)   for each parent j
#' sigma    ~ HalfNormal(1)
#' mu       = alpha + sum(beta_j * Pa_j)
#' node     ~ Normal(mu, sigma)
#' ```
#' Suitable for continuous, unbounded variables.
#'
#' @export
#'
#' @examples
#' mech <- LinearMechanism$new()
LinearMechanism <- R6::R6Class("LinearMechanism",
  inherit = Mechanism,
  public = list(
    #' @field prior_alpha Prior SD for the intercept (default 10).
    prior_alpha = 10,
    #' @field prior_beta Prior SD for coefficients (default 2).
    prior_beta = 2,
    #' @field prior_sigma Scale for the HalfNormal prior on sigma (default 1).
    prior_sigma = 1,

    #' @description Create a new LinearMechanism.
    #' @param prior_alpha Prior SD for the intercept.
    #' @param prior_beta Prior SD for regression coefficients.
    #' @param prior_sigma Scale for the HalfNormal prior on sigma.
    initialize = function(prior_alpha = 10, prior_beta = 2, prior_sigma = 1) {
      self$prior_alpha <- prior_alpha
      self$prior_beta  <- prior_beta
      self$prior_sigma <- prior_sigma
    },

    #' @description Build the brms formula.
    build_formula = function(node, parents) {
      if (length(parents) == 0) {
        rhs <- "1"
      } else {
        rhs <- paste(parents, collapse = " + ")
      }
      brms::bf(stats::as.formula(paste(node, "~", rhs)))
    },

    #' @description brms Gaussian family.
    family = function() gaussian(),

    #' @description Predict conditional mean using posterior draws.
    #'
    #' Evaluates `mu = alpha + sum(beta_j * parent_j)` for each
    #' (chain, draw) pair, returning a `[n_chains, n_draws]` matrix.
    predict_mean = function(node, parents, parent_values, fit) {
      draws    <- posterior::as_draws_array(fit)
      n_chains <- dim(draws)[[2]]
      n_draws  <- dim(draws)[[1]]

      # Extract intercept: draws[iter, chain, var] -> t() -> [n_chains, n_draws]
      alpha <- t(draws[, , "b_Intercept", drop = TRUE])   # [n_chains, n_draws]
      mu    <- alpha

      for (parent in parents) {
        beta_name <- paste0("b_", parent)
        beta      <- t(draws[, , beta_name, drop = TRUE])  # [n_chains, n_draws]
        mu        <- mu + beta * parent_values[[parent]]
      }

      mu
    }
  )
)
