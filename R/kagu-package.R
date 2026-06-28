#' kagu: Bayesian Graphical Causal Models
#'
#' @description
#' Fit Bayesian graphical causal models specified as directed acyclic graphs
#' (DAGs). Each node's conditional distribution `P(X | parents)` is fitted as
#' an independent Bayesian model via **brms**, exploiting the Markov blanket
#' factorisation. Causal effects are estimated by propagating interventions
#' forward through the structural model using the do-operator.
#'
#' ## Main entry point
#'
#' ```r
#' model <- KaguModel$new(dag = list(x = c(), y = c("x")))
#' model$fit(data)
#' effect <- model$effects("x", "y")
#' effect$summary()
#' ```
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import R6
#' @importFrom ggplot2 ggplot aes geom_density geom_ribbon geom_line
#'   geom_segment geom_point geom_text annotate labs theme_minimal theme_void
#'   theme element_blank coord_equal
#' @importFrom posterior as_draws_array as_draws_df summarise_draws rhat ess_bulk
#' @importFrom bayestestR hdi
#' @importFrom tibble as_tibble
#' @importFrom stats as.formula sd
## usethis namespace: end
NULL
