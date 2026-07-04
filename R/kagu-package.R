#' kagu: Bayesian Graphical Causal Models
#'
#' @description
#' Fit Bayesian graphical causal models specified as directed acyclic graphs
#' (DAGs). Each node's conditional distribution `P(X | parents)` is modelled as
#' a Gaussian process of its parents, exploiting the Markov blanket
#' factorisation. Causal effects are estimated by propagating interventions
#' forward through the structural model using the do-operator, and causal
#' structure can be discovered as a posterior distribution over DAGs.
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
#' @importFrom posterior summarise_draws
#' @importFrom bayestestR hdi
#' @importFrom tibble as_tibble
#' @importFrom stats sd
## usethis namespace: end
NULL
