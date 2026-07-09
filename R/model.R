#' KaguModel — Bayesian graphical causal model
#'
#' @description
#' The main user-facing class. Specify a DAG, optionally assign mechanisms
#' to nodes, fit to data, and extract causal effects.
#'
#' @examples
#' \dontrun{
#' library(kagu)
#'
#' model <- KaguModel$new(
#'   dag = list(
#'     age     = c(),
#'     smoking = c("age"),
#'     health  = c("smoking", "age")
#'   )
#' )
#'
#' model$fit(data)
#'
#' effect <- model$effects("smoking", "health")
#' effect$summary()
#' }
#'
#' @export
KaguModel <- R6::R6Class("KaguModel",
  public = list(
    #' @field dag Named list — the DAG specification.
    dag = NULL,
    #' @field mechanisms Named list of `Mechanism` instances, one per node.
    mechanisms = NULL,
    #' @field data The training `data.frame` (set after `$fit()`).
    data = NULL,
    #' @field traces Named list of mechanism fit objects, one per node (after fit).
    traces = NULL,
    #' @field .fitted Logical — whether `$fit()` has been called.
    .fitted = FALSE,

    # -------------------------------------------------------------------------
    # Constructor

    #' @description Create a new KaguModel.
    #' @param dag Named list mapping each node to a character vector of its
    #'   parent node names. Root nodes map to `c()`.
    #' @param mechanisms Optional named list of `Mechanism` instances. Any
    #'   node not specified receives a `GPMechanism` by default.
    initialize = function(dag, mechanisms = NULL) {
      validate_dag(dag)
      self$dag <- dag

      nodes <- names(dag)
      if (is.null(mechanisms)) mechanisms <- list()

      self$mechanisms <- setNames(
        lapply(nodes, function(n) {
          if (!is.null(mechanisms[[n]])) mechanisms[[n]] else GPMechanism$new()
        }),
        nodes
      )

      self$traces  <- list()
      self$.fitted <- FALSE
    },

    # -------------------------------------------------------------------------
    # Fitting

    #' @description Fit each node's conditional distribution in topological order.
    #' @param data A `data.frame` with one column per node.
    #' @param ... Additional arguments forwarded to each node's mechanism `$fit()`.
    #' @return `self` invisibly (for method chaining).
    fit = function(data, ...) {
      .check_data_nodes(self$dag, data)

      self$data <- data
      order     <- topological_sort(self$dag)
      n_nodes   <- length(order)

      cli::cli_alert_info("Fitting {n_nodes} node{?s} in topological order …")
      for (i in seq_along(order)) {
        node <- order[[i]]
        # Persistent per-node line (visible across consoles, RStudio, scripts);
        # printed *before* the fit so the user sees which node is running.
        cli::cli_alert("[{i}/{n_nodes}] fitting node {.field {node}} …")

        self$traces[[node]] <- fit_node(
          node      = node,
          parents   = self$dag[[node]],
          data      = data,
          mechanism = self$mechanisms[[node]],
          ...
        )
      }
      cli::cli_alert_success("Fitted {n_nodes} node{?s}.")

      self$.fitted <- TRUE
      invisible(self)
    },

    # -------------------------------------------------------------------------
    # Causal effects

    #' @description Estimate the causal effect of `source` on `target`.
    #' @param source Character scalar — the intervention (treatment) node.
    #' @param target Character scalar — the outcome node.
    #' @param values Optional numeric vector of length 2 `c(from, to)` giving
    #'   the explicit intervention contrast.
    #' @param std_units Logical — if `TRUE`, compute the effect of a 1-SD
    #'   increase centred at the mean.
    #' @param conditions Optional named list of node values to condition on
    #'   (fixes those nodes at the given values during propagation). Can also be
    #'   a list of such lists to compute and overlay multiple conditions.
    #' @param sweep Logical — if `TRUE`, compute the dose-response curve.
    #' @param sweep_n Integer — number of points in the sweep grid (default 50).
    #' @param sweep_range Numeric vector `c(min, max)` for the sweep grid.
    #'   Defaults to the observed range of `source`.
    #' @param hdi Numeric in (0, 1) — HDI probability (default 0.90).
    #' @return An `EffectResult` object.
    effects = function(source, target, values = NULL, std_units = FALSE,
                       conditions = NULL, sweep = FALSE, sweep_n = 50L,
                       sweep_range = NULL, hdi = 0.90) {
      if (!self$.fitted) stop("Call $fit() before $effects().")

      is_list_of_lists <- is.list(conditions) && length(conditions) > 0 && is.list(conditions[[1]]) && !is.data.frame(conditions[[1]])

      if (is_list_of_lists) {
        # Compute the first condition as the base result
        res <- compute_effect(
          self, source, target, values, std_units, conditions[[1]],
          sweep, sweep_n, sweep_range, hdi
        )
        
        label_fn <- function(cond) {
          paste(paste(names(cond), unlist(cond), sep="="), collapse=", ")
        }
        res$conditions_label <- label_fn(conditions[[1]])
        
        # Attach the remaining results for automatic comparison plotting
        res$compare_results <- list()
        for (i in seq_along(conditions)[-1]) {
          c_res <- compute_effect(
            self, source, target, values, std_units, conditions[[i]],
            sweep, sweep_n, sweep_range, hdi
          )
          lbl <- label_fn(conditions[[i]])
          c_res$conditions_label <- lbl
          res$compare_results[[lbl]] <- c_res
        }
        return(res)
      }

      compute_effect(
        model       = self,
        source      = source,
        target      = target,
        values      = values,
        std_units   = std_units,
        conditions  = conditions,
        sweep       = sweep,
        sweep_n     = sweep_n,
        sweep_range = sweep_range,
        hdi         = hdi
      )
    },

    # -------------------------------------------------------------------------
    # Summaries and diagnostics

    #' @description Per-node summary table.
    #'
    #' For each node, reports the direct local effect of each parent (the
    #' function's gradient at the parents' means — comparable to a regression
    #' coefficient) and the residual noise sd, each with posterior mean, sd and
    #' HDI.
    #' @param hdi_prob Numeric — HDI probability (default 0.90; currently fixed).
    #' @return A `tibble` with columns `node`, `term`, `mean`, `sd`,
    #'   `hdi_lower`, `hdi_upper`.
    summary = function(hdi_prob = 0.90) {
      if (!self$.fitted) stop("Call $fit() before $summary().")
      build_summary_table(self)
    },

    #' @description Fit diagnostics for each node.
    #'
    #' The Gaussian-process sampler produces a single chain, so r-hat / ESS do
    #' not apply; this reports the posterior draw count and the residual noise sd
    #' per node.
    #' @param node Optional character scalar. If `NULL`, runs for all nodes.
    #' @return A `tibble` with `node`, `n_chains`, `n_draws`, `sigma`, `sigma_sd`.
    diagnostics = function(node = NULL) {
      if (!self$.fitted) stop("Call $fit() before $diagnostics().")
      target_nodes <- if (is.null(node)) names(self$traces) else node
      results <- lapply(target_nodes, function(n) {
        mech  <- self$mechanisms[[n]]
        shape <- mech$posterior_shape(self$traces[[n]])
        terms <- mech$node_terms(n, self$dag[[n]], self$data, self$traces[[n]])
        sig   <- terms[terms$term == "sigma (noise)", ]
        tibble::tibble(
          node = n, n_chains = shape[["n_chains"]], n_draws = shape[["n_draws"]],
          sigma = sig$mean, sigma_sd = sig$sd
        )
      })
      do.call(rbind, results)
    },

    # -------------------------------------------------------------------------
    # Plots

    #' @description Visualise the DAG structure.
    #' @param node_pos Optional named list of `c(row, col)` grid coordinates
    #'   for manual node placement (overrides auto-layout for specified nodes).
    #' @return A `ggplot` object.
    plot_dag = function(node_pos = NULL) {
      kagu_plot_dag(self$dag, node_pos = node_pos)
    },

    #' @description Plot the posterior for a fitted node.
    #'
    #' Shows the node's direct local effects (each parent's gradient at the
    #' means) and residual noise, as posterior means with HDI intervals.
    #' @param node Character scalar — the node to plot.
    #' @return A `ggplot` object.
    plot_posterior = function(node) {
      if (!self$.fitted) stop("Call $fit() before $plot_posterior().")
      terms <- self$mechanisms[[node]]$node_terms(
        node, self$dag[[node]], self$data, self$traces[[node]]
      )
      kagu_plot_posterior(terms, node)
    },

    # -------------------------------------------------------------------------
    # Save / load

    #' @description Save the fitted model to disk.
    #' @param path File path (e.g. `"model.rds"`).
    #' @param include_data Logical — whether to save the training data alongside
    #'   the model (default `TRUE`, so `$effects()` works immediately on load).
    save = function(path, include_data = TRUE) {
      kagu_save(self, path, include_data = include_data)
    },

    # -------------------------------------------------------------------------
    # Print

    #' @description Print a concise model summary.
    #' @description Print summary of the model.
    #' @param ... Ignored.
    print = function(...) {
      status <- if (self$.fitted) "fitted" else "unfitted"
      nodes  <- names(self$dag)
      cat(sprintf("<KaguModel [%s]>\n", status))
      cat(sprintf("  Nodes (%d): %s\n", length(nodes),
                  paste(nodes, collapse = ", ")))
      if (self$.fitted) {
        n_draws <- tryCatch({
          s <- self$mechanisms[[nodes[[1]]]]$posterior_shape(self$traces[[nodes[[1]]]])
          sprintf("%d chain x %d draws", s[["n_chains"]], s[["n_draws"]])
        }, error = function(e) "?")
        cat(sprintf("  Posterior: %s\n", n_draws))
      }
      invisible(self)
    }
  )
)

# Static load method (called as KaguModel$load(path))
KaguModel$load <- function(path) kagu_load(path)

# Static structure-discovery method (called as KaguModel$discover(data, ...)).
# Delegates to kagu_discover(); see ?kagu_discover for full documentation.
KaguModel$discover <- function(data, nodes = NULL, dags = NULL, disallowed = NULL,
                               required = NULL, mechanisms = NULL, prior = "uniform",
                               allow_empty = FALSE, ...) {
  kagu_discover(
    data, nodes = nodes, dags = dags, disallowed = disallowed, required = required,
    mechanisms = mechanisms, prior = prior, allow_empty = allow_empty, ...
  )
}
