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
    #' @field traces Named list of `brmsfit` objects, one per node (after fit).
    traces = NULL,
    #' @field .fitted Logical — whether `$fit()` has been called.
    .fitted = FALSE,

    # -------------------------------------------------------------------------
    # Constructor

    #' @description Create a new KaguModel.
    #' @param dag Named list mapping each node to a character vector of its
    #'   parent node names. Root nodes map to `c()`.
    #' @param mechanisms Optional named list of `Mechanism` instances. Any
    #'   node not specified receives a `LinearMechanism` by default.
    initialize = function(dag, mechanisms = NULL) {
      validate_dag(dag)
      self$dag <- dag

      nodes <- names(dag)
      if (is.null(mechanisms)) mechanisms <- list()

      self$mechanisms <- setNames(
        lapply(nodes, function(n) {
          if (!is.null(mechanisms[[n]])) mechanisms[[n]] else LinearMechanism$new()
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
    #' @param draws Integer — post-warmup draws per chain (default 1000).
    #' @param tune Integer — warmup draws per chain (default 1000).
    #' @param chains Integer — number of MCMC chains (default 4).
    #' @param ... Additional arguments forwarded to `brms::brm()`.
    #' @return `self` invisibly (for method chaining).
    fit = function(data, draws = 1000L, tune = 1000L, chains = 4L, ...) {
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
          draws     = draws,
          tune      = tune,
          chains    = chains,
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
    #'   (fixes those nodes at the given values during propagation).
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

    #' @description Coefficient table across all node mechanisms.
    #' @param hdi_prob Numeric — HDI probability (default 0.90).
    #' @return A `tibble` with columns `node`, `variable`, `mean`, `sd`,
    #'   `hdi_lower`, `hdi_upper`, `rhat`, `ess_bulk`.
    summary = function(hdi_prob = 0.90) {
      if (!self$.fitted) stop("Call $fit() before $summary().")
      build_summary_table(self$traces, hdi_prob)
    },

    #' @description R-hat and ESS diagnostics for fitted nodes.
    #' @param node Optional character scalar. If `NULL`, runs for all nodes.
    #' @return A `tibble` from `posterior::summarise_draws()`.
    diagnostics = function(node = NULL) {
      if (!self$.fitted) stop("Call $fit() before $diagnostics().")
      target_nodes <- if (is.null(node)) names(self$traces) else node
      results <- lapply(target_nodes, function(n) {
        df       <- posterior::summarise_draws(self$traces[[n]])
        df$node  <- n
        df[, c("node", setdiff(names(df), "node"))]
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
    #' @param node Character scalar — the node to plot.
    #' @return A `ggplot` object.
    plot_posterior = function(node) {
      if (!self$.fitted) stop("Call $fit() before $plot_posterior().")
      kagu_plot_posterior(self$traces[[node]], node)
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
    print = function(...) {
      status <- if (self$.fitted) "fitted" else "unfitted"
      nodes  <- names(self$dag)
      cat(sprintf("<KaguModel [%s]>\n", status))
      cat(sprintf("  Nodes (%d): %s\n", length(nodes),
                  paste(nodes, collapse = ", ")))
      if (self$.fitted) {
        n_draws <- tryCatch({
          s <- .chain_draw_shape(self$traces[[nodes[[1]]]])
          sprintf("%d chains x %d draws", s[["n_chains"]], s[["n_draws"]])
        }, error = function(e) "?")
        cat(sprintf("  Posterior: %s\n", n_draws))
      }
      invisible(self)
    }
  )
)

# Static load method (called as KaguModel$load(path))
KaguModel$load <- function(path) kagu_load(path)
