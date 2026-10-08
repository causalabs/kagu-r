#' Prior / posterior predictive check results
#'
#' @description
#' Holds replicated datasets from a [KaguModel]'s `$prior_predictive()` or
#' `$posterior_predictive()` alongside the observed data. `$y()` and `$yrep()`
#' return data in the shapes the **bayesplot** package expects (a vector and a
#' `[draws, observations]` matrix), so any `bayesplot::ppc_*()` function works
#' directly:
#'
#' ```r
#' pp <- model$posterior_predictive()
#' bayesplot::ppc_dens_overlay(pp$y("y"), pp$yrep("y")[1:50, ])
#' ```
#'
#' `$summary()` reports predictive p-values for summary statistics, and
#' `$cor_check()` / `$plot_cor()` compare the correlation between every pair of
#' variables with what the model reproduces - a whole-graph check (for `joint`
#' and `prior` results) with no direct bayesplot equivalent.
#'
#' @export
PredictiveResult <- R6::R6Class("PredictiveResult",
  public = list(
    #' @field type Character - `"prior"`, `"conditional"` or `"joint"`.
    type = NULL,
    #' @field data The observed `data.frame`.
    data = NULL,
    #' @field replicates Named list of `[ndraws, n]` matrices, one per node.
    replicates = NULL,
    #' @field dag The model's DAG.
    dag = NULL,

    #' @description Create a PredictiveResult (normally via a [KaguModel]).
    #' @param type Character - `"prior"`, `"conditional"` or `"joint"`.
    #' @param data The observed `data.frame`.
    #' @param replicates Named list of `[ndraws, n]` matrices, one per node.
    #' @param dag The model's DAG.
    initialize = function(type, data, replicates, dag) {
      self$type       <- match.arg(type, c("prior", "conditional", "joint"))
      self$data       <- as.data.frame(data)[names(replicates)]
      self$replicates <- replicates
      self$dag        <- dag
    },

    #' @description Node names.
    #' @return Character vector.
    nodes = function() names(self$replicates),

    #' @description Number of replicated datasets.
    #' @return Integer.
    ndraws = function() nrow(self$replicates[[1]]),

    #' @description Observed values of a node (bayesplot's `y`).
    #' @param node Character scalar - the node name.
    #' @return Numeric vector of length `n`.
    y = function(node) {
      private$check_node(node)
      as.numeric(self$data[[node]])
    },

    #' @description Replicated values of a node (bayesplot's `yrep`).
    #' @param node Character scalar - the node name.
    #' @return Numeric `[ndraws, n]` matrix.
    yrep = function(node) {
      private$check_node(node)
      self$replicates[[node]]
    },

    #' @description Predictive p-values for summary statistics.
    #'
    #' For each node and statistic `T`, compares `T(y)` with the distribution of
    #' `T(yrep)` over the replicated datasets. `p_value = Pr(T(yrep) >= T(y))`:
    #' values near 0 or 1 mean the observed statistic is extreme relative to
    #' what the model generates.
    #' @param stats Character vector of function names (e.g. `"sd"`), or a named
    #'   list of functions each mapping a numeric vector to a scalar.
    #' @param nodes Optional character vector of nodes (default: all).
    #' @param prob Numeric - probability mass of the reported replicate interval.
    #' @return A `tibble` with columns `node`, `stat`, `observed`, `rep_mean`,
    #'   `rep_lower`, `rep_upper`, `p_value`.
    summary = function(stats = c("mean", "sd", "min", "max"), nodes = NULL,
                       prob = 0.90) {
      fns <- .resolve_stats(stats)
      if (is.null(nodes)) nodes <- self$nodes()
      a <- (1 - prob) / 2
      rows <- list()
      for (nd in nodes) {
        y <- self$y(nd); yrep <- self$yrep(nd)
        for (st in names(fns)) {
          t_obs <- fns[[st]](y)
          t_rep <- apply(yrep, 1, fns[[st]])
          q <- stats::quantile(t_rep, c(a, 1 - a), names = FALSE)
          rows[[length(rows) + 1L]] <- tibble::tibble(
            node = nd, stat = st, observed = t_obs, rep_mean = mean(t_rep),
            rep_lower = q[1], rep_upper = q[2], p_value = mean(t_rep >= t_obs)
          )
        }
      }
      do.call(rbind, rows)
    },

    #' @description Compare observed pairwise correlations with replicated ones.
    #'
    #' A whole-graph check: under a joint (or prior) predictive, each replicated
    #' dataset carries the dependence structure the DAG implies, so a pair whose
    #' observed correlation falls outside its replicated range points to a
    #' missing (or spurious) edge or a misspecified mechanism. Not available for
    #' `conditional` results, whose nodes are simulated separately.
    #' @param prob Numeric - probability mass of the replicate interval.
    #' @param method Correlation method passed to [stats::cor()].
    #' @return A `tibble` with columns `var1`, `var2`, `observed`, `rep_median`,
    #'   `rep_lower`, `rep_upper`, `p_value`, `outside`.
    cor_check = function(prob = 0.90, method = "pearson") {
      if (self$type == "conditional") {
        stop(paste(
          "Correlations between separately simulated nodes are not meaningful",
          "for a conditional check. Use $posterior_predictive(type = \"joint\")",
          "or $prior_predictive()."
        ))
      }
      nodes <- self$nodes()
      idx   <- which(upper.tri(diag(length(nodes))), arr.ind = TRUE)
      idx   <- idx[order(idx[, "row"], idx[, "col"]), , drop = FALSE]
      a     <- (1 - prob) / 2

      rows <- lapply(seq_len(nrow(idx)), function(k) {
        v1 <- nodes[[idx[k, "row"]]]; v2 <- nodes[[idx[k, "col"]]]
        obs <- stats::cor(self$y(v1), self$y(v2), method = method)
        r1  <- self$yrep(v1); r2 <- self$yrep(v2)
        rep <- suppressWarnings(vapply(seq_len(nrow(r1)), function(s) {
          stats::cor(r1[s, ], r2[s, ], method = method)
        }, numeric(1)))
        rep <- rep[is.finite(rep)]
        q   <- stats::quantile(rep, c(a, 0.5, 1 - a), names = FALSE)
        tibble::tibble(
          var1 = v1, var2 = v2, observed = obs, rep_median = q[2],
          rep_lower = q[1], rep_upper = q[3], p_value = mean(rep >= obs),
          outside = obs < q[1] | obs > q[3]
        )
      })
      do.call(rbind, rows)
    },

    #' @description Plot `$cor_check()`: each pair's replicated correlation
    #'   interval with the observed value; pairs outside their interval are
    #'   highlighted.
    #' @param prob Numeric - probability mass of the replicate interval.
    #' @param method Correlation method passed to [stats::cor()].
    #' @return A `ggplot` object.
    plot_cor = function(prob = 0.90, method = "pearson") {
      cc <- self$cor_check(prob = prob, method = method)
      cc$pair <- factor(paste(cc$var1, cc$var2, sep = " - "),
                        levels = rev(paste(cc$var1, cc$var2, sep = " - ")))
      cc$status <- ifelse(cc$outside, "outside replicated range", "consistent")

      ggplot(cc, aes(y = .data$pair)) +
        geom_segment(aes(x = .data$rep_lower, xend = .data$rep_upper,
                         yend = .data$pair),
                     linewidth = 2.2, colour = "#b2dfdb", lineend = "round") +
        geom_point(aes(x = .data$rep_median), colour = "#26a69a", size = 1.8) +
        geom_point(aes(x = .data$observed, colour = .data$status),
                   shape = 18, size = 3.6) +
        ggplot2::scale_colour_manual(
          values = c("consistent" = "#1c2826", "outside replicated range" = "#e65100"),
          drop = FALSE
        ) +
        labs(
          x = "Correlation", y = NULL, colour = "Observed",
          title = sprintf("Pairwise correlations: observed vs %s replicates",
                          if (self$type == "prior") "prior" else "posterior"),
          subtitle = sprintf("Bars: %d%% interval of replicated correlations",
                             round(100 * prob))
        ) +
        theme_minimal() +
        theme(legend.position = "bottom")
    },

    #' @description Print a short description.
    #' @param ... Ignored.
    print = function(...) {
      cat(sprintf("<PredictiveResult [%s]: %d nodes, %d replicated datasets of n = %d>\n",
                  self$type, length(self$nodes()), self$ndraws(), nrow(self$data)))
      cat("  Use $y(node) / $yrep(node) with bayesplot::ppc_*(), or $summary().\n")
      invisible(self)
    }
  ),

  private = list(
    check_node = function(node) {
      if (!is.character(node) || length(node) != 1L || !node %in% self$nodes()) {
        stop(sprintf("Unknown node '%s'. Available: %s.", format(node),
                     paste(self$nodes(), collapse = ", ")))
      }
    }
  )
)

#' Turn a `stats` argument into a named list of functions
#' @noRd
.resolve_stats <- function(stats) {
  if (is.character(stats)) {
    return(setNames(lapply(stats, match.fun), stats))
  }
  if (is.list(stats) && !is.null(names(stats)) && all(nzchar(names(stats))) &&
      all(vapply(stats, is.function, logical(1)))) {
    return(stats)
  }
  stop("`stats` must be function names or a named list of functions.")
}
