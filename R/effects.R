#' Causal effect estimation via the do-operator
#'
#' @name effects
NULL

# =============================================================================
# EffectResult R6 class
# =============================================================================

#' Result of a causal effect query
#'
#' @description
#' Returned by `KaguModel$effects()`. Holds the full posterior over the
#' causal effect along with metadata.
#'
#' - **Scalar effect**: `samples` is a `[n_chains, n_draws]` matrix.
#' - **Sweep effect**: `samples` is a `[n_sweep, n_chains, n_draws]` array.
#'
#' @export
EffectResult <- R6::R6Class("EffectResult",
  public = list(
    #' @field source Character — the intervention node.
    source = NULL,
    #' @field target Character — the outcome node.
    target = NULL,
    #' @field from_value Numeric — intervention baseline value (or sweep grid).
    from_value = NULL,
    #' @field to_value Numeric — intervention target value (or sweep grid).
    to_value = NULL,
    #' @field samples Numeric array — posterior effect samples.
    samples = NULL,
    #' @field hdi_prob Numeric — HDI probability used for summaries and plots.
    hdi_prob = NULL,
    #' @field std_units Logical — whether the effect is in SD units.
    std_units = NULL,
    #' @field conditions Named list of conditioned node values, or `NULL`.
    conditions = NULL,
    #' @field sweep_values Numeric vector of sweep grid values, or `NULL`.
    sweep_values = NULL,

    #' @description Create an EffectResult (normally called by `compute_effect`).
    initialize = function(source, target, from_value, to_value, samples,
                          hdi = 0.90, std_units = FALSE, conditions = NULL,
                          sweep_values = NULL) {
      self$source       <- source
      self$target       <- target
      self$from_value   <- from_value
      self$to_value     <- to_value
      self$samples      <- samples
      self$hdi_prob     <- hdi
      self$std_units    <- std_units
      self$conditions   <- conditions
      self$sweep_values <- sweep_values
    },

    #' @description Is this a sweep result?
    #' @return Logical scalar.
    is_sweep = function() length(dim(self$samples)) == 3L,

    # -------------------------------------------------------------------------
    # Summary

    #' @description Tabular summary of the effect posterior.
    #'
    #' Returns a `tibble`, consistent with `KaguModel$summary()`. The HDI
    #' columns use the same `hdi_lower` / `hdi_upper` naming as
    #' [build_summary_table()].
    #'
    #' - **Scalar**: a one-row tibble with `source`, `target`, `from`, `to`,
    #'   `mean`, `sd`, `hdi_lower`, `hdi_upper`.
    #' - **Sweep**: one row per grid point with `source`, `target`, `x`,
    #'   `mean`, `sd`, `hdi_lower`, `hdi_upper`.
    #'
    #' @return A `tibble`.
    summary = function() {
      if (self$is_sweep()) {
        n_sweep <- dim(self$samples)[[1]]
        flat    <- matrix(self$samples, nrow = n_sweep)
        h       <- t(vapply(seq_len(n_sweep), function(i) {
          hh <- bayestestR::hdi(flat[i, ], ci = self$hdi_prob)
          c(hh$CI_low, hh$CI_high)
        }, numeric(2)))

        tibble::tibble(
          source    = self$source,
          target    = self$target,
          x         = self$sweep_values,
          mean      = rowMeans(flat),
          sd        = apply(flat, 1, stats::sd),
          hdi_lower = h[, 1],
          hdi_upper = h[, 2]
        )
      } else {
        flat <- as.vector(self$samples)
        h    <- bayestestR::hdi(flat, ci = self$hdi_prob)

        tibble::tibble(
          source    = self$source,
          target    = self$target,
          from      = self$from_value[[1]],
          to        = self$to_value[[length(self$to_value)]],
          mean      = mean(flat),
          sd        = stats::sd(flat),
          hdi_lower = h$CI_low,
          hdi_upper = h$CI_high
        )
      }
    },

    # -------------------------------------------------------------------------
    # Plot

    #' @description Plot the posterior effect distribution.
    #'
    #' - **Scalar**: density plot with mean point and HDI bar.
    #' - **Sweep**: dose-response line with HDI ribbon.
    #'
    #' @return A `ggplot` object.
    plot = function() {
      if (self$is_sweep()) {
        .plot_sweep(self)
      } else {
        .plot_scalar(self)
      }
    },

    # -------------------------------------------------------------------------
    # Diagnostics

    #' @description R-hat and ESS diagnostics for the effect posterior.
    #'
    #' Chain structure is preserved, so r-hat values are genuine.
    #'
    #' @return A `tibble` from `posterior::summarise_draws()`.
    diagnostics = function() {
      if (self$is_sweep()) {
        stop("$diagnostics() is not supported for sweep effects.")
      }
      draws <- .to_draws_array(self$samples, varname = "effect")
      posterior::summarise_draws(draws)
    },

    #' @description Print method — shows the summary table.
    print = function(...) {
      kind <- if (self$is_sweep()) "sweep" else "scalar"
      cat(sprintf("<EffectResult [%s]: %s → %s>\n",
                  kind, self$source, self$target))
      print(self$summary())
      invisible(self)
    }
  )
)

# =============================================================================
# compute_effect
# =============================================================================

#' Compute the causal effect of source on target
#'
#' Called internally by `KaguModel$effects()`. Not usually called directly.
#'
#' @param model A fitted `KaguModel`.
#' @param source,target Character scalars — intervention and outcome nodes.
#' @param values Optional `c(from, to)` intervention contrast.
#' @param std_units Logical — 1-SD effect.
#' @param conditions Optional named list of fixed node values.
#' @param sweep Logical — compute dose-response curve.
#' @param sweep_n,sweep_range Sweep grid parameters.
#' @param hdi Numeric — HDI probability.
#' @return An `EffectResult`.
#' @export
compute_effect <- function(model, source, target,
                            values = NULL, std_units = FALSE,
                            conditions = NULL, sweep = FALSE,
                            sweep_n = 50L, sweep_range = NULL,
                            hdi = 0.90) {
  nodes <- names(model$dag)
  if (!source %in% nodes) stop(sprintf("Node '%s' not in DAG.", source))
  if (!target %in% nodes) stop(sprintf("Node '%s' not in DAG.", target))

  # Zero effect if source has no causal path to target
  if (!is_ancestor(model$dag, source, target)) {
    shape    <- .chain_draw_shape(model$traces[[target]])
    zero_mat <- matrix(0, nrow = shape[["n_chains"]], ncol = shape[["n_draws"]])
    return(EffectResult$new(
      source = source, target = target,
      from_value = 0, to_value = 0,
      samples = zero_mat, hdi = hdi,
      std_units = FALSE, conditions = conditions
    ))
  }

  src_data  <- model$data[[source]]
  src_mean  <- mean(src_data)
  src_sd    <- stats::sd(src_data)
  shape     <- .chain_draw_shape(model$traces[[target]])
  n_chains  <- shape[["n_chains"]]
  n_draws   <- shape[["n_draws"]]

  # ---- Sweep ----------------------------------------------------------------
  if (sweep) {
    if (is.null(sweep_range)) sweep_range <- range(src_data)
    grid     <- seq(sweep_range[[1]], sweep_range[[2]], length.out = sweep_n)
    samples  <- array(NA_real_, dim = c(sweep_n, n_chains, n_draws))

    for (i in seq_along(grid)) {
      samples[i, , ] <- .propagate(model, source, target, grid[[i]],
                                   conditions, n_chains, n_draws)
    }

    return(EffectResult$new(
      source = source, target = target,
      from_value = grid, to_value = grid,
      samples = samples, hdi = hdi,
      std_units = FALSE, conditions = conditions,
      sweep_values = grid
    ))
  }

  # ---- Scalar ---------------------------------------------------------------
  if (!is.null(values)) {
    from_val     <- values[[1]]
    to_val       <- values[[2]]
    scale        <- 1.0
    display_from <- from_val
    display_to   <- to_val
  } else if (std_units) {
    from_val     <- src_mean - 0.5 * src_sd
    to_val       <- src_mean + 0.5 * src_sd
    scale        <- 1.0
    display_from <- from_val
    display_to   <- to_val
  } else {
    # Epsilon central-difference gradient at the mean
    eps          <- src_sd * 1e-5
    from_val     <- src_mean - eps
    to_val       <- src_mean + eps
    scale        <- 1.0 / (2.0 * eps)
    display_from <- src_mean
    display_to   <- src_mean + 1.0
  }

  out_from <- .propagate(model, source, target, from_val, conditions, n_chains, n_draws)
  out_to   <- .propagate(model, source, target, to_val,   conditions, n_chains, n_draws)

  EffectResult$new(
    source     = source,
    target     = target,
    from_value = display_from,
    to_value   = display_to,
    samples    = (out_to - out_from) * scale,
    hdi        = hdi,
    std_units  = std_units,
    conditions = conditions
  )
}

# =============================================================================
# Forward propagation
# =============================================================================

#' Propagate an intervention forward through the DAG
#'
#' @param model Fitted KaguModel.
#' @param source,target Intervention and outcome node names.
#' @param value Numeric scalar — fixed value for the source node.
#' @param conditions Named list of additional fixed nodes.
#' @param n_chains,n_draws Posterior dimensions.
#' @return `[n_chains, n_draws]` matrix of target node values.
#' @noRd
.propagate <- function(model, source, target, value,
                       conditions, n_chains, n_draws) {
  order         <- topological_sort(model$dag)
  node_samples  <- list()

  for (node in order) {
    if (node == source) {
      node_samples[[node]] <- matrix(value, nrow = n_chains, ncol = n_draws)

    } else if (!is.null(conditions) && node %in% names(conditions)) {
      node_samples[[node]] <- matrix(conditions[[node]],
                                     nrow = n_chains, ncol = n_draws)

    } else {
      parents        <- model$dag[[node]]
      parent_values  <- setNames(
        lapply(parents, function(p) node_samples[[p]]),
        parents
      )
      node_samples[[node]] <- model$mechanisms[[node]]$predict_mean(
        node          = node,
        parents       = parents,
        parent_values = parent_values,
        fit           = model$traces[[node]]
      )
    }

    if (node == target) break
  }

  node_samples[[target]]
}

# =============================================================================
# Plot helpers (called by EffectResult$plot)
# =============================================================================

.plot_scalar <- function(result) {
  flat <- as.vector(result$samples)
  df   <- data.frame(effect = flat)
  h    <- bayestestR::hdi(flat, ci = result$hdi_prob)
  m    <- mean(flat)
  sd_  <- stats::sd(flat)

  x_label <- sprintf("Effect of %s on %s", result$source, result$target)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$effect)) +
    ggplot2::geom_density(fill = "#b2dfdb", colour = "#26a69a",
                          linewidth = 0.8, alpha = 0.85) +
    ggplot2::annotate("segment",
      x = h$CI_low, xend = h$CI_high, y = 0, yend = 0,
      colour = "#26a69a", linewidth = 2.5
    ) +
    ggplot2::annotate("point",
      x = m, y = 0, colour = "#26a69a", size = 3
    ) +
    ggplot2::labs(x = x_label, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.y  = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank()
    )
}

.plot_sweep <- function(result) {
  x       <- result$sweep_values
  n_sweep <- dim(result$samples)[[1]]
  # Flatten chains and draws: [n_sweep, n_chains * n_draws]
  flat    <- matrix(result$samples, nrow = n_sweep)

  means  <- rowMeans(flat)
  hdi_lo <- apply(flat, 1, function(s) bayestestR::hdi(s, ci = result$hdi_prob)$CI_low)
  hdi_hi <- apply(flat, 1, function(s) bayestestR::hdi(s, ci = result$hdi_prob)$CI_high)

  df <- data.frame(x = x, mean = means, lower = hdi_lo, upper = hdi_hi)

  x_label <- result$source
  y_label <- sprintf("E[%s | do(%s = x)]", result$target, result$source)
  hdi_pct <- sprintf("%.0f%% HDI", result$hdi_prob * 100)

  ggplot2::ggplot(df, ggplot2::aes(x = .data$x)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = .data$lower, ymax = .data$upper),
      fill = "#b2dfdb", alpha = 0.6
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data$mean),
      colour = "#26a69a", linewidth = 1
    ) +
    ggplot2::labs(x = x_label, y = y_label,
                  caption = hdi_pct) +
    ggplot2::theme_minimal(base_size = 12)
}
