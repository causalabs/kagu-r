#' Summary table for a fitted KaguModel
#'
#' @name summary
NULL

#' Build a parameter summary table across all nodes
#'
#' @param traces Named list of `brmsfit` objects, one per node.
#' @param hdi_prob Numeric — HDI probability (default 0.90).
#' @return A `tibble` with columns `node`, `variable`, `mean`, `sd`,
#'   `hdi_lower`, `hdi_upper`, `rhat`, `ess_bulk`.
#' @export
build_summary_table <- function(traces, hdi_prob = 0.90) {
  results <- lapply(names(traces), function(node) {
    summ <- posterior::summarise_draws(
      traces[[node]],
      mean = mean,
      sd   = stats::sd,
      # `summarise_draws()` hands the lambda a per-variable draws_array (2-D:
      # iterations x chains). bayestestR::hdi() dispatches to its draws method
      # which assumes a 3-D array and errors, so flatten to a plain vector.
      ~ bayestestR::hdi(as.numeric(.), ci = hdi_prob),
      rhat     = posterior::rhat,
      ess_bulk = posterior::ess_bulk
    )

    # Keep only population-level and sigma parameters
    summ <- summ[grepl("^b_|^sigma", summ$variable), ]

    df       <- as.data.frame(summ)
    df$node  <- node

    # bayestestR::hdi() adds CI / CI_low / CI_high; drop the redundant CI level
    # column and rename the bounds to the documented hdi_lower / hdi_upper.
    df$CI <- NULL
    if ("CI_low" %in% names(df))  names(df)[names(df) == "CI_low"]  <- "hdi_lower"
    if ("CI_high" %in% names(df)) names(df)[names(df) == "CI_high"] <- "hdi_upper"

    df[, c("node", "variable", "mean", "sd",
           "hdi_lower", "hdi_upper", "rhat", "ess_bulk")]
  })

  tibble::as_tibble(do.call(rbind, results))
}
