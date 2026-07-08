#' DAG and posterior visualisation
#'
#' @name plots
NULL

#' Plot the DAG structure
#'
#' Produces a `ggplot` showing nodes and directed edges. Layout is computed
#' automatically from topological depth; individual nodes can be repositioned
#' via `node_pos`.
#'
#' @param dag Named list — the DAG specification.
#' @param node_pos Optional named list of `c(row, col)` coordinates for manual
#'   node placement. Partial overrides are supported — unspecified nodes keep
#'   the auto-layout position.
#' @return A `ggplot` object.
#' @export
kagu_plot_dag <- function(dag, node_pos = NULL) {
  nodes  <- names(dag)
  depth  <- node_depth(dag)

  # --- Auto-layout: group nodes by depth layer, spread within each layer -----
  layers     <- split(nodes, depth[nodes])
  positions  <- list()

  for (d in names(layers)) {
    layer_nodes <- layers[[d]]
    n           <- length(layer_nodes)
    col_offset  <- seq_len(n) - (n + 1L) / 2.0
    for (i in seq_along(layer_nodes)) {
      positions[[layer_nodes[[i]]]] <- c(row = as.numeric(d), col = col_offset[[i]])
    }
  }

  # --- Apply manual overrides ------------------------------------------------
  if (!is.null(node_pos)) {
    for (n in names(node_pos)) {
      positions[[n]] <- c(row = node_pos[[n]][[1]], col = node_pos[[n]][[2]])
    }
  }

  # --- Build node data frame -------------------------------------------------
  node_df <- data.frame(
    name = nodes,
    x    = sapply(nodes, function(n) positions[[n]][["col"]]),
    y    = sapply(nodes, function(n) -positions[[n]][["row"]]),  # flip y so depth 0 is top
    stringsAsFactors = FALSE
  )

  # --- Build edge data frame -------------------------------------------------
  edges <- do.call(rbind, lapply(nodes, function(child) {
    parents <- dag[[child]]
    if (length(parents) == 0) return(NULL)
    data.frame(
      from   = parents,
      to     = rep(child, length(parents)),
      x      = sapply(parents, function(p) node_df$x[node_df$name == p]),
      y      = sapply(parents, function(p) node_df$y[node_df$name == p]),
      xend   = node_df$x[node_df$name == child],
      yend   = node_df$y[node_df$name == child],
      stringsAsFactors = FALSE
    )
  }))

  # Shorten edges so they start/end at the node circle boundary rather than
  # its centre (leaves a clean gap for the arrowhead, Pearl-diagram style).
  if (!is.null(edges) && nrow(edges) > 0) {
    SHRINK <- 0.13
    dx     <- edges$xend - edges$x
    dy     <- edges$yend - edges$y
    len    <- sqrt(dx^2 + dy^2)
    edges$x    <- edges$x    + SHRINK * dx / len
    edges$y    <- edges$y    + SHRINK * dy / len
    edges$xend <- edges$xend - SHRINK * dx / len
    edges$yend <- edges$yend - SHRINK * dy / len
  }

  # --- Plot (black & white, literature / Pearl style) ------------------------
  p <- ggplot2::ggplot(node_df, ggplot2::aes(x = .data$x, y = .data$y))

  if (!is.null(edges) && nrow(edges) > 0) {
    p <- p + ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = .data$x, y = .data$y,
                   xend = .data$xend, yend = .data$yend),
      arrow     = grid::arrow(length = grid::unit(6, "pt"), type = "closed"),
      colour    = "black",
      linewidth = 0.5
    )
  }

  p +
    # open node: white fill masks any crossing edge, thick black outline
    ggplot2::geom_point(shape = 21, size = 9, stroke = 1.2,
                        fill = "white", colour = "black") +
    # variable name set just below the node, on a translucent rounded white
    # label so it stays legible where an edge passes behind it
    ggplot2::geom_label(ggplot2::aes(label = .data$name),
                        nudge_y       = -0.17,
                        size          = 3.6,
                        colour        = "black",
                        fill          = "#FFFFFFB3",
                        label.size    = 0,
                        label.r       = grid::unit(0.12, "lines"),
                        label.padding = grid::unit(0.12, "lines")) +
    # generous expansion + clip = "off" so nodes/labels are never cut off
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.18)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.22)) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(14, 14, 14, 14))
}

#' Plot the posterior distribution over DAGs
#'
#' Ranked bar chart of the posterior probability of each candidate DAG from a
#' structure search. If the true (data-generating) DAG is supplied, its bar is
#' highlighted.
#'
#' @param result A [DiscoveryResult] (from `KaguModel$discover()`).
#' @param top_n Integer — number of top-ranked DAGs to display (default 20).
#' @param true_dag Optional DAG specification to highlight (matched by its set
#'   of directed edges).
#' @return A `ggplot` object.
#' @export
kagu_plot_discovery <- function(result, top_n = 20L, true_dag = NULL) {
  ord <- order(result$prob, decreasing = TRUE)
  ord <- ord[seq_len(min(top_n, length(ord)))]

  df <- data.frame(
    rank      = seq_along(ord),
    id        = result$labels[ord],
    prob      = result$prob[ord],
    highlight = FALSE
  )
  df$id <- factor(df$id, levels = df$id)   # keep posterior order on the axis

  if (!is.null(true_dag)) {
    true_edges  <- .dag_edges(true_dag)
    df$highlight <- vapply(ord, function(i) {
      identical(.dag_edges(result$dags[[i]]), true_edges)
    }, logical(1))
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$id, y = .data$prob,
                                        fill = .data$highlight)) +
    ggplot2::geom_col(width = 0.85, colour = "black", linewidth = 0.3) +
    ggplot2::scale_fill_manual(
      values = c(`FALSE` = "grey85", `TRUE` = "#26a69a"), guide = "none"
    ) +
    ggplot2::labs(x = "DAG (id, ranked by posterior)", y = "P(G | data)") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  if (!is.null(true_dag) && any(df$highlight)) {
    hl <- df[df$highlight, , drop = FALSE]
    p <- p + ggplot2::annotate(
      "text", x = hl$rank[[1]], y = hl$prob[[1]],
      label = "true DAG", vjust = -0.6, size = 3.4, colour = "#1f8e83"
    )
  }
  p
}

#' Plot a fitted node's summary terms
#'
#' Shows the node's direct local effects (each parent's gradient at the parents'
#' means) and residual noise as posterior means with HDI intervals.
#'
#' @param terms A `tibble` of node terms (from a mechanism's `$node_terms()`),
#'   with columns `term`, `mean`, `hdi_lower`, `hdi_upper`.
#' @param node Character scalar — node name (used for the plot title).
#' @return A `ggplot` object.
#' @export
kagu_plot_posterior <- function(terms, node) {
  terms$term <- factor(terms$term, levels = rev(terms$term))

  ggplot2::ggplot(terms, ggplot2::aes(y = .data$term)) +
    ggplot2::geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.4) +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = .data$hdi_lower, xmax = .data$hdi_upper),
      colour = "#26a69a", linewidth = 1.2
    ) +
    ggplot2::geom_point(ggplot2::aes(x = .data$mean), size = 2.4,
                        colour = "#26a69a") +
    ggplot2::labs(
      title = sprintf("Node '%s': direct local effects", node),
      x = "posterior (mean and 90% HDI)", y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank())
}
