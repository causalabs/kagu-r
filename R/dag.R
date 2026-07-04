#' DAG validation and graph utilities
#'
#' Pure-R implementations of the graph operations needed for causal inference.
#' No external graph library required.
#'
#' @name dag
NULL

#' Validate a DAG specification
#'
#' @param dag Named list mapping each node name to a character vector of its
#'   parent node names. Root nodes map to `character(0)` or `c()`.
#' @return Invisible `NULL`. Raises an error if the DAG is invalid.
#' @export
#'
#' @examples
#' dag <- list(x = c(), y = c("x"), z = c("x", "y"))
#' validate_dag(dag)  # silent — DAG is valid
validate_dag <- function(dag) {
  nodes <- names(dag)
  if (is.null(nodes) || length(nodes) == 0) {
    stop("DAG must have at least one node.")
  }

  for (node in nodes) {
    parents <- dag[[node]]
    unknown <- setdiff(parents, nodes)
    if (length(unknown) > 0) {
      stop(sprintf(
        "Node '%s' references unknown parent(s): %s",
        node, paste(unknown, collapse = ", ")
      ))
    }
  }

  # Cycle detection via topological sort — will error if cycle found
  tryCatch(
    topological_sort(dag),
    error = function(e) stop("DAG contains a cycle.")
  )

  invisible(NULL)
}

#' Topological sort (Kahn's algorithm)
#'
#' @param dag Named list as described in [validate_dag()].
#' @return Character vector of node names in a valid topological order
#'   (parents always before children).
#' @export
topological_sort <- function(dag) {
  nodes <- names(dag)

  # Build children map and in-degree count
  in_degree <- setNames(integer(length(nodes)), nodes)
  children  <- setNames(vector("list", length(nodes)), nodes)

  for (node in nodes) {
    for (parent in dag[[node]]) {
      in_degree[[node]] <- in_degree[[node]] + 1L
      children[[parent]] <- c(children[[parent]], node)
    }
  }

  queue  <- nodes[in_degree == 0L]
  sorted <- character(0)

  while (length(queue) > 0) {
    node   <- queue[[1]]
    queue  <- queue[-1]
    sorted <- c(sorted, node)

    for (child in children[[node]]) {
      in_degree[[child]] <- in_degree[[child]] - 1L
      if (in_degree[[child]] == 0L) {
        queue <- c(queue, child)
      }
    }
  }

  if (length(sorted) != length(nodes)) {
    stop("DAG contains a cycle.")
  }

  sorted
}

#' Descendants of a node
#'
#' @param dag Named list as described in [validate_dag()].
#' @param node Character scalar — the node to query.
#' @return Character vector of all nodes reachable from `node` following
#'   directed edges (i.e. the node's causal descendants).
#' @export
descendants <- function(dag, node) {
  children <- .children_map(dag)
  result   <- character(0)
  queue    <- children[[node]]

  while (length(queue) > 0) {
    current <- queue[[1]]
    queue   <- queue[-1]
    if (!(current %in% result)) {
      result <- c(result, current)
      queue  <- c(queue, children[[current]])
    }
  }
  result
}

#' Ancestors of a node
#'
#' @param dag Named list as described in [validate_dag()].
#' @param node Character scalar — the node to query.
#' @return Character vector of all nodes with a directed path to `node`
#'   (i.e. the node's causal ancestors).
#' @export
ancestors <- function(dag, node) {
  result <- character(0)
  queue  <- dag[[node]]

  while (length(queue) > 0) {
    current <- queue[[1]]
    queue   <- queue[-1]
    if (!(current %in% result)) {
      result <- c(result, current)
      queue  <- c(queue, dag[[current]])
    }
  }
  result
}

#' Test whether one node is an ancestor of another
#'
#' @param dag Named list as described in [validate_dag()].
#' @param source Character scalar — the potential ancestor.
#' @param target Character scalar — the potential descendant.
#' @return Logical scalar.
#' @export
is_ancestor <- function(dag, source, target) {
  source %in% ancestors(dag, target)
}

#' Depth of each node in the DAG
#'
#' Depth is defined as the length of the longest path from any root node.
#' Root nodes have depth 0.
#'
#' @param dag Named list as described in [validate_dag()].
#' @return Named integer vector of node depths.
#' @export
node_depth <- function(dag) {
  depth <- setNames(integer(length(dag)), names(dag))
  for (node in topological_sort(dag)) {
    if (length(dag[[node]]) > 0) {
      depth[[node]] <- max(depth[dag[[node]]]) + 1L
    }
  }
  depth
}

#' Enumerate all DAGs over a set of nodes
#'
#' Generates every directed acyclic graph over `nodes`, optionally excluding a
#' set of disallowed directed edges. Each returned DAG is in the named-list
#' format used throughout the package (node -> character vector of parents).
#'
#' The number of DAGs grows super-exponentially in the number of nodes, so the
#' search space is only tractable for small problems. Use `disallowed` (e.g. to
#' encode a known temporal ordering) to prune it.
#'
#' @param nodes Character vector of node names.
#' @param disallowed Optional list of length-2 character vectors `c(from, to)`,
#'   each forbidding the directed edge `from -> to`. Directions are independent:
#'   forbidding `a -> b` still permits `b -> a`.
#' @return A list of DAG specifications. Over 3 unconstrained nodes there are 25.
#' @export
#'
#' @examples
#' length(enumerate_dags(c("a", "b", "c")))            # 25
#' # forbid b -> a and c -> a (e.g. 'a' comes first in time)
#' enumerate_dags(c("a", "b", "c"),
#'                disallowed = list(c("b", "a"), c("c", "a")))
enumerate_dags <- function(nodes, disallowed = NULL) {
  if (length(nodes) < 1L) stop("`nodes` must contain at least one node.")
  if (anyDuplicated(nodes)) stop("`nodes` must be unique.")

  # All candidate directed edges (ordered pairs), minus the disallowed ones.
  grid  <- expand.grid(from = nodes, to = nodes, stringsAsFactors = FALSE)
  edges <- grid[grid$from != grid$to, , drop = FALSE]

  if (!is.null(disallowed)) {
    dis  <- vapply(disallowed, function(e) paste(e[[1]], e[[2]], sep = "\r"),
                   character(1))
    keep <- !(paste(edges$from, edges$to, sep = "\r") %in% dis)
    edges <- edges[keep, , drop = FALSE]
  }

  m <- nrow(edges)
  if (m > 20L) {
    cli::cli_abort(c(
      "The search space is too large: {m} candidate directed edges (2^{m} subsets).",
      "i" = "Add more {.arg disallowed} edges (e.g. a temporal ordering) or
             reduce {.arg nodes}."
    ))
  }

  # Each subset of candidate edges that is acyclic is exactly one DAG.
  powers <- bitwShiftL(1L, seq_len(m) - 1L)
  dags   <- vector("list", 0L)
  for (mask in seq.int(0L, 2L^m - 1L)) {
    incl <- bitwAnd(as.integer(mask), powers) != 0L
    dag  <- .edges_to_dag(edges[incl, , drop = FALSE], nodes)
    if (.is_acyclic(dag)) dags[[length(dags) + 1L]] <- dag
  }
  dags
}

# --- Internal helpers --------------------------------------------------------

#' Build a DAG (named parent-list) from an edge data.frame
#' @noRd
.edges_to_dag <- function(edges, nodes) {
  setNames(
    lapply(nodes, function(n) as.character(edges$from[edges$to == n])),
    nodes
  )
}

#' Is a DAG specification acyclic? (reuses Kahn's algorithm)
#' @noRd
.is_acyclic <- function(dag) {
  tryCatch({ topological_sort(dag); TRUE }, error = function(e) FALSE)
}

.children_map <- function(dag) {
  nodes    <- names(dag)
  children <- setNames(vector("list", length(nodes)), nodes)
  for (node in nodes) {
    for (parent in dag[[node]]) {
      children[[parent]] <- c(children[[parent]], node)
    }
  }
  children
}
