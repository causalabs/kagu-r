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
#' @param required Optional list of length-2 character vectors `c(from, to)`,
#'   each requiring the directed edge `from -> to` to be present in all candidate DAGs.
#' @param allow_empty Logical — whether to include the completely empty (edgeless)
#'   graph in the search space. Defaults to `FALSE`.
#' @return A list of DAG specifications. Over 3 unconstrained nodes there are 24
#'   (with `allow_empty = FALSE`).
#' @export
#'
#' @examples
#' length(enumerate_dags(c("a", "b", "c")))            # 24
#' # forbid b -> a and c -> a (e.g. 'a' comes first in time)
#' enumerate_dags(c("a", "b", "c"),
#'                disallowed = list(c("b", "a"), c("c", "a")))
enumerate_dags <- function(nodes, disallowed = NULL, required = NULL, allow_empty = FALSE) {
  if (length(nodes) < 1L) stop("`nodes` must contain at least one node.")
  if (anyDuplicated(nodes)) stop("`nodes` must be unique.")

  # All candidate directed edges (ordered pairs), minus the disallowed ones.
  grid  <- expand.grid(from = nodes, to = nodes, stringsAsFactors = FALSE)
  edges <- grid[grid$from != grid$to, , drop = FALSE]

  .fmt <- function(el) vapply(el, function(e) paste(e[[1]], e[[2]], sep = "\r"), character(1))
  edge_strs <- paste(edges$from, edges$to, sep = "\r")

  if (!is.null(disallowed)) {
    dis  <- .fmt(disallowed)
    keep <- !(edge_strs %in% dis)
    edges <- edges[keep, , drop = FALSE]
    edge_strs <- edge_strs[keep]
  }

  req_edges <- data.frame(from = character(0), to = character(0), stringsAsFactors = FALSE)
  if (!is.null(required)) {
    req <- .fmt(required)
    is_req <- edge_strs %in% req
    req_edges <- edges[is_req, , drop = FALSE]
    edges <- edges[!is_req, , drop = FALSE]
    
    if (!.is_acyclic(.edges_to_dag(req_edges, nodes))) {
       stop("The `required` edges contain a cycle.")
    }
  }

  m <- nrow(edges)
  if (m > 24L) {
    cli::cli_abort(c(
      "The search space is too large: {m} candidate 'maybe' directed edges (2^{m} subsets).",
      "i" = "Add more {.arg disallowed} or {.arg required} edges, or reduce {.arg nodes}."
    ))
  }

  # Each subset of candidate edges that is acyclic is exactly one DAG.
  powers <- bitwShiftL(1L, seq_len(m) - 1L)
  dags   <- vector("list", 0L)
  
  max_mask <- if (m == 0L) 0L else 2L^m - 1L
  for (mask in seq.int(0L, max_mask)) {
    if (!allow_empty && mask == 0L && nrow(req_edges) == 0L && (m > 0L || length(nodes) > 1L)) next
    
    if (m > 0L) {
      incl <- bitwAnd(as.integer(mask), powers) != 0L
      comb_edges <- rbind(req_edges, edges[incl, , drop = FALSE])
    } else {
      comb_edges <- req_edges
    }
    
    dag  <- .edges_to_dag(comb_edges, nodes)
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
  nodes <- names(dag)
  n_nodes <- length(nodes)
  
  # Map nodes to integer indices for speed
  node_idx <- setNames(seq_len(n_nodes), nodes)
  
  in_degree <- integer(n_nodes)
  # Pre-allocate adj list
  adj <- vector("list", n_nodes)
  
  for (i in seq_len(n_nodes)) {
    parents <- dag[[i]]
    in_degree[i] <- length(parents)
    for (pa in parents) {
      pa_idx <- node_idx[[pa]]
      adj[[pa_idx]] <- c(adj[[pa_idx]], i)
    }
  }
  
  queue <- which(in_degree == 0L)
  visited <- 0L
  
  while (length(queue) > 0L) {
    u <- queue[1]
    queue <- queue[-1]
    visited <- visited + 1L
    
    for (v in adj[[u]]) {
      in_degree[v] <- in_degree[v] - 1L
      if (in_degree[v] == 0L) {
        queue <- c(queue, v)
      }
    }
  }
  
  visited == n_nodes
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
