#' Causal structure discovery
#'
#' @name discover
NULL

# =============================================================================
# DiscoveryResult R6 class
# =============================================================================

#' Result of a causal structure search
#'
#' @description
#' Returned by `KaguModel$discover()` (equivalently [kagu_discover()]). Holds the
#' posterior probability distribution over candidate DAGs, together with the
#' per-graph log marginal likelihoods and the cache of unique local fits.
#'
#' @export
DiscoveryResult <- R6::R6Class("DiscoveryResult",
  public = list(
    #' @field dags List of candidate DAG specifications evaluated.
    dags = NULL,
    #' @field labels Character vector of short ids (one per DAG) — user-supplied
    #'   names if the `dags` list was named, otherwise spreadsheet-style ids
    #'   (`"A"`, `"B"`, ..., `"AA"`). Used to reference DAGs in summaries and in
    #'   `$get_dag()` / `$plot_dag()`.
    labels = NULL,
    #' @field nodes Character vector of node names.
    nodes = NULL,
    #' @field data The data used for the search (for refitting via `$as_model()`).
    data = NULL,
    #' @field log_marglik Numeric vector — log marginal likelihood per DAG.
    log_marglik = NULL,
    #' @field log_prior Numeric vector — log prior per DAG.
    log_prior = NULL,
    #' @field log_posterior Numeric vector — normalised log posterior per DAG.
    log_posterior = NULL,
    #' @field prob Numeric vector — posterior probability per DAG (sums to 1).
    prob = NULL,
    #' @field local_logml Named numeric — cached log marginal likelihood of each
    #'   unique local model, keyed by `"node|sorted,parents"`.
    local_logml = NULL,
    #' @field local_fits Named list — the fitted model for each unique local
    #'   model (reused for effect propagation; same keys as `local_logml`).
    local_fits = NULL,
    #' @field mechanisms Named list of `Mechanism` instances, one per node.
    mechanisms = NULL,
    #' @field n_models Integer — number of candidate DAGs.
    n_models = NULL,
    #' @field n_unique_fits Integer — number of unique local models fitted.
    n_unique_fits = NULL,

    #' @description Create a DiscoveryResult (normally called by `kagu_discover`).
    #' @param dags,labels,nodes,data,log_marglik,log_prior,log_posterior,prob,local_logml,local_fits,mechanisms,n_unique_fits
    #'   Internal fields (see the corresponding `$field` documentation).
    initialize = function(dags, labels, nodes, data, log_marglik, log_prior,
                          log_posterior, prob, local_logml, local_fits,
                          mechanisms, n_unique_fits) {
      self$dags          <- dags
      self$labels        <- labels
      self$nodes         <- nodes
      self$data          <- data
      self$log_marglik   <- log_marglik
      self$log_prior     <- log_prior
      self$log_posterior <- log_posterior
      self$prob          <- prob
      self$local_logml   <- local_logml
      self$local_fits    <- local_fits
      self$mechanisms    <- mechanisms
      self$n_models      <- length(dags)
      self$n_unique_fits <- n_unique_fits
    },

    #' @description Ranked summary table of the most probable DAGs.
    #'
    #' Each DAG is referenced by its short `id` (see `$labels`) rather than by a
    #' full edge listing — inspect a structure with `$get_dag(id)` or
    #' `$plot_dag(id)`.
    #' @param top_n Integer — number of top DAGs to return (default 10).
    #' @return A `tibble` with `rank`, `id`, `n_edges`, `log_marglik`,
    #'   `posterior_prob`.
    summary = function(top_n = 10L) {
      ord <- order(self$prob, decreasing = TRUE)
      ord <- ord[seq_len(min(top_n, length(ord)))]
      tibble::tibble(
        rank           = seq_along(ord),
        id             = self$labels[ord],
        n_edges        = vapply(ord, function(i) length(.dag_edges(self$dags[[i]])),
                                integer(1)),
        log_marglik    = self$log_marglik[ord],
        posterior_prob = self$prob[ord]
      )
    },

    #' @description Posterior probability of every candidate DAG.
    #' @return A `tibble` with `rank`, `id`, `n_edges`, `log_marglik`,
    #'   `posterior_prob` (all DAGs).
    probabilities = function() self$summary(top_n = self$n_models),

    #' @description Look up a candidate DAG by its `id` (see `$labels`).
    #' @param id Character scalar — a DAG id from `$summary()` / `$labels`.
    #' @return A DAG specification (named list of parent vectors).
    get_dag = function(id) {
      i <- match(id, self$labels)
      if (is.na(i)) {
        stop(sprintf("No DAG with id '%s'. Available ids: %s",
                     id, paste(self$labels, collapse = ", ")))
      }
      self$dags[[i]]
    },

    #' @description Plot a candidate DAG by its `id` (see `$labels`).
    #' @param id Character scalar — a DAG id from `$summary()` / `$labels`.
    #' @param node_pos Optional named list of `c(row, col)` node positions,
    #'   passed to [kagu_plot_dag()].
    #' @return A `ggplot` object.
    plot_dag = function(id, node_pos = NULL) {
      kagu_plot_dag(self$get_dag(id), node_pos = node_pos)
    },

    #' @description Marginal posterior probability of each directed edge.
    #'
    #' For each directed edge, sums the posterior probability of all DAGs that
    #' contain it: `P(from -> to | data) = sum_{G: from->to in G} P(G | data)`.
    #'
    #' @return A `tibble` with `from`, `to`, `prob`, sorted by `prob`.
    edge_probabilities = function() {
      acc <- list()
      for (i in seq_along(self$dags)) {
        for (e in .dag_edges(self$dags[[i]])) {
          acc[[e]] <- (if (is.null(acc[[e]])) 0 else acc[[e]]) + self$prob[[i]]
        }
      }
      if (!length(acc)) {
        return(tibble::tibble(from = character(0), to = character(0),
                              prob = numeric(0)))
      }
      parts <- strsplit(names(acc), "→", fixed = TRUE)
      out <- tibble::tibble(
        from = vapply(parts, `[[`, character(1), 1L),
        to   = vapply(parts, `[[`, character(1), 2L),
        prob = unlist(acc, use.names = FALSE)
      )
      out[order(out$prob, decreasing = TRUE), ]
    },

    #' @description Estimate a causal effect, averaged over the DAG posterior.
    #'
    #' Computes the **Bayesian model-averaged** causal effect,
    #' \eqn{P(Q \mid X) = \sum_G P(Q \mid X, G)\, P(G \mid X)}, by combining each
    #' DAG's effect posterior weighted by its posterior probability. This reflects
    #' structural uncertainty as well as parameter uncertainty, and is preferred
    #' over committing to a single DAG (e.g. the MAP) when the structure is not
    #' certain.
    #'
    #' The signature mirrors `KaguModel$effects()` and the return value is the
    #' same [EffectResult] type, so summaries, plots and downstream code are
    #' identical to the single-DAG workflow. DAGs in which `source` has no causal
    #' path to `target` contribute a zero effect, exactly as they should.
    #'
    #' @param source,target Character scalars — intervention and outcome nodes.
    #' @param values,std_units,conditions,sweep,sweep_n,sweep_range,hdi Passed to
    #'   `KaguModel$effects()` (identical meaning).
    #' @param prob_threshold Numeric — DAGs with posterior probability below this
    #'   are skipped for efficiency (default 1e-3).
    #' @param n_samples Integer — size of the pooled mixture posterior (default 4000).
    #' @return An [EffectResult].
    effects = function(source, target, values = NULL, std_units = FALSE,
                       conditions = NULL, sweep = FALSE, sweep_n = 50L,
                       sweep_range = NULL, hdi = 0.90,
                       prob_threshold = 1e-3, n_samples = 4000L) {
      if (!source %in% self$nodes) stop(sprintf("Node '%s' not in nodes.", source))
      if (!target %in% self$nodes) stop(sprintf("Node '%s' not in nodes.", target))

      ord  <- order(self$prob, decreasing = TRUE)
      keep <- ord[self$prob[ord] >= prob_threshold]

      ers      <- vector("list", length(keep))
      weights  <- numeric(length(keep))
      is_anc   <- logical(length(keep))
      for (j in seq_along(keep)) {
        i   <- keep[[j]]
        dag <- self$dags[[i]]
        model         <- KaguModel$new(dag, mechanisms = self$mechanisms[names(dag)])
        model$data    <- self$data
        model$traces  <- stats::setNames(
          lapply(names(dag), function(n) self$local_fits[[.local_key(n, dag[[n]])]]),
          names(dag)
        )
        model$.fitted <- TRUE
        ers[[j]]    <- model$effects(source, target, values = values,
                                     std_units = std_units, conditions = conditions,
                                     sweep = sweep, sweep_n = sweep_n,
                                     sweep_range = sweep_range, hdi = hdi)
        weights[[j]] <- self$prob[[i]]
        is_anc[[j]]  <- is_ancestor(dag, source, target)
      }

      .pool_effects(ers, weights / sum(weights), is_anc, n_samples)
    },

    #' @description The maximum a posteriori (most probable) DAG.
    #'
    #' Provided for inspection of the single most probable structure. For
    #' estimating causal effects under structural uncertainty, prefer
    #' `$effects()`, which averages over the whole posterior.
    #' @return A DAG specification (named list of parent vectors).
    map = function() self$dags[[which.max(self$prob)]],

    #' @description Alias for `$map()`.
    best = function() self$map(),

    #' @description Refit a chosen DAG as a full [KaguModel].
    #' @param which Either `"map"` (default) for the most probable DAG, or a DAG
    #'   specification to fit.
    #' @param ... Passed to `KaguModel$fit()`.
    #' @return A fitted `KaguModel`.
    as_model = function(which = "map", ...) {
      dag   <- if (identical(which, "map")) self$map() else which
      model <- KaguModel$new(dag)
      model$fit(self$data, ...)
      model
    },

    #' @description Plot the posterior over DAGs.
    #' @param top_n Integer — number of top DAGs to show (default 20).
    #' @param true_dag Optional DAG specification to highlight (e.g. the known
    #'   data-generating structure).
    #' @return A `ggplot` object.
    plot = function(top_n = 20L, true_dag = NULL) {
      kagu_plot_discovery(self, top_n = top_n, true_dag = true_dag)
    },

    #' @description Print method.
    #' @description Print summary of the discovery result.
    #' @param ... Ignored.
    print = function(...) {
      cat(sprintf("<DiscoveryResult: %d DAGs, %d unique local fits>\n",
                  self$n_models, self$n_unique_fits))
      print(self$summary(top_n = 5L))
      invisible(self)
    }
  )
)

# =============================================================================
# kagu_discover — the engine behind KaguModel$discover()
# =============================================================================

#' Discover causal structure from data
#'
#' Fits a GCM to every candidate DAG over the supplied variables and returns a
#' posterior probability distribution over those DAGs,
#' \deqn{P(G \mid X) \propto P(X \mid G)\, P(G),}
#' where each graph's marginal likelihood \eqn{P(X \mid G)} is obtained from the
#' per-node Gaussian-process fits and \eqn{P(G)} is the prior.
#'
#' This is normally called as the static method **`KaguModel$discover(data,
#' ...)`**; `kagu_discover()` is the underlying function.
#'
#' Two properties keep the computation tractable. First, by the Markov property
#' a DAG's marginal likelihood factorises over nodes,
#' \eqn{\log P(X \mid G) = \sum_i \log P(X_i \mid \mathrm{Pa}_G(X_i))}, so a
#' graph's score is a sum of independent per-node terms. Second, the same
#' `(node, parent-set)` local model recurs across many DAGs, so each unique
#' local model is fitted only once and cached — collapsing the work from
#' `#DAGs` fits to at most `p * 2^(p-1)`.
#'
#' @section Priors:
#' Only a uniform prior over DAGs is currently supported (`prior = "uniform"`),
#' so the posterior is driven entirely by the marginal likelihoods. The `prior`
#' argument is reserved for future options — e.g. sparsity-favouring priors,
#' edge/temporal constraints, or fully custom priors over structures.
#'
#' @section Marginal likelihoods and Markov equivalence:
#' Each node's marginal likelihood is the closed-form type-II (empirical-Bayes)
#' evidence of its Gaussian-process model, so discovery involves no bridge
#' sampling and no Stan compilation. Note that Markov-equivalent DAGs are
#' statistically indistinguishable from observational data and will therefore
#' receive (near-)equal posterior mass — a faithful representation of structural
#' uncertainty, not a defect.
#'
#' @param data A `data.frame` with one column per variable.
#' @param nodes Optional character vector of node names. Defaults to all columns
#'   of `data`.
#' @param dags Optional explicit list of candidate DAGs to score (each a named
#'   list of parent vectors). If provided, `disallowed` and `required` are ignored,
#'   and the search space is restricted exactly to these DAGs.
#' @param disallowed Optional list of length-2 character vectors `c(from, to)`,
#'   each forbidding the directed edge `from -> to`. Useful for encoding a known
#'   temporal ordering and for pruning the search space.
#' @param required Optional list of length-2 character vectors `c(from, to)`,
#'   each requiring the directed edge `from -> to` to be present in all candidate DAGs.
#' @param mechanisms Optional named list of [Mechanism] instances, one per node.
#'   Any node not specified receives a [GPMechanism].
#' @param prior Character — prior over DAGs. Only `"uniform"` is supported.
#' @param allow_empty Logical — whether to include the completely edgeless graph
#'   in the search space (default `FALSE`). Usually, researchers are looking
#'   for at least some causal structure, so the empty graph is omitted.
#' @param ... Additional arguments forwarded to each node's mechanism `$fit()`.
#' @return A [DiscoveryResult].
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 300
#' a <- rnorm(n)
#' b <- 0.8 * a + rnorm(n, sd = 0.5)
#' c <- 1.2 * b + rnorm(n, sd = 0.5)
#' df <- data.frame(a = a, b = b, c = c)
#'
#' # Equivalent to KaguModel$discover(df)
#' res <- kagu_discover(df)
#' res$summary()
#' res$edge_probabilities()
#' res$plot(true_dag = list(a = c(), b = "a", c = "b"))
#' }
kagu_discover <- function(data, nodes = NULL, dags = NULL, disallowed = NULL,
                          required = NULL, mechanisms = NULL, prior = "uniform",
                          allow_empty = FALSE, ...) {
  if (is.null(nodes)) nodes <- names(data)

  missing <- setdiff(nodes, names(data))
  if (length(missing)) {
    cli::cli_abort(
      "{cli::qty(missing)}Node{?s} {.val {missing}} {?is/are} not {?a column/columns} in {.arg data}."
    )
  }
  prior <- match.arg(prior, "uniform")

  if (is.null(mechanisms)) mechanisms <- list()
  mech_for <- function(n) {
    if (!is.null(mechanisms[[n]])) mechanisms[[n]] else GPMechanism$new()
  }
  resolved_mech <- stats::setNames(lapply(nodes, mech_for), nodes)

  # ---- Enumerate candidate DAGs --------------------------------------------
  if (!is.null(dags)) {
    for (i in seq_along(dags)) {
      if (!.is_acyclic(dags[[i]])) stop(sprintf("Provided DAG %d contains a cycle.", i))
      missing_dag_nodes <- setdiff(nodes, names(dags[[i]]))
      if (length(missing_dag_nodes)) stop(sprintf("Provided DAG %d is missing nodes: %s", i, paste(missing_dag_nodes, collapse=", ")))
    }
  } else {
    dags <- enumerate_dags(nodes, disallowed = disallowed, required = required, allow_empty = allow_empty)
  }
  n_models <- length(dags)

  # Stable short id per DAG: user-supplied names if the list was fully named
  # (unique, non-empty), otherwise spreadsheet-style ids "A", "B", ..., "AA".
  labels <- names(dags)
  if (is.null(labels) || any(!nzchar(labels)) || anyDuplicated(labels)) {
    labels <- vapply(seq_len(n_models), .letter_id, character(1))
  }

  # ---- Collect unique local models (node, sorted parent set) ---------------
  registry <- list()
  for (dag in dags) {
    for (n in nodes) {
      k <- .local_key(n, dag[[n]])
      if (is.null(registry[[k]])) registry[[k]] <- list(node = n, parents = dag[[n]])
    }
  }
  keys     <- names(registry)
  n_unique <- length(keys)

  # ---- Fit each unique local model once; keep both fit and logML -----------
  cli::cli_alert_info(
    "Scoring {n_models} DAG{?s} via {n_unique} unique local fit{?s} …"
  )
  local_logml <- stats::setNames(numeric(n_unique), keys)
  local_fits  <- stats::setNames(vector("list", n_unique), keys)
  t_start <- Sys.time()
  for (i in seq_along(keys)) {
    info <- registry[[keys[[i]]]]
    rhs  <- if (length(info$parents)) paste(info$parents, collapse = " + ") else "1"
    # Progress with a rolling time estimate.
    eta <- ""
    if (i > 1L) {
      per <- as.numeric(difftime(Sys.time(), t_start, units = "secs")) / (i - 1L)
      eta <- sprintf(" — ~%s remaining", .fmt_duration(per * (n_unique - i + 1L)))
    }
    cli::cli_alert("[{i}/{n_unique}] fitting {.field {info$node}} ~ {rhs}{eta}")

    res <- .fit_and_marglik(
      info$node, info$parents, data, mech_for(info$node), ...
    )
    local_fits[[keys[[i]]]]  <- res$fit
    local_logml[[keys[[i]]]] <- res$logml
  }
  if (anyNA(local_logml)) {
    warning("Some local fits failed; DAGs depending on them get probability 0.")
  }

  # ---- Combine via Markov factorisation, add prior, normalise --------------
  # Pre-calculate the logml for each unique local fit
  
  # For each DAG, extract the logml for each node, sum them, and return
  log_marglik <- vapply(dags, function(dag) {
    # .local_key(n, dag[[n]]) is fast, but we can just map and sum
    s <- 0
    for (n in nodes) {
      val <- local_logml[[.local_key(n, dag[[n]])]]
      if (is.na(val)) return(-Inf)
      s <- s + val
    }
    s
  }, numeric(1))

  log_prior <- rep(-log(n_models), n_models)          # uniform
  log_unnorm <- log_marglik + log_prior
  log_z      <- .logsumexp(log_unnorm)
  if (!is.finite(log_z)) {
    stop("All candidate DAGs have non-finite marginal likelihood; cannot normalise.")
  }
  log_post <- log_unnorm - log_z
  prob     <- exp(log_post)

  cli::cli_alert_success("Done — posterior over {n_models} DAG{?s}.")

  DiscoveryResult$new(
    dags = dags, labels = labels, nodes = nodes, data = data,
    log_marglik = log_marglik, log_prior = log_prior,
    log_posterior = log_post, prob = prob,
    local_logml = local_logml, local_fits = local_fits,
    mechanisms = resolved_mech, n_unique_fits = n_unique
  )
}

# =============================================================================
# Internal helpers
# =============================================================================

#' Column-letter id for an index: 1 -> "A", 26 -> "Z", 27 -> "AA"
#' @noRd
.letter_id <- function(i) {
  out <- ""
  while (i > 0) {
    r   <- (i - 1L) %% 26L
    out <- paste0(LETTERS[[r + 1L]], out)
    i   <- (i - 1L) %/% 26L
  }
  out
}

#' Cache key for a local model: "node|sorted,parents"
#' @noRd
.local_key <- function(node, parents) {
  paste0(node, "|", paste(sort(parents), collapse = ","))
}

#' Canonical (sorted) set of "parent -> child" edge labels for a DAG
#' @noRd
.dag_edges <- function(dag) {
  e <- character(0)
  for (child in names(dag)) {
    for (parent in dag[[child]]) {
      e <- c(e, paste0(parent, "→", child))
    }
  }
  sort(e)
}

#' Human-readable edge string, e.g. "a->c; b->c" (or "(empty graph)")
#' @noRd
.format_edges <- function(dag) {
  e <- .dag_edges(dag)
  if (length(e)) paste(e, collapse = "; ") else "(empty graph)"
}

#' Numerically stable log-sum-exp (ignoring non-finite entries)
#' @noRd
.logsumexp <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m)))
}

#' Human-readable duration, e.g. "45s" or "3m 05s"
#' @noRd
.fmt_duration <- function(secs) {
  secs <- max(0, round(secs))
  if (secs < 60) return(sprintf("%ds", secs))
  sprintf("%dm %02ds", secs %/% 60, secs %% 60)
}

#' Draw a weighted mixture sample from a list of numeric vectors
#'
#' Allocates `N` draws across the components in proportion to `weights` and
#' samples (with replacement) from each — the empirical realisation of a mixture
#' distribution. Used to combine per-DAG effect posteriors into one.
#' @noRd
.mixture_sample <- function(vecs, weights, N) {
  counts <- as.integer(round(N * weights))
  counts[[1]] <- counts[[1]] + (N - sum(counts))   # absorb rounding error
  out <- unlist(lapply(seq_along(vecs), function(j) {
    if (counts[[j]] <= 0 || !length(vecs[[j]])) return(numeric(0))
    sample(vecs[[j]], counts[[j]], replace = TRUE)
  }), use.names = FALSE)
  out
}

#' Pool a list of EffectResults into a single Bayesian model-averaged EffectResult
#'
#' Mixes the per-DAG effect posteriors with the given (normalised) weights into
#' one [EffectResult] of the same shape, so the result is indistinguishable in
#' type from a single-DAG effect query. Display metadata (`from`/`to`) is taken
#' from the most probable DAG in which `source` actually causes `target`.
#' @noRd
.pool_effects <- function(ers, weights, is_anc, n_samples) {
  # Representative for from/to display: most probable ancestor DAG (ers are
  # already in descending-probability order), else the first.
  rep_idx <- if (any(is_anc)) which(is_anc)[[1]] else 1L
  rep_er  <- ers[[rep_idx]]

  if (rep_er$is_sweep()) {
    n_sweep <- dim(rep_er$samples)[[1]]
    pooled  <- array(NA_real_, dim = c(n_sweep, 1L, n_samples))
    for (s in seq_len(n_sweep)) {
      vecs <- lapply(ers, function(r) as.vector(r$samples[s, , ]))
      pooled[s, 1, ] <- .mixture_sample(vecs, weights, n_samples)
    }
  } else {
    vecs   <- lapply(ers, function(r) as.vector(r$samples))
    pooled <- matrix(.mixture_sample(vecs, weights, n_samples), nrow = 1L)
  }

  EffectResult$new(
    source = rep_er$source, target = rep_er$target,
    from_value = rep_er$from_value, to_value = rep_er$to_value,
    samples = pooled, hdi = rep_er$hdi_prob, std_units = rep_er$std_units,
    conditions = rep_er$conditions, sweep_values = rep_er$sweep_values
  )
}
