# Result of a causal structure search

Returned by `KaguModel$discover()` (equivalently
[`kagu_discover()`](https://causalabs.github.io/kagu-r/reference/kagu_discover.md)).
Holds the posterior probability distribution over candidate DAGs,
together with the per-graph log marginal likelihoods and the cache of
unique local fits.

## Public fields

- `dags`:

  List of candidate DAG specifications evaluated.

- `labels`:

  Character vector of short ids (one per DAG) - user-supplied names if
  the `dags` list was named, otherwise spreadsheet-style ids (`"A"`,
  `"B"`, ..., `"AA"`). Used to reference DAGs in summaries and in
  `$get_dag()` / `$plot_dag()`.

- `nodes`:

  Character vector of node names.

- `data`:

  The data used for the search (for refitting via `$as_model()`).

- `log_marglik`:

  Numeric vector - log marginal likelihood per DAG.

- `log_prior`:

  Numeric vector - log prior per DAG.

- `log_posterior`:

  Numeric vector - normalised log posterior per DAG.

- `prob`:

  Numeric vector - posterior probability per DAG (sums to 1).

- `local_logml`:

  Named numeric - cached log marginal likelihood of each unique local
  model, keyed by `"node|sorted,parents"`.

- `local_fits`:

  Named list - the fitted model for each unique local model (reused for
  effect propagation; same keys as `local_logml`).

- `mechanisms`:

  Named list of `Mechanism` instances, one per node.

- `n_models`:

  Integer - number of candidate DAGs.

- `n_unique_fits`:

  Integer - number of unique local models fitted.

## Methods

### Public methods

- [`DiscoveryResult$new()`](#method-DiscoveryResult-initialize)

- [`DiscoveryResult$summary()`](#method-DiscoveryResult-summary)

- [`DiscoveryResult$probabilities()`](#method-DiscoveryResult-probabilities)

- [`DiscoveryResult$get_dag()`](#method-DiscoveryResult-get_dag)

- [`DiscoveryResult$plot_dag()`](#method-DiscoveryResult-plot_dag)

- [`DiscoveryResult$edge_probabilities()`](#method-DiscoveryResult-edge_probabilities)

- [`DiscoveryResult$effects()`](#method-DiscoveryResult-effects)

- [`DiscoveryResult$map()`](#method-DiscoveryResult-map)

- [`DiscoveryResult$best()`](#method-DiscoveryResult-best)

- [`DiscoveryResult$as_model()`](#method-DiscoveryResult-as_model)

- [`DiscoveryResult$plot()`](#method-DiscoveryResult-plot)

- [`DiscoveryResult$print()`](#method-DiscoveryResult-print)

- [`DiscoveryResult$clone()`](#method-DiscoveryResult-clone)

------------------------------------------------------------------------

### `DiscoveryResult$new()`

Create a DiscoveryResult (normally called by `kagu_discover`).

#### Usage

    DiscoveryResult$new(
      dags,
      labels,
      nodes,
      data,
      log_marglik,
      log_prior,
      log_posterior,
      prob,
      local_logml,
      local_fits,
      mechanisms,
      n_unique_fits
    )

#### Arguments

- `dags, labels, nodes, data, log_marglik, log_prior, log_posterior, prob, local_logml, local_fits, mechanisms, n_unique_fits`:

  Internal fields (see the corresponding `$field` documentation).

------------------------------------------------------------------------

### `DiscoveryResult$summary()`

Ranked summary table of the most probable DAGs.

Each DAG is referenced by its short `id` (see `$labels`) rather than by
a full edge listing - inspect a structure with `$get_dag(id)` or
`$plot_dag(id)`.

#### Usage

    DiscoveryResult$summary(top_n = 10L)

#### Arguments

- `top_n`:

  Integer - number of top DAGs to return (default 10).

#### Returns

A `tibble` with `rank`, `id`, `n_edges`, `log_marglik`,
`posterior_prob`.

------------------------------------------------------------------------

### `DiscoveryResult$probabilities()`

Posterior probability of every candidate DAG.

#### Usage

    DiscoveryResult$probabilities()

#### Returns

A `tibble` with `rank`, `id`, `n_edges`, `log_marglik`, `posterior_prob`
(all DAGs).

------------------------------------------------------------------------

### `DiscoveryResult$get_dag()`

Look up a candidate DAG by its `id` (see `$labels`).

#### Usage

    DiscoveryResult$get_dag(id)

#### Arguments

- `id`:

  Character scalar - a DAG id from `$summary()` / `$labels`.

#### Returns

A DAG specification (named list of parent vectors).

------------------------------------------------------------------------

### `DiscoveryResult$plot_dag()`

Plot a candidate DAG by its `id` (see `$labels`).

#### Usage

    DiscoveryResult$plot_dag(id, node_pos = NULL)

#### Arguments

- `id`:

  Character scalar - a DAG id from `$summary()` / `$labels`.

- `node_pos`:

  Optional named list of `c(row, col)` node positions, passed to
  [`kagu_plot_dag()`](https://causalabs.github.io/kagu-r/reference/kagu_plot_dag.md).

#### Returns

A `ggplot` object.

------------------------------------------------------------------------

### `DiscoveryResult$edge_probabilities()`

Marginal posterior probability of each directed edge.

For each directed edge, sums the posterior probability of all DAGs that
contain it: `P(from -> to | data) = sum_{G: from->to in G} P(G | data)`.

#### Usage

    DiscoveryResult$edge_probabilities()

#### Returns

A `tibble` with `from`, `to`, `prob`, sorted by `prob`.

------------------------------------------------------------------------

### `DiscoveryResult$effects()`

Estimate a causal effect, averaged over the DAG posterior.

Computes the **Bayesian model-averaged** causal effect, \\P(Q \mid X) =
\sum_G P(Q \mid X, G)\\ P(G \mid X)\\, by combining each DAG's effect
posterior weighted by its posterior probability. This reflects
structural uncertainty as well as parameter uncertainty, and is
preferred over committing to a single DAG (e.g. the MAP) when the
structure is not certain.

The signature mirrors `KaguModel$effects()` and the return value is the
same
[EffectResult](https://causalabs.github.io/kagu-r/reference/EffectResult.md)
type, so summaries, plots and downstream code are identical to the
single-DAG workflow. DAGs in which `source` has no causal path to
`target` contribute a zero effect, exactly as they should.

#### Usage

    DiscoveryResult$effects(
      source,
      target,
      values = NULL,
      std_units = FALSE,
      conditions = NULL,
      sweep = FALSE,
      sweep_n = 50L,
      sweep_range = NULL,
      hdi = 0.9,
      prob_threshold = 0.001,
      n_samples = 4000L
    )

#### Arguments

- `source, target`:

  Character scalars - intervention and outcome nodes.

- `values, std_units, conditions, sweep, sweep_n, sweep_range, hdi`:

  Passed to `KaguModel$effects()` (identical meaning).

- `prob_threshold`:

  Numeric - DAGs with posterior probability below this are skipped for
  efficiency (default 1e-3).

- `n_samples`:

  Integer - size of the pooled mixture posterior (default 4000).

#### Returns

An
[EffectResult](https://causalabs.github.io/kagu-r/reference/EffectResult.md).

------------------------------------------------------------------------

### `DiscoveryResult$map()`

The maximum a posteriori (most probable) DAG.

Provided for inspection of the single most probable structure. For
estimating causal effects under structural uncertainty, prefer
`$effects()`, which averages over the whole posterior.

#### Usage

    DiscoveryResult$map()

#### Returns

A DAG specification (named list of parent vectors).

------------------------------------------------------------------------

### `DiscoveryResult$best()`

Alias for `$map()`.

#### Usage

    DiscoveryResult$best()

------------------------------------------------------------------------

### `DiscoveryResult$as_model()`

Refit a chosen DAG as a full
[KaguModel](https://causalabs.github.io/kagu-r/reference/KaguModel.md).

#### Usage

    DiscoveryResult$as_model(which = "map", ...)

#### Arguments

- `which`:

  Either `"map"` (default) for the most probable DAG, or a DAG
  specification to fit.

- `...`:

  Passed to `KaguModel$fit()`.

#### Returns

A fitted `KaguModel`.

------------------------------------------------------------------------

### `DiscoveryResult$plot()`

Plot the posterior over DAGs.

#### Usage

    DiscoveryResult$plot(top_n = 20L, true_dag = NULL)

#### Arguments

- `top_n`:

  Integer - number of top DAGs to show (default 20).

- `true_dag`:

  Optional DAG specification to highlight (e.g. the known
  data-generating structure).

#### Returns

A `ggplot` object.

------------------------------------------------------------------------

### `DiscoveryResult$print()`

Print method.

Print summary of the discovery result.

#### Usage

    DiscoveryResult$print(...)

#### Arguments

- `...`:

  Ignored.

------------------------------------------------------------------------

### `DiscoveryResult$clone()`

The objects of this class are cloneable with this method.

#### Usage

    DiscoveryResult$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
