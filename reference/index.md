# Package index

## Model

Main user-facing class.

- [`KaguModel`](https://causalabs.github.io/kagu-r/reference/KaguModel.md)
  : KaguModel - Bayesian graphical causal model

## Mechanisms

Probabilistic models for node conditional distributions.

- [`Mechanism`](https://causalabs.github.io/kagu-r/reference/Mechanism.md)
  : Mechanism abstract base class
- [`GPMechanism`](https://causalabs.github.io/kagu-r/reference/GPMechanism.md)
  : Gaussian-process mechanism (default)

## Effects

Causal effect estimation and results.

- [`EffectResult`](https://causalabs.github.io/kagu-r/reference/EffectResult.md)
  : Result of a causal effect query
- [`compute_effect()`](https://causalabs.github.io/kagu-r/reference/compute_effect.md)
  : Compute the causal effect of source on target

## Structure discovery

Posterior over DAGs via marginal-likelihood model comparison.

- [`kagu_discover()`](https://causalabs.github.io/kagu-r/reference/kagu_discover.md)
  : Discover causal structure from data
- [`DiscoveryResult`](https://causalabs.github.io/kagu-r/reference/DiscoveryResult.md)
  : Result of a causal structure search
- [`enumerate_dags()`](https://causalabs.github.io/kagu-r/reference/enumerate_dags.md)
  : Enumerate all DAGs over a set of nodes
- [`kagu_plot_discovery()`](https://causalabs.github.io/kagu-r/reference/kagu_plot_discovery.md)
  : Plot the posterior distribution over DAGs

## DAG utilities

Graph operations on DAG specifications.

- [`validate_dag()`](https://causalabs.github.io/kagu-r/reference/validate_dag.md)
  : Validate a DAG specification
- [`topological_sort()`](https://causalabs.github.io/kagu-r/reference/topological_sort.md)
  : Topological sort (Kahn's algorithm)
- [`ancestors()`](https://causalabs.github.io/kagu-r/reference/ancestors.md)
  : Ancestors of a node
- [`descendants()`](https://causalabs.github.io/kagu-r/reference/descendants.md)
  : Descendants of a node
- [`is_ancestor()`](https://causalabs.github.io/kagu-r/reference/is_ancestor.md)
  : Test whether one node is an ancestor of another
- [`node_depth()`](https://causalabs.github.io/kagu-r/reference/node_depth.md)
  : Depth of each node in the DAG

## Plots

- [`kagu_plot_dag()`](https://causalabs.github.io/kagu-r/reference/kagu_plot_dag.md)
  : Plot the DAG structure
- [`kagu_plot_posterior()`](https://causalabs.github.io/kagu-r/reference/kagu_plot_posterior.md)
  : Plot a fitted node's summary terms

## Save / load

- [`kagu_save()`](https://causalabs.github.io/kagu-r/reference/kagu_save.md)
  : Save a fitted model to disk

- [`kagu_load()`](https://causalabs.github.io/kagu-r/reference/kagu_load.md)
  :

  Load a model saved with
  [`kagu_save()`](https://causalabs.github.io/kagu-r/reference/kagu_save.md)

## Fitting

- [`fit_node()`](https://causalabs.github.io/kagu-r/reference/fit_node.md)
  : Fit a single node's conditional distribution

## Summaries

- [`build_summary_table()`](https://causalabs.github.io/kagu-r/reference/build_summary_table.md)
  : Build a per-node summary table across all nodes

## Internal

File-level documentation pages.

- [`dag`](https://causalabs.github.io/kagu-r/reference/dag.md) : DAG
  validation and graph utilities
- [`discover`](https://causalabs.github.io/kagu-r/reference/discover.md)
  : Causal structure discovery
- [`effects`](https://causalabs.github.io/kagu-r/reference/effects.md) :
  Causal effect estimation via the do-operator
- [`inference`](https://causalabs.github.io/kagu-r/reference/inference.md)
  : Per-node model fitting
- [`io`](https://causalabs.github.io/kagu-r/reference/io.md) : Save and
  load fitted KaguModel objects
- [`mechanisms`](https://causalabs.github.io/kagu-r/reference/mechanisms.md)
  : Mechanism base class and the Gaussian-process implementation
- [`plots`](https://causalabs.github.io/kagu-r/reference/plots.md) : DAG
  and posterior visualisation
- [`summary`](https://causalabs.github.io/kagu-r/reference/summary.md) :
  Summary table for a fitted KaguModel
