# Plot the posterior distribution over DAGs

Ranked bar chart of the posterior probability of each candidate DAG from
a structure search. If the true (data-generating) DAG is supplied, its
bar is highlighted.

## Usage

``` r
kagu_plot_discovery(result, top_n = 20L, true_dag = NULL)
```

## Arguments

- result:

  A
  [DiscoveryResult](https://causalabs.github.io/kagu-r/reference/DiscoveryResult.md)
  (from `KaguModel$discover()`).

- top_n:

  Integer - number of top-ranked DAGs to display (default 20).

- true_dag:

  Optional DAG specification to highlight (matched by its set of
  directed edges).

## Value

A `ggplot` object.
