# Plot a fitted node's summary terms

Shows the node's direct local effects (each parent's gradient at the
parents' means) and residual noise as posterior means with HDI
intervals.

## Usage

``` r
kagu_plot_posterior(terms, node)
```

## Arguments

- terms:

  A `tibble` of node terms (from a mechanism's `$node_terms()`), with
  columns `term`, `mean`, `hdi_lower`, `hdi_upper`.

- node:

  Character scalar - node name (used for the plot title).

## Value

A `ggplot` object.
