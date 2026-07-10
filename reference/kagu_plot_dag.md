# Plot the DAG structure

Produces a `ggplot` showing nodes and directed edges. Layout is computed
automatically from topological depth; individual nodes can be
repositioned via `node_pos`.

## Usage

``` r
kagu_plot_dag(dag, node_pos = NULL)
```

## Arguments

- dag:

  Named list - the DAG specification.

- node_pos:

  Optional named list of `c(row, col)` coordinates for manual node
  placement. Partial overrides are supported - unspecified nodes keep
  the auto-layout position.

## Value

A `ggplot` object.
