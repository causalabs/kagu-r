# Build a per-node summary table across all nodes

For each node, the node's mechanism contributes summary rows - for a
Gaussian process, the direct local effect of each parent (the fitted
function's gradient at the parents' means, comparable to a regression
coefficient) and the residual noise sd.

## Usage

``` r
build_summary_table(model)
```

## Arguments

- model:

  A fitted `KaguModel`.

## Value

A `tibble` with columns `node`, `term`, `mean`, `sd`, `hdi_lower`,
`hdi_upper`.
