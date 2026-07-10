# Depth of each node in the DAG

Depth is defined as the length of the longest path from any root node.
Root nodes have depth 0.

## Usage

``` r
node_depth(dag)
```

## Arguments

- dag:

  Named list as described in
  [`validate_dag()`](https://causalabs.github.io/kagu-r/reference/validate_dag.md).

## Value

Named integer vector of node depths.
