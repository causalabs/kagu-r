# Topological sort (Kahn's algorithm)

Topological sort (Kahn's algorithm)

## Usage

``` r
topological_sort(dag)
```

## Arguments

- dag:

  Named list as described in
  [`validate_dag()`](https://causalabs.github.io/kagu-r/reference/validate_dag.md).

## Value

Character vector of node names in a valid topological order (parents
always before children).
