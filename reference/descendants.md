# Descendants of a node

Descendants of a node

## Usage

``` r
descendants(dag, node)
```

## Arguments

- dag:

  Named list as described in
  [`validate_dag()`](https://causalabs.github.io/kagu-r/reference/validate_dag.md).

- node:

  Character scalar - the node to query.

## Value

Character vector of all nodes reachable from `node` following directed
edges (i.e. the node's causal descendants).
