# Ancestors of a node

Ancestors of a node

## Usage

``` r
ancestors(dag, node)
```

## Arguments

- dag:

  Named list as described in
  [`validate_dag()`](https://causalabs.github.io/kagu-r/reference/validate_dag.md).

- node:

  Character scalar - the node to query.

## Value

Character vector of all nodes with a directed path to `node` (i.e. the
node's causal ancestors).
