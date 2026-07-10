# Test whether one node is an ancestor of another

Test whether one node is an ancestor of another

## Usage

``` r
is_ancestor(dag, source, target)
```

## Arguments

- dag:

  Named list as described in
  [`validate_dag()`](https://causalabs.github.io/kagu-r/reference/validate_dag.md).

- source:

  Character scalar - the potential ancestor.

- target:

  Character scalar - the potential descendant.

## Value

Logical scalar.
