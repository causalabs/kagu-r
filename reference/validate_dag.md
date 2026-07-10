# Validate a DAG specification

Validate a DAG specification

## Usage

``` r
validate_dag(dag)
```

## Arguments

- dag:

  Named list mapping each node name to a character vector of its parent
  node names. Root nodes map to `character(0)` or
  [`c()`](https://rdrr.io/r/base/c.html).

## Value

Invisible `NULL`. Raises an error if the DAG is invalid.

## Examples

``` r
dag <- list(x = c(), y = c("x"), z = c("x", "y"))
validate_dag(dag)  # silent - DAG is valid
```
