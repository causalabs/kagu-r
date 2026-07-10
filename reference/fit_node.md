# Fit a single node's conditional distribution

Thin dispatcher that delegates to the node's
[Mechanism](https://causalabs.github.io/kagu-r/reference/Mechanism.md).
Kept as a convenience/entry point; `KaguModel$fit()` calls it once per
node.

## Usage

``` r
fit_node(node, parents, data, mechanism, ...)
```

## Arguments

- node:

  Character scalar - name of the node to fit.

- parents:

  Character vector of parent node names (empty for root nodes).

- data:

  A `data.frame` with columns for all nodes.

- mechanism:

  A `Mechanism` instance (defaults to
  [GPMechanism](https://causalabs.github.io/kagu-r/reference/GPMechanism.md)
  in `KaguModel`).

- ...:

  Additional arguments forwarded to the mechanism's `$fit()`.

## Value

A mechanism-specific fit object.
