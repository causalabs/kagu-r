# Compute the causal effect of source on target

Called internally by `KaguModel$effects()`. Not usually called directly.

## Usage

``` r
compute_effect(
  model,
  source,
  target,
  values = NULL,
  std_units = FALSE,
  conditions = NULL,
  sweep = FALSE,
  sweep_n = 50L,
  sweep_range = NULL,
  hdi = 0.9
)
```

## Arguments

- model:

  A fitted `KaguModel`.

- source, target:

  Character scalars - intervention and outcome nodes.

- values:

  Optional `c(from, to)` intervention contrast.

- std_units:

  Logical - 1-SD effect.

- conditions:

  Optional named list of fixed node values.

- sweep:

  Logical - compute dose-response curve.

- sweep_n, sweep_range:

  Sweep grid parameters.

- hdi:

  Numeric - HDI probability.

## Value

An `EffectResult`.
