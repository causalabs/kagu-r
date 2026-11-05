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
  hdi = 0.9,
  noise = TRUE,
  n_mc = 50L
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

- noise:

  Logical - how exogenous variation is handled during propagation.

  - `TRUE` (default): the interventional mean
    `E[target | do(source = x)]` over the observed population. Every
    observed row is a unit: root nodes keep their observed values and
    every other node gets its posterior observation noise. For each
    posterior draw the target is averaged over all units, so the noise
    is integrated out rather than added to the effect's uncertainty, and
    the uncertainty is about the mechanisms only (the data are treated
    as the whole population, not a sample from a larger one). The same
    units are reused for both arms of a contrast and for every sweep
    point (common random numbers), so Monte Carlo error largely cancels.

  - `FALSE`: the noise-free plug-in. Root nodes are fixed at their mean
    and each node is set to its conditional mean. Identical to `TRUE`
    when every mechanism on the path is linear; otherwise it is the
    effect in a system with the exogenous variation switched off. See
    [`vignette("noise_propagation")`](https://causalabs.github.io/kagu-r/articles/noise_propagation.md).

- n_mc:

  Integer - minimum number of simulated units per posterior draw when
  `noise = TRUE` (default 50). The observed rows are replicated to reach
  it, and only when a noisy node lies between the roots and the target;
  with more rows than `n_mc` it has no effect.

## Value

An `EffectResult`.
