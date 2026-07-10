# KaguModel - Bayesian graphical causal model

The main user-facing class. Specify a DAG, optionally assign mechanisms
to nodes, fit to data, and extract causal effects.

## Public fields

- `dag`:

  Named list - the DAG specification.

- `mechanisms`:

  Named list of `Mechanism` instances, one per node.

- `data`:

  The training `data.frame` (set after `$fit()`).

- `traces`:

  Named list of mechanism fit objects, one per node (after fit).

- `.fitted`:

  Logical - whether `$fit()` has been called.

## Methods

### Public methods

- [`KaguModel$new()`](#method-KaguModel-initialize)

- [`KaguModel$fit()`](#method-KaguModel-fit)

- [`KaguModel$effects()`](#method-KaguModel-effects)

- [`KaguModel$summary()`](#method-KaguModel-summary)

- [`KaguModel$diagnostics()`](#method-KaguModel-diagnostics)

- [`KaguModel$plot_dag()`](#method-KaguModel-plot_dag)

- [`KaguModel$plot_posterior()`](#method-KaguModel-plot_posterior)

- [`KaguModel$save()`](#method-KaguModel-save)

- [`KaguModel$print()`](#method-KaguModel-print)

- [`KaguModel$clone()`](#method-KaguModel-clone)

------------------------------------------------------------------------

### `KaguModel$new()`

Create a new KaguModel.

#### Usage

    KaguModel$new(dag, mechanisms = NULL)

#### Arguments

- `dag`:

  Named list mapping each node to a character vector of its parent node
  names. Root nodes map to [`c()`](https://rdrr.io/r/base/c.html).

- `mechanisms`:

  Optional named list of `Mechanism` instances. Any node not specified
  receives a `GPMechanism` by default.

------------------------------------------------------------------------

### `KaguModel$fit()`

Fit each node's conditional distribution in topological order.

#### Usage

    KaguModel$fit(data, ...)

#### Arguments

- `data`:

  A `data.frame` with one column per node.

- `...`:

  Additional arguments forwarded to each node's mechanism `$fit()`.

#### Returns

`self` invisibly (for method chaining).

------------------------------------------------------------------------

### `KaguModel$effects()`

Estimate the causal effect of `source` on `target`.

#### Usage

    KaguModel$effects(
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

#### Arguments

- `source`:

  Character scalar - the intervention (treatment) node.

- `target`:

  Character scalar - the outcome node.

- `values`:

  Optional numeric vector of length 2 `c(from, to)` giving the explicit
  intervention contrast.

- `std_units`:

  Logical - if `TRUE`, compute the effect of a 1-SD increase centred at
  the mean.

- `conditions`:

  Optional named list of node values to condition on (fixes those nodes
  at the given values during propagation). Can also be a list of such
  lists to compute and overlay multiple conditions.

- `sweep`:

  Logical - if `TRUE`, compute the dose-response curve.

- `sweep_n`:

  Integer - number of points in the sweep grid (default 50).

- `sweep_range`:

  Numeric vector `c(min, max)` for the sweep grid. Defaults to the
  observed range of `source`.

- `hdi`:

  Numeric in (0, 1) - HDI probability (default 0.90).

#### Returns

An `EffectResult` object.

------------------------------------------------------------------------

### `KaguModel$summary()`

Per-node summary table.

For each node, reports the direct local effect of each parent (the
function's gradient at the parents' means - comparable to a regression
coefficient) and the residual noise sd, each with posterior mean, sd and
HDI.

#### Usage

    KaguModel$summary(hdi_prob = 0.9)

#### Arguments

- `hdi_prob`:

  Numeric - HDI probability (default 0.90; currently fixed).

#### Returns

A `tibble` with columns `node`, `term`, `mean`, `sd`, `hdi_lower`,
`hdi_upper`.

------------------------------------------------------------------------

### `KaguModel$diagnostics()`

Fit diagnostics for each node.

The Gaussian-process sampler produces a single chain, so r-hat / ESS do
not apply; this reports the posterior draw count and the residual noise
sd per node.

#### Usage

    KaguModel$diagnostics(node = NULL)

#### Arguments

- `node`:

  Optional character scalar. If `NULL`, runs for all nodes.

#### Returns

A `tibble` with `node`, `n_chains`, `n_draws`, `sigma`, `sigma_sd`.

------------------------------------------------------------------------

### `KaguModel$plot_dag()`

Visualise the DAG structure.

#### Usage

    KaguModel$plot_dag(node_pos = NULL)

#### Arguments

- `node_pos`:

  Optional named list of `c(row, col)` grid coordinates for manual node
  placement (overrides auto-layout for specified nodes).

#### Returns

A `ggplot` object.

------------------------------------------------------------------------

### `KaguModel$plot_posterior()`

Plot the posterior for a fitted node.

Shows the node's direct local effects (each parent's gradient at the
means) and residual noise, as posterior means with HDI intervals.

#### Usage

    KaguModel$plot_posterior(node)

#### Arguments

- `node`:

  Character scalar - the node to plot.

#### Returns

A `ggplot` object.

------------------------------------------------------------------------

### `KaguModel$save()`

Save the fitted model to disk.

#### Usage

    KaguModel$save(path, include_data = TRUE)

#### Arguments

- `path`:

  File path (e.g. `"model.rds"`).

- `include_data`:

  Logical - whether to save the training data alongside the model
  (default `TRUE`, so `$effects()` works immediately on load).

------------------------------------------------------------------------

### `KaguModel$print()`

Print a concise model summary.

Print summary of the model.

#### Usage

    KaguModel$print(...)

#### Arguments

- `...`:

  Ignored.

------------------------------------------------------------------------

### `KaguModel$clone()`

The objects of this class are cloneable with this method.

#### Usage

    KaguModel$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
if (FALSE) { # \dontrun{
library(kagu)

model <- KaguModel$new(
  dag = list(
    age     = c(),
    smoking = c("age"),
    health  = c("smoking", "age")
  )
)

model$fit(data)

effect <- model$effects("smoking", "health")
effect$summary()
} # }
```
