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

- [`KaguModel$prior_predictive()`](#method-KaguModel-prior_predictive)

- [`KaguModel$posterior_predictive()`](#method-KaguModel-posterior_predictive)

- [`KaguModel$function_draws()`](#method-KaguModel-function_draws)

- [`KaguModel$plot_dag()`](#method-KaguModel-plot_dag)

- [`KaguModel$plot_posterior()`](#method-KaguModel-plot_posterior)

- [`KaguModel$save()`](#method-KaguModel-save)

- [`KaguModel$print()`](#method-KaguModel-print)

- [`KaguModel$clone()`](#method-KaguModel-clone)

------------------------------------------------------------------------

### `KaguModel$new()`

Create a new KaguModel.

#### Usage

    KaguModel$new(dag, mechanisms = NULL, default_mechanism = NULL)

#### Arguments

- `dag`:

  Named list mapping each node to a character vector of its parent node
  names. Root nodes map to [`c()`](https://rdrr.io/r/base/c.html).

- `mechanisms`:

  Optional named list of `Mechanism` instances. Any node not specified
  receives a copy of `default_mechanism`.

- `default_mechanism`:

  Optional `Mechanism` used for every node not listed in `mechanisms`
  (each node gets its own copy). Defaults to `GPMechanism$new()`; pass
  e.g. `GPMechanism$new(prior = gp_prior(...))` to set one prior for the
  model.

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

### `KaguModel$prior_predictive()`

Simulate replicated datasets from the **prior** predictive distribution,
by ancestral sampling through the DAG: each draw samples every node's
hyperparameters, function and noise from its mechanism's prior, feeding
simulated parents into their children. Needs no fit.

#### Usage

    KaguModel$prior_predictive(data = NULL, ndraws = 100L)

#### Arguments

- `data`:

  A `data.frame`, used only for each variable's location and scale
  (priors are on the standardised scale) and the number of rows.
  Defaults to the fitted data.

- `ndraws`:

  Integer - number of replicated datasets.

#### Returns

A
[PredictiveResult](https://causalabs.github.io/kagu-r/reference/PredictiveResult.md)
with `type = "prior"`.

------------------------------------------------------------------------

### `KaguModel$posterior_predictive()`

Simulate replicated datasets from the **posterior** predictive
distribution.

- `type = "conditional"`: every node is simulated given its *observed*
  parents. This checks each mechanism (local fit) in isolation.

- `type = "joint"`: the whole graph is simulated ancestrally from the
  fitted model, with each node's posterior draw fed its simulated
  parents. This checks what the DAG implies globally - e.g. a missing
  edge shows up as a dependence the replicated data cannot reproduce.

Hyperparameters are point estimates, so the replicates do not carry
hyperparameter uncertainty and can be slightly too narrow.

#### Usage

    KaguModel$posterior_predictive(ndraws = 100L, type = c("conditional", "joint"))

#### Arguments

- `ndraws`:

  Integer - number of replicated datasets (at most the number of stored
  posterior draws).

- `type`:

  `"conditional"` or `"joint"`.

#### Returns

A
[PredictiveResult](https://causalabs.github.io/kagu-r/reference/PredictiveResult.md).

------------------------------------------------------------------------

### `KaguModel$function_draws()`

Prior or posterior draws of a node's conditional-mean function against
one parent, with the other parents held fixed - for plotting what
functions the prior allows or the posterior has learned.

#### Usage

    KaguModel$function_draws(
      node,
      parent = NULL,
      ndraws = 50L,
      n_grid = 100L,
      prior = FALSE,
      at = NULL,
      data = NULL
    )

#### Arguments

- `node`:

  Character scalar - the node whose function to draw.

- `parent`:

  Character scalar - the parent to vary (default: the first).

- `ndraws`:

  Integer - number of function draws.

- `n_grid`:

  Integer - number of grid points across the parent's range.

- `prior`:

  Logical - draw from the prior (`TRUE`) or posterior.

- `at`:

  Optional named list of values for the other parents (default: their
  means).

- `data`:

  Optional `data.frame` (defaults to the fitted data).

#### Returns

A `tibble` with columns `draw`, `parent_value`, `value`.

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
