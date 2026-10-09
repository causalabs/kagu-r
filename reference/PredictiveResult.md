# Prior / posterior predictive check results

Holds replicated datasets from a
[KaguModel](https://causalabs.github.io/kagu-r/reference/KaguModel.md)'s
`$prior_predictive()` or `$posterior_predictive()` alongside the
observed data. `$y()` and `$yrep()` return data in the shapes the
**bayesplot** package expects (a vector and a `[draws, observations]`
matrix), so any `bayesplot::ppc_*()` function works directly:

    pp <- model$posterior_predictive()
    bayesplot::ppc_dens_overlay(pp$y("y"), pp$yrep("y")[1:50, ])

`$summary()` reports predictive p-values for summary statistics, and
`$cor_check()` / `$plot_cor()` compare the correlation between every
pair of variables with what the model reproduces - a whole-graph check
(for `joint` and `prior` results) with no direct bayesplot equivalent.

## Public fields

- `type`:

  Character - `"prior"`, `"conditional"` or `"joint"`.

- `data`:

  The observed `data.frame`.

- `replicates`:

  Named list of `[ndraws, n]` matrices, one per node.

- `dag`:

  The model's DAG.

## Methods

### Public methods

- [`PredictiveResult$new()`](#method-PredictiveResult-initialize)

- [`PredictiveResult$nodes()`](#method-PredictiveResult-nodes)

- [`PredictiveResult$ndraws()`](#method-PredictiveResult-ndraws)

- [`PredictiveResult$y()`](#method-PredictiveResult-y)

- [`PredictiveResult$yrep()`](#method-PredictiveResult-yrep)

- [`PredictiveResult$summary()`](#method-PredictiveResult-summary)

- [`PredictiveResult$cor_check()`](#method-PredictiveResult-cor_check)

- [`PredictiveResult$plot_cor()`](#method-PredictiveResult-plot_cor)

- [`PredictiveResult$print()`](#method-PredictiveResult-print)

- [`PredictiveResult$clone()`](#method-PredictiveResult-clone)

------------------------------------------------------------------------

### `PredictiveResult$new()`

Create a PredictiveResult (normally via a
[KaguModel](https://causalabs.github.io/kagu-r/reference/KaguModel.md)).

#### Usage

    PredictiveResult$new(type, data, replicates, dag)

#### Arguments

- `type`:

  Character - `"prior"`, `"conditional"` or `"joint"`.

- `data`:

  The observed `data.frame`.

- `replicates`:

  Named list of `[ndraws, n]` matrices, one per node.

- `dag`:

  The model's DAG.

------------------------------------------------------------------------

### `PredictiveResult$nodes()`

Node names.

#### Usage

    PredictiveResult$nodes()

#### Returns

Character vector.

------------------------------------------------------------------------

### `PredictiveResult$ndraws()`

Number of replicated datasets.

#### Usage

    PredictiveResult$ndraws()

#### Returns

Integer.

------------------------------------------------------------------------

### `PredictiveResult$y()`

Observed values of a node (bayesplot's `y`).

#### Usage

    PredictiveResult$y(node)

#### Arguments

- `node`:

  Character scalar - the node name.

#### Returns

Numeric vector of length `n`.

------------------------------------------------------------------------

### `PredictiveResult$yrep()`

Replicated values of a node (bayesplot's `yrep`).

#### Usage

    PredictiveResult$yrep(node)

#### Arguments

- `node`:

  Character scalar - the node name.

#### Returns

Numeric `[ndraws, n]` matrix.

------------------------------------------------------------------------

### `PredictiveResult$summary()`

Predictive p-values for summary statistics.

For each node and statistic `T`, compares `T(y)` with the distribution
of `T(yrep)` over the replicated datasets.
`p_value = Pr(T(yrep) >= T(y))`: values near 0 or 1 mean the observed
statistic is extreme relative to what the model generates.

#### Usage

    PredictiveResult$summary(
      stats = c("mean", "sd", "min", "max"),
      nodes = NULL,
      prob = 0.9
    )

#### Arguments

- `stats`:

  Character vector of function names (e.g. `"sd"`), or a named list of
  functions each mapping a numeric vector to a scalar.

- `nodes`:

  Optional character vector of nodes (default: all).

- `prob`:

  Numeric - probability mass of the reported replicate interval.

#### Returns

A `tibble` with columns `node`, `stat`, `observed`, `rep_mean`,
`rep_lower`, `rep_upper`, `p_value`.

------------------------------------------------------------------------

### `PredictiveResult$cor_check()`

Compare observed pairwise correlations with replicated ones.

A whole-graph check: under a joint (or prior) predictive, each
replicated dataset carries the dependence structure the DAG implies, so
a pair whose observed correlation falls outside its replicated range
points to a missing (or spurious) edge or a misspecified mechanism. Not
available for `conditional` results, whose nodes are simulated
separately.

#### Usage

    PredictiveResult$cor_check(prob = 0.9, method = "pearson")

#### Arguments

- `prob`:

  Numeric - probability mass of the replicate interval.

- `method`:

  Correlation method passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html).

#### Returns

A `tibble` with columns `var1`, `var2`, `observed`, `rep_median`,
`rep_lower`, `rep_upper`, `p_value`, `outside`.

------------------------------------------------------------------------

### `PredictiveResult$plot_cor()`

Plot `$cor_check()`: each pair's replicated correlation interval with
the observed value; pairs outside their interval are highlighted.

#### Usage

    PredictiveResult$plot_cor(prob = 0.9, method = "pearson")

#### Arguments

- `prob`:

  Numeric - probability mass of the replicate interval.

- `method`:

  Correlation method passed to
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html).

#### Returns

A `ggplot` object.

------------------------------------------------------------------------

### `PredictiveResult$print()`

Print a short description.

#### Usage

    PredictiveResult$print(...)

#### Arguments

- `...`:

  Ignored.

------------------------------------------------------------------------

### `PredictiveResult$clone()`

The objects of this class are cloneable with this method.

#### Usage

    PredictiveResult$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
