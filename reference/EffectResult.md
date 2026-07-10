# Result of a causal effect query

Returned by `KaguModel$effects()`. Holds the full posterior over the
causal effect along with metadata.

- **Scalar effect**: `samples` is a `[n_chains, n_draws]` matrix.

- **Sweep effect**: `samples` is a `[n_sweep, n_chains, n_draws]` array.

## Public fields

- `source`:

  Character - the intervention node.

- `target`:

  Character - the outcome node.

- `from_value`:

  Numeric - intervention baseline value (or sweep grid).

- `to_value`:

  Numeric - intervention target value (or sweep grid).

- `samples`:

  Numeric array - posterior effect samples.

- `hdi_prob`:

  Numeric - HDI probability used for summaries and plots.

- `std_units`:

  Logical - whether the effect is in SD units.

- `conditions`:

  Named list of conditioned node values, or `NULL`.

- `sweep_values`:

  Numeric vector of sweep grid values, or `NULL`.

- `conditions_label`:

  String representation of conditions (for plotting).

- `compare_results`:

  List of other `EffectResult`s (for multi-condition plots).

## Methods

### Public methods

- [`EffectResult$new()`](#method-EffectResult-initialize)

- [`EffectResult$is_sweep()`](#method-EffectResult-is_sweep)

- [`EffectResult$summary()`](#method-EffectResult-summary)

- [`EffectResult$plot()`](#method-EffectResult-plot)

- [`EffectResult$diagnostics()`](#method-EffectResult-diagnostics)

- [`EffectResult$print()`](#method-EffectResult-print)

- [`EffectResult$clone()`](#method-EffectResult-clone)

------------------------------------------------------------------------

### `EffectResult$new()`

Create an EffectResult (normally called by `compute_effect`).

#### Usage

    EffectResult$new(
      source,
      target,
      from_value,
      to_value,
      samples,
      hdi = 0.9,
      std_units = FALSE,
      conditions = NULL,
      sweep_values = NULL
    )

#### Arguments

- `source`:

  Character - the intervention node.

- `target`:

  Character - the outcome node.

- `from_value`:

  Numeric - intervention baseline value.

- `to_value`:

  Numeric - intervention target value.

- `samples`:

  Numeric array - posterior effect samples.

- `hdi`:

  Numeric - HDI probability.

- `std_units`:

  Logical - whether the effect is in SD units.

- `conditions`:

  Named list of conditioned node values.

- `sweep_values`:

  Numeric vector of sweep grid values.

------------------------------------------------------------------------

### `EffectResult$is_sweep()`

Is this a sweep result?

#### Usage

    EffectResult$is_sweep()

#### Returns

Logical scalar.

------------------------------------------------------------------------

### `EffectResult$summary()`

Tabular summary of the effect posterior.

Returns a `tibble`, consistent with `KaguModel$summary()`. The HDI
columns use the same `hdi_lower` / `hdi_upper` naming as
[`build_summary_table()`](https://causalabs.github.io/kagu-r/reference/build_summary_table.md).

- **Scalar**: a one-row tibble with `source`, `target`, `from`, `to`,
  `mean`, `sd`, `hdi_lower`, `hdi_upper`.

- **Sweep**: one row per grid point with `source`, `target`, `x`,
  `mean`, `sd`, `hdi_lower`, `hdi_upper`.

#### Usage

    EffectResult$summary()

#### Returns

A `tibble`.

------------------------------------------------------------------------

### `EffectResult$plot()`

Plot the posterior effect distribution.

- **Scalar**: density plot with mean point and HDI bar.

- **Sweep**: dose-response line with HDI ribbon.

#### Usage

    EffectResult$plot(compare = NULL, base_label = NULL)

#### Arguments

- `compare`:

  Optional named list of other `EffectResult` instances to overlay
  (useful for comparing sweeps under different interactions/conditions).
  Only valid if this result and all comparison results are
  `sweep = TRUE`. If this effect was generated with a list of
  conditions, `compare` is populated automatically.

- `base_label`:

  Optional character string to label this `EffectResult` in the legend
  when plotting with `compare` (defaults to `"Base"` or the condition
  label).

#### Returns

A `ggplot` object.

------------------------------------------------------------------------

### `EffectResult$diagnostics()`

R-hat and ESS diagnostics for the effect posterior.

Chain structure is preserved, so r-hat values are genuine.

#### Usage

    EffectResult$diagnostics()

#### Returns

A `tibble` from
[`posterior::summarise_draws()`](https://mc-stan.org/posterior/reference/draws_summary.html).

------------------------------------------------------------------------

### `EffectResult$print()`

Print method - shows the summary table.

Print summary of the effect result.

#### Usage

    EffectResult$print(...)

#### Arguments

- `...`:

  Ignored.

------------------------------------------------------------------------

### `EffectResult$clone()`

The objects of this class are cloneable with this method.

#### Usage

    EffectResult$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
