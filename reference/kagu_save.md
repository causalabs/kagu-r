# Save a fitted model to disk

Save a fitted model to disk

## Usage

``` r
kagu_save(model, path, include_data = TRUE)
```

## Arguments

- model:

  A fitted `KaguModel`.

- path:

  File path to write to (e.g. `"model.rds"`).

- include_data:

  Logical - whether to save the training data alongside the model.
  Defaults to `TRUE` so that `$effects()` works immediately after
  loading without needing to reattach the data.

## Value

Invisible `NULL`.
