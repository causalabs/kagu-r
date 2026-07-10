# Load a model saved with [`kagu_save()`](https://causalabs.github.io/kagu-r/reference/kagu_save.md)

Load a model saved with
[`kagu_save()`](https://causalabs.github.io/kagu-r/reference/kagu_save.md)

## Usage

``` r
kagu_load(path)
```

## Arguments

- path:

  File path to read from.

## Value

A `KaguModel` with `dag`, `mechanisms`, `traces`, and (if saved) `data`
restored. The model is marked as fitted.
