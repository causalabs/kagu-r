# Prior distributions for Gaussian-process hyperparameters

Constructors for the prior families accepted by
[`gp_prior()`](https://causalabs.github.io/kagu-r/reference/gp_prior.md).
Each returns a small `kagu_prior` object (a family name plus its
parameters) that Kagu can both sample from (prior predictive simulation)
and evaluate (MAP fitting).

- `prior_loguniform(lower, upper)` - uniform on `log(x)` over
  `[lower, upper]`.

- `prior_lognormal(meanlog, sdlog)` - `log(x) ~ Normal(meanlog, sdlog)`.

- `prior_invgamma(shape, scale)` - inverse-gamma; a common lengthscale
  prior because it puts little mass near zero (very wiggly functions).

- `prior_halfnormal(scale)` - `|Normal(0, scale)|`.

- `prior_normal(mean, sd)` - unbounded; only valid for `root_mean`.

## Usage

``` r
prior_loguniform(lower, upper)

prior_lognormal(meanlog = 0, sdlog = 1)

prior_invgamma(shape, scale)

prior_halfnormal(scale = 1)

prior_normal(mean = 0, sd = 1)
```

## Arguments

- lower, upper:

  Positive bounds of the log-uniform prior.

- meanlog, sdlog:

  Mean and sd of `log(x)` for the log-normal prior.

- shape:

  Shape of the inverse-gamma prior.

- scale:

  Scale of the inverse-gamma or half-normal prior.

- mean, sd:

  Mean and sd of the normal prior.

## Value

A `kagu_prior` object.

## See also

[`gp_prior()`](https://causalabs.github.io/kagu-r/reference/gp_prior.md)

## Examples

``` r
prior_lognormal(0, 0.5)
#> <kagu_prior> lognormal(meanlog = 0, sdlog = 0.5) 
prior_invgamma(5, 5)
#> <kagu_prior> invgamma(shape = 5, scale = 5) 
```
