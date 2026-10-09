# Priors for a GPMechanism's hyperparameters

Bundles the priors a
[GPMechanism](https://causalabs.github.io/kagu-r/reference/GPMechanism.md)
places on its hyperparameters. All are on the **standardised** scale:
before fitting, each parent and the node itself are centred and scaled
to unit sd, so `signal_sd = 1` means "functions vary about as much as
the data", and a lengthscale of 1 means "the function changes over
roughly one sd of the parent".

The defaults are log-uniform over exactly the box the type-II maximum
likelihood search is bounded to. Maximising the likelihood inside that
box is the same as finding the posterior mode under this prior, so the
defaults make Kagu's implicit prior explicit without changing any fit.
They are deliberately vague; see
[`vignette("predictive_checks")`](https://causalabs.github.io/kagu-r/articles/predictive_checks.md)
for what they imply about the data and how to choose something more
informative.

Root nodes (no parents) are modelled by a Normal; `root_mean` and
`root_sd` are used for prior predictive simulation only - a root node's
posterior uses a non-informative reference prior, which these weak
priors barely differ from once there is any data.

## Usage

``` r
gp_prior(
  lengthscale = prior_loguniform(exp(-3), exp(4)),
  signal_sd = prior_loguniform(exp(-2), exp(2)),
  noise_sd = prior_loguniform(exp(-3), exp(1)),
  root_mean = prior_normal(0, 1),
  root_sd = prior_lognormal(0, 0.5)
)
```

## Arguments

- lengthscale:

  Prior on each parent's kernel lengthscale.

- signal_sd:

  Prior on the kernel signal sd (the function's scale).

- noise_sd:

  Prior on the observation noise sd.

- root_mean:

  Prior on a root node's (standardised) mean.

- root_sd:

  Prior on a root node's (standardised) sd.

## Value

A `kagu_gp_prior` object.

## See also

[prior_families](https://causalabs.github.io/kagu-r/reference/prior_families.md),
[GPMechanism](https://causalabs.github.io/kagu-r/reference/GPMechanism.md)

## Examples

``` r
gp_prior()   # the implicit default
#> <kagu_gp_prior> (standardised scale)
#>   lengthscale  loguniform(lower = 0.04979, upper = 54.6)
#>   signal_sd    loguniform(lower = 0.1353, upper = 7.389)
#>   noise_sd     loguniform(lower = 0.04979, upper = 2.718)
#>   root_mean    normal(mean = 0, sd = 1)
#>   root_sd      lognormal(meanlog = 0, sdlog = 0.5)

# Weakly informative: smooth functions, signal and noise on the data's scale
gp_prior(
  lengthscale = prior_lognormal(0, 0.5),
  signal_sd   = prior_halfnormal(1),
  noise_sd    = prior_halfnormal(1)
)
#> <kagu_gp_prior> (standardised scale)
#>   lengthscale  lognormal(meanlog = 0, sdlog = 0.5)
#>   signal_sd    halfnormal(scale = 1)
#>   noise_sd     halfnormal(scale = 1)
#>   root_mean    normal(mean = 0, sd = 1)
#>   root_sd      lognormal(meanlog = 0, sdlog = 0.5)
```
