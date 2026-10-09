# Gaussian-process mechanism (default)

Models each node's conditional mean as an **exact Gaussian process** of
its parents, with a squared-exponential (ARD) kernel - one lengthscale
per parent, which captures non-linearity and interactions automatically.
Kernel hyperparameters (lengthscales, signal and noise variance) are
point-estimated by type-II maximum likelihood (`hyper = "ml"`) or
maximum a posteriori (`hyper = "map"`) under the priors in
[`gp_prior()`](https://causalabs.github.io/kagu-r/reference/gp_prior.md);
there is no MCMC. A node with no parents is modelled by its marginal (a
Normal). This is Kagu's default and only mechanism.

The default
[`gp_prior()`](https://causalabs.github.io/kagu-r/reference/gp_prior.md)
is log-uniform over the box the likelihood search is bounded to, so the
default ML fit *is* the MAP fit under that prior. Supplying any other
prior switches the default to `hyper = "map"`, so the prior you check
with `$prior_predictive()` is the prior the fit actually uses.

Structure discovery uses the **exact** GP log marginal likelihood,
which - in contrast to fast basis/eigenfunction approximations - is well
calibrated: for a genuinely unidentifiable (e.g. linear-Gaussian) edge
it does not manufacture spurious confidence about direction.

Posterior function samples for effect propagation are drawn by
**pathwise / decoupled sampling** (a random-feature prior plus the exact
data update), so each draw is a coherent function evaluable at any
point - giving correctly correlated uncertainty for do-calculus
contrasts.

## Super class

[`Mechanism`](https://causalabs.github.io/kagu-r/reference/Mechanism.md)
-\> `GPMechanism`

## Public fields

- `num_results`:

  Integer - number of posterior draws to keep (default 1000).

- `n_features`:

  Integer - random Fourier features for pathwise sampling (default 300).

- `jitter`:

  Numeric - diagonal jitter for numerical stability (default 1e-6).

- `prior`:

  A
  [`gp_prior()`](https://causalabs.github.io/kagu-r/reference/gp_prior.md) -
  priors on the (standardised) hyperparameters.

- `hyper`:

  Character - `"ml"` (type-II maximum likelihood) or `"map"`.

## Methods

### Public methods

- [`GPMechanism$new()`](#method-GPMechanism-initialize)

- [`GPMechanism$fit()`](#method-GPMechanism-fit)

- [`GPMechanism$predict_mean()`](#method-GPMechanism-predict_mean)

- [`GPMechanism$log_marglik()`](#method-GPMechanism-log_marglik)

- [`GPMechanism$posterior_shape()`](#method-GPMechanism-posterior_shape)

- [`GPMechanism$node_terms()`](#method-GPMechanism-node_terms)

- [`GPMechanism$simulate_prior()`](#method-GPMechanism-simulate_prior)

- [`GPMechanism$simulate_posterior()`](#method-GPMechanism-simulate_posterior)

- [`GPMechanism$function_draws()`](#method-GPMechanism-function_draws)

- [`GPMechanism$clone()`](#method-GPMechanism-clone)

------------------------------------------------------------------------

### `GPMechanism$new()`

Create a new GPMechanism.

#### Usage

    GPMechanism$new(
      num_results = 1000L,
      n_features = 300L,
      jitter = 1e-06,
      prior = gp_prior(),
      hyper = NULL
    )

#### Arguments

- `num_results`:

  Integer - number of posterior draws to keep.

- `n_features`:

  Integer - number of random Fourier features.

- `jitter`:

  Numeric - diagonal jitter added to the kernel.

- `prior`:

  A
  [`gp_prior()`](https://causalabs.github.io/kagu-r/reference/gp_prior.md)
  object (default: the implicit log-uniform box).

- `hyper`:

  `"ml"` or `"map"`. `NULL` (default) picks `"ml"` when every kernel
  prior is log-uniform (where ML and MAP coincide) and `"map"`
  otherwise.

------------------------------------------------------------------------

### `GPMechanism$fit()`

Fit the node (exact GP if it has parents, marginal Normal if not).

#### Usage

    GPMechanism$fit(node, parents, data, ...)

#### Arguments

- `node`:

  Character scalar - the node name (response).

- `parents`:

  Character vector of parent node names.

- `data`:

  A `data.frame` with columns for the node and its parents.

- `...`:

  Backend-specific arguments.

------------------------------------------------------------------------

### `GPMechanism$predict_mean()`

Conditional mean draws (see
[Mechanism](https://causalabs.github.io/kagu-r/reference/Mechanism.md)).

#### Usage

    GPMechanism$predict_mean(node, parents, parent_values, fit)

#### Arguments

- `node`:

  Character scalar - the node name.

- `parents`:

  Character vector of parent node names.

- `parent_values`:

  Named list of `[n_chains, n_draws]` matrices, one per parent, giving
  the parent values for each posterior draw.

- `fit`:

  A fit object from `$fit()`.

------------------------------------------------------------------------

### `GPMechanism$log_marglik()`

Exact GP log marginal likelihood (see
[Mechanism](https://causalabs.github.io/kagu-r/reference/Mechanism.md)).

#### Usage

    GPMechanism$log_marglik(node, parents, data)

#### Arguments

- `node`:

  Character scalar - the node name.

- `parents`:

  Character vector of parent node names.

- `data`:

  A `data.frame`.

------------------------------------------------------------------------

### `GPMechanism$posterior_shape()`

Posterior shape (single chain).

#### Usage

    GPMechanism$posterior_shape(fit)

#### Arguments

- `fit`:

  A fit object from `$fit()`.

------------------------------------------------------------------------

### `GPMechanism$node_terms()`

Per-node summary rows: each parent's direct local effect (function
gradient at the parent means) plus the residual noise sd.

#### Usage

    GPMechanism$node_terms(node, parents, data, fit)

#### Arguments

- `node`:

  Character scalar - the node name.

- `parents`:

  Character vector of parent node names.

- `data`:

  A `data.frame`.

- `fit`:

  A fit object from `$fit()`.

------------------------------------------------------------------------

### `GPMechanism$simulate_prior()`

Prior predictive draws (see
[Mechanism](https://causalabs.github.io/kagu-r/reference/Mechanism.md)).
Each draw samples hyperparameters from `prior`, an exact GP function at
the parents' rows, and observation noise - on the standardised scale -
then maps back to the data's location and scale.

#### Usage

    GPMechanism$simulate_prior(node, parents, parent_values, data, ndraws)

#### Arguments

- `node`:

  Character scalar - the node name.

- `parents`:

  Character vector of parent node names.

- `parent_values`:

  Named list of `[ndraws, n]` matrices of simulated parent values, one
  per parent (empty for a root node).

- `data`:

  A `data.frame` - sets the node's location/scale and `n`.

- `ndraws`:

  Integer - number of prior predictive draws.

------------------------------------------------------------------------

### `GPMechanism$simulate_posterior()`

Posterior predictive draws (see
[Mechanism](https://causalabs.github.io/kagu-r/reference/Mechanism.md)):
the matched pathwise function draw at each row's parents plus noise,
with the noise sd drawn as in `$node_terms()`.

#### Usage

    GPMechanism$simulate_posterior(node, parents, parent_values, fit, draws, n)

#### Arguments

- `node`:

  Character scalar - the node name.

- `parents`:

  Character vector of parent node names.

- `parent_values`:

  Named list, one entry per parent: either a length-`n` vector shared by
  every draw (e.g. the observed parents) or a `[length(draws), n]`
  matrix with one row per draw.

- `fit`:

  A fit object from `$fit()`.

- `draws`:

  Integer vector - which posterior draws to use.

- `n`:

  Integer - number of observations to simulate per draw.

------------------------------------------------------------------------

### `GPMechanism$function_draws()`

Prior or posterior function draws on a grid (see
[Mechanism](https://causalabs.github.io/kagu-r/reference/Mechanism.md)).

#### Usage

    GPMechanism$function_draws(node, parents, grid, data, fit = NULL, ndraws = 50L)

#### Arguments

- `node`:

  Character scalar - the node name.

- `parents`:

  Character vector of parent node names.

- `grid`:

  A `data.frame` with one column per parent.

- `data`:

  A `data.frame` - sets the node's location/scale.

- `fit`:

  A fit object from `$fit()`, or `NULL` for prior draws.

- `ndraws`:

  Integer - number of function draws.

------------------------------------------------------------------------

### `GPMechanism$clone()`

The objects of this class are cloneable with this method.

#### Usage

    GPMechanism$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.

## Examples

``` r
mech <- GPMechanism$new()

# Weakly informative priors, fitted by MAP
mech <- GPMechanism$new(prior = gp_prior(lengthscale = prior_lognormal(0, 0.5)))
mech$hyper
#> [1] "map"
```
