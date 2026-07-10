# Gaussian-process mechanism (default)

Models each node's conditional mean as an **exact Gaussian process** of
its parents, with a squared-exponential (ARD) kernel - one lengthscale
per parent, which captures non-linearity and interactions automatically.
Kernel hyperparameters (lengthscales, signal and noise variance) are set
by type-II maximum likelihood; there is no MCMC. A node with no parents
is modelled by its marginal (a Normal). This is Kagu's default and only
mechanism.

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

## Methods

### Public methods

- [`GPMechanism$new()`](#method-GPMechanism-initialize)

- [`GPMechanism$fit()`](#method-GPMechanism-fit)

- [`GPMechanism$predict_mean()`](#method-GPMechanism-predict_mean)

- [`GPMechanism$log_marglik()`](#method-GPMechanism-log_marglik)

- [`GPMechanism$posterior_shape()`](#method-GPMechanism-posterior_shape)

- [`GPMechanism$node_terms()`](#method-GPMechanism-node_terms)

- [`GPMechanism$clone()`](#method-GPMechanism-clone)

------------------------------------------------------------------------

### `GPMechanism$new()`

Create a new GPMechanism.

#### Usage

    GPMechanism$new(num_results = 1000L, n_features = 300L, jitter = 1e-06)

#### Arguments

- `num_results`:

  Integer - number of posterior draws to keep.

- `n_features`:

  Integer - number of random Fourier features.

- `jitter`:

  Numeric - diagonal jitter added to the kernel.

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
```
