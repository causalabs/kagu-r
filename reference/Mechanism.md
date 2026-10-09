# Mechanism abstract base class

A mechanism owns the full lifecycle of a node's conditional model,
decoupled from any particular inference backend. Subclasses implement:

- `$fit(node, parents, data, ...)` - fit the local model, returning an
  opaque fit object.

- `$predict_mean(node, parents, parent_values, fit)` - the conditional
  mean `E[node | parents]` for each posterior draw, as a
  `[n_chains, n_draws]` matrix (Gaussian processes use a single chain,
  so `n_chains = 1`).

- `$log_marglik(node, parents, data)` - the log marginal likelihood
  `log P(node | parents)`, used by structure discovery.

- `$posterior_shape(fit)` - `c(n_chains, n_draws)` for the fit.

- `$node_terms(node, parents, data, fit)` - summary rows for
  `model$summary()`.

- `$simulate_prior(node, parents, parent_values, data, ndraws)` -
  replicated node values drawn from the prior predictive distribution.

- `$simulate_posterior(node, parents, parent_values, fit, draws, n)` -
  replicated node values drawn from the posterior predictive
  distribution.

- `$function_draws(node, parents, grid, data, fit, ndraws)` - prior or
  posterior draws of the node's conditional-mean function on a grid.

## Methods

### Public methods

- [`Mechanism$fit()`](#method-Mechanism-fit)

- [`Mechanism$predict_mean()`](#method-Mechanism-predict_mean)

- [`Mechanism$log_marglik()`](#method-Mechanism-log_marglik)

- [`Mechanism$posterior_shape()`](#method-Mechanism-posterior_shape)

- [`Mechanism$node_terms()`](#method-Mechanism-node_terms)

- [`Mechanism$simulate_prior()`](#method-Mechanism-simulate_prior)

- [`Mechanism$simulate_posterior()`](#method-Mechanism-simulate_posterior)

- [`Mechanism$function_draws()`](#method-Mechanism-function_draws)

- [`Mechanism$clone()`](#method-Mechanism-clone)

------------------------------------------------------------------------

### `Mechanism$fit()`

Fit the node's local model.

#### Usage

    Mechanism$fit(node, parents, data, ...)

#### Arguments

- `node`:

  Character scalar - the node name (response).

- `parents`:

  Character vector of parent node names.

- `data`:

  A `data.frame` with columns for the node and its parents.

- `...`:

  Backend-specific arguments.

#### Returns

An opaque fit object.

------------------------------------------------------------------------

### `Mechanism$predict_mean()`

Conditional mean for each posterior draw.

#### Usage

    Mechanism$predict_mean(node, parents, parent_values, fit)

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

#### Returns

A `[n_chains, n_draws]` numeric matrix.

------------------------------------------------------------------------

### `Mechanism$log_marglik()`

Log marginal likelihood `log P(node | parents)`.

#### Usage

    Mechanism$log_marglik(node, parents, data)

#### Arguments

- `node`:

  Character scalar - the node name.

- `parents`:

  Character vector of parent node names.

- `data`:

  A `data.frame`.

#### Returns

A single numeric.

------------------------------------------------------------------------

### `Mechanism$posterior_shape()`

Posterior shape of a fit.

#### Usage

    Mechanism$posterior_shape(fit)

#### Arguments

- `fit`:

  A fit object from `$fit()`.

#### Returns

Named integer vector `c(n_chains, n_draws)`.

------------------------------------------------------------------------

### `Mechanism$node_terms()`

Summary rows for this node (used by `model$summary()`).

#### Usage

    Mechanism$node_terms(node, parents, data, fit)

#### Arguments

- `node`:

  Character scalar - the node name.

- `parents`:

  Character vector of parent node names.

- `data`:

  A `data.frame`.

- `fit`:

  A fit object from `$fit()`.

#### Returns

A `tibble` with columns `node`, `term`, `mean`, `sd`, `hdi_lower`,
`hdi_upper`.

------------------------------------------------------------------------

### `Mechanism$simulate_prior()`

Simulate node values from the prior predictive distribution.

#### Usage

    Mechanism$simulate_prior(node, parents, parent_values, data, ndraws)

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

#### Returns

A `[ndraws, n]` numeric matrix.

------------------------------------------------------------------------

### `Mechanism$simulate_posterior()`

Simulate node values from the posterior predictive distribution
(conditional mean plus observation noise).

#### Usage

    Mechanism$simulate_posterior(node, parents, parent_values, fit, draws, n)

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

#### Returns

A `[length(draws), n]` numeric matrix.

------------------------------------------------------------------------

### `Mechanism$function_draws()`

Prior (`fit = NULL`) or posterior draws of the node's conditional-mean
function, evaluated on a grid of parent values.

#### Usage

    Mechanism$function_draws(node, parents, grid, data, fit = NULL, ndraws = 50L)

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

#### Returns

A `[ndraws, nrow(grid)]` numeric matrix.

------------------------------------------------------------------------

### `Mechanism$clone()`

The objects of this class are cloneable with this method.

#### Usage

    Mechanism$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
