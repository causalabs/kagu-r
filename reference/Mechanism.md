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

## Methods

### Public methods

- [`Mechanism$fit()`](#method-Mechanism-fit)

- [`Mechanism$predict_mean()`](#method-Mechanism-predict_mean)

- [`Mechanism$log_marglik()`](#method-Mechanism-log_marglik)

- [`Mechanism$posterior_shape()`](#method-Mechanism-posterior_shape)

- [`Mechanism$node_terms()`](#method-Mechanism-node_terms)

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

### `Mechanism$clone()`

The objects of this class are cloneable with this method.

#### Usage

    Mechanism$clone(deep = FALSE)

#### Arguments

- `deep`:

  Whether to make a deep clone.
