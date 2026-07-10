# Discover causal structure from data

Fits a GCM to every candidate DAG over the supplied variables and
returns a posterior probability distribution over those DAGs, \$\$P(G
\mid X) \propto P(X \mid G)\\ P(G),\$\$ where each graph's marginal
likelihood \\P(X \mid G)\\ is obtained from the per-node
Gaussian-process fits and \\P(G)\\ is the prior.

## Usage

``` r
kagu_discover(
  data,
  nodes = NULL,
  dags = NULL,
  disallowed = NULL,
  required = NULL,
  mechanisms = NULL,
  prior = "uniform",
  allow_empty = FALSE,
  ...
)
```

## Arguments

- data:

  A `data.frame` with one column per variable.

- nodes:

  Optional character vector of node names. Defaults to all columns of
  `data`.

- dags:

  Optional explicit list of candidate DAGs to score (each a named list
  of parent vectors). If provided, `disallowed` and `required` are
  ignored, and the search space is restricted exactly to these DAGs.

- disallowed:

  Optional list of length-2 character vectors `c(from, to)`, each
  forbidding the directed edge `from -> to`. Useful for encoding a known
  temporal ordering and for pruning the search space.

- required:

  Optional list of length-2 character vectors `c(from, to)`, each
  requiring the directed edge `from -> to` to be present in all
  candidate DAGs.

- mechanisms:

  Optional named list of
  [Mechanism](https://causalabs.github.io/kagu-r/reference/Mechanism.md)
  instances, one per node. Any node not specified receives a
  [GPMechanism](https://causalabs.github.io/kagu-r/reference/GPMechanism.md).

- prior:

  Character - prior over DAGs. Only `"uniform"` is supported.

- allow_empty:

  Logical - whether to include the completely edgeless graph in the
  search space (default `FALSE`). Usually, researchers are looking for
  at least some causal structure, so the empty graph is omitted.

- ...:

  Additional arguments forwarded to each node's mechanism `$fit()`.

## Value

A
[DiscoveryResult](https://causalabs.github.io/kagu-r/reference/DiscoveryResult.md).

## Details

This is normally called as the static method
**`KaguModel$discover(data, ...)`**; `kagu_discover()` is the underlying
function.

Two properties keep the computation tractable. First, by the Markov
property a DAG's marginal likelihood factorises over nodes, \\\log P(X
\mid G) = \sum_i \log P(X_i \mid \mathrm{Pa}\_G(X_i))\\, so a graph's
score is a sum of independent per-node terms. Second, the same
`(node, parent-set)` local model recurs across many DAGs, so each unique
local model is fitted only once and cached - collapsing the work from
`#DAGs` fits to at most `p * 2^(p-1)`.

## Priors

Only a uniform prior over DAGs is currently supported
(`prior = "uniform"`), so the posterior is driven entirely by the
marginal likelihoods. The `prior` argument is reserved for future
options - e.g. sparsity-favouring priors, edge/temporal constraints, or
fully custom priors over structures.

## Marginal likelihoods and Markov equivalence

Each node's marginal likelihood is the closed-form type-II
(empirical-Bayes) evidence of its Gaussian-process model, so discovery
involves no bridge sampling and no Stan compilation. Note that
Markov-equivalent DAGs are statistically indistinguishable from
observational data and will therefore receive (near-)equal posterior
mass - a faithful representation of structural uncertainty, not a
defect.

## Examples

``` r
if (FALSE) { # \dontrun{
set.seed(1)
n <- 300
a <- rnorm(n)
b <- 0.8 * a + rnorm(n, sd = 0.5)
c <- 1.2 * b + rnorm(n, sd = 0.5)
df <- data.frame(a = a, b = b, c = c)

# Equivalent to KaguModel$discover(df)
res <- kagu_discover(df)
res$summary()
res$edge_probabilities()
res$plot(true_dag = list(a = c(), b = "a", c = "b"))
} # }
```
