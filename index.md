<div class="kagu-hero">
<div class="kagu-hero__inner">
<h1 class="kagu-hero__title">Kagu</h1>
<p class="kagu-hero__tagline">Bayesian graphical causal models in R — specify a DAG, fit it with full Bayesian inference, and read off interventional effects.</p>
<div class="kagu-hero__install"><code>remotes::install_github("your-org/kagu-r")</code></div>
<div class="kagu-hero__actions">
<a class="kagu-btn kagu-btn--primary" href="articles/quickstart.html">Get started →</a>
<a class="kagu-btn kagu-btn--ghost" href="reference/index.html">Reference</a>
</div>
</div>
</div>

## Why Kagu

A causal system is specified as a directed acyclic graph (DAG). Each node is
fitted as an independent Bayesian model via
[brms](https://paul-buerkner.github.io/brms/), and causal effects are extracted
by propagating interventions forward through the structural model — carrying the
full posterior all the way through.

```r
library(kagu)

model <- KaguModel$new(
  dag = list(
    age      = c(),
    smoking  = c("age"),
    exercise = c("age"),
    fitness  = c("exercise", "age"),
    health   = c("smoking", "fitness", "age")
  )
)
model$fit(data)

# Causal effect via the do-operator, returned as a tidy tibble
model$effects("smoking", "health")$summary()
#> # A tibble: 1 × 8
#>   source  target  from    to   mean     sd hdi_lower hdi_upper
#>   <chr>   <chr>  <dbl> <dbl>  <dbl>  <dbl>     <dbl>     <dbl>
#> 1 smoking health   8.2   9.2  -0.31   0.09     -0.45     -0.16
```

## Core concepts

### Structural causal models

The system is represented as structural equations `X_i = f_i(Pa(X_i), ε_i)`.
Each `f_i` is a *mechanism* fitted as an independent Bayesian model, so the
joint distribution factorises over nodes:
`P(X_1, …, X_n) = ∏ P(X_i | Pa(X_i))`.

### Causal effect estimation

Effects are computed via the *do-operator*: fix the treatment node at `x`
(severing its incoming edges), propagate forward through the DAG in topological
order, and compare `E[Y | do(X = x)]` against `E[Y | do(X = x')]`.

## Mechanisms

| Class | Model | Use case |
|---|---|---|
| `LinearMechanism` | `X_i ~ Normal(α + Σ βⱼ Paⱼ, σ)` | Continuous, unbounded |

## Roadmap

- [x] Core DAG-based Bayesian fitting (brms)
- [x] Do-calculus effect estimation with full posterior
- [x] Sweep plots with HDI ribbon
- [x] R-hat diagnostics on effect posteriors
- [x] Save / load with data
- [x] Documentation site (pkgdown)
- [x] CI/CD (GitHub Actions)
- [ ] **GLM mechanisms** — LogNormal, Gamma, Poisson, NegBinom, Bernoulli, Beta, Ordered
- [ ] User-defined priors via mechanism configuration
- [ ] Model fit diagnostics per node (posterior predictive checks, LOO)
- [ ] Model comparison per node (WAIC / LOO)
