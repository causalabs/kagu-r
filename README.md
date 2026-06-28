# kagu

**kagu** is an R package for fitting Bayesian graphical causal models (GCMs).
Models are specified as directed acyclic graphs (DAGs), fitted with full
Bayesian inference via [brms](https://paul-buerkner.github.io/brms/), and
causal effects are extracted by propagating interventions forward through the
structural model.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("your-org/kagu-r")
```

---

## Example

```r
library(kagu)

# --- 1. Fit ---
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

# --- 2. Causal effects ---

# Unit effect at the mean — equivalent to a regression coefficient
effect <- model$effects("smoking", "health")
effect$summary()
#> # A tibble: 1 × 8
#>   source  target  from    to   mean     sd hdi_lower hdi_upper
#>   <chr>   <chr>  <dbl> <dbl>  <dbl>  <dbl>     <dbl>     <dbl>
#> 1 smoking health   8.2   9.2  -0.31   0.09     -0.45     -0.16

# Effect of a 1-SD increase centred at the mean
model$effects("smoking", "health", std_units = TRUE)$summary()

# Specific contrast: non-smoker (0) vs heavy smoker (20 cig/day)
model$effects("smoking", "health", values = c(0, 20))$summary()

# Conditional effect: effect of smoking for older adults only
model$effects("smoking", "health", values = c(0, 20),
              conditions = list(age = 65))$summary()

# --- 3. Dose-response curve ---
sweep <- model$effects("smoking", "health", sweep = TRUE)
sweep$plot()

# --- 4. Summaries and diagnostics ---
model$summary()             # coefficient table across all nodes
model$diagnostics("health") # r-hat and ESS

# --- 5. Plots ---
model$plot_dag()              # DAG visualisation
model$plot_posterior("health")

# --- 6. Save and load ---
model$save("health_model.rds")
loaded <- KaguModel$load("health_model.rds")
```

---

## Core concepts

### Structural causal models

Kagu represents a causal system as structural equations:

```
X_i = f_i(Pa(X_i), ε_i)
```

Each `f_i` is a *mechanism* fitted as an independent Bayesian model. The joint
distribution factorises over nodes:

```
P(X_1, ..., X_n) = ∏ P(X_i | Pa(X_i))
```

### Causal effect estimation

Effects are computed via the *do-operator*:

1. Fix the treatment node at `x` (severing its incoming edges).
2. Propagate forward through the DAG in topological order.
3. Compare `E[Y | do(X = x)]` to `E[Y | do(X = x')]`.

---

## Mechanisms

| Class | Model | Use case |
|---|---|---|
| `LinearMechanism` | `X_i ~ Normal(α + Σ βⱼ Paⱼ, σ)` | Continuous, unbounded |

---

## Roadmap

- [x] Core DAG-based Bayesian fitting (brms)
- [x] Do-calculus effect estimation with full posterior
- [x] Sweep plots with HDI ribbon
- [x] R-hat diagnostics on effect posteriors
- [x] Save/load with data
- [x] Documentation site (pkgdown)
- [x] CI/CD (GitHub Actions)
- [ ] **GLM mechanisms** — `LogNormalMechanism`, `GammaMechanism`,
  `PoissonMechanism`, `NegBinomialMechanism`, `BernoulliMechanism`,
  `BetaMechanism`, `OrderedMechanism`
- [ ] User-defined priors via mechanism configuration
- [ ] Model fit diagnostics per node (posterior predictive checks, LOO)
- [ ] Model comparison per node (WAIC / LOO)
