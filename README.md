# kagu

**kagu** is an R package for fitting Bayesian graphical causal models (GCMs).
Models are specified as directed acyclic graphs (DAGs); each node's conditional
distribution is modelled as a Gaussian process of its parents, and causal
effects are extracted by propagating interventions forward through the
structural model, carrying the full posterior throughout.

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

# Local effect at the mean (for a linear relationship, a regression coefficient)
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
model$summary()             # per-node direct local effects + residual noise
model$diagnostics("health") # draw count and residual noise per node

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
| `GPMechanism` | `X_i ~ GaussianProcess(Pa(X_i))` | Continuous; linear or nonlinear |

Each mechanism is a Gaussian process over a smooth eigen-basis of its parents,
fitted in closed form (conjugate Normal posterior; noise and amplitude set by
type-II maximum likelihood). This gives full posterior uncertainty on effects
and a closed-form marginal likelihood for structure discovery — no MCMC.

---

## Roadmap

- [x] Core DAG-based causal modelling (Gaussian-process mechanisms)
- [x] Do-calculus effect estimation with full posterior
- [x] Sweep plots with HDI ribbon
- [x] Causal structure discovery — posterior over DAGs
- [x] Bayesian model averaging of effects over the DAG posterior
- [x] Save/load with data
- [x] Documentation site (pkgdown)
- [x] CI/CD (GitHub Actions)
- [ ] Non-Gaussian outcome families (counts, binary, bounded)
- [ ] User-defined priors over DAGs (sparsity, edge/temporal constraints)
