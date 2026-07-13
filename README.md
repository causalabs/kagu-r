# kagu

**kagu** is an R package for fitting Bayesian graphical causal models (GCMs).
You describe a system as a directed acyclic graph (DAG) of which variables
cause which. Kagu fits each node as an independent Gaussian process of its
parents, and causal questions (total effects, conditional effects,
dose-response curves) are answered by propagating interventions through that
fitted graph, so different questions are queries against the same fitted
model.

Full documentation: <https://www.kagu-r.org>

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("causalabs/kagu-r")
```

---

## How it works

The system is represented as structural equations `X_i = f_i(Pa(X_i), ε_i)`.
Each `f_i` is a *mechanism* fitted as an independent Bayesian model, so the
joint distribution factorises over nodes:

```
P(X_1, ..., X_n) = ∏ P(X_i | Pa(X_i))
```

Effects are computed via the *do-operator*: fix the treatment node at `x`
(severing its incoming edges), propagate forward through the DAG in
topological order, and compare `E[Y | do(X = x)]` against `E[Y | do(X = x')]`.
Because each node's mechanism is a Gaussian process, non-linearities and
interactions between parents are picked up automatically, without being
specified in the model formula.

---

## A worked example

Consider a small ecological study of 50 animals. `age` confounds
`sociality` and `condition` (a body condition index). `sociality` also
drives `food_sharing`, which in turn improves `condition`, so
`food_sharing` is a mediator. And `sociality` has a direct effect on
`condition` that runs in opposite directions for each sex.

```r
library(kagu)

# sex is -1 for female, 1 for male
n   <- 50
sex <- sample(c(-1, 1), n, replace = TRUE)
age <- rnorm(n)
sociality    <- 0.5 * age + rnorm(n, sd = 1.5)
food_sharing <- 0.8 * sociality + rnorm(n, sd = 1.5)
condition    <- 0.5 * age + 0.5 * sex +
  0.8 * food_sharing - 1.0 * sex * sociality +
  rnorm(n, sd = 0.5)
data <- data.frame(
  age, sex, sociality, food_sharing, condition
)

model <- KaguModel$new(
  dag = list(
    age = c(), sex = c(),
    sociality    = "age",
    food_sharing = "sociality",
    condition    = c("age", "sex",
                      "sociality", "food_sharing")
  )
)
model$fit(data)
```

A single regression of `condition` on all the other variables gives a
statistically significant *negative* coefficient for sociality:

```r
fit <- lm(condition ~ sociality + food_sharing + age + sex,
          data = data)
summary(fit)$coefficients
#>              Estimate Std. Error t value Pr(>|t|)
#> (Intercept)    0.0822     0.2111  0.3895   0.6987
#> sociality     -0.3939     0.1786 -2.2052   0.0326
#> food_sharing   0.8268     0.1659  4.9829   0.0000
#> age            0.3257     0.2352  1.3847   0.1730
#> sex            0.3511     0.2109  1.6651   0.1029
```

That coefficient does not correspond to sociality's causal effect.
`food_sharing` is a mediator, so including it in the regression removes
sociality's indirect (positive) route to condition. And with no
`sex * sociality` interaction term, the fit averages two large, opposite
direct effects into a single small coefficient with the wrong sign. Querying
the fitted Kagu model separately for each sex keeps them apart:

```r
# Female
model$effects("sociality", "condition",
              conditions = list(sex = -1))$summary()
#> # A tibble: 1 × 8
#>   source    target      from    to  mean    sd hdi_lower hdi_upper
#>   <chr>     <chr>      <dbl> <dbl> <dbl> <dbl>     <dbl>     <dbl>
#> 1 sociality condition -0.091 0.909  1.44 0.150      1.18      1.74

# Male
model$effects("sociality", "condition",
              conditions = list(sex = 1))$summary()
#> # A tibble: 1 × 8
#>   source    target      from    to   mean    sd hdi_lower hdi_upper
#>   <chr>     <chr>      <dbl> <dbl>  <dbl> <dbl>     <dbl>     <dbl>
#> 1 sociality condition -0.091 0.909 -0.573 0.128    -0.839    -0.338
```

For females, a one-unit increase in sociality raises condition by about 1.4;
for males it lowers it by about 0.6. Averaged over both sexes these largely
cancel, which is what the single regression coefficient reflects. The
interaction term was not specified anywhere: because each node is a Gaussian
process, it is captured during fitting. The
[ecology case study](https://causalabs.github.io/kagu-r/articles/ecology_case_study.html)
vignette works through this example in full, including causal discovery and
the population-level effect.

Other things you can do with the same fitted model:

```r
# Effect of a 1-SD increase in sociality, centred at the mean
model$effects("sociality", "condition", std_units = TRUE)$summary()

# Dose-response curve
sweep <- model$effects("sociality", "condition", sweep = TRUE)
sweep$plot()

# Summaries, diagnostics, and plots
model$summary()
model$plot_dag()

# Save and reload a fitted model
model$save("model.rds")
loaded <- KaguModel$load("model.rds")
```

---

## Relationship to regression

Regression tools such as `lm()`, `glm()`, and `brms` fit one conditional
distribution: an outcome given a chosen set of predictors. Reading a causal
effect off that fit requires you to choose the adjustment set (which variables
to condition on) and to specify any interaction terms, based on your
assumptions about the causal structure. `brms` additionally returns a full
posterior rather than a point estimate and standard error.

Kagu takes the causal graph as its input rather than a formula. Each node is
modelled separately, and an effect is computed by intervening on the fitted
graph, so the adjustment set for a query follows from the graph. A different
causal question, including one that integrates over a mediator, is another
query against the same fitted model.

Given a correctly specified regression, with the right adjustment set and
interaction terms, the two agree. The difference is whether the causal
assumptions are expressed as a formula or as a graph.

| | `lm()` / `glm()` | `brms` | Kagu |
|---|---|---|---|
| Input | a regression formula | a formula and priors | a causal DAG |
| Adjustment set | chosen when writing the formula | chosen when writing the formula | follows from the graph |
| Uncertainty | point estimate and standard error | full posterior | full posterior |
| Another causal question | new formula, refit | new formula, refit | new query, same fit |

---

## Mechanisms

| Class | Model | Use case |
|---|---|---|
| `GPMechanism` | `X_i ~ GaussianProcess(Pa(X_i))` | Continuous; linear or nonlinear |

Each mechanism is an exact Gaussian process with a squared exponential (ARD)
kernel, one lengthscale per parent, so non-linearity and interactions are
captured automatically. Hyperparameters are set by type II maximum
likelihood, with no MCMC required, and pathwise (decoupled) sampling provides
coherent posterior draws for do-calculus effect propagation.

---

## Learn more

- [Quickstart](https://causalabs.github.io/kagu-r/articles/quickstart.html):
  a first model, fit and queried in a few lines.
- [Confounder, mediator, collider, M-bias](https://causalabs.github.io/kagu-r/articles/comparison.html):
  four classic DAG patterns compared against OLS.
- [Case study: sociality and fitness in ecology](https://causalabs.github.io/kagu-r/articles/ecology_case_study.html):
  the full version of the example above, including causal discovery and the
  Table II fallacy.
- [Modelling interactions](https://causalabs.github.io/kagu-r/articles/interactions.html):
  how the Gaussian process mechanism recovers interactions automatically.
- [Causal structure discovery](https://causalabs.github.io/kagu-r/articles/discovery.html):
  treating the DAG itself as uncertain and estimating a posterior over
  structures.

---

## Roadmap

- [x] Core DAG-based causal modelling (Gaussian-process mechanisms)
- [x] Do-calculus effect estimation with full posterior
- [x] Sweep plots with HDI ribbon
- [x] Causal structure discovery via a posterior over DAGs
- [x] Bayesian model averaging of effects over the DAG posterior
- [x] Save/load with data
- [x] Documentation site (pkgdown)
- [x] CI/CD (GitHub Actions)
- [ ] Non-Gaussian outcome families (counts, binary, bounded)
- [ ] User-defined priors over DAGs (sparsity, edge/temporal constraints)
