<div class="kagu-hero">
<div class="kagu-hero__inner">
<h1 class="kagu-hero__title">Kagu</h1>
<p class="kagu-hero__tagline">Bayesian graphical causal models in R. Specify a causal DAG, fit each node with a Gaussian process, and estimate interventional effects with full posterior uncertainty.</p>
<div class="kagu-hero__install">
<code id="kagu-install-cmd" data-cmd='remotes::install_github("causalabs/kagu-r")'>remotes::install_github("causalabs/kagu-r")</code>
<button class="kagu-hero__copy" type="button" data-copy-target="kagu-install-cmd" aria-label="Copy install command">
  <svg class="kagu-hero__copy-icon" viewBox="0 0 16 16" width="16" height="16" fill="currentColor" aria-hidden="true">
    <path fill-rule="evenodd" d="M0 6.75C0 5.784.784 5 1.75 5h1.5a.75.75 0 010 1.5h-1.5a.25.25 0 00-.25.25v7.5c0 .138.112.25.25.25h7.5a.25.25 0 00.25-.25v-1.5a.75.75 0 011.5 0v1.5A1.75 1.75 0 019.25 16h-7.5A1.75 1.75 0 010 14.25v-7.5z"></path>
    <path fill-rule="evenodd" d="M5 1.75C5 .784 5.784 0 6.75 0h7.5C15.216 0 16 .784 16 1.75v7.5A1.75 1.75 0 0114.25 11h-7.5A1.75 1.75 0 015 9.25v-7.5zm1.75-.25a.25.25 0 00-.25.25v7.5c0 .138.112.25.25.25h7.5a.25.25 0 00.25-.25v-7.5a.25.25 0 00-.25-.25h-7.5z"></path>
  </svg>
  <svg class="kagu-hero__copy-check" viewBox="0 0 16 16" width="16" height="16" fill="currentColor" aria-hidden="true">
    <path fill-rule="evenodd" d="M13.78 4.22a.75.75 0 010 1.06l-7.25 7.25a.75.75 0 01-1.06 0L2.22 9.28a.75.75 0 011.06-1.06L6 10.94l6.72-6.72a.75.75 0 011.06 0z"></path>
  </svg>
</button>
</div>
<div class="kagu-hero__actions">
<a class="kagu-btn kagu-btn--primary" href="articles/quickstart.html">Get started →</a>
<a class="kagu-btn kagu-btn--ghost" href="reference/index.html">Reference</a>
</div>
</div>
</div>

## Overview

Kagu fits Bayesian graphical causal models (GCMs). You describe a system as a
directed acyclic graph (DAG) of which variables cause which, and Kagu fits
each node as an independent Gaussian process of its parents. Causal questions
(total effects, conditional effects, dose-response curves) are answered by
propagating interventions through that fitted graph, so different questions
are queries against the same fitted model.

## How it works

The system is represented as structural equations `X_i = f_i(Pa(X_i), ε_i)`.
Each `f_i` is a *mechanism* fitted as an independent Bayesian model, so the
joint distribution factorises over nodes:
`P(X_1, …, X_n) = ∏ P(X_i | Pa(X_i))`.

Effects are computed via the *do-operator*: fix the treatment node at `x`
(severing its incoming edges), propagate forward through the DAG in
topological order, and compare `E[Y | do(X = x)]` against `E[Y | do(X = x')]`.
Because each node's mechanism is a Gaussian process, non-linearities and
interactions between parents are picked up automatically, without being
specified in the model formula.

## A worked example

Consider a small ecological study of 50 animals. `age` confounds `sociality`
and `condition` (a body condition index). `sociality` also drives
`food_sharing`, which in turn improves `condition`, so `food_sharing` is a
mediator. And `sociality` has a direct effect on `condition` that runs in
opposite directions for each sex.

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
[ecology case study](articles/ecology_case_study.html) vignette works through
this example in full, including causal discovery and the population-level
effect.

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

## Mechanisms

| Class | Model | Use case |
|---|---|---|
| `GPMechanism` | `X_i ~ GaussianProcess(Pa(X_i))` | Continuous; linear or nonlinear |

Each mechanism is an exact Gaussian process with a squared exponential (ARD)
kernel, one lengthscale per parent, so non-linearity and interactions are
captured automatically. Hyperparameters are set by type II maximum
likelihood, with no MCMC required, and pathwise (decoupled) sampling provides
coherent posterior draws for do-calculus effect propagation.

## Learn more

- [Quickstart](articles/quickstart.html): a first model, fit and queried in a
  few lines.
- [Confounder, mediator, collider, M-bias](articles/comparison.html): four
  classic DAG patterns compared against OLS.
- [Case study: sociality and fitness in ecology](articles/ecology_case_study.html):
  the full version of the example above, including causal discovery and the
  Table II fallacy.
- [Modelling interactions](articles/interactions.html): how the Gaussian
  process mechanism recovers interactions automatically.
- [Causal structure discovery](articles/discovery.html): treating the DAG
  itself as uncertain and estimating a posterior over structures.

## Roadmap

- [x] Core DAG-based causal modelling (Gaussian-process mechanisms)
- [x] Do-calculus effect estimation with full posterior
- [x] Sweep plots with HDI ribbon
- [x] Causal structure discovery via a posterior over DAGs
- [x] Bayesian model averaging of effects over the DAG posterior
- [x] Save / load with data
- [x] Documentation site (pkgdown)
- [x] CI/CD (GitHub Actions)
- [ ] Non-Gaussian outcome families (counts, binary, bounded)
- [ ] User-defined priors over DAGs (sparsity, edge constraints)
