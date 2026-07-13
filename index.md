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

The system is represented as structural equations
$X_i = f_i\bigl(\mathrm{Pa}(X_i),\, \varepsilon_i\bigr)$, where $\mathrm{Pa}(X_i)$
are the parents of node $i$. Each $f_i$ is a *mechanism* fitted as an independent
Bayesian model, so the joint distribution factorises over nodes:

$$P(X_1, \ldots, X_n) = \prod_{i=1}^{n} P\bigl(X_i \mid \mathrm{Pa}(X_i)\bigr).$$

Effects are computed via the *do-operator*: fix the treatment node at $x$
(severing its incoming edges), propagate forward through the DAG in topological
order, and compare $\mathbb{E}\!\left[Y \mid \mathrm{do}(X = x)\right]$ against
$\mathbb{E}\!\left[Y \mid \mathrm{do}(X = x')\right]$. Because each node's
mechanism is a Gaussian process, non-linearities and interactions between parents
are picked up automatically, without being specified in the model formula.

## A worked example

Say we study 50 animals and want to know whether being more **social** improves
body **condition**. The causal story has three wrinkles that trip up a naive
regression:

- **age** is a common cause of both sociality and condition (a *confounder*),
- sociality drives **food sharing**, which itself improves condition (a *mediator*),
- and **sex** changes how sociality pays off (an *effect modifier*).

We simulate data with exactly that structure, so we know the ground truth: here,
sociality *helps females and harms males*.

```r
library(kagu)

set.seed(42)
n   <- 50
age <- rnorm(n)
sex <- sample(c(-1, 1), n, replace = TRUE)          # -1 = female, 1 = male
sociality    <- 0.6 * age + rnorm(n)
food_sharing <- 0.8 * sociality + rnorm(n)
condition    <- 0.5 * age + 0.4 * food_sharing -
                sex * sociality +                    # the sex * sociality interaction
                rnorm(n)
data <- data.frame(age, sex, sociality, food_sharing, condition)
```

We hand Kagu the DAG and fit every node as a Gaussian process of its parents:

```r
model <- KaguModel$new(
  dag = list(
    age = c(), sex = c(),
    sociality    = "age",
    food_sharing = "sociality",
    condition    = c("age", "sex", "sociality", "food_sharing")
  )
)
model$fit(data)
```

Now ask the causal question. The **total effect** of sociality on condition is a
do-calculus query: Kagu adjusts for the confounder `age` and propagates through
the mediator `food_sharing` automatically, returning a full posterior rather than
a point estimate.

```r
model$effects("sociality", "condition")$summary()
#> # A tibble: 1 x 8
#>   source    target       from    to  mean    sd hdi_lower hdi_upper
#>   <chr>     <chr>       <dbl> <dbl> <dbl> <dbl>     <dbl>     <dbl>
#> 1 sociality condition -0.0627 0.937 0.220 0.465    -0.460      1.05
```

On average the effect looks small and its interval brushes zero, but that is *not*
"no effect": the positive effect in females and the negative effect in males
cancel out. The **same fitted model** answers the follow-up without refitting,
just by conditioning on `sex`:

```r
model$effects("sociality", "condition", conditions = list(sex = -1))$summary()$mean  # females
#> [1] 1.32
model$effects("sociality", "condition", conditions = list(sex =  1))$summary()$mean  # males
#> [1] -0.95
```

Strongly positive for females ($\approx 1.3$), negative for males ($\approx -1.0$).
That interaction was never written into a formula — each node is a Gaussian
process, so it is learned during fitting and recovered on query. The
[ecology case study](articles/ecology_case_study.html) works through this same
example end to end, including causal discovery and the Table II fallacy.

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
