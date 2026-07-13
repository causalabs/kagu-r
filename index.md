# Kagu

Bayesian graphical causal models in R. Specify a causal DAG, fit each
node with a Gaussian process, and estimate interventional effects with
full posterior uncertainty.

`remotes::install_github(“causalabs/kagu-r”)`
![](data:image/svg+xml;base64,PHN2ZyBjbGFzcz0ia2FndS1oZXJvX19jb3B5LWljb24iIHZpZXdib3g9IjAgMCAxNiAxNiIgd2lkdGg9IjE2IiBoZWlnaHQ9IjE2IiBmaWxsPSJjdXJyZW50Q29sb3IiIGFyaWEtaGlkZGVuPSJ0cnVlIj48cGF0aCBmaWxsLXJ1bGU9ImV2ZW5vZGQiIGQ9Ik0wIDYuNzVDMCA1Ljc4NC43ODQgNSAxLjc1IDVoMS41YS43NS43NSAwIDAxMCAxLjVoLTEuNWEuMjUuMjUgMCAwMC0uMjUuMjV2Ny41YzAgLjEzOC4xMTIuMjUuMjUuMjVoNy41YS4yNS4yNSAwIDAwLjI1LS4yNXYtMS41YS43NS43NSAwIDAxMS41IDB2MS41QTEuNzUgMS43NSAwIDAxOS4yNSAxNmgtNy41QTEuNzUgMS43NSAwIDAxMCAxNC4yNXYtNy41eiIgLz48cGF0aCBmaWxsLXJ1bGU9ImV2ZW5vZGQiIGQ9Ik01IDEuNzVDNSAuNzg0IDUuNzg0IDAgNi43NSAwaDcuNUMxNS4yMTYgMCAxNiAuNzg0IDE2IDEuNzV2Ny41QTEuNzUgMS43NSAwIDAxMTQuMjUgMTFoLTcuNUExLjc1IDEuNzUgMCAwMTUgOS4yNXYtNy41em0xLjc1LS4yNWEuMjUuMjUgMCAwMC0uMjUuMjV2Ny41YzAgLjEzOC4xMTIuMjUuMjUuMjVoNy41YS4yNS4yNSAwIDAwLjI1LS4yNXYtNy41YS4yNS4yNSAwIDAwLS4yNS0uMjVoLTcuNXoiIC8+PC9zdmc+)![](data:image/svg+xml;base64,PHN2ZyBjbGFzcz0ia2FndS1oZXJvX19jb3B5LWNoZWNrIiB2aWV3Ym94PSIwIDAgMTYgMTYiIHdpZHRoPSIxNiIgaGVpZ2h0PSIxNiIgZmlsbD0iY3VycmVudENvbG9yIiBhcmlhLWhpZGRlbj0idHJ1ZSI+PHBhdGggZmlsbC1ydWxlPSJldmVub2RkIiBkPSJNMTMuNzggNC4yMmEuNzUuNzUgMCAwMTAgMS4wNmwtNy4yNSA3LjI1YS43NS43NSAwIDAxLTEuMDYgMEwyLjIyIDkuMjhhLjc1Ljc1IDAgMDExLjA2LTEuMDZMNiAxMC45NGw2LjcyLTYuNzJhLjc1Ljc1IDAgMDExLjA2IDB6IiAvPjwvc3ZnPg==)

[Get started
→](https://causalabs.github.io/kagu-r/articles/quickstart.md)
[Reference](https://causalabs.github.io/kagu-r/reference/index.md)

## Overview

Kagu fits Bayesian graphical causal models (GCMs). You describe a system
as a directed acyclic graph (DAG) of which variables cause which, and
Kagu fits each node as an independent Gaussian process of its parents.
Causal questions (total effects, conditional effects, dose-response
curves) are answered by propagating interventions through that fitted
graph, so different questions are queries against the same fitted model.

## How it works

The system is represented as structural equations \\X_i =
f_i\bigl(\mathrm{Pa}(X_i),\\ \varepsilon_i\bigr)\\, where
\\\mathrm{Pa}(X_i)\\ are the parents of node \\i\\. Each \\f_i\\ is a
*mechanism* fitted as an independent Bayesian model, so the joint
distribution factorises over nodes:

\\P(X_1, \ldots, X_n) = \prod\_{i=1}^{n} P\bigl(X_i \mid
\mathrm{Pa}(X_i)\bigr).\\

Effects are computed via the *do-operator*: fix the treatment node at
\\x\\ (severing its incoming edges), propagate forward through the DAG
in topological order, and compare \\\mathbb{E}\\\left\[Y \mid
\mathrm{do}(X = x)\right\]\\ against \\\mathbb{E}\\\left\[Y \mid
\mathrm{do}(X = x')\right\]\\. Because each node’s mechanism is a
Gaussian process, non-linearities and interactions between parents are
picked up automatically, without being specified in the model formula.

## A worked example

Say we study 50 animals and want to know whether being more **social**
improves body **condition**. The causal story has three wrinkles that
trip up a naive regression:

- **age** is a common cause of both sociality and condition (a
  *confounder*),
- sociality drives **food sharing**, which itself improves condition (a
  *mediator*),
- and **sex** changes how sociality pays off (an *effect modifier*).

We simulate data with exactly that structure, so we know the ground
truth: here, sociality *helps females and harms males*.

``` r

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

We hand Kagu the DAG and fit every node as a Gaussian process of its
parents:

``` r

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

Now ask the causal question. The **total effect** of sociality on
condition is a do-calculus query: Kagu adjusts for the confounder `age`
and propagates through the mediator `food_sharing` automatically,
returning a full posterior rather than a point estimate.

``` r

model$effects("sociality", "condition")$summary()
#> # A tibble: 1 x 8
#>   source    target       from    to  mean    sd hdi_lower hdi_upper
#>   <chr>     <chr>       <dbl> <dbl> <dbl> <dbl>     <dbl>     <dbl>
#> 1 sociality condition -0.0627 0.937 0.220 0.465    -0.460      1.05
```

On average the effect looks small and its interval brushes zero, but
that is *not* “no effect”: the positive effect in females and the
negative effect in males cancel out. The **same fitted model** answers
the follow-up without refitting, just by conditioning on `sex`:

``` r

model$effects("sociality", "condition", conditions = list(sex = -1))$summary()$mean  # females
#> [1] 1.32
model$effects("sociality", "condition", conditions = list(sex =  1))$summary()$mean  # males
#> [1] -0.95
```

Strongly positive for females (\\\approx 1.3\\), negative for males
(\\\approx -1.0\\). That interaction was never written into a formula —
each node is a Gaussian process, so it is learned during fitting and
recovered on query. The [ecology case
study](https://causalabs.github.io/kagu-r/articles/ecology_case_study.md)
works through this same example end to end, including causal discovery
and the Table II fallacy.

## Learn more

- [Quickstart](https://causalabs.github.io/kagu-r/articles/quickstart.md):
  a first model, fit and queried in a few lines.
- [Confounder, mediator, collider,
  M-bias](https://causalabs.github.io/kagu-r/articles/comparison.md):
  four classic DAG patterns compared against OLS.
- [Case study: sociality and fitness in
  ecology](https://causalabs.github.io/kagu-r/articles/ecology_case_study.md):
  the full version of the example above, including causal discovery and
  the Table II fallacy.
- [Modelling
  interactions](https://causalabs.github.io/kagu-r/articles/interactions.md):
  how the Gaussian process mechanism recovers interactions
  automatically.
- [Causal structure
  discovery](https://causalabs.github.io/kagu-r/articles/discovery.md):
  treating the DAG itself as uncertain and estimating a posterior over
  structures.
