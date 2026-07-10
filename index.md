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

The system is represented as structural equations
`X_i = f_i(Pa(X_i), ε_i)`. Each `f_i` is a *mechanism* fitted as an
independent Bayesian model, so the joint distribution factorises over
nodes: `P(X_1, …, X_n) = ∏ P(X_i | Pa(X_i))`.

Effects are computed via the *do-operator*: fix the treatment node at
`x` (severing its incoming edges), propagate forward through the DAG in
topological order, and compare `E[Y | do(X = x)]` against
`E[Y | do(X = x')]`. Because each node’s mechanism is a Gaussian
process, non-linearities and interactions between parents are picked up
automatically, without being specified in the model formula.

## A worked example

``` r

library(kagu)

model <- KaguModel$new(
  dag = list(
    age = c(), sex = c(),
    sociality    = "age",
    food_sharing = "sociality",
    condition    = c("age", "sex", "sociality", "food_sharing")
  )
)
model$fit(data)

# Extract causal effect distributions
model$effects("sociality", "condition")$summary()
```

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
