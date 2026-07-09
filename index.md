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

```r
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
