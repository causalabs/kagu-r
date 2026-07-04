# Kagu (R) — Developer Guide

This is the R port of Kagu. The Python version lives in a separate repo.

---

## Project identity

- **Package name**: `kagu`
- **Import convention**: `library(kagu)` — all classes exported at top level
- **CRAN**: not yet submitted
- **Companion**: Python version at `../pangolin/` (or `kagu` on GitHub)

---

## Stack

- **R ≥ 4.2**, managed with `renv` or direct `remotes`
- **BayesGPfit** — supplies the smooth eigen-basis; each node is a Gaussian
  process fitted as a conjugate Bayesian linear model on that basis (no MCMC)
- **posterior** — chain/draw structure, r-hat, ESS
- **bayestestR** — HDI computation
- **ggplot2** — all plots
- **R6** — OO class system (mirrors Python's class design)
- **pkgdown** — documentation site
- **testthat 3** — testing

---

## Package layout

```
R/
  kagu-package.R   # package-level docs and namespace imports
  dag.R            # validate_dag, topological_sort, ancestors, descendants, node_depth
  utils.R          # .chain_draw_shape, .to_draws_array, .quietly (internal only)
  mechanisms.R     # Mechanism R6 ABC + GPMechanism (+ conjugate GP internals)
  inference.R      # fit_node() dispatcher + .fit_and_marglik()
  model.R          # KaguModel R6 class
  effects.R        # EffectResult R6 class + compute_effect() + .propagate()
  discover.R       # kagu_discover() + DiscoveryResult (posterior over DAGs)
  plots.R          # kagu_plot_dag(), kagu_plot_posterior(), kagu_plot_discovery()
  summary.R        # build_summary_table()
  io.R             # kagu_save(), kagu_load()
tests/testthat/
  helper.R         # shared fixtures (make_chain_data etc.)
  test-dag.R
  test-mechanisms.R
  test-effects.R
  test-discover.R
  test-model.R
vignettes/
  quickstart.Rmd   # whistle-stop tour
  comparison.Rmd   # confounder/mediator/collider/M-bias vs OLS
  discovery.Rmd    # structure discovery — posterior over DAGs
_pkgdown.yml       # pkgdown site config
.github/workflows/
  ci.yml           # R-CMD-check on push/PR
  docs.yml         # pkgdown build + GitHub Pages deploy
```

---

## Key classes and design

### `KaguModel` (`R/model.R`)

Main R6 class. Constructor:
```r
KaguModel$new(dag, mechanisms = NULL)
```
- `dag`: named list, each node -> character vector of parent names.
- `mechanisms`: optional named list of `Mechanism` instances; defaults to
  `GPMechanism` for unspecified nodes.

Key methods:
```r
model$fit(data, ...)                       # ... forwarded to mechanism$fit()
model$effects(source, target, values=NULL, std_units=FALSE, conditions=NULL,
              sweep=FALSE, sweep_n=50, sweep_range=NULL, hdi=0.90)
model$summary(hdi_prob=0.90)
model$diagnostics(node=NULL)
model$plot_dag(node_pos=NULL)
model$plot_posterior(node)
model$save(path, include_data=TRUE)
KaguModel$load(path)                        # static method on the generator
KaguModel$discover(data, ...)               # static method — structure discovery
```

### `Mechanism` ABC + `GPMechanism` (`R/mechanisms.R`)

Backend-agnostic interface (a mechanism owns its full lifecycle):
- `$fit(node, parents, data, ...)` → opaque fit object
- `$predict_mean(node, parents, parent_values, fit)` → `[1, n_draws]` matrix
- `$log_marglik(node, parents, data)` → scalar (for discovery)
- `$posterior_shape(fit)` → `c(n_chains=1, n_draws)`
- `$node_terms(node, parents, data, fit)` → summary rows

`GPMechanism` (default and only): each node is a Gaussian process over the
`BayesGPfit` eigen-basis, fitted as a **conjugate Bayesian linear regression**
(`.gp_conj_fit`) with the noise/amplitude hyperparameters set by type-II ML.
Root nodes (no parents) use a marginal Normal. Key internals:
- `.gp_basis()` — standardise inputs + build eigen-basis (`GP.std.grids`,
  `GP.eigen.funcs.fast`, `GP.eigen.value`).
- `.gp_conj_fit()` — exact Gaussian posterior over basis coefficients, pre-drawn
  coefficient samples (for matched-draw propagation), residual-sd draws, and the
  log marginal likelihood.
- `.gp_type2_evidence()` — closed-form type-II marginal likelihood (Woodbury /
  matrix-determinant lemma); used by discovery.

**Important:** BayesGPfit's own `GP.fast.Bayes.fit` sampler is **not** used for
fitting — its noise/uncertainty quantification is designed for dense imaging
grids and is wrong for scattered regression (sigma scales with n; 1D degenerate).
We use only its (exported) basis functions and do the Bayesian regression
ourselves.

### `EffectResult` (`R/effects.R`)

Fields: `source`, `target`, `from_value`, `to_value`, `samples`,
`hdi_prob`, `std_units`, `conditions`, `sweep_values`.

Sample shapes:
- Scalar: `[n_chains, n_draws]` matrix
- Sweep: `[n_sweep, n_chains, n_draws]` array

Methods: `$summary()`, `$plot()`, `$diagnostics()`, `$is_sweep()`.

### `compute_effect()` (`R/effects.R`)

Intervention modes (same as Python):
- `values = c(a, b)`: explicit contrast
- `std_units = TRUE`: 1-SD effect centred at mean
- Default: epsilon central-difference gradient — `eps = sd * 1e-5`,
  `scale = 1/(2*eps)`. Display shows `mean → mean+1`.

### `.propagate()` (`R/effects.R`)

Internal. Iterates topological order, fixes source and any `conditions` nodes,
calls `mechanism$predict_mean()` for all others, stops at target.

---

## Posterior draw convention

All posterior samples use `(n_chains, n_draws)` matrices throughout. The GP is a
single "chain" (`n_chains = 1`); get the shape via
`mechanism$posterior_shape(fit)`. `predict_mean()` draws the GP function at
matched draw indices — `f^(d)(parent^(d))` — from the pre-drawn coefficient
samples, preserving coherent per-draw propagation.

`EffectResult$diagnostics()` still uses `.to_draws_array()` (in `utils.R`), but
with a single chain r-hat is not meaningful — treat the effect `sd`/HDI as the
uncertainty summary.

---

## Documentation

- **pkgdown** with Bootstrap 5 — dark navbar (`#212529`), teal accent (`#26a69a`)
- IBM Plex Sans body font, IBM Plex Mono code font (both via Google Fonts in `_pkgdown.yml`)
- Vignettes in `vignettes/` — rendered by knitr, shown under "Articles" in pkgdown

**Deploy locally:**
```r
pkgdown::build_site()
```

---

## CI/CD

- **`ci.yml`**: `r-lib/actions/check-r-package` on ubuntu + macos, release R
- **`docs.yml`**: `pkgdown::build_site_github_pages()` + `JamesIves/github-pages-deploy-action`

---

## Testing

```r
devtools::test()
# or
testthat::test_package("kagu")
```

- Fitting-dependent tests use `skip_on_cran()`.
- Fitting is fast (closed-form, no MCMC), so the whole suite runs in ~1s.
- `helper.R` provides `make_chain_data()`, `CHAIN_DAG`, `CONFOUNDED_DAG`
  (`SAMPLE_KWARGS` is now an empty list — the GP mechanism has no sampling knobs).

---

## Design decisions (R-specific)

**GP via conjugate Bayesian linear regression, not the Gibbs sampler**
BayesGPfit's `GP.fast.Bayes.fit` is built for dense imaging grids; its noise and
uncertainty are wrong for scattered regression (sigma scales with n, 1D is
degenerate). We use only its exported basis functions and do a conjugate Normal
regression ourselves — exact posterior, correct uncertainty, closed-form
evidence, and it fits in milliseconds.

**Type-II marginal likelihood for discovery**
Each node's evidence is the closed-form empirical-Bayes marginal likelihood of
the basis-linear GP model (`.gp_type2_evidence`). No bridge sampling, no MCMC —
discovery over a few variables runs in a fraction of a second.

**`posterior` package for effect diagnostics only**
`.to_draws_array()` reconstructs a `draws_array` from the `[1, n_draws]` effect
samples for `posterior::summarise_draws()`. With a single chain r-hat is `NA`;
the `sd`/HDI is the uncertainty summary.

**R6 over S3/S4**
Mirrors the Python class design closely. Users coming from Python or the R
`keras`/`tidymodels` ecosystem find it natural. All classes exported directly
so users write `KaguModel$new(...)` after `library(kagu)`.

**No `igraph` / `tidygraph` dependency**
DAG utilities are ~80 lines of pure R — Kahn's algorithm, BFS ancestors/
descendants. Keeps the dependency footprint small.

---

## Roadmap

- [x] Core DAG-based causal modelling (Gaussian-process mechanisms)
- [x] Do-calculus effect estimation with full posterior
- [x] Sweep plots with HDI ribbon (ggplot2)
- [x] Save/load with data (saveRDS/readRDS)
- [x] Documentation site (pkgdown)
- [x] CI/CD (GitHub Actions)
- [x] Causal structure discovery — posterior over DAGs via closed-form GP
  marginal likelihoods (`KaguModel$discover` / `kagu_discover`)
- [x] Bayesian model averaging of effects over the DAG posterior
  (`DiscoveryResult$effects`)
- [ ] Non-Gaussian outcome families (counts, binary, bounded)
- [ ] User-defined priors over DAGs (sparsity, edge/temporal constraints)
- [ ] Model fit diagnostics (posterior predictive checks, LOO)
