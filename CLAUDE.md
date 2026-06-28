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
- **brms** — Bayesian inference backend (one `brm()` call per node)
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
  mechanisms.R     # Mechanism R6 ABC + LinearMechanism
  inference.R      # fit_node() — wraps brms::brm()
  model.R          # KaguModel R6 class
  effects.R        # EffectResult R6 class + compute_effect() + .propagate()
  plots.R          # kagu_plot_dag(), kagu_plot_posterior()
  summary.R        # build_summary_table()
  io.R             # kagu_save(), kagu_load()
tests/testthat/
  helper.R         # shared fixtures (make_chain_data etc.)
  test-dag.R
  test-mechanisms.R
  test-effects.R
  test-model.R
vignettes/
  quickstart.Rmd   # whistle-stop tour
  comparison.Rmd   # confounder/mediator/collider/M-bias vs OLS
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
  `LinearMechanism` for unspecified nodes.

Key methods:
```r
model$fit(data, draws=1000, tune=1000, chains=4, ...)
model$effects(source, target, values=NULL, std_units=FALSE, conditions=NULL,
              sweep=FALSE, sweep_n=50, sweep_range=NULL, hdi=0.90)
model$summary(hdi_prob=0.90)
model$diagnostics(node=NULL)
model$plot_dag(node_pos=NULL)
model$plot_posterior(node)
model$save(path, include_data=TRUE)
KaguModel$load(path)   # static method on the generator
```

### `Mechanism` ABC + `LinearMechanism` (`R/mechanisms.R`)

Two methods to implement in subclasses:
- `$build_formula(node, parents)` → `brmsformula`
- `$predict_mean(node, parents, parent_values, fit)` → `[n_chains, n_draws]` matrix
- `$family()` → `brmsfamily`

`LinearMechanism` priors:
- `b_Intercept ~ Normal(0, prior_alpha)` (default 10)
- `b_<parent> ~ Normal(0, prior_beta)` (default 2)
- `sigma ~ HalfNormal(prior_sigma)` (default 1)

**Important:** brms always back-transforms and reports the uncentered
`b_Intercept` in the posterior, so `predict_mean` arithmetic is
straightforward regardless of internal sampling parameterisation.

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

All posterior samples use `(n_chains, n_draws)` matrices throughout:
- `n_chains = dim(as_draws_array(fit))[[2]]`
- `n_draws  = dim(as_draws_array(fit))[[1]]`

`t(draws[,, "b_Intercept"])` converts from `[n_draws, n_chains]` to
`[n_chains, n_draws]`.

`EffectResult$diagnostics()` uses `.to_draws_array()` (in `utils.R`) to
reconstruct a `posterior::draws_array` from `[n_chains, n_draws]` and passes
it to `posterior::summarise_draws()` for genuine r-hat.

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

- Brms-dependent tests use `skip_on_cran()`.
- Use `SAMPLE_KWARGS = list(draws=500, tune=500, chains=2, silent=2)` for speed.
- `helper.R` provides `make_chain_data()`, `CHAIN_DAG`, `CONFOUNDED_DAG`.

---

## Design decisions (R-specific)

**cmdstanr default backend**
`fit_node()` defaults to `backend = "cmdstanr"`. Users can override per-call
or via `model$fit(..., backend = "rstan")`. cmdstanr is in `Suggests` (not
`Imports`) so the package installs without it — brms will error at fit time
with a clear message if cmdstanr isn't available.

**`b_Intercept` is always the raw intercept**
brms may centre predictors internally for sampling efficiency but always
back-transforms to the uncentered `b_Intercept` in the posterior, so
`predict_mean = b_Intercept + Σ b_j * parent_j` is exact.

**`posterior` package for chain/draw structure**
Use `posterior::as_draws_array(fit)` (dims `[n_iter, n_chains, n_vars]`) then
transpose individual parameter slices to `[n_chains, n_draws]` matrices. This
preserves chain structure for genuine r-hat via `posterior::summarise_draws()`.

**R6 over S3/S4**
Mirrors the Python class design closely. Users coming from Python or the R
`keras`/`tidymodels` ecosystem find it natural. All classes exported directly
so users write `KaguModel$new(...)` after `library(kagu)`.

**No `igraph` / `tidygraph` dependency**
DAG utilities are ~80 lines of pure R — Kahn's algorithm, BFS ancestors/
descendants. Keeps the dependency footprint small.

---

## Roadmap

- [x] Core DAG-based Bayesian fitting (brms)
- [x] Do-calculus effect estimation with full posterior
- [x] Sweep plots with HDI ribbon (ggplot2)
- [x] R-hat diagnostics on effect posteriors (posterior package)
- [x] Save/load with data (saveRDS/readRDS)
- [x] Documentation site (pkgdown)
- [x] CI/CD (GitHub Actions)
- [ ] GLM mechanisms (LogNormal, Gamma, Poisson, NegBinom, Bernoulli, Beta, Ordered)
- [ ] User-defined priors via mechanism configuration
- [ ] Model fit diagnostics (posterior predictive checks, LOO)
- [ ] Model comparison per node (WAIC / LOO)
