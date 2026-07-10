# Contributing to kagu

Thanks for your interest in contributing. This document explains how to
report problems and propose changes.

## Reporting bugs and requesting features

Please open an issue at <https://github.com/causalabs/kagu-r/issues>.

For a bug report, the most useful thing you can include is a small,
reproducible example (a [reprex](https://reprex.tidyverse.org)) together
with the output of
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html). A minimal
DAG, a snippet of simulated data, and the exact call that misbehaves are
usually enough.

For a feature request, describe the causal question or workflow you are
trying to support, not only the API you have in mind, so we can weigh
alternatives.

## Branching and pull requests

We follow a GitHub-flow model. `main` is the live branch: it is
protected, and changes only reach it through a reviewed pull request
that passes all checks.

1.  Create a branch off `main` with a short, descriptive name for the
    work (e.g. `fix-sweep-hdi`, `discovery-edge-priors`).
2.  Make your change, keeping commits focused and their messages
    descriptive.
3.  **Bump the package version** in `DESCRIPTION` (every pull request
    must raise it above the current `main`).
4.  Open a pull request against `main`. Small, self-contained PRs are
    much easier to review and merge than large ones; if you are planning
    something big, open an issue first so we can agree on the approach.

Continuous integration runs `R CMD check` (including the test suite) on
Linux, macOS and Windows, and verifies the version bump, for every pull
request. All checks must pass before a PR can be merged.

## Development setup

The package is developed with [devtools](https://devtools.r-lib.org):

``` r

# install.packages("devtools")
devtools::load_all()      # load the package from source
devtools::test()          # run the testthat suite
devtools::document()      # regenerate NAMESPACE and man/ from roxygen
devtools::check()         # a local R CMD check
```

Please:

- **Add tests** for any behaviour you add or fix, under
  `tests/testthat/`. Fitting-dependent tests should `skip_on_cran()`.
- **Document exported functions** with roxygen2 and run
  `devtools::document()` so `man/` and `NAMESPACE` stay in sync with the
  source.
- **Match the surrounding style.** Follow the conventions already in the
  file you are editing (naming, comment density, the R6 mechanism
  interface, and the `[n_chains, n_draws]` posterior draw convention).
  `AGENTS.md` documents the package’s architecture and design decisions.

## Vignettes and the website

Documentation lives in `vignettes/` (rendered by pkgdown) and the
landing page in `index.md`. Several vignettes fit Gaussian processes and
take a little while to knit. The website is built and deployed
automatically from `main`; you do not need to build or commit it
yourself.
