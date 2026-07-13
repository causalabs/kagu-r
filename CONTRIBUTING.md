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
    must raise it above the current `main`). See [Choosing the
    version](#choosing-the-version) below.
4.  Open a pull request against `main`. Small, self-contained PRs are
    much easier to review and merge than large ones; if you are planning
    something big, open an issue first so we can agree on the approach.

Continuous integration runs `R CMD check` (including the test suite) on
Linux, macOS and Windows, and verifies the version bump, for every pull
request. All checks must pass before a PR can be merged.

## Choosing the version

Kagu uses `MAJOR.MINOR.PATCH` versions. Choose the component according
to the user-visible scope of the change, not the number of lines
changed. In particular, a patch is a small, backward-compatible release;
it does not have to be a bug fix.

### Patch: `x.y.Z`

Increment the patch version when a change refines existing behaviour
without materially expanding the public API or requiring users to change
their code. Examples include:

- small usability or visual improvements, such as better DAG spacing,
  labels, or arrow placement;
- bug fixes and numerical-stability improvements;
- performance improvements that preserve results and semantics;
- clearer documentation, examples, warnings, or error messages; and
- tests, maintenance, and internal refactoring with no intended API
  change.

A small enhancement to an existing feature can therefore be a patch. The
test is whether users receive a better version of something Kagu already
does, rather than a meaningfully new capability.

### Minor: `x.Y.0`

Increment the minor version for a backward-compatible release that
meaningfully expands what users can do. Examples include:

- a new exported function, class, or substantial plotting method;
- a new outcome family, mechanism type, diagnostic, or modelling
  workflow;
- new arguments or return data that form a significant supported
  extension to the public API; or
- a deprecation that gives users advance notice of a future breaking
  change.

Small supporting fixes and documentation may be included in the same
minor release without separate version increments.

### Major: `X.0.0`

Increment the major version when users may need to change existing code
or reconsider existing results. Examples include:

- removing or renaming exported functions, classes, methods, or
  arguments;
- changing established defaults or causal-effect semantics in a way that
  can materially alter results;
- changing documented return types or posterior sample shapes; or
- otherwise breaking compatibility with code written for the previous
  major version.

When a pull request sits near a boundary, describe the user impact in
the PR. Maintainers may adjust the version during review, especially
when several pull requests target `main` concurrently.

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
