chain_model <- local({
  data  <- make_chain_data()
  model <- KaguModel$new(dag = CHAIN_DAG)
  do.call(model$fit, c(list(data = data), SAMPLE_KWARGS))
  list(model = model, data = data)
})

test_that("unit effect on chain X->M->Y has correct sign", {
  skip_on_cran()
  effect <- chain_model$model$effects("x", "y")
  # True effect is -6 per unit of x
  expect_lt(mean(effect$samples), -3)
})

test_that("explicit values contrast recovers known effect magnitude", {
  skip_on_cran()
  effect <- chain_model$model$effects("x", "y", values = c(0, 1))
  expect_lt(mean(effect$samples), -3)
})

test_that("sweep returns (n_sweep, n_chains, n_draws) array", {
  skip_on_cran()
  effect <- chain_model$model$effects("x", "y", sweep = TRUE, sweep_n = 30L)
  expect_true(effect$is_sweep())
  dims <- dim(effect$samples)
  expect_length(dims, 3L)
  expect_equal(dims[[1]], 30L)      # n_sweep
  expect_equal(dims[[2]], 1L)       # GP mechanism is a single chain
  expect_gt(dims[[3]], 0L)          # n_draws
})

test_that("sweep is monotonically decreasing (x->m->y with pos then neg coef)", {
  skip_on_cran()
  effect <- chain_model$model$effects(
    "x", "y", sweep = TRUE, sweep_n = 20L, sweep_range = c(-2, 2)
  )
  n_sweep <- dim(effect$samples)[[1]]
  flat    <- matrix(effect$samples, nrow = n_sweep)
  means   <- rowMeans(flat)
  expect_true(all(diff(means) < 0), "Sweep should be monotonically decreasing")
})

test_that("conditions on mediator blocks causal path", {
  skip_on_cran()
  effect <- chain_model$model$effects("x", "y", values = c(0, 1),
                                       conditions = list(m = 0))
  # With M fixed, X has near-zero effect on Y
  expect_lt(abs(mean(effect$samples)), 1.5)
})

test_that("zero effect returned when source is not ancestor of target", {
  skip_on_cran()
  effect <- chain_model$model$effects("y", "x")
  expect_true(all(effect$samples == 0))
})

test_that("scalar summary returns a tibble with the expected columns", {
  skip_on_cran()
  effect <- chain_model$model$effects("x", "y", values = c(0, 1))
  df     <- effect$summary()
  expect_true(tibble::is_tibble(df))
  expect_equal(nrow(df), 1L)
  expect_identical(
    names(df),
    c("source", "target", "from", "to", "mean", "sd", "hdi_lower", "hdi_upper")
  )
  expect_identical(df$source, "x")
  expect_identical(df$target, "y")
})

test_that("sweep summary returns one tibble row per grid point", {
  skip_on_cran()
  effect <- chain_model$model$effects("x", "y", sweep = TRUE, sweep_n = 10L)
  df     <- effect$summary()
  expect_true(tibble::is_tibble(df))
  expect_equal(nrow(df), 10L)
  expect_identical(
    names(df),
    c("source", "target", "x", "mean", "sd", "hdi_lower", "hdi_upper")
  )
})

test_that("diagnostics returns a data frame with rhat column", {
  skip_on_cran()
  effect <- chain_model$model$effects("x", "y")
  diag   <- effect$diagnostics()
  expect_true(is.data.frame(diag) || tibble::is_tibble(diag))
  expect_true("rhat" %in% names(diag))
})
