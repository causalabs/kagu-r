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

# ---- Noise-averaged vs plug-in propagation ----------------------------------

saturating_model <- local({
  set.seed(RNG_SEED)
  x <- rnorm(N)
  m <- x + rnorm(N)
  y <- 3 * tanh(2 * m) + rnorm(N, sd = 0.2)
  model <- KaguModel$new(dag = CHAIN_DAG)
  do.call(model$fit, c(list(data = data.frame(x = x, m = m, y = y)), SAMPLE_KWARGS))
  model
})

test_that("noise-averaged and plug-in effects agree on a linear chain", {
  skip_on_cran()
  avg  <- chain_model$model$effects("x", "y", values = c(0, 1))
  plug <- chain_model$model$effects("x", "y", values = c(0, 1), noise = FALSE)
  expect_equal(mean(avg$samples), mean(plug$samples), tolerance = 0.05)
})

test_that("both arms share exogenous units, so a null contrast is exactly zero", {
  skip_on_cran()
  effect <- chain_model$model$effects("x", "y", values = c(0.5, 0.5))
  expect_true(all(effect$samples == 0))
})

test_that("noise averaging recovers the interventional mean through a saturating mechanism", {
  skip_on_cran()
  # True E[y | do(x = 1)] - E[y | do(x = 0)] with m = x + N(0, 1)
  z     <- stats::qnorm(stats::ppoints(1e5))
  truth <- mean(3 * tanh(2 * (1 + z))) - mean(3 * tanh(2 * z))   # ~1.92
  avg   <- mean(saturating_model$effects("x", "y", values = c(0, 1))$samples)
  plug  <- mean(saturating_model$effects("x", "y", values = c(0, 1), noise = FALSE)$samples)
  expect_lt(abs(avg - truth), 0.35)
  expect_lt(abs(avg - truth), abs(plug - truth))
})

test_that("effects() is reproducible under set.seed and leaves the RNG stream usable", {
  skip_on_cran()
  set.seed(7); a <- chain_model$model$effects("x", "y", values = c(0, 1))$samples
  set.seed(7); b <- chain_model$model$effects("x", "y", values = c(0, 1))$samples
  expect_identical(a, b)
  expect_false(identical(a, chain_model$model$effects("x", "y", values = c(0, 1))$samples))
})

test_that("noise = TRUE without training data gives an informative error", {
  skip_on_cran()
  model      <- chain_model$model$clone()
  model$data <- NULL
  expect_error(model$effects("x", "y", values = c(0, 1)), "noise = FALSE")
})

test_that("noise averaging treats the observed rows as the population", {
  skip_on_cran()
  set.seed(RNG_SEED)
  a <- rnorm(80L); b <- rnorm(80L)
  y <- a * b + rnorm(80L, sd = 0.3)
  model <- KaguModel$new(dag = list(a = c(), b = c(), y = c("a", "b")))
  do.call(model$fit, c(list(data = data.frame(a = a, b = b, y = y)), SAMPLE_KWARGS))

  # Per draw, the effect is the plain average over every observed b
  fit    <- model$traces$y
  manual <- colMeans(kagu:::.gp_eval(fit, cbind(1, b), matched = FALSE) -
                     kagu:::.gp_eval(fit, cbind(0, b), matched = FALSE))
  effect <- model$effects("a", "y", values = c(0, 1))
  expect_equal(as.vector(effect$samples), manual)

  # With no noisy node on the path there is no Monte Carlo error to reduce
  expect_equal(model$effects("a", "y", values = c(0, 1), n_mc = 1000L)$samples,
               effect$samples)
})
