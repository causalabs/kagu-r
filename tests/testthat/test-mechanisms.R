# Tests for GPMechanism (the default and only mechanism)

test_that("GPMechanism fits a node with parents and predicts mean draws", {
  set.seed(1L)
  n  <- 300L
  x  <- rnorm(n)
  y  <- 1.5 * x + rnorm(n, sd = 0.3)
  df <- data.frame(x = x, y = y)

  mech <- GPMechanism$new(num_results = 400L)
  fit  <- mech$fit("y", "x", df)

  expect_identical(fit$kind, "gp")
  expect_equal(mech$posterior_shape(fit)[["n_draws"]], 400L)

  pvals <- list(x = matrix(rnorm(400L), nrow = 1L))
  preds <- mech$predict_mean("y", "x", pvals, fit)
  expect_true(is.matrix(preds))
  expect_equal(dim(preds), c(1L, 400L))
})

test_that("GPMechanism models a root node as a marginal Normal", {
  set.seed(1L)
  df   <- data.frame(x = rnorm(500L, mean = 3, sd = 2))
  mech <- GPMechanism$new(num_results = 500L)
  fit  <- mech$fit("x", character(0), df)

  expect_identical(fit$kind, "root")
  preds <- mech$predict_mean("x", character(0), list(), fit)
  expect_equal(dim(preds), c(1L, 500L))
  # posterior of the mean is centred near the sample mean
  expect_lt(abs(mean(preds) - 3), 0.3)
})

test_that("GPMechanism recovers a known linear slope via node_terms", {
  set.seed(1L)
  n  <- 500L
  x  <- rnorm(n)
  y  <- 2.5 * x + rnorm(n, sd = 0.4)
  df <- data.frame(x = x, y = y)

  mech  <- GPMechanism$new()
  fit   <- mech$fit("y", "x", df)
  terms <- mech$node_terms("y", "x", df, fit)

  slope <- terms[terms$term == "x", ]
  expect_lt(abs(slope$mean - 2.5), 0.2)      # correct point estimate
  expect_lt(slope$sd, 0.15)                  # and tight uncertainty

  sigma <- terms[terms$term == "sigma (noise)", ]
  expect_lt(abs(sigma$mean - 0.4), 0.1)
})

test_that("GPMechanism marginal likelihood favours the true parent set", {
  set.seed(1L)
  n  <- 400L
  a  <- rnorm(n); b <- rnorm(n)
  cc <- 1.0 * a + 0.8 * b + rnorm(n, sd = 0.3)
  noise <- rnorm(n)
  df <- data.frame(a = a, b = b, c = cc, noise = noise)

  mech <- GPMechanism$new()
  ll_true  <- mech$log_marglik("c", c("a", "b"), df)
  ll_under <- mech$log_marglik("c", "a", df)
  ll_wrong <- mech$log_marglik("c", "noise", df)
  ll_over  <- mech$log_marglik("c", c("a", "b", "noise"), df)

  expect_gt(ll_true, ll_under)   # missing a real parent is worse
  expect_gt(ll_true, ll_wrong)   # an irrelevant parent is worse
  expect_gt(ll_true, ll_over)    # an extra noise parent is penalised (Occam)
})
