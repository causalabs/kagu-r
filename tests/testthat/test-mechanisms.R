test_that("LinearMechanism builds a valid brms formula with parents", {
  mech    <- LinearMechanism$new()
  formula <- mech$build_formula("y", c("x1", "x2"))
  expect_s3_class(formula, "brmsformula")
  expect_match(deparse(formula$formula), "y ~ x1 \\+ x2")
})

test_that("LinearMechanism builds intercept-only formula for root node", {
  mech    <- LinearMechanism$new()
  formula <- mech$build_formula("x", c())
  expect_s3_class(formula, "brmsformula")
  expect_match(deparse(formula$formula), "x ~ 1")
})

test_that("LinearMechanism returns correct family", {
  mech <- LinearMechanism$new()
  expect_s3_class(mech$family(), "brmsfamily")
})

test_that("LinearMechanism$predict_mean returns (n_chains, n_draws) matrix", {
  skip_on_cran()
  data    <- make_chain_data()
  mech    <- LinearMechanism$new()
  formula <- mech$build_formula("m", c("x"))
  fit     <- do.call(brms::brm, c(
    list(formula = formula, data = data, family = mech$family(),
         center = FALSE, refresh = 0),
    SAMPLE_KWARGS
  ))

  n_chains <- 2L
  n_draws  <- 500L
  pvals    <- list(x = matrix(0, nrow = n_chains, ncol = n_draws))
  preds    <- mech$predict_mean("m", c("x"), pvals, fit)

  expect_true(is.matrix(preds))
  expect_equal(dim(preds), c(n_chains, n_draws))
})

test_that("LinearMechanism recovers known coefficients", {
  skip_on_cran()
  set.seed(1L)
  n  <- 300L
  x  <- rnorm(n)
  y  <- 2.0 + 1.5 * x + rnorm(n, sd = 0.5)
  df <- data.frame(x = x, y = y)

  mech <- LinearMechanism$new()
  fit  <- do.call(brms::brm, c(
    list(formula = mech$build_formula("y", c("x")),
         data = df, family = mech$family(),
         center = FALSE, refresh = 0),
    SAMPLE_KWARGS
  ))

  draws <- posterior::as_draws_df(fit)
  expect_lt(abs(mean(draws$b_Intercept) - 2.0), 0.4)
  expect_lt(abs(mean(draws$b_x)         - 1.5), 0.3)
})
