# --- Prior families and gp_prior() -------------------------------------------

test_that("prior samplers respect their supports", {
  set.seed(1)
  lu <- .prior_sample(prior_loguniform(0.5, 2), 1000)
  expect_true(all(lu >= 0.5 & lu <= 2))
  expect_true(all(.prior_sample(prior_halfnormal(1), 1000) >= 0))
  expect_true(all(.prior_sample(prior_invgamma(3, 3), 1000) > 0))
  expect_true(all(.prior_sample(prior_lognormal(0, 1), 1000) > 0))
})

test_that("prior log densities are -Inf outside the support", {
  expect_equal(.prior_logdens(prior_loguniform(0.5, 2), 3), -Inf)
  expect_equal(.prior_logdens(prior_halfnormal(1), -1), -Inf)
  expect_true(is.finite(.prior_logdens(prior_invgamma(2, 2), 1)))
})

test_that("gp_prior() validates its components", {
  expect_s3_class(gp_prior(), "kagu_gp_prior")
  expect_error(gp_prior(lengthscale = 1), "must be a prior")
  expect_error(gp_prior(noise_sd = prior_normal(0, 1)), "positive")
  expect_s3_class(gp_prior(root_mean = prior_normal(2, 1)), "kagu_gp_prior")
})

test_that("default gp_prior() reproduces the historical search box exactly", {
  p <- gp_prior()
  expect_identical(.prior_log_bounds(p$lengthscale), c(-3, 4))
  expect_identical(2 * .prior_log_bounds(p$signal_sd), c(-4, 4))
  expect_identical(2 * .prior_log_bounds(p$noise_sd), c(-6, 2))
})

test_that("GPMechanism picks ML for log-uniform priors and MAP otherwise", {
  expect_equal(GPMechanism$new()$hyper, "ml")
  mech <- GPMechanism$new(prior = gp_prior(lengthscale = prior_lognormal(0, 0.5)))
  expect_equal(mech$hyper, "map")
  expect_equal(GPMechanism$new(prior = mech$prior, hyper = "ml")$hyper, "ml")
  expect_error(GPMechanism$new(prior = list()), "gp_prior")
})

test_that("MAP under the default log-uniform prior coincides with ML", {
  df <- make_chain_data(n = 80)
  ml  <- .gp_core(df$y, df["m"], hyper = "ml")
  map <- .gp_core(df$y, df["m"], hyper = "map")
  expect_equal(log(c(map$ls, map$sf2, map$sn2)), log(c(ml$ls, ml$sf2, ml$sn2)),
               tolerance = 1e-6)
})

test_that("an informative prior changes the MAP fit but logml stays the likelihood", {
  df  <- make_chain_data(n = 80)
  tight <- gp_prior(lengthscale = prior_lognormal(log(0.2), 0.05))
  fit <- .gp_core(df$y, df["m"], prior = tight, hyper = "map")
  expect_lt(abs(log(fit$ls) - log(0.2)), 0.5)   # pulled towards the prior
  ml <- .gp_core(df$y, df["m"])
  expect_lte(fit$logml, ml$logml + 1e-8)        # ML maximises the likelihood
})

test_that("default_mechanism is copied to every unlisted node", {
  mech  <- GPMechanism$new(prior = gp_prior(noise_sd = prior_halfnormal(1)))
  model <- KaguModel$new(CHAIN_DAG, default_mechanism = mech)
  expect_equal(model$mechanisms$y$hyper, "map")
  expect_false(identical(model$mechanisms$x, model$mechanisms$y))  # separate copies
  expect_error(KaguModel$new(CHAIN_DAG, default_mechanism = list()), "Mechanism")
})

# --- Prior predictive ----------------------------------------------------------

test_that("prior_predictive() needs no fit and returns [ndraws, n] per node", {
  df    <- make_chain_data(n = 60)
  model <- KaguModel$new(CHAIN_DAG)
  pp    <- model$prior_predictive(df, ndraws = 25)
  expect_s3_class(pp, "PredictiveResult")
  expect_equal(pp$type, "prior")
  expect_setequal(pp$nodes(), names(CHAIN_DAG))
  for (nd in pp$nodes()) expect_equal(dim(pp$yrep(nd)), c(25L, 60L))
  expect_error(KaguModel$new(CHAIN_DAG)$prior_predictive(), "Supply `data`")
})

test_that("a tighter prior gives a tighter prior predictive", {
  set.seed(3)
  df <- make_chain_data(n = 60)
  vague <- KaguModel$new(CHAIN_DAG)$prior_predictive(df, ndraws = 200)
  tight <- KaguModel$new(CHAIN_DAG, default_mechanism = GPMechanism$new(prior = gp_prior(
    signal_sd = prior_halfnormal(0.5), noise_sd = prior_halfnormal(0.5)
  )))$prior_predictive(df, ndraws = 200)
  sd_v <- apply(vague$yrep("y"), 1, sd); sd_t <- apply(tight$yrep("y"), 1, sd)
  expect_lt(stats::median(sd_t), stats::median(sd_v))
})

test_that("prior function draws handle duplicated (discrete) parent values", {
  set.seed(4)
  df <- data.frame(g = rep(c(-1, 1), 30), y = rnorm(60))
  model <- KaguModel$new(list(g = c(), y = "g"))
  pp <- model$prior_predictive(df, ndraws = 10)   # simulated g is continuous
  expect_true(all(is.finite(pp$yrep("y"))))
  fd <- model$function_draws("y", "g", ndraws = 5, n_grid = 10, prior = TRUE, data = df)
  expect_equal(nrow(fd), 50L)
})

# --- Posterior predictive -------------------------------------------------------

test_that("posterior_predictive() returns conditional and joint replicates", {
  skip_on_cran()
  df    <- make_chain_data(n = 80)
  model <- KaguModel$new(CHAIN_DAG)
  suppressMessages(model$fit(df))

  expect_error(KaguModel$new(CHAIN_DAG)$posterior_predictive(), "fit")
  for (type in c("conditional", "joint")) {
    pp <- model$posterior_predictive(ndraws = 40, type = type)
    expect_equal(pp$type, type)
    for (nd in pp$nodes()) expect_equal(dim(pp$yrep(nd)), c(40L, 80L))
  }
  expect_error(model$posterior_predictive(ndraws = 5000), "exceeds")
})

test_that("conditional replicates reproduce the observed mean and sd", {
  skip_on_cran()
  df    <- make_chain_data(n = 150)
  model <- KaguModel$new(CHAIN_DAG)
  suppressMessages(model$fit(df))
  s <- model$posterior_predictive(ndraws = 100)$summary(stats = c("mean", "sd"))
  expect_true(all(s$p_value > 0.01 & s$p_value < 0.99))
})

test_that("the joint correlation check flags a missing edge", {
  skip_on_cran()
  df  <- make_chain_data(n = 150)          # truth: x -> m -> y
  bad <- KaguModel$new(list(x = c(), m = c(), y = "m"))
  suppressMessages(bad$fit(df))
  cc <- bad$posterior_predictive(ndraws = 100, type = "joint")$cor_check()
  xm <- cc[cc$var1 == "x" & cc$var2 == "m", ]
  expect_true(xm$outside)
  expect_equal(xm$p_value, 0)              # observed corr ~0.97, replicated ~0
})

test_that("cor_check() refuses conditional results; plot_cor() returns a ggplot", {
  skip_on_cran()
  df    <- make_chain_data(n = 60)
  model <- KaguModel$new(CHAIN_DAG)
  suppressMessages(model$fit(df))
  expect_error(model$posterior_predictive(20)$cor_check(), "joint")
  expect_s3_class(model$posterior_predictive(20, "joint")$plot_cor(), "ggplot")
})

test_that("summary() accepts named custom statistics and validates nodes", {
  skip_on_cran()
  df    <- make_chain_data(n = 60)
  model <- KaguModel$new(CHAIN_DAG)
  suppressMessages(model$fit(df))
  pp <- model$posterior_predictive(20)
  s  <- pp$summary(stats = list(iqr = stats::IQR), nodes = "y")
  expect_equal(s$stat, "iqr")
  expect_true(s$p_value >= 0 && s$p_value <= 1)
  expect_error(pp$yrep("nope"), "Unknown node")
  expect_error(pp$summary(stats = list(stats::IQR)), "named list")
})

test_that("posterior function draws are returned in long format", {
  skip_on_cran()
  df    <- make_chain_data(n = 60)
  model <- KaguModel$new(CHAIN_DAG)
  suppressMessages(model$fit(df))
  fd <- model$function_draws("m", ndraws = 7, n_grid = 12)
  expect_equal(names(fd), c("draw", "parent_value", "value"))
  expect_equal(nrow(fd), 7L * 12L)
  expect_error(model$function_draws("x"), "no parents")
  expect_error(model$function_draws("y", "x"), "not a parent")
})
