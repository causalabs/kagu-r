# Tests for causal structure discovery (KaguModel$discover / kagu_discover)

test_that("enumerate_dags counts match known values", {
  expect_equal(length(enumerate_dags(c("a", "b", "c"))), 25L)
  expect_equal(length(enumerate_dags(c("a", "b", "c", "d"))), 543L)
  expect_equal(length(enumerate_dags("a")), 1L)
})

test_that("enumerate_dags respects disallowed edges and stays acyclic", {
  # forbid everything pointing into 'a' -> 'a' is forced to be a root
  dags <- enumerate_dags(c("a", "b", "c"),
                         disallowed = list(c("b", "a"), c("c", "a")))
  expect_true(length(dags) < 25L)

  for (dag in dags) {
    expect_length(dag[["a"]], 0L)                 # 'a' never has parents
    expect_silent(topological_sort(dag))         # acyclic
  }
})

test_that("enumerate_dags returns each DAG exactly once", {
  dags <- enumerate_dags(c("a", "b", "c", "d"))
  keys <- vapply(dags, function(d) paste(kagu:::.dag_edges(d), collapse = "|"),
                 character(1))
  expect_false(as.logical(anyDuplicated(keys)))
})

test_that("enumerate_dags errors on an over-large search space", {
  expect_error(enumerate_dags(letters[1:6]), "too large")
})

test_that(".logsumexp is numerically correct and handles -Inf", {
  expect_equal(kagu:::.logsumexp(log(c(1, 2, 3))), log(6))
  expect_equal(kagu:::.logsumexp(c(-Inf, log(2), log(3))), log(5))
  expect_equal(kagu:::.logsumexp(rep(-Inf, 3)), -Inf)
})

test_that(".local_key is parent-order invariant", {
  expect_identical(kagu:::.local_key("y", c("a", "b")),
                   kagu:::.local_key("y", c("b", "a")))
})

test_that(".dag_edges produces canonical parent->child labels", {
  dag <- list(a = c(), b = "a", c = c("a", "b"))
  expect_identical(kagu:::.dag_edges(dag),
                   sort(c("a→b", "a→c", "b→c")))
})

# --- brms-dependent integration tests --------------------------------------

test_that("discovery yields a valid posterior over DAGs", {
  skip_on_cran()
  set.seed(42)
  n <- 250
  a <- rnorm(n)
  b <- 0.8 * a + rnorm(n, sd = 0.5)
  cc <- 1.2 * b + rnorm(n, sd = 0.5)        # true: a -> b -> c
  df <- data.frame(a = a, b = b, c = cc)

  res <- KaguModel$discover(df, draws = 500L, tune = 500L, chains = 2L)

  expect_s3_class(res, "DiscoveryResult")
  expect_equal(sum(res$prob), 1, tolerance = 1e-8)
  expect_true(all(is.finite(res$prob)))
  expect_equal(res$n_models, 25L)
  expect_equal(res$n_unique_fits, 12L)          # 3 nodes: 3 * 2^2

  ep <- res$edge_probabilities()
  expect_true(all(ep$prob >= 0 & ep$prob <= 1))

  # The a-b and b-c skeleton edges should carry far more mass than a-c
  # (a and c are conditionally independent given b).
  edge_prob <- function(from, to) {
    row <- ep[ep$from == from & ep$to == to, ]
    if (nrow(row)) row$prob else 0
  }
  ab <- edge_prob("a", "b") + edge_prob("b", "a")
  bc <- edge_prob("b", "c") + edge_prob("c", "b")
  ac <- edge_prob("a", "c") + edge_prob("c", "a")
  expect_gt(ab, ac)
  expect_gt(bc, ac)
})

test_that("model-averaged effects return an EffectResult like the single-DAG API", {
  skip_on_cran()
  set.seed(42)
  n <- 250
  a <- rnorm(n)
  b <- 0.8 * a + rnorm(n, sd = 0.5)
  cc <- 1.2 * b + rnorm(n, sd = 0.5)               # true total a -> c ~ 0.96
  df <- data.frame(a = a, b = b, c = cc)

  # Temporal ordering identifies the chain (a < b < c).
  res <- KaguModel$discover(
    df, disallowed = list(c("b", "a"), c("c", "a"), c("c", "b")),
    draws = 500L, tune = 500L, chains = 2L
  )

  eff <- res$effects("a", "c")
  expect_s3_class(eff, "EffectResult")

  # Same tibble contract as KaguModel$effects().
  df_eff <- eff$summary()
  expect_identical(
    names(df_eff),
    c("source", "target", "from", "to", "mean", "sd", "hdi_lower", "hdi_upper")
  )
  expect_identical(df_eff$source, "a")
  expect_identical(df_eff$target, "c")

  # Sweep path also yields an EffectResult of the right shape.
  sw <- res$effects("a", "c", sweep = TRUE, sweep_n = 8L)
  expect_s3_class(sw, "EffectResult")
  expect_true(sw$is_sweep())
  expect_equal(nrow(sw$summary()), 8L)
})
