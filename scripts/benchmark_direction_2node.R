# =============================================================================
# Calibration benchmark: 2-node edge orientation
#
# For a 2-node linear-Gaussian system a -> b, the two orientations (a -> b and
# b -> a) are Markov-equivalent: they imply the identical observational
# distribution and cannot be told apart from observational data. A *calibrated*
# structure score must therefore pick the true direction only ~50% of the time
# and should rarely be confident. This checks that kagu's exact-GP marginal
# likelihood behaves that way, and — crucially — that its confidence does NOT
# grow with n (the failure mode of truncated-basis GP approximations).
#
# Run from the package root:  Rscript scripts/benchmark_direction_2node.R
# =============================================================================
suppressMessages(devtools::load_all(".", quiet = TRUE))

mech <- GPMechanism$new()

# P(a -> b) = softmax of the marginal-likelihood difference (root terms cancel
# on the standardised scale).
p_ab <- function(a, b) {
  df <- data.frame(a = a, b = b)
  d  <- mech$log_marglik("b", "a", df) - mech$log_marglik("a", "b", df)
  1 / (1 + exp(-d))
}

run <- function(n, beta, reps = 25L) {
  ps <- vapply(seq_len(reps), function(s) {
    set.seed(s + 1000L * n + 7L * round(beta * 10))
    a <- rnorm(n)
    b <- beta * a + rnorm(n, sd = 1)                  # TRUE direction is a -> b
    p_ab(a, b)
  }, numeric(1))
  data.frame(
    n = n, beta = beta,
    mean_P_ab   = round(mean(ps), 3),                       # want ~0.50
    pct_pick_ab = round(mean(ps > 0.5) * 100, 0),           # want ~50%
    pct_conf    = round(mean(ps > 0.9 | ps < 0.1) * 100, 0) # overconfident? want low
  )
}

cat("2-node a -> b (Markov-equivalent orientations).\n")
cat("Calibrated => mean_P_ab ~ 0.50, pct_pick_ab ~ 50, pct_conf low & flat in n.\n\n")

grid <- expand.grid(n = c(50L, 100L, 150L, 200L), beta = c(0.5, 1.0, 1.5))
t0   <- Sys.time()
out  <- do.call(rbind, Map(function(n, b) run(n, b), grid$n, grid$beta))
print(out[order(out$beta, out$n), ], row.names = FALSE)
cat(sprintf("\nelapsed: %.0fs\n", as.numeric(Sys.time() - t0, units = "secs")))
