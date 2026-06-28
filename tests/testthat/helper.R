# Shared fixtures for kagu tests

RNG_SEED <- 42L
N        <- 200L

SAMPLE_KWARGS <- list(draws = 500L, tune = 500L, chains = 2L, silent = 2L)

make_chain_data <- function(seed = RNG_SEED, n = N) {
  set.seed(seed)
  x <- rnorm(n)
  m <- 2.0 * x + rnorm(n, sd = 0.5)
  y <- -3.0 * m + rnorm(n, sd = 0.5)
  data.frame(x = x, m = m, y = y)
}

make_confounded_data <- function(seed = RNG_SEED, n = N) {
  set.seed(seed)
  age     <- rnorm(n, 50, 10)
  smoking <- 0.3 * age + rnorm(n, sd = 2)
  health  <- -0.5 * smoking + 0.2 * age + rnorm(n, sd = 2)
  data.frame(age = age, smoking = smoking, health = health)
}

CHAIN_DAG <- list(x = c(), m = c("x"), y = c("m"))

CONFOUNDED_DAG <- list(
  age     = c(),
  smoking = c("age"),
  health  = c("smoking", "age")
)
