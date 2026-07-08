suppressMessages(devtools::load_all(".", quiet = TRUE))
library(dplyr)
library(ggplot2)

# Ensure reproducibility
# set.seed(1)

# --- Configuration ---
n_iters <- 50
n_samples_range <- 50:150   # exact GP is O(n^3): draw a modest sample per scenario
nodes_range <- 5:20
# Competitors are fully random DAGs (not near-miss mutations), so every candidate
# has a different parent set per node and the unique-fit cache barely helps —
# hence a modest candidate pool.
n_candidate_dags <- 20
edge_prob <- 0.15 # Keep sparsity reasonable for 20 nodes
noise_sd <- 0.5

cli::cli_alert_info("Starting discovery benchmark: {n_iters} random scenarios")
cli::cli_bullets(c(
  "*" = "Samples per scenario: {min(n_samples_range)}-{max(n_samples_range)} (random)",
  "*" = "Nodes: {min(nodes_range)} to {max(nodes_range)}",
  "*" = "Candidate DAGs per scenario: {n_candidate_dags} (including the true one)"
))

# --- Helpers ---

# Non-linear transformations for generating data
funcs <- list(
  function(x) sin(2 * x),
  function(x) cos(2 * x),
  function(x) tanh(x) * 2,
  function(x) sign(x) * sqrt(abs(x)) * 1.5,
  function(x) x^2 / 2 - 1,
  function(x) 0.8 * x
)

generate_dag_and_data <- function(p, n, edge_prob, noise_sd) {
  # Generate a random topological sort
  nodes <- paste0("x", 1:p)
  dag <- setNames(vector("list", p), nodes)

  # Populate edges based on topological sort to ensure acyclicity
  for (j in 2:p) {
    for (i in 1:(j-1)) {
      if (runif(1) < edge_prob) {
        dag[[nodes[j]]] <- c(dag[[nodes[j]]], nodes[i])
      }
    }
  }

  # Generate data
  data <- as.data.frame(matrix(0, nrow = n, ncol = p))
  names(data) <- nodes

  for (node in nodes) {
    parents <- dag[[node]]
    if (length(parents) == 0) {
      data[[node]] <- rnorm(n, sd = 1.0)
    } else {
      # Sum non-linear contributions from each parent
      val <- rep(0, n)
      for (pa in parents) {
        f <- sample(funcs, 1)[[1]]
        sign_factor <- sample(c(-1, 1), 1)
        val <- val + sign_factor * f(data[[pa]])
      }
      # Add noise
      data[[node]] <- val + rnorm(n, sd = noise_sd)
    }
  }

  dag_to_string <- function(d) {
    e <- character(0)
    for (n in names(d)) {
      for (pa in d[[n]]) e <- c(e, paste0(pa, "->", n))
    }
    if (length(e) == 0) return("<empty>")
    paste(sort(e), collapse = ";")
  }

  true_str <- dag_to_string(dag)
  seen_strs <- new.env(parent = emptyenv())
  seen_strs[[true_str]] <- TRUE

  # Generate (K-1) competitors as *fully random* DAGs over the same nodes: a
  # random topological order plus random edges. These are structurally unrelated
  # to the true DAG, so this asks whether the true structure stands out from an
  # arbitrary field rather than from near-identical twins.
  random_dag <- function() {
    perm <- sample(nodes)
    d <- setNames(rep(list(character(0)), length(nodes)), nodes)
    for (j in 2:length(perm)) {
      for (i in 1:(j - 1)) {
        if (runif(1) < edge_prob) d[[perm[j]]] <- c(d[[perm[j]]], perm[i])
      }
    }
    d
  }
  alt_dags <- list()
  attempts <- 0L
  while (length(alt_dags) < (n_candidate_dags - 1) && attempts < 10000L) {
    attempts <- attempts + 1L
    new_dag  <- random_dag()
    new_str  <- dag_to_string(new_dag)
    if (is.null(seen_strs[[new_str]])) {          # random_dag is acyclic by construction
      seen_strs[[new_str]] <- TRUE
      alt_dags[[length(alt_dags) + 1]] <- new_dag
    }
  }

  candidate_dags <- c(list(dag), alt_dags)

  list(dag = dag, data = data, p = p, candidates = candidate_dags)
}

# --- Run Benchmark ---

results <- list()
pb <- progress::progress_bar$new(total = n_iters, format = "  Simulating [:bar] :percent eta: :eta")

for (i in 1:n_iters) {
  pb$tick()

  p <- sample(nodes_range, 1)
  n_samples <- sample(n_samples_range, 1)
  sim <- generate_dag_and_data(p, n_samples, edge_prob, noise_sd)

  # Run discovery by passing the explicit list of DAGs (local package version).
  res <- suppressMessages(kagu_discover(sim$data, dags = sim$candidates))

  # Locate the true DAG by its edge set (the summary references DAGs by id now,
  # so we match on structure rather than a formatted edge string) and read off
  # its rank and posterior probability.
  true_edges <- kagu:::.dag_edges(sim$dag)
  true_i     <- which(vapply(res$dags, function(d)
                      identical(kagu:::.dag_edges(d), true_edges), logical(1)))[1]
  if (is.na(true_i)) {
    rank <- NA
    prob <- 0
  } else {
    ord  <- order(res$prob, decreasing = TRUE)
    rank <- which(ord == true_i)
    prob <- res$prob[true_i]
  }

  results[[i]] <- data.frame(
    iter = i,
    nodes = p,
    n_samples = n_samples,
    n_edges = length(unlist(sim$dag)),
    true_rank = rank,
    true_prob = prob,
    top_1_match = (rank == 1),
    top_10_pct_match = (rank <= ceiling(n_candidate_dags * 0.10)),
    top_20_pct_match = (rank <= ceiling(n_candidate_dags * 0.20))
  )
}

res_df <- bind_rows(results)

# --- Summary ---
cat("\n\n=== Benchmark Results ===\n")
res_df %>%
  group_by(nodes) %>%
  summarise(
    n_runs = n(),
    top_1_accuracy = mean(top_1_match, na.rm = TRUE),
    top_10_pct_accuracy = mean(top_10_pct_match, na.rm = TRUE),
    top_20_pct_accuracy = mean(top_20_pct_match, na.rm = TRUE),
    avg_rank = mean(true_rank, na.rm = TRUE),
    median_rank = median(true_rank, na.rm = TRUE),
    avg_prob = mean(true_prob, na.rm = TRUE)
  ) %>%
  print()

cat("\nOverall Top-1 Accuracy:", mean(res_df$top_1_match, na.rm = TRUE), "\n")
cat("Overall Top-10% Accuracy:", mean(res_df$top_10_pct_match, na.rm = TRUE), "\n")
cat("Overall Top-20% Accuracy:", mean(res_df$top_20_pct_match, na.rm = TRUE), "\n")
cat("Overall Average Rank:", mean(res_df$true_rank, na.rm = TRUE), "\n")

# Save detailed results
if (!dir.exists("outputs")) dir.create("outputs")
write.csv(res_df, "outputs/benchmark_discovery_results.csv", row.names = FALSE)
cli::cli_alert_success("Benchmark complete. Results saved to outputs/benchmark_discovery_results.csv")
