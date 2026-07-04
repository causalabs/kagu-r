library(kagu)
library(dplyr)
library(ggplot2)

# Ensure reproducibility
# set.seed(1)

# --- Configuration ---
n_iters <- 50
n_samples <- 300
nodes_range <- 10:20
n_candidate_dags <- 100
edge_prob <- 0.15 # Keep sparsity reasonable for 20 nodes
noise_sd <- 0.5

cli::cli_alert_info("Starting discovery benchmark: {n_iters} random scenarios")
cli::cli_bullets(c(
  "*" = "Samples per scenario: {n_samples}",
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

  # Format the true edges string to match Kagu's formatting
  # Kagu sorts parents alphabetically and joins by semicolon
  # e.g., "x1->x2; x1->x3"
  edges <- character(0)
  for (node in nodes) {
    for (pa in dag[[node]]) {
      edges <- c(edges, paste0(pa, "→", node)) # using the right arrow char
    }
  }
  edges <- sort(edges)
  true_edges_str <- if (length(edges) > 0) paste(edges, collapse = "; ") else "(empty graph)"
  
  dag_to_string <- function(d) {
    e <- character(0)
    for (n in names(d)) {
      for (pa in d[[n]]) e <- c(e, paste0(pa, "->", n))
    }
    if (length(e) == 0) return("")
    paste(sort(e), collapse = ";")
  }

  true_str <- dag_to_string(dag)
  seen_strs <- new.env(parent = emptyenv())
  seen_strs[[true_str]] <- TRUE

  # Generate (K-1) alternative maximally difficult DAGs
  # We do this by randomly perturbing the true DAG by exactly 1 or 2 edges (add, remove, reverse)
  alt_dags <- list()
  while (length(alt_dags) < (n_candidate_dags - 1)) {
    new_dag <- dag
    
    n_changes <- sample(1:2, 1)
    for (c in 1:n_changes) {
      from <- sample(nodes, 1)
      to <- sample(nodes, 1)
      if (from != to) {
        if (from %in% new_dag[[to]]) {
          if (runif(1) < 0.5) {
            new_dag[[to]] <- setdiff(new_dag[[to]], from) # Remove
          } else {
            new_dag[[to]] <- setdiff(new_dag[[to]], from) # Reverse
            new_dag[[from]] <- c(new_dag[[from]], to)
          }
        } else if (to %in% new_dag[[from]]) {
          new_dag[[from]] <- setdiff(new_dag[[from]], to) # Reverse
          new_dag[[to]] <- c(new_dag[[to]], from)
        } else {
          new_dag[[to]] <- c(new_dag[[to]], from)         # Add
        }
      }
    }
    
    # Check cyclicity and uniqueness
    if (kagu:::.is_acyclic(new_dag)) {
      new_str <- dag_to_string(new_dag)
      if (is.null(seen_strs[[new_str]])) {
        seen_strs[[new_str]] <- TRUE
        alt_dags[[length(alt_dags) + 1]] <- new_dag
      }
    }
  }
  
  candidate_dags <- c(list(dag), alt_dags)
  
  list(dag = dag, data = data, true_edges_str = true_edges_str, p = p, candidates = candidate_dags)
}

# --- Run Benchmark ---

results <- list()
pb <- progress::progress_bar$new(total = n_iters, format = "  Simulating [:bar] :percent eta: :eta")

for (i in 1:n_iters) {
  pb$tick()

  p <- sample(nodes_range, 1)
  sim <- generate_dag_and_data(p, n_samples, edge_prob, noise_sd)

  # Run discovery by passing the explicit list of Dags
  # We must use devtools::load_all() internally since the script relies on the package namespace
  # Actually, we already loaded the library at the top, but the new code might not be installed.
  # Let's ensure it calls the local version
  res <- withr::with_options(list(kagu.quiet = TRUE), kagu_discover(sim$data, dags = sim$candidates))

  summ <- res$summary(top_n = res$n_models)

  true_edges <- sim$true_edges_str

  # Find where the true DAG ranks
  match_idx <- which(summ$edges == true_edges)

  if (length(match_idx) == 0) {
    rank <- NA
    prob <- 0
  } else {
    rank <- summ$rank[match_idx]
    prob <- summ$posterior_prob[match_idx]
  }

  results[[i]] <- data.frame(
    iter = i,
    nodes = p,
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
