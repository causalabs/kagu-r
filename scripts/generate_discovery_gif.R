library(kagu)
library(dplyr)
library(ggplot2)
library(gganimate)

# Ensure reproducibility
set.seed(42)

# We generate a dataset from a non-linear 4-node system.
# w causes x and y (non-linearly). x and y jointly cause z (linearly).
# Non-linear relationships break Markov equivalence, so Kagu's GP mechanism
# can uniquely identify the true causal direction as N grows.
max_n <- 500
w <- rnorm(max_n)
x <- sin(2 * w) + rnorm(max_n, sd = 0.5)
y <- cos(2 * w) + rnorm(max_n, sd = 0.5)
z <- 0.8 * x + 0.8 * y + rnorm(max_n, sd = 0.5)
df_full <- data.frame(w, x, y, z)

# Evaluate over an increasing sequence of sample sizes
ns <- seq(20, 500, by = 20)

cli::cli_alert_info("Fitting {length(ns)} subsets of data to track posterior convergence...")

# To keep the space manageable (543 -> 200 DAGs) we encode the light domain knowledge
# that w is temporally first (nothing causes it). 
disallowed <- list(c("x", "w"), c("y", "w"), c("z", "w"))

# Loop over N and collect all posterior probabilities
results <- do.call(rbind, lapply(ns, function(n) {
  df <- df_full[1:n, ]
  
  # Suppress fitting output for the animation loop
  res <- withr::with_options(list(kagu.quiet = TRUE), kagu_discover(df, disallowed = disallowed))
  
  probs <- res$probabilities()
  probs$N <- n
  probs
}))

# We want a stable x-axis and consistent bars. 
# Let's find the true DAG and the next best 4 DAGs across the entire sequence.
true_dag <- "w→x; w→y; x→z; y→z"

# Find highest scoring DAGs overall (other than the true one)
top_other <- results %>%
  filter(edges != true_dag) %>%
  group_by(edges) %>%
  summarize(max_prob = max(posterior_prob), .groups = "drop") %>%
  arrange(desc(max_prob)) %>%
  head(4) %>%
  pull(edges)

keep_dags <- c(true_dag, top_other)

# Group the rest into "Other"
plot_df <- results %>%
  mutate(
    dag_label = if_else(edges %in% keep_dags, edges, "Other"),
    is_true = if_else(edges == true_dag, "True DAG", "Other DAGs")
  ) %>%
  group_by(N, dag_label, is_true) %>%
  summarise(probability = sum(posterior_prob), .groups = "drop")

# Ensure ordering of factors so true DAG is always first, then the specific ones, then "Other"
plot_df$dag_label <- factor(plot_df$dag_label, levels = c(true_dag, top_other, "Other"))

cli::cli_alert_info("Animating the sequence using gganimate...")

# Create the plot
p <- ggplot(plot_df, aes(x = dag_label, y = probability, fill = is_true)) +
  geom_col() +
  scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
  scale_fill_manual(values = c("True DAG" = "#26a69a", "Other DAGs" = "gray70")) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16)
  ) +
  labs(
    title = "Causal Structure Discovery: N = {closest_state}",
    subtitle = "Posterior probability over DAGs converging to the true non-linear structure",
    x = "Candidate Structure",
    y = "Posterior Probability P(G | X)"
  ) +
  transition_states(N, transition_length = 2, state_length = 1) +
  ease_aes('cubic-in-out')

# Render and save the gif
if (!dir.exists("outputs")) dir.create("outputs")
animate(p, nframes = 150, fps = 20, width = 800, height = 500, res = 100,
        renderer = gifski_renderer("outputs/discovery_convergence.gif"))

cli::cli_alert_success("Saved animation to outputs/discovery_convergence.gif")
