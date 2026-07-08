# Publication figure: the three candidate DAGs from the ecology case study.
# Run from the package root:  Rscript figures/candidate_dags_figure.R
suppressMessages({ library(ggplot2); library(patchwork) })

# ---- The three candidate structures (shared node layout) --------------------
dag_true <- list(age = c(), sex = c(), sociality = "age",
                 food_sharing = "sociality",
                 condition = c("age", "sex", "sociality", "food_sharing"))
dag_alt1 <- list(age = c(), sex = c(), food_sharing = c(),
                 sociality = c("age", "food_sharing"),
                 condition = c("age", "sex", "sociality", "food_sharing"))
dag_alt2 <- list(age = c(), sex = c(), condition = c("age", "sex"),
                 food_sharing = "condition",
                 sociality = c("age", "food_sharing", "condition"))

pos <- list(age = c(0, 0), sex = c(0, 2), sociality = c(1, 0),
            food_sharing = c(2, 0), condition = c(1, 2))
labels <- c(age = "Age", sex = "Sex", sociality = "Sociality",
            food_sharing = "Food sharing", condition = "Condition")

draw_dag <- function(dag, title) {
  nd <- data.frame(name = names(pos),
                   x =  vapply(pos, `[`, numeric(1), 2),
                   y = -vapply(pos, `[`, numeric(1), 1))
  nd$lab <- labels[nd$name]
  ed <- do.call(rbind, lapply(names(dag), function(ch) {
    pa <- dag[[ch]]; if (!length(pa)) return(NULL)
    data.frame(from = pa, to = ch, stringsAsFactors = FALSE)
  }))
  ed$x  <- nd$x[match(ed$from, nd$name)]; ed$y  <- nd$y[match(ed$from, nd$name)]
  ed$xe <- nd$x[match(ed$to,   nd$name)]; ed$ye <- nd$y[match(ed$to,   nd$name)]
  dx <- ed$xe - ed$x; dy <- ed$ye - ed$y; L <- sqrt(dx^2 + dy^2); sh <- 0.34
  ed$x  <- ed$x  + sh * dx / L; ed$y  <- ed$y  + sh * dy / L
  ed$xe <- ed$xe - sh * dx / L; ed$ye <- ed$ye - sh * dy / L

  # Labels beside each node (left column -> left, right nodes -> right) so they
  # never sit on the vertical chain arrows.
  nd$xlab <- ifelse(nd$x <= 1, nd$x - 0.24, nd$x + 0.24)
  nd$hj   <- ifelse(nd$x <= 1, 1, 0)

  ggplot() +
    geom_segment(data = ed, aes(x, y, xend = xe, yend = ye),
                 arrow = grid::arrow(length = grid::unit(7, "pt"), type = "closed"),
                 linewidth = 0.5, colour = "grey15") +
    geom_point(data = nd, aes(x, y), size = 9, shape = 21,
               fill = "white", colour = "black", stroke = 1.1) +
    geom_text(data = nd, aes(x = xlab, y = y, label = lab), hjust = nd$hj,
              size = 3.1) +
    labs(title = title) +
    coord_equal(clip = "off") +
    scale_x_continuous(expand = expansion(mult = c(0.42, 0.3))) +
    scale_y_continuous(expand = expansion(mult = 0.16)) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13,
                                    margin = margin(b = 4)),
          plot.margin = margin(4, 10, 4, 10))
}

fig <- draw_dag(dag_true, "True") +
       draw_dag(dag_alt1, "Alternate 1") +
       draw_dag(dag_alt2, "Alternate 2") +
       plot_layout(nrow = 1)

ggsave("figures/candidate_dags.png", fig, width = 10, height = 2.9, dpi = 300, bg = "white")
ggsave("figures/candidate_dags.pdf", fig, width = 10, height = 2.9, bg = "white")
cat("saved to figures/\n")
