Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/x86_64")
suppressMessages({
  devtools::load_all("/Users/jordan/Documents/kagu-r", quiet = TRUE)
  library(ggplot2); library(patchwork)
})

# ---------------------------------------------------------------------------
# Ecology case-study data. n = 50: small and noisy enough that the weak
# age -> sociality edge cannot be pinned down, so discovery stays honestly
# uncertain while still preferring the true structure.
# ---------------------------------------------------------------------------
set.seed(42)
n            <- 50
sex          <- sample(c(-1, 1), n, replace = TRUE)
age          <- rnorm(n)
sociality    <- 0.5 * age + rnorm(n, sd = 1.5)
food_sharing <- 0.8 * sociality + rnorm(n, sd = 1.5)
condition    <- 0.5 * age + 0.5 * sex + 0.8 * food_sharing -
                1.0 * sex * sociality + rnorm(n, sd = 0.5)
df <- data.frame(age, sex, sociality, food_sharing, condition)

fc <- c("age", "sex", "sociality", "food_sharing")
dag_true <- list(age = c(), sex = c(), sociality = "age",
                 food_sharing = "sociality", condition = fc)
dag_alt1 <- list(age = c(), sex = c(), sociality = c(),           # sociality has no cause
                 food_sharing = "sociality", condition = fc)
dag_alt2 <- list(age = c(), sex = c(), condition = c("age", "sex"),  # condition drives sociality
                 food_sharing = "condition",
                 sociality = c("age", "food_sharing", "condition"))
cands <- list(dag_true, dag_alt1, dag_alt2)

res  <- KaguModel$discover(df, dags = cands)
prob <- vapply(seq_along(cands), function(i) {
  L <- res$log_marglik; exp(L[i] - max(L)) / sum(exp(L - max(L)))
}, numeric(1))

# ---------------------------------------------------------------------------
# (a) three candidate DAGs, each with a posterior-probability bar beneath
# ---------------------------------------------------------------------------
pos <- list(age = c(0, 0), sex = c(0, 2), sociality = c(1, 0),
            food_sharing = c(2, 0), condition = c(1, 2))
labs <- c(age = "Age", sex = "Sex", sociality = "Sociality",
          food_sharing = "Food sharing", condition = "Condition")
teal <- "#26a69a"; slate <- "#90a4ae"

draw_dag <- function(dag, title) {
  nd <- data.frame(name = names(pos),
                   x =  vapply(pos, `[`, numeric(1), 2),
                   y = -vapply(pos, `[`, numeric(1), 1))
  nd$lab <- labs[nd$name]
  ed <- do.call(rbind, lapply(names(dag), function(ch) {
    pa <- dag[[ch]]; if (!length(pa)) return(NULL)
    data.frame(from = pa, to = ch, stringsAsFactors = FALSE)
  }))
  ed$x  <- nd$x[match(ed$from, nd$name)]; ed$y  <- nd$y[match(ed$from, nd$name)]
  ed$xe <- nd$x[match(ed$to,   nd$name)]; ed$ye <- nd$y[match(ed$to,   nd$name)]
  dx <- ed$xe - ed$x; dy <- ed$ye - ed$y; L <- sqrt(dx^2 + dy^2); sh <- 0.34
  ed$x  <- ed$x  + sh * dx / L; ed$y  <- ed$y  + sh * dy / L
  ed$xe <- ed$xe - sh * dx / L; ed$ye <- ed$ye - sh * dy / L
  nd$xlab <- ifelse(nd$x <= 1, nd$x - 0.45, nd$x + 0.45)
  nd$hj   <- ifelse(nd$x <= 1, 1, 0)
  ggplot() +
    geom_segment(data = ed, aes(x, y, xend = xe, yend = ye),
                 arrow = grid::arrow(length = grid::unit(6, "pt"), type = "closed"),
                 linewidth = 0.45, colour = "grey15") +
    geom_point(data = nd, aes(x, y), size = 7.5, shape = 21,
               fill = "white", colour = "black", stroke = 1) +
    geom_text(data = nd, aes(x = xlab, y = y, label = lab), hjust = nd$hj, size = 2.8) +
    labs(title = title) +
    coord_equal(clip = "off") +
    scale_x_continuous(expand = expansion(mult = c(0.62, 0.46))) +
    scale_y_continuous(expand = expansion(mult = 0.16)) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12,
                                    margin = margin(b = 2)),
          plot.margin = margin(2, 8, 0, 8))
}

prob_bar <- function(p, fill) {
  ggplot() +
    geom_rect(aes(xmin = 0, xmax = 1, ymin = 0, ymax = 1), fill = "grey90") +
    geom_rect(aes(xmin = 0, xmax = max(p, 0.0025), ymin = 0, ymax = 1), fill = fill) +
    annotate("text", x = 0.5, y = 1.9, label = sprintf("%.0f%%", 100 * p),
             size = 4.1, fontface = "bold", colour = "grey20") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 2.6), clip = "off") +
    theme_void() +
    theme(plot.margin = margin(1, 22, 4, 22))
}

col <- function(dag, title, p, fill)
  (draw_dag(dag, title) / prob_bar(p, fill)) + plot_layout(heights = c(4.2, 1))

panel_a <- col(dag_true, "True structure", prob[1], teal) |
           col(dag_alt1, "Alternative 1",  prob[2], slate) |
           col(dag_alt2, "Alternative 2",  prob[3], slate)
panel_a <- panel_a + plot_annotation(
  title = "Posterior probability over candidate structures") &
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))

# ---------------------------------------------------------------------------
# (b) sweep: total effect of sociality on condition, stratified by sex
# ---------------------------------------------------------------------------
model <- KaguModel$new(dag_true); model$fit(df)
pal_b <- c(Female = "#D55E00", Male = "#0072B2")
sweep_by_sex <- function(s) {
  d <- model$effects("sociality", "condition", sweep = TRUE, sweep_n = 40,
                     sweep_range = c(-2, 2), conditions = list(sex = s))$summary()
  d$Sex <- if (s == -1) "Female" else "Male"; d
}
sw   <- rbind(sweep_by_sex(-1), sweep_by_sex(1))
ends <- do.call(rbind, lapply(split(sw, sw$Sex), function(d) d[which.max(d$x), ]))
femY <- ends$mean[ends$Sex == "Female"]; malY <- ends$mean[ends$Sex == "Male"]

panel_b <- ggplot(sw, aes(x, colour = Sex, fill = Sex)) +
  geom_ribbon(aes(ymin = hdi_lower, ymax = hdi_upper), alpha = 0.16, colour = NA) +
  geom_line(aes(y = mean), linewidth = 1.2) +
  annotate("text", x = 2, y = femY, label = "Female", colour = pal_b[["Female"]],
           hjust = 1, vjust = -0.7, fontface = "bold", size = 4.2) +
  annotate("text", x = 2, y = malY, label = "Male", colour = pal_b[["Male"]],
           hjust = 1, vjust = 1.6, fontface = "bold", size = 4.2) +
  scale_colour_manual(values = pal_b, guide = "none") +
  scale_fill_manual(values = pal_b, guide = "none") +
  labs(x = "Intervention: do(Sociality)", y = "Expected condition",
       title = "Sex-dependent effect recovered by the GP") +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
        plot.margin = margin(6, 14, 6, 10))

# ---------------------------------------------------------------------------
fig <- wrap_elements(panel_a) / panel_b +
  plot_layout(heights = c(1, 1.15)) +
  plot_annotation(tag_levels = list(c("(a)", "(b)"))) &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("figures/discovery_sweep.png", fig, width = 10, height = 7.4, dpi = 200, bg = "white")
ggsave("figures/discovery_sweep.pdf", fig, width = 10, height = 7.4, bg = "white")
message(sprintf("probs: True=%.1f%%  Alt1=%.1f%%  Alt2=%.2f%%",
                100*prob[1], 100*prob[2], 100*prob[3]))
message("wrote figures/discovery_sweep.{png,pdf}")
