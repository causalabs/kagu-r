Sys.setenv(RSTUDIO_PANDOC = "/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/x86_64")
suppressMessages({
  devtools::load_all("/Users/jordan/Documents/kagu-r", quiet = TRUE)
  library(ggplot2); library(patchwork)
})

# ---- Simulate: sociality -> condition, effect diverging by sex --------------
set.seed(11)
n         <- 500
sex       <- rbinom(n, 1, 0.5)                     # 0 = female, 1 = male
sociality <- 0.4 * sex + rnorm(n)
slope     <- ifelse(sex == 0, 1.5, -1.5)
condition <- slope * sociality + 0.5 * sex + rnorm(n, sd = 0.6)
df <- data.frame(sex = sex, sociality = sociality, condition = condition)

model <- KaguModel$new(list(sex = c(), sociality = "sex",
                            condition = c("sociality", "sex")))
model$fit(df)

pal <- c(Female = "#D55E00", Male = "#0072B2")

# ---- Right panel: sweeps stratified by sex ----------------------------------
sweep_by_sex <- function(s) {
  d <- model$effects("sociality", "condition", sweep = TRUE, sweep_n = 40,
                     sweep_range = c(-2.5, 2.5), conditions = list(sex = s))$summary()
  d$Sex <- if (s == 0) "Female" else "Male"; d
}
sw <- rbind(sweep_by_sex(0), sweep_by_sex(1))

ends <- do.call(rbind, lapply(split(sw, sw$Sex), function(d) d[which.max(d$x), ]))
femY <- ends$mean[ends$Sex == "Female"]; malY <- ends$mean[ends$Sex == "Male"]

p_sweep <- ggplot(sw, aes(x, colour = Sex, fill = Sex)) +
  geom_ribbon(aes(ymin = hdi_lower, ymax = hdi_upper), alpha = 0.16, colour = NA) +
  geom_line(aes(y = mean), linewidth = 1.2) +
  annotate("text", x = 2.5, y = femY, label = "Female", colour = pal[["Female"]],
           hjust = 1, vjust = -0.6, fontface = "bold", size = 4.4) +
  annotate("text", x = 2.5, y = malY, label = "Male", colour = pal[["Male"]],
           hjust = 1, vjust = 1.5, fontface = "bold", size = 4.4) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_fill_manual(values = pal, guide = "none") +
  labs(x = "Sociality", y = "Condition") +
  theme_classic(base_size = 13) +
  theme(axis.title = element_text(size = 13),
        plot.margin = margin(6, 12, 6, 28))

# ---- Left panels: three candidate DAGs (custom, crisp arrows) ---------------
draw_dag <- function(dag, pos, title) {
  nd <- data.frame(name = names(pos),
                   x =  vapply(pos, `[`, numeric(1), 2),
                   y = -vapply(pos, `[`, numeric(1), 1))
  ed <- do.call(rbind, lapply(names(dag), function(ch) {
    pa <- dag[[ch]]; if (!length(pa)) return(NULL)
    data.frame(from = pa, to = ch, stringsAsFactors = FALSE)
  }))
  ed$x  <- nd$x[match(ed$from, nd$name)]; ed$y  <- nd$y[match(ed$from, nd$name)]
  ed$xe <- nd$x[match(ed$to,   nd$name)]; ed$ye <- nd$y[match(ed$to,   nd$name)]
  dx <- ed$xe - ed$x; dy <- ed$ye - ed$y; L <- sqrt(dx^2 + dy^2); sh <- 0.31
  ed$x  <- ed$x  + sh * dx / L; ed$y  <- ed$y  + sh * dy / L
  ed$xe <- ed$xe - sh * dx / L; ed$ye <- ed$ye - sh * dy / L

  ggplot() +
    geom_segment(data = ed, aes(x, y, xend = xe, yend = ye),
                 arrow = grid::arrow(length = grid::unit(8, "pt"), type = "closed"),
                 linewidth = 0.55, colour = "grey15") +
    geom_point(data = nd, aes(x, y), size = 10, shape = 21,
               fill = "white", colour = "black", stroke = 1.1) +
    geom_text(data = nd, aes(x, y - 0.46, label = name), size = 3.3) +
    labs(title = title) +
    coord_equal(clip = "off") +
    scale_x_continuous(expand = expansion(mult = 0.28)) +
    scale_y_continuous(expand = expansion(mult = 0.34)) +
    theme_void(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12,
                                    margin = margin(b = 1)),
          plot.margin = margin(2, 6, 2, 6))
}

pos <- list(sex = c(0, 0), sociality = c(1, -1), condition = c(1, 1))
p_true <- draw_dag(list(sex = c(), sociality = "sex",
                        condition = c("sociality", "sex")), pos, "True")
p_a1   <- draw_dag(list(sex = c(), sociality = c("sex", "condition"),
                        condition = "sex"), pos, "Alternate 1")
p_a2   <- draw_dag(list(sex = c(), sociality = "sex",
                        condition = "sociality"), pos, "Alternate 2")

dag_col <- p_true / p_a1 / p_a2

fig <- (wrap_elements(dag_col) | p_sweep) +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 15))

ggsave("figures/kagu_pubfig.png", fig, width = 9.5, height = 6.3, dpi = 300, bg = "white")
ggsave("figures/kagu_pubfig.pdf", fig, width = 9.5, height = 6.3, bg = "white")
cat("saved to figures/\n")
