# Generates the static figures used in paper.qmd
# Run from the `paper/` directory.
#
# The rest of the pipeline is Python, run after this script:
#   make_fig_approaches.py   figures/fig_approaches.png  (reads the csv below)
#   make_fig_overview.py     figures/fig_overview.png    (embeds the logo)
#   add_figure_titles.py     puts headings on the two vignette figures

library(ggplot2)
library(patchwork)
library(lme4)
library(cogmod)

dir.create("figures", showWarnings = FALSE)

source("wagenmakers.R") # reconstructs the data from rtdists::speed_acc

# Intermediate consumed by make_fig_approaches.py
write.csv(wagenmakers(), "wagenmakers2008.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# Figure 1: two conditions with identical means but radically different shapes.
# Uses the bundled `badlm` dataset, as does the opening of the RT models
# vignette (see ?badlm and data-raw/badlm.R).
# ---------------------------------------------------------------------------

sim <- cogmod::badlm

m <- lmer(RT ~ Condition + (1 | Participant), data = sim)
print(summary(m))

means <- aggregate(RT ~ Condition, data = sim, FUN = mean)
print(means)

p1 <- ggplot(sim, aes(x = RT, fill = Condition)) +
  geom_histogram(bins = 100, alpha = 0.75, position = "identity") +
  geom_vline(
    data = means, aes(xintercept = RT, color = Condition),
    linetype = "dashed", linewidth = 0.8, show.legend = FALSE
  ) +
  scale_fill_manual(values = c(A = "#FF9800", B = "#3F51B5")) +
  scale_color_manual(values = c(A = "#E65100", B = "#1A237E")) +
  coord_cartesian(xlim = c(0, 2.5)) +
  labs(x = "Reaction Time (s)", y = "Count") +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = c(0.9, 0.8),
    panel.grid.minor = element_blank()
  )

# ---------------------------------------------------------------------------
# The same data under an ex-Gaussian likelihood. Only the family and the two
# extra formulas differ from the linear mixed model above. Feeds panel C of
# the figure and Table 1 (@tbl-badlm).
# ---------------------------------------------------------------------------

library(brms)

f_eg <- bf(
  RT ~ Condition + (1 | Participant),
  sigma ~ Condition,
  tau ~ Condition,
  family = rt_exgaussian()
)

# Condition B has essentially no exponential component, so `tau` sits at the
# boundary of its support and the default adapt_delta produces ~1% divergent
# transitions. Raising it to 0.999 brings that down to ~0.1% and roughly
# doubles the sampling time; the estimates are unchanged. Note that running
# *longer* chains does not help, since the divergences come from rare
# excursions into the tau -> 0 tail rather than from too few draws.
m_eg <- brm(
  f_eg,
  data = sim,
  stanvars = rt_exgaussian_stanvars(),
  init = 0,
  chains = 4, cores = 4, iter = 2000, seed = 1,
  control = list(adapt_delta = 0.999, max_treedepth = 12),
  backend = "cmdstanr",
  # Cached so that re-rendering the figure does not refit; `on_change`
  # invalidates the cache if the formula, data or family are edited.
  file = "m_eg_badlm", file_refit = "on_change"
)
print(m_eg)

# Response scale (softplus link on all three parameters), in milliseconds.
d <- brms::as_draws_df(m_eg)
softplus <- function(x) log1p(exp(x))
cell <- function(v) sprintf("%.3f [%.3f, %.3f]", mean(v), quantile(v, .025), quantile(v, .975))

pars <- list(
  mu = list(
    A = softplus(d$b_Intercept),
    B = softplus(d$b_Intercept + d$b_ConditionB)
  ),
  sigma = list(
    A = softplus(d$b_sigma_Intercept),
    B = softplus(d$b_sigma_Intercept + d$b_sigma_ConditionB)
  ),
  tau = list(
    A = softplus(d$b_tau_Intercept),
    B = softplus(d$b_tau_Intercept + d$b_tau_ConditionB)
  )
)
pars$`mu + tau` <- list(A = pars$mu$A + pars$tau$A, B = pars$mu$B + pars$tau$B)

for (nm in names(pars)) {
  a <- pars[[nm]]$A
  b <- pars[[nm]]$B
  cat(sprintf(
    "%-9s A = %-24s B = %-24s diff = %s\n",
    nm, cell(a), cell(b), cell(b - a)
  ))
}

# ---------------------------------------------------------------------------
# Panels B and C: the condition effect as each model reports it. Green = the
# interval excludes zero, red = it does not. The linear model has a single red
# row because the mean is the only quantity it estimates, and the mean is the
# one quantity in which these two conditions genuinely do not differ.
# ---------------------------------------------------------------------------

# Everything below is in seconds, matching the x axis of panels A and B.

# Linear mixed model: Wald interval on the condition slope.
lmm_est <- lme4::fixef(m)[["ConditionB"]]
lmm_se <- sqrt(diag(vcov(m)))[["ConditionB"]]
lmm_ci <- lmm_est + c(-1, 1) * 1.96 * lmm_se

# Ex-Gaussian: posterior of the B - A difference on each parameter.
eg_diff <- lapply(pars, function(p) p$B - p$A)

fmt <- function(est, lo, hi) {
  sprintf("%+.3f [%+.3f, %+.3f]", est, lo, hi)
}
row_from_draws <- function(label, draws) {
  ci <- unname(quantile(draws, c(.025, .975)))
  data.frame(
    label = label,
    value = fmt(mean(draws), ci[1], ci[2]),
    excludes_zero = ci[1] > 0 || ci[2] < 0
  )
}

tbl_lmm <- data.frame(
  label = "Condition B vs. A",
  value = fmt(lmm_est, lmm_ci[1], lmm_ci[2]),
  excludes_zero = lmm_ci[1] > 0 || lmm_ci[2] < 0
)

# Only the three actual parameters of the family; `mu + tau` is the implied
# mean rather than a parameter, so it is reported in @tbl-badlm instead.
tbl_eg <- rbind(
  row_from_draws("μ  (Gaussian bulk)", eg_diff$mu),
  row_from_draws("σ  (Gaussian SD)", eg_diff$sigma),
  row_from_draws("τ  (exponential tail)", eg_diff$tau)
)

# Both panels are laid out on the same row grid so their rows line up.
n_rows <- max(nrow(tbl_lmm), nrow(tbl_eg))

panel_table <- function(df, title, subtitle, note = NULL) {
  df$row <- seq_len(nrow(df))
  p <- ggplot(df) +
    geom_rect(
      aes(
        xmin = 0, xmax = 1,
        ymin = -row - 0.40, ymax = -row + 0.40,
        fill = excludes_zero
      ),
      alpha = 0.18, colour = NA
    ) +
    geom_text(aes(x = 0.035, y = -row, label = label),
      hjust = 0, size = 3.2, colour = "#22252E"
    ) +
    geom_text(aes(x = 0.965, y = -row, label = value, colour = excludes_zero),
      hjust = 1, size = 3.2, fontface = "bold"
    ) +
    scale_fill_manual(values = c(`TRUE` = "#2E7D32", `FALSE` = "#C62828"), guide = "none") +
    scale_colour_manual(values = c(`TRUE` = "#1B5E20", `FALSE` = "#B71C1C"), guide = "none") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-n_rows - 0.6, -0.3), expand = c(0, 0)) +
    labs(title = title, subtitle = subtitle) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(
        face = "bold", size = 10.5, colour = "#22252E",
        margin = margin(b = 1)
      ),
      plot.subtitle = element_text(
        size = 8.6, colour = "#6B7280",
        margin = margin(b = 5)
      ),
      plot.margin = margin(2, 6, 2, 6)
    )

  # The linear model's panel is mostly empty, which is the point; a short note
  # says so rather than leaving the reader to infer it from the white space.
  if (!is.null(note)) {
    p <- p + annotate("text",
      x = 0.035, y = -(nrow(df) + n_rows) / 2 - 0.35, label = note,
      hjust = 0, vjust = 0.5, size = 2.9, fontface = "italic",
      colour = "#6B7280", lineheight = 1.15
    )
  }
  p
}

t1 <- panel_table(tbl_lmm, "Linear mixed model", "Effect of Condition (s, 95% CI)",
  note = paste("The mean is the only parameter in which",
    "an effect could have appeared.",
    sep = "\n"
  )
)
t2 <- panel_table(tbl_eg, "Ex-Gaussian model", "Effect of Condition (s, 95% CI)")

# Sizes are set so the heading prints at ~10.5 pt once the 6.5 in figure is
# scaled to the width it occupies in the PDF.
p <- (p1 / (t1 | t2)) +
  plot_layout(heights = c(1, 0.52)) +
  plot_annotation(
    tag_levels = "A",
    title = "Identical means, different distributions",
    subtitle = "The linear model finds nothing; the ex-Gaussian finds three significant effects",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, colour = "#22252E"),
      plot.subtitle = element_text(face = "italic", size = 11.4, colour = "#6B7280")
    )
  )

ggsave("figures/fig_motivation.png", p, width = 6.5, height = 4.7, dpi = 300, bg = "white")

cat("\nDone.\n")
