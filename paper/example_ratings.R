# Example 4 of the paper: bounded ratings with cogmod_betagate().
# Run from the `paper/` directory. Produces figures/fig_ratings.png and the
# numbers quoted in the text.
#
# The CHOCO (Choice-Confidence) family is deliberately not covered by the
# paper - it is the subject of a separate paper of its own - so nothing here
# should grow to include it.
#
# THE SCENARIO
# ------------
# A taste test. Participants rate how much they enjoyed each of two dishes on a
# visual analogue slider from "hated it" (0) to "loved it" (1).
#
# One dish is a margherita pizza: nobody has strong feelings about it, and the
# ratings pile up in the middle. The other is heavy on coriander, which a
# substantial minority of people experience as tasting of soap - a real and
# partly heritable split. Its ratings are bimodal: half the room loves it, half
# hates it, and few people are indifferent.
#
# Both dishes have the *same mean rating*. A t-test therefore reports that the
# dishes are equally liked, which is true of the average and false of everyone.
# This is the opening figure's argument (@fig-motivation) transplanted to a
# bounded scale, and it is the reason a mean is a poor summary of a rating.

library(cogmod)
library(brms)
library(ggplot2)
library(patchwork)
library(ggnewscale)

# The two dishes have *identical* generating means (0.50), so any difference in
# the sample means is sampling error: with n = 250 per dish the SE of that
# difference is about 0.027, and across 400 seeds only 4% give a significant
# t-test. This seed is one whose realised means coincide closely, so the figure
# shows the typical case rather than a 1-in-25 draw.
set.seed(19)

n <- 250
true <- list(
  Margherita = list(mu = 0.50, phi = 10.0, pex = 0.02, bex = 0.5),
  Coriander  = list(mu = 0.50, phi = 0.6, pex = 0.40, bex = 0.5)
)

df <- do.call(rbind, lapply(names(true), function(k) {
  p <- true[[k]]
  data.frame(
    Dish = k,
    score = rcogmod_betagate(n, mu = p$mu, phi = p$phi, pex = p$pex, bex = p$bex)
  )
}))
df$Dish <- factor(df$Dish, levels = c("Margherita", "Coriander"))

# ---- What a mean comparison sees --------------------------------------------

cat("--- observed ---\n")
print(aggregate(score ~ Dish, df, function(v) {
  round(c(
    mean = mean(v), sd = sd(v), prop_0 = mean(v == 0),
    prop_1 = mean(v == 1), prop_mid = mean(v > 0.35 & v < 0.65)
  ), 3)
}))
tt <- t.test(score ~ Dish, data = df)
cat(sprintf(
  "\nt-test: t(%.0f) = %.2f, p = %.3f, diff = %+.3f [%+.3f, %+.3f]\n",
  tt$parameter, tt$statistic, tt$p.value,
  diff(rev(tt$estimate)), tt$conf.int[1], tt$conf.int[2]
))

# ---- The Gaussian model -----------------------------------------------------

m_norm <- brm(
  score ~ Dish,
  data = df,
  chains = 4, cores = 4, iter = 2000, seed = 1,
  backend = "cmdstanr",
  file = "m_ratings_normal", file_refit = "on_change"
)
print(m_norm)

# ---- The Beta-Gate model ----------------------------------------------------

# cogmod_priors() has nothing to add for the rating families yet - it covers the
# RT and choice families, whose parameters are in seconds - so the brms defaults
# stand here, and only the Stan code has to be attached.
f <- bf(
  score ~ Dish,
  phi ~ Dish,
  pex ~ Dish,
  bex ~ Dish,
  family = cogmod_betagate()
)

m <- brm(
  f,
  data = df,
  stanvars = cogmod_stanvars(f),
  init = 0,
  chains = 4, cores = 4, iter = 2000, seed = 1,
  control = list(adapt_delta = 0.95),
  backend = "cmdstanr",
  file = "m_betagate_ratings", file_refit = "on_change"
)
print(m)

# ---- Response-scale parameters ----------------------------------------------

d <- brms::as_draws_df(m)
softplus <- function(x) log1p(exp(x))
cell <- function(v) {
  sprintf("%.2f [%.2f, %.2f]", mean(v), quantile(v, .025), quantile(v, .975))
}

pars <- list(
  mu  = list(link = plogis, int = "b_Intercept", slope = "b_DishCoriander"),
  phi = list(link = softplus, int = "b_phi_Intercept", slope = "b_phi_DishCoriander"),
  pex = list(link = plogis, int = "b_pex_Intercept", slope = "b_pex_DishCoriander"),
  bex = list(link = plogis, int = "b_bex_Intercept", slope = "b_bex_DishCoriander")
)

cat("\n--- Beta-Gate, response scale ---\n")
est <- list()
for (nm in names(pars)) {
  p <- pars[[nm]]
  a <- p$link(d[[p$int]])
  b <- p$link(d[[p$int]] + d[[p$slope]])
  est[[nm]] <- list(A = a, B = b)
  cat(sprintf(
    "%-4s Margherita = %-20s Coriander = %-20s  (generated: %.2f -> %.2f)\n",
    nm, cell(a), cell(b), true$Margherita[[nm]], true$Coriander[[nm]]
  ))
}

# ---- Figure -----------------------------------------------------------------

means <- aggregate(score ~ Dish, df, mean)

pA <- ggplot() +

  # Margherita: yellow -> red
  geom_histogram(
    data = subset(df, Dish == "Margherita"),
    aes(x = score, fill = after_stat(x)),
    bins = 40,
    alpha = 0.8,
    colour = NA
  ) +
  scale_fill_gradient(
    low = "#FFF176",
    high = "#E53935",
    guide = "none"
  ) +
  ggnewscale::new_scale_fill() +

  # Coriander: cyan -> green
  geom_histogram(
    data = subset(df, Dish == "Coriander"),
    aes(x = score, fill = after_stat(x)),
    bins = 40,
    alpha = 0.8,
    colour = NA
  ) +
  scale_fill_gradient(
    low = "#00BCD4",
    high = "#43A047",
    guide = "none"
  ) +

  # Means
  geom_vline(
    data = means,
    aes(xintercept = score, colour = Dish),
    linetype = "dashed",
    linewidth = 0.8
  ) +
  scale_colour_manual(
    values = c(
      Margherita = "#FF9800",
      Coriander = "#20BFA9"
    ),
    guide = "none"
  ) +

  # Emoji legend
  annotate(
    "text",
    x = 0.64, y = Inf,
    label = "🍕 Margherita",
    colour = "#FF9800",
    hjust = 0,
    vjust = 1.5,
    size = 3
  ) +
  annotate(
    "text",
    x = 0.82, y = Inf,
    label = "🌿 Coriander",
    colour = "#43A047",
    hjust = 0,
    vjust = 1.5,
    size = 3
  ) +
  labs(
    x = "Enjoyment rating",
    y = "Count"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

fmt <- function(est, lo, hi) sprintf("%+.2f [%+.2f, %+.2f]", est, lo, hi)
row_from_draws <- function(label, draws) {
  ci <- unname(quantile(draws, c(.025, .975)))
  data.frame(
    label = label, value = fmt(mean(draws), ci[1], ci[2]),
    excludes_zero = ci[1] > 0 || ci[2] < 0
  )
}

norm_draws <- brms::as_draws_df(m_norm)$b_DishCoriander
tbl_norm <- row_from_draws("Mean rating", norm_draws)
tbl_bg <- rbind(
  row_from_draws("μ  (underlying tendency)", est$mu$B - est$mu$A),
  row_from_draws("φ  (consistency)", est$phi$B - est$phi$A),
  row_from_draws("pₑₓ  (extreme responding)", est$pex$B - est$pex$A),
  row_from_draws("bₑₓ  (direction of extremes)", est$bex$B - est$bex$A)
)

n_rows <- max(nrow(tbl_norm), nrow(tbl_bg))
panel_table <- function(df, title, subtitle, note = NULL) {
  df$row <- seq_len(nrow(df))
  p <- ggplot(df) +
    geom_rect(aes(
      xmin = 0, xmax = 1, ymin = -row - 0.40, ymax = -row + 0.40,
      fill = excludes_zero
    ), alpha = 0.18, colour = NA) +
    geom_text(aes(x = 0.02, y = -row, label = label),
      hjust = 0, size = 2.95,
      colour = "#22252E"
    ) +
    geom_text(aes(x = 0.985, y = -row, label = value, colour = excludes_zero),
      hjust = 1, size = 2.95, fontface = "bold"
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
      plot.subtitle = element_text(size = 8.6, colour = "#6B7280", margin = margin(b = 5)),
      plot.margin = margin(2, 6, 2, 6)
    )
  if (!is.null(note)) {
    p <- p + annotate("text",
      x = 0.035, y = -(nrow(df) + n_rows) / 2 - 0.35,
      label = note, hjust = 0, vjust = 0.5, size = 2.9,
      fontface = "italic", colour = "#6B7280", lineheight = 1.15
    )
  }
  p
}

t1 <- panel_table(tbl_norm, "Gaussian model", "Coriander − Margherita (95% CI)",
  note = paste("The mean is the only parameter",
    "in which an effect could appear.",
    sep = "\n"
  )
)
t2 <- panel_table(tbl_bg, "Beta-Gate model", "Coriander − Margherita (95% CI)")

# Build the second row first: `plot_layout()` binds to the composition it is
# added to, so writing `t1 | t2 + plot_layout(...)` would attach the widths to
# t2 alone and squeeze the four-row panel rather than widen it.
row2 <- (t1 | t2) + plot_layout(widths = c(0.75, 1.25))

p <- (pA / row2) +
  plot_layout(heights = c(1, 0.62)) +
  plot_annotation(
    tag_levels = "A",
    title = "Focusing on average ratings can be misleading",
    subtitle = "The margherita divides nobody; the coriander divides the room",
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, colour = "#22252E"),
      plot.subtitle = element_text(face = "italic", size = 11.4, colour = "#6B7280")
    )
  )
p
ggsave("figures/fig_ratings.png", p, width = 6.5, height = 5, dpi = 300, bg = "white")
cat("\nwrote figures/fig_ratings.png\n")
