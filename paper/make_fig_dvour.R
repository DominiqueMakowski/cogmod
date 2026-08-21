# Generates figures/fig_dvour.png, used in Example 2 of paper.qmd.
#
# Reads the mixed shifted-LogNormal fitted by refit_reliability.R (and by the
# `Assessing Reliability` vignette), so it does not refit anything. Run from the
# `paper/` directory.
#
#   Rscript make_fig_dvour.R
#
# WHAT THE FIGURE HAS TO SHOW
# ---------------------------
# D-vour is a ratio of two spreads, and a table of ratios hides both of them.
# Each panel therefore draws the two quantities the ratio is made of:
#
#   - BETWEEN participants: how far apart the 30 estimates are. This is the
#     spread of the points along the x axis, and the marginal density drawn
#     above them.
#   - WITHIN participants: how sure we are of any one of them. This is the
#     width of each participant's interval.
#
# A parameter is usable as an individual score when the first is large relative
# to the second - when the points spread further than the bars are wide. That
# comparison is visible at a glance here and invisible in a table.

library(dplyr)
library(ggplot2)
library(brms)
library(cogmod)
library(modelbased)
library(performance)
library(stringr)


m <- readRDS(file.path("..", "vignettes", "models", "m_reliability.rds"))

random <- estimate_grouplevel(m)
dv <- as.data.frame(performance_dvour(random))

d <- as.data.frame(random)

# Readable panel names, ordered so the story reads left to right and top to
# bottom: the response formula first, then the two dpars that act as negative
# controls.
comp_lab <- c(conditional = "mu", sigma = "sigma", ndt = "ndt")
par_lab <- c(
  Intercept = "Intercept (Cond. A)", ConditionB = "Condition B (vs. A)",
  ConditionC = "Condition C (vs. A)"
)

d$comp <- factor(comp_lab[d$Component], levels = comp_lab)
d$par <- factor(par_lab[d$Parameter], levels = par_lab)
dv$comp <- factor(comp_lab[dv$Component], levels = comp_lab)
dv$par <- factor(par_lab[dv$Parameter], levels = par_lab)

# Order participants by their estimate *within each panel*, so the between-
# participant spread reads as a slope rather than as noise.
d <- d |>
  group_by(comp, par) |>
  mutate(rank = rank(Median, ties.method = "first")) |>
  ungroup()

n_p <- length(unique(d$Level))

# ---------------------------------------------------------------------------
# The marginal density of the point estimates, drawn as a band above the
# participants. Scaled to a fixed height per panel so that panels with very
# different x ranges stay comparable - the shape is the message, not the area.
# ---------------------------------------------------------------------------

BAND_LO <- n_p + 2
BAND_HI <- n_p + 11

dens <- d |>
  group_by(comp, par) |>
  group_modify(function(g, key) {
    dd <- stats::density(g$Median, adjust = 1.1)
    data.frame(x = dd$x, h = dd$y / max(dd$y))
  }) |>
  ungroup() |>
  mutate(y = BAND_LO + h * (BAND_HI - BAND_LO))

# D-vour label, printed in the top-left of each panel.
lab <- dv |>
  transmute(comp, par,
    txt = sprintf("D-vour = %.2f", D_vour),
    strong = D_vour > 0.5
  )

p <- d |>
  ggplot(aes(x = Median, y = rank)) +
  geom_vline(
    xintercept = 0, linetype = "dashed", colour = "grey",
    linewidth = 0.3
  ) +
  # within-participant uncertainty
  geom_pointrange(aes(xmin = CI_low, xmax = CI_high, color = par),
    size = 0.3,
    linewidth = 0.32
  ) +
  # between-participant spread
  geom_ribbon(
    data = dens, aes(x = x, ymin = BAND_LO, ymax = y, fill = par),
    inherit.aes = FALSE, alpha = 0.45
  ) +
  # geom_hline(yintercept = BAND_LO - 0.8, colour = "grey80", linewidth = 0.25) +
  geom_text(
    data = lab, aes(label = txt, colour = strong),
    x = -Inf, y = BAND_HI + 2.5, hjust = -0.06, vjust = 1,
    size = 2.6, fontface = "bold", inherit.aes = FALSE
  ) +
  scale_colour_manual(
    values = c(
      `TRUE` = "#4CAF50", `FALSE` = "red",
      "Intercept (Cond. A)" = "#9C27B0", "Condition B (vs. A)" = "#2196F3", "Condition C (vs. A)" = "#FF5722"
    ),
    guide = "none"
  ) +
  scale_fill_manual(
    values = c("Intercept (Cond. A)" = "#9C27B0", "Condition B (vs. A)" = "#2196F3", "Condition C (vs. A)" = "#FF5722"),
    guide = "none"
  ) +
  scale_y_continuous(limits = c(0, BAND_HI + 4), expand = c(0, 0)) +
  facet_grid(comp ~ par, scales = "free_x", switch = "y") +
  labs(
    x = "Participant deviation from the group",
    y = NULL,
    title = "Interindividual vs. Intraindividual Variability",
    subtitle = paste(
      "Each point is one participants, with its 95% interval (intraindividual variability);",
      "the marginal distribution is the spread of those points (interindividual variability)"
    )
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold", size = 9.5),
    strip.placement = "outside",
    plot.title = element_text(face = "bold", size = 11, colour = "#22252E"),
    plot.subtitle = element_text(face = "italic", size = 8.4, colour = "#6B7280"),
    panel.spacing.x = unit(0.7, "lines"),
    panel.spacing.y = unit(0.5, "lines")
  )
p
ggsave("figures/fig_dvour.png", p,
  width = 6.5, height = 6.5, dpi = 300,
  bg = "white"
)
cat("wrote figures/fig_dvour.png\n")

# The numbers the text quotes, so they cannot drift from the figure.
cat("\n--- D-vour ---\n")
print(dv[, c("Component", "Parameter", "D_vour")], digits = 3, row.names = FALSE)
cat("\n--- between-participant SD of the point estimates ---\n")
print(
  as.data.frame(d |> group_by(comp, par) |>
    summarise(
      sd_between = sd(Median),
      mean_ci_width = mean(CI_high - CI_low), .groups = "drop"
    )),
  digits = 3, row.names = FALSE
)
