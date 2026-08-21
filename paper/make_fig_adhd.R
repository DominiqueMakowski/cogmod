# Generates figures/fig_adhd.png, used in Example 2 of paper.qmd.
#
# Reads the objects saved by example_cnp.R (which must be run first, after
# download_cnp.R). Run from the `paper/` directory.
#
# Panel A: per-participant switch costs, defined as mean RT on switch trials
#          minus mean RT on non-switch trials.
# Panel B: response-scale decomposition from the shifted LogNormal model.

library(dplyr)
library(ggplot2)
library(brms)
library(cogmod)
library(patchwork)
library(ggdist)

DIR <- "cnp_data"

m <- readRDS(file.path(DIR, "m_ln.rds"))
d <- readRDS(file.path(DIR, "d.rds"))
sc <- readRDS(file.path(DIR, "sc.rds"))

cols <- c(CONTROL = "#00BCD4", ADHD = "#E91E63")

# -------------------------------------------------------------------------
# Panel A: conventional switch-cost index
# -------------------------------------------------------------------------

pA <- ggplot(
  sc,
  aes(x = 1000 * cost, fill = Diagnosis)
) +
  geom_density(alpha = .45, colour = NA) +
  geom_vline(
    data = sc |>
      group_by(Diagnosis) |>
      summarise(m = 1000 * mean(cost)),
    aes(
      xintercept = m,
      colour = Diagnosis
    ),
    linetype = "dashed",
    linewidth = .6
  ) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols) +
  labs(
    x = "Switch cost (ms)",
    y = NULL,
    title = "A. The conventional index",
    subtitle = "mean switch RT − mean non-switch RT"
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "italic"),
    legend.position = "top",
    legend.title = element_blank(),
    legend.key.size = unit(.3, "cm"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# -------------------------------------------------------------------------
# Response-scale decomposition
# -------------------------------------------------------------------------

p <- posterior::as_draws_df(m)

sp <- function(x) {
  log1p(exp(x))
}

# Non-decision time is on a log link and expressed directly in seconds.
ndt_c <- exp(p$b_ndt_Intercept)
ndt_a <- exp(
  p$b_ndt_Intercept +
    p$b_ndt_DiagnosisADHD
)

# Sigma is on a softplus link.
s_c <- sp(p$b_sigma_Intercept)
s_a <- sp(
  p$b_sigma_Intercept +
    p$b_sigma_DiagnosisADHD
)

# For a LogNormal distribution, the mean of the decision component
# is exp(mu + sigma^2 / 2).
dec_c <- exp(
  p$b_Intercept +
    s_c^2 / 2
)

dec_a <- exp(
  p$b_Intercept +
    p$b_DiagnosisADHD +
    s_a^2 / 2
)

eff <- bind_rows(
  tibble(
    par = "Total\nmean RT",
    v = 1000 * (
      (ndt_a + dec_a) -
        (ndt_c + dec_c)
    )
  ),
  tibble(
    par = "Non-decision\ntime",
    v = 1000 * (
      ndt_a - ndt_c
    )
  ),
  tibble(
    par = "Decision\ncomponent",
    v = 1000 * (
      dec_a - dec_c
    )
  )
) |>
  mutate(
    par = factor(
      par,
      levels = c(
        "Total\nmean RT",
        "Non-decision\ntime",
        "Decision\ncomponent"
      )
    )
  )

# -------------------------------------------------------------------------
# Panel B: model-based decomposition
# -------------------------------------------------------------------------

pB <- ggplot(
  eff,
  aes(x = v, y = par)
) +
  stat_halfeye(
    aes(fill = after_stat(x > 0)),
    .width = .95,
    normalize = "groups",
    slab_alpha = .75,
    point_size = 1.2,
    interval_size = .6
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey35"
  ) +
  scale_fill_manual(
    values = c(
      `TRUE` = "#E91E63",
      `FALSE` = "#00BCD4"
    ),
    guide = "none"
  ) +
  labs(
    x = "ADHD − control (ms)",
    y = NULL,
    title = "B. What the model separates",
    subtitle = "ADHD − control"
  ) +
  theme_classic(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "italic"),
  )

# -------------------------------------------------------------------------
# Combine and save
# -------------------------------------------------------------------------

g <- pA / pB

g

ggsave("figures/fig_adhd.png", g, width = 7.0, height = 5, dpi = 300, bg = "white")
cat("wrote figures/fig_adhd.png\n")
