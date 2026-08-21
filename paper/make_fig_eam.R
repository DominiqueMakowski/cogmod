# Generates figures/fig_eam.png, used in Example 1 of paper.qmd.
#
# Reads the five choice+RT models fitted in the `Decision Making Models`
# vignette (vignettes/models/*.rds, see vignettes/articles/decision_making.qmd), so it
# does not refit anything. Run from the `paper/` directory.
#
#   Rscript make_fig_eam.R
#
# Panel A: posterior predictive check, correct responses drawn upwards and
#          errors downwards, each scaled by its own predicted frequency.
# Panel B: what each architecture costs against what it buys.

library(brms)
library(cogmod)
library(ggplot2)
library(patchwork)

set.seed(1)

NDRAWS <- 100      # posterior predictive draws per model
NLINES <- 50       # of which this many are drawn as individual faint curves

VIG <- file.path("..", "vignettes", "models")

models <- list(
  LBA     = readRDS(file.path(VIG, "m_lba2.rds")),
  LNR     = readRDS(file.path(VIG, "m_lnr.rds")),
  `DDM-5` = readRDS(file.path(VIG, "m_ddm5.rds")),
  DDM     = readRDS(file.path(VIG, "m_ddm.rds")),
  RDM     = readRDS(file.path(VIG, "m_rdm.rds"))
)

# The data are taken from the fitted object rather than rebuilt here, so the
# figure cannot drift away from the models the vignette actually fitted
# (participants 1-3 of Wagenmakers et al., RT <= 2 s, correct and error trials).
df <- models$DDM$data
n_obs <- nrow(df)
pct_correct <- 100 * mean(df$Error == 0)
pct_error <- 100 * mean(df$Error == 1)

# Ordered by ELPD (best first); the levels drive every legend and axis below.
LEVELS <- names(models)
COLORS <- c(
  LBA = "#4CAF50", LNR = "#FF9800", `DDM-5` = "#2196F3",
  DDM = "#673AB7", RDM = "#E91E63"
)

# ---------------------------------------------------------------------------
# Posterior predictive draws. `posterior_predict()` returns the two outcomes
# interleaved (odd columns = RT, even columns = response), one pair per trial.
# ---------------------------------------------------------------------------

# `with_outliers()` matters here. By default `posterior_predict()` returns the
# decision process alone, while `log_lik()` - and so `loo()` - is the full
# mixture. Panel B ranks the models by LOO, so panel A has to check the same
# distribution that ranking was computed from.
pred <- lapply(LEVELS, function(nm) {
  pp <- brms::posterior_predict(cogmod::with_outliers(models[[nm]]), ndraws = NDRAWS)
  rt <- pp[, seq(1, 2 * n_obs, by = 2), drop = FALSE]
  resp <- pp[, seq(2, 2 * n_obs, by = 2), drop = FALSE]
  data.frame(
    Model = nm,
    draw = rep(seq_len(NDRAWS), times = n_obs),
    rt = as.vector(rt),
    response = as.vector(resp)
  )
})
names(pred) <- LEVELS
pred <- do.call(rbind, pred)
pred$Model <- factor(pred$Model, levels = LEVELS)

# A near-zero drift rate takes an arbitrarily long time to reach the threshold,
# so the races occasionally predict enormous RTs. The observed data were
# trimmed at 2 s, so predictions beyond 3 s are discarded throughout: in panel
# A a single extreme draw would flatten every curve (`stat_density()` spreads
# its evaluation grid over the whole x range), and it would move a whole draw's
# error mean by several hundred milliseconds in the statistic printed below.
cat(sprintf("dropped %d predicted RTs > 3 s (%.3f%%)\n",
            sum(pred$rt > 3), 100 * mean(pred$rt > 3)))
pred <- pred[pred$rt < 3, ]
pred_plot <- pred

# ---------------------------------------------------------------------------
# Panel A: correct upwards, errors downwards, both scaled by the *proportion*
# of that response, so the area under each curve is its predicted frequency.
# ---------------------------------------------------------------------------

correct <- pred_plot[pred_plot$response == 0, ]
error <- pred_plot[pred_plot$response == 1, ]
lines_draws <- seq_len(NLINES)
bw <- 0.02

pA <- ggplot(df, aes(x = RT)) +
  geom_histogram(
    data = df[df$Error == 0, ], aes(y = after_stat(count) / (n_obs * bw)),
    binwidth = bw, fill = "#2E7D32", alpha = 0.25
  ) +
  geom_histogram(
    data = df[df$Error == 1, ], aes(y = -after_stat(count) / (n_obs * bw)),
    binwidth = bw, fill = "#B71C1C", alpha = 0.25
  ) +
  # One faint curve per posterior draw
  geom_line(
    data = correct[correct$draw %in% lines_draws, ],
    aes(x = rt, y = after_stat(count) / n_obs, color = Model,
        group = interaction(Model, draw)),
    stat = "density", alpha = 0.045, linewidth = 0.25
  ) +
  geom_line(
    data = error[error$draw %in% lines_draws, ],
    aes(x = rt, y = -after_stat(count) / n_obs, color = Model,
        group = interaction(Model, draw)),
    stat = "density", alpha = 0.045, linewidth = 0.25
  ) +
  # Pooled over draws
  geom_line(
    data = correct, aes(x = rt, y = after_stat(count) / (n_obs * NDRAWS), color = Model),
    stat = "density", linewidth = 0.75
  ) +
  geom_line(
    data = error, aes(x = rt, y = -after_stat(count) / (n_obs * NDRAWS), color = Model),
    stat = "density", linewidth = 0.75
  ) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.3) +
  annotate("text", x = 1.29, y = 2.6, hjust = 1, size = 2.9, fontface = "bold",
           color = "#2E7D32", label = sprintf("Correct (%.1f%%)", pct_correct)) +
  annotate("text", x = 1.29, y = -0.42, hjust = 1, size = 2.9, fontface = "bold",
           color = "#B71C1C", label = sprintf("Errors (%.1f%%)", pct_error)) +
  scale_color_manual(values = COLORS, breaks = LEVELS) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.1), nrow = 1)) +
  coord_cartesian(xlim = c(0.28, 1.3)) +
  labs(
    title = "Posterior predictive check",
    subtitle = "Observed trials and 100 draws per model",
    x = "RT (s)", y = "Density  (up = correct, down = error)", color = NULL
  ) +
  theme_minimal(base_size = 9) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(t = -4, b = 0),
    legend.key.height = unit(0.75, "lines"),
    panel.grid.minor = element_blank()
  )

# ---------------------------------------------------------------------------
# The error-timing statistic, printed rather than plotted.
#
# This used to be panel B. The paper no longer discusses which architecture
# matches this particular statistic - on three pooled participants that is a
# partial result, not a finding - so the figure is down to the two panels that
# stand on their own. The numbers are still computed and printed, because they
# are what tells you whether the predictive check in panel A is hiding a
# systematic error-timing miss.
# ---------------------------------------------------------------------------

observed <- mean(df$RT[df$Error == 1]) - mean(df$RT[df$Error == 0])
cat(sprintf("observed mean(error) - mean(correct) = %+.3f s\n", observed))
for (nm in LEVELS) {
  d <- pred[pred$Model == nm, ]
  s <- tapply(seq_len(nrow(d)), d$draw, function(i) {
    mean(d$rt[i][d$response[i] == 1]) - mean(d$rt[i][d$response[i] == 0])
  })
  v <- as.numeric(s)
  cat(sprintf("%-6s predicted %+.3f [%+.3f, %+.3f]\n", nm, mean(v),
              quantile(v, .025), quantile(v, .975)))
}

# ---------------------------------------------------------------------------
# Panel B: predictive accuracy against sampling cost.
# ---------------------------------------------------------------------------

elpd <- vapply(models, function(m) m$criteria$loo$estimates["elpd_loo", "Estimate"], numeric(1))
best <- names(which.max(elpd))
cmp <- loo::loo_compare(lapply(models, function(m) m$criteria$loo))
diffs <- data.frame(
  Model = LEVELS,
  elpd_diff = elpd[LEVELS] - max(elpd),
  se_diff = vapply(LEVELS, function(nm) {
    if (nm == best) return(0)
    # SE of the pairwise difference against the best model
    a <- models[[best]]$criteria$loo$pointwise[, "elpd_loo"]
    b <- models[[nm]]$criteria$loo$pointwise[, "elpd_loo"]
    sqrt(length(a)) * sd(b - a)
  }, numeric(1)),
  time = vapply(models[LEVELS], function(m) {
    median(attributes(m$fit)$metadata$time$chain$total)
  }, numeric(1))
)
diffs$Model <- factor(diffs$Model, levels = LEVELS)
print(diffs)

pB <- ggplot(diffs, aes(x = time, y = elpd_diff, color = Model)) +
  geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
  geom_linerange(aes(ymin = elpd_diff - se_diff, ymax = elpd_diff + se_diff), linewidth = 0.5) +
  geom_point(size = 2.2) +
  ggrepel::geom_text_repel(aes(label = Model), size = 2.8, fontface = "bold",
                           seed = 1, min.segment.length = 0.4,
                           box.padding = 0.35, show.legend = FALSE) +
  scale_x_log10(breaks = c(10, 30, 100, 300), labels = c("10", "30", "100", "300")) +
  scale_color_manual(values = COLORS, guide = "none") +
  labs(
    title = "What each architecture costs",
    subtitle = "Predictive accuracy against sampling time",
    x = "Median sampling time per chain (s, log scale)",
    y = "ELPD vs. best (±1 SE)"
  ) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank())

# ---------------------------------------------------------------------------

p <- (pA | pB) +
  plot_layout(widths = c(1.15, 1)) +
  plot_annotation(
    tag_levels = "A",
    title = "Five evidence accumulation models, one interface",
    subtitle = sprintf(
      paste("Choices and reaction times from three participants of",
            "Wagenmakers et al. (2008), %s trials"),
      format(n_obs, big.mark = ",")
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 11, colour = "#22252E"),
      plot.subtitle = element_text(face = "italic", size = 8.6, colour = "#6B7280")
    )
  ) &
  theme(plot.tag = element_text(face = "bold", size = 10))

# Sized so that, once scaled to the width it occupies in the two-column PDF,
# the heading prints at ~10.5 pt like the other figures.
ggsave("figures/fig_eam.png", p, width = 7.2, height = 3.9, dpi = 300, bg = "white")
cat("\nwrote figures/fig_eam.png\n")
