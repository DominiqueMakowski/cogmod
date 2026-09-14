# Figure for warmstart_ablation.R: the effective draws per second of every
# run, per target data set and family, drawn from the results CSV that script
# writes. warmstart_ablation.R sources this once its runs are done; it can
# also be run on its own to redraw the figure from the CSV without refitting
# anything (from the repository root):
#
#   Rscript benchmarks/warmstart_ablation_figure.R
#
# The figure is written to man/figures/, which is where the articles read
# their static images from, so that vignettes/articles/performance.qmd can
# include it with
#
#   ![](../man/figures/warmstart_ablation.png)
#
# What the figure shows. One row of panels per family, one column per target
# data set (the pilot's participants included; other participants of the same
# task; other participants with shifted RTs - see the design notes at the top
# of warmstart_ablation.R). Each panel: the smallest bulk ESS over a run's
# parameters per second of its wall time (slowest chain, warmup and sampling),
# which is the number that matters. Each family row has its own y axis - the
# families differ by an order of magnitude in cost, and the point is the
# comparison of runs within a family, across the three targets. Bars are the
# mean over seeds, points the seeds. A cell where a seed failed to converge is
# faded and labelled with its Rhat: a stuck chain gives an ESS/s near zero,
# which is a failure, not a slow run. The colours group the runs by what they
# are: grey for the two cold runs everything is measured against, blues for
# the runs that change only the sampler (the more of the warm start they use,
# the darker), and two warm colours for the runs that also change the model by
# re-centring the priors - red for the full warm start plus the pilot's
# priors, orange for priors built from the target's own reference fit. Those
# two sample a different posterior from every other bar, so their ESS/s is not
# a like-for-like speedup, and the eye should not read them as the best
# sampler setting. Whether the different posterior is a legitimate one depends
# on the column: the caption says how.

suppressPackageStartupMessages(library(ggplot2))

# Drawing order: cold runs, then the sampler-only cells in order of how much
# of the warm start they use, then the two cells that change the model.
ABLATION_RUNS <- c("reference", "base", "step_size", "inv_metric", "both",
                   "both + inits", "all + pilot priors", "oracle priors")
ABLATION_COLOURS <- c(
  reference = "grey65", base = "grey40",
  step_size = "#90CAF9", inv_metric = "#42A5F5", both = "#1E88E5", "both + inits" = "#0D47A1",
  "all + pilot priors" = "#C62828", "oracle priors" = "#FB8C00"
)
ABLATION_LABELS <- c(
  reference = "reference (cold, full warmup)",
  base = "base (cold, short warmup)",
  step_size = "base + step_size",
  inv_metric = "base + inv_metric",
  both = "base + inv_metric + step_size",
  "both + inits" = "base + inv_metric + step_size + inits (full warm start)",
  "all + pilot priors" = "full warm start + priors from the pilot",
  "oracle priors" = "base + priors from the target's own reference fit (oracle)"
)
# The target data sets, as the columns of the figure.
ABLATION_DATASETS <- c(
  pilot_included = "Target includes the pilot's participants",
  new_participants = "Other participants, same task",
  shifted = "Other participants, RTs shifted by 0.3 s (ndt doubled)"
)

# `results`: the rows of warmstart_ablation.csv (one per fit). The pilot rows
# are left out - the pilot is fitted to fewer participants, so its speed is
# not on the same footing as the runs. `families` gives the row order and
# `datasets` the column order; anything not named comes after, in order of
# appearance, and a run or data set the labels do not know is drawn under its
# own name (a run in light grey).
warmstart_ablation_figure <- function(results, families = NULL, datasets = NULL) {
  if (is.null(results$dataset)) results$dataset <- "pilot_included"
  full <- results[results$run != "pilot" & is.finite(results$ess_per_s), , drop = FALSE]
  if (!nrow(full)) stop("No full-data runs in the results; nothing to draw.", call. = FALSE)

  ordered_levels <- function(x, preferred) c(intersect(preferred, unique(x)), setdiff(unique(x), preferred))
  fam_levels <- ordered_levels(full$family, families)
  ds_levels <- ordered_levels(full$dataset, if (is.null(datasets)) names(ABLATION_DATASETS) else datasets)
  run_levels <- ordered_levels(full$run, ABLATION_RUNS)
  full$family <- factor(toupper(full$family), levels = toupper(fam_levels))
  full$dataset <- factor(full$dataset, levels = ds_levels)
  full$run <- factor(full$run, levels = run_levels)
  colours <- ABLATION_COLOURS[run_levels]
  colours[is.na(colours)] <- "grey80"
  names(colours) <- run_levels

  means <- aggregate(ess_per_s ~ dataset + family + run, full, mean)
  # the label sits above whatever is tallest in the bar's column: the mean, or
  # a seed's point when one lands above it
  tops <- aggregate(ess_per_s ~ dataset + family + run, full, max)
  # A cell where any seed failed to converge is flagged with its worst Rhat:
  # an ESS/s of nearly zero from a stuck chain is a failure, not a slow run,
  # and the two must not be read from the same axis without being told apart.
  worst <- aggregate(rhat_max ~ dataset + family + run, full, max)
  tops <- merge(tops, worst, by = c("dataset", "family", "run"))
  # 1.1 rather than the textbook 1.05: with 1200 draws a healthy fit of these
  # models sits at 1.01-1.07, and a chain that is actually stuck shows 1.3-3.
  tops$failed <- tops$rhat_max > 1.1
  tops$label <- ifelse(tops$failed,
                       sprintf("%s\nRhat %.2f", signif(tops$ess_per_s, 2), tops$rhat_max),
                       as.character(signif(tops$ess_per_s, 2)))
  means$failed <- tops$failed[match(paste(means$dataset, means$family, means$run),
                                    paste(tops$dataset, tops$family, tops$run))]
  n_seeds <- length(unique(full$seed))

  run_label <- function(x) ifelse(x %in% names(ABLATION_LABELS), ABLATION_LABELS[x], x)
  ds_label <- function(x) ifelse(x %in% names(ABLATION_DATASETS), ABLATION_DATASETS[x], x)
  first <- full[1, ]
  subtitle <- sprintf(
    paste0("Smallest bulk ESS over all parameters, per second of wall time ",
           "(slowest chain, warmup and sampling); higher is better.\n",
           "One pilot per family, fitted to 5 participants; each column a target of ",
           "%d participants, %d trials; %d chains with %d sampling iterations each; ",
           "warmup %d for the reference, %d for every other run.\n",
           "%d seed%s per cell: bars are their mean%s."),
    first$participants, first$n_obs, first$chains, first$sample_iter,
    max(full$warmup), min(full$warmup),
    n_seeds, if (n_seeds == 1) "" else "s", if (n_seeds > 1) ", points the seeds" else ""
  )

  p <- ggplot(means, aes(x = run, y = ess_per_s, fill = run)) +
    geom_col(aes(alpha = failed), width = 0.75) +
    scale_alpha_manual(values = c(`FALSE` = 1, `TRUE` = 0.35), guide = "none") +
    geom_text(data = tops, aes(label = label, colour = failed), vjust = -0.3, size = 2.7,
              lineheight = 0.85) +
    scale_colour_manual(values = c(`FALSE` = "grey20", `TRUE` = "#B71C1C"), guide = "none") +
    facet_grid(family ~ dataset, scales = "free_y", switch = "y",
               labeller = labeller(dataset = ds_label)) +
    scale_fill_manual(values = colours, labels = run_label,
                      guide = guide_legend(nrow = 2, byrow = TRUE)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(
      x = NULL, y = NULL, fill = NULL,
      title = "What each half of a warm start buys, and what happens when the source no longer fits",
      subtitle = subtitle,
      caption = paste0(
        "Blue runs sample the same posterior as the grey ones, only with a different metric, ",
        "step size or starting values; the two priors runs also change the model.\n",
        "Orange: priors built from a full fit of the very same target data, so the prior already ",
        "knows the answer - an upper bound on what informative priors could buy, not something ",
        "one can do.\n",
        "Red: the full warm start plus priors from the pilot. That double-counts the pilot's ",
        "participants in the first column, is the legitimate use of an independent sample in the ",
        "second,\nand a confidently wrong prior in the third - compare the answers, not only the speed. ",
        "A faded bar labelled with an Rhat did not converge on at least one seed (a stuck chain)."
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      strip.placement = "outside",
      strip.text.x = element_text(size = 10),
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 11),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.spacing.x = unit(1, "lines"),
      panel.spacing.y = unit(1, "lines"),
      legend.position = "bottom",
      legend.key.size = unit(0.9, "lines"),
      legend.text = element_text(size = 9),
      plot.subtitle = element_text(size = 9.5, colour = "grey30"),
      plot.caption = element_text(hjust = 0, colour = "grey30"),
      plot.title.position = "plot",
      plot.caption.position = "plot"
    )
  if (n_seeds > 1) {
    p <- p + geom_point(data = full, shape = 21, fill = "white", size = 1.5, stroke = 0.5,
                        position = position_jitter(width = 0.12, height = 0, seed = 1))
  }
  p
}

# Where the figure goes: where the articles look for it.
ablation_figure_file <- function() file.path("man", "figures", "warmstart_ablation.png")

save_warmstart_ablation_figure <- function(results, files = ablation_figure_file(),
                                           families = NULL, datasets = NULL,
                                           width = 12, height = NULL, dpi = 150) {
  p <- warmstart_ablation_figure(results, families, datasets)
  if (is.null(height)) {
    n_rows <- length(unique(results$family[results$run != "pilot"]))
    height <- 3 + 2.3 * max(n_rows, 1)
  }
  for (f in files) {
    dir.create(dirname(f), showWarnings = FALSE, recursive = TRUE)
    ggsave(f, p, width = width, height = height, dpi = dpi, bg = "white")
  }
  invisible(p)
}

# Run directly (not sourced): redraw from the CSV on disk.
if (sys.nframe() == 0L) {
  results <- read.csv(file.path("benchmarks", "results", "warmstart_ablation.csv"),
                      stringsAsFactors = FALSE)
  save_warmstart_ablation_figure(results)
  cat("Figure written to:", ablation_figure_file(), "\n")
}
