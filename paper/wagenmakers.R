# Experiment 1 of Wagenmakers, Ratcliff, Gomez & McKoon (2008), as distributed
# in the `rtdists` package (courtesy of the original authors). The data are not
# redistributed with `cogmod`; this script reconstructs the subset used in the
# paper and in the vignettes.
#
# The only exclusion is `RT <= 2` s. This deliberately matches the vignettes
# exactly - the `censor` column and the "error" response level are *not*
# filtered out here, because the models the paper reports were fitted on the
# unfiltered subset. Applying the original authors' censoring instead drops a
# further five trials and puts every n in the paper three or four out from the
# fitted models, which is the sort of drift this script exists to prevent.
#
# Sourced by make_figures.R (which writes wagenmakers2008.csv for
# make_fig_approaches.py).

wagenmakers <- function() {
  if (!requireNamespace("rtdists", quietly = TRUE)) {
    stop("The 'rtdists' package is required to reconstruct these data.")
  }
  e <- new.env()
  utils::data("speed_acc", package = "rtdists", envir = e)
  speed_acc <- e$speed_acc

  d <- speed_acc[speed_acc$rt <= 2, ]

  data.frame(
    Participant = as.integer(as.character(d$id)),
    Condition = unname(c(accuracy = "Accuracy", speed = "Speed")[
      as.character(d$condition)
    ]),
    RT = d$rt,
    Error = as.integer(as.character(d$response) != as.character(d$stim_cat)),
    Frequency = unname(c(high = "High", low = "Low", very_low = "Very Low")[
      sub("^nw_", "", as.character(d$frequency))
    ])
  )
}
