# Experiment 1 of Wagenmakers, Ratcliff, Gomez & McKoon (2008), as distributed
# in the `rtdists` package (courtesy of the original authors). The data are not
# redistributed with `cogmod`; this script reconstructs the subset used in the
# paper and in the vignettes.
#
# `censor` flags the trials the original authors excluded: uninterpretable
# responses, RT < 180 ms and RT > 3 s. Following them, we additionally cap RT at
# 2 s.
#
# Sourced by make_figures.R and make_fig_eam.R.

wagenmakers <- function() {
  if (!requireNamespace("rtdists", quietly = TRUE)) {
    stop("The 'rtdists' package is required to reconstruct these data.")
  }
  e <- new.env()
  utils::data("speed_acc", package = "rtdists", envir = e)
  speed_acc <- e$speed_acc

  d <- speed_acc[
    !speed_acc$censor & speed_acc$response != "error" & speed_acc$rt <= 2,
  ]

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
