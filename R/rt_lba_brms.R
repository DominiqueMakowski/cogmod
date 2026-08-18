#' @rdname rrt_lba
#' @param link_mu,link_sigma,link_sigmabias,link_bs,link_ndt,link_poutlier Link
#'   functions for the parameters.
#' @param predict_outliers Logical; whether `posterior_predict()` should include
#'   the outlier component. `FALSE` (the default) fixes `poutlier` to zero for
#'   prediction, so predictions describe the decision process alone; the
#'   likelihood is always the full mixture either way. See [with_outliers()].
#' @inheritParams rrt_lognormal
#' @export
rt_lba <- function(
  link_mu = "softplus",
  link_sigma = "softplus",
  link_sigmabias = "softplus",
  link_bs = "softplus",
  link_ndt = "log",
  link_poutlier = "logit",
  predict_outliers = FALSE,
  minrt = 0.3
) {
  .rt_shifted_family(
    "rt_lba",
    links = c(link_mu, link_sigma, link_sigmabias, link_bs),
    predict_outliers = predict_outliers,
    minrt = minrt
  )
}


#' @keywords internal
.rt_lba_lpdf <- function(minrt = 0.3) {
  .rt_shifted_lpdf("rt_lba", minrt = minrt)
}


#' @rdname rrt_lba
#' @export
rt_lba_lpdf_expose <- function(minrt = 0.3) {
  .rt_shifted_expose("rt_lba", minrt)
}


#' @rdname rrt_lba
#' @export
rt_lba_stanvars <- function(minrt = 0.3) {
  brms::stanvar(scode = .rt_lba_lpdf(.as_minrt(minrt)), block = "functions")
}


# brms Post-processing Functions ------------------------------------------

#' @rdname rrt_lba
#' @inheritParams lnr
#' @export
log_lik_rt_lba <- function(i, prep) {
  .log_lik_shifted("rt_lba", i, prep)
}


#' @rdname rrt_lba
#' @export
posterior_predict_rt_lba <- function(i, prep, predict_outliers = NULL, ...) {
  .posterior_predict_shifted("rt_lba", i, prep, predict_outliers)
}


#' @rdname rrt_lba
#' @export
posterior_epred_rt_lba <- function(prep, predict_outliers = NULL) {
  # Not merely expensive: the expectation does not exist. The decision time is
  # (b - U(0, A)) / v with v a normal truncated at zero, and since that density
  # is positive at v = 0, E[1/v] diverges. Simulating and summarising is the
  # only option, and even the sample mean will not settle down.
  stop(
    "The single-accumulator LBA decision time has no finite mean - E[1/drift] ",
    "diverges because the truncated-normal drift has positive density at 0 - ",
    "so there is no expectation to return.\n",
    "Use `posterior_predict()` and summarise the draws with a median or ",
    "another quantile instead.",
    call. = FALSE
  )
}
