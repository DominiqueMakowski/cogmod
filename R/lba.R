#' @title Linear Ballistic Accumulator (LBA) Model Simulation
#'
#' @description
#' Simulates random draws (reaction times and choices) from a two-choice Linear Ballistic Accumulator (LBA) model.
#'
#' In this parametrization, each accumulator has its own independent drift rate distribution:
#' - Accumulator 0 has drift drawn from N(driftzero, sigmazero^2).
#' - Accumulator 1 has drift drawn from N(driftone, sigmaone^2).
#'
#' For each trial, drift rates are sampled on an individual basis until at least one of the two is positive.
#' The starting point for each accumulator is sampled uniformly from (0, sigmabias). The decision threshold is defined as sigmabias + bs.
#' The decision time for an accumulator is calculated as (b - start)/drift, and if its drift is not positive, its decision time is set to Inf.
#' The winning accumulator (the one whose decision time is minimal) determines the response, and the final reaction time is the sum of its decision time
#' and a fixed non-decision time (ndt).
#'
#' @param n Number of simulated trials. Must be a positive integer.
#' @param driftzero Mean drift rate for the first accumulator (accumulator 0). Range: (-Inf, Inf).
#' @param driftone Mean drift rate for the second accumulator (accumulator 1). Range: (-Inf, Inf).
#' @param sigmazero Standard deviation of the drift rate for the first accumulator. Must be positive. Range: (0, Inf).
#' @param sigmaone Standard deviation of the drift rate for the second accumulator. Must be positive. Range: (0, Inf).
#' @param sigmabias Maximum starting point for the uniform distribution of starting evidence (0, sigmabias). Must be positive.
#'   Range: (0, Inf). Default: 0.5.
#' @param bs Additional amount beyond `sigmabias` to set the decision threshold (b = sigmabias + bs). Must be positive.
#'   Range: (0, Inf). Default: 0.5.
#' @param ndt Non-decision time, representing processes such as encoding and motor response. Must be non-negative.
#'   Range: [0, Inf). Default: 0.3.
#' @param max_iter Maximum iterations allowed (per trial) for resampling drift rates if both
#'   are non-positive. Default: 100.
#'
#' @details
#' **Psychological Interpretation:**
#' - **Drift Rate (`driftzero`, `driftone`)**: Reflects the rate at which evidence accumulates for each choice. Higher drift rates indicate faster
#'   evidence accumulation and a higher likelihood of selecting the corresponding choice. Differences in drift rates between accumulators
#'   can represent differences in preference, difficulty, or bias between the two options.
#' - **Drift Rate Variability (`sigmazero`, `sigmaone`)**: Captures trial-to-trial variability in the evidence accumulation process.
#'   Higher variability indicates less consistent evidence accumulation, leading to greater variability in reaction times and choices.
#' - **Start Point Variability (`sigmabias`)**: Represents the range of initial evidence levels for each accumulator. Larger values of `sigmabias` introduce
#'   more variability in reaction times, as the starting point can vary more widely between trials.
#' - **Threshold (`b = sigmabias + bs`)**: Boundary separation (`bs`). Represents the amount of evidence required to make a decision. Higher thresholds lead to longer reaction times
#'   but more accurate decisions, as more evidence is required before a choice is made.
#' - **Non-Decision Time (`ndt`)**: Accounts for processes unrelated to evidence accumulation, such as sensory encoding and motor response.
#'   This parameter shifts all reaction times by a constant amount.
#'
#' @references
#' - Brown, S. D., & Heathcote, A. (2008). The simplest complete model of choice response time: Linear ballistic accumulation.
#'     *Cognitive Psychology*, *57*(3), 153-178. \doi{10.1016/j.cogpsych.2007.12.002}
#'
#' @examples
#' df <- rlba(n = 1000, driftzero = 3, driftone = 2,
#'            sigmazero = 0.5, sigmaone = 0.5,
#'            sigmabias = 0.5, bs = 0.5, ndt = 0.3)
#' hist(df$rt[df$response == 0], breaks = 50, col = rgb(0,0,1,0.5))
#' hist(df$rt[df$response == 1], breaks = 50, col = rgb(1,0,0,0.5), add = TRUE)
#'
#' @export
rlba <- function(
  n,
  driftzero = 3,
  driftone = 3,
  sigmazero = 1,
  sigmaone = 1,
  sigmabias = 0.5,
  bs = 0.5,
  ndt = 0.3,
  max_iter = 100
) {
  # --- Input Validation ---
  if (length(n) != 1 || n <= 0 || n != floor(n)) {
    stop("n must be a single positive integer.")
  }

  # Recycle every parameter to length n *before* validating or drawing.
  # `posterior_predict_lba()` calls this with one value per posterior draw, so
  # the parameters have to be honoured trial by trial: scalar `rnorm(1, mean =
  # driftzero, ...)` inside a loop would silently reuse `driftzero[1]` for
  # every trial, and scalar `||` on a vector is an error from R 4.3 onwards.
  driftzero <- rep_len(driftzero, n)
  driftone <- rep_len(driftone, n)
  sigmazero <- rep_len(sigmazero, n)
  sigmaone <- rep_len(sigmaone, n)
  sigmabias <- rep_len(sigmabias, n)
  bs <- rep_len(bs, n)
  ndt <- rep_len(ndt, n)

  if (
    any(sigmabias <= 0) ||
      any(bs <= 0) ||
      any(sigmazero <= 0) ||
      any(sigmaone <= 0) ||
      any(ndt < 0)
  ) {
    stop(
      "sigmabias, bs, sigmazero, sigmaone must be positive; ndt must be non-negative."
    )
  }

  # --- Derived Parameter ---
  b <- sigmabias + bs # Decision threshold

  # Sample drift rates, resampling only those trials where neither accumulator
  # has a positive drift (such a trial would never produce a response).
  v0 <- stats::rnorm(n, mean = driftzero, sd = sigmazero)
  v1 <- stats::rnorm(n, mean = driftone, sd = sigmaone)

  bad <- v0 <= 0 & v1 <= 0
  iter <- 0
  while (any(bad) && iter < max_iter) {
    iter <- iter + 1
    nb <- sum(bad)
    v0[bad] <- stats::rnorm(nb, mean = driftzero[bad], sd = sigmazero[bad])
    v1[bad] <- stats::rnorm(nb, mean = driftone[bad], sd = sigmaone[bad])
    bad <- v0 <= 0 & v1 <= 0
  }
  if (any(bad)) {
    warning(sprintf(
      "%d trial(s) reached max_iter; forcing accumulator 0 positive.",
      sum(bad)
    ))
    v0[bad] <- abs(stats::rnorm(
      sum(bad),
      mean = driftzero[bad],
      sd = sigmazero[bad]
    ))
  }

  # Starting points from U(0, sigmabias)
  start0 <- stats::runif(n, min = 0, max = sigmabias)
  start1 <- stats::runif(n, min = 0, max = sigmabias)

  # Decision times, treating a non-positive drift as never finishing.
  time0 <- ifelse(v0 > 0, (b - start0) / v0, Inf)
  time1 <- ifelse(v1 > 0, (b - start1) / v1, Inf)

  # The accumulator that reaches threshold first wins (ties go to 1, as before).
  data.frame(
    rt = ndt + pmin(time0, time1),
    response = ifelse(time0 < time1, 0L, 1L)
  )
}


#' The density function `dlba` calculates the likelihood of observing a specific
#' reaction time `rt` and response `response`, given the LBA parameters. It is
#' based on the formulation by Brown & Heathcote (2008), where the likelihood
#' is the product of the probability density of the winning accumulator finishing
#' at time `t = rt - ndt` and the probability (survival function) that the losing
#' accumulator has not finished by time `t`. This implementation assumes that
#' the `rlba` function ensures at least one positive drift per trial, so no
#' additional normalization by `(1 - pnegative)` is required.
#' @rdname rlba
#' @inheritParams rlnr
#' @export
dlba <- function(
  x,
  driftzero = 3,
  driftone = 3,
  sigmazero = 1,
  sigmaone = 1,
  sigmabias = 0.5,
  bs = 0.5,
  ndt = 0.3,
  response,
  log = FALSE
) {
  # Recycle every argument to a common length. This is essential because `dlba()`
  # is called under two different conventions: directly, with a vector of
  # observations and scalar parameters; and from `log_lik_lba()`, with a *single*
  # observation and one parameter value per posterior draw. Indexing `x` with a
  # logical vector longer than itself (as an earlier version did) silently
  # produces NAs, which then get floored to `.Machine$double.eps` - yielding a
  # constant log-density of about -36 for every draw.
  n_out <- max(
    length(x),
    length(response),
    length(driftzero),
    length(driftone),
    length(sigmazero),
    length(sigmaone),
    length(sigmabias),
    length(bs),
    length(ndt)
  )
  x <- rep(x, length.out = n_out)
  response <- rep(response, length.out = n_out)
  driftzero <- rep(driftzero, length.out = n_out)
  driftone <- rep(driftone, length.out = n_out)
  sigmazero <- rep(sigmazero, length.out = n_out)
  sigmaone <- rep(sigmaone, length.out = n_out)
  sigmabias <- rep(sigmabias, length.out = n_out)
  bs <- rep(bs, length.out = n_out)
  ndt <- rep(ndt, length.out = n_out)

  A <- sigmabias
  b <- sigmabias + bs
  dt <- x - ndt

  dens <- rep(.Machine$double.eps, n_out)

  valid <- !is.na(dt) & dt > 0
  idx0 <- which(valid & response == 0)
  if (length(idx0) > 0) {
    f0 <- .lba_defectivedensity(
      dt[idx0], driftzero[idx0], sigmazero[idx0], A[idx0], b[idx0]
    )
    F1 <- .lba_cumulative(
      dt[idx0], driftone[idx0], sigmaone[idx0], A[idx0], b[idx0]
    )
    dens[idx0] <- f0 * (1 - F1)
  }
  idx1 <- which(valid & response == 1)
  if (length(idx1) > 0) {
    f1 <- .lba_defectivedensity(
      dt[idx1], driftone[idx1], sigmaone[idx1], A[idx1], b[idx1]
    )
    F0 <- .lba_cumulative(
      dt[idx1], driftzero[idx1], sigmazero[idx1], A[idx1], b[idx1]
    )
    dens[idx1] <- f1 * (1 - F0)
  }

  # Guard against NA/NaN/Inf arising from extreme parameter draws (e.g. during
  # warmup or in poorly-mixing chains) so downstream code (e.g. loo) doesn't crash.
  dens[!is.finite(dens)] <- .Machine$double.eps
  dens[dens < .Machine$double.eps] <- .Machine$double.eps

  out <- if (log) log(dens) else dens

  # Return -Inf (or a very small positive number) for RT at or below ndt,
  # matching the `if (Y <= ndt) return negative_infinity();` guard in the Stan
  # implementation (see `.lba_lpdf()`).
  below_ndt <- !is.na(x) & !is.na(ndt) & x <= ndt
  out[below_ndt] <- if (log) -Inf else .Machine$double.eps

  out
}


# Internal helper functions --------------------------------------------------

# Helper function: defective density for an accumulator
#' @keywords internal
.lba_defectivedensity <- function(dt, v, s, A, b) {
  # dt: decision time(s) (must be positive)
  # v: mean drift for the accumulator
  # s: standard deviation for the accumulator's drift
  # A: start point range (sigmabias)
  # b: decision threshold (sigmabias + bs)
  n1 <- (b - A - v * dt) / (dt * s)
  n2 <- (b - v * dt) / (dt * s)
  f_val <- (1 / A) *
    (-v *
      stats::pnorm(n1) +
      s * stats::dnorm(n1) +
      v * stats::pnorm(n2) -
      s * stats::dnorm(n2))
  f_val
}

# Helper function: cumulative density function for an accumulator
#' @keywords internal
.lba_cumulative <- function(dt, v, s, A, b) {
  n1 <- (b - A - v * dt) / (dt * s)
  n2 <- (b - v * dt) / (dt * s)
  F_val <- 1 +
    ((b - A - v * dt) / A) * stats::pnorm(n1) -
    ((b - v * dt) / A) * stats::pnorm(n2) +
    ((dt * s) / A) * (stats::dnorm(n1) - stats::dnorm(n2))
  F_val
}
