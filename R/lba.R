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
rlba <- function(n,
                 driftzero = 3,
                 driftone = 3,
                 sigmazero = 1,
                 sigmaone = 1,
                 sigmabias = 0.5,
                 bs = 0.5,
                 ndt = 0.3,
                 max_iter = 100) {

  # --- Input Validation ---
  if (length(n) != 1 || n <= 0 || n != floor(n))
    stop("n must be a single positive integer.")
  if (sigmabias <= 0 || bs <= 0 || sigmazero <= 0 || sigmaone <= 0 || ndt < 0)
    stop("sigmabias, bs, sigmazero, sigmaone must be positive; ndt must be non-negative.")

  # --- Derived Parameter ---
  b <- sigmabias + bs  # Decision threshold

  # --- Prepare Output ---
  rates0 <- numeric(n)
  rates1 <- numeric(n)
  choices <- integer(n)
  rts <- numeric(n)

  for (i in seq_len(n)) {
    # Sample both drift rates until at least one is positive.
    iter <- 0
    repeat {
      iter <- iter + 1
      v0 <- stats::rnorm(1, mean = driftzero, sd = sigmazero)
      v1 <- stats::rnorm(1, mean = driftone, sd = sigmaone)
      if (v0 > 0 || v1 > 0) break
      if (iter >= max_iter) {
        warning(sprintf("Trial %d reached max_iter; forcing accumulator 0 positive.", i))
        v0 <- abs(stats::rnorm(1, mean = driftzero, sd = sigmazero))
        break
      }
    }
    rates0[i] <- v0
    rates1[i] <- v1

    # Compute starting points from U(0, sigmabias)
    start0 <- stats::runif(1, min = 0, max = sigmabias)
    start1 <- stats::runif(1, min = 0, max = sigmabias)
    
    # Compute decision times using only positive drifts, treat non-positive as Inf.
    time0 <- if (v0 > 0) { (b - start0) / v0 } else { Inf }
    time1 <- if (v1 > 0) { (b - start1) / v1 } else { Inf }
    
    # Choose the accumulator that reached threshold first.
    if (time0 < time1) {
      choices[i] <- 0
      rts[i] <- ndt + time0
    } else {
      choices[i] <- 1
      rts[i] <- ndt + time1
    }
  }
  
  # Return Results
  data.frame(rt = rts, response = choices)
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
dlba <- function(x, driftzero = 3, driftone = 3, sigmazero = 1, sigmaone = 1,
                 sigmabias = 0.5, bs = 0.5, ndt = 0.3, response, log = FALSE) {
  # Return -Inf (or a very small positive number) for RT below ndt.
  below_ndt <- x < ndt
  out <- if (log) rep(-Inf, length(x)) else rep(.Machine$double.eps, length(x))
  
  keep <- !below_ndt
  if (any(keep)) {
    A <- sigmabias
    b <- sigmabias + bs
    dt <- x[keep] - ndt
    dens <- rep(.Machine$double.eps, length(dt))
    
    valid <- dt > 0
    if (any(valid)) {
      idx0 <- valid & (response[keep] == 0)
      if (any(idx0)) {
        f0 <- .lba_defectivedensity(dt[idx0], driftzero, sigmazero, A, b)
        F1 <- .lba_cumulative(dt[idx0], driftone, sigmaone, A, b)
        dens[idx0] <- f0 * (1 - F1)
      }
      idx1 <- valid & (response[keep] == 1)
      if (any(idx1)) {
        f1 <- .lba_defectivedensity(dt[idx1], driftone, sigmaone, A, b)
        F0 <- .lba_cumulative(dt[idx1], driftzero, sigmazero, A, b)
        dens[idx1] <- f1 * (1 - F0)
      }
    }

    dens[dens < .Machine$double.eps] <- .Machine$double.eps
    
    if (log)
      out[keep] <- log(dens)
    else
      out[keep] <- dens
  }
  
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
  f_val <- (1/A) * (-v * stats::pnorm(n1) + s * stats::dnorm(n1) +
                     v * stats::pnorm(n2) - s * stats::dnorm(n2))
  f_val
}

# Helper function: cumulative density function for an accumulator
#' @keywords internal
.lba_cumulative <- function(dt, v, s, A, b) {
  n1 <- (b - A - v * dt) / (dt * s)
  n2 <- (b - v * dt) / (dt * s)
  F_val <- 1 + ((b - A - v * dt) / A) * stats::pnorm(n1) -
    ((b - v * dt) / A) * stats::pnorm(n2) +
    ((dt * s) / A) * (stats::dnorm(n1) - stats::dnorm(n2))
  F_val
}