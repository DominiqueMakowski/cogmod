#' Single-Accumulator LBA model
#'
#' @description
#' This function simulates reaction times using a single-accumulator
#' version of the Linear Ballistic Accumulator (LBA) model.
#'
#' @param n A single positive integer giving the number of trials.
#' @param drift Mean drift rate.
#' @param sigma Standard deviation of the drift rate.
#' @param sigmabias The starting-point range (A); must be positive.
#' @param bs The additional offset such that b = sigmabias + bs; must be positive.
#' @param ndt Non-decision time.
#' 
#' @examples
#' # Simulate 1000 trials with specified parameters
#' rts <- rrt_lba(n = 1000, drift = 3, sigma = 1, sigmabias = 0.5, bs = 0.5, ndt = 0.3)
#' hist(rts, breaks = 100, xlab = "RT (s)", col = "lightblue")
#' 
#' @export
rrt_lba <- function(n,
                    drift   = 3,
                    sigma   = 1,
                    sigmabias = 0.5,
                    bs      = 0.5,
                    ndt     = 0.3) {

  # — Input checks (as before) —
  if (length(n)!=1 || n<=0 || n!=floor(n)) {
    stop("`n` must be a single positive integer.")
  }
  if (any(c(sigma, sigmabias, bs) <= 0) || ndt < 0) {
    stop("`sigma`, `sigmabias`, `bs` > 0 and `ndt` >= 0.")
  }

  # — Derived threshold —
  b <- sigmabias + bs

  # — Vectorized drift sampling: truncated Normal(mu=drift, sigma=sigma) on (0,Inf) —
  # Uses the internal function .rnorm_truncated() copied from msm::rtnorm (https://github.com/chjackson/msm/blob/master/R/tnorm.R)
  v <- .rnorm_truncated(n, mean = drift, sd = sigma, lower = 0, upper = Inf)

  # — Vectorized uniform start‐point sampling over [0, A] —
  sp <- stats::runif(n, min = 0, max = sigmabias)

  # — Decision times for all trials at once —
  dt <- (b - sp) / v

  # — Final RTs: add non‐decision time —
  ndt + dt
}



#' @rdname rrt_lba
#' @inheritParams lnr
#' @export
drt_lba <- function(x, drift = 1, sigma = 1, sigmabias = 0.5, bs = 0.5, ndt = 0.3, log = FALSE) {
  out <- numeric(length(x))
  T <- x - ndt
  ok <- T > 0
  
  if (!any(ok)) {
    return(if (log) rep(-Inf, length(x)) else out)
  }
  
  A <- sigmabias
  b <- sigmabias + bs
  Tp <- T[ok]
  st <- sigma * Tp
  
  # Compute the defective density using the formula from Brown & Heathcote (2008)
  z1 <- (b - A - drift * Tp) / st
  z2 <- (b - drift * Tp) / st
  
  fval <- (1 / A) * (
    drift * (stats::pnorm(z2) - stats::pnorm(z1)) +
    sigma * (stats::dnorm(z1) - stats::dnorm(z2))
  )
  
  # Normalize by the probability that drift rate is positive
  prob_positive_drift <- 1 - stats::pnorm(0, mean = drift, sd = sigma)
  fval <- fval / prob_positive_drift
  
  fval <- pmax(fval, 1e-10)  # Avoid log(0)
  out[ok] <- if (log) log(fval) else fval
  out
}