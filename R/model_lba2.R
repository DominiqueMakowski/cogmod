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
#' The starting point for each accumulator is sampled uniformly from (0, sigmabias). The decision threshold is defined as sigmabias + boundary.
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
#' @param boundary Additional amount beyond `sigmabias` to set the decision threshold (b = sigmabias + boundary). Must be positive.
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
#' - **Threshold (`b = sigmabias + boundary`)**: Boundary separation (`boundary`). Represents the amount of evidence required to make a decision. Higher thresholds lead to longer reaction times
#'   but more accurate decisions, as more evidence is required before a choice is made.
#' - **Non-Decision Time (`ndt`)**: Accounts for processes unrelated to evidence accumulation, such as sensory encoding and motor response.
#'   This parameter shifts all reaction times by a constant amount.
#'
#' @references
#' - Brown, S. D., & Heathcote, A. (2008). The simplest complete model of choice response time: Linear ballistic accumulation.
#'     *Cognitive Psychology*, *57*(3), 153-178. \doi{10.1016/j.cogpsych.2007.12.002}
#'
#' @examples
#' df <- rcogmod_lba2(n = 1000, driftzero = 3, driftone = 2,
#'            sigmazero = 0.5, sigmaone = 0.5,
#'            sigmabias = 0.5, boundary = 0.5, ndt = 0.3)
#' hist(df$rt[df$response == 0], breaks = 50, col = rgb(0,0,1,0.5))
#' hist(df$rt[df$response == 1], breaks = 50, col = rgb(1,0,0,0.5), add = TRUE)
#'
#' @export
rcogmod_lba2 <- function(
  n,
  driftzero = 3,
  driftone = 3,
  sigmazero = 1,
  sigmaone = 1,
  sigmabias = 0.5,
  boundary = 0.5,
  ndt = 0.3,
  max_iter = 100
) {
  # --- Input Validation ---
  if (length(n) != 1 || n <= 0 || n != floor(n)) {
    stop("n must be a single positive integer.")
  }

  # Recycle every parameter to length n *before* validating or drawing.
  # `posterior_predict_cogmod_lba2()` calls this with one value per posterior draw, so
  # the parameters have to be honoured trial by trial: scalar `rnorm(1, mean =
  # driftzero, ...)` inside a loop would silently reuse `driftzero[1]` for
  # every trial, and scalar `||` on a vector is an error from R 4.3 onwards.
  driftzero <- rep_len(driftzero, n)
  driftone <- rep_len(driftone, n)
  sigmazero <- rep_len(sigmazero, n)
  sigmaone <- rep_len(sigmaone, n)
  sigmabias <- rep_len(sigmabias, n)
  boundary <- rep_len(boundary, n)
  ndt <- rep_len(ndt, n)

  if (
    any(sigmabias <= 0) ||
      any(boundary <= 0) ||
      any(sigmazero <= 0) ||
      any(sigmaone <= 0) ||
      any(ndt < 0)
  ) {
    stop(
      "sigmabias, boundary, sigmazero, sigmaone must be positive; ndt must be non-negative."
    )
  }

  # --- Derived Parameter ---
  b <- sigmabias + boundary # Decision threshold

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


#' The density function `dcogmod_lba2` calculates the likelihood of observing a specific
#' reaction time `rt` and response `response`, given the LBA parameters. It is
#' based on the formulation by Brown & Heathcote (2008), where the likelihood
#' is the product of the probability density of the winning accumulator finishing
#' at time `t = rt - ndt` and the probability (survival function) that the losing
#' accumulator has not finished by time `t`. This implementation assumes that
#' the `rcogmod_lba2` function ensures at least one positive drift per trial, so no
#' additional normalization by `(1 - pnegative)` is required.
#' @rdname rcogmod_lba2
#' @inheritParams rcogmod_lnr
#' @export
dcogmod_lba2 <- function(
  x,
  driftzero = 3,
  driftone = 3,
  sigmazero = 1,
  sigmaone = 1,
  sigmabias = 0.5,
  boundary = 0.5,
  ndt = 0.3,
  response,
  log = FALSE
) {
  # Recycle every argument to a common length. This is essential because `dcogmod_lba2()`
  # is called under two different conventions: directly, with a vector of
  # observations and scalar parameters; and from `log_lik_cogmod_lba2()`, with a *single*
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
    length(boundary),
    length(ndt)
  )
  x <- rep(x, length.out = n_out)
  response <- rep(response, length.out = n_out)
  driftzero <- rep(driftzero, length.out = n_out)
  driftone <- rep(driftone, length.out = n_out)
  sigmazero <- rep(sigmazero, length.out = n_out)
  sigmaone <- rep(sigmaone, length.out = n_out)
  sigmabias <- rep(sigmabias, length.out = n_out)
  boundary <- rep(boundary, length.out = n_out)
  ndt <- rep(ndt, length.out = n_out)

  A <- sigmabias
  b <- sigmabias + boundary
  dt <- x - ndt

  dens <- rep(.Machine$double.eps, n_out)

  valid <- !is.na(dt) & dt > 0
  idx0 <- which(valid & response == 0)
  if (length(idx0) > 0) {
    f0 <- .cogmod_lba2_defectivedensity(
      dt[idx0], driftzero[idx0], sigmazero[idx0], A[idx0], b[idx0]
    )
    F1 <- .cogmod_lba2_cumulative(
      dt[idx0], driftone[idx0], sigmaone[idx0], A[idx0], b[idx0]
    )
    dens[idx0] <- f0 * (1 - F1)
  }
  idx1 <- which(valid & response == 1)
  if (length(idx1) > 0) {
    f1 <- .cogmod_lba2_defectivedensity(
      dt[idx1], driftone[idx1], sigmaone[idx1], A[idx1], b[idx1]
    )
    F0 <- .cogmod_lba2_cumulative(
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
  # implementation (see `.cogmod_lba2_lpdf()`).
  below_ndt <- !is.na(x) & !is.na(ndt) & x <= ndt
  out[below_ndt] <- if (log) -Inf else .Machine$double.eps

  out
}


# Internal helper functions --------------------------------------------------

# Helper function: defective density for an accumulator
#' @keywords internal
.cogmod_lba2_defectivedensity <- function(dt, v, s, A, b) {
  # dt: decision time(s) (must be positive)
  # v: mean drift for the accumulator
  # s: standard deviation for the accumulator's drift
  # A: start point range (sigmabias)
  # b: decision threshold (sigmabias + boundary)
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
.cogmod_lba2_cumulative <- function(dt, v, s, A, b) {
  n1 <- (b - A - v * dt) / (dt * s)
  n2 <- (b - v * dt) / (dt * s)
  F_val <- 1 +
    ((b - A - v * dt) / A) * stats::pnorm(n1) -
    ((b - v * dt) / A) * stats::pnorm(n2) +
    ((dt * s) / A) * (stats::dnorm(n1) - stats::dnorm(n2))
  F_val
}


# Stanvars ----------------------------------------------------------------

#' @keywords internal
.cogmod_lba2_lpdf <- function() {
  "
// Log probability density function for the LBA model using the new parametrization.
// Y: observed reaction time.
// mu: mean drift for accumulator 0.
// driftone: mean drift for accumulator 1.
// sigmazero, sigmaone: standard deviations of drift for accumulators 0 and 1.
// sigmabias: starting-point range, A = sigmabias.
// boundary: threshold offset such that b = sigmabias + boundary.
// tau: scale factor for non-decision time (0-1, scaled by minrt).
// minrt: minimum possible reaction time (used to scale tau).
// dec: decision indicator (0 if accumulator 0 wins; 1 if accumulator 1 wins).
real cogmod_lba2_lpdf(real Y,
              real mu,
              real driftone,
              real sigmazero,
              real sigmaone,
              real sigmabias,
              real boundary,
              real tau,
              real minrt,
              int dec) {
  // Derived parameters.
  real A = sigmabias;             // Starting-point range.
  real ndt = tau * minrt;         // Non-decision time.

  // Check if the observed RT is less than or equal to ndt.
  if (Y <= ndt) return negative_infinity();

  // Now compute the adjusted decision time.
  real t = Y - ndt;

  // --- Early Parameter Validity Checks ---
  if (sigmazero <= 0 || sigmaone <= 0) return negative_infinity();
  if (dec != 0 && dec != 1) return negative_infinity();

  // --- Standard Case: Both drifts are finite.
  // Define winning and losing accumulator parameters based on dec.
  real v_win;
  real s_win;
  real v_loss;
  real s_loss;
  if (dec == 0) {
    v_win  = mu;
    s_win  = sigmazero;
    v_loss = driftone;
    s_loss = sigmaone;
  } else {  // dec == 1.
    v_win  = driftone;
    s_win  = sigmaone;
    v_loss = mu;
    s_loss = sigmazero;
  }

  real inv_A = inv(A);
  real INV_SQRT_2PI = 0.3989422804014327;

  // --- Winning accumulator: defective density ---
  real st_win = fmax(s_win * t, 1e-10);
  real inv_st_win = inv(st_win);
  real z1w = (boundary - v_win * t) * inv_st_win;        // = (b - A - v_win*t) / st_win
  real z2w = z1w + A * inv_st_win;                 // = (b - v_win*t) / st_win
  real Phi1w = Phi(z1w);
  real Phi2w = Phi(z2w);
  real phi1w = exp(-0.5 * square(z1w)) * INV_SQRT_2PI;
  real phi2w = exp(-0.5 * square(z2w)) * INV_SQRT_2PI;
  real f_win = inv_A * (v_win * (Phi2w - Phi1w) + s_win * (phi1w - phi2w));

  // --- Losing accumulator: survival probability ---
  real st_loss = fmax(s_loss * t, 1e-10);
  real inv_st_loss = inv(st_loss);
  real z1l = (boundary - v_loss * t) * inv_st_loss;      // = (b - A - v_loss*t) / st_loss
  real z2l = z1l + A * inv_st_loss;                // = (b - v_loss*t) / st_loss
  real Phi1l = Phi(z1l);
  real Phi2l = Phi(z2l);
  real phi1l = exp(-0.5 * square(z1l)) * INV_SQRT_2PI;
  real phi2l = exp(-0.5 * square(z2l)) * INV_SQRT_2PI;
  real g1 = z1l * Phi1l + phi1l;
  real g2 = z2l * Phi2l + phi2l;
  real S_loss = (st_loss * inv_A) * (g2 - g1);     // = 1 - CDF_loss(t)

  // --- Combine into the joint (defective) density ---
  // Floor the *joint* density (not each factor separately) at double
  // precision epsilon, mirroring the R reference implementation `dcogmod_lba2()`.
  // Flooring each factor individually (e.g. at 1e-10) is much more
  // aggressive and can bias the log-density in the tails relative to R.
  real EPS = 2.220446049250313e-16;  // .Machine$double.eps
  real dens = fmax(f_win, 0) * fmin(fmax(S_loss, 0), 1);
  real log_lik = log(fmax(dens, EPS));

  if (is_nan(log_lik) || is_inf(log_lik)) return negative_infinity();
  return log_lik;
}
"
}


#' @rdname rcogmod_lba2
#' @examples
#' # You can expose the lpdf function as follows:
#' # cogmod_lba2_lpdf <- cogmod_lba2_lpdf_expose()
#' # cogmod_lba2_lpdf(...)
#'
#' @export
cogmod_lba2_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # This wraps the (new) Stan function cogmod_lba2_lpdf using our new parametrization.
  stancode <- paste0(
    "functions {
",
    .cogmod_lba2_lpdf(),
    "
}"
  )

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$cogmod_lba2_lpdf
}

#' @rdname rcogmod_lba2
#' @export
cogmod_lba2_stanvars <- function() {
  brms::stanvar(scode = .cogmod_lba2_lpdf(), block = "functions")
}


#' @rdname rcogmod_lba2
#' @param link_mu,link_driftone,link_sigmazero,link_sigmaone,link_sigmabias,link_boundary,link_tau,link_minrt Link functions for the parameters.
#' @export
cogmod_lba2 <- function(
  link_mu = "identity",
  link_driftone = "identity",
  link_sigmazero = "softplus",
  link_sigmaone = "softplus",
  link_sigmabias = "softplus",
  link_boundary = "softplus",
  link_tau = "logit",
  link_minrt = "identity"
) {
  brms::custom_family(
    name = "cogmod_lba2",
    dpars = c(
      "mu",
      "driftone",
      "sigmazero",
      "sigmaone",
      "sigmabias",
      "boundary",
      "tau",
      "minrt"
    ),
    links = c(
      link_mu,
      link_driftone,
      link_sigmazero,
      link_sigmaone,
      link_sigmabias,
      link_boundary,
      link_tau,
      link_minrt
    ),
    lb = c(NA, NA, 0, 0, 0, 0, 0, 0),
    ub = c(NA, NA, NA, NA, NA, NA, 1, NA),
    type = "real",
    vars = c("dec[n]")
  )
}

# brms Post-processing Functions ------------------------------------------

#' @rdname rcogmod_lba2
#' @inheritParams rcogmod_lnr
#' @export
log_lik_cogmod_lba2 <- function(i, prep) {
  # Extract observed reaction time for observation i
  y <- prep$data$Y[i]
  if (is.na(y)) {
    return(NA_real_)
  } # Handle missing RTs

  # Get the parameter draws for observation i:
  driftzero <- brms::get_dpar(prep, "mu", i = i)
  driftone <- brms::get_dpar(prep, "driftone", i = i)
  sigmazero <- brms::get_dpar(prep, "sigmazero", i = i)
  sigmaone <- brms::get_dpar(prep, "sigmaone", i = i)
  sigmabias <- brms::get_dpar(prep, "sigmabias", i = i)
  boundary <- brms::get_dpar(prep, "boundary", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i) # Minimum possible reaction time

  # Calculate non-decision time from tau and minrt
  ndt <- tau * minrt

  # Get decision indicator for observation i (should be 0 or 1).
  response <- prep$data[["dec"]][i]
  if (!response %in% c(0, 1)) {
    warning(
      "Response (dec) must be 0 or 1. Got: ",
      response,
      " for observation ",
      i
    )
    return(rep(-Inf, length(driftzero)))
  }

  # Compute log-likelihood using our R density function 'dcogmod_lba2'
  ll <- dcogmod_lba2(
    x = y,
    response = response,
    driftzero = driftzero,
    driftone = driftone,
    sigmazero = sigmazero,
    sigmaone = sigmaone,
    sigmabias = sigmabias,
    boundary = boundary,
    ndt = ndt,
    log = TRUE
  )
  ll # Return vector of log-likelihoods (one per draw)
}

#' @rdname rcogmod_lba2
#' @export
posterior_predict_cogmod_lba2 <- function(i, prep, ...) {
  # Get the parameter draws for observation i:
  driftzero <- brms::get_dpar(prep, "mu", i = i)
  driftone <- brms::get_dpar(prep, "driftone", i = i)
  sigmazero <- brms::get_dpar(prep, "sigmazero", i = i)
  sigmaone <- brms::get_dpar(prep, "sigmaone", i = i)
  sigmabias <- brms::get_dpar(prep, "sigmabias", i = i)
  boundary <- brms::get_dpar(prep, "boundary", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  # Calculate non-decision time from tau and minrt
  ndt <- tau * minrt

  # Number of posterior draws:
  n_draws <- length(driftzero)

  # Generate predictions using the R simulation function 'rcogmod_lba2'
  # It is assumed that 'rcogmod_lba2' is vectorized over parameters.
  sim_data <- rcogmod_lba2(
    n = n_draws,
    driftzero = driftzero,
    driftone = driftone,
    sigmazero = sigmazero,
    sigmaone = sigmaone,
    sigmabias = sigmabias,
    boundary = boundary,
    ndt = ndt
  )
  as.matrix(sim_data)
}

#' @rdname rcogmod_lba2
#' @export
posterior_epred_cogmod_lba2 <- function(prep) {
  stop(
    "Calculating the posterior expected prediction (epred) for the LBA model ",
    "is computationally prohibitive within this framework.\n",
    "Please use `posterior_predict()` to obtain draws from the posterior ",
    "predictive distribution and calculate summaries manually if needed."
  )
}
