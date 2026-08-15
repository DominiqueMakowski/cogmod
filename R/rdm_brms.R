# Stanvars ----------------------------------------------------------------

#' @keywords internal
.rdm_lpdf <- function() {
  "
// ---------------------------------------------------------------------------
// Two-accumulator Racing Diffusion Model (Tillman, Van Zandt & Logan, 2020).
//
// Two Wald accumulators race to a common threshold b = bs + sigmabias, each
// starting from z ~ Uniform(0, sigmabias). The observed trial contributes the
// winner's density times the loser's survival.
//
// Everything is carried in log space. That is not stylistic: for moderate
// drift the survival underflows to exactly zero (at drift 6 it does so by
// t = 4), and `log(0)` hands the sampler a zero gradient, which stalls it
// silently rather than erroring. The same applies to the density, whose
// individual terms underflow long before the density itself does.
// ---------------------------------------------------------------------------

// Below this ratio of start-point range to sqrt(t), the threshold difference
// quotient cancels; the midpoint (plain Wald) limit is then both cheaper and
// more accurate, being second-order in the range.
// Below RDM_EPS_V the survival's 1 / (2 * drift) factor takes its driftless
// limit instead.

// log(u * Phi(u) + phi(u)) = log of the antiderivative of Phi.
// For u < -10 both terms underflow *and* they cancel to order u^-2, so factor
// out phi(u) and expand 1 - |u| * MillsRatio(|u|) asymptotically.
real rdm_log_g(real u) {
  if (u > -10) {
    return log(u * Phi(u) + exp(-0.5 * square(u)) * 0.3989422804014327);
  }
  real x2 = square(u);
  real s = 0;
  real term = 1;
  for (j in 1:12) {
    term *= (2.0 * j - 1) / x2;
    s += (j % 2 == 1) ? term : -term;
  }
  return -0.5 * square(u) - 0.9189385332046727 + log(s);
}

// log(Phi(b) - Phi(a)) for b >= a, using whichever tail keeps both arguments
// away from a saturating normal CDF.
//
// The upper tail is written as std_normal_lcdf(-a) rather than
// std_normal_lccdf(a): Stan's lccdf collapses to -inf once its argument passes
// about 8.3 (and is already wrong in the 3rd decimal at 8), whereas its lcdf
// stays accurate past -30. The two are mathematically identical, and the
// distinction is not academic here -- alpha reaches 10 for a reaction time only
// a few milliseconds above the non-decision time.
real rdm_log_diff_Phi(real a, real b) {
  if (b <= a) return negative_infinity();
  if (a >= 0) return log_diff_exp(std_normal_lcdf(-a | ), std_normal_lcdf(-b | ));
  if (b <= 0) return log_diff_exp(std_normal_lcdf(b | ), std_normal_lcdf(a | ));
  return log(Phi(b) - Phi(a));
}

// Add the signed term sign2 * exp(log2) to the running total (log1, sign1).
// Returns [log|total|, sign]. Mirrors the R helper `.log_sum_signed()`.
vector rdm_signed_add(real log1, real sign1, real log2, real sign2) {
  vector[2] out;
  if (sign2 == 0 || (is_inf(log2) && log2 < 0)) {
    out[1] = log1; out[2] = sign1; return out;
  }
  if (sign1 == 0 || (is_inf(log1) && log1 < 0)) {
    out[1] = log2; out[2] = sign2; return out;
  }
  if (sign1 == sign2) {
    out[1] = log_sum_exp(log1, log2); out[2] = sign1; return out;
  }
  if (log1 > log2) {
    out[1] = log_diff_exp(log1, log2); out[2] = sign1; return out;
  }
  if (log2 > log1) {
    out[1] = log_diff_exp(log2, log1); out[2] = sign2; return out;
  }
  out[1] = negative_infinity(); out[2] = 0; return out;
}

// log density of the first passage time t (net of non-decision time) for a
// diffusion with drift nu, threshold offset k, start-point range A.
real rdm_wald_ldens(real t, real nu, real k, real A) {
  if (t <= 0) return negative_infinity();
  real st = sqrt(t);

  if (A / st < 1e-4) {              // midpoint plain Wald
    real bm = k + 0.5 * A;
    return log(bm) - 0.5 * (1.8378770664093453 + 3 * log(t))
           - square(bm - nu * t) / (2 * t);
  }

  real alpha = (k - nu * t) / st;
  real beta  = (k + A - nu * t) / st;

  // T1 = drift * (Phi(beta) - Phi(alpha))
  real l1 = log(abs(nu)) + rdm_log_diff_Phi(alpha, beta);
  real s1 = nu > 0 ? 1 : (nu < 0 ? -1 : 0);

  // T2 = (phi(alpha) - phi(beta)) / sqrt(t); phi is even, so the closer
  // argument to zero wins and decides the sign.
  real lpa = -0.5 * square(alpha) - 0.9189385332046727;
  real lpb = -0.5 * square(beta)  - 0.9189385332046727;
  real s2; real l2;
  if (abs(alpha) <= abs(beta)) {
    s2 = 1;  l2 = log_diff_exp(lpa, lpb) - log(st);
  } else {
    s2 = -1; l2 = log_diff_exp(lpb, lpa) - log(st);
  }

  vector[2] r = rdm_signed_add(l1, s1, l2, s2);
  return r[1] - log(A);
}

// log survival of the first passage time t. Computed from the closed-form
// antiderivative of the Wald survival in the threshold, NOT as log(1 - CDF):
// the latter underflows to -inf at ordinary parameter values.
real rdm_wald_lsurv(real t, real nu, real k, real A) {
  if (t <= 0) return 0;             // log(1): nothing finishes before the ndt
  real st = sqrt(t);
  real b  = k + A;

  if (A / st < 1e-4) {              // midpoint plain Wald survival
    real bm = k + 0.5 * A;
    real m1 = std_normal_lcdf((bm - nu * t) / st | );
    real m2 = 2 * nu * bm + std_normal_lcdf(-(bm + nu * t) / st | );
    return m1 > m2 ? log_diff_exp(m1, m2) : negative_infinity();
  }

  real alpha = (k - nu * t) / st;
  real beta  = (b - nu * t) / st;

  if (abs(nu) < 1e-7) {             // driftless limit
    // S = 1 + (2 sqrt(t) / A) * (g(-b/st) - g(-k/st)); already carries its 1/A.
    real lgap = log_diff_exp(rdm_log_g(-k / st), rdm_log_g(-b / st));
    vector[2] rz = rdm_signed_add(0, 1,
                                  0.6931471805599453 + log(st) - log(A) + lgap,
                                  -1);
    return fmin(rz[1], 0);
  }

  real sv   = nu > 0 ? 1 : -1;
  real linv = -log(2 * abs(nu));
  real wa = -(k + nu * t) / st;
  real wb = -(b + nu * t) / st;

  // Six signed terms; positive overall. exp(2*nu*b) * Phi(w) is kept as a
  // single exponent so it stays finite when the factors separately overflow
  // and underflow.
  vector[2] r = rdm_signed_add(0.5 * log(t) + rdm_log_g(beta),  1,
                               0.5 * log(t) + rdm_log_g(alpha), -1);
  r = rdm_signed_add(r[1], r[2], linv + 2 * nu * b + std_normal_lcdf(wb | ), -sv);
  r = rdm_signed_add(r[1], r[2], linv + 2 * nu * k + std_normal_lcdf(wa | ),  sv);
  r = rdm_signed_add(r[1], r[2], linv + std_normal_lcdf(beta | ),  -sv);
  r = rdm_signed_add(r[1], r[2], linv + std_normal_lcdf(alpha | ), sv);

  return fmin(r[1] - log(A), 0);
}

// Log probability density for a single RDM observation.
// Y: observed reaction time.
// mu: drift rate of accumulator 0.
// driftone: drift rate of accumulator 1.
// sigmabias: start-point range A; z ~ Uniform(0, sigmabias).
// bs: threshold offset, so the threshold is b = bs + sigmabias.
// tau: non-decision time as a proportion (0-1) of minrt.
// minrt: minimum possible reaction time, used to scale tau.
// dec: which accumulator won (0 for mu, 1 for driftone).
real rdm_lpdf(real Y, real mu, real driftone, real sigmabias, real bs,
              real tau, real minrt, int dec) {
  if (bs <= 0 || sigmabias <= 0) return negative_infinity();
  if (minrt < 0 || tau < 0 || tau > 1) return negative_infinity();
  if (dec != 0 && dec != 1) return negative_infinity();

  real ndt = tau * minrt;
  real t = Y - ndt;
  if (t <= 0) return negative_infinity();

  real v_win  = dec == 0 ? mu : driftone;
  real v_loss = dec == 0 ? driftone : mu;

  return rdm_wald_ldens(t, v_win, bs, sigmabias)
       + rdm_wald_lsurv(t, v_loss, bs, sigmabias);
}
"
}


#' @rdname rrdm
#' @examples
#' # You can expose the lpdf function as follows:
#' # rdm_lpdf <- rdm_lpdf_expose()
#' # rdm_lpdf(0.5, 2, 1.5, 0.2, 0.5, 0.5, 0.2, 0)
#'
#' @export
rdm_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Wrap the function Stan block (done normally by brms on model compilation)
  stancode <- paste0(
    "functions {
",
    .rdm_lpdf(),
    "
}"
  )

  # `force_recompile` so that repeated calls in one session do not hit
  # cmdstanr's cache, which refuses to expose functions from a pre-compiled
  # model (same reason as `ddm_lpdf_expose()`).
  mod <- cmdstanr::cmdstan_model(
    cmdstanr::write_stan_file(stancode),
    force_recompile = TRUE
  )
  mod$expose_functions()
  mod$functions$rdm_lpdf
}


#' @rdname rrdm
#' @export
rdm_stanvars <- function() {
  brms::stanvar(scode = .rdm_lpdf(), block = "functions")
}


#' @rdname rrdm
#' @param link_mu,link_driftone,link_sigmabias,link_bs,link_tau,link_minrt Link
#'   functions for the parameters.
#' @details
#' The `brms` family names the drift of the first accumulator `mu` (as `brms`
#' requires) and that of the second `driftone`, and calls the start-point range
#' `sigmabias` to match [lba()], where it denotes the same quantity. Note that
#' this is *not* the same thing as `bias` in [ddm()], which is a relative
#' starting point in `[0, 1]`. Both drifts use a softplus link with a lower
#' bound of zero, following [rt_invgaussian()]: a Wald drift must be
#' non-negative for the accumulator to be a proper first passage time.
#'
#' Note that `sigmabias` and `bs` are only weakly identified from each other,
#' because they enter the threshold as the sum `b = bs + sigmabias` and trade
#' off almost freely: on simulated data with 4000 trials the profile
#' log-likelihood varies by only about 3 units as `sigmabias` ranges from 0 to
#' half the threshold, while `bs` slides to compensate. With flat priors the
#' sampler tends to wander down the `sigmabias -> 0` ridge (the plain Wald race)
#' and produce divergent transitions. A weakly informative prior on
#' `sigmabias` fixes this, for example
#' `brms::prior(normal(-1, 1), class = "Intercept", dpar = "sigmabias")`, which
#' under the softplus link is centred near `0.31`. The sum `bs + sigmabias` is
#' well identified either way, so it is the more trustworthy quantity to
#' interpret and to compare across conditions. The same caveat applies to
#' [lba()], which shares this parameterisation.
#' @export
rdm <- function(
  link_mu = "softplus",
  link_driftone = "softplus",
  link_sigmabias = "softplus",
  link_bs = "softplus",
  link_tau = "logit",
  link_minrt = "identity"
) {
  brms::custom_family(
    name = "rdm",
    dpars = c(
      "mu",
      "driftone",
      "sigmabias",
      "bs",
      "tau",
      "minrt"
    ),
    links = c(
      link_mu,
      link_driftone,
      link_sigmabias,
      link_bs,
      link_tau,
      link_minrt
    ),
    lb = c(0, 0, 0, 0, 0, 0),
    ub = c(NA, NA, NA, NA, 1, NA),
    type = "real", # Continuous outcome variable (RT)
    vars = "dec[n]" # Required additional data variable for the decision
  )
}


# brms Post-processing Functions ------------------------------------------

#' @rdname rrdm
#' @inheritParams rlnr
#' @export
log_lik_rdm <- function(i, prep, ...) {
  y <- prep$data$Y[i]
  dec <- prep$data[["dec"]][i]

  if (is.na(y) || is.na(dec)) {
    return(rep(NA_real_, prep$ndraws))
  }
  if (!dec %in% c(0, 1)) {
    warning(
      "Response ('dec') must be 0 or 1. Found: ",
      dec,
      " for observation ",
      i,
      ". Returning -Inf log-likelihood."
    )
    return(rep(-Inf, prep$ndraws))
  }

  mu <- brms::get_dpar(prep, "mu", i = i)
  driftone <- brms::get_dpar(prep, "driftone", i = i)
  sigmabias <- brms::get_dpar(prep, "sigmabias", i = i)
  bs <- brms::get_dpar(prep, "bs", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  drdm(
    x = rep(y, length.out = length(mu)),
    response = dec,
    vzero = mu,
    vone = driftone,
    bs = bs,
    bias = sigmabias,
    ndt = tau * minrt,
    log = TRUE
  )
}


#' @rdname rrdm
#' @inheritParams rlnr
#' @export
posterior_predict_rdm <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  driftone <- brms::get_dpar(prep, "driftone", i = i)
  sigmabias <- brms::get_dpar(prep, "sigmabias", i = i)
  bs <- brms::get_dpar(prep, "bs", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  sim_data <- rrdm(
    n = prep$ndraws,
    vzero = mu,
    vone = driftone,
    bs = bs,
    bias = sigmabias,
    ndt = tau * minrt
  )
  as.matrix(sim_data)
}


#' @rdname rrdm
#' @inheritParams rlnr
#' @export
posterior_epred_rdm <- function(prep) {
  stop(
    "Calculating the posterior expected prediction (epred) for the RDM model ",
    "is computationally prohibitive within this framework.\n",
    "Please use `posterior_predict()` to obtain draws from the posterior ",
    "predictive distribution and calculate summaries manually if needed."
  )
}
