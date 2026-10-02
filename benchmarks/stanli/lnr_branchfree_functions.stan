functions {
// The LNR likelihood at sigmabias = 0 with no control flow that depends on a
// parameter, which is what stanli 0.18 needs (README.md, "Why it refuses").
// Same arguments as cogmod_lnr_lpdf() minus sigmabias. Benchmark use only.

// log Phi(x) on the erfc route alone. The package's cogmod_log_Phi() switches
// to an asymptotic series below x = -25, a branch on a parameter that stanli
// cannot compile; without it the value underflows to log(0) near x = -38 and
// the gradient goes NaN (the LNR bug fixed in 0.3.3). Fine on the vignette's
// data, which never gets there; not fine in general.
real bf_log_Phi(real x) {
  return log(0.5 * erfc(-x * 0.7071067811865476));
}

real cogmod_lnr_lpdf(real Y, real mu, real nuone, real sigmazero, real sigmaone, real ndt, real poutlier, int dec) {
  real lp_out = 0.6904993792294275 - 12.5 * square(Y);
  // Winner and loser picked by data, which stanli folds when it builds the
  // model; a branch on `dec` is not a runtime-control region.
  real nu_w = dec == 0 ? mu : nuone;
  real s_w = dec == 0 ? sigmazero : sigmaone;
  real nu_l = dec == 0 ? nuone : mu;
  real s_l = dec == 0 ? sigmaone : sigmazero;
  // Replaces `if (t_adj <= 0) return log(poutlier) + lp_out;`. Comparisons are
  // unsupported in stanli even outside a region, so there is no indicator to
  // multiply by. At t = 1e-300 s the decision component is around -1e6 for
  // any plausible sigma, so it drops out of log_mix() exactly in double
  // precision, with a zero adjoint. It stops being exact only for sigmas
  // above ~20 on the log scale, which the priors put nowhere.
  real t = fmax(Y - ndt, 1e-300);
  real lp_dec = lognormal_lpdf(t | -nu_w, s_w)
              + bf_log_Phi((-nu_l - log(t)) / s_l);
  return log_mix(poutlier, lp_out, lp_dec);
}
}
