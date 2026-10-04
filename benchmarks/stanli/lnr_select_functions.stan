functions {
// The LNR likelihood with the package's numerics and no control flow at all,
// which is what stanli 0.18.1 runs fast (README.md). Each branch of
// cogmod_lnr_lpdf() becomes a "select" made of fmin()/fmax() and step():
// every arm is evaluated on an input clamped into its own domain, so it is
// finite everywhere, and step() picks one. The arm not picked sits at a
// constant (its clamp edge) with weight 0, so it passes back an exact zero -
// never 0 * inf - and the picked arm gets the whole gradient. Same arguments
// as cogmod_lnr_lpdf() minus sigmabias (sigmabias = 0 only).

// cogmod_log_Phi() as a select: the asymptotic series below x = -25, erfc
// from -25 up, as the package's `x < -25` splits it. The two meet at -25 to
// the last bit in value (see the package's comment), not in derivative (186
// ULP apart), so the split point has to be the same: Stan's step(0) is 1, so
// `step(-25 - x)` would pick the series at -25 itself, and `1 - step(x + 25)`
// does not. One difference remains: above x = 0 this takes
// log(0.5 * erfc(-x / sqrt2)) rather than log1p(-0.5 * erfc(x / sqrt2)). Same
// gradient, and the absolute error stays below 1e-16, but the relative
// accuracy goes: 41-49 ULP at x = 3, and 0 in place of -1.13e-19 at x = 9
// (stanli#422). A third arm would restore it at the price of a second erfc.
// Measured at 1.1 ms against 0.8-0.9 ms for the bare erfc route on the
// vignette model (stanli 0.18.1), so the select is nearly free.
real sel_log_Phi(real x) {
  real xl = fmin(x, -25);
  real xh = fmax(x, -25);
  real z = inv_square(xl);
  real series = 1 + z * (-1 + z * (3 + z * (-15 + z * (105 - 945 * z))));
  real lo = -0.5 * square(xl) - log(-xl) - 0.91893853320467274 + log(series);
  real hi = log(0.5 * erfc(-xh * 0.7071067811865476));
  real w = 1 - step(x + 25);  // 1 strictly below -25
  return w * lo + (1 - w) * hi;
}

// No checks, which is not the same as needing none. dec in {0, 1} is
// guaranteed only by .cogmod_checkdata(), which runs from cogmod_priors()
// alone, and Y > 0 likewise: brms declares Y with no lower bound. A softplus
// link underflows to exactly 0 below about -745 on the link scale; then
// lognormal_lpdf() rejects a zero s_w, but a zero s_l divides by 0 unchecked.
// Data-only branches cost 4x in stanli 0.18.1 (4.8 ms against 1.1), although
// it has the data when it builds the model.
real cogmod_lnr_lpdf(real Y, real mu, real nuone, real sigmazero, real sigmaone, real ndt, real poutlier, int dec) {
  real lp_out = 0.6904993792294275 - 12.5 * square(Y);
  // Winner and loser picked by data; a ternary on data costs nothing.
  real nu_w = dec == 0 ? mu : nuone;
  real s_w = dec == 0 ? sigmazero : sigmaone;
  real nu_l = dec == 0 ? nuone : mu;
  real s_l = dec == 0 ? sigmaone : sigmazero;
  // `if (t_adj <= 0) return log(poutlier) + lp_out;` as a mask inside
  // log_mix(): at or below ndt the decision component is evaluated at a
  // clamped time, where it is finite, and pushed to -1e300, where log_mix()
  // gives it an exact zero weight and adjoint. step(ndt - Y) is 1 at Y = ndt
  // itself, as `t_adj <= 0` is; step(Y - ndt) - 1 would keep the decision
  // component there. Not exact at poutlier = 0: `orig` returns -inf, this
  // -1e300 with a NaN gradient (log_mix()'s partials take 0 * inf); either
  // way the proposal is rejected. Written as a blend instead -
  // w * log_mix(...) + (1 - w) * (...) - or with log_sum_exp() in place of
  // log_mix(), the same value cost 5-6 ms in stanli 0.18.1: the mixture had
  // to stay one log_mix() for stanli to keep it fast.
  real t = fmax(Y - ndt, 1e-300);
  real lp_dec = lognormal_lpdf(t | -nu_w, s_w) + sel_log_Phi((-nu_l - log(t)) / s_l);
  return log_mix(poutlier, lp_out, lp_dec - step(ndt - Y) * 1e300);
}
}
