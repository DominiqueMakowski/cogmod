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
// above, meeting at -25 to the last bit (see the package's comment), so the
// blend has no step. One difference: above x = 0 this takes
// log(0.5 * erfc(-x / sqrt2)) rather than log1p(-0.5 * erfc(x / sqrt2)), which
// loses relative accuracy only in a value already within 1e-16 of 0 (absolute
// error <= 1e-16) and gives the same gradient. A third arm would restore it
// at the price of a second erfc. Measured at 1.1 ms against 0.8-0.9 ms for
// the bare erfc route on the vignette model, so the select is nearly free.
real sel_log_Phi(real x) {
  real xl = fmin(x, -25);
  real xh = fmax(x, -25);
  real z = inv_square(xl);
  real series = 1 + z * (-1 + z * (3 + z * (-15 + z * (105 - 945 * z))));
  real lo = -0.5 * square(xl) - log(-xl) - 0.91893853320467274 + log(series);
  real hi = log(0.5 * erfc(-xh * 0.7071067811865476));
  real w = step(-25 - x);  // 1 below -25
  return w * lo + (1 - w) * hi;
}

// No checks. The parameter ones are guaranteed by brms's links and bounds;
// the data ones (dec in {0, 1}, Y > 0) by the data block and
// .cogmod_checkdata(). Even branches on data alone cost 4x in stanli (4.8 ms
// against 1.1), although it has the data when it builds the model.
real cogmod_lnr_lpdf(real Y, real mu, real nuone, real sigmazero, real sigmaone, real ndt, real poutlier, int dec) {
  real lp_out = 0.6904993792294275 - 12.5 * square(Y);
  // Winner and loser picked by data; a ternary on data costs nothing.
  real nu_w = dec == 0 ? mu : nuone;
  real s_w = dec == 0 ? sigmazero : sigmaone;
  real nu_l = dec == 0 ? nuone : mu;
  real s_l = dec == 0 ? sigmaone : sigmazero;
  // `if (t_adj <= 0) return log(poutlier) + lp_out;` as a mask inside
  // log_mix(): below ndt the decision component is evaluated at a clamped
  // time, where it is finite, and pushed to -1e300, where log_mix() gives it
  // an exact zero weight and adjoint. Written as a blend instead -
  // w * log_mix(...) + (1 - w) * (...) - or with log_sum_exp() in place of
  // log_mix(), the same value cost 5-6 ms: the mixture has to stay one
  // log_mix() for stanli to keep it fast.
  real t = fmax(Y - ndt, 1e-300);
  real lp_dec = lognormal_lpdf(t | -nu_w, s_w) + sel_log_Phi((-nu_l - log(t)) / s_l);
  return log_mix(poutlier, lp_out, lp_dec + (step(Y - ndt) - 1) * 1e300);
}
}
