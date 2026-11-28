// Variant direct: the R side's .ddm_lfpt0() / .ddm_lfpt() transcribed.
//
// Log first-passage density of the driftless, unit-boundary diffusion at its
// lower boundary, started at wl, at rescaled time u (Navarro & Fuss 2009,
// 1e-12 absolute on the density). Each series is summed with its leading
// term's exponent factored out, so the result stays finite far into both
// tails. The term counts are the R side's rule (small-time k = -m..m with
// m = max(2, ceil(ks / 2)); large-time K = max(2, ceil(kl))), read off
// thresholds on u computed once in R, so nothing is taped to choose them.
// `w` = 1 - wl is unused here; the signature is shared with wfpt_rec.stan.
real cogmod_wfpt0(real u, real wl, real w) {
  if (u < 0.20276100128382937) {
    int m = u < 0.073189601146901873 ? 2 : 3;
    real s = 0;
    for (k in -m:m) {
      real wk = wl + 2 * k;
      s += wk * exp(-(square(wk) - square(wl)) / (2 * u));
    }
    if (!(s > 0)) return negative_infinity();
    return log(s) - square(wl) / (2 * u) - 0.5 * log(2 * pi() * u^3);
  }
  int K = u < 0.22671899719224897 ? 6
          : (u < 0.34879271234654191 ? 5
             : (u < 0.60757953464589354 ? 4
                : (u < 1.3274606720345303 ? 3 : 2)));
  real s = 0;
  for (k in 1:K) {
    s += k * exp(-(k^2 - 1) * square(pi()) * u / 2) * sin(k * pi() * wl);
  }
  if (!(s > 0)) return negative_infinity();
  return log(pi()) + log(s) - square(pi()) * u / 2;
}

// Stan's wiener_lpdf() convention: the density at the upper boundary, with w
// measured from the lower one. That is the lower boundary of the reflected
// process started at 1 - w with drift -v; sv integrates out in closed form.
//
// Outside 1e-140 < t / a^2 < 1e140 the log density is below -1e100 or so,
// and the partials of t / a^2 itself overflow on the tape (u / a^2 = 4e380
// at a = 8e-96) while the value is still finite: reverse mode multiplies a
// zero adjoint into the infinite partial, and the whole gradient is NaN
// under any mixture weight. So the answer there is the -inf it rounds to,
// as a constant, tested on products whose partials cannot overflow.
real cogmod_wfpt_lpdf(real t, real a, real w, real v, real sv) {
  real a2 = square(a);
  if (t > 1e140 * a2 || t < 1e-140 * a2) return negative_infinity();
  real wl = 1 - w;
  real lf = cogmod_wfpt0(t / a2, wl, w) - 2 * log(a);
  if (sv == 0) return lf + v * a * wl - square(v) * t / 2;
  real s2 = square(sv);
  return lf + (s2 * square(a * wl) + 2 * v * a * wl - square(v) * t) / (2 * (1 + s2 * t))
         - 0.5 * log1p(s2 * t);
}
