// Variant unrolled: wfpt_rec.stan with both series written out at fixed
// length and no loop - the small-time one at m = 3 (k = -3..3) and the
// large-time one at K = 6, the most either needs at 1e-12 on its side of the
// switch. The terms the shorter rules would skip are below 1e-12 of the sum,
// so the value moves in the last bits only. The point is stanli: a loop
// whose trip count depends on a parameter keeps the trial loop out of
// stanli's region map, which then re-records it at every new point.
//
// Log first-passage density of the driftless, unit-boundary diffusion at its
// lower boundary, started at wl, at rescaled time u; w = 1 - wl is passed
// alongside because the caller has it exactly, where 1 - wl would round.
real cogmod_wfpt0(real u, real wl, real w) {
  if (u < 0.20276100128382937) {
    // Small time: sum_k (wl + 2k) exp(-2k(k + wl) / u), with every term a
    // product of e1 (the k = -1 factor), a1 (k = 1) and q = exp(-4 / u); the
    // k = 0 and k = -1 terms through expm1(), which does not cancel as w -> 0.
    real r = 2 / u;
    real e1 = exp(-w * r);
    real a1 = exp(-(1 + wl) * r);
    real q = a1 * e1;
    real q2 = square(q);
    real b2 = square(e1) * q;
    real a2 = square(a1) * q;
    real s = -expm1(-w * r) - w * (1 + e1)
             + (wl + 2) * a1 + (wl - 4) * b2 + (wl + 4) * a2
             + (wl - 6) * (b2 * e1 * q2) + (wl + 6) * (a2 * a1 * q2);
    if (!(s > 0)) return negative_infinity();
    return log(s) - square(wl) / (2 * u) - 0.5 * log(2 * pi()) - 1.5 * log(u);
  }
  // Large time: sum_k k E^(k^2 - 1) sin(k pi wl), k = 1..6, E = exp(-pi^2 u / 2);
  // sin(k theta) by rotation, theta at whichever of wl and w is nearer 0, and
  // sin(k pi wl) = (-1)^(k + 1) sin(k pi w) when it is w.
  real E = exp(-square(pi()) * u / 2);
  real E2 = square(E);
  real g2 = E2 * E;
  real g3 = g2 * g2 * E2;
  real g4 = g3 * g2 * E2 * E2;
  real g5 = g4 * g3 * E;
  real g6 = g5 * g3 * E2 * E;
  int flip = wl > 0.5;
  int f = flip ? -1 : 1;
  real theta = pi() * (flip ? w : wl);
  real s1 = sin(theta);
  real c1 = cos(theta);
  real s2 = 2 * s1 * c1;
  real c2 = 1 - 2 * square(s1);
  real s3 = s2 * c1 + c2 * s1;
  real c3 = c2 * c1 - s2 * s1;
  real s4 = 2 * s2 * c2;
  real c4 = 1 - 2 * square(s2);
  real s5 = s4 * c1 + c4 * s1;
  real c5 = c4 * c1 - s4 * s1;
  real s6 = s5 * c1 + c5 * s1;
  real s = s1 + f * 2 * g2 * s2 + 3 * g3 * s3 + f * 4 * g4 * s4 + 5 * g5 * s5 + f * 6 * g6 * s6;
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
