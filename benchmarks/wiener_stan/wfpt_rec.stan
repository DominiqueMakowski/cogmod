// Variant rec: the same series and term counts as wfpt_direct.stan, with
// the per-term transcendental calls replaced by products.
//
// Log first-passage density of the driftless, unit-boundary diffusion at its
// lower boundary, started at wl, at rescaled time u; w = 1 - wl is passed
// alongside because the caller has it exactly, where 1 - wl would round.
real cogmod_wfpt0(real u, real wl, real w) {
  if (u < 0.20276100128382937) {
    // Small time: sum_k (wl + 2k) exp(-2k(k + wl) / u), k = -m..m. The terms
    // are products of e1 (the k = -1 factor), a1 (k = 1) and q = exp(-4 / u):
    // a[k+1] = a[k] a1 q^k, and b[j+1] = b[j] e1 q^j for k = -j. The k = 0
    // and k = -1 terms are summed as -expm1(-2w/u) - w (1 + e1): written as
    // wl + (wl - 2) e1 they cancel as w -> 0, both near 1, sum near 2w/u.
    real r = 2 / u;
    real e1 = exp(-w * r);
    real a1 = exp(-(1 + wl) * r);
    real q = a1 * e1;
    real s = -expm1(-w * r) - w * (1 + e1);
    real ak = a1;
    real bj = e1;
    real qk = 1;
    int m = u < 0.073189601146901873 ? 2 : 3;
    for (k in 1:m) {
      s += (wl + 2 * k) * ak;
      qk *= q;
      if (k < m) {
        bj *= e1 * qk;
        s += (wl - 2 * (k + 1)) * bj;
        ak *= a1 * qk;
      }
    }
    if (!(s > 0)) return negative_infinity();
    return log(s) - square(wl) / (2 * u) - 0.5 * log(2 * pi()) - 1.5 * log(u);
  }
  // Large time: sum_k k exp(-(k^2 - 1) pi^2 u / 2) sin(k pi wl), k = 1..K.
  // exp(-(k^2 - 1) pi^2 u / 2) = E^(k^2 - 1), by g[k+1] = g[k] E^(2k + 1);
  // sin(k theta) by rotation from one sin() and one cos(). theta is pi times
  // whichever of wl and w is nearer 0, so a start point near either boundary
  // keeps its relative precision: sin(k pi wl) = (-1)^(k + 1) sin(k pi w).
  int K = u < 0.22671899719224897 ? 6
          : (u < 0.34879271234654191 ? 5
             : (u < 0.60757953464589354 ? 4
                : (u < 1.3274606720345303 ? 3 : 2)));
  real E = exp(-square(pi()) * u / 2);
  real E2 = square(E);
  real h = E2 * E;
  real g = 1;
  int flip = wl > 0.5;
  real theta = pi() * (flip ? w : wl);
  real st = sin(theta);
  real ct = cos(theta);
  real sk = st;
  real ck = ct;
  real s = st;
  int sg = 1;
  for (k in 2:K) {
    real sn = sk * ct + ck * st;
    ck = ck * ct - sk * st;
    sk = sn;
    g *= h;
    h *= E2;
    if (flip) sg = -sg;
    s += sg * k * g * sk;
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
