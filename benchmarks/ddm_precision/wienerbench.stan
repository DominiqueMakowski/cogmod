data {
  int<lower=1> N;
  vector<lower=0>[N] y;
  int<lower=1,upper=5> branch;
  real sw_lo; real sw_hi; real st0_lo; real st0_hi;
  real<lower=0> prec;
}
parameters {
  real v;
  real<lower=0.3, upper=3> a;
  real<lower=0.25, upper=0.75> w;
  real<lower=0.001, upper=2> sv;
  real<lower=sw_lo, upper=sw_hi> sw;
  real<lower=st0_lo, upper=st0_hi> st0;
}
model {
  v ~ normal(1, 1); a ~ normal(1.5, 0.5); w ~ normal(0.5, 0.1); sv ~ normal(0.5, 0.5);
  sw ~ normal(0.5 * (sw_lo + sw_hi), 1);
  st0 ~ normal(0.5 * (st0_lo + st0_hi), 1);
  if (branch == 1)      { for (n in 1:N) target += wiener_lpdf(y[n] | a, 0.1, w, v); }
  else if (branch == 2) { for (n in 1:N) target += wiener_lpdf(y[n] | a, 0.1, w, v, sv); }
  else if (branch == 3) { for (n in 1:N) target += wiener_lpdf(y[n] | a, 0.1, w, v, sv, sw, st0, prec); }
  else if (branch == 4) { for (n in 1:N) target += wiener_lpdf(y[n] | a, 0.1, w, v, sv, sw, 0, prec); }
  else                  { for (n in 1:N) target += wiener_lpdf(y[n] | a, 0.1, w, v, sv, 0, st0, prec); }
}
