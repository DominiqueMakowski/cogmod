// From https://github.com/DominiqueMakowski/cogmod/issues/5 (Sean Talts,
// 2026-10-04), copied verbatim below this header: the brms program for the
// vignette LNR (lnr_orig.stan) with the likelihood rewritten by hand - rows
// split by `dec` in transformed data, vector arithmetic per group, each
// parameter branch a 0/1 blend of clamped arms, scalar sigmas. Data,
// parameters, priors and generated quantities are brms's, unchanged
// (bench.R checks). Benchmark variant `rw`; `rw_v` is derived from it.
// The LNR benchmark program from brms 2.23.0, rewritten by hand (see above).
functions {
  vector lnr_log_Phi(vector x, int n) {
    vector[n] wh = ceil(fmin(fmax(x, 0.0), 1.0));
    vector[n] he = 0.5 * erfc(((2 * wh - 1) .* fmax(x, -25.0)) * 0.7071067811865476);
    vector[n] res = (1 - wh) .* log(he + wh) + wh .* log1p(-he);
    if (min(x) < -25.0) {
      vector[n] wl = ceil(fmin(fmax(-25.0 - x, 0.0), 1.0));
      vector[n] xl = fmin(x, -25.0);
      vector[n] z = inv_square(xl);
      vector[n] series = 1 + z .* (-1 + z .* (3 + z .* (-15 + z .* (105 - 945 * z))));
      res = wl .* (-0.5 * square(xl) - log(-xl) - 0.91893853320467274 + log(series))
            + (1 - wl) .* res;
    }
    return res;
  }

  vector lnr_group(int n, vector Yg, vector lpout, vector nu_w, real s_w, vector nu_l,
                   real s_l, vector ndt, real lp, real l1mp, real tfloor,
                   vector ones) {
    vector[n] t = Yg - ndt;
    vector[n] w = ceil(fmin(fmax(t, 0.0), 1.0));
    vector[n] lt = log(fmax(t, tfloor));
    vector[n] sw = s_w * ones;
    vector[n] isq = square(inv(sw));
    vector[n] lognorm = -0.91893853320467274 - 0.5 * (square(lt + nu_w) .* isq)
                        - log(sw) - lt;
    vector[n] lpdec = lognorm + lnr_log_Phi((-nu_l - lt) ./ (s_l * ones), n);
    vector[n] a = lp * ones + lpout;
    vector[n] mix = log_sum_exp(a, l1mp * ones + lpdec);
    return w .* mix + (1 - w) .* a;
  }
}
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  array[N] int<lower=0,upper=1> dec;  // decisions
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  int<lower=1> Kc;  // number of population-level effects after centering
  int<lower=1> K_nuone;  // number of population-level effects
  matrix[N, K_nuone] X_nuone;  // population-level design matrix
  int<lower=1> Kc_nuone;  // number of population-level effects after centering
  int<lower=1> K_ndt;  // number of population-level effects
  matrix[N, K_ndt] X_ndt;  // population-level design matrix
  int<lower=1> Kc_ndt;  // number of population-level effects after centering
  int prior_only;  // should the likelihood be ignored?
}

transformed data {
  real min_Y = min(Y);
  matrix[N, Kc] Xc;  // centered version of X without an intercept
  vector[Kc] means_X;  // column means of X before centering
  matrix[N, Kc_nuone] Xc_nuone;  // centered version of X_nuone without an intercept
  vector[Kc_nuone] means_X_nuone;  // column means of X_nuone before centering
  matrix[N, Kc_ndt] Xc_ndt;  // centered version of X_ndt without an intercept
  vector[Kc_ndt] means_X_ndt;  // column means of X_ndt before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
  for (i in 2:K_nuone) {
    means_X_nuone[i - 1] = mean(X_nuone[, i]);
    Xc_nuone[, i - 1] = X_nuone[, i] - means_X_nuone[i - 1];
  }
  for (i in 2:K_ndt) {
    means_X_ndt[i - 1] = mean(X_ndt[, i]);
    Xc_ndt[, i - 1] = X_ndt[, i] - means_X_ndt[i - 1];
  }
  int n_bad = 0;
  int N1 = 0;
  for (n in 1:N) {
    if (Y[n] <= 0) n_bad += 1;
    else N1 += dec[n];
  }
  int N0 = N - n_bad - N1;
  array[N0] int idx0;
  array[N1] int idx1;
  {
    int k0 = 0;
    int k1 = 0;
    for (n in 1:N) {
      if (Y[n] > 0) {
        if (dec[n] == 0) {
          k0 += 1;
          idx0[k0] = n;
        } else {
          k1 += 1;
          idx1[k1] = n;
        }
      }
    }
  }
  vector[N0] Y0 = Y[idx0];
  vector[N1] Y1 = Y[idx1];
  vector[N0] lpout0 = 0.6904993792294275 - 12.5 * square(Y0);
  vector[N1] lpout1 = 0.6904993792294275 - 12.5 * square(Y1);
  vector[N0] ones0 = rep_vector(1.0, N0);
  vector[N1] ones1 = rep_vector(1.0, N1);
  real tfloor = 1e-300;
  if (N0 + N1 > 0) tfloor = fmax(1e-300, min(append_row(Y0, Y1)) * 8.673617379884035e-19);
}
parameters {
  vector[Kc] b;  // regression coefficients
  real Intercept;  // temporary intercept for centered predictors
  vector[Kc_nuone] b_nuone;  // regression coefficients
  real Intercept_nuone;  // temporary intercept for centered predictors
  real Intercept_sigmazero;  // temporary intercept for centered predictors
  real Intercept_sigmaone;  // temporary intercept for centered predictors
  vector[Kc_ndt] b_ndt;  // regression coefficients
  real Intercept_ndt;  // temporary intercept for centered predictors
  real<lower=0,upper=1> poutlier;
}
transformed parameters {
  real sigmabias = 0;
  // prior contributions to the log posterior
  real lprior = 0;
  lprior += student_t_lpdf(Intercept | 3, 0.6, 2.5);
  lprior += normal_lpdf(b_nuone | 0, 0.5);
  lprior += normal_lpdf(Intercept_nuone | 0.7, 1.5);
  lprior += normal_lpdf(Intercept_sigmazero | 0, 1);
  lprior += normal_lpdf(Intercept_sigmaone | 0, 1);
  lprior += normal_lpdf(b_ndt | 0, 0.2);
  lprior += normal_lpdf(Intercept_ndt | -1.2, 0.5);
  lprior += exponential_lpdf(poutlier | 100)
    - 1 * exponential_lcdf(1 | 100);
}
model {
  // likelihood including constants
  if (!prior_only) {
    vector[N] mu = Intercept + Xc * b;
    vector[N] nuone = Intercept_nuone + Xc_nuone * b_nuone;
    vector[N] ndt = exp(Intercept_ndt + Xc_ndt * b_ndt);
    real sigmazero = log1p_exp(Intercept_sigmazero);
    real sigmaone = log1p_exp(Intercept_sigmaone);
    if (n_bad > 0) {
      target += negative_infinity();
    } else {
      int bad = sigmazero <= 0 || sigmaone <= 0 || sigmabias < 0 || poutlier < 0
                || poutlier > 1 || (poutlier <= 0 && max(ndt - Y) >= 0);
      real ok = 1 - bad;
      real sz = sigmazero * ok + bad;
      real so = sigmaone * ok + bad;
      real pout = (poutlier + (poutlier <= 0) * 4.9406564584124654e-324) * ok
                  + 0.5 * bad;
      real lp = log(pout);
      real l1mp = log1m(pout);
      vector[N] per;
      per[idx0] = lnr_group(N0, Y0, lpout0, mu[idx0], sz, nuone[idx0], so,
                            ndt[idx0], lp, l1mp, tfloor, ones0);
      per[idx1] = lnr_group(N1, Y1, lpout1, nuone[idx1], so, mu[idx1], sz,
                            ndt[idx1], lp, l1mp, tfloor, ones1);
      target += ok * sum(per) + log(ok);
    }
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept - dot_product(means_X, b);
  // actual population-level intercept
  real b_nuone_Intercept = Intercept_nuone - dot_product(means_X_nuone, b_nuone);
  // actual population-level intercept
  real b_sigmazero_Intercept = Intercept_sigmazero;
  // actual population-level intercept
  real b_sigmaone_Intercept = Intercept_sigmaone;
  // actual population-level intercept
  real b_ndt_Intercept = Intercept_ndt - dot_product(means_X_ndt, b_ndt);
}


