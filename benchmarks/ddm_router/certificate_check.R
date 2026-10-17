# Does the proposed certificate imply that today's router sends a trial to
# the classic wiener_lpdf()? (README.md, "The proof".) Pure R, no Stan.
#
#   Rscript benchmarks/ddm_router/certificate_check.R
#
# router() is an R transcription of cogmod_ddm_log_density_scale() at
# sv = sw = st0 = 0 (R/model_ddm.R, .DDM_STAN_PRELUDE). 2026-10-16: bound
# -560.57; 2e6 points, 605002 certified, none with router <= -600, the
# lowest -532.44; every typical posterior trial certified.
router <- function(t, v, a, w) {
  w_lo <- pmax(w, 1e-12); w_hi <- pmin(w, 1 - 1e-12)
  t_lo <- pmax(t, 1e-12); tau_hi <- t / a^2; tau_lo <- t_lo / a^2
  ws <- pmin(pmax(sqrt(tau_hi), w_lo), w_hi)
  small <- log(ws) - ws^2 / (2 * tau_hi) - 0.5 * log(2 * pi * tau_hi^3)
  wl <- pmin(pmax(0.5, w_lo), w_hi)
  large <- log(pi) - pi^2 * tau_lo / 2 + log(sin(pi * wl))
  num <- -v^2 * t_lo - 2 * v * a * ifelse(v > 0, w_lo, w_hi)
  -2 * log(a) + pmin(small, large) + num / 2
}
T_ <- 50; A <- 100; w0 <- 1e-3; K <- 300
cert <- function(t, v, a, w) {
  a2 <- a^2
  t > 0.01 * a2 & t < T_ * a2 & a < A & w > w0 & w < 1 - w0 & v * (v * t + 2 * a * w) < 2 * K
}
lb <- -2 * log(A) + min(log(w0) - 50 - 0.5 * log(2 * pi) - 1.5 * log(T_),
                        log(pi) - pi^2 * T_ / 2 + log(sin(pi * w0))) - K
cat("analytic lower bound inside the certificate:", round(lb, 2), "(must be > -600)\n")

set.seed(1)
n <- 2e6
# Log-uniform over wide ranges, so that the certificate's edges are hit often.
a <- exp(runif(n, log(1e-3), log(200)))
t <- exp(runif(n, log(1e-4), log(50))) * a^2 * exp(runif(n, -1, 1))
w <- plogis(runif(n, -12, 12))
v <- sign(runif(n, -1, 1)) * exp(runif(n, log(1e-3), log(60)))
c_ <- cert(t, v, a, w)
r <- router(t, v, a, w)
cat("points:", n, " certified:", sum(c_), " lowest router value among them:",
    round(min(r[c_]), 2), " violations (certified, router <= -600):", sum(c_ & r <= -600), "\n")

# Trials like those of a DDM-4 posterior (multilevel.R's ddm4 truth).
ty <- data.frame(t = runif(1e5, 0.05, 2.5), v = rnorm(1e5, 1.4, 0.6) * sample(c(-1, 1), 1e5, TRUE),
                 a = log1p(exp(rnorm(1e5, 1.1, 0.4))), w = plogis(rnorm(1e5, 0, 0.3)))
reach <- ty$t > 0.01 * ty$a^2
cat("typical trials that reach the router today and are certified:",
    mean(with(ty[reach, ], cert(t, v, a, w))), "\n")
