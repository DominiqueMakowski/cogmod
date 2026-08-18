context("RDM")

test_that("rcogmod_rdm behaves correctly under various conditions", {
  n <- 10000 # For proportion checks
  set.seed(456) # Use a single seed for the combined test

  # --- Parameters ---
  boundary_val <- 0.5
  bias_val <- 0.2
  ndt_val <- 0.15

  # --- Test 1: Basic Structure and Properties ---
  sim_data_t1 <- rcogmod_rdm(n = 100, vzero = 0.8, vone = 0.6, boundary = boundary_val, bias = bias_val, ndt = ndt_val)

  # Check structure
  expect_named(sim_data_t1, c("rt", "response"), label = "Output has correct names")
  expect_equal(nrow(sim_data_t1), 100, label = "Output has correct number of rows")

  # Check values
  expect_true(all(sim_data_t1$rt >= ndt_val), label = "All RTs >= NDT")
  expect_true(all(is.finite(sim_data_t1$rt)), label = "All RTs are finite")
  expect_true(all(sim_data_t1$response %in% c(0L, 1L)), label = "Responses are 0 or 1")

  # --- Test 2: Choice Proportions Align with Drift Rates ---

  # Case 2.1: vzero > vone
  sim_21 <- rcogmod_rdm(n = n, vzero = 1.0, vone = 0.7, boundary = boundary_val, bias = bias_val, ndt = ndt_val)
  expect_gt(mean(sim_21$response == 0), 0.5, label = "vzero > vone: expect P(response=0) > 0.5")

  # Case 2.2: vone > vzero
  sim_22 <- rcogmod_rdm(n = n, vzero = 0.6, vone = 0.9, boundary = boundary_val, bias = bias_val, ndt = ndt_val)
  expect_gt(mean(sim_22$response == 1), 0.5, label = "vone > vzero: expect P(response=1) > 0.5")

  # Case 2.3: a zero drift is slow, not absent. Driftless Brownian motion still
  # reaches the threshold with probability one, so that accumulator sometimes
  # wins -- and the RTs stay finite.
  sim_23 <- rcogmod_rdm(n = n, vzero = 0.8, vone = 0, boundary = boundary_val, bias = bias_val, ndt = ndt_val)
  expect_gt(mean(sim_23$response == 0), 0.5, label = "vone = 0: vzero usually wins")
  expect_gt(mean(sim_23$response == 1), 0, label = "vone = 0: still wins occasionally")
  expect_true(all(sim_23$rt >= ndt_val), label = "All RTs >= NDT (vone = 0)")
  expect_true(all(is.finite(sim_23$rt)), label = "All RTs are finite (vone = 0)")
})


test_that("dcogmod_rdm matches the simulated distribution and log=TRUE works", {
  set.seed(2025)
  pars <- list(vzero = 0.9, vone = 0.7, boundary = 0.4, bias = 0.15, ndt = 0.2)

  dat <- do.call(rcogmod_rdm, c(list(n = 20000), pars))

  # Marginal density against a kernel estimate of the simulated RTs
  emp <- density(dat$rt, bw = "SJ", n = 64, from = pars$ndt + 1e-8)
  keep <- emp$x > pars$ndt + 1e-3
  y_th <- do.call(dcogmod_rdm, c(list(x = emp$x[keep]), pars))
  expect_gt(cor(emp$y[keep], y_th), 0.99)

  # log argument
  pts <- c(0.25, 0.6, 1.0)
  log_d <- do.call(dcogmod_rdm, c(list(x = pts, log = TRUE), pars))
  raw_d <- do.call(dcogmod_rdm, c(list(x = pts, log = FALSE), pars))
  expect_equal(log_d, log(raw_d), tolerance = 1e-8)
})


test_that("dcogmod_rdm defective densities decompose the marginal and match simulation", {
  set.seed(11)
  pars <- list(vzero = 0.9, vone = 0.7, boundary = 0.4, bias = 0.15, ndt = 0.2)
  xs <- c(0.25, 0.4, 0.8, 1.5, 3)

  marg <- do.call(dcogmod_rdm, c(list(x = xs), pars))
  d0 <- do.call(dcogmod_rdm, c(list(x = xs, response = 0), pars))
  d1 <- do.call(dcogmod_rdm, c(list(x = xs, response = 1), pars))

  # The per-response densities are defective and must partition the marginal
  expect_equal(d0 + d1, marg, tolerance = 1e-12)

  # Each integrates to that response's probability, and they sum to one
  i0 <- integrate(function(z) do.call(dcogmod_rdm, c(list(x = z, response = 0), pars)),
                  0, Inf, rel.tol = 1e-9, subdivisions = 2000)$value
  i1 <- integrate(function(z) do.call(dcogmod_rdm, c(list(x = z, response = 1), pars)),
                  0, Inf, rel.tol = 1e-9, subdivisions = 2000)$value
  expect_equal(i0 + i1, 1, tolerance = 1e-6)

  sim <- do.call(rcogmod_rdm, c(list(n = 50000), pars))
  expect_equal(mean(sim$response == 0), i0, tolerance = 0.01)

  expect_error(do.call(dcogmod_rdm, c(list(x = 0.5, response = 2), pars)))
})


test_that("pcogmod_rdm agrees with integrating dcogmod_rdm and with simulation", {
  set.seed(99)
  pars <- list(vzero = 0.9, vone = 0.7, boundary = 0.4, bias = 0.15, ndt = 0.2)
  sim <- do.call(rcogmod_rdm, c(list(n = 50000), pars))

  for (q in c(0.3, 0.5, 1.0, 2.0)) {
    surv <- do.call(pcogmod_rdm, c(list(q = q, lower.tail = FALSE), pars))
    num <- integrate(function(z) do.call(dcogmod_rdm, c(list(x = z), pars)),
                     q, Inf, rel.tol = 1e-10, subdivisions = 2000)$value
    expect_equal(surv, num, tolerance = 1e-7,
                 label = sprintf("pcogmod_rdm survival vs integrated density at q=%.1f", q))
    expect_equal(surv, mean(sim$rt > q), tolerance = 0.02,
                 label = sprintf("pcogmod_rdm survival vs empirical at q=%.1f", q))
  }

  # CDF and survival are complementary, and log.p is consistent
  q <- c(0.3, 0.8, 2.0)
  lo <- do.call(pcogmod_rdm, c(list(q = q), pars))
  hi <- do.call(pcogmod_rdm, c(list(q = q, lower.tail = FALSE), pars))
  expect_equal(lo + hi, rep(1, length(q)), tolerance = 1e-12)
  expect_equal(do.call(pcogmod_rdm, c(list(q = q, log.p = TRUE), pars)), log(lo),
               tolerance = 1e-10)
})


test_that("the Wald survival stays finite where 1 - CDF would underflow", {
  # This is the regression that motivated the log-space rewrite. Computing the
  # survival as 1 - CDF returns exactly 0 here, which would give the sampler a
  # -Inf log-likelihood with a zero gradient.
  n1 <- function(v) rep(v, 1)
  for (drift in c(3, 6, 10)) {
    for (t in c(4, 8, 20)) {
      ls <- .swald(t, n1(drift), n1(0.5), n1(0.2), n1(0), log.p = TRUE)
      expect_true(is.finite(ls),
                  label = sprintf("log survival finite at drift=%g, t=%g", drift, t))
      expect_lt(ls, 0)
    }
  }
  # ... and it is still the right number, matching numerical integration.
  dens <- function(z) .dwald(z, rep(3, length(z)), rep(0.5, length(z)),
                             rep(0.2, length(z)), rep(0, length(z)))
  num <- integrate(dens, 3, Inf, rel.tol = 1e-11, subdivisions = 2000)$value
  expect_equal(.swald(3, 3, 0.5, 0.2, 0), num, tolerance = 1e-8)
})


test_that("the Wald density and CDF survive extreme drift x threshold", {
  # exp(2 * drift * b) overflows on its own here; folding it into a single
  # exponent with the log-CDF keeps the product representable.
  for (prm in list(c(80, 4, 1), c(200, 1.6, 0.4), c(500, 4, 1))) {
    ld <- .dwald(0.5, prm[1], prm[2], prm[3], 0.1, log = TRUE)
    ls <- .swald(0.5, prm[1], prm[2], prm[3], 0.1, log.p = TRUE)
    expect_true(is.finite(ld) && !is.nan(ld))
    expect_true(is.finite(ls) && !is.nan(ls))
  }
})


test_that("the driftless Wald CDF is correct (regression: was off by a factor of 2)", {
  # `.pwald()`'s drift == 0 branch used to return half the right answer. It was
  # unreachable from dcogmod_rdm() at the time, but the Stan version reaches it
  # whenever the drift approaches zero.
  k <- 0.5
  A <- 0.2
  for (t in c(0.3, 1.0, 3.0)) {
    num <- integrate(function(z) .dwald(z, rep(0, length(z)), rep(k, length(z)),
                                        rep(A, length(z)), rep(0, length(z))),
                     0, t, rel.tol = 1e-11, subdivisions = 2000)$value
    expect_equal(.pwald(t, 0, k, A, 0), num, tolerance = 1e-8,
                 label = sprintf("driftless CDF at t=%g", t))
  }
  # And it is continuous as the drift crosses zero
  s <- vapply(c(1e-5, 1e-7, 1e-9, 0, -1e-9), function(v) {
    .swald(1.1, v, 1, 0.3, 0.1, log.p = TRUE)
  }, numeric(1))
  expect_lt(max(abs(diff(s))), 1e-5)
})


test_that("the Wald density integrates to one across the parameter grid", {
  for (prm in list(c(3, 0.5, 0.2), c(1, 1, 0.5), c(0.5, 2, 0.1),
                   c(6, 0.3, 0.05), c(0, 1, 0.3))) {
    nu <- prm[1]; k <- prm[2]; A <- prm[3]
    tot <- integrate(function(z) .dwald(z, rep(nu, length(z)), rep(k, length(z)),
                                        rep(A, length(z)), rep(0, length(z))),
                     0, Inf, rel.tol = 1e-10, subdivisions = 3000)$value
    expect_equal(tot, 1, tolerance = 1e-6,
                 label = sprintf("total mass at drift=%g, boundary=%g, bias=%g", nu, k, A))
  }
})


test_that("dcogmod_rdm edge-cases x<=ndt give zero density", {
  p <- list(vzero = 0.8, vone = 0.6, boundary = 0.5, bias = 0.2, ndt = 0.15)
  xs <- c(p$ndt - 0.01, p$ndt, p$ndt + 0.05)

  dens <- dcogmod_rdm(xs, vzero = p$vzero, vone = p$vone, boundary = p$boundary, bias = p$bias, ndt = p$ndt)

  expect_equal(dens[1:2], c(0, 0), label = "density = 0 for x <= ndt")
  expect_gt(dens[3], 0, label = "density > 0 for x > ndt")

  # Immediately above the ndt the density really is below the smallest double,
  # so the natural scale is 0 -- but the log density stays finite and usable,
  # which is the property the sampler depends on. (The previous implementation
  # substituted .Machine$double.xmin here, which is a fabricated value.)
  lg <- dcogmod_rdm(p$ndt + c(1e-8, 1e-6, 1e-4), vzero = p$vzero, vone = p$vone,
             boundary = p$boundary, bias = p$bias, ndt = p$ndt, log = TRUE)
  expect_true(all(is.finite(lg)), label = "log density finite just above the ndt")
  expect_true(all(diff(lg) > 0), label = "log density increases away from the ndt")

  # dcogmod_rdm errors on invalid inputs
  expect_error(dcogmod_rdm(0.5, vzero = -0.1, vone = 0.5, boundary = 0.5, bias = 0.2, ndt = 0.1))
  expect_error(dcogmod_rdm(0.5, vzero = 0.5, vone = -0.1, boundary = 0.5, bias = 0.2, ndt = 0.1))
  expect_error(dcogmod_rdm(0.5, vzero = 0.5, vone = 0.5, boundary = 0, bias = 0.2, ndt = 0.1))
  expect_error(dcogmod_rdm(0.5, vzero = 0.5, vone = 0.5, boundary = -1, bias = 0.2, ndt = 0.1))
  expect_error(dcogmod_rdm(0.5, vzero = 0.5, vone = 0.5, boundary = 0.5, bias = 0, ndt = 0.1))
  expect_error(dcogmod_rdm(0.5, vzero = 0.5, vone = 0.5, boundary = 0.5, bias = -1, ndt = 0.1))
  expect_error(dcogmod_rdm(0.5, vzero = 0.5, vone = 0.5, boundary = 0.5, bias = 0.2, ndt = -0.1))
})


context("RDM - brms")

test_that("Stan cogmod_rdm_lpdf matches R dcogmod_rdm", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  cogmod_rdm_lpdf_stan <- stan_fun("cogmod_rdm")

  grid <- expand.grid(
    Y = c(0.25, 0.4, 0.8, 1.5, 3.0, 8.0),
    mu = c(0, 0.5, 2, 6),
    driftone = c(0.3, 1.5, 5),
    sigmabias = c(0.05, 0.3, 1.0),
    boundary = c(0.2, 0.8, 2.0),
    tau = c(0.3, 0.9),
    dec = c(0, 1),
    KEEP.OUT.ATTRS = FALSE
  )
  minrt <- 0.2
  grid <- grid[grid$Y > grid$tau * minrt + 1e-9, ]

  stan_ll <- mapply(
    cogmod_rdm_lpdf_stan,
    grid$Y, grid$mu, grid$driftone, grid$sigmabias, grid$boundary,
    grid$tau, minrt, grid$dec
  )
  r_ll <- mapply(
    function(Y, mu, d1, sb, boundary, tau, dec) {
      dcogmod_rdm(x = Y, response = dec, vzero = mu, vone = d1, boundary = boundary,
           bias = sb, ndt = tau * minrt, log = TRUE)
    },
    grid$Y, grid$mu, grid$driftone, grid$sigmabias, grid$boundary, grid$tau, grid$dec
  )

  expect_equal(is.finite(stan_ll), is.finite(r_ll),
               label = "Stan and R agree on which points are finite")
  expect_true(all(is.finite(r_ll)), label = "no finite-grid point returns -Inf")
  expect_equal(stan_ll, r_ll, tolerance = 1e-6)
})


test_that("Stan cogmod_rdm_lpdf is finite in the tails and for fast responses", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  cogmod_rdm_lpdf_stan <- stan_fun("cogmod_rdm")

  # Slow responses with a fast losing accumulator: the survival term is what
  # collapses here if it is computed as 1 - CDF.
  for (prm in list(c(5, 10), c(8, 10), c(20, 6))) {
    ll <- cogmod_rdm_lpdf_stan(prm[1], prm[2], prm[2], 0.2, 0.5, 0.5, 0.2, 0)
    expect_true(is.finite(ll),
                label = sprintf("finite at Y=%g, drift=%g", prm[1], prm[2]))
  }

  # Responses just above the non-decision time, where (boundary - drift*t)/sqrt(t)
  # grows past 8 and Stan's own std_normal_lccdf would return -Inf.
  for (dy in c(1e-2, 1e-3, 1e-4)) {
    ll <- cogmod_rdm_lpdf_stan(0.1 + dy, 2, 1.5, 0.3, 1.0, 0.5, 0.2, 0)
    expect_true(is.finite(ll),
                label = sprintf("finite at %g above the ndt", dy))
  }

  # Y at or below the non-decision time is impossible
  expect_equal(cogmod_rdm_lpdf_stan(0.05, 2, 1, 0.3, 0.5, 0.9, 0.2, 0), -Inf)
  # Invalid inputs
  expect_equal(cogmod_rdm_lpdf_stan(0.5, 2, 1, 0.3, 0.5, 0.5, 0.2, 2), -Inf)
  expect_equal(cogmod_rdm_lpdf_stan(0.5, 2, 1, 0, 0.5, 0.5, 0.2, 0), -Inf)
  expect_equal(cogmod_rdm_lpdf_stan(0.5, 2, 1, 0.3, 0, 0.5, 0.2, 0), -Inf)
})


# End-to-end parameter recovery. Left commented out (as in test-model_lba2.R) because
# it takes several minutes, but run it before a release.
#
# test_that("cogmod_rdm model recovers its parameters with brms", {
#   skip_on_cran()
#   skip_if_not_installed("cmdstanr")
#
#   set.seed(2024)
#   true <- list(vzero = 2.5, vone = 1.6, boundary = 0.8, bias = 0.3,
#                tau = 0.6, minrt = 0.25)
#   df <- rcogmod_rdm(n = 4000, vzero = true$vzero, vone = true$vone, boundary = true$boundary,
#              bias = true$bias, ndt = true$tau * true$minrt)
#
#   f <- brms::bf(
#     rt | dec(response) ~ 1,
#     driftone ~ 1,
#     sigmabias ~ 1,
#     boundary ~ 1,
#     tau ~ 1,
#     minrt = true$minrt,
#     family = cogmod_rdm()
#   )
#
#   # `sigmabias` and `boundary` trade off along a very flat ridge (see ?cogmod_rdm), so a
#   # weakly informative prior is needed to keep the sampler off the
#   # sigmabias -> 0 edge. Without it this fit produces hundreds of divergences
#   # and sigmabias collapses to zero.
#   pr <- brms::prior(normal(-1, 1), class = "Intercept", dpar = "sigmabias")
#
#   fit <- brms::brm(f, data = df, stanvars = cogmod_rdm_stanvars(), prior = pr,
#                    init = 0.5, chains = 2, iter = 1000, warmup = 500,
#                    control = list(adapt_delta = 0.95), backend = "cmdstanr")
#
#   sp <- function(x) log1p(exp(x))
#   po <- brms::posterior_summary(fit)[, "Estimate"]
#   expect_equal(sp(po[["b_Intercept"]]), true$vzero, tolerance = 0.15)
#   expect_equal(sp(po[["b_driftone_Intercept"]]), true$vone, tolerance = 0.15)
#   expect_equal(plogis(po[["b_tau_Intercept"]]), true$tau, tolerance = 0.15)
#   # The threshold sum is the well-identified combination
#   expect_equal(sp(po[["b_boundary_Intercept"]]) + sp(po[["b_sigmabias_Intercept"]]),
#                true$boundary + true$bias, tolerance = 0.15)
#
#   ll <- brms::log_lik(fit, ndraws = 20)
#   expect_true(all(is.finite(ll)))
#   pp <- brms::posterior_predict(fit, ndraws = 20, newdata = df[1:5, ])
#   expect_true(all(is.finite(pp)))
# })


test_that("cogmod_rdm() family is registered consistently with the Stan signature", {
  fam <- cogmod_rdm()
  expect_equal(fam$name, "cogmod_rdm")
  # Order matters: brms passes the dpars to cogmod_rdm_lpdf in exactly this order.
  expect_equal(fam$dpars, c("mu", "driftone", "sigmabias", "boundary", "tau", "minrt"))
  # brms stores the bounds as character
  expect_equal(fam$ub[["tau"]], "1")
  expect_true(all(unlist(fam$lb) == "0"))
  expect_true(all(is.na(unlist(fam$ub)[c("mu", "driftone", "sigmabias", "boundary")])))
  expect_equal(fam$vars, "dec[n]")
})
