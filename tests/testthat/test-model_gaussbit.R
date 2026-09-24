context("Gaussian-probit - joint Gaussian RT and probit choice")

# Reference density, written out from the bivariate Normal it is built on: the
# RT is Normal(mu, sigma), the latent w is Normal(mudec, 1), the two correlate
# r = tanh(rho), and dec = 1 when w > 0. Given the standardized RT u, w is
# Normal(mudec + r u, 1 - r^2). The family evaluates the same probability as
# Phi(mudec cosh(rho) + u sinh(rho)); this keeps the textbook form, division
# and square root included, so the two are independent derivations.
#
# `response` is a scalar here, so the tail is picked with `if`, not ifelse() -
# which returns the length of its test and would recycle the first time's
# probability over every other.
ref_ldens <- function(y, mu, sigma, mudec, rho, response) {
  r <- tanh(rho)
  u <- (y - mu) / sigma
  z <- (mudec + r * u) / sqrt(1 - r^2)
  stats::dnorm(y, mu, sigma, log = TRUE) +
    stats::pnorm(z, lower.tail = response == 1, log.p = TRUE)
}

make_prep <- function(y, dec, mu, sigma, mudec, rho, n_draws = 10) {
  structure(
    list(
      data = list(Y = y, dec = dec),
      family = cogmod_gaussbit(),
      dpars = list(
        mu = rep(mu, n_draws), sigma = rep(sigma, n_draws),
        mudec = rep(mudec, n_draws), rho = rep(rho, n_draws)
      )
    ),
    class = "brmsprep"
  )
}

gp_pars <- list(
  list(mu = 0.6, sigma = 0.15, mudec = -1.3, rho = 0.4),
  list(mu = 0.5, sigma = 0.3, mudec = 0.7, rho = -1.2),
  list(mu = 1, sigma = 0.1, mudec = 2.5, rho = 2.5)
)


# dcogmod_gaussbit ------------------------------------------------------

test_that("dcogmod_gaussbit matches the bivariate Normal construction", {
  y <- c(0.1, 0.45, 0.6, 0.9, 1.6)
  for (p in gp_pars) {
    for (k in 0:1) {
      expect_equal(
        do.call(dcogmod_gaussbit, c(list(x = y, response = k, log = TRUE), p)),
        do.call(ref_ldens, c(list(y = y, response = k), p)),
        # the reference's sqrt(1 - r^2) loses a few digits at rho = 2.5
        tolerance = 1e-10
      )
    }
  }
})


test_that("the density integrates to one and both marginals are exact", {
  # The two properties the family's whole case rests on: whatever rho is, the
  # RT is exactly Normal(mu, sigma) and P(dec = 1) is exactly pnorm(mudec).
  for (p in gp_pars) {
    f <- function(t, k) do.call(dcogmod_gaussbit, c(list(x = t, response = k), p))
    lo <- p$mu - 12 * p$sigma
    hi <- p$mu + 12 * p$sigma
    i0 <- stats::integrate(f, lo, hi, k = 0, rel.tol = 1e-10)$value
    i1 <- stats::integrate(f, lo, hi, k = 1, rel.tol = 1e-10)$value
    expect_equal(i0 + i1, 1, tolerance = 1e-9)
    expect_equal(i1, stats::pnorm(p$mudec), tolerance = 1e-9)
    tt <- p$mu + p$sigma * c(-2, -0.3, 0, 1.1, 3)
    expect_equal(f(tt, 0) + f(tt, 1), stats::dnorm(tt, p$mu, p$sigma),
                 tolerance = 1e-12)
  }
})


test_that("rho = 0 is exactly a Gaussian times a probit", {
  # The identity behind "rho = 0 is the default analysis": the same likelihood
  # as gaussian() on the RT plus bernoulli("probit") on the choice.
  y <- c(0.2, 0.55, 1.3)
  for (k in 0:1) {
    expect_equal(
      dcogmod_gaussbit(y, mu = 0.6, sigma = 0.2, mudec = -0.8, rho = 0,
                       response = k, log = TRUE),
      stats::dnorm(y, 0.6, 0.2, log = TRUE) +
        stats::pnorm(if (k == 1) -0.8 else 0.8, log.p = TRUE),
      tolerance = 1e-14
    )
  }
})


test_that("recoding dec flips mudec and rho and nothing else", {
  y <- c(0.3, 0.6, 1.1)
  for (k in 0:1) {
    expect_equal(
      dcogmod_gaussbit(y, 0.6, 0.2, -1, 0.7, response = k, log = TRUE),
      dcogmod_gaussbit(y, 0.6, 0.2, 1, -0.7, response = 1 - k, log = TRUE),
      tolerance = 1e-14
    )
  }
})


test_that("rho > 0 makes slow trials more likely to be dec = 1", {
  # With dec coding errors, that is slow errors - the sign the docs promise.
  cond1 <- function(t, rho) {
    exp(dcogmod_gaussbit(t, 0.6, 0.15, -1.3, rho, response = 1, log = TRUE) -
          stats::dnorm(t, 0.6, 0.15, log = TRUE))
  }
  expect_gt(cond1(0.9, 0.5), cond1(0.3, 0.5))
  expect_lt(cond1(0.9, -0.5), cond1(0.3, -0.5))
  expect_equal(cond1(0.9, 0), cond1(0.3, 0), tolerance = 1e-14)
})


test_that("the far tails stay finite", {
  # A 12 s trial at a large rho puts the probit argument ~ 540 units out, which
  # cogmod_log_Phi() and pnorm(log.p = TRUE) both carry.
  for (k in 0:1) {
    expect_true(is.finite(dcogmod_gaussbit(12, response = k, rho = 1.5, log = TRUE)))
    expect_true(is.finite(dcogmod_gaussbit(-3, response = k, rho = -1.5, log = TRUE)))
  }
})


test_that("missing, infinite and empty input", {
  expect_equal(
    dcogmod_gaussbit(c(NA, Inf, -Inf), response = c(1, 1, 0), log = TRUE),
    rep(-Inf, 3)
  )
  expect_equal(dcogmod_gaussbit(c(NA, Inf), response = 0), c(0, 0))
  expect_length(dcogmod_gaussbit(numeric(0), response = numeric(0)), 0)
  expect_equal(nrow(rcogmod_gaussbit(0)), 0)
})


test_that("invalid parameters are rejected", {
  expect_error(dcogmod_gaussbit(0.5, sigma = 0, response = 1), "sigma")
  expect_error(dcogmod_gaussbit(0.5, response = 2), "response")
  expect_error(dcogmod_gaussbit(0.5), "response")
  expect_error(rcogmod_gaussbit(10, sigma = -1), "sigma")
  # the locations are unbounded
  expect_true(is.finite(dcogmod_gaussbit(0.5, mu = -2, mudec = -9, rho = -4,
                                         response = 0, log = TRUE)))
})


test_that("dcogmod_gaussbit is vectorized over every argument", {
  args <- list(x = c(0.4, 0.7), mu = c(0.5, 0.6), sigma = c(0.1, 0.2),
               mudec = c(-1, 0.5), rho = c(0.3, -0.2), response = c(0, 1))
  vec <- do.call(dcogmod_gaussbit, args)
  one <- vapply(1:2, function(i) {
    do.call(dcogmod_gaussbit, lapply(args, `[`, i))
  }, numeric(1))
  expect_equal(vec, one)
})


# rcogmod_gaussbit ------------------------------------------------------

test_that("rcogmod_gaussbit reproduces its own density", {
  set.seed(11)
  d <- rcogmod_gaussbit(2e5, mu = 0.6, sigma = 0.15, mudec = -1.3, rho = 0.4)
  expect_named(d, c("rt", "response"))
  expect_setequal(unique(d$response), c(0, 1))
  expect_equal(mean(d$response), stats::pnorm(-1.3), tolerance = 0.02)
  expect_equal(mean(d$rt), 0.6, tolerance = 0.002)
  expect_equal(stats::sd(d$rt), 0.15, tolerance = 0.005)
  # the conditional mean RT of each response, by integrating the density
  for (k in 0:1) {
    f <- function(t) dcogmod_gaussbit(t, 0.6, 0.15, -1.3, 0.4, response = k)
    e <- stats::integrate(function(t) t * f(t), -1, 2)$value /
      stats::integrate(f, -1, 2)$value
    expect_equal(mean(d$rt[d$response == k]), e, tolerance = 0.005)
  }
})


# brms methods ------------------------------------------------------------

test_that("log_lik_cogmod_gaussbit matches dcogmod_gaussbit", {
  prep <- make_prep(y = c(0.5, 0.9), dec = c(0, 1), mu = 0.6, sigma = 0.15,
                    mudec = -1.3, rho = 0.4)
  for (i in 1:2) {
    expect_equal(
      log_lik_cogmod_gaussbit(i, prep),
      rep(dcogmod_gaussbit(prep$data$Y[i], 0.6, 0.15, -1.3, 0.4,
                           response = prep$data$dec[i], log = TRUE), 10)
    )
  }
  prep$data$Y[1] <- NA
  expect_true(is.na(log_lik_cogmod_gaussbit(1, prep)))
  prep$data$dec <- NULL
  expect_error(log_lik_cogmod_gaussbit(2, prep), "dec")
})


test_that("posterior_predict_cogmod_gaussbit simulates both columns jointly", {
  set.seed(12)
  prep <- make_prep(y = 0.5, dec = 0, mu = 0.6, sigma = 0.15, mudec = -1.3,
                    rho = 0.6, n_draws = 4e4)
  out <- posterior_predict_cogmod_gaussbit(1, prep)
  expect_equal(dim(out), c(4e4, 2))
  expect_equal(colnames(out), c("rt", "response"))
  expect_equal(mean(out[, "response"]), stats::pnorm(-1.3), tolerance = 0.05)
  # Joint, not conditional on the observed RT: the simulated errors are slow.
  expect_gt(mean(out[out[, 2] == 1, 1]), mean(out[out[, 2] == 0, 1]) + 0.05)
})


test_that("posterior_epred_cogmod_gaussbit is the mean RT", {
  prep <- structure(list(dpars = list(mu = matrix(c(0.5, 0.7), 1))),
                    class = "brmsprep")
  expect_equal(posterior_epred_cogmod_gaussbit(prep), matrix(c(0.5, 0.7), 1))
})


# family ------------------------------------------------------------------

test_that("cogmod_gaussbit() builds a valid brms custom family", {
  fam <- cogmod_gaussbit()
  expect_s3_class(fam, "customfamily")
  # Order matters: brms passes the dpars to cogmod_gaussbit_lpdf in this one.
  expect_equal(fam$dpars, c("mu", "sigma", "mudec", "rho"))
  expect_equal(unname(cogmod:::.family_links(fam)),
               c("identity", "log", "identity", "identity"))
  expect_equal(fam$vars, "dec[n]")
  expect_true("cogmod_gaussbit" %in% cogmod:::.cogmod_families())
  expect_false("cogmod_gaussbit" %in% cogmod:::.OUTLIER_FAMILIES)
  expect_equal(cogmod:::.checkdata_class("cogmod_gaussbit"), "choice")
})


test_that("cens() is refused, the choice being modelled already", {
  f <- brms::bf(RT | cens(Error) ~ 1, family = cogmod_gaussbit())
  expect_error(cogmod_stanvars(f), "already models the errors")
})


# Stan code ---------------------------------------------------------------

test_that("stanvars carry the likelihood through cogmod_log_Phi()", {
  code <- cogmod_gaussbit_stanvars()[[1]]$scode
  expect_true(grepl("real cogmod_gaussbit_lpdf", code))
  expect_true(grepl("int dec", code))
  expect_true(grepl("real cogmod_log_Phi", code))
  expect_true(grepl("mudec * cosh(rho) + u * sinh(rho)", code, fixed = TRUE))
  # No normal tail other than cogmod_log_Phi(), comments aside (they name the
  # built-ins to say why they are not used).
  body <- gsub("cogmod_log_Phi", "", gsub("//[^\n]*", "", code))
  expect_false(grepl("normal_lcdf|normal_lccdf|Phi", body))
})


test_that("Stan cogmod_gaussbit_lpdf matches dcogmod_gaussbit", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  lpdf <- stan_fun("cogmod_gaussbit")
  grid <- covering_grid(
    Y = c(-0.2, 0.05, 0.4, 0.6, 1.2, 4, 12),
    mu = c(0.3, 0.6, 1.5),
    sigma = c(0.02, 0.15, 0.6),
    mudec = c(-3, -1.3, 0, 2),
    rho = c(-2.5, -0.4, 0, 0.3, 1.5),
    dec = 0:1,
    # the far tail at a strong correlation sends the probit argument through
    # cogmod_log_Phi()'s asymptotic branch, so that slice is swept in full
    always = function(g) g$sigma == 0.02 & g$mu == 0.6
  )
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    stan <- lpdf(g$Y, g$mu, g$sigma, g$mudec, g$rho, as.integer(g$dec))
    r <- dcogmod_gaussbit(g$Y, g$mu, g$sigma, g$mudec, g$rho,
                          response = g$dec, log = TRUE)
    # Relative: the far-tail log-densities run to 1e5 and beyond.
    expect_lt(abs(stan - r) / max(1, abs(r)), 1e-12)
  }

  expect_equal(lpdf(0.5, 0.6, 0, -1, 0, 0L), -Inf)
  expect_equal(lpdf(0.5, 0.6, 0.15, -1, 0, 2L), -Inf)
})


# priors and inits --------------------------------------------------------

test_that("at rho = 0 the priors are those of gaussian() + bernoulli('probit')", {
  set.seed(13)
  sim <- rcogmod_gaussbit(200, rho = 0.3)
  d <- data.frame(RT = sim$rt, Error = sim$response,
                  Condition = rep(c("a", "b"), length.out = 200))
  f_gp <- brms::bf(RT | dec(Error) ~ Condition, mudec ~ Condition, rho = 0,
                   family = cogmod_gaussbit())
  f_mv <- brms::bf(RT ~ Condition) +
    brms::bf(Error ~ Condition, family = brms::bernoulli("probit")) +
    brms::set_rescor(FALSE)
  gp <- cogmod_priors(f_gp, d)
  mv <- brms::get_prior(f_mv, data = d)

  # Rename the multivariate model's rows into the family's terms: the RT
  # response is the family's own linear predictor, the choice is `mudec`.
  mv$dpar[mv$resp == "Error"] <- "mudec"
  key <- function(p) {
    p <- as.data.frame(p)[, c("prior", "class", "coef", "dpar")]
    p <- p[order(p$dpar, p$class, p$coef), ]
    rownames(p) <- NULL
    p
  }
  expect_equal(key(gp), key(mv))
})


test_that("cogmod_priors adds a prior on rho, on one scale in both forms", {
  sim <- rcogmod_gaussbit(200, rho = 0.3)
  d <- data.frame(RT = sim$rt, Error = sim$response, x = rnorm(200))
  modelled <- as.data.frame(cogmod_priors(
    brms::bf(RT | dec(Error) ~ 1, rho ~ x, family = cogmod_gaussbit()), d
  ))
  omitted <- as.data.frame(cogmod_priors(
    brms::bf(RT | dec(Error) ~ 1, family = cogmod_gaussbit()), d
  ))
  expect_equal(modelled$prior[modelled$dpar == "rho" & modelled$class == "Intercept"],
               "normal(0, 0.5)")
  expect_true(all(modelled$prior[modelled$dpar == "rho" & modelled$class == "b" &
                                   !nzchar(modelled$coef)] == "normal(0, 0.5)"))
  expect_equal(omitted$prior[omitted$class == "rho"], "normal(0, 0.5)")
  expect_equal(omitted$prior[omitted$class == "mudec"], "student_t(3, 0, 2.5)")
})


test_that("cogmod_inits covers the declared parameters", {
  sim <- rcogmod_gaussbit(150)
  d <- data.frame(RT = sim$rt, Error = sim$response,
                  Condition = rep(c("a", "b"), length.out = 150))
  f <- brms::bf(RT | dec(Error) ~ Condition, sigma ~ 1, mudec ~ Condition,
                rho ~ 1, family = cogmod_gaussbit())
  vals <- cogmod_inits(f, d)(1)
  code <- brms::make_stancode(f, data = d, stanvars = cogmod_stanvars(f))
  declared <- vapply(cogmod:::.stan_param_decls(code), `[[`, character(1), "name")
  expect_setequal(names(vals), declared)

  v0 <- cogmod_inits(f, d, jitter = 0)(1)
  expect_equal(v0$Intercept, 0.6, tolerance = 1e-8)
  expect_equal(exp(v0$Intercept_sigma), 0.2, tolerance = 1e-8)
  expect_equal(v0$Intercept_rho, 0, tolerance = 1e-8)
})
