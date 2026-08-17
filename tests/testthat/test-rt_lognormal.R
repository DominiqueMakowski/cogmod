context("Shifted LogNormal - brms")

make_prep <- function(y, mu, sigma, ndt, poutlier, n_draws = 10) {
  structure(
    list(
      data = list(Y = y),
      dpars = list(
        mu = rep(mu, n_draws),
        sigma = rep(sigma, n_draws),
        ndt = rep(ndt, n_draws),
        poutlier = rep(poutlier, n_draws)
      )
    ),
    class = "brmsprep"
  )
}

# Reference implementation of the mixture, written out longhand
ref_dens <- function(y, mu, sigma, ndt, poutlier) {
  d_out <- dlnorm(y, log(0.15), 1.5)
  d_dec <- if (y > ndt) dlnorm(y - ndt, mu, sigma) else 0
  poutlier * d_out + (1 - poutlier) * d_dec
}


# drt_lognormal -----------------------------------------------------------

test_that("drt_lognormal matches the mixture density", {
  mu <- -0.7
  sigma <- 0.5
  ndt <- 0.3
  poutlier <- 0.02

  for (y in c(0.05, 0.31, 0.5, 0.9, 1.5, 3)) {
    expect_equal(
      drt_lognormal(y, mu, sigma, ndt, poutlier),
      ref_dens(y, mu, sigma, ndt, poutlier),
      tolerance = 1e-12,
      label = sprintf("density at y = %.2f", y)
    )
    expect_equal(
      drt_lognormal(y, mu, sigma, ndt, poutlier, log = TRUE),
      log(ref_dens(y, mu, sigma, ndt, poutlier)),
      tolerance = 1e-12,
      label = sprintf("log-density at y = %.2f", y)
    )
  }
})

test_that("responses faster than ndt keep positive density (the whole point)", {
  # Under the old tau * minrt parameterization these were 0, which is what
  # forced ndt below the fastest observed RT.
  d <- drt_lognormal(0.1, ndt = 0.4, poutlier = 0.02)
  expect_gt(d, 0)
  expect_equal(d, 0.02 * dlnorm(0.1, log(0.15), 1.5), tolerance = 1e-12)

  # With no outlier component they are impossible again
  expect_equal(drt_lognormal(0.1, ndt = 0.4, poutlier = 0), 0)
  expect_true(
    is.infinite(drt_lognormal(0.1, ndt = 0.4, poutlier = 0, log = TRUE))
  )
})

test_that("poutlier = 0 recovers the plain shifted LogNormal", {
  expect_equal(
    drt_lognormal(0.9, mu = -0.7, sigma = 0.5, ndt = 0.3, poutlier = 0),
    dlnorm(0.6, -0.7, 0.5),
    tolerance = 1e-12
  )
})

test_that("drt_lognormal is vectorized and integrates to one", {
  d <- drt_lognormal(c(0.5, 0.9), mu = c(-0.7, -0.5), sigma = 0.5, ndt = 0.3)
  expect_length(d, 2)
  expect_equal(
    d,
    c(dlnorm(0.2, -0.7, 0.5), dlnorm(0.6, -0.5, 0.5)),
    tolerance = 1e-12
  )

  total <- integrate(
    function(z) drt_lognormal(z, ndt = 0.3, poutlier = 0.05),
    lower = 0, upper = Inf
  )$value
  expect_equal(total, 1, tolerance = 1e-4)
})

test_that("drt_lognormal returns 0 density for invalid parameters", {
  expect_warning(d <- drt_lognormal(0.5, sigma = -0.5))
  expect_equal(d, 0)
  expect_warning(d <- drt_lognormal(0.5, ndt = -0.3))
  expect_equal(d, 0)
  expect_warning(d <- drt_lognormal(0.5, poutlier = 1.5))
  expect_equal(d, 0)
})


# rrt_lognormal -----------------------------------------------------------

test_that("rrt_lognormal recovers the mixture mean", {
  set.seed(123)
  mu <- -0.7
  sigma <- 0.5
  ndt <- 0.3
  poutlier <- 0.05

  rts <- rrt_lognormal(2e5, mu, sigma, ndt, poutlier)
  theo <- (1 - poutlier) * (exp(mu + sigma^2 / 2) + ndt) +
    poutlier * exp(log(0.15) + 1.5^2 / 2)

  expect_equal(mean(rts), theo, tolerance = 0.05)
  expect_true(all(rts > 0))
  expect_length(rts, 2e5)
})

test_that("rrt_lognormal produces responses below ndt only via outliers", {
  set.seed(1)
  expect_true(all(rrt_lognormal(5000, ndt = 0.3, poutlier = 0) > 0.3))
  expect_true(any(rrt_lognormal(5000, ndt = 0.3, poutlier = 0.5) < 0.3))
})

test_that("rrt_lognormal errors on invalid parameters", {
  expect_error(rrt_lognormal(10, sigma = -1), "sigma")
  expect_error(rrt_lognormal(10, ndt = -1), "ndt")
  expect_error(rrt_lognormal(10, poutlier = 2), "poutlier")
  expect_error(rrt_lognormal(-5), "non-negative integer")
})


# brms methods ------------------------------------------------------------

test_that("log_lik_rt_lognormal matches drt_lognormal", {
  for (y in c(0.1, 0.35, 0.9, 2)) {
    prep <- make_prep(y, mu = -0.7, sigma = 0.5, ndt = 0.3, poutlier = 0.02)
    expect_equal(
      log_lik_rt_lognormal(1, prep),
      rep(log(ref_dens(y, -0.7, 0.5, 0.3, 0.02)), 10),
      tolerance = 1e-10,
      label = sprintf("log_lik at y = %.2f", y)
    )
  }
})

test_that("log_lik_rt_lognormal returns -Inf for invalid parameters", {
  prep <- structure(
    list(
      data = list(Y = 0.5),
      dpars = list(
        mu = c(-0.7, -0.7, -0.7),
        sigma = c(-0.5, 0.5, 0.5),
        ndt = c(0.3, -0.3, 0.3),
        poutlier = c(0.02, 0.02, 1.5)
      )
    ),
    class = "brmsprep"
  )
  expect_true(all(suppressWarnings(log_lik_rt_lognormal(1, prep)) == -Inf))
})

test_that("posterior_predict_rt_lognormal recovers the mixture mean", {
  set.seed(123)
  n_draws <- 1e5
  mu <- -0.7
  sigma <- 0.5
  ndt <- 0.3
  poutlier <- 0.05

  prep <- structure(
    list(dpars = list(
      mu = rep(mu, n_draws), sigma = rep(sigma, n_draws),
      ndt = rep(ndt, n_draws), poutlier = rep(poutlier, n_draws)
    )),
    class = "brmsprep"
  )

  rts <- posterior_predict_rt_lognormal(1, prep)
  theo <- (1 - poutlier) * (exp(mu + sigma^2 / 2) + ndt) +
    poutlier * exp(log(0.15) + 1.5^2 / 2)

  expect_equal(mean(rts), theo, tolerance = 0.05)
  expect_true(all(rts > 0))
})

test_that("posterior_epred_rt_lognormal equals the mixture of component means", {
  mu <- matrix(c(-0.7, -0.5), nrow = 2, ncol = 2)
  sigma <- matrix(0.5, nrow = 2, ncol = 2)
  ndt <- matrix(0.3, nrow = 2, ncol = 2)
  poutlier <- matrix(0.02, nrow = 2, ncol = 2)

  prep <- structure(
    list(dpars = list(mu = mu, sigma = sigma, ndt = ndt, poutlier = poutlier)),
    class = "brmsprep"
  )

  expect_equal(
    posterior_epred_rt_lognormal(prep),
    0.98 * (exp(mu + sigma^2 / 2) + ndt) + 0.02 * exp(log(0.15) + 1.5^2 / 2)
  )
})


# family ------------------------------------------------------------------

test_that("rt_lognormal() builds a valid brms custom family", {
  fam <- rt_lognormal()

  expect_s3_class(fam, "customfamily")
  expect_identical(fam$dpars, c("mu", "sigma", "ndt", "poutlier"))
  expect_identical(
    unname(c(fam$link, fam$link_sigma, fam$link_ndt, fam$link_poutlier)),
    c("identity", "softplus", "log", "logit")
  )
  # tau and minrt are gone
  expect_false(any(c("tau", "minrt") %in% fam$dpars))
})

test_that("stanvars carry the likelihood with the hard-coded outlier component", {
  sv <- rt_lognormal_stanvars()
  code <- paste(vapply(sv, function(x) x$scode, character(1)), collapse = "\n")

  expect_true(grepl("rt_lognormal_lpdf", code, fixed = TRUE))
  expect_true(grepl("log_mix", code, fixed = TRUE))
  expect_true(grepl("log(0.15), 1.5", code, fixed = TRUE))
})
