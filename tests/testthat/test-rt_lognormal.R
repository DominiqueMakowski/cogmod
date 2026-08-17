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
  d_out <- 2 * dt(y / 0.4, df = 3) / 0.4
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
  expect_equal(d, 0.02 * 2 * dt(0.1 / 0.4, df = 3) / 0.4, tolerance = 1e-12)

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
    poutlier * (2 * 0.4 * sqrt(3) * gamma(2) / (sqrt(pi) * 2 * gamma(1.5)))

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

test_that("posterior_predict_rt_lognormal excludes outliers by default", {
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

  # No flag on the family: the decision process alone, so nothing below ndt
  rts <- posterior_predict_rt_lognormal(1, prep)
  expect_equal(mean(rts), exp(mu + sigma^2 / 2) + ndt, tolerance = 0.02)
  expect_true(all(rts > ndt))

  # Asking for the mixture recovers the mixture mean
  rts_mix <- posterior_predict_rt_lognormal(1, prep, predict_outliers = TRUE)
  theo <- (1 - poutlier) * (exp(mu + sigma^2 / 2) + ndt) +
    poutlier * (2 * 0.4 * sqrt(3) * gamma(2) / (sqrt(pi) * 2 * gamma(1.5)))

  expect_equal(mean(rts_mix), theo, tolerance = 0.05)
  expect_true(all(rts_mix > 0))
  expect_true(any(rts_mix < ndt))
})

test_that("posterior_epred_rt_lognormal excludes outliers by default", {
  mu <- matrix(c(-0.7, -0.5), nrow = 2, ncol = 2)
  sigma <- matrix(0.5, nrow = 2, ncol = 2)
  ndt <- matrix(0.3, nrow = 2, ncol = 2)
  poutlier <- matrix(0.02, nrow = 2, ncol = 2)

  prep <- structure(
    list(dpars = list(mu = mu, sigma = sigma, ndt = ndt, poutlier = poutlier)),
    class = "brmsprep"
  )

  # Default: the decision-component mean exactly
  expect_equal(
    posterior_epred_rt_lognormal(prep),
    exp(mu + sigma^2 / 2) + ndt
  )

  # With the outlier component: the mixture of component means
  expect_equal(
    posterior_epred_rt_lognormal(prep, predict_outliers = TRUE),
    0.98 * (exp(mu + sigma^2 / 2) + ndt) + 0.02 * (2 * 0.4 * sqrt(3) * gamma(2) / (sqrt(pi) * 2 * gamma(1.5)))
  )

  # and the two differ
  expect_false(isTRUE(all.equal(
    posterior_epred_rt_lognormal(prep),
    posterior_epred_rt_lognormal(prep, predict_outliers = TRUE)
  )))
})

test_that("the family flag drives predictions when no argument is given", {
  set.seed(123)
  n_draws <- 5e4
  mu <- -0.7
  sigma <- 0.5
  ndt <- 0.3
  poutlier <- 0.2 # exaggerated so the two differ unmistakably

  dpars <- list(
    mu = rep(mu, n_draws), sigma = rep(sigma, n_draws),
    ndt = rep(ndt, n_draws), poutlier = rep(poutlier, n_draws)
  )
  on <- structure(
    list(dpars = dpars,
         family = list(name = "rt_lognormal", predict_outliers = TRUE)),
    class = "brmsprep"
  )
  off <- structure(
    list(dpars = dpars,
         family = list(name = "rt_lognormal", predict_outliers = FALSE)),
    class = "brmsprep"
  )

  expect_true(any(posterior_predict_rt_lognormal(1, on) < ndt))
  expect_true(all(posterior_predict_rt_lognormal(1, off) > ndt))

  # an explicit argument overrides the flag either way
  expect_true(all(
    posterior_predict_rt_lognormal(1, on, predict_outliers = FALSE) > ndt
  ))
})

test_that(".predict_outliers resolves the argument, then the family flag", {
  bare <- structure(list(family = list(name = "rt_lognormal")),
                    class = "brmsprep")
  on <- structure(list(family = list(name = "rt_lognormal",
                                     predict_outliers = TRUE)),
                  class = "brmsprep")

  # no flag anywhere: exclude them
  expect_false(cogmod:::.predict_outliers(NULL, bare))
  # flag on the family is honoured
  expect_true(cogmod:::.predict_outliers(NULL, on))
  # an explicit argument overrides the flag either way
  expect_false(cogmod:::.predict_outliers(FALSE, on))
  expect_true(cogmod:::.predict_outliers(TRUE, bare))
})

test_that("with_outliers sets the flag where prepare_predictions can see it", {
  # brms rebuilds prep$family from object$formula$family, so the flag must live
  # there; object$family alone is dropped. Verified against brms 2.23.1. If this
  # changes upstream, with_outliers() and its docs need revisiting.
  fake <- structure(
    list(family = list(name = "rt_lognormal"),
         formula = list(family = list(name = "rt_lognormal"))),
    class = "brmsfit"
  )
  m <- with_outliers(fake)
  expect_true(m$formula$family$predict_outliers)
  expect_true(m$family$predict_outliers)

  m2 <- without_outliers(m)
  expect_false(m2$formula$family$predict_outliers)
  expect_false(m2$family$predict_outliers)

  for (f in list(with_outliers, without_outliers)) {
    expect_error(f(list()), "must be a brmsfit")
    wrong <- structure(list(family = list(name = "lnr")), class = "brmsfit")
    expect_error(f(wrong), "rt_lognormal")
  }
})

test_that("brms does not forward ... to custom family prediction methods", {
  # This is why the flag rides on the model instead of being an argument.
  flat <- function(f) gsub("\\s+", " ", paste(deparse(body(f)), collapse = " "))

  # posterior_predict.brmsfit sends ... to prepare_predictions(), then calls the
  # brmsprep method with named arguments only
  fit_body <- flat(brms:::posterior_predict.brmsfit)
  expect_true(grepl("prepare_predictions(object", fit_body, fixed = TRUE))
  expect_false(grepl("...", sub(".*posterior_predict\\(prep", "", fit_body),
                     fixed = TRUE))

  # epred and log_lik have no ... at all in their custom dispatch
  expect_false("..." %in% names(formals(brms:::posterior_epred_custom)))
  expect_false("..." %in% names(formals(brms:::log_lik_custom)))
})


# p_outlier ---------------------------------------------------------------

test_that("p_outlier returns the mixture responsibility", {
  y <- c(0.05, 0.4, 0.9, 3.0) # below ndt / just after / bulk / far tail
  nd <- 4
  prep <- list(
    data = list(Y = y),
    dpars = list(
      mu = matrix(-0.9, nd, length(y)), sigma = matrix(0.5, nd, length(y)),
      ndt = matrix(0.25, nd, length(y)), poutlier = matrix(0.02, nd, length(y))
    )
  )

  out <- with_mocked_bindings(
    p_outlier(structure(list(), class = "brmsfit")),
    prepare_predictions = function(object, ...) prep,
    get_dpar = function(prep, dpar, ...) prep$dpars[[dpar]],
    .package = "brms"
  )

  # Reference responsibility, written out longhand
  ref <- vapply(y, function(v) {
    g <- 2 * dt(v / 0.4, df = 3) / 0.4
    f <- if (v > 0.25) dlnorm(v - 0.25, -0.9, 0.5) else 0
    0.02 * g / (0.02 * g + 0.98 * f)
  }, numeric(1))

  expect_equal(out$p_outlier, ref, tolerance = 1e-12)
  expect_equal(out$rt, y)
  expect_equal(out$fast, y < median(y))

  # The shape the article plots: certainty below ndt, near zero in the bulk,
  # and rising again in the far slow tail because the half-t has heavier tails
  # than the LogNormal.
  expect_equal(out$p_outlier[1], 1)
  expect_lt(out$p_outlier[3], out$p_outlier[2])
  expect_gt(out$p_outlier[4], out$p_outlier[3])
})

test_that("p_outlier rejects non-brmsfit input", {
  expect_error(p_outlier(list()), "must be a brmsfit")
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

test_that("predict_outliers survives the family validation brm() performs", {
  # The flag has to be an ordinary field on the family and reach
  # object$formula$family untouched, or prepare_predictions() never sees it.
  # Verified against brms 2.23.1.
  expect_false(rt_lognormal()$predict_outliers)
  expect_true(rt_lognormal(predict_outliers = TRUE)$predict_outliers)

  fam <- rt_lognormal(predict_outliers = TRUE)
  expect_true(brms:::validate_family(fam)$predict_outliers)

  f <- brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1, family = fam)
  expect_true(f$family$predict_outliers)

  d <- data.frame(RT = c(0.4, 0.5, 0.6, 0.9, 1.2))
  expect_true(brms:::validate_formula(f, data = d)$family$predict_outliers)
})

test_that("stanvars carry the likelihood with the hard-coded outlier component", {
  sv <- rt_lognormal_stanvars()
  code <- paste(vapply(sv, function(x) x$scode, character(1)), collapse = "\n")

  expect_true(grepl("rt_lognormal_lpdf", code, fixed = TRUE))
  expect_true(grepl("log_mix", code, fixed = TRUE))
  expect_true(grepl("student_t_lpdf(Y | 3, 0, 0.4)", code, fixed = TRUE))
})
