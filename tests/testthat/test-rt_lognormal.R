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
  d_out <- 2 * dt(y / 0.3, df = 3) / 0.3
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
  # Without the outlier component these are exactly 0, which is what forces ndt
  # below the fastest observed RT.
  d <- drt_lognormal(0.1, ndt = 0.4, poutlier = 0.02)
  expect_gt(d, 0)
  expect_equal(d, 0.02 * 2 * dt(0.1 / 0.3, df = 3) / 0.3, tolerance = 1e-12)

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
    poutlier * (2 * 0.3 * sqrt(3) * gamma(2) / (sqrt(pi) * 2 * gamma(1.5)))

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
    poutlier * (2 * 0.3 * sqrt(3) * gamma(2) / (sqrt(pi) * 2 * gamma(1.5)))

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
    0.98 * (exp(mu + sigma^2 / 2) + ndt) + 0.02 * (2 * 0.3 * sqrt(3) * gamma(2) / (sqrt(pi) * 2 * gamma(1.5)))
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
    g <- 2 * dt(v / 0.3, df = 3) / 0.3
    f <- if (v > 0.25) dlnorm(v - 0.25, -0.9, 0.5) else 0
    0.02 * g / (0.02 * g + 0.98 * f)
  }, numeric(1))

  expect_equal(out$p_outlier, ref, tolerance = 1e-12)
  expect_equal(out$rt, y)
  expect_named(out, c("rt", "p_outlier"))

  # The shape the article plots: certainty below ndt, near zero in the bulk,
  # and rising again in the far slow tail because the half-t has heavier tails
  # than the LogNormal.
  expect_equal(out$p_outlier[1], 1)
  expect_lt(out$p_outlier[3], out$p_outlier[2])
  expect_gt(out$p_outlier[4], out$p_outlier[3])
})

test_that("p_outlier stays finite where the natural-scale ratio would be 0/0", {
  # Responsibility is computed in logs. On the natural scale both densities
  # underflow to exactly 0 in the far tails, turning the ratio into 0/0 before
  # either component is genuinely negligible.
  # sigma = 0.05 puts the decision process at 0.62-0.70 s, so 0.66 is the bulk
  # and everything else here is a far tail on one side or the other.
  y <- c(1e-8, 0.05, 0.66, 50, 1e6)
  nd <- 3
  prep <- list(
    data = list(Y = y),
    dpars = list(
      mu = matrix(-0.9, nd, length(y)),
      sigma = matrix(0.05, nd, length(y)), # tight enough to underflow dlnorm
      ndt = matrix(0.25, nd, length(y)),
      poutlier = matrix(1e-8, nd, length(y)) # tiny rate, as in early warmup
    )
  )

  out <- with_mocked_bindings(
    p_outlier(structure(list(), class = "brmsfit")),
    prepare_predictions = function(object, ...) prep,
    get_dpar = function(prep, dpar, ...) prep$dpars[[dpar]],
    .package = "brms"
  )

  expect_true(all(is.finite(out$p_outlier)))
  expect_true(all(out$p_outlier >= 0 & out$p_outlier <= 1))
  expect_equal(out$p_outlier[1], 1) # below ndt: only the outlier could produce it
  expect_lt(out$p_outlier[3], 1e-3) # in the bulk: the decision process owns it
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
  # minrt is a constant on the family, never a dpar: brms would estimate it
  expect_false("minrt" %in% fam$dpars)
  expect_equal(fam$minrt, 0.3)
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

test_that("stanvars carry the likelihood with the outlier component", {
  sv <- rt_lognormal_stanvars()
  code <- paste(vapply(sv, function(x) x$scode, character(1)), collapse = "\n")

  expect_true(grepl("rt_lognormal_lpdf", code, fixed = TRUE))
  expect_true(grepl("log_mix", code, fixed = TRUE))
  expect_true(grepl("student_t_lpdf(Y | 3, 0, 0.3)", code, fixed = TRUE))
})


# priors ------------------------------------------------------------------

test_that("brms leaves the ndt and poutlier intercepts flat", {
  # Both are improper, and both directions of each are flat in the likelihood
  # (poutlier -> 1, ndt -> 0), so a model fitted without priors has an improper
  # posterior. The article says so; this pins the upstream behaviour it rests
  # on. Verified against brms 2.23.1.
  d <- data.frame(RT = c(0.4, 0.5, 0.6, 0.9, 1.2))
  p <- brms::get_prior(
    brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1),
    family = rt_lognormal(), data = d
  )
  flat <- p$prior[p$class == "Intercept" & p$dpar %in% c("ndt", "poutlier")]
  expect_true(all(flat == ""))
})

test_that("cogmod_priors survives arbitrary formula shapes", {
  # The whole point of reading get_prior() rather than guessing: every row comes
  # from the model brms is actually going to build, so a prior matching no
  # parameter is impossible by construction.
  set.seed(1)
  d <- data.frame(
    RT = rrt_lognormal(200, ndt = 0.3, poutlier = 0.02),
    x = rnorm(200), g = factor(sample(letters[1:5], 200, TRUE)),
    cond = factor(sample(c("A", "B"), 200, TRUE))
  )
  fam <- rt_lognormal()
  forms <- list(
    plain = brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1, family = fam),
    no_intercept = brms::bf(RT ~ 0 + Intercept + x, sigma ~ 1,
                            ndt ~ 0 + Intercept + x, poutlier ~ 1, family = fam),
    interaction = brms::bf(RT ~ cond * x, sigma ~ 1, ndt ~ cond,
                           poutlier ~ 1 + (1 | g), family = fam),
    ranef_both = brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1 + (1 | g),
                          poutlier ~ 1 + (1 + x | g), family = fam),
    smooth = brms::bf(RT ~ s(x), sigma ~ 1, ndt ~ s(x), poutlier ~ 1,
                      family = fam)
  )

  for (nm in names(forms)) {
    p <- cogmod_priors(forms[[nm]], d)
    expect_s3_class(p, "brmsprior")
    # No ndt or poutlier row is left flat unless a blanket row covers it -
    # brms reports a blanket row plus one per coefficient, and filling the
    # blanket is what makes the per-coefficient ones proper.
    np <- p[p$dpar %in% c("ndt", "poutlier") & !nzchar(p$prior), ]
    covered <- vapply(seq_len(nrow(np)), function(i) {
      any(!nzchar(p$coef) & nzchar(p$prior) & p$class == np$class[i] &
            p$dpar == np$dpar[i] & p$group == np$group[i])
    }, logical(1))
    expect_true(all(covered), label = nm)
    expect_gt(sum(p$dpar %in% c("ndt", "poutlier") & nzchar(p$prior)), 0)
    # brms accepts it, silently: no unmatched prior, no unused blanket row
    expect_silent(
      brms::make_stancode(forms[[nm]], data = d, family = rt_lognormal(),
                          prior = p, stanvars = rt_lognormal_stanvars())
    )
  }

  # `0 + Intercept` has no Intercept class, so the location has to land on the
  # b coefficient named "Intercept" instead
  p <- cogmod_priors(forms$no_intercept, d)
  b <- p[p$dpar == "ndt" & p$class == "b", ]
  expect_equal(b$prior[b$coef == "Intercept"], "normal(-1.5, 0.15)")
  expect_equal(b$prior[!nzchar(b$coef)], "normal(0, 1)")
})

test_that("cogmod_priors passes other families straight through", {
  d <- data.frame(y = rnorm(30), x = rnorm(30))

  f <- brms::bf(y ~ x, family = stats::gaussian())
  # the message names the family it found, not "<none found>"
  expect_message(p <- cogmod_priors(f, d), "gaussian")
  expect_message(cogmod_priors(f, d), "rt_lognormal")
  # exactly the brms defaults, untouched
  expect_equal(p, brms::validate_prior(brms::empty_prior(), f, d,
                                       family = stats::gaussian()))

  # a formula carrying no family at all is a message, not an error
  expect_message(q <- cogmod_priors(brms::bf(y ~ x), d), "none found")
  expect_s3_class(q, "brmsprior")
})

test_that("cogmod_priors moves the ndt location with minrt", {
  d <- data.frame(RT = rrt_lognormal(50, ndt = 0.3, poutlier = 0.02))
  mk <- function(fam) {
    brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1, family = fam)
  }

  p <- cogmod_priors(mk(rt_lognormal()), d)
  expect_equal(p$prior[p$dpar == "ndt"], "normal(-1.5, 0.15)")

  # ndt is a location in time and shifts with the timescale; poutlier does not
  ms <- cogmod_priors(mk(rt_lognormal(minrt = 300)),
                      transform(d, RT = RT * 1000))
  expect_equal(
    ms$prior[ms$dpar == "ndt"],
    sprintf("normal(%s, 0.15)", formatC(-1.5 + log(1000), format = "g",
                                       digits = 4, width = 1))
  )
  expect_equal(ms$prior[ms$dpar == "poutlier"], "normal(-5, 1)")
})

test_that("returned priors compose and replace with c()", {
  d <- data.frame(RT = rrt_lognormal(50, ndt = 0.3, poutlier = 0.02))
  f <- brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1,
                family = rt_lognormal())

  p <- c(
    cogmod_priors(f, d),
    brms::prior(normal(-2, 0.1), class = "Intercept", dpar = "ndt"),
    replace = TRUE
  )
  expect_s3_class(p, "brmsprior")
  code <- brms::make_stancode(f, data = d, family = rt_lognormal(), prior = p,
                              stanvars = rt_lognormal_stanvars())
  expect_true(grepl("| -2, 0.1", code, fixed = TRUE))
})

test_that("the likelihood is flat in both directions a prior has to cover", {
  y <- rrt_lognormal(400, mu = -0.7, sigma = 0.5, ndt = 0.3, poutlier = 0.02)
  ll <- function(nd, p) sum(drt_lognormal(y, -0.7, 0.5, nd, p, log = TRUE))

  # poutlier -> 1: mu, sigma and ndt drop out of the density entirely
  expect_equal(ll(1e-9, 1 - 1e-12), ll(0.1, 1 - 1e-12), tolerance = 1e-8)
  # and the plateau is flat in poutlier itself once it saturates
  expect_equal(ll(0.3, plogis(40)), ll(0.3, plogis(60)), tolerance = 1e-8)

  # ndt -> 0 on a log link: flat in log(ndt), which has nothing to do with the
  # mixture and is why a prior on poutlier alone is not enough
  expect_equal(ll(1e-8, 0.02), ll(1e-12, 0.02), tolerance = 1e-6)
})


# minrt ----------------------------------------------------------------

test_that("minrt defaults to 0.5 and reproduces the fixed 0.4 scale", {
  expect_equal(rt_lognormal()$minrt, 0.3)
  expect_equal(cogmod:::.validate_minrt(), 0.3)
  expect_equal(cogmod:::.validate_minrt(0.3), 0.3)
  # the shipped article numbers
  expect_equal(cogmod:::.mcontam(), 0.3307, tolerance = 1e-3)
  expect_equal(
    drt_lognormal(0.1, ndt = 0.3, poutlier = 0.02),
    0.02 * 2 * dt(0.1 / 0.3, df = 3) / 0.3,
    tolerance = 1e-12
  )
})

test_that("minrt makes the density equivariant to the unit of measurement", {
  # The unmixed shifted LogNormal is exactly scale-equivariant; a fixed outlier
  # scale breaks that, and scaling it with minrt restores it.
  s <- drt_lognormal(1, mu = -0.7, sigma = 0.5, ndt = 0.3, poutlier = 0.02)
  ms <- 1000 * drt_lognormal(1000,
    mu = -0.7 + log(1000), sigma = 0.5, ndt = 300,
    poutlier = 0.02, minrt = 300
  )
  expect_equal(s, ms, tolerance = 1e-10)

  # and the fixed default does not have that property
  ms_fixed <- 1000 * drt_lognormal(1000,
    mu = -0.7 + log(1000), sigma = 0.5, ndt = 300, poutlier = 0.02
  )
  expect_false(isTRUE(all.equal(s, ms_fixed)))

  # simulation follows the same scale
  set.seed(1)
  expect_gt(
    median(rrt_lognormal(2e4, ndt = 300, poutlier = 1, minrt = 300)),
    100
  )
})

test_that("minrt reaches the Stan code and survives brms validation", {
  code <- rt_lognormal_stanvars(minrt = 300)[[1]]$scode
  expect_true(grepl("student_t_lpdf(Y | 3, 0, 300)", code, fixed = TRUE))

  # a family can be passed instead of a number, so the two cannot drift apart
  fam <- rt_lognormal(minrt = 300)
  expect_equal(fam$minrt, 300)
  expect_identical(rt_lognormal_stanvars(fam)[[1]]$scode, code)

  # it has to reach prep$family, like predict_outliers
  expect_equal(brms:::validate_family(fam)$minrt, 300)
  d <- data.frame(RT = c(400, 500, 600, 900, 1200))
  f <- brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1, family = fam)
  expect_equal(brms:::validate_formula(f, data = d)$family$minrt, 300)

  # minrt is NOT a dpar: brms would estimate it whenever a user omitted it
  expect_false("minrt" %in% fam$dpars)

  expect_error(rt_lognormal(minrt = -1), "positive")
  expect_error(rt_lognormal(minrt = c(1, 2)), "single")
})

test_that("brms methods pick minrt up off the family", {
  n_draws <- 8
  dpars <- list(
    mu = rep(-0.7 + log(1000), n_draws), sigma = rep(0.5, n_draws),
    ndt = rep(300, n_draws), poutlier = rep(0.02, n_draws)
  )
  prep <- structure(
    list(data = list(Y = 1000), dpars = dpars,
         family = list(name = "rt_lognormal", minrt = 300)),
    class = "brmsprep"
  )

  expect_equal(
    log_lik_rt_lognormal(1, prep),
    rep(log(drt_lognormal(1000, -0.7 + log(1000), 0.5, 300, 0.02,
                          minrt = 300)), n_draws),
    tolerance = 1e-10
  )

  # a family with no minrt falls back to the reference value, so models
  # fitted before it was adjustable still predict correctly
  bare <- prep
  bare$family <- list(name = "rt_lognormal")
  expect_equal(
    log_lik_rt_lognormal(1, bare),
    rep(log(drt_lognormal(1000, -0.7 + log(1000), 0.5, 300, 0.02)), n_draws),
    tolerance = 1e-10
  )

  # epred mixes in the component mean at the right scale
  mat <- function(v) matrix(v, 2, 2)
  p2 <- structure(
    list(dpars = list(mu = mat(-0.7), sigma = mat(0.5), ndt = mat(0.3),
                      poutlier = mat(0.02)),
         family = list(name = "rt_lognormal", minrt = 1,
                       predict_outliers = TRUE)),
    class = "brmsprep"
  )
  expect_equal(
    posterior_epred_rt_lognormal(p2),
    mat(0.98 * (exp(-0.7 + 0.5^2 / 2) + 0.3) + 0.02 * cogmod:::.mcontam(1))
  )
})


# omitted dpars -----------------------------------------------------------

test_that("omitting ndt/poutlier from bf() still yields proper priors", {
  # A dpar left out of bf() becomes a plain auxiliary parameter on the natural
  # scale, not a linear predictor on the link scale, so it needs a different
  # prior. brms's own defaults there are uniform(0, min_Y) for ndt - the very
  # min-RT bound this parameterization removes - and flat over [0, 1] for
  # poutlier, which puts half its mass above 0.5.
  set.seed(1)
  d <- data.frame(RT = rrt_lognormal(100, ndt = 0.3, poutlier = 0.02))
  f <- brms::bf(RT ~ 1, family = rt_lognormal())

  raw <- brms::get_prior(f, data = d, family = rt_lognormal())
  expect_match(raw$prior[raw$class == "ndt"], "min_Y")
  expect_equal(raw$prior[raw$class == "poutlier"], "")

  p <- cogmod_priors(f, d)
  expect_equal(p$prior[p$class == "ndt"], "lognormal(-1.2, 0.2)")
  expect_equal(p$prior[p$class == "poutlier"], "exponential(100)")

  code <- brms::make_stancode(f, data = d, family = rt_lognormal(), prior = p,
                              stanvars = rt_lognormal_stanvars())
  expect_true(grepl("lognormal_lpdf(ndt | -1.2, 0.2)", code, fixed = TRUE))
  expect_true(grepl("exponential_lpdf(poutlier | 100)", code, fixed = TRUE))
  expect_false(grepl("uniform_lpdf(ndt", code, fixed = TRUE))
})

test_that("the natural-scale priors describe the same belief as the link ones", {
  # lognormal(m, s) on ndt is exactly normal(m, s) on log(ndt)
  set.seed(1)
  x <- rlnorm(2e5, -1.2, 0.2)
  expect_equal(mean(log(x)), -1.2, tolerance = 0.01)
  expect_equal(sd(log(x)), 0.2, tolerance = 0.02)

  # exponential(100) keeps the centre of logit-normal(-5, 1) but moves the mode
  # to zero: omitting poutlier from the formula says you expect no outliers.
  # both medians sit at about 0.7%: 0.00693 vs 0.00669, i.e. within 0.03 of a
  # percentage point of each other (a relative check would be the wrong test
  # here - what matters is that neither says "expect a percent or two")
  expect_lt(abs(qexp(0.5, 100) - plogis(-5)), 5e-4)
  expect_lt(qexp(0.95, 100), 0.05)
  # essentially all of it lies inside the [0, 1] support, so truncation is moot
  expect_gt(pexp(1, 100), 1 - 1e-12)
})
