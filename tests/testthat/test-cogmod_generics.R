context("cogmod_priors / cogmod_inits / cogmod_stanvars")

# The three generics all read the family off the model rather than being told
# it, so most of what is worth testing is that they agree with the Stan program
# brms is actually going to build.

set.seed(1)
d_ig <- data.frame(
  RT = rcogmod_invgaussian(150, ndt = 0.2, poutlier = 0.02),
  Condition = rep(c("a", "b"), 75),
  id = rep(1:10, 15)
)
d_ex <- data.frame(
  RT = rcogmod_exgaussian(150),
  Condition = rep(c("a", "b"), 75),
  id = rep(1:10, 15)
)

# The names brms declares in the `parameters` block, for comparison against
# what cogmod_inits() supplies.
declared_params <- function(formula, data) {
  code <- suppressWarnings(
    brms::make_stancode(formula, data = data, family = formula$family)
  )
  vapply(cogmod:::.stan_param_decls(code), `[[`, character(1), "name")
}


# cogmod_inits: completeness ----------------------------------------------

test_that("cogmod_inits covers every declared parameter", {
  models <- list(
    simple = brms::bf(RT ~ 1, family = cogmod_invgaussian()),
    dpars = brms::bf(RT ~ Condition, boundary ~ Condition, ndt ~ Condition,
                     family = cogmod_invgaussian()),
    nonlinear_intercept = brms::bf(
      RT ~ 0 + Intercept + Condition, ndt ~ 0 + Intercept,
      family = cogmod_invgaussian()
    ),
    grouped = brms::bf(RT ~ Condition + (Condition | id), ndt ~ (1 | id),
                       family = cogmod_lognormal())
  )
  for (nm in names(models)) {
    f <- models[[nm]]
    vals <- cogmod_inits(f, d_ig)(1)
    expect_setequal(names(vals), declared_params(f, d_ig))
    expect_true(all(vapply(vals, function(v) all(is.finite(v)), logical(1))),
                info = nm)
  }
})


test_that("cogmod_inits matches the declared dimensions", {
  f <- brms::bf(RT ~ Condition + (Condition | id), ndt ~ 1,
                family = cogmod_lognormal())
  vals <- cogmod_inits(f, d_ig)(1)
  sdata <- suppressWarnings(
    brms::make_standata(f, data = d_ig, family = f$family)
  )
  expect_length(vals$b, sdata$Kc)
  expect_length(vals$sd_1, sdata$M_1)
  expect_equal(dim(vals$z_1), c(sdata$M_1, sdata$N_1))
  # A Cholesky factor has to start as one; the identity is the safe choice.
  expect_equal(vals$L_1, diag(sdata$M_1))
})


# cogmod_inits: values ----------------------------------------------------

test_that("cogmod_inits puts the targets on the right scale", {
  f <- brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1,
                family = cogmod_lognormal(minrt = 0.3))
  vals <- cogmod_inits(f, d_ig, jitter = 0)(1)
  # ndt starts at a third of minrt, on the log link
  expect_equal(vals$Intercept_ndt, log(0.1))
  expect_equal(vals$Intercept_poutlier, qlogis(0.02))
  # mu has an identity link here, sigma a softplus one
  expect_equal(vals$Intercept, -0.7)
  expect_equal(vals$Intercept_sigma, log(expm1(0.5)))
})


test_that("cogmod_inits leaves an omitted dpar on the natural scale", {
  f <- brms::bf(RT ~ 1, family = cogmod_lognormal(minrt = 0.3))
  vals <- cogmod_inits(f, d_ig, jitter = 0)(1)
  # Omitted from the formula, so brms declares it as a plain auxiliary
  # parameter with no link applied.
  expect_equal(vals$ndt, 0.1)
  expect_equal(vals$poutlier, 0.02)
})


test_that("cogmod_inits follows minrt and the links on the family", {
  f <- brms::bf(RT ~ 1, ndt ~ 1, family = cogmod_invgaussian(minrt = 0.24))
  expect_equal(cogmod_inits(f, d_ig, jitter = 0)(1)$Intercept_ndt, log(0.08))

  # mu's target is 3; softplus by default, log if asked for.
  g <- brms::bf(RT ~ 1, family = cogmod_invgaussian())
  h <- brms::bf(RT ~ 1, family = cogmod_invgaussian(link_mu = "log"))
  expect_equal(cogmod_inits(g, d_ig, jitter = 0)(1)$Intercept, log(expm1(3)))
  expect_equal(cogmod_inits(h, d_ig, jitter = 0)(1)$Intercept, log(3))
})


test_that("cogmod_inits gives the target to `Intercept` under 0 + Intercept", {
  f <- brms::bf(RT ~ Condition, ndt ~ 0 + Intercept + Condition,
                family = cogmod_lognormal(minrt = 0.3))
  vals <- cogmod_inits(f, d_ig, jitter = 0)(1)
  # There is no Intercept_ndt to set: the coefficient named Intercept inside
  # b_ndt is the intercept, and the slope beside it is a slope.
  expect_null(vals$Intercept_ndt)
  expect_equal(vals$b_ndt, c(log(0.1), 0))
})


test_that("cogmod_inits does not mistake a spline coefficient for a dpar", {
  # A smooth on mu declares `vector[Ks] bs`, brms' name for the unpenalized
  # spline coefficients. No cogmod dpar is called `bs` any more (that is why
  # the threshold is `boundary`), but only a scalar `real` is ever a dpar, so
  # the guard is tested rather than assumed.
  set.seed(2)
  dd <- cbind(d_ig, x = rnorm(nrow(d_ig)))
  f <- brms::bf(RT ~ s(x), boundary ~ 1, ndt ~ 1, family = cogmod_invgaussian())
  vals <- cogmod_inits(f, dd, jitter = 0)(1)
  expect_equal(vals$bs, 0)
  expect_equal(vals$Intercept_boundary, log(expm1(0.5)))
})


test_that("cogmod_inits jitters without leaving the bounds", {
  f <- brms::bf(RT ~ 1, family = cogmod_lognormal())
  inits <- cogmod_inits(f, d_ig, jitter = 1)
  draws <- replicate(300, unlist(inits(1)))
  expect_true(all(draws["ndt", ] > 0))
  expect_true(all(draws["poutlier", ] > 0 & draws["poutlier", ] < 1))
  expect_true(all(draws["sigma", ] > 0))
  # jitter = 0 is the deterministic case
  fixed <- cogmod_inits(f, d_ig, jitter = 0)
  expect_identical(fixed(1), fixed(2))
})


# cogmod_inits: families --------------------------------------------------

test_that("cogmod_inits supports cogmod_exgaussian", {
  f <- brms::bf(RT ~ Condition, sigma ~ Condition, tau ~ Condition + (1 | id),
                family = cogmod_exgaussian())
  vals <- cogmod_inits(f, d_ex, jitter = 0)(1)
  expect_setequal(names(vals), declared_params(f, d_ex))
  # All three are softplus by default and all three are targeted.
  expect_equal(vals$Intercept, log(expm1(0.4)))
  expect_equal(vals$Intercept_sigma, log(expm1(0.1)))
  expect_equal(vals$Intercept_tau, log(expm1(0.2)))

  # Omitted from the formula: natural scale, no link.
  g <- brms::bf(RT ~ 1, family = cogmod_exgaussian())
  omitted <- cogmod_inits(g, d_ex, jitter = 0)(1)
  expect_equal(omitted$sigma, 0.1)
  expect_equal(omitted$tau, 0.2)
})


test_that("cogmod_inits refuses a family it has nothing to offer", {
  expect_error(
    cogmod_inits(brms::bf(RT ~ 1, family = brms::lognormal()), d_ig),
    "nothing to offer"
  )
  expect_error(cogmod_inits(RT ~ 1, d_ig), "none found on the formula")
})


# cogmod_stanvars ---------------------------------------------------------

test_that("cogmod_stanvars dispatches on the family", {
  f <- brms::bf(RT ~ 1, family = cogmod_lognormal())
  expect_identical(cogmod_stanvars(f), cogmod_lognormal_stanvars())
  expect_identical(cogmod_stanvars(cogmod_exgaussian()), cogmod_exgaussian_stanvars())
  expect_identical(cogmod_stanvars(cogmod_choco()), cogmod_choco_stanvars())
  expect_identical(cogmod_stanvars(cogmod_lnr()), cogmod_lnr_stanvars())
})


test_that("cogmod_stanvars takes minrt off the family", {
  f <- brms::bf(RT ~ 1, family = cogmod_lognormal(minrt = 0.25))
  expect_identical(cogmod_stanvars(f), cogmod_lognormal_stanvars(minrt = 0.25))
  # ...which is the point: the family's minrt is baked into the Stan code, so
  # the default would silently disagree with it.
  expect_false(identical(cogmod_stanvars(f), cogmod_lognormal_stanvars()))
  # An explicit argument still wins.
  expect_identical(cogmod_stanvars(f, minrt = 0.4),
                   cogmod_lognormal_stanvars(minrt = 0.4))
})


test_that("cogmod_stanvars refuses a family it has no code for", {
  expect_error(
    cogmod_stanvars(brms::bf(RT ~ 1, family = brms::lognormal())),
    "no Stan code"
  )
  expect_error(cogmod_stanvars(RT ~ 1), "none found on the formula")
})


# cogmod_priors -----------------------------------------------------------

test_that("cogmod_priors returns rows that match real parameters", {
  f <- brms::bf(RT ~ Condition, ndt ~ 1, poutlier ~ 1,
                family = cogmod_lognormal())
  p <- cogmod_priors(f, d_ig)
  expect_s3_class(p, "brmsprior")
  expect_true(any(p$dpar == "ndt" & p$class == "Intercept" & nzchar(p$prior)))
})


# Stan declaration parsing ------------------------------------------------

test_that(".parse_stan_decl handles the shapes brms emits", {
  parse1 <- function(s) cogmod:::.parse_stan_decl(s)
  expect_equal(parse1("vector[Kc] b"),
               list(name = "b", type = "vector", dims = "Kc", bounds = ""))
  expect_equal(parse1("real<lower=0,upper=1> poutlier"),
               list(name = "poutlier", type = "real", dims = character(0),
                    bounds = "lower=0,upper=1"))
  # Dimensions nest, which is why this is not a single regex.
  expect_equal(parse1("vector[knots_1[1]] zs_1_1")$dims, "knots_1[1]")
  expect_equal(parse1("matrix[M_1, N_1] z_1")$dims, c("M_1", "N_1"))
  expect_equal(parse1("array[M_2] vector[N_2] z_2")$dims, c("M_2", "N_2"))
  expect_null(parse1("int not_a_parameter_type"))
})


# Complex formula shapes --------------------------------------------------

# The three generics all claim to survive whatever brms builds. This sweeps the
# formula features that change the *shape* of the Stan program - extra parameter
# blocks, extra dimensions, extra design matrices - and asserts three things for
# each: inits name exactly the declared parameters, priors validate without a
# complaint, and brms builds Stan code from the result silently.
test_that("the generics survive complex model specifications", {
  # `me()` is the one term brms cannot resolve from its own namespace: while
  # collecting grouping variables it eval()s the term on the search path. That
  # is a brms quirk rather than anything to do with cogmod, so the search path
  # is arranged for the duration rather than the spec being dropped.
  if (!"package:brms" %in% search()) {
    suppressMessages(library(brms))
    on.exit(detach("package:brms"), add = TRUE)
  }

  set.seed(7)
  n <- 200
  dc <- data.frame(
    RT = rcogmod_lognormal(n, mu = -0.7, sigma = 0.4, ndt = 0.2, poutlier = 0.02),
    Condition = factor(rep(c("a", "b"), length.out = n)),
    Block = factor(rep(c("p", "q", "r"), length.out = n)),
    x = rnorm(n), z = rnorm(n),
    ord = factor(rep(1:5, length.out = n), ordered = TRUE),
    id = factor(rep(1:15, length.out = n)),
    item = factor(rep(1:6, length.out = n))
  )
  dc$xse <- runif(n, 0.05, 0.2)
  fam <- cogmod_lognormal()

  specs <- list(
    smooth = brms::bf(RT ~ s(x), ndt ~ 1, family = fam),
    smooth_on_dpar = brms::bf(RT ~ Condition, ndt ~ s(x), family = fam),
    tensor = brms::bf(RT ~ t2(x, z), ndt ~ 1, family = fam),
    smooth_by = brms::bf(RT ~ s(x, by = Condition), ndt ~ 1, family = fam),
    re_correlated = brms::bf(RT ~ Condition + (1 + Condition | id), ndt ~ 1,
                             family = fam),
    re_uncorrelated = brms::bf(RT ~ Condition + (1 + Condition || id), ndt ~ 1,
                               family = fam),
    re_crossed = brms::bf(RT ~ (1 | id) + (1 | item), ndt ~ 1, family = fam),
    re_nested = brms::bf(RT ~ (1 | id / item), ndt ~ 1, family = fam),
    re_shared = brms::bf(RT ~ (1 | g | id), ndt ~ (1 | g | id), family = fam),
    re_several_dpars = brms::bf(RT ~ (1 | id), sigma ~ (1 | id),
                                ndt ~ (1 | id), family = fam),
    monotonic = brms::bf(RT ~ mo(ord), ndt ~ 1, family = fam),
    monotonic_on_dpar = brms::bf(RT ~ Condition, ndt ~ mo(ord), family = fam),
    interaction_poly = brms::bf(RT ~ Condition * Block + poly(x, 2), ndt ~ 1,
                                family = fam),
    no_intercept = brms::bf(RT ~ 0 + Intercept + Condition,
                            ndt ~ 0 + Intercept, family = fam),
    measurement_error = brms::bf(RT ~ me(x, xse), ndt ~ 1, family = fam),
    kitchen_sink = brms::bf(
      RT ~ s(x) + mo(ord) + (1 + Condition | id),
      sigma ~ Condition, ndt ~ 0 + Intercept + (1 | id), poutlier ~ 1,
      family = fam
    ),
    invgaussian = brms::bf(
      RT ~ s(x) + Condition + (1 + Condition | id) + (1 | item),
      boundary ~ mo(ord), ndt ~ 0 + Intercept + (1 | id),
      family = cogmod_invgaussian()
    ),
    loggamma = brms::bf(RT ~ Condition + (1 | id), shape ~ Condition,
                        ndt ~ s(x), family = cogmod_loggamma()),
    exgaussian = brms::bf(RT ~ s(x) + (1 + Condition | id), tau ~ mo(ord),
                          family = cogmod_exgaussian())
  )

  for (nm in names(specs)) {
    f <- specs[[nm]]
    vals <- cogmod_inits(f, dc)(1)
    expect_setequal(names(vals), declared_params(f, dc))
    expect_true(all(vapply(vals, function(v) all(is.finite(v)), logical(1))),
                label = paste("finite inits for", nm))

    # cogmod_exgaussian is passed through with a message by design; everything else
    # has to be silent, including the validate_prior() call inside.
    p <- if (identical(nm, "exgaussian")) {
      suppressMessages(cogmod_priors(f, dc))
    } else {
      expect_silent(cogmod_priors(f, dc))
    }
    # ...and brms must not warn that one of those rows goes unused.
    expect_silent(
      brms::make_stancode(f, data = dc, family = f$family, prior = p,
                          stanvars = cogmod_stanvars(f))
    )
  }
})


# Family invariants -------------------------------------------------------

# `cogmod_invgaussian(link_bs = ...)` survived the rename of the `bs` dpar to
# `boundary` for a while, because a word-boundary search does not see `bs`
# inside `link_bs`. The argument and the dpar it sets have to agree, or the
# argument silently sets nothing.
test_that("every link_<dpar> argument names a real dpar of its family", {
  for (nm in cogmod:::.cogmod_families()) {
    fam <- get(nm)()
    args <- sub("^link_", "",
                grep("^link_", names(formals(get(nm))), value = TRUE))
    expect_setequal(args, fam$dpars)
  }
})


test_that("the shifted registry agrees with the families it builds", {
  # .mixture_spec() rather than .shifted_spec(): .OUTLIER_FAMILIES now spans
  # both the RT-only registry and the choice registry (cogmod_lnr, so far),
  # and the shape checked below - dpars/links/lb/ub/init lining up - holds for
  # either kind.
  for (nm in cogmod:::.OUTLIER_FAMILIES) {
    spec <- cogmod:::.mixture_spec(nm)
    fam <- get(nm)()
    # the registry's decision dpars, then the two the parameterization adds
    expect_equal(fam$dpars, c(spec$dpars, "ndt", "poutlier"), label = nm)
    expect_length(spec$links, length(spec$dpars))
    expect_length(spec$lb, length(spec$dpars))
    expect_length(spec$ub, length(spec$dpars))
    # every decision dpar has a starting value
    expect_setequal(names(spec$init), spec$dpars)
    # the Stan parameter check names only real parameters
    for (tok in unlist(strsplit(spec$stan_check, "[^A-Za-z0-9_]+"))) {
      if (grepl("^[a-z][a-z0-9_]*$", tok)) {
        expect_true(tok %in% spec$dpars, label = paste(nm, tok))
      }
    }
  }
})


test_that("the R parameter checks reject exactly what Stan rejects", {
  # Both come from the registry's lb/ub now, so this pins them together.
  expect_error(rcogmod_lognormal(5, sigma = -1), "sigma")
  expect_error(rcogmod_lognormal(5, sigma = 0), "sigma")
  expect_error(rcogmod_gamma(5, mu = 0), "mu")
  expect_error(rcogmod_invgaussian(5, boundary = -0.1), "boundary")
  expect_error(rcogmod_lba1(5, sigmabias = 0), "sigmabias")
  # shape is unconstrained for the log-gamma, so it must not be rejected
  expect_silent(rcogmod_loggamma(5, shape = -3))
  # ndt = 0 and poutlier = 0 are legitimate, not boundary violations
  expect_silent(rcogmod_lognormal(5, ndt = 0, poutlier = 0))
  # a density gets a warning and zero rather than an error, so that one bad
  # posterior draw cannot abort a whole call
  expect_warning(d <- dcogmod_lognormal(0.5, sigma = -1), "must be greater")
  expect_equal(d, 0)
})


# The shared Stan model ---------------------------------------------------

test_that("the suite's shared Stan model stands in for *_lpdf_expose()", {
  # helper-stan.R compiles every family's lpdf into one model, because nine
  # separate compilations cost far more than the tests they serve. That is only
  # legitimate if the shared model returns what the per-family exposer returns,
  # so the substitution is checked here - once, against the one family whose
  # own compilation is being paid for anyway.
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  shared <- stan_fun("cogmod_lognormal")
  own <- cogmod_lognormal_lpdf_expose()
  args <- list(0.9, -0.7, 0.5, 0.3, 0.02)
  expect_equal(do.call(shared, args), do.call(own, args), tolerance = 1e-12)

  # and every family the helper claims to carry is actually in there, so a new
  # family added to the package without a line in .STAN_LPDF_FAMILIES fails
  # here rather than in whichever test happens to ask for it first
  for (nm in .STAN_LPDF_FAMILIES) {
    suffix <- if (identical(nm, "cogmod_betadiscrete")) "_lpmf" else "_lpdf"
    expect_true(is.function(stan_fun(nm, suffix)), label = nm)
  }
})
