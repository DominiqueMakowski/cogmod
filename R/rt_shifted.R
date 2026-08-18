# Shared machinery for the shifted RT families
# =============================================
#
# Every family in this group has the same structure: a decision-time
# distribution shifted by `ndt`, mixed with a half Student-t outlier component
# of weight `poutlier` and scale `minrt`. Only the decision distribution
# differs. Keeping the mixture in one place is what stops the eight families
# drifting apart - the Stan code, the R density, the RNG, the likelihood, the
# predictions and the priors are all generated from the one registry below, so a
# family cannot end up with, say, an outlier component in Stan but not in
# `log_lik()`.
#
# See ?rrt_lognormal for the account of why `ndt` is expressed directly rather
# than as a fraction of an observed minimum, what the outlier component is for,
# and why `minrt` is a constant on the family rather than a dpar.


# The registry ------------------------------------------------------------

# Each entry describes the *decision* component only.
#
#  dpars/links/lb/ub : the distributional parameters, before ndt and poutlier
#  stan_check        : Stan expression that is TRUE for invalid parameters
#  stan_dens         : Stan expression for the decision log-density at `t_adj`
#  ldens             : R log-density of the decision component, at t > 0
#  rng               : R sampler for the decision component
#  mean              : E[decision time]; Inf where it does not exist
#  label             : human-readable name, used in the generated Stan comments
#' @keywords internal
.RT_SHIFTED <- list(
  rt_lognormal = list(
    dpars = c("mu", "sigma"),
    links = c("identity", "softplus"),
    lb = c(NA, 0), ub = c(NA, NA),
    stan_check = "sigma <= 0",
    stan_dens = "lognormal_lpdf(t_adj | mu, sigma)",
    ldens = function(t, p) stats::dlnorm(t, p$mu, p$sigma, log = TRUE),
    rng = function(n, p) stats::rlnorm(n, p$mu, p$sigma),
    mean = function(p) exp(p$mu + p$sigma^2 / 2),
    label = "LogNormal"
  ),
  rt_loggamma = list(
    dpars = c("mu", "sigma", "shape"),
    links = c("identity", "softplus", "identity"),
    lb = c(NA, 0, NA), ub = c(NA, NA, NA),
    stan_check = "sigma <= 0",
    stan_dens = paste0(
      "-log(sigma) - log(t_adj) + rt_loggamma_lconst(shape)",
      " + rt_loggamma_lkernel(shape, (log(t_adj) - mu) / sigma)"
    ),
    ldens = function(t, p) .dloggamma(t, p$mu, p$sigma, p$shape, log = TRUE),
    rng = function(n, p) .rloggamma(n, p$mu, p$sigma, p$shape),
    mean = function(p) .mloggamma(p$mu, p$sigma, p$shape),
    label = "Log-Gamma"
  ),
  rt_invgaussian = list(
    # mu is the drift rate and bs the decision threshold, so the underlying
    # inverse Gaussian has mean bs / mu and shape bs^2.
    dpars = c("mu", "bs"),
    links = c("softplus", "softplus"),
    lb = c(0, 0), ub = c(NA, NA),
    stan_check = "mu <= 0 || bs <= 0",
    stan_dens = paste0(
      "0.5 * (2 * log(bs) - (log(2) + log(pi())) - 3 * log(t_adj))",
      " - square(bs) * square(t_adj - bs / mu) / (2 * square(bs / mu) * t_adj)"
    ),
    ldens = function(t, p) .dwald_raw(t, p$mu, p$bs),
    rng = function(n, p) .rwald_raw(n, p$mu, p$bs),
    mean = function(p) p$bs / p$mu,
    label = "Wald (inverse Gaussian)"
  ),
  rt_gamma = list(
    dpars = c("mu", "sigma"), # shape, scale
    links = c("softplus", "softplus"),
    lb = c(0, 0), ub = c(NA, NA),
    stan_check = "mu <= 0 || sigma <= 0",
    stan_dens = "gamma_lpdf(t_adj | mu, inv(sigma))",
    ldens = function(t, p) {
      stats::dgamma(t, shape = p$mu, scale = p$sigma, log = TRUE)
    },
    rng = function(n, p) stats::rgamma(n, shape = p$mu, scale = p$sigma),
    mean = function(p) p$mu * p$sigma,
    label = "Gamma"
  ),
  rt_invgamma = list(
    dpars = c("mu", "sigma"), # shape, scale
    links = c("softplus", "softplus"),
    lb = c(0, 0), ub = c(NA, NA),
    stan_check = "mu <= 0 || sigma <= 0",
    stan_dens = "inv_gamma_lpdf(t_adj | mu, sigma)",
    ldens = function(t, p) {
      p$mu * log(p$sigma) - lgamma(p$mu) - (p$mu + 1) * log(t) - p$sigma / t
    },
    rng = function(n, p) 1 / stats::rgamma(n, shape = p$mu, rate = p$sigma),
    # The mean of an inverse Gamma exists only for shape > 1.
    mean = function(p) ifelse(p$mu > 1, p$sigma / (p$mu - 1), Inf),
    label = "inverse Gamma"
  ),
  rt_weibull = list(
    dpars = c("mu", "sigma"), # shape, scale
    links = c("softplus", "softplus"),
    lb = c(0, 0), ub = c(NA, NA),
    stan_check = "mu <= 0 || sigma <= 0",
    stan_dens = "weibull_lpdf(t_adj | mu, sigma)",
    ldens = function(t, p) {
      stats::dweibull(t, shape = p$mu, scale = p$sigma, log = TRUE)
    },
    rng = function(n, p) stats::rweibull(n, shape = p$mu, scale = p$sigma),
    mean = function(p) p$sigma * gamma(1 + 1 / p$mu),
    label = "Weibull"
  ),
  rt_invweibull = list(
    dpars = c("mu", "sigma"), # shape, scale
    links = c("softplus", "softplus"),
    lb = c(0, 0), ub = c(NA, NA),
    stan_check = "mu <= 0 || sigma <= 0",
    stan_dens = "frechet_lpdf(t_adj | mu, sigma)",
    ldens = function(t, p) {
      log(p$mu) - log(p$sigma) - (1 + p$mu) * (log(t) - log(p$sigma)) -
        (t / p$sigma)^(-p$mu)
    },
    rng = function(n, p) p$sigma * (-log(stats::runif(n)))^(-1 / p$mu),
    # E[Frechet] is finite only for shape > 1.
    mean = function(p) ifelse(p$mu > 1, p$sigma * gamma(1 - 1 / p$mu), Inf),
    label = "inverse Weibull (Frechet)"
  ),
  rt_lba = list(
    # Single-accumulator LBA: mu is the drift rate, sigma its across-trial SD
    # (fixed to 1 by convention - the evidence scale is arbitrary), sigmabias
    # the start-point range A, and bs the threshold offset, so b = sigmabias + bs.
    dpars = c("mu", "sigma", "sigmabias", "bs"),
    links = c("softplus", "softplus", "softplus", "softplus"),
    lb = c(NA, 0, 0, 0), ub = c(NA, NA, NA, NA),
    stan_check = "sigma <= 0 || sigmabias <= 0 || bs <= 0",
    stan_dens = "rt_lba_decision_lpdf(t_adj | mu, sigma, sigmabias, bs)",
    prelude = ".LBA_STAN_PRELUDE",
    ldens = function(t, p) .dlba_raw(t, p$mu, p$sigma, p$sigmabias, p$bs),
    rng = function(n, p) .rlba_raw(n, p$mu, p$sigma, p$sigmabias, p$bs),
    # E[(b - U(0, A)) / v] with v a normal truncated at 0 diverges: the drift
    # density is positive at v = 0, so E[1/v] does not exist. posterior_epred()
    # therefore has nothing to return - see posterior_epred_rt_lba().
    mean = NULL,
    label = "single-accumulator LBA"
  ),
  rt_logweibull = list(
    # log(decision time) ~ Gumbel(mu, sigma), i.e. the log-Weibull.
    dpars = c("mu", "sigma"),
    links = c("identity", "softplus"),
    lb = c(NA, 0), ub = c(NA, NA),
    stan_check = "sigma <= 0",
    stan_dens = "gumbel_lpdf(log(t_adj) | mu, sigma) - log(t_adj)",
    ldens = function(t, p) {
      lt <- log(t)
      z <- (lt - p$mu) / p$sigma
      -log(p$sigma) - z - exp(-z) - lt
    },
    rng = function(n, p) exp(p$mu - p$sigma * log(-log(stats::runif(n)))),
    # E[exp(Gumbel)] = exp(mu) * Gamma(1 - sigma), and only exists for
    # sigma < 1. Note this is NOT exp(mu + sigma * gamma_euler), which is the
    # geometric mean rather than the mean.
    mean = function(p) ifelse(p$sigma < 1, exp(p$mu) * gamma(1 - p$sigma), Inf),
    label = "Log-Weibull (Gumbel)"
  )
)


# The families built on this machinery. with_outliers(), without_outliers(),
# p_outlier() and cogmod_priors() all work on every one of them.
#' @keywords internal
.OUTLIER_FAMILIES <- names(.RT_SHIFTED)


#' @keywords internal
.rt_spec <- function(name) {
  spec <- .RT_SHIFTED[[name]]
  if (is.null(spec)) {
    stop(
      "'", name, "' is not one of the shifted RT families. Supported: ",
      paste0(.OUTLIER_FAMILIES, "()", collapse = ", "), ".",
      call. = FALSE
    )
  }
  spec
}


# Family constructor ------------------------------------------------------

# Builds the brms custom family for `name`, appending `ndt` and `poutlier` to
# whatever decision dpars the registry lists.
#' @keywords internal
.rt_shifted_family <- function(name, links = NULL, predict_outliers = FALSE,
                               minrt = 0.3) {
  spec <- .rt_spec(name)
  dec_links <- if (is.null(links)) spec$links else links
  if (length(dec_links) != length(spec$dpars)) {
    stop("Expected ", length(spec$dpars), " link functions for ", name, "().",
         call. = FALSE)
  }
  fam <- brms::custom_family(
    name = name,
    dpars = c(spec$dpars, "ndt", "poutlier"),
    links = c(dec_links, "log", "logit"),
    lb = c(spec$lb, 0, 0),
    ub = c(spec$ub, NA, 1),
    type = "real"
  )
  # Both ride on the family because that is the only thing brms carries down to
  # a custom family's prediction methods, and because a dpar left out of the
  # formula would be estimated rather than defaulted. See ?rrt_lognormal.
  fam$predict_outliers <- isTRUE(predict_outliers)
  invisible(.validate_minrt(minrt))
  fam$minrt <- minrt
  fam
}


# Stan code ---------------------------------------------------------------

# Generates the `<name>_lpdf` Stan function: the same mixture skeleton for every
# family, with only the decision density swapped in. `minrt` is baked in as a
# literal because Stan functions cannot see the data block, and because a dpar
# would be estimated whenever the user left it out of the formula.
#' @keywords internal
.rt_shifted_lpdf <- function(name, minrt = 0.3, prelude = "") {
  spec <- .rt_spec(name)
  # Families whose decision density needs helper functions declare them here.
  if (!is.null(spec$prelude)) prelude <- get(spec$prelude)
  # width = 1 so formatC does not pad the literal out with spaces
  scale <- formatC(.validate_minrt(minrt), format = "g", digits = 15, width = 1)
  args <- paste(
    sprintf("real %s", c("Y", spec$dpars, "ndt", "poutlier")),
    collapse = ", "
  )
  sprintf(
    "%s
// Log-likelihood for one observation from the shifted %s model.
// Y: observed reaction time.
// ndt: non-decision time, same unit as Y (> 0).
// poutlier: proportion of responses from the outlier process, in [0, 1].
//
// The outlier component is a half Student-t with 3 degrees of freedom and scale
// %s (= minrt). It keeps the density strictly positive below `ndt`, where the
// shifted decision component has none. That is what removes the hard min-RT
// boundary and lets `ndt` be estimated directly rather than as a fraction of an
// observed minimum. Because the scale is tied to `minrt`, the likelihood is
// equivariant to the unit Y is measured in.
real %s_lpdf(%s) {
    // Parameter checks
    if (%s || ndt < 0 || poutlier < 0 || poutlier > 1) {
      return negative_infinity();
    }
    if (Y <= 0) return negative_infinity();

    // log(2) folds the symmetric Student-t onto [0, Inf).
    real lp_out = log(2) + student_t_lpdf(Y | 3, 0, %s);
    real t_adj  = Y - ndt;

    // Faster than the non-decision time: only the outlier component can have
    // produced this response.
    if (t_adj <= 0) return log(poutlier) + lp_out;

    return log_mix(poutlier, lp_out, %s);
}
",
    prelude, spec$label, scale, name, args, spec$stan_check, scale,
    spec$stan_dens
  )
}


# Density, RNG and mean ---------------------------------------------------

# Recycle the decision dpars plus ndt/poutlier to a common length. `...` holds
# the decision dpars, by name.
#' @keywords internal
.prepare_rt_shifted <- function(name, x = NULL, n = NULL, ndt, poutlier, ...) {
  spec <- .rt_spec(name)
  dec <- list(...)
  missing <- setdiff(spec$dpars, names(dec))
  if (length(missing)) {
    stop("Missing parameter(s) for ", name, "(): ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  dec <- dec[spec$dpars]

  if (any(ndt < 0, na.rm = TRUE)) stop("`ndt` must be non-negative.")
  if (any(poutlier < 0 | poutlier > 1, na.rm = TRUE)) {
    stop("`poutlier` must be in [0, 1].")
  }

  lens <- c(vapply(dec, length, integer(1)), length(ndt), length(poutlier))
  if (!is.null(x)) {
    m <- max(length(x), lens)
    if (m == 0) stop("At least one input vector must have non-zero length.")
  } else if (!is.null(n)) {
    if (length(n) > 1) n <- length(n)
    if (length(n) != 1 || n < 0 || n != floor(n)) {
      stop("n must be a single non-negative integer.")
    }
    m <- max(n, lens)
    if (n == 0) m <- 0
  } else {
    stop("Either 'x' or 'n' must be provided.")
  }

  params <- lapply(dec, rep_len, m)
  params$ndt <- rep_len(ndt, m)
  params$poutlier <- rep_len(poutlier, m)
  if (!is.null(x)) params$x <- rep_len(x, m)
  params$ndraws <- m
  params
}


# Log-density of the decision component alone, at t (which may be <= 0, giving
# -Inf). Shape is preserved, so this works on draws x observations matrices.
#' @keywords internal
.ldec <- function(name, t, p) {
  spec <- .rt_spec(name)
  ok <- is.finite(t) & t > 0
  # pmax() only keeps the density functions from being handed a non-positive
  # number; those entries are overwritten with -Inf immediately afterwards.
  ld <- spec$ldens(pmax(t, 1e-300), p)
  ld[!ok] <- -Inf
  ld[is.na(ld)] <- -Inf
  dim(ld) <- dim(t)
  ld
}


# Mixture log-density, shared by every family's d*() function.
#' @keywords internal
.drt_shifted <- function(name, x, ndt, poutlier, minrt = 0.3, log = FALSE,
                         ...) {
  params <- tryCatch(
    .prepare_rt_shifted(name, x = x, ndt = ndt, poutlier = poutlier, ...),
    error = function(e) {
      warning(conditionMessage(e), ". Returning 0 density / -Inf log-density.")
      list(ndraws = length(x), error = TRUE)
    }
  )
  if (!is.null(params$error)) {
    return(rep(ifelse(log, -Inf, 0), params$ndraws))
  }

  lp_out <- .dcontam(params$x, minrt = minrt, log = TRUE)
  lp_dec <- .ldec(name, params$x - params$ndt, params)

  ld <- .log_mix(params$poutlier, lp_out, lp_dec)
  ld[is.na(ld)] <- -Inf
  if (log) ld else exp(ld)
}


# Mixture RNG, shared by every family's r*() function.
#' @keywords internal
.rrt_shifted <- function(name, n, ndt, poutlier, minrt = 0.3, ...) {
  spec <- .rt_spec(name)
  params <- .prepare_rt_shifted(name, n = n, ndt = ndt, poutlier = poutlier,
                                ...)
  m <- params$ndraws

  is_out <- stats::runif(m) < params$poutlier
  out <- numeric(m)

  if (any(is_out)) {
    out[is_out] <- abs(
      .validate_minrt(minrt) * stats::rt(sum(is_out), df = .POUTLIER_DF)
    )
  }
  if (any(!is_out)) {
    keep <- !is_out
    p <- lapply(params[spec$dpars], function(v) v[keep])
    out[keep] <- spec$rng(sum(keep), p) + params$ndt[keep]
  }
  out
}


# brms methods ------------------------------------------------------------

# Pull the decision dpars for observation i (or all of them) off a brmsprep.
#' @keywords internal
.dpars_from_prep <- function(name, prep, i = NULL) {
  spec <- .rt_spec(name)
  out <- lapply(spec$dpars, function(d) {
    if (is.null(i)) brms::get_dpar(prep, d) else brms::get_dpar(prep, d, i = i)
  })
  stats::setNames(out, spec$dpars)
}


#' @keywords internal
.log_lik_shifted <- function(name, i, prep) {
  if (!"Y" %in% names(prep$data)) {
    stop("Outcome variable 'Y' not found in prep$data.")
  }
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_)

  dec <- .dpars_from_prep(name, prep, i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  poutlier <- brms::get_dpar(prep, "poutlier", i = i)

  n_draws <- max(vapply(c(dec, list(ndt, poutlier)), length, integer(1)))
  if (n_draws == 0) return(numeric(0))

  ll <- do.call(.drt_shifted, c(
    list(name = name, x = rep(y, length.out = n_draws), ndt = ndt,
         poutlier = poutlier, minrt = .minrt(prep), log = TRUE),
    dec
  ))
  ll[is.na(ll)] <- -Inf
  ll
}


#' @keywords internal
.posterior_predict_shifted <- function(name, i, prep, predict_outliers = NULL) {
  dec <- .dpars_from_prep(name, prep, i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  poutlier <- if (.predict_outliers(predict_outliers, prep)) {
    brms::get_dpar(prep, "poutlier", i = i)
  } else {
    0
  }
  n_draws <- max(vapply(c(dec, list(ndt)), length, integer(1)))

  out <- do.call(.rrt_shifted, c(
    list(name = name, n = n_draws, ndt = ndt, poutlier = poutlier,
         minrt = .minrt(prep)),
    dec
  ))
  as.matrix(out)
}


#' @keywords internal
.posterior_epred_shifted <- function(name, prep, predict_outliers = NULL) {
  spec <- .rt_spec(name)
  dec <- .dpars_from_prep(name, prep)
  ndt <- brms::get_dpar(prep, "ndt")

  # get_dpar() returns a bare scalar for a dpar that is constant across draws
  # and observations, so the draws x observations shape is taken from whichever
  # one is actually a matrix.
  dims <- NULL
  for (d in c(dec, list(ndt))) {
    if (is.matrix(d)) {
      dims <- dim(d)
      break
    }
  }
  n <- if (is.null(dims)) max(vapply(dec, length, integer(1))) else prod(dims)
  dec <- lapply(dec, rep_len, n)

  if (is.null(spec$mean)) {
    stop(
      "The decision time of ", name, "() has no finite mean, so there is no ",
      "expectation to return. Use posterior_predict() and summarise the draws ",
      "instead.",
      call. = FALSE
    )
  }
  mean_dec <- spec$mean(dec)
  if (!is.null(dims) && length(mean_dec) == prod(dims)) dim(mean_dec) <- dims
  mean_dec <- mean_dec + ndt

  if (!.predict_outliers(predict_outliers, prep)) {
    return(mean_dec)
  }
  poutlier <- brms::get_dpar(prep, "poutlier")
  (1 - poutlier) * mean_dec + poutlier * .mcontam(.minrt(prep))
}


# Unshifted Wald ----------------------------------------------------------

# The decision component of rt_invgaussian(), written out here rather than
# delegating to drt_invgaussian()/rrt_invgaussian(): those are the *mixture*
# functions and route back through .drt_shifted(), so calling them from the
# registry would recurse.
#' @keywords internal
.dwald_raw <- function(t, drift, bs) {
  ig_mu <- bs / drift
  lambda <- bs^2
  0.5 * (log(lambda) - log(2 * pi) - 3 * log(t)) -
    lambda * (t - ig_mu)^2 / (2 * ig_mu^2 * t)
}

# Michael-Schucany-Haas two-root method.
#' @keywords internal
.rwald_raw <- function(n, drift, bs) {
  ig_mu <- bs / drift
  lambda <- bs^2
  y <- stats::rnorm(n)^2
  z <- y * (ig_mu / lambda)
  x1_over_mu <- 1 + z / 2 * (1 - sqrt(1 + 4 / z))
  u <- stats::runif(n)
  ig_mu * ifelse(u < 1 / (1 + x1_over_mu), x1_over_mu, 1 / x1_over_mu)
}


# CDF of the half Student-t outlier component.
#' @keywords internal
.pcontam <- function(x, minrt = 0.3) {
  s <- .validate_minrt(minrt)
  out <- 2 * stats::pt(x / s, df = .POUTLIER_DF) - 1
  out[x <= 0] <- 0
  out
}


# Shared expose helper: compiles the generated lpdf and hands it back as an R
# function, for checking the Stan code against the R density.
#' @keywords internal
.rt_shifted_expose <- function(name, minrt = 0.3) {
  insight::check_if_installed("cmdstanr")
  stancode <- paste0(
    "functions {
",
    .rt_shifted_lpdf(name, .as_minrt(minrt)),
    "}"
  )
  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions[[paste0(name, "_lpdf")]]
}


# Unshifted single-accumulator LBA ----------------------------------------

# Decision component of rt_lba(), written out here for the same reason as the
# Wald: drt_lba() is the mixture and would recurse.
#' @keywords internal
.dlba_raw <- function(t, drift, sigma, sigmabias, bs) {
  st <- pmax(sigma * t, 1e-10)
  z1 <- (bs - drift * t) / st
  n2 <- .lba_dens_over_A(drift, sigma, st, z1, sigmabias / st)
  log_norm <- log1p(-exp(stats::pnorm(-drift / sigma, log.p = TRUE)))
  ifelse(n2 > 0, log(n2) - log_norm, -Inf)
}

#' @keywords internal
.rlba_raw <- function(n, drift, sigma, sigmabias, bs) {
  v <- .rnorm_truncated(n, mean = drift, sd = sigma, lower = 0)
  sp <- stats::runif(n, min = 0, max = sigmabias)
  (sigmabias + bs - sp) / v
}

# Stan counterpart of .lba_dens_over_A() and .dlba_raw(). Same two branches, and
# the same reason for them: both differences cancel as the start-point range
# goes to zero, and the result is then divided by that range.
#' @keywords internal
.LBA_STAN_PRELUDE <- "
// LBA defective density divided by the start-point range A, with
// z2 = z1 + delta and delta = A / (sigma * t). Both differences below vanish
// linearly in delta, so evaluating them directly and dividing by A loses every
// significant digit for a small start-point range. The Taylor expansion is used
// there instead; its truncation error is ~1e-12 at the switch, where the direct
// form still has ~12 good digits.
real rt_lba_dens_over_A(real drift, real sigma, real st, real z1, real delta) {
  real phi1 = exp(-0.5 * square(z1)) * 0.3989422804014327;
  if (delta < 1e-4) {
    real series = (drift + sigma * z1)
      - (delta / 2) * (drift * z1 + sigma * (square(z1) - 1))
      + (square(delta) / 6) * (drift * (square(z1) - 1)
                               + sigma * (z1 * square(z1) - 3 * z1));
    return phi1 * series / st;
  }
  real z2 = z1 + delta;
  // upper tail when both are positive, so Phi(z2) - Phi(z1) does not cancel
  real dPhi = z1 > 0 ? (Phi(-z1) - Phi(-z2)) : (Phi(z2) - Phi(z1));
  real dphi = phi1 * -expm1(-delta * (z1 + z2) / 2);
  return (drift * dPhi + sigma * dphi) / (delta * st);
}

// Log density of the single-accumulator LBA decision time (no shift).
real rt_lba_decision_lpdf(real t, real drift, real sigma, real sigmabias, real bs) {
  real st = fmax(sigma * t, 1e-10);
  real z1 = (bs - drift * t) / st;
  real f = rt_lba_dens_over_A(drift, sigma, st, z1, sigmabias / st);
  if (f <= 0) return negative_infinity();
  return log(f) - log1m_exp(std_normal_lcdf(-drift / sigma));
}
"
