# Shared machinery for the shifted RT families
# =============================================
#
# Every family in this group has the same structure: a decision-time
# distribution shifted by `ndt`, mixed with a half Student-t outlier component
# of weight `poutlier` and scale `minrt`. Only the decision distribution
# differs. Keeping the mixture in one place is what stops the nine families
# drifting apart - the Stan code, the R density, the RNG, the likelihood, the
# predictions, the parameter checks and the priors are all generated from the
# one registry below, so a family cannot end up with, say, an outlier component
# in Stan but not in `log_lik()`, or an R check that rejects what Stan accepts.
#
# See ?rcogmod_lognormal for the account of why `ndt` is expressed directly rather
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
.SHIFTED <- list(
  cogmod_lognormal = list(
    dpars = c("mu", "sigma"),
    links = c("identity", "softplus"),
    lb = c(NA, 0), ub = c(NA, NA),
    stan_check = "sigma <= 0",
    stan_dens = "lognormal_lpdf(t_adj | mu, sigma)",
    ldens = function(t, p) stats::dlnorm(t, p$mu, p$sigma, log = TRUE),
    rng = function(n, p) stats::rlnorm(n, p$mu, p$sigma),
    mean = function(p) exp(p$mu + p$sigma^2 / 2),
    init = list(mu = -0.7, sigma = 0.5),
    dpar_doc = c(
      "mu: mean of the decision time on the log scale (meanlog).",
      "sigma: SD of the decision time on the log scale (> 0)."
    ),
    label = "LogNormal"
  ),
  cogmod_loggamma = list(
    dpars = c("mu", "sigma", "shape"),
    links = c("identity", "softplus", "identity"),
    lb = c(NA, 0, NA), ub = c(NA, NA, NA),
    stan_check = "sigma <= 0",
    stan_dens = paste0(
      "-log(sigma) - log(t_adj) + cogmod_loggamma_lconst(shape)",
      " + cogmod_loggamma_lkernel(shape, (log(t_adj) - mu) / sigma)"
    ),
    prelude = ".LOGGAMMA_STAN_PRELUDE",
    ldens = function(t, p) .dloggamma(t, p$mu, p$sigma, p$shape, log = TRUE),
    rng = function(n, p) .rloggamma(n, p$mu, p$sigma, p$shape),
    mean = function(p) .mloggamma(p$mu, p$sigma, p$shape),
    init = list(mu = -0.7, sigma = 0.5, shape = 0),
    dpar_doc = c(
      "mu: location of the decision time on the log scale.",
      "sigma: scale of the decision time on the log scale (> 0).",
      paste("shape: shape of the log-gamma on the log-RT scale",
            "(unconstrained). shape = 0 is the"),
      "   LogNormal, shape = sigma the Gamma, shape = 1 the Weibull."
    ),
    note = c(
      paste("Note that when sigma * shape >= 1 the decision density is",
            "unbounded at `ndt`, so"),
      paste("the likelihood is unbounded too; the prior on shape is what",
            "keeps the sampler out"),
      "of that region. See ?rcogmod_loggamma."
    ),
    label = "Log-Gamma"
  ),
  cogmod_invgaussian = list(
    # mu is the drift rate and boundary the decision threshold, so the underlying
    # inverse Gaussian has mean boundary / mu and shape boundary^2.
    dpars = c("mu", "boundary"),
    links = c("softplus", "softplus"),
    lb = c(0, 0), ub = c(NA, NA),
    stan_check = "mu <= 0 || boundary <= 0",
    stan_dens = paste0(
      "0.5 * (2 * log(boundary) - (log(2) + log(pi())) - 3 * log(t_adj))",
      " - square(boundary) * square(t_adj - boundary / mu) / (2 * square(boundary / mu) * t_adj)"
    ),
    ldens = function(t, p) .dwald_raw(t, p$mu, p$boundary),
    rng = function(n, p) .rwald_raw(n, p$mu, p$boundary),
    mean = function(p) p$boundary / p$mu,
    init = list(mu = 3, boundary = 0.5),
    label = "Wald (inverse Gaussian)"
  ),
  cogmod_gamma = list(
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
    init = list(mu = 2, sigma = 0.2),
    label = "Gamma"
  ),
  cogmod_invgamma = list(
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
    init = list(mu = 4, sigma = 1.5),
    label = "inverse Gamma"
  ),
  cogmod_weibull = list(
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
    init = list(mu = 2, sigma = 0.5),
    label = "Weibull"
  ),
  cogmod_invweibull = list(
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
    init = list(mu = 3, sigma = 0.4),
    label = "inverse Weibull (Frechet)"
  ),
  cogmod_lba1 = list(
    # Single-accumulator LBA: mu is the drift rate, sigma its across-trial SD
    # (fixed to 1 by convention - the evidence scale is arbitrary), sigmabias
    # the start-point range A, and boundary the threshold offset, so b = sigmabias + boundary.
    dpars = c("mu", "sigma", "sigmabias", "boundary"),
    links = c("softplus", "softplus", "softplus", "softplus"),
    lb = c(NA, 0, 0, 0), ub = c(NA, NA, NA, NA),
    stan_check = "sigma <= 0 || sigmabias <= 0 || boundary <= 0",
    stan_dens = "cogmod_lba1_decision_lpdf(t_adj | mu, sigma, sigmabias, boundary)",
    prelude = ".LBA_STAN_PRELUDE",
    ldens = function(t, p) .dlba1_raw(t, p$mu, p$sigma, p$sigmabias, p$boundary),
    rng = function(n, p) .rlba1_raw(n, p$mu, p$sigma, p$sigmabias, p$boundary),
    # E[(b - U(0, A)) / v] with v a normal truncated at 0 diverges: the drift
    # density is positive at v = 0, so E[1/v] does not exist. posterior_epred()
    # therefore has nothing to return - see posterior_epred_cogmod_lba1().
    mean = NULL,
    init = list(mu = 3, sigma = 1, sigmabias = 0.5, boundary = 0.5),
    label = "single-accumulator LBA"
  ),
  cogmod_logweibull = list(
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
    init = list(mu = -0.8, sigma = 0.3),
    label = "Log-Weibull (Gumbel)"
  )
)


# The families built on the outlier mixture. with_outliers(),
# without_outliers(), p_outlier(), cogmod_priors() and cogmod_inits() all work
# on every one of them. The choice families in core_choice.R share the mixture
# but not the univariate density, so they have a registry of their own and are
# folded in here.
#' @keywords internal
.OUTLIER_FAMILIES <- c(names(.SHIFTED), .CHOICE_FAMILIES)


#' @keywords internal
.shifted_spec <- function(name) {
  spec <- .SHIFTED[[name]]
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
.shifted_family <- function(name, links = NULL, predict_outliers = FALSE,
                               minrt = 0.3) {
  spec <- .shifted_spec(name)
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
  # formula would be estimated rather than defaulted. See ?rcogmod_lognormal.
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
.shifted_lpdf <- function(name, minrt = 0.3, prelude = "") {
  spec <- .shifted_spec(name)
  # Families whose decision density needs helper functions declare them here.
  if (!is.null(spec$prelude)) prelude <- get(spec$prelude)
  # width = 1 so formatC does not pad the literal out with spaces
  scale <- formatC(.validate_minrt(minrt), format = "g", digits = 15, width = 1)

  # The outlier term log(2) + student_t_lpdf(Y | nu, 0, minrt) has *no*
  # parameter in it - nu, the location and the scale are all constants, and Y is
  # data. Called as written, Stan nevertheless recomputes its normalising
  # constant (two lgammas and a log) for every observation on every leapfrog
  # step. Folding that constant into a literal here and keeping only the
  # Y-dependent part measured ~1.4x faster per gradient evaluation on a
  # 4000-observation LogNormal fit, with the posterior unchanged.
  nu <- .POUTLIER_DF
  lc <- log(2) + lgamma((nu + 1) / 2) - lgamma(nu / 2) -
    0.5 * log(nu * pi) - log(.validate_minrt(minrt))
  lp_out <- sprintf(
    "%s - %s * log1p(square(Y / %s) / %s)",
    formatC(lc, format = "g", digits = 17, width = 1),
    formatC((nu + 1) / 2, format = "g", digits = 15, width = 1),
    scale,
    formatC(nu, format = "g", digits = 15, width = 1)
  )
  args <- paste(
    sprintf("real %s", c("Y", spec$dpars, "ndt", "poutlier")),
    collapse = ", "
  )
  # Families that have something of their own to say about their parameters, or
  # about a boundary in their parameter space, carry it in the registry so the
  # generated code keeps it.
  dpar_doc <- if (is.null(spec$dpar_doc)) "" else {
    paste0(paste0("// ", spec$dpar_doc, collapse = "\n"), "\n")
  }
  note <- if (is.null(spec$note)) "" else {
    paste0("//\n", paste0("// ", spec$note, collapse = "\n"), "\n")
  }
  sprintf(
    "%s
// Log-likelihood for one observation from the shifted %s model.
// Y: observed reaction time.
%s// ndt: non-decision time, same unit as Y (> 0).
// poutlier: proportion of responses from the outlier process, in [0, 1].
//
// The outlier component is a half Student-t with 3 degrees of freedom and scale
// %s (= minrt). It keeps the density strictly positive below `ndt`, where the
// shifted decision component has none. That is what removes the hard min-RT
// boundary and lets `ndt` be estimated directly rather than as a fraction of an
// observed minimum. Because the scale is tied to `minrt`, the likelihood is
// equivariant to the unit Y is measured in.
//
// It is written out rather than called as student_t_lpdf() because every one of
// its arguments is constant: written that way, Stan recomputes the normalising
// constant for every observation on every leapfrog step.
%sreal %s_lpdf(%s) {
    // Parameter checks
    if (%s || ndt < 0 || poutlier < 0 || poutlier > 1) {
      return negative_infinity();
    }
    if (Y <= 0) return negative_infinity();

    // The leading constant includes the log(2) that folds the symmetric
    // Student-t onto [0, Inf).
    real lp_out = %s;
    real t_adj  = Y - ndt;

    // Faster than the non-decision time: only the outlier component can have
    // produced this response.
    if (t_adj <= 0) return log(poutlier) + lp_out;

    return log_mix(poutlier, lp_out, %s);
}
",
    prelude, spec$label, dpar_doc, scale, note, name, args, spec$stan_check,
    lp_out, spec$stan_dens
  )
}


# Density, RNG and mean ---------------------------------------------------

# Recycle the decision dpars plus ndt/poutlier to a common length. `...` holds
# the decision dpars, by name.
#' @keywords internal
.prepare_shifted <- function(name, x = NULL, n = NULL, ndt, poutlier, ...) {
  spec <- .shifted_spec(name)
  dec <- list(...)
  missing <- setdiff(spec$dpars, names(dec))
  if (length(missing)) {
    stop("Missing parameter(s) for ", name, "(): ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  dec <- dec[spec$dpars]

  # The decision parameters are checked against the registry's own bounds, so
  # the R functions reject exactly what `stan_check` rejects and the two cannot
  # disagree. `lb` is strict throughout the registry - `sigma <= 0` rather than
  # `sigma < 0` - which is why the comparison is `<=`.
  for (j in seq_along(spec$dpars)) {
    d <- spec$dpars[j]
    if (!is.na(spec$lb[j]) && any(dec[[d]] <= spec$lb[j], na.rm = TRUE)) {
      stop("`", d, "` must be greater than ", spec$lb[j], ".", call. = FALSE)
    }
    if (!is.na(spec$ub[j]) && any(dec[[d]] >= spec$ub[j], na.rm = TRUE)) {
      stop("`", d, "` must be less than ", spec$ub[j], ".", call. = FALSE)
    }
  }
  # ndt and poutlier are not in the registry: every family shares them, and
  # both admit their boundary values (ndt = 0 is an unshifted model, and
  # poutlier = 0 switches the outlier component off).
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
  spec <- .shifted_spec(name)
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
.dshifted <- function(name, x, ndt, poutlier, minrt = 0.3, log = FALSE,
                         ...) {
  params <- tryCatch(
    .prepare_shifted(name, x = x, ndt = ndt, poutlier = poutlier, ...),
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
.rshifted <- function(name, n, ndt, poutlier, minrt = 0.3, ...) {
  spec <- .shifted_spec(name)
  params <- .prepare_shifted(name, n = n, ndt = ndt, poutlier = poutlier,
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
  spec <- .shifted_spec(name)
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

  ll <- do.call(.dshifted, c(
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

  out <- do.call(.rshifted, c(
    list(name = name, n = n_draws, ndt = ndt, poutlier = poutlier,
         minrt = .minrt(prep)),
    dec
  ))
  as.matrix(out)
}


#' @keywords internal
.posterior_epred_shifted <- function(name, prep, predict_outliers = NULL) {
  spec <- .shifted_spec(name)
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

# The decision component of cogmod_invgaussian(), written out here rather than
# delegating to dcogmod_invgaussian()/rcogmod_invgaussian(): those are the *mixture*
# functions and route back through .dshifted(), so calling them from the
# registry would recurse.
#' @keywords internal
.dwald_raw <- function(t, drift, boundary) {
  ig_mu <- boundary / drift
  lambda <- boundary^2
  0.5 * (log(lambda) - log(2 * pi) - 3 * log(t)) -
    lambda * (t - ig_mu)^2 / (2 * ig_mu^2 * t)
}

# Michael-Schucany-Haas two-root method.
#' @keywords internal
.rwald_raw <- function(n, drift, boundary) {
  ig_mu <- boundary / drift
  lambda <- boundary^2
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
.shifted_expose <- function(name, minrt = 0.3) {
  insight::check_if_installed("cmdstanr")
  stancode <- paste0(
    "functions {
",
    .shifted_lpdf(name, .as_minrt(minrt)),
    "}"
  )
  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions[[paste0(name, "_lpdf")]]
}


# Unshifted single-accumulator LBA ----------------------------------------

# Decision component of cogmod_lba1(), written out here for the same reason as the
# Wald: dcogmod_lba1() is the mixture and would recurse.
#' @keywords internal
.dlba1_raw <- function(t, drift, sigma, sigmabias, boundary) {
  st <- pmax(sigma * t, 1e-10)
  z1 <- (boundary - drift * t) / st
  n2 <- .lba1_dens_over_A(drift, sigma, st, z1, sigmabias / st)
  log_norm <- log1p(-exp(stats::pnorm(-drift / sigma, log.p = TRUE)))
  ifelse(n2 > 0, log(n2) - log_norm, -Inf)
}

#' @keywords internal
.rlba1_raw <- function(n, drift, sigma, sigmabias, boundary) {
  v <- .rnorm_truncated(n, mean = drift, sd = sigma, lower = 0)
  sp <- stats::runif(n, min = 0, max = sigmabias)
  (sigmabias + boundary - sp) / v
}

# Stan counterpart of .lba1_dens_over_A() and .dlba1_raw(). Same two branches, and
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
real cogmod_lba1_dens_over_A(real drift, real sigma, real st, real z1, real delta) {
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
real cogmod_lba1_decision_lpdf(real t, real drift, real sigma, real sigmabias, real boundary) {
  real st = fmax(sigma * t, 1e-10);
  real z1 = (boundary - drift * t) / st;
  real f = cogmod_lba1_dens_over_A(drift, sigma, st, z1, sigmabias / st);
  if (f <= 0) return negative_infinity();
  return log(f) - log1m_exp(std_normal_lcdf(-drift / sigma));
}
"


# Log-Gamma helpers -------------------------------------------------------

# The two functions the Log-Gamma decision density needs. They are pulled out
# of the lpdf as a prelude for the same reason as the LBA's: the density itself
# is one expression in the registry, and the maths that makes it stable near
# shape = 0 does not fit in one.
#' @keywords internal
.LOGGAMMA_STAN_PRELUDE <- "
// Normalizing constant of the standardized log-gamma with shape k = 1 / shape^2,
// with its leading k already cancelled against the equal and opposite k hidden
// in the kernel below. Doing that cancellation on paper is what keeps the
// density accurate near shape = 0, where both terms reach 1e8 by |shape| = 1e-4.
// Below |shape| = 1e-3 the Stirling series is used, where 1 / (12k) = shape^2 / 12; the
// two branches agree to about 1e-9 where they meet, so the density is smooth
// through the switch and finite at shape = 0 exactly.
real cogmod_loggamma_lconst(real shape) {
  if (abs(shape) < 1e-3) {
    return -0.5 * log(2 * pi()) - square(shape) / 12 + pow(shape, 6) / 360;
  }
  real k = inv_square(shape);
  return (k - 0.5) * log(k) - lgamma(k) - k;
}

// (shape * w - expm1(shape * w)) / shape^2, the log-gamma kernel with the same k removed.
// At shape = 0 this is exactly -w^2 / 2, the Normal kernel, which is why shape = 0
// gives back the LogNormal rather than approaching it.
real cogmod_loggamma_lkernel(real shape, real w) {
  real u = shape * w;
  if (abs(u) < 1e-6) {
    return -0.5 * square(w) - shape * pow(w, 3) / 6 - square(shape) * pow(w, 4) / 24;
  }
  return (u - expm1(u)) / square(shape);
}
"


# The outlier component ---------------------------------------------------

# Shared by every shifted family, so it lives here rather than with any one of
# them. `minrt` is the half-t's scale and also the unit the whole
# parameterization is equivariant to, which is why validating and reading it
# are shared too.

# The outlier component is a half Student-t on [0, Inf) with `.POUTLIER_DF`
# degrees of freedom and scale `minrt`.
#
# The half-t is flat at the origin (zero derivative), which matters: the fastest
# responses are the ones least likely to be decisions, so the component must not
# thin out there. A LogNormal or Gamma vanishes at 0 and an Exponential peaks
# there with maximal slope; all three get this wrong. df = 3 keeps the tails
# heavy while keeping the mean finite for posterior_epred().
.POUTLIER_DF <- 3

# `minrt` is used directly as the half-t scale, with no conversion factor. That
# fixes the component's shape relative to it: 61% of the mass falls below
# `minrt`, and the density there is still 66% of its peak, whatever value is
# supplied. Tying the scale to a quantity in the same unit as the data is what
# makes the likelihood equivariant to that unit.
#
# Raising it makes `poutlier` claim more of the body and the slow tail; lowering
# it much below 0.25 s stops the component covering contaminants that land near
# `ndt`, and `ndt` gets noisier. `ndt` itself is insensitive across the whole
# usable range.
#' @keywords internal
.validate_minrt <- function(minrt = 0.3) {
  if (!is.numeric(minrt) || length(minrt) != 1 || is.na(minrt) ||
    minrt <= 0) {
    stop("`minrt` must be a single positive number.", call. = FALSE)
  }
  minrt
}

# Read `minrt` off a family, a brmsprep or a brmsfit, falling back to the
# default for models fitted before it was adjustable.
#' @keywords internal
.minrt <- function(x) {
  m <- x$family$minrt
  if (is.null(m)) eval(formals(.validate_minrt)$minrt) else m
}

# Density of the half-t component.
#' @keywords internal
.dcontam <- function(x, minrt = 0.3, log = FALSE) {
  s <- .validate_minrt(minrt)
  ld <- log(2) + stats::dt(x / s, df = .POUTLIER_DF, log = TRUE) - log(s)
  ld[x <= 0] <- -Inf
  if (log) ld else exp(ld)
}

# Mean of the half-t; finite for df > 1.
#' @keywords internal
.mcontam <- function(minrt = 0.3) {
  nu <- .POUTLIER_DF
  2 * .validate_minrt(minrt) * sqrt(nu) * gamma((nu + 1) / 2) /
    (sqrt(pi) * (nu - 1) * gamma(nu / 2))
}


# Accept either a number or the family itself, so that
# `cogmod_lognormal_stanvars(fam)` cannot drift out of step with the family the
# model was fitted with.
#' @keywords internal
.as_minrt <- function(x) {
  if (inherits(x, "brmsfamily") || inherits(x, "customfamily") || is.list(x)) {
    m <- x$minrt
    if (is.null(m)) {
      stop(
        "`minrt` was given a family that carries no `minrt`. ",
        "Build it with cogmod_lognormal().",
        call. = FALSE
      )
    }
    return(m)
  }
  x
}


# Resolve whether predictions should include the outlier component. An explicit
# argument wins; otherwise the flag carried on the family is used (set by
# cogmod_lognormal() or with_outliers()), defaulting to excluding it.
#' @keywords internal
.predict_outliers <- function(predict_outliers, prep) {
  if (!is.null(predict_outliers)) {
    return(isTRUE(predict_outliers))
  }
  isTRUE(prep$family$predict_outliers)
}


# Switching the outlier component on and off ------------------------------

# User-facing, and shared by every family built on the mixture, so they belong
# with the mixture rather than with the LogNormal they were first written for.

#' Include or exclude the outlier component in predictions
#'
#' @description
#' Switches the `predict_outliers` flag on a model fitted with [cogmod_lognormal()]
#' or [cogmod_loggamma()], controlling whether `posterior_predict()` and `posterior_epred()` describe the
#' fitted mixture or the decision process alone.
#'
#' Predictions **exclude** the outlier component by default, because for almost
#' every downstream use it is a nuisance: it pulls expected values toward its own
#' mean (0.331 s at the default `minrt`, and proportional to it) and adds a
#' spike of implausibly fast draws to posterior predictive samples. It is also a deliberately fixed regularizer rather than a
#' claim about how guesses are distributed, so simulating from it means
#' simulating from something the model does not assert.
#'
#' `with_outliers()` restores the mixture. The main reason to want it is
#' `brms::pp_check()`: on untrimmed data the decision-only predictive has no fast
#' spike to match the one in the data, which reads as misfit. Use
#' `pp_check(with_outliers(m))` for a like-for-like check.
#'
#' The flag is stored **on the model** rather than passed as an argument, because
#' `brms` and the packages built on it (`insight`, `modelbased`,
#' `marginaleffects`, `emmeans`) do not forward extra arguments down to a custom
#' family's prediction methods - `posterior_epred()` reaches the family method
#' with `prep` and nothing else. Carrying it on the object is what makes it work
#' through all of them. The same flag can be set up front with
#' `cogmod_lognormal(predict_outliers = TRUE)`.
#'
#' `log_lik()` is unaffected and has no equivalent switch: the likelihood *is*
#' the mixture, and dropping a component from it would not be a different summary
#' of the same model but a different model. One consequence worth knowing is that
#' `posterior_predict()` and `log_lik()` do not describe the same distribution by
#' default. This also desyncs `loo_pit()`, `loo_predict()` and `bayes_R2()` from
#' `loo()`, not just hand-rolled checks - anything that compares a simulated
#' replicate against the likelihood should be run on `with_outliers()`.
#'
#' @param object A `brmsfit` fitted with [cogmod_lognormal()], [cogmod_loggamma()]
#'   or any other family built on the outlier mixture - see the *Supported
#'   families* section of [cogmod_priors()] for the full list, which includes
#'   the choice-and-RT families such as [cogmod_lnr()].
#'
#' @return The model, with the flag set. The fit itself is untouched - only how
#'   predictions are summarised changes.
#'
#' @examples
#' # f <- bf(RT ~ Condition, ndt ~ 1, poutlier ~ 1, family = cogmod_lognormal())
#' # m <- brms::brm(f, data = df, stanvars = cogmod_stanvars(f))
#' #
#' # # the decision process alone - the default, everywhere downstream
#' # brms::posterior_epred(m)
#' # modelbased::estimate_means(m, by = "Condition")
#' # marginaleffects::avg_predictions(m, by = "Condition")
#' #
#' # # the fitted mixture, e.g. for a like-for-like predictive check
#' # brms::pp_check(with_outliers(m))
#' # without_outliers(m)  # back to the default
#'
#' @export
with_outliers <- function(object) {
  .set_predict_outliers(object, TRUE)
}


#' @rdname with_outliers
#' @export
without_outliers <- function(object) {
  .set_predict_outliers(object, FALSE)
}


#' @keywords internal
.set_predict_outliers <- function(object, value) {
  if (!inherits(object, "brmsfit")) {
    stop("`object` must be a brmsfit.", call. = FALSE)
  }
  if (!isTRUE(object$family$name %in% .OUTLIER_FAMILIES)) {
    stop("`object` must have been fitted with one of the families carrying a ",
         "`poutlier` parameter: ",
         paste0(.OUTLIER_FAMILIES, "()", collapse = ", "), ".",
         call. = FALSE)
  }
  # brms rebuilds prep$family from object$formula$family, so the flag has to
  # live there to survive prepare_predictions(). object$family is set too, so
  # that inspecting the fit tells the same story.
  object$formula$family$predict_outliers <- value
  object$family$predict_outliers <- value
  object
}


#' Per-trial outlier probabilities
#'
#' @description
#' Posterior probability that each response was generated by the outlier
#' component rather than by the decision process, for a model fitted with any
#' family built on the outlier mixture - [cogmod_lognormal()], [cogmod_loggamma()],
#' [cogmod_lnr()] and the rest listed in the *Supported families* section of
#' [cogmod_priors()]. This is the mixture *responsibility*
#' `poutlier * g(rt) / (poutlier * g(rt) + (1 - poutlier) * f(rt - ndt))`,
#' averaged over posterior draws.
#'
#' Responses faster than `ndt` come out at 1, responses in the heart of the
#' distribution near 0, and responses in either tail somewhere in between - the
#' model discriminates by evidence rather than by a cutoff. A response in the
#' middle can still be an outlier; a low probability means the data cannot tell,
#' not that the trial is clean.
#'
#' Averaging the responsibility over draws gives `P(trial i came from the
#' outlier component | data)` directly, so the posterior mean *is* the quantity
#' of interest and there is no interval to report alongside it. Pass
#' `summary = FALSE` for the raw draws if you need the spread.
#'
#' @param object A `brmsfit` fitted with [cogmod_lognormal()], [cogmod_loggamma()]
#'   or any other family built on the outlier mixture - see the *Supported
#'   families* section of [cogmod_priors()].
#' @param summary Logical; if `TRUE` (default) returns a data frame with one row
#'   per observation. If `FALSE`, returns the full draws x observations matrix.
#'
#' @return A data frame with columns `rt` and `p_outlier`, in the order the
#'   observations appear in the model frame, or a draws x observations matrix if
#'   `summary = FALSE`.
#'
#' @examples
#' # f <- bf(RT ~ 1, ndt ~ 1, poutlier ~ 1, family = cogmod_lognormal())
#' # m <- brms::brm(f, data = df, stanvars = cogmod_stanvars(f))
#' # head(p_outlier(m))
#'
#' @export
p_outlier <- function(object, summary = TRUE) {
  if (!inherits(object, "brmsfit")) {
    stop("`object` must be a brmsfit.")
  }
  prep <- brms::prepare_predictions(object)
  if (!"Y" %in% names(prep$data)) {
    stop("Outcome variable 'Y' not found in the fitted model.")
  }
  y <- prep$data$Y
  n <- length(y)

  # Models fitted before the family carried a name fall back to the LogNormal,
  # which is the only family this ever supported.
  fam <- prep$family$name
  if (is.null(fam)) fam <- "cogmod_lognormal"
  spec <- .mixture_spec(fam)

  mu <- brms::get_dpar(prep, "mu")
  n_draws <- if (is.matrix(mu)) nrow(mu) else brms::ndraws(object)
  as_mat <- function(x) {
    if (is.matrix(x)) x else matrix(x, nrow = n_draws, ncol = n)
  }

  # Whatever decision dpars this family happens to have - the registry knows,
  # so nothing here needs to name them.
  dpars <- lapply(spec$dpars, function(d) as_mat(brms::get_dpar(prep, d)))
  names(dpars) <- spec$dpars
  ndt <- as_mat(brms::get_dpar(prep, "ndt"))
  poutlier <- as_mat(brms::get_dpar(prep, "poutlier"))
  Y <- matrix(y, nrow = n_draws, ncol = n, byrow = TRUE)

  # Done in logs throughout: both components underflow to exactly 0 in the tails
  # on the natural scale, which turns the ratio into 0/0 well before either
  # component is genuinely negligible.
  #
  # A choice family's contaminant has to produce a choice as well as a time, so
  # the outlier density carries the 1 / K that spreads the guess over the
  # response options, and the decision density needs the observed choice.
  if (is.null(spec$K)) {
    lp_out <- .dcontam(Y, minrt = .minrt(object), log = TRUE)
    lp_dec <- .ldec(fam, Y - ndt, dpars)
  } else {
    resp <- matrix(.dec_from_prep(prep), nrow = n_draws, ncol = n,
                   byrow = TRUE)
    lp_out <- .lout_choice(Y, spec$K, minrt = .minrt(object))
    lp_dec <- .ldec_choice(fam, Y - ndt, resp, dpars)
  }

  log_num <- log(poutlier) + lp_out
  log_den <- .log_mix(poutlier, lp_out, lp_dec)
  r <- exp(log_num - log_den)
  # Both components vanish (log_den = -Inf): nothing but the outlier could have
  # produced the response.
  r[!is.finite(r)] <- 1

  if (!summary) return(r)

  data.frame(rt = y, p_outlier = colMeans(r))
}
