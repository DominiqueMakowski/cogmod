#' Ex-Gaussian Model (Classical Parameterization) for brms
#'
#' Provides the necessary functions to use the "classical" parameterization
#' of the Ex-Gaussian distribution as a custom family in `brms`. Unlike
#' `brms`'s built-in `exgaussian()` family - in which `mu` indexes the mean of
#' the *entire* distribution (Gaussian + exponential components combined) -
#' this parameterization follows the convention familiar to experimental
#' psychologists, in which `mu` and `sigma` are the mean and SD of the
#' Gaussian component alone, and `tau` is the mean of the exponential
#' component (the tail). The mean of the full distribution is `mu + tau`.
#'
#' This distinction matters because changes in the Gaussian location (`mu`)
#' and changes in the exponential tail (`tau`) can offset one another at the
#' level of the overall mean, so effects estimated on `brms`'s default `mu`
#' can lead to different (and potentially incorrect) inferences than effects
#' estimated on this classical `mu`.
#'
#' By default, all three parameters use a `"softplus"` link
#' (`log(1 + exp(x))`), rather than `"identity"` or `"log"`. `mu` must remain
#' strictly positive since it represents the (unobserved) center of the
#' Gaussian component. An `"identity"` link would allow the linear predictor
#' to cross zero or go negative, which is invalid for all three parameters.
#' A `"log"` link would enforce positivity too, but its curvature explodes
#' as the linear predictor departs from zero, producing extreme gradients
#' and making priors/sampling harder to calibrate, especially for `mu` and
#' `tau`, which are already both on the RT scale (seconds) and can take
#' comparatively large values. `"softplus"` is positive-constrained like
#' `"log"` but behaves almost linearly (`softplus(x) ~ x`) away from zero,
#' making it easier to specify weakly-informative priors directly on the
#' RT scale while still guaranteeing valid, strictly positive parameters.
#'
#' @param link_mu,link_sigma,link_tau Character of the type of link used to
#'   model the ex-Gaussian parameters. Defaults to `"softplus"` for all three
#'   (see Details).
#'
#' @return A `brms::custom_family` object.
#' @export
rt_exgaussian <- function(
    link_mu = "softplus",
    link_sigma = "softplus",
    link_tau = "softplus"
) {
    brms::custom_family(
        name = "rt_exgaussian",
        dpars = c("mu", "sigma", "tau"), # mu/sigma = Gaussian location/SD, tau = exponential mean
        links = c(link_mu, link_sigma, link_tau),
        lb = c(0, 0, 0), # Lower bounds: mu > 0 (RT scale), sigma > 0, tau > 0
        ub = c(NA, NA, NA),
        type = "real" # Continuous outcome variable (RT)
    )
}

#' @keywords internal
.rt_exgaussian_lpdf <- function() {
    "
// Log-likelihood for a single observation from the classical Ex-Gaussian distribution.
// Y: observed reaction time.
// mu: mean of the Gaussian component (> 0).
// sigma: SD of the Gaussian component (> 0).
// tau: mean of the exponential component, i.e., the tail (> 0).
real rt_exgaussian_lpdf(real Y, real mu, real sigma, real tau) {
    // Parameter checks
    if (sigma <= 0 || tau <= 0) return negative_infinity();

    // Stan's built-in exp_mod_normal is parameterized with the rate of the
    // exponential component (beta = 1 / tau)
    return exp_mod_normal_lpdf(Y | mu, sigma, inv(tau));
}
"
}

#' @rdname rt_exgaussian
#' @export
rt_exgaussian_lpdf_expose <- function() {
    insight::check_if_installed("cmdstanr")

    # Wrap the function Stan block
    stancode <- paste0(
        "functions {
",
        .rt_exgaussian_lpdf(),
        "
}"
    )

    mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
    mod$expose_functions()
    mod$functions$rt_exgaussian_lpdf
}

#' @rdname rt_exgaussian
#' @export
rt_exgaussian_stanvars <- function() {
    brms::stanvar(scode = .rt_exgaussian_lpdf(), block = "functions")
}


# brms methods ------------------------------------------------------------

#' @rdname rt_exgaussian
#' @inheritParams lnr
#' @export
log_lik_rt_exgaussian <- function(i, prep) {
    # Extract observation
    if (!"Y" %in% names(prep$data)) {
        stop("Outcome variable 'Y' not found in prep$data.")
    }
    y <- prep$data$Y[i]
    if (is.na(y)) {
        return(NA_real_)
    }

    # Get parameters for observation i across all draws
    mu <- brms::get_dpar(prep, "mu", i = i)
    sigma <- brms::get_dpar(prep, "sigma", i = i)
    tau <- brms::get_dpar(prep, "tau", i = i)

    # Determine number of draws
    n_draws <- length(mu)
    if (n_draws == 0) {
        return(numeric(0))
    }

    # Replicate the scalar y to match the number of draws
    y_vec <- rep(y, length.out = n_draws)

    # Log-density of the classical Ex-Gaussian (Gaussian convolved with exponential),
    # computed on the log scale for numerical stability. Equivalent to Stan's
    # exp_mod_normal_lpdf(y | mu, sigma, 1 / tau).
    ll <- -log(tau) +
        (mu / tau) +
        (sigma^2) / (2 * tau^2) -
        (y_vec / tau) +
        stats::pnorm(
            y_vec,
            mean = mu + (sigma^2) / tau,
            sd = sigma,
            log.p = TRUE
        )

    # Ensure correct handling for invalid parameters
    ll[sigma <= 0 | tau <= 0] <- -Inf

    # Ensure no other NaN/NA values
    ll[is.nan(ll) | is.na(ll)] <- -Inf

    ll
}

#' @rdname rt_exgaussian
#' @inheritParams lnr
#' @export
posterior_predict_rt_exgaussian <- function(i, prep, ...) {
    # Get parameters for observation i across all draws
    mu <- brms::get_dpar(prep, "mu", i = i)
    sigma <- brms::get_dpar(prep, "sigma", i = i)
    tau <- brms::get_dpar(prep, "tau", i = i)

    # Number of posterior draws
    n_draws <- length(mu)

    # Simulate as Gaussian + Exponential (vectorized)
    stats::rnorm(n = n_draws, mean = mu, sd = sigma) +
        stats::rexp(n = n_draws, rate = 1 / tau)
}

#' @rdname rt_exgaussian
#' @export
posterior_epred_rt_exgaussian <- function(prep) {
    # Extract draws for the necessary parameters (matrices: draws x observations)
    mu <- brms::get_dpar(prep, "mu")
    tau <- brms::get_dpar(prep, "tau")

    # E[ExGaussian] = E[Gaussian] + E[Exponential] = mu + tau
    mu + tau
}
