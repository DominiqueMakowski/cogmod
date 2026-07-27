#' @rdname rrt_exgaussian
#' @param link_mu,link_sigma,link_tau Character of the type of link used to
#'   model the ex-Gaussian parameters. Defaults to `"softplus"` for all three
#'   (see Details).
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

#' @rdname rrt_exgaussian
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

#' @rdname rrt_exgaussian
#' @export
rt_exgaussian_stanvars <- function() {
    brms::stanvar(scode = .rt_exgaussian_lpdf(), block = "functions")
}


# brms methods ------------------------------------------------------------

#' @rdname rrt_exgaussian
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

#' @rdname rrt_exgaussian
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

#' @rdname rrt_exgaussian
#' @export
posterior_epred_rt_exgaussian <- function(prep) {
    # Extract draws for the necessary parameters (matrices: draws x observations)
    mu <- brms::get_dpar(prep, "mu")
    tau <- brms::get_dpar(prep, "tau")

    # E[ExGaussian] = E[Gaussian] + E[Exponential] = mu + tau
    mu + tau
}
