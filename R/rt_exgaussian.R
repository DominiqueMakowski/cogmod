#' @title Ex-Gaussian Model (Classical Parameterization)
#'
#' @description
#' Density, random generation, and `brms` custom family for the Ex-Gaussian
#' distribution, using the "classical" parameterization familiar to
#' experimental psychologists, in which `mu` and `sigma` are the mean and SD
#' of the Gaussian component alone, and `tau` is the mean of the exponential
#' component (the tail). This is unlike `brms`'s built-in `exgaussian()`
#' family, in which `mu` indexes the mean of the *entire* distribution
#' (Gaussian + exponential components combined).
#'
#' Functions:
#' - `rrt_exgaussian()`: Simulates random draws from the Ex-Gaussian distribution.
#' - `drt_exgaussian()`: Computes the density (likelihood) of the Ex-Gaussian distribution.
#' - `rt_exgaussian()`: Creates a `brms::custom_family()` for use in `brms` models.
#'
#' @details
#' The Ex-Gaussian distribution is the sum of an independent Normal (Gaussian)
#' random variable with mean `mu` and SD `sigma`, and an Exponential random
#' variable with mean `tau` (rate `1 / tau`). Unlike `brms`'s built-in
#' `exgaussian()` family - in which `mu` indexes the mean of the *entire*
#' distribution (Gaussian + exponential components combined) - here `mu` is
#' the mean of the Gaussian component alone, so that the mean of the full
#' distribution is `mu + tau`.
#'
#' This distinction matters because changes in the Gaussian location (`mu`)
#' and changes in the exponential tail (`tau`) can offset one another at the
#' level of the overall mean, so effects estimated on `brms`'s default `mu`
#' can lead to different (and potentially incorrect) inferences than effects
#' estimated on this classical `mu`.
#'
#' In the `brms` custom family (`rt_exgaussian()`), all three parameters use
#' a `"softplus"` link (`log(1 + exp(x))`) by default, rather than
#' `"identity"` or `"log"`. `mu` must remain strictly positive since it
#' represents the (unobserved) center of the Gaussian component. An
#' `"identity"` link would allow the linear predictor to cross zero or go
#' negative, which is invalid for all three parameters. A `"log"` link would
#' enforce positivity too, but its curvature explodes as the linear predictor
#' departs from zero, producing extreme gradients and making priors/sampling
#' harder to calibrate, especially for `mu` and `tau`, which are already both
#' on the RT scale (seconds) and can take comparatively large values.
#' `"softplus"` is positive-constrained like `"log"` but behaves almost
#' linearly (`softplus(x) ~ x`) away from zero, making it easier to specify
#' weakly-informative priors directly on the RT scale while still guaranteeing
#' valid, strictly positive parameters.
#'
#' @param n Number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param mu Mean of the Gaussian component. Must be positive. Range: (0, Inf).
#' @param sigma SD of the Gaussian component. Must be positive. Range: (0, Inf).
#' @param tau Mean of the exponential component (the tail). Must be positive. Range: (0, Inf).
#'
#' @references
#' - Matzke, D., & Wagenmakers, E. J. (2009). Psychological interpretation of the ex-Gaussian
#'     and shifted Wald parameters: A diffusion model analysis. *Psychonomic Bulletin & Review*,
#'     *16*(5), 798–817. \doi{10.3758/PBR.16.5.798}
#'
#' @examples
#' # Simulate 1000 RTs
#' rts <- rrt_exgaussian(1000, mu = 0.5, sigma = 0.1, tau = 0.2)
#' # hist(rts, breaks = 50, main = "Simulated Ex-Gaussian RTs", xlab = "Reaction Time")
#'
#' @export
rrt_exgaussian <- function(n, mu = 0.5, sigma = 0.1, tau = 0.2) {
    # Prepare and validate all inputs for RNG
    params <- .prepare_exgaussian(
        x = NULL,
        n = n,
        mu = mu,
        sigma = sigma,
        tau = tau
    )

    stats::rnorm(params$ndraws, mean = params$mu, sd = params$sigma) +
        stats::rexp(params$ndraws, rate = 1 / params$tau)
}


#' @rdname rrt_exgaussian
#' @param x Vector of quantiles (observed reaction times).
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @export
drt_exgaussian <- function(x, mu = 0.5, sigma = 0.1, tau = 0.2, log = FALSE) {
    # Prepare and validate inputs for density
    params <- .prepare_exgaussian(
        x = x,
        n = NULL,
        mu = mu,
        sigma = sigma,
        tau = tau
    )

    # Log-density of the classical Ex-Gaussian (Gaussian convolved with exponential),
    # computed on the log scale for numerical stability. Equivalent to Stan's
    # exp_mod_normal_lpdf(x | mu, sigma, 1 / tau).
    log_density <- -log(params$tau) +
        (params$mu / params$tau) +
        (params$sigma^2) / (2 * params$tau^2) -
        (params$x / params$tau) +
        stats::pnorm(
            params$x,
            mean = params$mu + (params$sigma^2) / params$tau,
            sd = params$sigma,
            log.p = TRUE
        )

    # Guard against numerical over/underflow
    log_density[is.nan(log_density)] <- -Inf

    if (log) {
        log_density
    } else {
        exp(log_density)
    }
}


# Internals ---------------------------------------------------------------

#' @keywords internal
.prepare_exgaussian <- function(x = NULL, n = NULL, mu, sigma, tau) {
    # Validate parameters once
    if (any(mu <= 0)) {
        stop("`mu` must be positive.")
    }
    if (any(sigma <= 0)) {
        stop("`sigma` must be positive.")
    }
    if (any(tau <= 0)) {
        stop("`tau` must be positive.")
    }

    # Determine target length:
    if (!is.null(x)) {
        # For PDF: based on x length
        m <- length(x)
    } else if (!is.null(n)) {
        # For RNG: based on n draw count
        if (length(n) != 1 || n < 0 || n != floor(n)) {
            stop("n must be a single non-negative integer.")
        }
        m <- n
    } else {
        stop("Either 'x' or 'n' must be provided.")
    }

    # Recycle vectors to length m
    params <- list(
        mu = rep(mu, length.out = m),
        sigma = rep(sigma, length.out = m),
        tau = rep(tau, length.out = m)
    )

    # Incorporate x if provided
    if (!is.null(x)) {
        params$x <- rep(x, length.out = m)
    }
    params$ndraws <- m
    params
}
