#' @title Drift Diffusion Model (DDM)
#'
#' @description
#' The Drift Diffusion Model (DDM) is a widely used model for decision-making tasks.
#' It assumes that evidence accumulates over time until it reaches one of two decision boundaries.
#'
#' The 'cogmod_ddm' family in this package is not a full re-implementation of the DDM, but uses of the
#' `brms` wiener family. It is simply a reparametrized version for consistency with the other
#' race models in the package.
#'
#' @param drift Drift rate. Can take any real value.
#' @param boundary Decision threshold (boundary separation). Must be positive.
#' @param bias Starting point bias, as a proportion of the boundary separation
#'   *measured from the boundary coded `0`*. Must be in (0, 1). See Details.
#' @param ndt Non-decision time. Must be non-negative.
#' @param ... Other arguments to be passed to [brms::rwiener()] or [brms::dwiener()],
#' @details
#' # Response coding
#'
#' The response coded `1` (`response` here, `dec()` in a `brms` formula) is the
#' **upper** boundary and the response coded `0` is the **lower** one, following
#' `brms`'s own `wiener()` family. Two consequences are worth keeping in mind
#' when reading a fitted model:
#'
#' * `bias` is measured *from the lower boundary*, so `bias > 0.5` places the
#'   starting point closer to the response coded `1`, and `bias < 0.5` closer to
#'   the response coded `0`. Since first-passage times are shorter for the nearer
#'   boundary, `bias > 0.5` makes responses coded `1` *faster* than responses
#'   coded `0`, and `bias < 0.5` makes them slower. At exactly `bias = 0.5` the
#'   two conditional RT distributions are identical, whatever the drift rate.
#' * A positive `drift` pushes the accumulator towards the response coded `1`.
#'   When `0` codes correct responses and `1` codes errors (a common choice when
#'   the task has no natural stimulus-to-boundary mapping), good performance
#'   therefore corresponds to a *negative* drift rate.
#'
#' # Implementation
#'
#' The underlying 4-parameter process is simulated and evaluated with
#' [brms::rwiener()]/[brms::dwiener()] (which require the `RWiener` package).
#' The full 7-parameter model is built on top of these rather than delegated to
#' another package: between-trial variability is simulated by drawing the
#' per-trial parameters, and evaluated by combining a closed-form drift
#' correction with Gauss-Legendre quadrature over the starting point and
#' non-decision time.
#' @param sigmadrift Inter-trial variability in drift rate (`sv` in the usual notation). Must be
#'   non-negative. Default `0` (no variability, i.e. the classic 4-parameter DDM).
#' @param sigmabias Inter-trial variability in the starting point, expressed as a fraction (in
#'   `[0, 1)`) of the maximum allowed range, i.e. `sw = sigmabias * min(2*bias, 2*(1-bias))`
#'   (`sw`/`sz` in the usual notation). Default `0`.
#' @param sigmatau Inter-trial variability in non-decision time, expressed as a fraction of
#'   `minrt`, i.e. `st0 = sigmatau * minrt` (`st0` in the usual notation), with `ndt` the lower
#'   bound of the resulting Uniform distribution. Default `0`.
#' @param minrt Minimum reaction time. Only required when `sigmatau > 0` (used to scale
#'   `sigmatau` into `st0`).
#'
#' @examples
#' # Simulate data
#' # data <- rcogmod_ddm(1000, drift = 0.2, boundary = 1, bias = 0.5, ndt = 0.3)
#' # hist(data[data$response == 0, "rt"], breaks = 50, main = "Reaction Times", xlab = "RT")
#' # hist(data[data$response == 1, "rt"], breaks = 50, add = TRUE, col = rgb(1, 0, 0, 0.5))
#'
#' # Compute density
#' # dcogmod_ddm(x = c(0.5, 0.7), drift = 0.2, boundary = 1, bias = 0.5, resp = c(0, 1), ndt = 0.3)
#'
#' @inheritParams rcogmod_lnr
#' @export
rcogmod_ddm <- function(
  n,
  drift,
  boundary,
  bias,
  ndt,
  sigmadrift = 0,
  sigmabias = 0,
  sigmatau = 0,
  minrt = NULL,
  ...
) {
  # Prepare and validate parameters
  params <- .prepare_cogmod_ddm(
    n = n,
    drift = drift,
    boundary = boundary,
    bias = bias,
    ndt = ndt,
    sigmadrift = sigmadrift,
    sigmabias = sigmabias,
    sigmatau = sigmatau,
    minrt = minrt
  )

  # Between-trial variability needs no special sampler: it *is* a hierarchical
  # draw. Sampling the per-trial drift, starting point and non-decision time and
  # then running the ordinary 4-parameter process is exact, and vastly cheaper
  # than inverting a numerically integrated CDF (which is what the `rtdists`
  # route did). See `.cogmod_ddm_draw_trialwise()` for the parameterization.
  if (.cogmod_ddm_has_variability(params)) {
    params <- .cogmod_ddm_draw_trialwise(params)
  }

  # Simulate data using rwiener
  sim_data <- brms::rwiener(
    params$n,
    alpha = params$boundary,
    beta = params$bias,
    delta = params$drift,
    tau = params$ndt,
    ...
  )

  # Return as a data frame
  data.frame(rt = sim_data$q, response = sim_data$resp)
}

#' @rdname rcogmod_ddm
#' @inheritParams rcogmod_lnr
#' @export
dcogmod_ddm <- function(
  x,
  drift,
  boundary,
  bias,
  ndt,
  response,
  log = FALSE,
  sigmadrift = 0,
  sigmabias = 0,
  sigmatau = 0,
  minrt = NULL,
  ...
) {
  # Prepare and validate parameters
  params <- .prepare_cogmod_ddm(
    x = x,
    drift = drift,
    boundary = boundary,
    bias = bias,
    ndt = ndt,
    response = response,
    sigmadrift = sigmadrift,
    sigmabias = sigmabias,
    sigmatau = sigmatau,
    minrt = minrt
  )

  if (.cogmod_ddm_has_variability(params)) {
    dens <- .cogmod_ddm_density_var(params, ...)
    return(if (log) log(dens) else dens)
  }

  # Compute density using dwiener
  brms::dwiener(
    x = params$x,
    alpha = params$boundary,
    beta = params$bias,
    delta = params$drift,
    resp = params$response,
    tau = params$ndt,
    log = log,
    ...
  )
}


# Internals ---------------------------------------------------------------

#' @keywords internal
#' @noRd
.prepare_cogmod_ddm <- function(
  n = NULL,
  x = NULL,
  drift,
  boundary,
  bias,
  ndt,
  response = NULL,
  sigmadrift = 0,
  sigmabias = 0,
  sigmatau = 0,
  minrt = NULL
) {
  # --- Basic Validation ---
  if (any(boundary <= 0, na.rm = TRUE)) {
    stop("boundary must be positive.")
  }
  if (any(bias <= 0 | bias >= 1, na.rm = TRUE)) {
    stop("bias must be in (0, 1).")
  }
  if (any(ndt < 0, na.rm = TRUE)) {
    stop("ndt must be non-negative.")
  }
  if (any(sigmadrift < 0, na.rm = TRUE)) {
    stop("sigmadrift must be non-negative.")
  }
  if (any(sigmabias < 0 | sigmabias >= 1, na.rm = TRUE)) {
    stop("sigmabias must be in [0, 1).")
  }
  if (any(sigmatau < 0, na.rm = TRUE)) {
    stop("sigmatau must be non-negative.")
  }
  if (any(sigmatau > 0, na.rm = TRUE) && is.null(minrt)) {
    stop("minrt must be provided when sigmatau > 0.")
  }
  if (!is.null(minrt) && any(minrt <= 0, na.rm = TRUE)) {
    stop("minrt must be positive.")
  }

  # --- Determine Target Length ---
  if (!is.null(n)) {
    # For simulation (rcogmod_ddm)
    if (length(n) != 1 || n <= 0 || n != floor(n)) {
      stop("n must be a single positive integer.")
    }
    m <- n # Keep n as a single value for rwiener
  } else if (!is.null(x)) {
    # For density computation (dcogmod_ddm)
    if (is.null(response)) {
      stop("response must be provided for dcogmod_ddm.")
    }
    if (any(!response %in% c(0, 1), na.rm = TRUE)) {
      stop("response must contain only 0 or 1.")
    }
    param_lengths <- c(
      length(x),
      length(drift),
      length(boundary),
      length(bias),
      length(ndt),
      length(response),
      length(sigmadrift),
      length(sigmabias),
      length(sigmatau)
    )
    m <- max(param_lengths)
  } else {
    stop("Internal error: Either 'n' or 'x' must be provided.")
  }

  # --- Recycle Parameters ---
  params <- list(
    drift = rep_len(drift, m),
    boundary = rep_len(boundary, m),
    bias = rep_len(bias, m),
    ndt = rep_len(ndt, m),
    sigmadrift = rep_len(sigmadrift, m),
    sigmabias = rep_len(sigmabias, m),
    sigmatau = rep_len(sigmatau, m)
  )

  # Add x and response if provided
  if (!is.null(x)) {
    params$x <- rep_len(x, m)
  }
  if (!is.null(response)) {
    params$response <- rep_len(response, m)
  }
  if (!is.null(minrt)) {
    params$minrt <- rep_len(minrt, m)
  }

  # Add n for rcogmod_ddm
  if (!is.null(n)) {
    params$n <- n # Keep n as a single value
  }

  params$ndraws <- m
  params
}

# Internal helpers for the extended (7-parameter) DDM ------------------------

#' @keywords internal
#' @noRd
.cogmod_ddm_has_variability <- function(params) {
  any(params$sigmadrift != 0, na.rm = TRUE) ||
    any(params$sigmabias != 0, na.rm = TRUE) ||
    any(params$sigmatau != 0, na.rm = TRUE)
}

#' Draw per-trial DDM parameters from their between-trial distributions
#'
#' Between-trial variability is, by definition, a hierarchical draw: the drift
#' is Normal, the starting point Uniform and the non-decision time Uniform
#' across trials. Drawing them and then running the plain 4-parameter process
#' reproduces the 7-parameter model exactly, so no dedicated sampler is needed.
#'
#' The parameterization matches the Stan likelihood in [cogmod_ddm_stanvars()]:
#' `sw = sigmabias * min(2*bias, 2*(1-bias))` with the starting point centred
#' on `bias`, and `st0 = sigmatau * minrt` with `ndt` the *lower bound* of the
#' non-decision time (not its midpoint).
#'
#' @keywords internal
#' @noRd
.cogmod_ddm_draw_trialwise <- function(params) {
  n <- params$n

  params$drift <- stats::rnorm(n, params$drift, params$sigmadrift)

  sw <- params$sigmabias * pmin(2 * params$bias, 2 * (1 - params$bias))
  params$bias <- stats::runif(n, params$bias - sw / 2, params$bias + sw / 2)

  st0 <- if (is.null(params$minrt)) rep(0, n) else params$sigmatau * params$minrt
  params$ndt <- pmax(stats::runif(n, params$ndt, params$ndt + st0), 0)

  params
}


#' Gauss-Legendre nodes and weights on `[-1, 1]`
#'
#' Computed with the Golub-Welsch algorithm (eigendecomposition of the Jacobi
#' matrix), so that no additional dependency is needed.
#'
#' @keywords internal
#' @noRd
.gauss_legendre <- function(n) {
  i <- seq_len(n - 1)
  b <- i / sqrt(4 * i^2 - 1)
  j <- matrix(0, n, n)
  j[cbind(i, i + 1)] <- b
  j[cbind(i + 1, i)] <- b
  e <- eigen(j, symmetric = TRUE)
  list(nodes = rev(e$values), weights = rev(2 * e$vectors[1, ]^2))
}


#' DDM density with between-trial variability in drift rate
#'
#' Integrating the 4-parameter Wiener density over a Normal drift has a closed
#' form: the density with zero drift, times a correction factor. No numerical
#' integration is required, so this costs one [brms::dwiener()] call.
#'
#' @keywords internal
#' @noRd
.cogmod_ddm_density_sv <- function(x, drift, boundary, bias, ndt, response, sigmadrift, ...) {
  dt <- x - ndt
  out <- rep(0, length(dt))
  ok <- !is.na(dt) & dt > 0
  if (!any(ok)) {
    return(out)
  }

  # Density of the driftless process; the drift enters through the correction.
  out[ok] <- brms::dwiener(
    x = x[ok],
    alpha = boundary[ok],
    beta = bias[ok],
    delta = 0,
    resp = response[ok],
    tau = ndt[ok],
    ...
  )

  # The upper boundary is the lower one with the drift and start point flipped.
  v <- ifelse(response[ok] == 1, -drift[ok], drift[ok])
  w <- ifelse(response[ok] == 1, 1 - bias[ok], bias[ok])
  sv2 <- sigmadrift[ok]^2
  num <- -v^2 * dt[ok] - 2 * v * boundary[ok] * w + sv2 * boundary[ok]^2 * w^2

  out[ok] <- out[ok] * exp(num / (2 * (1 + sv2 * dt[ok]))) / sqrt(1 + sv2 * dt[ok])
  out
}


#' Full (7-parameter) DDM density
#'
#' Drift variability is handled analytically by `.cogmod_ddm_density_sv()`; the
#' starting point and non-decision time are integrated out by Gauss-Legendre
#' quadrature over their Uniform distributions. The non-decision time is
#' integrated only up to `x`, since later start times contribute nothing -
#' this keeps the integrand smooth and the quadrature accurate.
#'
#' Validated against the Stan likelihood in [cogmod_ddm_stanvars()]: the maximum
#' relative error is at machine precision when only drift and starting-point
#' variability are present, and below 1e-5 with all three.
#'
#' @param nodes Number of quadrature nodes per integrated dimension.
#' @keywords internal
#' @noRd
.cogmod_ddm_density_var <- function(params, nodes = 25, ...) {
  n <- length(params$x)
  rec <- function(v) rep_len(v, n)

  x <- rec(params$x)
  drift <- rec(params$drift)
  boundary <- rec(params$boundary)
  bias <- rec(params$bias)
  ndt <- rec(params$ndt)
  response <- rec(params$response)
  sigmadrift <- rec(params$sigmadrift)

  sw <- rec(params$sigmabias) * pmin(2 * bias, 2 * (1 - bias))
  st0 <- if (is.null(params$minrt)) rep(0, n) else rec(params$sigmatau) * rec(params$minrt)

  # A degenerate "node" of weight 2 reproduces the point value once the 1/2
  # factor of the uniform average is applied.
  degenerate <- list(nodes = 0, weights = 2)
  quad <- .gauss_legendre(nodes)
  w_quad <- if (any(sw > 0)) quad else degenerate
  t_quad <- if (any(st0 > 0)) quad else degenerate

  span <- pmax(pmin(ndt + st0, x) - ndt, 0)
  scale <- ifelse(st0 > 0, span / (2 * st0), 1 / 2)

  out <- numeric(n)
  for (a in seq_along(w_quad$nodes)) {
    bias_a <- bias + (sw / 2) * w_quad$nodes[a]
    inner <- numeric(n)
    for (b in seq_along(t_quad$nodes)) {
      ndt_b <- ndt + (span / 2) * (t_quad$nodes[b] + 1)
      inner <- inner + t_quad$weights[b] *
        .cogmod_ddm_density_sv(x, drift, boundary, bias_a, ndt_b, response, sigmadrift, ...)
    }
    out <- out + (w_quad$weights[a] / 2) * inner * scale
  }
  out
}


#' Translate `cogmod`'s parameterization into `rtdists`' one
#'
#' Not used by the package itself - the 7-parameter model is implemented
#' natively - but kept as the reference mapping for cross-checking against
#' [rtdists::ddiffusion()]/[rtdists::rdiffusion()]. Note that `rtdists` takes
#' the starting point and its variability on the absolute (`a`-scaled) scale,
#' and treats `t0` as the *lower* bound of the non-decision time distribution.
#'
#' @keywords internal
#' @noRd
.cogmod_ddm_to_rtdists <- function(params) {
  sw <- params$sigmabias * pmin(2 * params$bias, 2 * (1 - params$bias))
  st0 <- if (!is.null(params$minrt)) {
    params$sigmatau * params$minrt
  } else {
    rep(0, length(params$sigmatau))
  }

  list(
    a = params$boundary,
    v = params$drift,
    t0 = params$ndt,
    z = params$bias * params$boundary,
    sv = params$sigmadrift,
    sz = sw * params$boundary,
    st0 = st0
  )
}

# CODE to BENCHMARK the backends

# # Load required packages
# library(brms)
# library(microbenchmark)

# # Function to randomly generate parameters
# generate_random_params <- function() {
#   list(
#     boundary = runif(1, 0.5, 2),  # Boundary separation: [0.5, 2]
#     bias = runif(1, 0.3, 0.7), # Starting point bias: [0.3, 0.7]
#     drift = runif(1, -1, 1),      # Drift rate: [-1, 1]
#     ndt = runif(1, 0.1, 0.5)   # Non-decision time: [0.1, 0.5]
#   )
# }

# # Benchmark settings
# n <- 5000  # Number of trials for simulation
# x <- seq(0.05, 2.5, length.out = 1000)  # Reaction times for density computation
# resp <- sample(c(1, 2), size = length(x), replace = TRUE)  # Random responses

# # Function to set backend and run rwiener with random parameters
# benchmark_rwiener <- function(backend) {
#   options(wiener_backend = backend)
#   params <- generate_random_params()
#   rwiener(n, boundary = params$boundary, bias = params$bias, delta = params$drift, tau = params$ndt)
# }

# # Function to set backend and run dwiener with random parameters
# benchmark_dwiener <- function(backend) {
#   options(wiener_backend = backend)
#   params <- generate_random_params()
#   dwiener(x, boundary = params$boundary, bias = params$bias, delta = params$drift, tau = params$ndt, resp = resp)
# }

# # Benchmark rwiener for both backends
# rwiener_benchmark <- microbenchmark(
#   Rwiener = benchmark_rwiener("Rwiener"),
#   rtdists = benchmark_rwiener("rtdists"),
#   times = 100
# )
# print("rwiener Benchmark:")
# print(rwiener_benchmark)
# autoplot(rwiener_benchmark) + ggtitle("rwiener Benchmark")

# # Benchmark dwiener for both backends
# dwiener_benchmark <- microbenchmark(
#   Rwiener = benchmark_dwiener("Rwiener"),
#   rtdists = benchmark_dwiener("rtdists"),
#   times = 100
# )

# # Print results
# print("dwiener Benchmark:")
# print(dwiener_benchmark)
# autoplot(dwiener_benchmark) + ggtitle("dwiener Benchmark")


#' @keywords internal
.cogmod_ddm_lpdf <- function() {
  "
// Reparameterized Wiener diffusion log-PDF for a single response, using
// Stan's `wiener_lpdf()` (7-parameter form, aka `wiener_full_lpdf()`).
// mu: drift rate parameter
// boundary: boundary separation parameter > 0
// bias: initial bias parameter in [0, 1]
// tau: non-decision time proportion in [0, 1]
// minrt: minimum reaction time > 0
// sigmadrift: inter-trial variability in drift (sv), >= 0
// sigmabias: fraction (in [0, 1)) of the maximum allowed starting-point range;
//            sw = sigmabias * fmin(2*bias, 2*(1-bias)) keeps bias +/- sw/2 inside (0, 1)
// sigmatau: scales minrt to yield the non-decision time range st0 (convenience only,
//           not a hard constraint on validity)
//
// Performance note: Stan's own `sw == 0 && st0 == 0` fast path inside
// `wiener_lpdf()`/`wiener_full_lpdf()` still delegates to the newer 'wiener5'
// algorithm (which also handles sv analytically), which is a structurally
// different and empirically much costlier implementation than the classic
// 4-parameter `wiener_lpdf()`, even when sv = 0. We therefore short-circuit
// ourselves whenever sw and st0 both vanish: to the classic 4-parameter
// density when sv is also 0 (fastest path, identical to the original
// `cogmod_ddm()`), or to the dedicated 5-parameter (sv-only) form otherwise, which
// is still much cheaper than the general 7-parameter form (the latter falls
// back to adaptive numerical quadrature whenever sw or st0 is nonzero).
real cogmod_ddm_lpdf(real y, real mu, real boundary, real bias,
              real tau, real minrt,
              real sigmadrift, real sigmabias, real sigmatau,
              int dec) {
    real ndt = tau * minrt;
    real sw  = sigmabias * fmin(2 * bias, 2 * (1 - bias));
    real st0 = sigmatau * minrt;

    if (sw == 0 && st0 == 0) {
        if (sigmadrift == 0) {
            // Fastest path: classic 4-parameter Navarro & Fuss density.
            if (dec == 1) {
                return wiener_lpdf(y | boundary, ndt, bias, mu);
            } else {
                return wiener_lpdf(y | boundary, ndt, 1 - bias, -mu);
            }
        } else {
            // Drift-variability-only: dedicated 5-parameter density.
            if (dec == 1) {
                return wiener_lpdf(y | boundary, ndt, bias, mu, sigmadrift);
            } else {
                return wiener_lpdf(y | boundary, ndt, 1 - bias, -mu, sigmadrift);
            }
        }
    }

    // General case: full 7-parameter Wiener density (adaptive quadrature).
    if (dec == 1) {
        return wiener_lpdf(y | boundary, ndt, bias, mu, sigmadrift, sw, st0);
    } else {
        return wiener_lpdf(y | boundary, ndt, 1 - bias, -mu, sigmadrift, sw, st0);
    }
}
"
}


#' @rdname rcogmod_ddm
#' @export
cogmod_ddm_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Wrap the function Stan block (done normally by brms on model compilation)
  stancode <- paste0(
    "functions {
",
    .cogmod_ddm_lpdf(),
    "
}"
  )

  mod <- cmdstanr::cmdstan_model(
    cmdstanr::write_stan_file(stancode),
    force_recompile = TRUE
  )
  mod$expose_functions()
  mod$functions$cogmod_ddm_lpdf
}


#' @rdname rcogmod_ddm
#' @export
cogmod_ddm_stanvars <- function() {
  brms::stanvar(scode = .cogmod_ddm_lpdf(), block = "functions")
}


#' @rdname rcogmod_ddm
#' @param link_mu,link_boundary,link_bias,link_tau,link_minrt Link functions for the parameters.
#' @param link_sigmadrift,link_sigmabias,link_sigmatau Link functions for the extra
#'   inter-trial variability parameters (fix them to 0, e.g. `sigmadrift = 0`, in the
#'   `brms::bf()` formula to recover the classic 4-parameter DDM).
#' @export
cogmod_ddm <- function(
  link_mu = "identity",
  link_boundary = "softplus",
  link_bias = "logit",
  link_tau = "logit",
  link_minrt = "identity",
  link_sigmadrift = "softplus",
  link_sigmabias = "logit",
  link_sigmatau = "logit"
) {
  brms::custom_family(
    name = "cogmod_ddm",
    dpars = c(
      "mu",
      "boundary",
      "bias",
      "tau",
      "minrt",
      "sigmadrift",
      "sigmabias",
      "sigmatau"
    ),
    links = c(
      link_mu,
      link_boundary,
      link_bias,
      link_tau,
      link_minrt,
      link_sigmadrift,
      link_sigmabias,
      link_sigmatau
    ),
    lb = c(NA, 0, 0, 0, 0, 0, 0, 0),
    ub = c(NA, NA, 1, 1, NA, NA, 1, 1),
    type = "real", # Response variable type
    vars = "dec[n]" # Required additional data variable for the decision
  )
}


# brms --------------------------------------------------------------------

#' @rdname rcogmod_ddm
#' @inheritParams rcogmod_lnr
#' @export
log_lik_cogmod_ddm <- function(i, prep, ...) {
  # Extract parameters for observation i across all posterior draws (vectors)
  drift <- brms::get_dpar(prep, "mu", i = i)
  boundary <- brms::get_dpar(prep, "boundary", i = i)
  bias <- brms::get_dpar(prep, "bias", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i) # Vector if varying, scalar if fixed
  sigmadrift <- brms::get_dpar(prep, "sigmadrift", i = i)
  sigmabias <- brms::get_dpar(prep, "sigmabias", i = i)
  sigmatau <- brms::get_dpar(prep, "sigmatau", i = i)

  # Extract observed data for observation i (scalars)
  y <- prep$data$Y[i] # Observed RT
  dec <- prep$data$dec[i] # Observed decision

  # Handle missing data for the observation
  if (is.na(y) || is.na(dec)) {
    # Return vector of NAs, one for each draw
    return(rep(NA_real_, prep$ndraws))
  }

  # Basic check for valid response coding (0 or 1)
  if (!dec %in% c(0, 1)) {
    warning(
      "Response ('dec') must be 0 or 1. Found: ",
      dec,
      " for observation ",
      i,
      ". Returning -Inf log-likelihood."
    )
    # Return -Inf for all draws for this observation
    return(rep(-Inf, prep$ndraws))
  }

  # Calculate ndt (vector, one element per draw)
  # Element-wise multiplication handles if minrt is scalar or vector
  ndt <- tau * minrt

  # Compute log-likelihood using dcogmod_ddm(), which automatically dispatches to
  # the full 7-parameter density whenever sigmadrift, sigmabias, or sigmatau
  # is non-zero, and to the fast 4-parameter brms::dwiener() otherwise.
  ll <- dcogmod_ddm(
    x = y,
    drift = drift,
    boundary = boundary,
    bias = bias,
    ndt = ndt,
    response = dec,
    log = TRUE,
    sigmadrift = sigmadrift,
    sigmabias = sigmabias,
    sigmatau = sigmatau,
    minrt = minrt,
    ...
  )

  # Return vector of log-likelihoods (one per posterior draw)
  ll
}


#' @rdname rcogmod_ddm
#' @inheritParams rcogmod_lnr
#' @export
posterior_epred_cogmod_ddm <- function(prep) {
  # Extract parameters
  mu <- prep$dpars$mu
  boundary <- prep$dpars$boundary
  bias <- prep$dpars$bias
  tau <- prep$dpars$tau
  minrt <- prep$dpars$minrt
  sigmadrift <- prep$dpars$sigmadrift
  sigmabias <- prep$dpars$sigmabias
  sigmatau <- prep$dpars$sigmatau

  # The closed-form E[RT] formula below is only exact for the classic
  # 4-parameter DDM. Warn (once) when inter-trial variability is present,
  # since there is no simple closed form for E[RT] in that case: users who
  # need exact predictions should use posterior_predict() and average the
  # simulated draws instead.
  if (any(sigmadrift != 0) || any(sigmabias != 0) || any(sigmatau != 0)) {
    warning(
      "posterior_epred() for `cogmod_ddm()` uses a closed-form approximation that is only exact ",
      "for the classic 4-parameter DDM (sigmadrift = sigmabias = sigmatau = 0). Since ",
      "inter-trial variability parameters are non-zero here, these expected-RT predictions ",
      "are only approximate. For exact predictions, use `posterior_predict()` and average ",
      "the simulated draws instead.",
      call. = FALSE
    )
  }

  # Compute ndt from tau and minrt
  ndt <- tau * minrt

  # Compute expected reaction time (E[RT]), i.e. the *unconditional* mean
  # first-passage time plus the non-decision time.
  # Formula adapted from https://doi.org/10.1016/j.jmp.2009.01.006, which takes
  # the starting point on the absolute scale: `bias` is a proportion of `boundary`,
  # so the quantity entering the formula is z = bias * boundary, not `bias` itself.
  z <- bias * boundary
  epred <- ndt - z / mu + boundary / mu * (exp(-2 * mu * z) - 1) / (exp(-2 * mu * boundary) - 1)

  # The expression above has a removable singularity at mu = 0 (and loses
  # precision to cancellation near it). The limit is ndt + z * (boundary - z).
  near_zero <- abs(mu) < 1e-6
  if (any(near_zero)) {
    # `+ 0 * mu` recycles the limit to the shape of the draws x observations
    # matrix, whatever mix of scalar and matrix dpars was supplied.
    lim <- ndt + z * (boundary - z) + 0 * mu
    epred[near_zero] <- lim[near_zero]
  }

  # Return expected reaction time
  epred
}


#' @rdname rcogmod_ddm
#' @inheritParams rcogmod_lnr
#' @export
posterior_predict_cogmod_ddm <- function(i, prep, ...) {
  # Extract parameters for observation i
  mu <- brms::get_dpar(prep, "mu", i = i)
  boundary <- brms::get_dpar(prep, "boundary", i = i)
  bias <- brms::get_dpar(prep, "bias", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)
  sigmadrift <- brms::get_dpar(prep, "sigmadrift", i = i)
  sigmabias <- brms::get_dpar(prep, "sigmabias", i = i)
  sigmatau <- brms::get_dpar(prep, "sigmatau", i = i)

  # Compute ndt from tau and minrt
  ndt <- tau * minrt

  # Simulate reaction times and responses using rcogmod_ddm(), which automatically
  # dispatches to the full 7-parameter simulation whenever sigmadrift,
  # sigmabias, or sigmatau is non-zero.
  sim_data <- rcogmod_ddm(
    n = prep$ndraws,
    drift = mu,
    boundary = boundary,
    bias = bias,
    ndt = ndt,
    sigmadrift = sigmadrift,
    sigmabias = sigmabias,
    sigmatau = sigmatau,
    minrt = minrt
  )

  # Return simulated reaction times
  as.matrix(sim_data)
}
