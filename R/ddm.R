#' @title Drift Diffusion Model (DDM)
#'
#' @description
#' The Drift Diffusion Model (DDM) is a widely used model for decision-making tasks.
#' It assumes that evidence accumulates over time until it reaches one of two decision boundaries.
#'
#' The 'ddm' family in this package is not a full re-implementation of the DDM, but uses of the
#' `brms` wiener family. It is simply a reparametrized version for consistency with the other
#' race models in the package.
#'
#' @param drift Drift rate. Can take any real value.
#' @param bs Decision threshold (boundary separation). Must be positive.
#' @param bias Starting point bias (proportion of boundary separation). Must be in (0, 1).
#' @param ndt Non-decision time. Must be non-negative.
#' @param ... Other arguments to be passed to [brms::rwiener()] or [brms::dwiener()],
#' @details
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
#' # data <- rddm(1000, drift = 0.2, bs = 1, bias = 0.5, ndt = 0.3)
#' # hist(data[data$response == 0, "rt"], breaks = 50, main = "Reaction Times", xlab = "RT")
#' # hist(data[data$response == 1, "rt"], breaks = 50, add = TRUE, col = rgb(1, 0, 0, 0.5))
#'
#' # Compute density
#' # dddm(x = c(0.5, 0.7), drift = 0.2, bs = 1, bias = 0.5, resp = c(0, 1), ndt = 0.3)
#'
#' @inheritParams rlnr
#' @export
rddm <- function(
  n,
  drift,
  bs,
  bias,
  ndt,
  sigmadrift = 0,
  sigmabias = 0,
  sigmatau = 0,
  minrt = NULL,
  ...
) {
  # Prepare and validate parameters
  params <- .prepare_ddm(
    n = n,
    drift = drift,
    bs = bs,
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
  # route did). See `.ddm_draw_trialwise()` for the parameterization.
  if (.ddm_has_variability(params)) {
    params <- .ddm_draw_trialwise(params)
  }

  # Simulate data using rwiener
  sim_data <- brms::rwiener(
    params$n,
    alpha = params$bs,
    beta = params$bias,
    delta = params$drift,
    tau = params$ndt,
    ...
  )

  # Return as a data frame
  data.frame(rt = sim_data$q, response = sim_data$resp)
}

#' @rdname rddm
#' @inheritParams rlnr
#' @export
dddm <- function(
  x,
  drift,
  bs,
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
  params <- .prepare_ddm(
    x = x,
    drift = drift,
    bs = bs,
    bias = bias,
    ndt = ndt,
    response = response,
    sigmadrift = sigmadrift,
    sigmabias = sigmabias,
    sigmatau = sigmatau,
    minrt = minrt
  )

  if (.ddm_has_variability(params)) {
    dens <- .ddm_density_var(params, ...)
    return(if (log) log(dens) else dens)
  }

  # Compute density using dwiener
  brms::dwiener(
    x = params$x,
    alpha = params$bs,
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
.prepare_ddm <- function(
  n = NULL,
  x = NULL,
  drift,
  bs,
  bias,
  ndt,
  response = NULL,
  sigmadrift = 0,
  sigmabias = 0,
  sigmatau = 0,
  minrt = NULL
) {
  # --- Basic Validation ---
  if (any(bs <= 0, na.rm = TRUE)) {
    stop("bs must be positive.")
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
    # For simulation (rddm)
    if (length(n) != 1 || n <= 0 || n != floor(n)) {
      stop("n must be a single positive integer.")
    }
    m <- n # Keep n as a single value for rwiener
  } else if (!is.null(x)) {
    # For density computation (dddm)
    if (is.null(response)) {
      stop("response must be provided for dddm.")
    }
    if (any(!response %in% c(0, 1), na.rm = TRUE)) {
      stop("response must contain only 0 or 1.")
    }
    param_lengths <- c(
      length(x),
      length(drift),
      length(bs),
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
    bs = rep_len(bs, m),
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

  # Add n for rddm
  if (!is.null(n)) {
    params$n <- n # Keep n as a single value
  }

  params$ndraws <- m
  params
}

# Internal helpers for the extended (7-parameter) DDM ------------------------

#' @keywords internal
#' @noRd
.ddm_has_variability <- function(params) {
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
#' The parameterization matches the Stan likelihood in [ddm_stanvars()]:
#' `sw = sigmabias * min(2*bias, 2*(1-bias))` with the starting point centred
#' on `bias`, and `st0 = sigmatau * minrt` with `ndt` the *lower bound* of the
#' non-decision time (not its midpoint).
#'
#' @keywords internal
#' @noRd
.ddm_draw_trialwise <- function(params) {
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
.ddm_density_sv <- function(x, drift, bs, bias, ndt, response, sigmadrift, ...) {
  dt <- x - ndt
  out <- rep(0, length(dt))
  ok <- !is.na(dt) & dt > 0
  if (!any(ok)) {
    return(out)
  }

  # Density of the driftless process; the drift enters through the correction.
  out[ok] <- brms::dwiener(
    x = x[ok],
    alpha = bs[ok],
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
  num <- -v^2 * dt[ok] - 2 * v * bs[ok] * w + sv2 * bs[ok]^2 * w^2

  out[ok] <- out[ok] * exp(num / (2 * (1 + sv2 * dt[ok]))) / sqrt(1 + sv2 * dt[ok])
  out
}


#' Full (7-parameter) DDM density
#'
#' Drift variability is handled analytically by [.ddm_density_sv()]; the
#' starting point and non-decision time are integrated out by Gauss-Legendre
#' quadrature over their Uniform distributions. The non-decision time is
#' integrated only up to `x`, since later start times contribute nothing -
#' this keeps the integrand smooth and the quadrature accurate.
#'
#' Validated against the Stan likelihood in [ddm_stanvars()]: the maximum
#' relative error is at machine precision when only drift and starting-point
#' variability are present, and below 1e-5 with all three.
#'
#' @param nodes Number of quadrature nodes per integrated dimension.
#' @keywords internal
#' @noRd
.ddm_density_var <- function(params, nodes = 25, ...) {
  n <- length(params$x)
  rec <- function(v) rep_len(v, n)

  x <- rec(params$x)
  drift <- rec(params$drift)
  bs <- rec(params$bs)
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
        .ddm_density_sv(x, drift, bs, bias_a, ndt_b, response, sigmadrift, ...)
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
.ddm_to_rtdists <- function(params) {
  sw <- params$sigmabias * pmin(2 * params$bias, 2 * (1 - params$bias))
  st0 <- if (!is.null(params$minrt)) {
    params$sigmatau * params$minrt
  } else {
    rep(0, length(params$sigmatau))
  }

  list(
    a = params$bs,
    v = params$drift,
    t0 = params$ndt,
    z = params$bias * params$bs,
    sv = params$sigmadrift,
    sz = sw * params$bs,
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
#     bs = runif(1, 0.5, 2),  # Boundary separation: [0.5, 2]
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
#   rwiener(n, bs = params$bs, bias = params$bias, delta = params$drift, tau = params$ndt)
# }

# # Function to set backend and run dwiener with random parameters
# benchmark_dwiener <- function(backend) {
#   options(wiener_backend = backend)
#   params <- generate_random_params()
#   dwiener(x, bs = params$bs, bias = params$bias, delta = params$drift, tau = params$ndt, resp = resp)
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
