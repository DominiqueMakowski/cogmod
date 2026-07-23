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
#' @param backend The backend to use for the simulation. Based on local benchmarks, `"rtdists"`
#'   is much faster for simulation but `"Rwiener"` is faster for density computation. Note that
#'   the `"Rwiener"` backend (used by default for the density, e.g., in `log_lik()`/`loo()`)
#'   requires the `RWiener` package to be installed. Ignored (forced to `rtdists`, via
#'   [rtdists::rdiffusion()]/[rtdists::ddiffusion()]) whenever `sigmadrift`, `sigmabias`, or
#'   `sigmatau` is non-zero, since only `rtdists` implements the full 7-parameter model.
#' @param sigmadrift Inter-trial variability in drift rate (`sv` in `rtdists` terms). Must be
#'   non-negative. Default `0` (no variability, i.e. the classic 4-parameter DDM).
#' @param sigmabias Inter-trial variability in the starting point, expressed as a fraction (in
#'   `[0, 1)`) of the maximum allowed range, i.e. `sw = sigmabias * min(2*bias, 2*(1-bias))`
#'   (`sw`/`sz` in `rtdists` terms). Default `0`.
#' @param sigmatau Inter-trial variability in non-decision time, expressed as a fraction of
#'   `minrt`, i.e. `st0 = sigmatau * minrt` (`st0` in `rtdists` terms). Default `0`.
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
  backend = "rtdists",
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

  if (.ddm_has_variability(params)) {
    insight::check_if_installed("rtdists")
    rt_params <- .ddm_to_rtdists(params)

    sim_data <- rtdists::rdiffusion(
      params$n,
      a = rt_params$a,
      v = rt_params$v,
      t0 = rt_params$t0,
      z = rt_params$z,
      sv = rt_params$sv,
      sz = rt_params$sz,
      st0 = rt_params$st0
    )

    return(data.frame(
      rt = sim_data$rt,
      response = ifelse(sim_data$response == "upper", 1, 0)
    ))
  }

  # Simulate data using rwiener
  sim_data <- brms::rwiener(
    params$n,
    alpha = params$bs,
    beta = params$bias,
    delta = params$drift,
    tau = params$ndt,
    backend = backend,
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
  backend = "Rwiener",
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
    insight::check_if_installed("rtdists")
    rt_params <- .ddm_to_rtdists(params)

    dens <- rtdists::ddiffusion(
      rt = params$x,
      response = params$response + 1L, # 0/1 -> 1/2 (lower/upper)
      a = rt_params$a,
      v = rt_params$v,
      t0 = rt_params$t0,
      z = rt_params$z,
      sv = rt_params$sv,
      sz = rt_params$sz,
      st0 = rt_params$st0
    )

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
    backend = backend,
    ...
  )
}


# Internals ---------------------------------------------------------------

#' @keywords internal
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

# Internal helpers for the extended (7-parameter) DDM via `rtdists` ---------

#' @keywords internal
.ddm_has_variability <- function(params) {
  any(params$sigmadrift != 0, na.rm = TRUE) ||
    any(params$sigmabias != 0, na.rm = TRUE) ||
    any(params$sigmatau != 0, na.rm = TRUE)
}

#' @keywords internal
.ddm_to_rtdists <- function(params) {
  # sw: our Stan/rtdists-relative starting-point variability (fraction of the
  # maximum allowed range), converted to rtdists' absolute (a-scaled) `sz`.
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
