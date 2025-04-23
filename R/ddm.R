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
#'   is much faster for simulation but `"Rwiener"` is faster for density computation.
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
rddm <- function(n, drift, bs, bias, ndt, backend = "rtdists", ...) {
  # Prepare and validate parameters
  params <- .prepare_ddm(n = n, drift = drift, bs = bs, bias = bias, ndt = ndt)

  # Simulate data using rwiener
  sim_data <- brms::rwiener(params$n, alpha = params$bs, beta = params$bias,
                            delta = params$drift, tau = params$ndt, backend = backend, ...)

  # Return as a data frame
  data.frame(rt = sim_data$q, response = sim_data$resp)
}

#' @rdname rddm
#' @inheritParams rlnr
#' @export
dddm <- function(x, drift, bs, bias, ndt, response, log = FALSE, backend = "Rwiener", ...) {
  # Prepare and validate parameters
  params <- .prepare_ddm(x = x, drift = drift, bs = bs, bias = bias, ndt = ndt, response = response)

  # Compute density using dwiener
  brms::dwiener(x = params$x, alpha = params$bs, beta = params$bias,
                delta = params$drift, resp = params$response, tau = params$ndt, log = log, backend = backend, ...)
}


# Internals ---------------------------------------------------------------

#' @keywords internal
.prepare_ddm <- function(n = NULL, x = NULL, drift, bs, bias, ndt, response = NULL) {
  # --- Basic Validation ---
  if (any(bs <= 0, na.rm = TRUE)) stop("bs must be positive.")
  if (any(bias <= 0 | bias >= 1, na.rm = TRUE)) stop("bias must be in (0, 1).")
  if (any(ndt < 0, na.rm = TRUE)) stop("ndt must be non-negative.")

  # --- Determine Target Length ---
  if (!is.null(n)) {
    # For simulation (rddm)
    if (length(n) != 1 || n <= 0 || n != floor(n)) {
      stop("n must be a single positive integer.")
    }
    m <- n  # Keep n as a single value for rwiener
  } else if (!is.null(x)) {
    # For density computation (dddm)
    if (is.null(response)) stop("response must be provided for dddm.")
    if (any(!response %in% c(0, 1), na.rm = TRUE)) stop("response must contain only 0 or 1.")
    param_lengths <- c(length(x), length(drift), length(bs), length(bias), length(ndt), length(response))
    m <- max(param_lengths)
  } else {
    stop("Internal error: Either 'n' or 'x' must be provided.")
  }

  # --- Recycle Parameters ---
  params <- list(
    drift    = rep_len(drift, m),
    bs = rep_len(bs, m),
    bias  = rep_len(bias, m),
    ndt   = rep_len(ndt, m)
  )

  # Add x and response if provided
  if (!is.null(x)) {
    params$x <- rep_len(x, m)
  }
  if (!is.null(response)) {
    params$response <- rep_len(response, m)
  }

  # Add n for rddm
  if (!is.null(n)) {
    params$n <- n  # Keep n as a single value
  }

  params$ndraws <- m
  params
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
