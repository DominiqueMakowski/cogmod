library(testthat)

test_that("rlnr generates valid outcomes", {
  set.seed(123)
  n <- 1000
  data <- rlnr(n, mu = 1, mudelta = 0.5, sigmazero = 1, sigmadelta = -0.5, tau = 0.2)

  # Check that the output is a data frame with the correct columns
  expect_true(is.data.frame(data))
  expect_true(all(c("rt", "response") %in% colnames(data)))

  # Check that reaction times are positive
  expect_true(all(data$rt > 0))

  # Check that responses are either 0 or 1
  expect_true(all(data$response %in% c(0, 1)))
})

test_that("dlnr computes valid densities", {
  set.seed(123)
  y <- 1.5
  mu <- 1
  mudelta <- 0.5
  sigmazero <- 1
  sigmadelta <- -0.5
  tau <- 0.2

  # Compute densities for both responses
  density_0 <- dlnr(y, mu, mudelta, sigmazero, sigmadelta, tau, response = 0)
  density_1 <- dlnr(y, mu, mudelta, sigmazero, sigmadelta, tau, response = 1)

  # Check that densities are finite and non-negative
  expect_true(is.finite(density_0))
  expect_true(is.finite(density_1))
  expect_true(density_0 >= 0)
  expect_true(density_1 >= 0)
})

test_that("dlnr computes valid log-densities", {
  set.seed(123)
  y <- 1.5
  mu <- 1
  mudelta <- 0.5
  sigmazero <- 1
  sigmadelta <- -0.5
  tau <- 0.2

  # Compute log-densities for both responses
  log_density_0 <- dlnr(y, mu, mudelta, sigmazero, sigmadelta, tau, response = 0, log = TRUE)
  log_density_1 <- dlnr(y, mu, mudelta, sigmazero, sigmadelta, tau, response = 1, log = TRUE)

  # Check that log-densities are finite and negative
  expect_true(is.finite(log_density_0))
  expect_true(is.finite(log_density_1))
  expect_true(log_density_0 < 0)
  expect_true(log_density_1 < 0)
})

test_that("dlnr joint density integrates to 1", {
  set.seed(123)
  mu <- 1
  mudelta <- 0.5
  sigmazero <- 1
  sigmadelta <- -0.5
  tau <- 0.2

  # Define the range of reaction times for integration
  y_range <- seq(tau + 1e-6, 10, length.out = 1000)

  # Compute the joint density for response = 0
  densities_0 <- sapply(y_range, function(y) {
    dlnr(y, mu, mudelta, sigmazero, sigmadelta, tau, response = 0)
  })

  # Compute the joint density for response = 1
  densities_1 <- sapply(y_range, function(y) {
    dlnr(y, mu, mudelta, sigmazero, sigmadelta, tau, response = 1)
  })

  # Compute the total (marginal) density
  total_density <- densities_0 + densities_1

  # Integrate the total density
  integral <- sum(total_density) * diff(y_range)[1]

  # Check that the integral is close to 1
  expect_equal(integral, 1, tolerance = 0.01)
})
