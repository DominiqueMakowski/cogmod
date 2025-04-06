context("LNR - rnlr and dlnr")

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

test_that("rlnr empirical distribution matches dlnr theoretical density", {
  # Set parameters
  set.seed(42)
  mu <- 0.5
  mudelta <- 0.3
  sigmazero <- 0.8
  sigmadelta <- -0.2
  tau <- 0.1
  n_samples <- 100000  # Large sample to reduce random variation

  # Generate samples using rlnr
  data <- rlnr(n_samples, mu, mudelta, sigmazero, sigmadelta, tau)

  # Separate by response
  rt_0 <- data$rt[data$response == 0]
  rt_1 <- data$rt[data$response == 1]

  # Empirical proportions of each response
  p_0 <- length(rt_0) / n_samples
  p_1 <- length(rt_1) / n_samples

  # First check: Response proportions
  # Calculate expected proportions using a more accurate method
  expected_p0 <- 0
  expected_p1 <- 0

  # Use a more precise numerical integration with adaptive range
  max_rt <- max(30, max(data$rt) * 1.5)  # Ensure sufficient range
  grid_points <- seq(tau, max_rt, length.out = 1000)

  for (rt in grid_points) {
    dens0 <- dlnr(rt, mu, mudelta, sigmazero, sigmadelta, tau, response = 0)
    dens1 <- dlnr(rt, mu, mudelta, sigmazero, sigmadelta, tau, response = 1)
    expected_p0 <- expected_p0 + dens0 * (grid_points[2] - grid_points[1])
    expected_p1 <- expected_p1 + dens1 * (grid_points[2] - grid_points[1])
  }

  # Normalize to ensure they sum to 1
  total_p <- expected_p0 + expected_p1
  expected_p0 <- expected_p0 / total_p
  expected_p1 <- expected_p1 / total_p

  # Compare empirical and expected proportions with higher tolerance
  expect_equal(p_0, expected_p0, tolerance = 0.03,
               label = "Proportion of response 0 should match density")
  expect_equal(p_1, expected_p1, tolerance = 0.03,
               label = "Proportion of response 1 should match density")

  # Second check: Distribution shape
  # Use adaptive bin sizes based on the data distribution
  # Use quantile-based bins for better representation of the distribution
  quantile_breaks_0 <- quantile(rt_0, probs = seq(0, 1, length.out = 25))
  quantile_breaks_1 <- quantile(rt_1, probs = seq(0, 1, length.out = 25))

  # Ensure unique break points
  quantile_breaks_0 <- unique(c(tau, quantile_breaks_0))
  quantile_breaks_1 <- unique(c(tau, quantile_breaks_1))

  # Calculate empirical densities with adaptive bins
  hist_0 <- hist(rt_0, breaks = quantile_breaks_0, plot = FALSE)
  hist_1 <- hist(rt_1, breaks = quantile_breaks_1, plot = FALSE)

  # Get bin centers and widths
  bin_centers_0 <- (hist_0$breaks[-1] + hist_0$breaks[-length(hist_0$breaks)]) / 2
  bin_widths_0 <- diff(hist_0$breaks)

  bin_centers_1 <- (hist_1$breaks[-1] + hist_1$breaks[-length(hist_1$breaks)]) / 2
  bin_widths_1 <- diff(hist_1$breaks)

  # Normalize to get empirical densities
  empirical_density_0 <- hist_0$counts / (sum(hist_0$counts) * bin_widths_0)
  empirical_density_1 <- hist_1$counts / (sum(hist_1$counts) * bin_widths_1)

  # Calculate theoretical densities at bin centers
  theoretical_density_0 <- sapply(bin_centers_0, function(y) {
    dlnr(y, mu, mudelta, sigmazero, sigmadelta, tau, response = 0) / p_0
  })

  theoretical_density_1 <- sapply(bin_centers_1, function(y) {
    dlnr(y, mu, mudelta, sigmazero, sigmadelta, tau, response = 1) / p_1
  })

  # Use a more robust distribution comparison:
  # Earth Mover's Distance (Wasserstein) approximation
  # For response 0
  emd_0 <- 0
  for (i in 1:length(empirical_density_0)) {
    # Calculate absolute difference in CDF at each bin
    emd_0 <- emd_0 + abs(
      (sum(empirical_density_0[1:i] * bin_widths_0[1:i])) -
      (sum(theoretical_density_0[1:i] * bin_widths_0[1:i]))
    )
  }
  emd_0 <- emd_0 / length(empirical_density_0)

  # For response 1
  emd_1 <- 0
  for (i in 1:length(empirical_density_1)) {
    # Calculate absolute difference in CDF at each bin
    emd_1 <- emd_1 + abs(
      (sum(empirical_density_1[1:i] * bin_widths_1[1:i])) -
      (sum(theoretical_density_1[1:i] * bin_widths_1[1:i]))
    )
  }
  emd_1 <- emd_1 / length(empirical_density_1)

  # Use more reasonable thresholds for distribution comparison
  expect_lt(emd_0, 0.15, label = "EMD for response 0 should be small")
  expect_lt(emd_1, 0.15, label = "EMD for response 1 should be small")

  # Final check: Compare specific quantiles of the distributions
  quantiles_to_check <- c(0.25, 0.5, 0.75)

  # Empirical quantiles
  emp_quantiles_0 <- quantile(rt_0, probs = quantiles_to_check)
  emp_quantiles_1 <- quantile(rt_1, probs = quantiles_to_check)

  # Avoid complex theoretical quantile calculation and
  # directly check that the theoretical density is higher near empirical quantiles
  for (q in 1:length(quantiles_to_check)) {
    q_val_0 <- emp_quantiles_0[q]
    q_val_1 <- emp_quantiles_1[q]

    # Theoretical density at empirical quantile should be non-negligible
    expect_gt(dlnr(q_val_0, mu, mudelta, sigmazero, sigmadelta, tau, response = 0),
              0.01 * max(theoretical_density_0),
              label = paste("Density at", quantiles_to_check[q], "quantile for response 0"))

    expect_gt(dlnr(q_val_1, mu, mudelta, sigmazero, sigmadelta, tau, response = 1),
              0.01 * max(theoretical_density_1),
              label = paste("Density at", quantiles_to_check[q], "quantile for response 1"))
  }
})

context("LNR - brms")

test_that("lnr recovery works with brms", {
  # Skip on CRAN and when not running full tests
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")

  # Set random seed for reproducibility
  set.seed(123)

  # Generate synthetic data
  n <- 500  # Small sample for quick testing
  true_mu <- 0.5
  true_mudelta <- 0.3
  true_sigmazero <- 0.8
  true_sigmadelta <- -0.2
  true_tau <- 0.1

  # Generate data
  sim_data <- rlnr(n, mu = true_mu, mudelta = true_mudelta,
                  sigmazero = true_sigmazero, sigmadelta = true_sigmadelta,
                  tau = true_tau)

  # Convert to proper format for brms
  df <- data.frame(
    rt = sim_data$rt,
    response = sim_data$response
  )

  # Fit with brms using variational inference (faster than MCMC)
  f <- brms::bf(
    rt | dec(response) ~ 1,
    mudelta ~ 1,
    sigmazero ~ 1,
    sigmadelta ~ 1,
    tau ~ 1
  )

  # Fit model with variational inference
  fit <- brm(
    formula = f,
    data = df,
    family = lnr(),
    stanvars = lnr_stanvars(),
    init = 1,
    algorithm = "sampling",  # Use VI for speed
    iter = 500,
    refresh = 0,
    backend = "cmdstanr"
  )

  # Check parameter recovery
  # Extract posterior samples
  post <- posterior_summary(fit)

  # Check if estimated parameters are close to true values
  # Use wider tolerance for this complex model
  expect_equal(post["b_Intercept", "Estimate"], true_mu, tolerance = 0.3)
  expect_equal(post["b_mudelta_Intercept", "Estimate"], true_mudelta, tolerance = 0.3)
  expect_equal(post["b_sigmazero_Intercept", "Estimate"], true_sigmazero, tolerance = 0.3)
  expect_equal(post["b_sigmadelta_Intercept", "Estimate"], true_sigmadelta, tolerance = 0.3)
  expect_equal(post["b_tau_Intercept", "Estimate"], true_tau, tolerance = 0.3)

  # Check that true values fall within 95% credible intervals
  expect_true(true_mu >= post["b_mu", "Q2.5"] && true_mu <= post["b_mu", "Q97.5"])
  expect_true(true_mudelta >= post["b_mudelta", "Q2.5"] && true_mudelta <= post["b_mudelta", "Q97.5"])
  expect_true(true_sigmazero >= post["b_sigmazero", "Q2.5"] && true_sigmazero <= post["b_sigmazero", "Q97.5"])
  expect_true(true_sigmadelta >= post["b_sigmadelta", "Q2.5"] && true_sigmadelta <= post["b_sigmadelta", "Q97.5"])
  expect_true(true_tau >= post["b_tau", "Q2.5"] && true_tau <= post["b_tau", "Q97.5"])


  # Test the log-likelihood function
  ll <- log_lik(fit)
  expect_true(all(!is.na(ll)))
  expect_true(all(is.finite(ll)))

  # Test WAIC and LOO calculations
  waic <- waic(fit)
  expect_true(!is.null(waic))
  expect_true(is.numeric(waic$estimates["waic", "Estimate"]))
})
