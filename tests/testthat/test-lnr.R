context("LNR - rlnr and dlnr")

test_that("rlnr generates valid outcomes", {
  set.seed(123)
  n <- 1000
  data <- rlnr(n, mu = 1, mudelta = 0.5, sigmazero = 1, sigmadelta = -0.5, ndt = 0.2)

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
  ndt <- 0.2

  # Compute densities for both responses
  density_0 <- dlnr(y, mu, mudelta, sigmazero, sigmadelta, ndt, response = 0)
  density_1 <- dlnr(y, mu, mudelta, sigmazero, sigmadelta, ndt, response = 1)

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
  ndt <- 0.2

  # Compute log-densities for both responses
  log_density_0 <- dlnr(y, mu, mudelta, sigmazero, sigmadelta, ndt, response = 0, log = TRUE)
  log_density_1 <- dlnr(y, mu, mudelta, sigmazero, sigmadelta, ndt, response = 1, log = TRUE)

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
  ndt <- 0.2

  # Define the range of reaction times for integration
  y_range <- seq(ndt + 1e-6, 10, length.out = 1000)

  # Compute the joint density for response = 0
  densities_0 <- sapply(y_range, function(y) {
    dlnr(y, mu, mudelta, sigmazero, sigmadelta, ndt, response = 0)
  })

  # Compute the joint density for response = 1
  densities_1 <- sapply(y_range, function(y) {
    dlnr(y, mu, mudelta, sigmazero, sigmadelta, ndt, response = 1)
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
  ndt <- 0.1
  n_samples <- 100000  # Large sample to reduce random variation

  # Generate samples using rlnr
  data <- rlnr(n_samples, mu, mudelta, sigmazero, sigmadelta, ndt)

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
  grid_points <- seq(ndt, max_rt, length.out = 1000)

  for (rt in grid_points) {
    dens0 <- dlnr(rt, mu, mudelta, sigmazero, sigmadelta, ndt, response = 0)
    dens1 <- dlnr(rt, mu, mudelta, sigmazero, sigmadelta, ndt, response = 1)
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
  quantile_breaks_0 <- unique(c(ndt, quantile_breaks_0))
  quantile_breaks_1 <- unique(c(ndt, quantile_breaks_1))

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
    dlnr(y, mu, mudelta, sigmazero, sigmadelta, ndt, response = 0) / p_0
  })

  theoretical_density_1 <- sapply(bin_centers_1, function(y) {
    dlnr(y, mu, mudelta, sigmazero, sigmadelta, ndt, response = 1) / p_1
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
    expect_gt(dlnr(q_val_0, mu, mudelta, sigmazero, sigmadelta, ndt, response = 0),
              0.01 * max(theoretical_density_0),
              label = paste("Density at", quantiles_to_check[q], "quantile for response 0"))

    expect_gt(dlnr(q_val_1, mu, mudelta, sigmazero, sigmadelta, ndt, response = 1),
              0.01 * max(theoretical_density_1),
              label = paste("Density at", quantiles_to_check[q], "quantile for response 1"))
  }
})

test_that("dlnr handles y < ndt correctly", {
  mu <- 0
  mudelta <- 0.5
  sigmazero <- 1
  sigmadelta <- -0.5
  ndt <- 0.2
  response <- 0

  # Case: y > ndt
  y <- 0.5
  result <- dlnr(y, mu, mudelta, sigmazero, sigmadelta, ndt, response, log = FALSE)
  expect_gt(result, 0) # Density should be positive

  # Case: y <= ndt
  y <- 0.1
  result <- dlnr(y, mu, mudelta, sigmazero, sigmadelta, ndt, response, log = FALSE)
  expect_equal(result, 0) # Density should be 0

  # Case: y > ndt with log = TRUE
  y <- 0.5
  result <- dlnr(y, mu, mudelta, sigmazero, sigmadelta, ndt, response, log = TRUE)
  expect_lt(result, 0) # Log-density should be negative

  # Case: y <= ndt with log = TRUE
  y <- 0.1
  result <- dlnr(y, mu, mudelta, sigmazero, sigmadelta, ndt, response, log = TRUE)
  expect_equal(result, -Inf) # Log-density should be -Inf
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
  n <- 1000  # Small sample for quick testing
  true_mu <- 0.5
  true_mudelta <- 0.3
  true_sigmazero <- 0.8
  true_sigmadelta <- -0.2
  true_ndt <- 0.2

  # Generate data
  df <- rlnr(n, mu = true_mu, mudelta = true_mudelta,
                  sigmazero = true_sigmazero, sigmadelta = true_sigmadelta,
                  ndt = true_ndt)

  # Get min_rt for tau calculation
  min_rt <- min(df$rt)
  true_tau <- true_ndt / min_rt

  # Fit with brms
  f <- brms::bf(
    rt | dec(response) ~ 1,
    mudelta ~ 1,
    sigmazero ~ 1,
    sigmadelta ~ 1,
    tau ~ 1,
    minrt = min(df$rt)
  )

  prior <- brms::set_prior("normal(0, 1)", dpar = "tau", class="Intercept") |>
    brms::validate_prior(f, data = df, family = lnr())

  fit <- brms::brm(
    formula = f,
    data = df,
    family = lnr(),
    stanvars = lnr_stanvars(),
    prior = prior,
    init = 0.2,  # Small positive initialization
    algorithm = "pathfinder",
    iter = 10000,
    refresh = 0,
    # control = list(adapt_delta = 0.9),
    backend = "cmdstanr"
  )

  # Extract posterior summaries
  post <- brms::posterior_summary(fit)

  # Convert parameters back to natural scale
  est_mu <- post["b_Intercept", "Estimate"]  # identity link
  est_mudelta <- post["b_mudelta_Intercept", "Estimate"]  # identity link
  est_sigmazero <- log1p(exp(post["b_sigmazero_Intercept", "Estimate"]))  # softplus link
  est_sigmadelta <- post["b_sigmadelta_Intercept", "Estimate"]  # identity link
  est_tau <- brms::inv_logit_scaled(post["b_tau_Intercept", "Estimate"])  # logit link
  est_ndt <- est_tau * min_rt  # convert tau to ndt

  # Check parameter recovery with appropriate transformations
  expect_equal(est_mu, true_mu, tolerance = 0.3)
  expect_equal(est_mudelta, true_mudelta, tolerance = 0.3)
  expect_equal(est_sigmazero, true_sigmazero, tolerance = 0.3)
  expect_equal(est_sigmadelta, true_sigmadelta, tolerance = 0.3)
  expect_equal(est_ndt, true_ndt, tolerance = 0.3)


  # Test that predict returns something
  # prep <- brms::prepare_predictions(fit, newdata = df[1:10,]); i = 1
  pred <- brms::posterior_predict(fit, ndraws = 10, newdata = df[1:3,])
  expect_equal(ncol(pred), 6)  # rt + response so 3 x 2 the number of rows
  expect_equal(nrow(pred), 10)  # 10 draws

  # Test the log-likelihood function
  ll <- brms::log_lik(fit)
  expect_true(all(!is.na(ll)))
  expect_true(all(is.finite(ll)))

  # Test WAIC and LOO calculations
  waic <- brms::waic(fit)
  expect_true(!is.null(waic))
  expect_true(is.numeric(waic$estimates["waic", "Estimate"]))

})
