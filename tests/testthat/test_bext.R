context("BEXT - rbext")

test_that("rbext generates valid values", {
  # Basic functionality
  set.seed(123)
  x <- rbext(1000, mu = 0.5, phi = 3, pex = 0.1, bex = 0.5)
  expect_true(all(x >= 0 & x <= 1))

  # Correct number of elements
  expect_equal(length(rbext(50, mu = 0.7, phi = 5, pex = 0.2, bex = 0.3)), 50)

  # Extreme values as expected
  set.seed(123)
  x_high_pex <- rbext(10000, mu = 0.5, phi = 3, pex = 0.8, bex = 0.5)
  zeros <- sum(x_high_pex == 0)
  ones <- sum(x_high_pex == 1)
  expect_gt(zeros + ones, 7000) # At least 70% extreme values
  expect_gt(zeros, 3000) # About 40% zeros (pex * (1-bex))
  expect_gt(ones, 3000)  # About 40% ones (pex * bex)

  # All zeros case
  all_zeros <- rbext(100, mu = 0.5, phi = 3, pex = 1, bex = 0)
  expect_equal(sum(all_zeros == 0), 100)

  # All ones case
  all_ones <- rbext(100, mu = 0.5, phi = 3, pex = 1, bex = 1)
  expect_equal(sum(all_ones == 1), 100)
})

test_that("rbext validates inputs correctly", {
  expect_error(rbext(10, mu = -0.1, phi = 3, pex = 0.1, bex = 0.5), "mu must be between 0 and 1")
  expect_error(rbext(10, mu = 1.1, phi = 3, pex = 0.1, bex = 0.5), "mu must be between 0 and 1")
  expect_error(rbext(10, mu = 0.5, phi = -1, pex = 0.1, bex = 0.5), "phi must be positive")
  expect_error(rbext(10, mu = 0.5, phi = 3, pex = -0.1, bex = 0.5), "pex must be between 0 and 1")
  expect_error(rbext(10, mu = 0.5, phi = 3, pex = 1.1, bex = 0.5), "pex must be between 0 and 1")
  expect_error(rbext(10, mu = 0.5, phi = 3, pex = 0.1, bex = -0.1), "bex must be between 0 and 1")
  expect_error(rbext(10, mu = 0.5, phi = 3, pex = 0.1, bex = 1.1), "bex must be between 0 and 1")
})

context("BEXT - dbext")

test_that("dbext calculates correct densities", {
  # Basic functionality
  x <- seq(0, 1, by = 0.1)
  d <- dbext(x, mu = 0.5, phi = 3, pex = 0.2, bex = 0.5)
  expect_equal(length(d), length(x))
  expect_true(all(d >= 0))

  # Point masses at 0 and 1
  expect_equal(dbext(0, mu = 0.5, phi = 3, pex = 0.2, bex = 0.5), 0.2 * (1 - 0.5)) # kleft
  expect_equal(dbext(1, mu = 0.5, phi = 3, pex = 0.2, bex = 0.5), 0.2 * 0.5)      # 1 - kright

  # Values outside [0,1] should have density 0
  expect_equal(dbext(-0.1, mu = 0.5, phi = 3, pex = 0.2, bex = 0.5), 0)
  expect_equal(dbext(1.1, mu = 0.5, phi = 3, pex = 0.2, bex = 0.5), 0)

  # Log densities
  d_log <- dbext(x, mu = 0.5, phi = 3, pex = 0.2, bex = 0.5, log = TRUE)
  expect_equal(d_log[d > 0], log(d[d > 0]))
  expect_true(all(d_log[d == 0] == -Inf))

  # Vectorization for mu and phi
  x_vec <- c(0.2, 0.4, 0.6, 0.8)
  mu_vec <- c(0.3, 0.4, 0.5, 0.6)
  phi_vec <- c(2, 3, 4, 5)

  # Test with vectorized mu
  d_mu_vec <- dbext(x_vec, mu = mu_vec, phi = 3, pex = 0.2, bex = 0.5)
  expect_equal(length(d_mu_vec), length(x_vec))

  # Test with vectorized phi
  d_phi_vec <- dbext(x_vec, mu = 0.5, phi = phi_vec, pex = 0.2, bex = 0.5)
  expect_equal(length(d_phi_vec), length(x_vec))
})

test_that("dbext validates inputs correctly", {
  expect_error(dbext(0.5, mu = -0.1, phi = 3, pex = 0.1, bex = 0.5), "mu must be between 0 and 1")
  expect_error(dbext(0.5, mu = 1.1, phi = 3, pex = 0.1, bex = 0.5), "mu must be between 0 and 1")
  expect_error(dbext(0.5, mu = 0.5, phi = -1, pex = 0.1, bex = 0.5), "phi must be positive")
  expect_error(dbext(0.5, mu = 0.5, phi = 3, pex = -0.1, bex = 0.5), "pex must be between 0 and 1")
  expect_error(dbext(0.5, mu = 0.5, phi = 3, pex = 1.1, bex = 0.5), "pex must be between 0 and 1")
  expect_error(dbext(0.5, mu = 0.5, phi = 3, pex = 0.1, bex = -0.1), "bex must be between 0 and 1")
  expect_error(dbext(0.5, mu = 0.5, phi = 3, pex = 0.1, bex = 1.1), "bex must be between 0 and 1")

  # Test invalid parameterization (kleft >= kright) for non-special case
  # Find pex and bex values that produce kleft >= kright
  pex_test <- 0.7
  bex_test <- 0.6  # This gives kleft = 0.28 and kright = 0.58, so kleft < kright (valid)

  # Adjust bex to create an invalid case: if bex = 0.9, then:
  # kleft = 0.7 * (1 - 0.9) = 0.07
  # kright = 1 - (0.7 * 0.9) = 0.37
  # So kleft < kright (still valid)

  # Note: pex = 1 is now a valid special case (Bernoulli distribution)
  # So this should not throw an error:
  expect_silent(dbext(0.5, mu = 0.5, phi = 3, pex = 1, bex = 0.5))
})

test_that("Density function integrates to approximately 1", {
  # Create a fine grid for numerical integration
  x_grid <- seq(0, 1, length.out = 10001)

  # Test different parameter configurations
  params <- list(
    list(mu = 0.5, phi = 3, pex = 0.2, bex = 0.5),
    list(mu = 0.7, phi = 10, pex = 0.1, bex = 0.3),
    list(mu = 0.3, phi = 5, pex = 0.5, bex = 0.8)
  )

  for (p in params) {
    d <- dbext(x_grid, mu = p$mu, phi = p$phi, pex = p$pex, bex = p$bex)

    # Adjust for point masses at 0 and 1 in numerical integration
    # The first and last points include the point masses
    integral <- sum(d[-c(1, length(d))]) * (x_grid[2] - x_grid[1])
    integral <- integral + d[1] + d[length(d)]

    expect_equal(integral, 1, tolerance = 0.01)
  }
})

test_that("pex = 1 special case works correctly", {
  # Test that rbext with pex = 1 generates only 0s and 1s
  n <- 1000
  set.seed(42)
  x1 <- rbext(n, mu = 0.5, phi = 3, pex = 1, bex = 0.3)

  # Should only contain 0s and 1s
  expect_equal(sort(unique(x1)), c(0, 1))

  # Check correct proportions: with bex = 0.3, should be 30% 1s and 70% 0s
  prop_ones <- sum(x1 == 1) / n
  expect_equal(prop_ones, 0.3, tolerance = 0.05)

  # Different bex values should result in different proportions
  set.seed(42)
  x2 <- rbext(n, mu = 0.5, phi = 3, pex = 1, bex = 0.7)
  prop_ones_2 <- sum(x2 == 1) / n
  expect_equal(prop_ones_2, 0.7, tolerance = 0.05)

  # Density function should have correct mass at endpoints
  expect_equal(dbext(0, mu = 0.5, phi = 3, pex = 1, bex = 0.3), 0.7)
  expect_equal(dbext(1, mu = 0.5, phi = 3, pex = 1, bex = 0.3), 0.3)

  # All values between 0 and 1 (exclusive) should have zero density
  x_mid <- seq(0.1, 0.9, by = 0.1)
  d_mid <- dbext(x_mid, mu = 0.5, phi = 3, pex = 1, bex = 0.3)
  expect_equal(sum(d_mid), 0)

  # Log density should work correctly
  d_log <- dbext(c(0, 0.5, 1), mu = 0.5, phi = 3, pex = 1, bex = 0.3, log = TRUE)
  expect_equal(d_log[1], log(0.7))
  expect_equal(d_log[2], -Inf)
  expect_equal(d_log[3], log(0.3))
})

test_that("dbext integrates to 1 for all parameterizations", {
  # Define a set of parameter combinations to test
  test_params <- list(
    list(mu = 0.5, phi = 3, pex = 0.1, bex = 0.5),  # Low extremes
    list(mu = 0.7, phi = 5, pex = 0.4, bex = 0.3),  # Moderate extremes, unbalanced
    list(mu = 0.3, phi = 8, pex = 0.8, bex = 0.6),  # High extremes, unbalanced
    list(mu = 0.5, phi = 3, pex = 1.0, bex = 0.5)   # Special case: all extremes
  )

  for (params in test_params) {
    # Create a fine grid for numerical integration of the continuous part
    x_cont <- seq(0.00001, 0.99999, length.out = 10000)

    # Calculate densities
    d_cont <- dbext(x_cont, mu = params$mu, phi = params$phi,
                   pex = params$pex, bex = params$bex)

    # Get masses at endpoints
    d_0 <- dbext(0, mu = params$mu, phi = params$phi,
                pex = params$pex, bex = params$bex)
    d_1 <- dbext(1, mu = params$mu, phi = params$phi,
                pex = params$pex, bex = params$bex)

    # Integrate the continuous part using trapezoidal rule
    int_cont <- sum(d_cont) * (x_cont[2] - x_cont[1])

    # Total probability should sum to approximately 1
    total_prob <- d_0 + int_cont + d_1

    expect_equal(total_prob, 1, tolerance = 0.01,
                 info = paste("Failed with params:",
                             "mu =", params$mu,
                             "phi =", params$phi,
                             "pex =", params$pex,
                             "bex =", params$bex))
  }
})

test_that("rbext and dbext are consistent with each other", {
  # Test with moderate extremes
  params <- list(mu = 0.6, phi = 4, pex = 0.3, bex = 0.4)
  set.seed(123)
  samples <- rbext(50000, mu = params$mu, phi = params$phi,
                  pex = params$pex, bex = params$bex)

  # Count zeros and ones
  zeros <- sum(samples == 0) / length(samples)
  ones <- sum(samples == 1) / length(samples)

  # Expected values from density function
  exp_zeros <- dbext(0, mu = params$mu, phi = params$phi,
                    pex = params$pex, bex = params$bex)
  exp_ones <- dbext(1, mu = params$mu, phi = params$phi,
                   pex = params$pex, bex = params$bex)

  # Check proportion of extremes
  expect_equal(zeros, exp_zeros, tolerance = 0.01)
  expect_equal(ones, exp_ones, tolerance = 0.01)

  # Test with all extremes (pex = 1)
  params2 <- list(mu = 0.5, phi = 3, pex = 1, bex = 0.25)
  set.seed(456)
  samples2 <- rbext(10000, mu = params2$mu, phi = params2$phi,
                   pex = params2$pex, bex = params2$bex)

  zeros2 <- sum(samples2 == 0) / length(samples2)
  ones2 <- sum(samples2 == 1) / length(samples2)

  # For pex = 1, zeros should be 1-bex and ones should be bex
  expect_equal(zeros2, 1 - params2$bex, tolerance = 0.01)
  expect_equal(ones2, params2$bex, tolerance = 0.01)
})

test_that("Edge cases are handled correctly", {
  # pex = 0 (pure beta, no extremes)
  x_values <- seq(0.1, 0.9, by = 0.1)
  pure_beta_dens <- dbext(x_values, mu = 0.5, phi = 3, pex = 0, bex = 0.5)
  beta_dens <- stats::dbeta(x_values, shape1 = 0.5 * 3, shape2 = 0.5 * 3)
  expect_equal(pure_beta_dens, beta_dens, tolerance = 1e-5)

  # No mass should be at extremes when pex = 0
  expect_equal(dbext(0, mu = 0.5, phi = 3, pex = 0, bex = 0.5), 0)
  expect_equal(dbext(1, mu = 0.5, phi = 3, pex = 0, bex = 0.5), 0)

  # bex = 0 (only zeros, no ones when pex = 1)
  all_zeros <- rbext(100, mu = 0.5, phi = 3, pex = 1, bex = 0)
  expect_equal(all_zeros, rep(0, 100))
  expect_equal(dbext(0, mu = 0.5, phi = 3, pex = 1, bex = 0), 1)
  expect_equal(dbext(1, mu = 0.5, phi = 3, pex = 1, bex = 0), 0)

  # bex = 1 (only ones, no zeros when pex = 1)
  all_ones <- rbext(100, mu = 0.5, phi = 3, pex = 1, bex = 1)
  expect_equal(all_ones, rep(1, 100))
  expect_equal(dbext(0, mu = 0.5, phi = 3, pex = 1, bex = 1), 0)
  expect_equal(dbext(1, mu = 0.5, phi = 3, pex = 1, bex = 1), 1)
})

test_that("rbext produces distribution consistent with dbext", {
  set.seed(123)
  # Generate random samples
  params <- list(mu = 0.6, phi = 4, pex = 0.3, bex = 0.4)
  samples <- rbext(100000, mu = params$mu, phi = params$phi,
                  pex = params$pex, bex = params$bex)

  # Count zeros and ones (point masses)
  zeros <- sum(samples == 0) / length(samples)
  ones <- sum(samples == 1) / length(samples)

  # Calculate expected values from dbext
  expected_zeros <- dbext(0, mu = params$mu, phi = params$phi,
                         pex = params$pex, bex = params$bex)
  expected_ones <- dbext(1, mu = params$mu, phi = params$phi,
                        pex = params$pex, bex = params$bex)

  # Test point masses match
  expect_equal(zeros, expected_zeros, tolerance = 0.01)
  expect_equal(ones, expected_ones, tolerance = 0.01)

  # Test the continuous part using histogram
  cont_samples <- samples[samples > 0 & samples < 1]
  if (length(cont_samples) > 1000) {
    # Create bins and count frequencies
    bins <- seq(0, 1, by = 0.05)
    hist_counts <- hist(cont_samples, breaks = bins, plot = FALSE)$counts
    hist_probs <- hist_counts / length(samples)

    # Calculate expected probabilities in each bin
    bin_centers <- (bins[-1] + bins[-length(bins)]) / 2
    expected_densities <- dbext(bin_centers, mu = params$mu, phi = params$phi,
                               pex = params$pex, bex = params$bex)
    expected_probs <- expected_densities * 0.05  # bin width

    # Test histogram shape matches density function (correlate rather than exact match)
    corr <- cor(hist_probs, expected_probs)
    expect_gt(corr, 0.9)  # Strong correlation between empirical and theoretical distributions
  }
})


test_that("BEXT functions match ordbetareg reference implementation", {
  skip_on_cran()

  # Source the ordbetareg distribution functions from GitHub
  source("https://raw.githubusercontent.com/saudiwin/ordbetareg_pack/master/R/distribution.R")

# Test several parameter configurations
  test_cases <- list(
    list(mu = 0.5, phi = 3, pex = 0.2, bex = 0.5),
    list(mu = 0.7, phi = 8, pex = 0.3, bex = 0.4),
    list(mu = 0.3, phi = 5, pex = 0.1, bex = 0.7)
  )

  for (params in test_cases) {
    # Calculate the equivalent kleft and kright from BEXT parameterization
    kleft <- params$pex * (1 - params$bex)  # probability mass at 0
    kright <- 1 - (params$pex * params$bex) # threshold above which outcomes are set to 1

    # Calculate equivalent cutpoints for ordbeta functions
    # Based on dordbeta implementation:
    # low = 1 - plogis(mu_ql - cutpoints[1]) should equal kleft
    # high = plogis(mu_ql - cutpoints[2]) should equal 1-kright
    mu_ql <- qlogis(params$mu)

    # Solving for cutpoints:
    # 1 - plogis(mu_ql - cutpoints[1]) = kleft
    # plogis(mu_ql - cutpoints[1]) = 1 - kleft
    # mu_ql - cutpoints[1] = qlogis(1 - kleft)
    # cutpoints[1] = mu_ql - qlogis(1 - kleft)
    cutpoint1 <- mu_ql - qlogis(1 - kleft)

    # plogis(mu_ql - cutpoints[2]) = 1-kright
    # mu_ql - cutpoints[2] = qlogis(1-kright)
    # cutpoints[2] = mu_ql - qlogis(1-kright)
    cutpoint2 <- mu_ql - qlogis(1-kright)

    # Test that our conversion is correct by checking the resulting probabilities
    low <- 1 - plogis(mu_ql - cutpoint1)
    high <- plogis(mu_ql - cutpoint2)
    middle <- plogis(mu_ql - cutpoint1) - plogis(mu_ql - cutpoint2)

    # Verify probabilities match expected values
    expect_equal(low, kleft, tolerance = 1e-5)
    expect_equal(high, 1-kright, tolerance = 1e-5)
    expect_equal(middle, kright-kleft, tolerance = 1e-5)

    # Test density function equivalence
    x_values <- seq(0.1, 0.9, by = 0.1)  # Avoid 0 and 1 for continuous comparison

    bext_dens <- dbext(x_values, mu = params$mu, phi = params$phi,
                        pex = params$pex, bex = params$bex)

    ordbeta_dens <- dordbeta(x_values, mu = params$mu, phi = params$phi,
                              cutpoints = c(cutpoint1, cutpoint2))

    expect_equal(bext_dens, ordbeta_dens, tolerance = 1e-5)

    # Test point masses at 0 and 1
    expect_equal(dbext(0, mu = params$mu, phi = params$phi,
                        pex = params$pex, bex = params$bex),
                  dordbeta(0, mu = params$mu, phi = params$phi,
                          cutpoints = c(cutpoint1, cutpoint2)))

    expect_equal(dbext(1, mu = params$mu, phi = params$phi,
                        pex = params$pex, bex = params$bex),
                  dordbeta(1, mu = params$mu, phi = params$phi,
                          cutpoints = c(cutpoint1, cutpoint2)))

    # Test random generation
    set.seed(42)
    bext_sample <- rbext(1000, mu = params$mu, phi = params$phi,
                          pex = params$pex, bex = params$bex)

    set.seed(42)
    ordbeta_sample <- rordbeta(1000, mu = params$mu, phi = params$phi,
                                cutpoints = c(cutpoint1, cutpoint2))

    # Compare distributions using statistical tests
    expect_equal(mean(bext_sample), mean(ordbeta_sample), tolerance = 0.05)
    expect_equal(sd(bext_sample), sd(ordbeta_sample), tolerance = 0.05)

    # Compare proportions of exact zeros and ones
    expect_equal(sum(bext_sample == 0) / 1000, sum(ordbeta_sample == 0) / 1000,
                  tolerance = 0.05)
    expect_equal(sum(bext_sample == 1) / 1000, sum(ordbeta_sample == 1) / 1000,
                  tolerance = 0.05)
  }
})

context("BEXT - brms")

test_that("BEXT model can recover parameters with brms using variational inference", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")

  # Set parameters for data generation
  true_mu <- 0.6
  true_phi <- 4
  true_pex <- 0.15
  true_bex <- 0.4

  # Generate synthetic data with known parameters
  set.seed(123)
  df <- data.frame(
    y = rbext(n = 1000, mu = true_mu, phi = true_phi,
              pex = true_pex, bex = true_bex)
  )

  # Formula for intercept-only model
  f <- brms::bf(
    y ~ 1,
    phi ~ 1,
    pex ~ 1,
    bex ~ 1
  )

  # Fit model with variational inference
  m <- brms::brm(f,
    data = df, family = bext(), stanvars = bext_stanvars(),
    iter = 10000, # More iterations for VI
    algorithm = "pathfinder", # Variational inference
    backend = "cmdstanr", seed = 123, refresh = 0
  )

  # Extract posterior means
  post <- brms::posterior_summary(m)
  est_mu <- brms::inv_logit_scaled(post["b_Intercept", "Estimate"])
  est_phi <- log(1 + exp(post["b_phi_Intercept", "Estimate"]))  # softplus(x) = log(1 + exp(x))
  est_pex <- brms::inv_logit_scaled(post["b_pex_Intercept", "Estimate"])
  est_bex <- brms::inv_logit_scaled(post["b_bex_Intercept", "Estimate"])

  # Check that estimates are reasonably close to true values
  # Using wider tolerance for VI since it's an approximation
  expect_equal(est_mu, true_mu, tolerance = 0.15)
  expect_equal(est_phi, true_phi, tolerance = 0.3)
  expect_equal(est_pex, true_pex, tolerance = 0.15)
  expect_equal(est_bex, true_bex, tolerance = 0.3)

  # Test posterior prediction
  pred <- brms::posterior_predict(m, ndraws = 10)
  expect_equal(nrow(pred), 10)
  expect_equal(ncol(pred), nrow(df))
  expect_true(all(pred >= 0 & pred <= 1))

  # Test with newdata
  pred <- brms::posterior_predict(m, ndraws = 10, newdata = df[1:5, ])
  expect_equal(nrow(pred), 10)
  expect_equal(ncol(pred), 5)

  # Also test log-likelihood calculation
  ll <- brms::log_lik(m, ndraws = 5)
  expect_equal(nrow(ll), 5)
  expect_equal(ncol(ll), nrow(df))
})


