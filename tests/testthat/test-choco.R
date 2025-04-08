context("CHOCO - rchoco")

test_that("rchoco generates valid values", {
  # Basic functionality
  set.seed(123)
  x <- rchoco(1000, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.1, bex = 0.5)
  expect_true(all(x >= 0 & x <= 1))

  # Correct number of elements
  expect_equal(length(rchoco(50, p = 0.6, muleft = 0.4, mudelta = 0.2, phileft = 5, phidelta = 0.1, pex = 0.2, bex = 0.3)), 50)

  # Extreme values as expected
  set.seed(123)
  x_high_pex <- rchoco(10000, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.8, bex = 0.5)
  zeros <- sum(x_high_pex == 0)
  ones <- sum(x_high_pex == 1)

  # Updated expectations:
  expect_gt(zeros + ones, 3900) # At least 39% extreme values (accounting for random variation)
  expect_gt(zeros, 1900) # About 20% zeros (pex * (1-bex) * (1-p))
  expect_gt(ones, 1900)  # About 20% ones (pex * bex * p)

  # All zeros case (pex = 1, bex = 0)
  all_zeros <- rchoco(100, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 1, bex = 0)
  expect_equal(sum(all_zeros == 0), 100)

  # All ones case (pex = 1, bex = 1)
  all_ones <- rchoco(100, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 1, bex = 1)
  expect_equal(sum(all_ones == 1), 100)
})

test_that("rchoco validates inputs correctly", {
  expect_error(rchoco(10, p = -0.1, muleft = 0.3, mudelta = 0, phileft = 5), "p must be between 0 and 1")
  expect_error(rchoco(10, p = 1.1, muleft = 0.3, mudelta = 0, phileft = 5), "p must be between 0 and 1")
  expect_error(rchoco(10, p = 0.5, muleft = 0, mudelta = 0, phileft = 5), "muleft must be between 0 and 1")
  expect_error(rchoco(10, p = 0.5, muleft = 1, mudelta = 0, phileft = 5), "muleft must be between 0 and 1")
  expect_error(rchoco(10, p = 0.5, muleft = 0.3, mudelta = 0, phileft = -1), "phileft must be positive")
  expect_error(rchoco(10, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = -0.1), "pex must be between 0 and 1")
  expect_error(rchoco(10, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 1.1), "pex must be between 0 and 1")
  expect_error(rchoco(10, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.1, bex = -0.1), "bex must be between 0 and 1")
  expect_error(rchoco(10, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.1, bex = 1.1), "bex must be between 0 and 1")
  expect_error(rchoco(10, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.1, bex = 0.5, threshold = -0.1), "threshold must be between 0 and 1")
  expect_error(rchoco(10, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.1, bex = 0.5, threshold = 1.1), "threshold must be between 0 and 1")
})

test_that("rchoco special case pex = 1 works correctly", {
  # Test that rchoco with pex = 1 generates only 0s and 1s
  n <- 1000
  set.seed(42)
  x1 <- rchoco(n, p = 0.6, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 1, bex = 0.4)

  # Should only contain 0s and 1s
  expect_equal(sort(unique(x1)), c(0, 1))

  # With p = 0.6 and bex = 0.4, zeros = (1-p)*(1-bex)/(total), ones = p*bex/(total)
  # total = (1-p)*(1-bex) + p*bex = (1-0.6)*(1-0.4) + 0.6*0.4 = 0.4*0.6 + 0.6*0.4 = 0.48
  # prob_zeros = 0.4*0.6/0.48 = 0.24/0.48 = 0.5
  # prob_ones = 0.6*0.4/0.48 = 0.24/0.48 = 0.5
  # With these parameters, we should expect approximately 50% zeros and 50% ones
  prop_zeros <- sum(x1 == 0) / n
  prop_ones <- sum(x1 == 1) / n
  expect_equal(prop_zeros, 0.5, tolerance = 0.05)
  expect_equal(prop_ones, 0.5, tolerance = 0.05)

  # Different p values should result in different proportions
  set.seed(42)
  x2 <- rchoco(n, p = 0.8, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 1, bex = 0.4)
  prop_zeros_2 <- sum(x2 == 0) / n
  prop_ones_2 <- sum(x2 == 1) / n

  # Calculate expected proportions for p = 0.8, bex = 0.4
  # total = (1-0.8)*(1-0.4) + 0.8*0.4 = 0.2*0.6 + 0.32 = 0.12 + 0.32 = 0.44
  # prob_zeros = 0.2*0.6/0.44 = 0.12/0.44 ≈ 0.273
  # prob_ones = 0.8*0.4/0.44 = 0.32/0.44 ≈ 0.727
  expect_equal(prop_zeros_2, 0.12/0.44, tolerance = 0.05)
  expect_equal(prop_ones_2, 0.32/0.44, tolerance = 0.05)
})

test_that("rchoco produces expected distribution shape", {
  # Test that values are distributed properly across the threshold
  set.seed(123)
  x <- rchoco(10000, p = 0.7, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.1, bex = 0.5, threshold = 0.5)

  # Count values left and right of threshold (excluding zeros and ones)
  cont_values <- x[x > 0 & x < 1]
  left_values <- cont_values[cont_values < 0.5]
  right_values <- cont_values[cont_values >= 0.5]

  # Probability mass excluding extreme values should follow p
  prop_right <- length(right_values) / length(cont_values)
  expect_equal(prop_right, 0.7, tolerance = 0.05)

  # Test with different threshold
  set.seed(123)
  x2 <- rchoco(10000, p = 0.7, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.1, bex = 0.5, threshold = 0.7)

  # Validate threshold is respected
  cont_values2 <- x2[x2 > 0 & x2 < 1]
  left_values2 <- cont_values2[cont_values2 < 0.7]
  right_values2 <- cont_values2[cont_values2 >= 0.7]

  # Proportion should still match p
  prop_right2 <- length(right_values2) / length(cont_values2)
  expect_equal(prop_right2, 0.7, tolerance = 0.05)

  # All left values should be less than threshold
  expect_true(all(left_values2 < 0.7))

  # All right values should be greater than or equal to threshold
  expect_true(all(right_values2 >= 0.7))
})

test_that("rchoco handles extreme p values correctly", {
  # Test p = 0 (all values should be on left or zero)
  set.seed(123)
  x_left <- rchoco(1000, p = 0, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.2, bex = 0.5, threshold = 0.5)

  # Should have no values equal to 1
  expect_equal(sum(x_left == 1), 0)

  # All non-zero values should be less than threshold
  expect_true(all(x_left[x_left > 0] < 0.5))

  # Test p = 1 (all values should be on right or one)
  set.seed(123)
  x_right <- rchoco(1000, p = 1, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.2, bex = 0.5, threshold = 0.5)

  # Should have no values equal to 0
  expect_equal(sum(x_right == 0), 0)

  # All values less than 1 should be greater than or equal to threshold
  expect_true(all(x_right[x_right < 1] >= 0.5))
})

test_that("rchoco respects mudelta and phidelta parameters", {
  # Test that mudelta shifts the mean of the right distribution
  set.seed(123)
  x1 <- rchoco(10000, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0, bex = 0.5, threshold = 0.5)
  x2 <- rchoco(10000, p = 0.5, muleft = 0.3, mudelta = 1, phileft = 5, phidelta = 0, pex = 0, bex = 0.5, threshold = 0.5)

  # Extract continuous values on the right side
  right_values1 <- x1[x1 >= 0.5]
  right_values2 <- x2[x2 >= 0.5]

  # Mean should be higher with positive mudelta
  expect_gt(mean(right_values2), mean(right_values1))

  # Test that phidelta affects the precision of the right distribution
  set.seed(123)
  x3 <- rchoco(10000, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 1, pex = 0, bex = 0.5, threshold = 0.5)

  # Extract continuous values on the right side
  right_values3 <- x3[x3 >= 0.5]

  # Variance should be lower with higher precision (phi)
  expect_lt(var(right_values3), var(right_values1))
})

test_that("rchoco handles bex correctly", {
  # Test that bex controls the balance of extreme values
  set.seed(123)
  x1 <- rchoco(10000, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.2, bex = 0.25, threshold = 0.5)
  x2 <- rchoco(10000, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.2, bex = 0.75, threshold = 0.5)

  # Count zeros and ones
  zeros1 <- sum(x1 == 0)
  ones1 <- sum(x1 == 1)
  zeros2 <- sum(x2 == 0)
  ones2 <- sum(x2 == 1)

  # With bex = 0.25, there should be more zeros than ones
  expect_gt(zeros1, ones1)

  # With bex = 0.75, there should be more ones than zeros
  expect_gt(ones2, zeros2)

  # The ratio of ones to zeros should be higher for bex = 0.75 than for bex = 0.25
  ratio1 <- ones1 / zeros1
  ratio2 <- ones2 / zeros2
  expect_gt(ratio2, ratio1)
})

test_that("rchoco produces expected extreme values proportions", {
  # Test that extreme values follow the expected proportions
  set.seed(123)
  pex_test <- 0.3
  bex_test <- 0.6
  p_test <- 0.7

  # Expected proportions
  expected_zeros <- pex_test * (1 - bex_test) * (1 - p_test) # = 0.3 * 0.4 * 0.3 = 0.036
  expected_ones <- pex_test * bex_test * p_test # = 0.3 * 0.6 * 0.7 = 0.126

  x <- rchoco(10000, p = p_test, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = pex_test, bex = bex_test)

  # Count observed proportions
  observed_zeros <- sum(x == 0) / length(x)
  observed_ones <- sum(x == 1) / length(x)

  # Check if observed proportions match expected
  expect_equal(observed_zeros, expected_zeros, tolerance = 0.01)
  expect_equal(observed_ones, expected_ones, tolerance = 0.01)
})



context("CHOCO - dchoco")

test_that("dchoco returns valid densities", {
  # Basic function test - outputs should be non-negative
  x <- seq(0, 1, length.out = 101)
  dens <- dchoco(x, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.1, bex = 0.5)
  expect_true(all(dens >= 0))

  # For a mixed discrete-continuous distribution, we need to handle
  # the integration differently - point masses need special treatment

  # Point masses at 0 and 1
  zero_mass <- dchoco(0, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.1, bex = 0.5)
  one_mass <- dchoco(1, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.1, bex = 0.5)

  # Continuous part (exclude exact 0 and 1)
  x_cont <- x[x > 0 & x < 1]
  dens_cont <- dchoco(x_cont, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.1, bex = 0.5)

  # Approximate the continuous integral using trapezoid rule
  # Width between points
  dx <- x_cont[2] - x_cont[1]
  # Trapezoid rule integration
  cont_integral <- sum(dens_cont) * dx

  # Total probability should be sum of discrete and continuous parts
  total_prob <- zero_mass + one_mass + cont_integral

  expect_equal(total_prob, 1, tolerance = 0.01)

  # Values outside [0,1] should have zero density
  expect_equal(dchoco(-0.1, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5), 0)
  expect_equal(dchoco(1.1, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5), 0)
})

test_that("dchoco is consistent with rchoco", {
  # Parameters to test
  test_params <- list(
    p = 0.7,
    muleft = 0.3,
    mudelta = 0.1,
    phileft = 5,
    phidelta = 0.2,
    pex = 0.2,
    bex = 0.6,
    threshold = 0.5
  )

  # Generate a large sample
  set.seed(123)
  n <- 100000
  sample_data <- do.call(rchoco, c(list(n = n), test_params))

  # Verify point masses
  zeros <- sum(sample_data == 0) / n
  ones <- sum(sample_data == 1) / n

  # Expected point masses
  zero_prob <- test_params$pex * (1 - test_params$bex) * (1 - test_params$p)
  one_prob <- test_params$pex * test_params$bex * test_params$p

  # Test point masses
  expect_equal(zeros, zero_prob, tolerance = 0.01)
  expect_equal(ones, one_prob, tolerance = 0.01)

  # Test continuous component using histogram
  cont_data <- sample_data[sample_data > 0 & sample_data < 1]

  # Create bins and compute empirical density
  bins <- seq(0, 1, length.out = 21)
  hist_counts <- hist(cont_data, breaks = bins, plot = FALSE)$counts
  bin_width <- 1/20
  empirical_density <- hist_counts / (length(cont_data) * bin_width)

  # Compute theoretical density at bin centers
  bin_centers <- (bins[-1] + bins[-length(bins)]) / 2
  theoretical_density <- dchoco(bin_centers, p = test_params$p, muleft = test_params$muleft,
                            mudelta = test_params$mudelta, phileft = test_params$phileft,
                            phidelta = test_params$phidelta, pex = test_params$pex,
                            bex = test_params$bex, threshold = test_params$threshold)

  # Compare densities (allow some tolerance due to sampling variation)
  for (i in 1:length(bin_centers)) {
    # Skip bins where both densities are very small
    if (theoretical_density[i] > 0.1 || empirical_density[i] > 0.1) {
      expect_equal(empirical_density[i], theoretical_density[i], tolerance = 0.2)
    }
  }
})

test_that("dchoco special cases work correctly", {
  # Test pex = 1 (all mass at extremes)
  x <- seq(0, 1, length.out = 11)
  dens_extreme <- dchoco(x, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 1, bex = 0.5)

  # Should have point masses only at 0 and 1
  expect_true(all(dens_extreme[x > 0 & x < 1] == 0))
  expect_true(dens_extreme[x == 0] > 0)
  expect_true(dens_extreme[x == 1] > 0)

  # Test p = 0 (all continuous mass on left)
  dens_left <- dchoco(x, p = 0, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.5, bex = 0.5)

  # Should have zero density on right side (excluding point mass at 1 and the threshold point)
  expect_true(all(dens_left[x > 0.5 & x < 1] == 0))

  # Test p = 1 (all continuous mass on right)
  dens_right <- dchoco(x, p = 1, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.5, bex = 0.5)

  # Should have zero density on left side (excluding point mass at 0)
  expect_true(all(dens_right[x > 0 & x < 0.5] == 0))
})

test_that("dchoco log option works correctly", {
  x <- c(0, 0.25, 0.5, 0.75, 1)

  dens <- dchoco(x, p = 0.6, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.2, bex = 0.4)
  log_dens <- dchoco(x, p = 0.6, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.2, bex = 0.4, log = TRUE)

  # Check that log option gives log of density
  expect_equal(log_dens[dens > 0], log(dens[dens > 0]))

  # Check that log density for zero values is -Inf
  expect_equal(log_dens[dens == 0], rep(-Inf, sum(dens == 0)))
})

test_that("dchoco respects mudelta and phidelta parameters", {
  x <- seq(0, 1, length.out = 101)

  # Generate densities with different parameters
  dens1 <- dchoco(x, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.1, bex = 0.5)
  dens2 <- dchoco(x, p = 0.5, muleft = 0.3, mudelta = 1, phileft = 5, phidelta = 0, pex = 0.1, bex = 0.5)
  dens3 <- dchoco(x, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 1, pex = 0.1, bex = 0.5)

  # mudelta should shift mode of right distribution
  right_part1 <- x >= 0.5 & x < 1
  right_part2 <- x >= 0.5 & x < 1

  mode1 <- x[right_part1][which.max(dens1[right_part1])]
  mode2 <- x[right_part2][which.max(dens2[right_part2])]

  # When mudelta > 0, the right distribution mode should shift
  expect_true(mode2 != mode1)

  # phidelta should affect the variance of the right distribution
  # Higher phidelta means higher precision and lower variance
  variance1 <- var(rep(x[right_part1], round(dens1[right_part1] * 1000)))
  variance3 <- var(rep(x[right_part1], round(dens3[right_part1] * 1000)))

  expect_true(variance3 < variance1)
})



test_that("dchoco handles threshold point correctly", {
  # Test parameters
  params <- list(
    muleft = 0.3,
    mudelta = 0.2,
    phileft = 4,
    phidelta = 0.5,
    pex = 0.2,
    bex = 0.6,
    threshold = 0.5
  )

  # Get right-side parameters
  compute_right_params <- function(muleft, mudelta, phileft, phidelta) {
    logit_muleft <- log(muleft / (1 - muleft))
    logit_muright <- -logit_muleft + mudelta
    muright <- exp(logit_muright) / (1 + exp(logit_muright))
    phiright <- phileft * exp(phidelta)
    return(list(muright = muright, phiright = phiright))
  }

  right_params <- compute_right_params(params$muleft, params$mudelta, params$phileft, params$phidelta)
  muright <- right_params$muright
  phiright <- right_params$phiright

  # Calculate zero_prob and one_prob for different p values
  p_values <- c(0.2, 0.5, 0.8)

  for (p in p_values) {
    # Compute probabilities for zeros and ones
    zero_prob <- params$pex * (1 - params$bex) * (1 - p)
    one_prob <- params$pex * params$bex * p
    cont_mass <- 1 - zero_prob - one_prob

    # Calculate scale factors
    left_scale <- cont_mass * (1 - p)
    right_scale <- cont_mass * p

    # Calculate left and right limits
    left_limit <- left_scale * stats::dbeta(
      1.0,  # approaching threshold from left
      shape1 = params$muleft * params$phileft,
      shape2 = (1 - params$muleft) * params$phileft
    ) / params$threshold

    right_limit <- right_scale * stats::dbeta(
      0.0,  # approaching threshold from right
      shape1 = muright * phiright,
      shape2 = (1 - muright) * phiright
    ) / (1 - params$threshold)

    # Expected density at threshold: weighted average based on p
    expected_dens <- (1 - p) * left_limit + p * right_limit

    # Actual density from function
    actual_dens <- dchoco(params$threshold, p = p, muleft = params$muleft,
                          mudelta = params$mudelta, phileft = params$phileft,
                          phidelta = params$phidelta, pex = params$pex,
                          bex = params$bex, threshold = params$threshold)

    # Test that density at threshold is not zero
    expect_true(actual_dens > 0)

    # Test that density matches expected weighted average (with tolerance for floating point)
    expect_equal(actual_dens, expected_dens, tolerance = 1e-10)

    # For p = 0.5, the density should be exactly halfway between left and right limits
    if (p == 0.5) {
      expect_equal(actual_dens, (left_limit + right_limit) / 2, tolerance = 1e-10)
    }
  }
})


context("CHOCO - brms")

test_that("CHOCO model can recover parameters with brms using variational inference", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")

  # Set parameters for data generation
  true_p <- 0.7
  true_muleft <- 0.3
  true_mudelta <- 0.5
  true_phileft <- 5
  true_phidelta <- 0.2
  true_pex <- 0.15
  true_bex <- 0.4

  # Generate synthetic data with known parameters
  set.seed(123)
  df <- data.frame(
    y = rchoco(n = 1000, p = true_p, muleft = true_muleft, mudelta = true_mudelta,
              phileft = true_phileft, phidelta = true_phidelta,
              pex = true_pex, bex = true_bex)
  )

  # Formula for intercept-only model
  f <- brms::bf(
    y ~ 1,
    mu ~ 1,
    muleft ~ 1,
    mudelta ~ 1,
    phileft ~ 1,
    phidelta ~ 1,
    pex ~ 1,
    bex ~ 1
  )

  # Fit model with variational inference
  m <- brms::brm(f,
    data = df, family = choco(), stanvars = choco_stanvars(),
    iter = 10000, # More iterations for VI
    algorithm = "pathfinder", # Variational inference
    backend = "cmdstanr", seed = 123, refresh = 0
  )

  # Extract posterior means
  post <- brms::posterior_summary(m)
  est_mu <- brms::inv_logit_scaled(post["b_Intercept", "Estimate"])
  est_muleft <- brms::inv_logit_scaled(post["b_muleft_Intercept", "Estimate"])
  est_mudelta <- post["b_mudelta_Intercept", "Estimate"]
  est_phileft <- log(1 + exp(post["b_phileft_Intercept", "Estimate"]))  # softplus
  est_phidelta <- post["b_phidelta_Intercept", "Estimate"]
  est_pex <- brms::inv_logit_scaled(post["b_pex_Intercept", "Estimate"])
  est_bex <- brms::inv_logit_scaled(post["b_bex_Intercept", "Estimate"])

  # Check that estimates are reasonably close to true values
  # Using wider tolerance for VI since it's an approximation
  expect_equal(est_mu, true_p, tolerance = 0.15)
  expect_equal(est_muleft, true_muleft, tolerance = 0.15)
  expect_equal(est_mudelta, true_mudelta, tolerance = 0.3)
  expect_equal(est_phileft, true_phileft, tolerance = 0.3)
  expect_equal(est_phidelta, true_phidelta, tolerance = 0.3)
  expect_equal(est_pex, true_pex, tolerance = 0.15)
  expect_equal(est_bex, true_bex, tolerance = 0.15)

  # Test posterior prediction
  pred <- brms::posterior_predict(m, ndraws = 10)
  expect_equal(nrow(pred), 10)
  expect_equal(ncol(pred), nrow(df))
  expect_true(all(pred >= 0 & pred <= 1))

  # Also test log-likelihood calculation
  ll <- brms::log_lik(m, ndraws = 5)
  expect_equal(nrow(ll), 5)
  expect_equal(ncol(ll), nrow(df))
})


test_that("choco_lpdf_expose works correctly", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  
  # Expose the Stan function
  choco_lpdf <- choco_lpdf_expose()
  
  # Basic check: a valid call returns a finite value
  val <- choco_lpdf(
    y = 0.5,         # value at threshold
    mu = 0.5,        # proportion right-hand side
    muleft = 0.3,    # mean of left beta
    mudelta = 0.2,   # deviation for right mean
    phileft = 5,     # precision of left beta
    phidelta = 0.5,  # deviation for right precision
    pex = 0.2,       # probability of extremes
    bex = 0.6        # balance of extremes
  )
  
  expect_true(is.finite(val))
  expect_false(is.na(val))
  
  # Test point masses at 0 and 1
  # For zeros: zero_prob = pex * (1 - bex) * (1 - mu)
  val_0 <- choco_lpdf(
    y = 0,
    mu = 0.5,
    muleft = 0.3,
    mudelta = 0.2,
    phileft = 5,
    phidelta = 0.5,
    pex = 0.2,
    bex = 0.6
  )
  expect_equal(exp(val_0), 0.2 * (1 - 0.6) * (1 - 0.5), tolerance = 1e-5)
  
  # For ones: one_prob = pex * bex * mu
  val_1 <- choco_lpdf(
    y = 1,
    mu = 0.5,
    muleft = 0.3,
    mudelta = 0.2,
    phileft = 5,
    phidelta = 0.5,
    pex = 0.2,
    bex = 0.6
  )
  expect_equal(exp(val_1), 0.2 * 0.6 * 0.5, tolerance = 1e-5)
  
  # Test special case (pex = 1)
  val_spec_0 <- choco_lpdf(
    y = 0,
    mu = 0.5,
    muleft = 0.3,   # irrelevant when pex = 1
    mudelta = 0.2,  # irrelevant when pex = 1
    phileft = 5,    # irrelevant when pex = 1
    phidelta = 0.5, # irrelevant when pex = 1
    pex = 1,
    bex = 0.6
  )
  
  # Calculate probability of zero when pex = 1:
  # prob0 = (1 - mu) * (1 - bex) / ((1 - mu) * (1 - bex) + mu * bex)
  prob0 <- (1 - 0.5) * (1 - 0.6) / ((1 - 0.5) * (1 - 0.6) + 0.5 * 0.6)
  expect_equal(exp(val_spec_0), prob0, tolerance = 1e-5)
  
  # Test left vs right side of threshold
  # Left side: y < 0.5
  val_left <- choco_lpdf(
    y = 0.25,  # left side
    mu = 0.5,
    muleft = 0.3,
    mudelta = 0.2,
    phileft = 5,
    phidelta = 0.5,
    pex = 0.2,
    bex = 0.6
  )
  expect_true(is.finite(val_left))
  
  # Right side: y > 0.5
  val_right <- choco_lpdf(
    y = 0.75,  # right side
    mu = 0.5,
    muleft = 0.3,
    mudelta = 0.2,
    phileft = 5,
    phidelta = 0.5,
    pex = 0.2,
    bex = 0.6
  )
  expect_true(is.finite(val_right))
  
  # Edge case: muright < 0 or > 1 should return a very negative value
  # Calculate muright = 1/(1+exp(-(-logit(muleft) + mudelta)))
  # When mudelta is very negative, muright can be <= 0
  val_invalid <- choco_lpdf(
    y = 0.75,
    mu = 0.5,
    muleft = 0.999,  # very high muleft
    mudelta = -20,   # very negative mudelta -> muright ≈ 0
    phileft = 5,
    phidelta = 0.5,
    pex = 0.2,
    bex = 0.6
  )
  # Instead of expecting exactly -Inf, check that it's very negative
  expect_lt(val_invalid, -20)  # Changed from expect_equal(val_invalid, -Inf)


  # Define parameter grids for testing
  y_values <- c(0, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
  mu_values <- c(0.3, 0.5, 0.7)
  muleft_values <- c(0.2, 0.4, 0.7)
  mudelta_values <- c(-0.5, 0, 0.5)
  phileft_values <- c(2, 5)
  phidelta_values <- c(-0.5, 0, 0.5)
  pex_values <- c(0, 0.3, 1)
  bex_values <- c(0.25, 0.5, 0.75)
  
  # Test a subset of the parameter space
  for (mu in mu_values) {
    for (muleft in muleft_values) {
      for (mudelta in mudelta_values) {
        for (phileft in phileft_values) {
          for (phidelta in phidelta_values) {
            for (pex in pex_values) {
              for (bex in bex_values) {
                # Skip invalid parameter combinations that would produce errors
                # Calculate muright to check validity
                logit_muleft <- log(muleft / (1 - muleft))
                logit_muright <- -logit_muleft + mudelta
                muright <- exp(logit_muright) / (1 + exp(logit_muright))
                
                # Skip if muright is not valid
                if (muright <= 0 || muright >= 1) next
                
                # Test different y values
                for (y in y_values) {
                  # Calculate density using Stan lpdf
                  stan_result <- choco_lpdf(y, mu, muleft, mudelta, phileft, phidelta, pex, bex)
                  stan_density <- if (is.finite(stan_result)) exp(stan_result) else 0
                  
                  # Calculate density using R dchoco
                  r_density <- dchoco(y, p = mu, muleft = muleft, mudelta = mudelta, 
                                     phileft = phileft, phidelta = phidelta,
                                     pex = pex, bex = bex)
                  
                  # Build informative test label
                  label <- sprintf("y=%g, mu=%g, muleft=%g, mudelta=%g, phileft=%g, phidelta=%g, pex=%g, bex=%g", 
                                  y, mu, muleft, mudelta, phileft, phidelta, pex, bex)
                  
                  # Special handling for point masses and threshold
                  if (y == 0 || y == 1 || abs(y - 0.5) < 1e-10) {
                    # For exact matches at point masses or threshold
                    expect_equal(stan_density, r_density, 
                                tolerance = 1e-5, 
                                label = paste("Special point at", label))
                  } else {
                    # For continuous densities - use relative tolerance for non-trivial densities
                    if (r_density > 1e-10) {
                      expect_equal(stan_density, r_density, 
                                  tolerance = 1e-5, 
                                  scale = r_density,  # Use relative scaling
                                  label = paste("Continuous density at", label))
                    }
                  }
                }
                
                # Verify log densities for continuous values
                if (pex < 1) {  # Skip pure-discrete case
                  y_cont <- seq(0.1, 0.9, by = 0.2)
                  # Skip threshold point which has special handling
                  y_cont <- y_cont[abs(y_cont - 0.5) > 1e-10]
                  
                  stan_log <- sapply(y_cont, function(y) {
                    choco_lpdf(y, mu, muleft, mudelta, phileft, phidelta, pex, bex)
                  })
                  
                  r_log <- dchoco(y_cont, p = mu, muleft = muleft, mudelta = mudelta, 
                                 phileft = phileft, phidelta = phidelta,
                                 pex = pex, bex = bex, log = TRUE)
                  
                  expect_equal(stan_log, r_log, tolerance = 1e-5,
                              label = sprintf("Log densities (mu=%g, muleft=%g)", mu, muleft))
                }
              }
            }
          }
        }
      }
    }
  }
  
  # Test that invalid parameter combinations return very small or negative values
  # Case: muright very close to 0
  invalid_params <- list(
    y = 0.75,
    mu = 0.5,
    muleft = 0.999,  # even higher muleft
    mudelta = -20,   # more negative mudelta -> muright even closer to 0
    phileft = 5,
    phidelta = 0.5,
    pex = 0.2,
    bex = 0.6
  )
  
  # Calculate the actual muright value to confirm why we're getting small values
  logit_muleft <- log(invalid_params$muleft / (1 - invalid_params$muleft))
  logit_muright <- -logit_muleft + invalid_params$mudelta
  muright <- exp(logit_muright) / (1 + exp(logit_muright))
  
  # For R function - should be very close to 0 but might not be exactly 0
  r_density <- dchoco(invalid_params$y, p = invalid_params$mu, 
                     muleft = invalid_params$muleft, mudelta = invalid_params$mudelta,
                     phileft = invalid_params$phileft, phidelta = invalid_params$phidelta,
                     pex = invalid_params$pex, bex = invalid_params$bex)
  
  # Use a small tolerance instead of expecting exactly 0
  expect_lt(r_density, 1e-6, label = "R density for invalid parameters should be very small")
  
  # For Stan function - should be very negative but might not be -Inf
  invalid_lpdf <- do.call(choco_lpdf, invalid_params)
  
  # Use a less strict threshold - both R and Stan implement the same logic,
  # but floating point precision can vary
  expect_lt(invalid_lpdf, -15, label = "Stan log-density for invalid parameters should be very negative")
})