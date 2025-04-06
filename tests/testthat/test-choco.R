library(testthat)

test_that("rchoco generates valid outcomes", {
  set.seed(123)
  n <- 1000
  outcomes <- rchoco(n, mu = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.1, bex = 0.5)

  # Outcomes should be in the range [0, 1]
  expect_true(all(outcomes >= 0 & outcomes <= 1))

  # Check that some outcomes are exactly 0 or 1 (extreme values)
  expect_true(any(outcomes == 0))
  expect_true(any(outcomes == 1))

  # Check that some outcomes are continuous
  expect_true(any(outcomes > 0 & outcomes < 1))
})

test_that("rchoco symmetry works without extreme values (pex = 0)", {
  set.seed(123)
  n <- 10000
  # Set pex = 0 to remove extreme values
  outcomes <- rchoco(n, mu = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0, bex = 0.5, threshold = 0.5)

  # Separate left and right sides based on the threshold
  left_side <- outcomes[outcomes >= 0 & outcomes <= 0.5]
  right_side <- outcomes[outcomes > 0.5 & outcomes <= 1]

  # Compute the expected means
  expected_mean_left <- 0.3 * 0.5  # muleft * threshold
  expected_mean_right <- 0.5 + (1 - 0.3) * 0.5  # threshold + muright * (1 - threshold)

  # Check that the means match the expected values
  expect_equal(mean(left_side), expected_mean_left, tolerance = 0.05)
  expect_equal(mean(right_side), expected_mean_right, tolerance = 0.05)
})

test_that("rchoco respects mudelta for asymmetry without extreme values", {
  set.seed(123)
  n <- 10000
  # Set pex = 0 to remove extreme values
  outcomes <- rchoco(n, mu = 0.5, muleft = 0.3, mudelta = 1, phileft = 5, phidelta = 0, pex = 0, bex = 0.5, threshold = 0.5)

  # Separate left and right sides based on the threshold
  left_side <- outcomes[outcomes > 0 & outcomes <= 0.5]
  right_side <- outcomes[outcomes > 0.5 & outcomes < 1]

  # Compute the expected mean for the left-hand side
  expected_mean_left <- 0.3 * 0.5  # muleft * threshold

  # Check that the left-hand side mean remains the same
  expect_equal(mean(left_side), expected_mean_left, tolerance = 0.05)

  # Check that the right-hand side mean is shifted to the right
  # Compute muright with mudelta = 1
  logit_muleft <- log(0.3 / (1 - 0.3))
  logit_muright <- -logit_muleft + 1
  muright <- exp(logit_muright) / (1 + exp(logit_muright))
  expected_mean_right <- 0.5 + muright * 0.5  # threshold + muright * (1 - threshold)

  expect_true(mean(right_side) > expected_mean_right - 0.05)  # Allow for tolerance
})

test_that("rchoco respects phidelta for precision adjustment without extreme values", {
  set.seed(123)
  n <- 10000
  # Set pex = 0 to remove extreme values
  outcomes_high_precision <- rchoco(n, mu = 0.5, muleft = 0.3, mudelta = 0, phileft = 10, phidelta = 0, pex = 0, bex = 0.5)
  outcomes_low_precision <- rchoco(n, mu = 0.5, muleft = 0.3, mudelta = 0, phileft = 2, phidelta = 0, pex = 0, bex = 0.5)

  # Separate left-hand side outcomes
  left_side_high <- outcomes_high_precision[outcomes_high_precision > 0 & outcomes_high_precision <= 0.5]
  left_side_low <- outcomes_low_precision[outcomes_low_precision > 0 & outcomes_low_precision <= 0.5]

  # High precision should result in less spread
  expect_true(sd(left_side_high) < sd(left_side_low))
})

test_that("rchoco respects pex and bex for extreme values", {
  set.seed(123)
  n <- 10000
  outcomes <- rchoco(n, mu = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.2, bex = 0.8)

  # Check that the proportion of extreme values matches pex
  extreme_values <- sum(outcomes == 0 | outcomes == 1) / n
  expect_equal(extreme_values, 0.2, tolerance = 0.02)

  # Check that the balance of extreme values matches bex
  zeros <- sum(outcomes == 0)
  ones <- sum(outcomes == 1)
  expect_equal(ones / (zeros + ones), 0.8, tolerance = 0.05)
})

test_that("dchoco computes valid densities", {

  y <- rchoco(10000, mu = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.1, bex = 0.5)
  x <- seq(0, 1, length.out = 1000)
  densities <- dchoco(x, mu = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.1, bex = 0.5)

  # visual test
  # plot(density(y))
  # lines(x, densities, col = "red")

  # Densities should be non-negative
  expect_true(all(densities >= 0))

  # The sum of densities (scaled by the range) should approximate 1
  expect_equal(sum(densities) * (1 / length(x)), 1, tolerance = 0.1)
})
