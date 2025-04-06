library(testthat)

test_that("rbext generates values in the range [0, 1]", {
  set.seed(123)
  n <- 1000
  results <- rbext(n, mu = 0.5, phi = 3, pex = 0.1, bex = 0.5)
  expect_true(all(results >= 0 & results <= 1))
})

test_that("rbext generates approximately correct proportions of zeros and ones", {
  set.seed(123)
  n <- 10000
  pex <- 0.2
  bex <- 0.75
  results <- rbext(n, mu = 0.5, phi = 3, pex = pex, bex = bex)
  
  # Expected proportions
  expected_zeros <- pex * (1 - bex)
  expected_ones <- pex * bex
  
  # Observed proportions
  observed_zeros <- mean(results == 0)
  observed_ones <- mean(results == 1)
  
  expect_equal(observed_zeros, expected_zeros, tolerance = 0.02)
  expect_equal(observed_ones, expected_ones, tolerance = 0.02)
})

test_that("rbext generates continuous values when pex = 0", {
  set.seed(123)
  n <- 1000
  results <- rbext(n, mu = 0.5, phi = 3, pex = 0, bex = 0.5)
  expect_true(all(results > 0 & results < 1))
})

test_that("rbext generates only zeros when pex = 1 and bex = 0", {
  set.seed(123)
  n <- 1000
  results <- rbext(n, mu = 0.5, phi = 3, pex = 1, bex = 0)
  expect_true(all(results == 0))
})

test_that("rbext generates only ones when pex = 1 and bex = 1", {
  set.seed(123)
  n <- 1000
  results <- rbext(n, mu = 0.5, phi = 3, pex = 1, bex = 1)
  expect_true(all(results == 1))
})

test_that("rbext handles edge cases for mu and phi", {
  set.seed(123)
  n <- 1000
  
  # Test with mu close to 0
  results <- rbext(n, mu = 0.01, phi = 3, pex = 0.1, bex = 0.5)
  expect_true(all(results >= 0 & results <= 1))
  
  # Test with mu close to 1
  results <- rbext(n, mu = 0.99, phi = 3, pex = 0.1, bex = 0.5)
  expect_true(all(results >= 0 & results <= 1))
  
  # Test with large phi
  results <- rbext(n, mu = 0.5, phi = 100, pex = 0.1, bex = 0.5)
  expect_true(all(results >= 0 & results <= 1))
})

test_that("rbext throws an error for invalid inputs", {
  expect_error(rbext(100, mu = -0.1, phi = 3, pex = 0.1, bex = 0.5), "mu must be between 0 and 1")
  expect_error(rbext(100, mu = 1.1, phi = 3, pex = 0.1, bex = 0.5), "mu must be between 0 and 1")
  expect_error(rbext(100, mu = 0.5, phi = -1, pex = 0.1, bex = 0.5), "phi must be positive")
  expect_error(rbext(100, mu = 0.5, phi = 3, pex = -0.1, bex = 0.5), "pex must be between 0 and 1")
  expect_error(rbext(100, mu = 0.5, phi = 3, pex = 1.1, bex = 0.5), "pex must be between 0 and 1")
  expect_error(rbext(100, mu = 0.5, phi = 3, pex = 0.1, bex = -0.1), "bex must be between 0 and 1")
  expect_error(rbext(100, mu = 0.5, phi = 3, pex = 0.1, bex = 1.1), "bex must be between 0 and 1")
})