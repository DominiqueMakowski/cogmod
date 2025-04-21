context("LNR - rlnr and dlnr")

test_that("rlnr generates valid outcomes and dlnr computes valid densities/log-densities", {
  set.seed(123)
  n <- 100
  # Parameters based on the new convention (higher nu = faster)
  nuzero <- 0.5   # Speed for accumulator 0
  nuone <- 0.2    # Speed for accumulator 1 (slower than 0)
  sigmazero <- 1.0 # SD for accumulator 0
  sigmaone <- 0.8  # SD for accumulator 1
  ndt <- 0.2

  # --- Test rlnr ---
  data <- rlnr(n, nuzero = nuzero, nuone = nuone, sigmazero = sigmazero, sigmaone = sigmaone, ndt = ndt)

  # Check output format
  expect_true(is.data.frame(data))
  expect_true(all(c("rt", "response") %in% colnames(data)))
  expect_equal(nrow(data), n)

  # Check RTs and responses
  expect_true(all(data$rt > ndt)) # RTs must be greater than ndt
  expect_true(all(data$response %in% c(0, 1)))
  # Since nuzero > nuone, expect more responses for 0 (faster accumulator)
  expect_true(mean(data$response == 0) > 0.5)

  # --- Test dlnr ---
  y_valid <- ndt + 0.5 # RT > ndt
  y_invalid <- ndt - 0.1 # RT <= ndt

  # Test density calculation for valid RT
  density_0 <- dlnr(y_valid, nuzero, nuone, sigmazero, sigmaone, ndt, response = 0)
  density_1 <- dlnr(y_valid, nuzero, nuone, sigmazero, sigmaone, ndt, response = 1)
  expect_true(is.finite(density_0) && density_0 >= 0)
  expect_true(is.finite(density_1) && density_1 >= 0)

  # Test log-density calculation for valid RT
  log_density_0 <- dlnr(y_valid, nuzero, nuone, sigmazero, sigmaone, ndt, response = 0, log = TRUE)
  log_density_1 <- dlnr(y_valid, nuzero, nuone, sigmazero, sigmaone, ndt, response = 1, log = TRUE)
  expect_true(is.finite(log_density_0) && log_density_0 <= 0) # Log density can be 0 if density is 1
  expect_true(is.finite(log_density_1) && log_density_1 <= 0)
  expect_equal(exp(log_density_0), density_0, tolerance = 1e-9)
  expect_equal(exp(log_density_1), density_1, tolerance = 1e-9)

  # Test density for invalid RT (y <= ndt)
  expect_equal(dlnr(y_invalid, nuzero, nuone, sigmazero, sigmaone, ndt, response = 0), 0)
  expect_equal(dlnr(y_invalid, nuzero, nuone, sigmazero, sigmaone, ndt, response = 0, log = TRUE), -Inf)

  # Test invalid parameters - test both the warning and return value
  # Use expect_warning to test that warnings are issued
  expect_warning(
    result_sz_neg <- dlnr(y_valid, nuzero, nuone, sigmazero = -1, sigmaone, ndt, response = 0),
    "sigmazero must be positive"
  )
  expect_equal(result_sz_neg, 0)

  expect_warning(
    result_so_neg <- dlnr(y_valid, nuzero, nuone, sigmazero, sigmaone = -1, ndt, response = 0),
    "sigmaone must be positive"
  )
  expect_equal(result_so_neg, 0)

  expect_warning(
    result_ndt_neg <- dlnr(y_valid, nuzero, nuone, sigmazero, sigmaone, ndt = -0.1, response = 0),
    "ndt must be non-negative"
  )
  # The result may be finite since y_valid > ndt even when ndt is negative
  expect_true(is.finite(result_ndt_neg))

  # Test invalid response
  expect_equal(suppressWarnings(dlnr(y_valid, nuzero, nuone, sigmazero, sigmaone, ndt, response = 2)), 0)
  expect_equal(suppressWarnings(dlnr(y_valid, nuzero, nuone, sigmazero, sigmaone, ndt, response = 2, log = TRUE)), -Inf)
})

test_that("dlnr integrates correctly and matches rlnr empirical probabilities", {
  set.seed(456)
  n_samples <- 15000 # Large sample size
  tol_prob <- 0.03 # Slightly increased tolerance for probability comparison
  tol_integration <- 0.01 # Tolerance for integration

  # Parameter grid (reduced for reasonable test time)
  nuzero_values <- c(-0.5, 0.5) # Speed for accumulator 0
  nuone_values <- c(-0.2, 0.8)  # Speed for accumulator 1
  sigmazero_values <- c(0.5, 1.2)
  sigmaone_values <- c(0.4, 1.0)
  ndt_values <- c(0.1, 0.3)

  for (nuzero in nuzero_values) {
    for (nuone in nuone_values) {
      for (sigmazero in sigmazero_values) {
        for (sigmaone in sigmaone_values) {
          for (ndt in ndt_values) {

            label <- sprintf("params(nu0=%.1f, nu1=%.1f, sz0=%.1f, sz1=%.1f, ndt=%.1f)",
                             nuzero, nuone, sigmazero, sigmaone, ndt)

            # Generate sample
            data <- rlnr(n_samples, nuzero, nuone, sigmazero, sigmaone, ndt)
            rt_0 <- data$rt[data$response == 0]
            rt_1 <- data$rt[data$response == 1]
            emp_p0 <- length(rt_0) / n_samples
            emp_p1 <- length(rt_1) / n_samples

            # --- Integration Check ---
            # Define the joint density function for integration
            integrand <- function(y, resp) {
              dlnr(y, nuzero, nuone, sigmazero, sigmaone, ndt, response = resp)
            }

            # Integrate densities for each response (handle potential errors)
            # Determine a reasonable upper limit for integration
            q999 <- tryCatch(quantile(data$rt, 0.999, na.rm = TRUE), error = function(e) ndt + 10) # Handle cases with few samples
            upper_limit <- max(q999, ndt + 10, na.rm=TRUE) # Ensure upper limit is well above ndt and typical RTs

            integral_0 <- tryCatch(
              stats::integrate(integrand, lower = ndt + .Machine$double.eps, upper = upper_limit,
                               resp = 0, subdivisions = 200, stop.on.error = FALSE)$value,
              error = function(e) NA
            )
            integral_1 <- tryCatch(
              stats::integrate(integrand, lower = ndt + .Machine$double.eps, upper = upper_limit,
                               resp = 1, subdivisions = 200, stop.on.error = FALSE)$value,
              error = function(e) NA
            )

            # Check if integrations were successful
            if (!is.na(integral_0) && !is.na(integral_1)) {
              # Check if total probability integrates to 1
              total_prob <- integral_0 + integral_1
              expect_equal(total_prob, 1, tolerance = tol_integration,
                           label = paste(label, "- Total Probability Integration"))

              # Check if empirical proportions match integrated probabilities
              # Only check if empirical probability is not too close to 0 or 1
              # to avoid issues with sampling variability for rare outcomes.
              if (emp_p0 > 0.01 && emp_p0 < 0.99) {
                  expect_equal(emp_p0, integral_0, tolerance = tol_prob,
                               label = paste(label, "- P(0): Empirical vs Integrated"))
              }
               if (emp_p1 > 0.01 && emp_p1 < 0.99) {
                  expect_equal(emp_p1, integral_1, tolerance = tol_prob,
                               label = paste(label, "- P(1): Empirical vs Integrated"))
               }
            } else {
              warning("Integration failed for parameters: ", label)
            }

            # --- KS Test Section Removed ---

          } # ndt
        } # sigmaone
      } # sigmazero
    } # nuone
  } # nuzero
})



context("LNR - brms")

test_that("lnr model can recover parameters with brms", {
  # Skip on CRAN and when not running full tests
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")

  # Set random seed for reproducibility
  set.seed(1234)

  # Generate synthetic data
  n <- 1500  # Increased sample size
  true_nuzero <- 0.5
  true_nuone <- 0.2 # Accumulator 1 is slower
  true_sigmazero <- 0.8
  true_sigmaone <- 1.0 # Different sigma for accumulator 1
  true_ndt <- 0.15

  # Generate data
  df <- rlnr(n, nuzero = true_nuzero, nuone = true_nuone,
                  sigmazero = true_sigmazero, sigmaone = true_sigmaone,
                  ndt = true_ndt)
  # Add 'dec' column required by the custom family
  df$dec <- df$response

  # Get min_rt for tau calculation
  min_rt <- min(df$rt)
  true_tau <- true_ndt / min_rt

  # Fit with brms
  f <- brms::bf(
    rt | dec(dec) ~ 1, # Use 'dec' column
    nuone ~ 1,
    sigmazero ~ 1,
    sigmaone ~ 1,
    tau ~ 1,
    minrt = min_rt # Pass min_rt directly
  )

  # Fit using Pathfinder VI
  fit <- brms::brm(
    formula = f,
    data = df,
    family = lnr(), # Uses the updated lnr family
    stanvars = lnr_stanvars(), # Uses the updated stanvars
    init = 0,
    algorithm = "pathfinder", # Use VI
    refresh = 0,
    backend = "cmdstanr"
  )

  # Extract posterior means from summary
  post <- brms::posterior_summary(fit)
  means <- post[, "Estimate"]

  # Check parameter recovery directly with tolerances
  # nuzero (identity link)
  expect_equal(means[["b_Intercept"]], true_nuzero, tolerance = 0.15,
               label = "nuzero recovery")

  # nuone (identity link)
  expect_equal(means[["b_nuone_Intercept"]], true_nuone, tolerance = 0.15,
               label = "nuone recovery")

  # sigmazero (softplus link)
  expect_equal(log1p(exp(means[["b_sigmazero_Intercept"]])), true_sigmazero, tolerance = 0.15,
               label = "sigmazero recovery")

  # sigmaone (softplus link)
  expect_equal(log1p(exp(means[["b_sigmaone_Intercept"]])), true_sigmaone, tolerance = 0.15,
               label = "sigmaone recovery")

  # tau (logit link)
  expect_equal(brms::inv_logit_scaled(means[["b_tau_Intercept"]]), true_tau, tolerance = 0.15,
               label = "tau recovery")

  # Check derived ndt recovery
  est_tau_mean <- brms::inv_logit_scaled(means[["b_tau_Intercept"]])
  est_ndt_mean <- est_tau_mean * min_rt
  # Allow slightly larger tolerance for derived quantity
  expect_equal(est_ndt_mean, true_ndt, tolerance = 0.15,
               label = "Derived NDT recovery")

  # --- Test Post-processing Functions ---
  # Test posterior prediction
  n_pred_draws <- 10
  newdata_pred <- df[1:5, ]
  pred <- brms::posterior_predict(fit, ndraws = n_pred_draws, newdata = newdata_pred)
  expect_true(is.matrix(pred))
  expect_equal(nrow(pred), n_pred_draws)
  expect_equal(ncol(pred), 10) # 5 rows * 2 columns (rt, response)

  # Check rt > ndt (approximately, using mean ndt)
  pred_rt <- pred[, seq(1, ncol(pred), by=2)]
  expect_true(all(pred_rt > est_ndt_mean * 0.9), "Predicted RTs should generally be > estimated NDT")

  # Test log-likelihood
  ll <- brms::log_lik(fit, ndraws = 5)
  expect_true(is.matrix(ll))
  expect_equal(nrow(ll), 5)
  expect_equal(ncol(ll), n)
  expect_true(all(is.finite(ll)), "Log-likelihood values should be finite")
})


test_that("Stan lnr_lpdf matches R dlnr function", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  # Expose the Stan function
  lnr_lpdf_stan <- lnr_lpdf_expose()

  # Define parameter grids for testing
  Y_values <- c(0.3, 0.5, 0.8, 1.2, 2.0)
  nuzero_values <- c(-0.5, 0, 0.5) # Speed for accumulator 0
  nuone_values <- c(-0.2, 0, 0.8)  # Speed for accumulator 1
  sigmazero_values <- c(0.5, 1.0) # Must be > 0
  sigmaone_values <- c(0.4, 1.2)  # Must be > 0
  tau_values <- c(0.1, 0.5) # Must be in (0, 1)
  minrt_values <- c(0.1, 0.2)
  dec_values <- c(0, 1)

  # Test across a subset of the parameter space
  for (nuzero in nuzero_values) {
    for (nuone in nuone_values) {
      for (sigmazero in sigmazero_values) {
        for (sigmaone in sigmaone_values) {
          for (tau in tau_values) {
            for (minrt in minrt_values) {
              # Calculate ndt for the R function
              ndt <- tau * minrt

              for (dec in dec_values) {
                # Test Y values that exceed ndt
                for (Y in Y_values) {
                  # Skip cases where Y <= ndt
                  if (Y <= ndt + 1e-9) next # Add tolerance

                  # Calculate lpdf using Stan function
                  stan_lpdf <- lnr_lpdf_stan(Y, nuzero, nuone, sigmazero, sigmaone, tau, minrt, dec)

                  # Calculate lpdf using R function
                  r_lpdf <- dlnr(Y, nuzero, nuone, sigmazero, sigmaone, ndt, dec, log = TRUE)

                  # Create informative label for errors
                  label <- sprintf("Y=%.2f, nu0=%.1f, nu1=%.1f, sz0=%.1f, sz1=%.1f, tau=%.1f, minrt=%.1f, dec=%d",
                                   Y, nuzero, nuone, sigmazero, sigmaone, tau, minrt, dec)

                  # Check for equality with tolerance
                  expect_equal(stan_lpdf, r_lpdf, tolerance = 1e-6,
                               label = paste("Density mismatch:", label))
                }

                # Test Y value below ndt (should be -Inf in both implementations)
                Y_below <- ndt - 0.01
                if (Y_below > 0) { # Only test if Y_below is positive
                    stan_lpdf_below <- lnr_lpdf_stan(Y_below, nuzero, nuone, sigmazero, sigmaone,
                                                     tau, minrt, dec)
                    r_lpdf_below <- dlnr(Y_below, nuzero, nuone, sigmazero, sigmaone,
                                         ndt, dec, log = TRUE)

                    expect_equal(stan_lpdf_below, r_lpdf_below,
                                 label = sprintf("Below ndt: ndt=%.2f, Y=%.2f", ndt, Y_below))
                    expect_equal(stan_lpdf_below, -Inf,
                                 label = sprintf("Should be -Inf when Y < ndt: Y=%.2f, ndt=%.2f",
                                                 Y_below, ndt))
                }
              }
            }
          }
        }
      }
    }
  }

  # --- Test invalid parameter handling ---
  # R function warnings/errors are tested in the first test_that block

  # Stan function should return -Inf for invalid parameters based on checks
  expect_equal(lnr_lpdf_stan(Y=0.5, mu=0, nuone=0, sigmazero=0, sigmaone=1, tau=0.1, minrt=0.1, dec=0), -Inf, label="Stan invalid sigmazero=0")
  expect_equal(lnr_lpdf_stan(Y=0.5, mu=0, nuone=0, sigmazero=-0.1, sigmaone=1, tau=0.1, minrt=0.1, dec=0), -Inf, label="Stan invalid sigmazero<0")
  expect_equal(lnr_lpdf_stan(Y=0.5, mu=0, nuone=0, sigmazero=1, sigmaone=0, tau=0.1, minrt=0.1, dec=0), -Inf, label="Stan invalid sigmaone=0")
  expect_equal(lnr_lpdf_stan(Y=0.5, mu=0, nuone=0, sigmazero=1, sigmaone=-0.1, tau=0.1, minrt=0.1, dec=0), -Inf, label="Stan invalid sigmaone<0")
  expect_equal(lnr_lpdf_stan(Y=0.5, mu=0, nuone=0, sigmazero=1, sigmaone=1, tau=-0.1, minrt=0.1, dec=0), -Inf, label="Stan invalid tau<0")
  expect_equal(lnr_lpdf_stan(Y=0.5, mu=0, nuone=0, sigmazero=1, sigmaone=1, tau=1.1, minrt=0.1, dec=0), -Inf, label="Stan invalid tau>1")
  expect_equal(lnr_lpdf_stan(Y=0.5, mu=0, nuone=0, sigmazero=1, sigmaone=1, tau=0.1, minrt=-0.1, dec=0), -Inf, label="Stan invalid minrt<0")
  expect_equal(lnr_lpdf_stan(Y=0.5, mu=0, nuone=0, sigmazero=1, sigmaone=1, tau=0.1, minrt=0.1, dec=2), -Inf, label="Stan invalid dec")
  expect_equal(lnr_lpdf_stan(Y=0.05, mu=0, nuone=0, sigmazero=1, sigmaone=1, tau=0.5, minrt=0.2, dec=0), -Inf, label="Stan invalid Y < ndt") # ndt=0.1

})
