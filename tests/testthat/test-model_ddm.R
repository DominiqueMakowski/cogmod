context("DDM - rcogmod_ddm and dcogmod_ddm")


test_that("rcogmod_ddm and dcogmod_ddm produces correct output shape", {
  n <- 1000
  drift <- 0.2
  boundary <- 1
  bias <- 0.5
  ndt <- 0.3

  # Run rcogmod_ddm
  sim_data <- rcogmod_ddm(n, drift = drift, boundary = boundary, bias = bias, ndt = ndt)

  # Check that the output is a data frame
  expect_s3_class(sim_data, "data.frame")

  # Check that the data frame has the correct number of rows
  expect_equal(nrow(sim_data), n)

  # Check that the data frame has the correct columns
  expect_named(sim_data, c("rt", "response"))

  # DDM
  x <- c(0.5, 0.7, 1.2)
  drift <- 0.2
  boundary <- 1
  bias <- 0.5
  ndt <- 0.3
  resp <- c(0, 1, 0)

  # Run dcogmod_ddm
  densities <- dcogmod_ddm(
    x = x,
    drift = drift,
    boundary = boundary,
    bias = bias,
    ndt = ndt,
    response = resp
  )

  # Check that the output is a numeric vector
  expect_type(densities, "double")

  # Check that the output has the correct length
  expect_equal(length(densities), length(x))
})


test_that("dcogmod_ddm matches densities of rcogmod_ddm simulated data", {
  set.seed(123)

  # Parameters for the test
  n <- 1000
  drift <- 0.5
  boundary <- 1.2
  bias <- 0.6
  ndt <- 0.3

  # Simulate data using rcogmod_ddm
  sim_data <- rcogmod_ddm(n, drift = drift, boundary = boundary, bias = bias, ndt = ndt)

  # Extract reaction times and responses
  rt <- sim_data$rt
  response <- sim_data$response

  # Compute densities using dcogmod_ddm
  densities <- dcogmod_ddm(
    x = rt,
    drift = drift,
    boundary = boundary,
    bias = bias,
    ndt = ndt,
    response = response,
    log = FALSE
  )

  # Check that densities are finite and non-negative
  expect_true(all(is.finite(densities)), "All densities should be finite")
  expect_true(all(densities >= 0), "All densities should be non-negative")

  # Define plausible RTs as those within a reasonable central range, excluding potential tails
  # Original upper bound (ndt + boundary / abs(drift)) might include all samples for some seeds/params.
  # Use a tighter upper bound to better isolate the central tendency.
  plausible_upper_bound <- ndt + boundary / (2 * abs(drift)) # e.g., ndt + half the mean decision time approx.
  plausible_rt_indices <- rt > ndt & rt < plausible_upper_bound
  plausible_densities <- densities[plausible_rt_indices]

  # Check that plausible RTs have higher average density than the overall average density
  expect_gt(
    mean(plausible_densities),
    mean(densities),
    "Mean density of plausible RTs should be higher than overall mean density"
  )
})

context("DDM - brms")

test_that("DDM model can recover parameters with brms", {
  skip_if_not_slow()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")

  set.seed(1234)
  n_obs <- 3000

  # True parameters
  true_nu <- 0.5
  true_boundary <- 1.2
  true_bias <- 0.6
  true_tau <- 0.8 # Proportion of minrt
  true_minrt <- 0.3
  true_ndt <- true_tau * true_minrt
  true_minrt <- true_minrt

  # Simulate data
  df <- rcogmod_ddm(
    n_obs,
    drift = true_nu,
    boundary = true_boundary,
    bias = true_bias,
    ndt = true_ndt
  )

  # Define brms formula
  f <- brms::bf(
    rt | dec(response) ~ 1,
    boundary ~ 1,
    bias ~ 1,
    tau ~ 1,
    minrt = true_minrt,
    sigmadrift = 0,
    sigmabias = 0,
    sigmatau = 0
  )

  # Fit model using brms
  fit <- brms::brm(
    formula = f,
    data = df,
    family = cogmod_ddm(),
    stanvars = cogmod_ddm_stanvars(),
    algorithm = "pathfinder", # Use VI
    backend = "cmdstanr",
    refresh = 0,
    init = 0
  )

  # Extract posterior means from summary
  post <- brms::posterior_summary(fit)
  means <- post[, "Estimate"]

  # Check parameter recovery
  expect_equal(
    means[["b_Intercept"]],
    true_nu,
    tolerance = 0.1,
    label = "drift recovery"
  )
  expect_equal(
    log1p(exp(means[["b_boundary_Intercept"]])),
    true_boundary,
    tolerance = 0.1,
    label = "boundary recovery"
  )
  expect_equal(
    brms::inv_logit_scaled(means[["b_bias_Intercept"]]),
    true_bias,
    tolerance = 0.1,
    label = "bias recovery"
  )
  expect_equal(
    brms::inv_logit_scaled(means[["b_tau_Intercept"]]),
    true_tau,
    tolerance = 0.1,
    label = "tau recovery"
  )

  # Check derived ndt recovery
  est_tau_mean <- brms::inv_logit_scaled(means[["b_tau_Intercept"]])
  est_ndt_mean <- est_tau_mean * min(df$rt)
  expect_equal(
    est_ndt_mean,
    true_ndt,
    tolerance = 0.2,
    label = "Derived NDT recovery"
  )

  # --- Test Post-processing Functions ---
  # Test posterior prediction dimensions
  n_pred_draws <- 10
  n_pred_obs <- 5
  # Ensure newdata has the 'response' column if needed by posterior_predict_cogmod_ddm
  # (It's not strictly needed for prediction itself, but good practice)
  newdata_pred <- df[1:n_pred_obs, , drop = FALSE]

  # Generate predictions
  pred <- brms::posterior_predict(
    fit,
    ndraws = n_pred_draws,
    newdata = newdata_pred
  )

  # Check dimensions: posterior_predict for wiener returns draws x (obs*2) matrix
  # Columns alternate between RT (q) and response (resp)
  expect_true(
    is.matrix(pred),
    label = "posterior_predict output should be a matrix"
  )
  expect_equal(
    nrow(pred),
    n_pred_draws,
    label = "posterior_predict rows should match ndraws"
  )
  expect_equal(
    ncol(pred),
    n_pred_obs * 2,
    label = "posterior_predict columns should be 2 * nrow(newdata)"
  )

  # Optional: Check predicted RTs are plausible (e.g., > estimated NDT)
  pred_rt_indices <- seq(1, ncol(pred), by = 2) # Indices for RT columns
  pred_rt <- pred[, pred_rt_indices]
  expect_true(
    all(pred_rt > est_ndt_mean * 0.9, na.rm = TRUE), # Allow some buffer
    label = "Predicted RTs should generally be > estimated NDT"
  )

  # Optional: Check predicted responses are plausible (0 or 1)
  pred_resp_indices <- seq(2, ncol(pred), by = 2) # Indices for response columns
  pred_resp <- pred[, pred_resp_indices]
  expect_true(
    all(pred_resp %in% c(0, 1)),
    label = "Predicted responses should be 0 or 1"
  )

  # Test log-likelihood dimensions
  ll <- brms::log_lik(fit, ndraws = 5) # Use a small number of draws for testing
  expect_true(is.matrix(ll), label = "log_lik output should be a matrix")
  expect_equal(nrow(ll), 5, label = "log_lik rows should match ndraws")
  expect_equal(
    ncol(ll),
    n_obs,
    label = "log_lik columns should match number of original observations"
  )
  expect_true(all(is.finite(ll)), "Log-likelihood values should be finite")
})


test_that("Stan DDM lpdf matches R dcogmod_ddm function", {
  skip_if_not_installed("cmdstanr")

  # Expose the Stan function
  cogmod_ddm_lpdf_stan <- stan_fun("cogmod_ddm")

  # Define parameter grids for testing
  Y_values <- c(0.5, 0.7, 1.0, 1.5)
  drift_values <- c(0.2, 0.5)
  boundary_values <- c(1.0, 1.5)
  bias_values <- c(0.3, 0.6)
  tau_values <- c(0.5, 0.8)
  minrt_values <- c(0.2, 0.3)
  dec_values <- c(0, 1)

  for (drift in drift_values) {
    for (boundary in boundary_values) {
      for (bias in bias_values) {
        for (tau in tau_values) {
          for (minrt in minrt_values) {
            ndt <- tau * minrt

            for (dec in dec_values) {
              for (Y in Y_values) {
                if (Y <= ndt) {
                  next
                } # Skip invalid cases where Y <= ndt

                # Calculate lpdf using Stan function (extra variability
                # params fixed to 0, should match the classic 4-param DDM)
                stan_lpdf <- cogmod_ddm_lpdf_stan(
                  Y,
                  drift,
                  boundary,
                  bias,
                  tau,
                  minrt,
                  0,
                  0,
                  0,
                  dec
                )

                # Calculate lpdf using R function
                r_lpdf <- dcogmod_ddm(
                  x = Y,
                  drift = drift,
                  boundary = boundary,
                  bias = bias,
                  ndt = ndt,
                  response = dec,
                  log = TRUE
                )

                # Compare results
                label <- sprintf(
                  "Y=%.2f, drift=%.2f, boundary=%.2f, bias=%.2f, tau=%.2f, minrt=%.2f, dec=%d",
                  Y,
                  drift,
                  boundary,
                  bias,
                  tau,
                  minrt,
                  dec
                )
                expect_equal(stan_lpdf, r_lpdf, tolerance = 1e-6, label = label)
              }
            }
          }
        }
      }
    }
  }
})


context("DDM - full (rtdists-backed) 7-parameter model")

test_that("dcogmod_ddm() with variability matches rtdists::ddiffusion()", {
  skip_on_cran()
  skip_if_not_installed("rtdists")

  grid <- expand.grid(
    drift = c(0.2, -0.3),
    boundary = c(1.0, 1.5),
    bias = c(0.4, 0.6),
    ndt = c(0.15, 0.25),
    sigmadrift = c(0, 0.3),
    sigmabias = c(0, 0.4),
    sigmatau = c(0, 0.3),
    minrt = c(0.2, 0.3),
    dec = c(0, 1),
    Y = c(0.6, 0.9, 1.3)
  )

  for (i in seq_len(nrow(grid))) {
    row <- grid[i, ]
    if (row$Y <= row$ndt) {
      next
    }

    r_lpdf <- dcogmod_ddm(
      x = row$Y,
      drift = row$drift,
      boundary = row$boundary,
      bias = row$bias,
      ndt = row$ndt,
      response = row$dec,
      log = TRUE,
      sigmadrift = row$sigmadrift,
      sigmabias = row$sigmabias,
      sigmatau = row$sigmatau,
      minrt = row$minrt
    )

    sw <- row$sigmabias * min(2 * row$bias, 2 * (1 - row$bias))
    ref_density <- rtdists::ddiffusion(
      rt = row$Y,
      response = row$dec + 1,
      a = row$boundary,
      v = row$drift,
      t0 = row$ndt,
      z = row$bias * row$boundary,
      sv = row$sigmadrift,
      sz = sw * row$boundary,
      st0 = row$sigmatau * row$minrt
    )

    label <- sprintf(
      "drift=%.2f, boundary=%.2f, bias=%.2f, ndt=%.2f, sigmadrift=%.2f, sigmabias=%.2f, sigmatau=%.2f, dec=%d, Y=%.2f",
      row$drift,
      row$boundary,
      row$bias,
      row$ndt,
      row$sigmadrift,
      row$sigmabias,
      row$sigmatau,
      row$dec,
      row$Y
    )
    # Tolerance accounts for rtdists' default numerical-integration precision
    # (precision = 3), which introduces an error on the order of ~3e-4.
    expect_equal(r_lpdf, log(ref_density), tolerance = 5e-4, label = label)
  }
})


test_that("Stan cogmod_ddm_lpdf matches dcogmod_ddm() when variability parameters are non-zero", {
  skip_if_not_installed("cmdstanr")
  skip_if_not_installed("rtdists")

  cogmod_ddm_lpdf_stan <- stan_fun("cogmod_ddm")

  grid <- expand.grid(
    drift = c(0.2, -0.3),
    boundary = c(1.0, 1.5),
    bias = c(0.4, 0.6),
    tau = c(0.5, 0.8),
    minrt = c(0.2, 0.3),
    sigmadrift = c(0, 0.3),
    sigmabias = c(0, 0.4),
    sigmatau = c(0, 0.3),
    dec = c(0, 1),
    Y = c(0.6, 0.9, 1.3)
  )

  for (i in seq_len(nrow(grid))) {
    row <- grid[i, ]
    ndt <- row$tau * row$minrt
    if (row$Y <= ndt) {
      next
    }

    stan_lpdf <- cogmod_ddm_lpdf_stan(
      row$Y,
      row$drift,
      row$boundary,
      row$bias,
      row$tau,
      row$minrt,
      row$sigmadrift,
      row$sigmabias,
      row$sigmatau,
      row$dec
    )
    r_lpdf <- dcogmod_ddm(
      x = row$Y,
      drift = row$drift,
      boundary = row$boundary,
      bias = row$bias,
      ndt = ndt,
      response = row$dec,
      log = TRUE,
      sigmadrift = row$sigmadrift,
      sigmabias = row$sigmabias,
      sigmatau = row$sigmatau,
      minrt = row$minrt
    )

    label <- sprintf(
      "drift=%.2f, boundary=%.2f, bias=%.2f, tau=%.2f, minrt=%.2f, sigmadrift=%.2f, sigmabias=%.2f, sigmatau=%.2f, dec=%d, Y=%.2f",
      row$drift,
      row$boundary,
      row$bias,
      row$tau,
      row$minrt,
      row$sigmadrift,
      row$sigmabias,
      row$sigmatau,
      row$dec,
      row$Y
    )
    # Tolerance accounts for rtdists' default numerical-integration precision
    # (precision = 3), which introduces an error on the order of ~3e-4.
    expect_equal(stan_lpdf, r_lpdf, tolerance = 5e-4, label = label)
  }
})


test_that("rcogmod_ddm() with variability produces valid simulated data", {
  skip_if_not_installed("rtdists")

  sim <- rcogmod_ddm(
    n = 500,
    drift = 0.4,
    boundary = 1.2,
    bias = 0.5,
    ndt = 0.2,
    sigmadrift = 0.3,
    sigmabias = 0.3,
    sigmatau = 0.2,
    minrt = 0.25
  )

  expect_s3_class(sim, "data.frame")
  expect_equal(nrow(sim), 500)
  expect_named(sim, c("rt", "response"))
  expect_true(all(sim$response %in% c(0, 1)))
  expect_true(all(is.finite(sim$rt)))
  expect_true(all(sim$rt > 0))
})


test_that("rcogmod_ddm() with variability samples the distribution dcogmod_ddm() describes", {
  # The simulator draws the per-trial parameters while the density integrates
  # them out analytically/by quadrature, so the two take genuinely different
  # routes and agreeing is informative.
  skip_on_cran()
  set.seed(404)

  p <- list(drift = -1.2, boundary = 1.3, bias = 0.4, ndt = 0.25, minrt = 0.30)
  grid <- seq(p$ndt + 1e-4, 12, length.out = 300)
  h <- diff(grid)
  trap <- function(y) sum((y[-1] + y[-length(y)]) / 2 * h)

  for (v in list(c(0.8, 0, 0), c(0, 0.5, 0), c(0, 0, 0.5), c(0.8, 0.4, 0.4))) {
    dens <- lapply(c(0, 1), function(r) {
      dcogmod_ddm(
        x = grid, drift = p$drift, boundary = p$boundary, bias = p$bias, ndt = p$ndt,
        response = r, sigmadrift = v[1], sigmabias = v[2], sigmatau = v[3],
        minrt = p$minrt
      )
    })
    sim <- rcogmod_ddm(
      n = 15000, drift = p$drift, boundary = p$boundary, bias = p$bias, ndt = p$ndt,
      sigmadrift = v[1], sigmabias = v[2], sigmatau = v[3], minrt = p$minrt
    )

    lab <- sprintf("sv=%.1f, sigmabias=%.1f, sigmatau=%.1f", v[1], v[2], v[3])
    p_lo <- trap(dens[[1]])
    p_up <- trap(dens[[2]])

    # A proper density over both responses.
    expect_equal(p_lo + p_up, 1, tolerance = 0.005, label = paste("mass:", lab))
    expect_equal(
      p_up, mean(sim$response == 1),
      tolerance = 0.02, label = paste("P(response = 1):", lab)
    )
    for (r in c(0, 1)) {
      expect_equal(
        trap(grid * dens[[r + 1]]) / c(p_lo, p_up)[r + 1],
        mean(sim$rt[sim$response == r]),
        tolerance = 0.02, label = sprintf("E[RT | response = %d]: %s", r, lab)
      )
    }
  }
})


test_that("posterior_epred_cogmod_ddm() matches the simulated mean RT", {
  skip_on_cran()
  # The closed-form mean first-passage time takes the starting point on the
  # absolute scale (z = bias * boundary), not as the proportion `bias`. Using the
  # proportion silently biases E[RT] whenever boundary != 1.
  set.seed(99)

  grid <- expand.grid(
    mu = c(-2.5, -1, -0.4, 0.4, 1, 2.5),
    boundary = c(0.8, 1.4, 2.0),
    bias = c(0.35, 0.5, 0.65)
  )

  for (i in seq_len(nrow(grid))) {
    g <- grid[i, ]
    prep <- list(dpars = list(
      mu = g$mu, boundary = g$boundary, bias = g$bias, tau = 0.8, minrt = 0.3125,
      sigmadrift = 0, sigmabias = 0, sigmatau = 0
    ))
    sim <- rcogmod_ddm(15000, drift = g$mu, boundary = g$boundary, bias = g$bias, ndt = 0.25)
    expect_equal(
      posterior_epred_cogmod_ddm(prep),
      mean(sim$rt),
      tolerance = 0.02,
      label = sprintf("mu=%.1f, boundary=%.1f, bias=%.2f", g$mu, g$boundary, g$bias)
    )
  }

  # The mu = 0 limit (removable singularity) is ndt + z * (boundary - z).
  prep0 <- list(dpars = list(
    mu = 0, boundary = 1.4, bias = 0.35, tau = 0.8, minrt = 0.3125,
    sigmadrift = 0, sigmabias = 0, sigmatau = 0
  ))
  z <- 0.35 * 1.4
  expect_equal(posterior_epred_cogmod_ddm(prep0), 0.25 + z * (1.4 - z), tolerance = 1e-8)
})


test_that("posterior_epred_cogmod_ddm() warns only when variability parameters are non-zero", {
  mock_prep <- function(sigmadrift, sigmabias, sigmatau) {
    list(
      dpars = list(
        mu = c(0.5, 0.5),
        boundary = c(1.2, 1.2),
        bias = c(0.5, 0.5),
        tau = c(0.7, 0.7),
        minrt = c(0.25, 0.25),
        sigmadrift = sigmadrift,
        sigmabias = sigmabias,
        sigmatau = sigmatau
      )
    )
  }

  expect_no_warning(posterior_epred_cogmod_ddm(mock_prep(c(0, 0), c(0, 0), c(0, 0))))
  expect_warning(
    posterior_epred_cogmod_ddm(mock_prep(c(0.3, 0.3), c(0, 0), c(0, 0))),
    "closed-form approximation"
  )
  expect_warning(
    posterior_epred_cogmod_ddm(mock_prep(c(0, 0), c(0.2, 0.2), c(0, 0))),
    "closed-form approximation"
  )
  expect_warning(
    posterior_epred_cogmod_ddm(mock_prep(c(0, 0), c(0, 0), c(0.1, 0.1))),
    "closed-form approximation"
  )
})
