context("RDM")

test_that("rrdm behaves correctly under various conditions", {
  n <- 10000 # For proportion checks
  set.seed(456) # Use a single seed for the combined test

  # --- Parameters ---
  bs_val <- 0.5
  bias_val <- 0.2
  ndt_val <- 0.15
  tol_prop <- 0.05 # Tolerance for proportion checks

  # --- Test 1: Basic Structure and Properties ---
  vzero_t1 <- 0.8
  vone_t1 <- 0.6
  sim_data_t1 <- rrdm(n = 100, vzero = vzero_t1, vone = vone_t1, bs = bs_val, bias = bias_val, ndt = ndt_val)

  # Check structure
  expect_named(sim_data_t1, c("rt", "choice"), label = "Output has correct names")
  expect_equal(nrow(sim_data_t1), 100, label = "Output has correct number of rows")

  # Check values
  expect_true(all(sim_data_t1$rt >= ndt_val), label = "All RTs >= NDT")
  expect_true(all(is.finite(sim_data_t1$rt)), label = "All RTs are finite")
  expect_true(all(sim_data_t1$choice %in% c(1L, 2L)), label = "Choices are 1 or 2")

  # --- Test 2: Choice Proportions Align with Drift Rates ---

  # Case 2.1: vzero > vone
  vzero_21 <- 1.0
  vone_21 <- 0.7
  sim_21 <- rrdm(n = n, vzero = vzero_21, vone = vone_21, bs = bs_val, bias = bias_val, ndt = ndt_val)
  prop_21 <- mean(sim_21$choice == 1)
  expect_gt(prop_21, 0.5, label = sprintf("vzero=%.1f > vone=%.1f: Expect P(choice=1) > 0.5", vzero_21, vone_21))

  # Case 2.2: vone > vzero
  vzero_22 <- 0.6
  vone_22 <- 0.9
  sim_22 <- rrdm(n = n, vzero = vzero_22, vone = vone_22, bs = bs_val, bias = bias_val, ndt = ndt_val)
  prop_22 <- mean(sim_22$choice == 2)
  expect_gt(prop_22, 0.5, label = sprintf("vone=%.1f > vzero=%.1f: Expect P(choice=2) > 0.5", vone_22, vzero_22))

  # Case 2.3: vone = 0
  vzero_23 <- 0.8
  vone_23 <- 0
  sim_23 <- rrdm(n = n, vzero = vzero_23, vone = vone_23, bs = bs_val, bias = bias_val, ndt = ndt_val)
  prop_23 <- mean(sim_23$choice == 1)
  expect_equal(prop_23, 1.0, tolerance = 1e-9, label = sprintf("vone=%.1f: Expect P(choice=1) == 1", vone_23))
  # Check RTs are still finite and >= ndt when one drift is zero
  expect_true(all(sim_23$rt >= ndt_val), label = "All RTs >= NDT (vone=0)")
  expect_true(all(is.finite(sim_23$rt)), label = "All RTs are finite (vone=0)")


  # Case 2.4: vzero = 0
  vzero_24 <- 0
  vone_24 <- 0.8
  sim_24 <- rrdm(n = n, vzero = vzero_24, vone = vone_24, bs = bs_val, bias = bias_val, ndt = ndt_val)
  prop_24 <- mean(sim_24$choice == 2)
  expect_equal(prop_24, 1.0, tolerance = 1e-9, label = sprintf("vzero=%.1f: Expect P(choice=2) == 1", vzero_24))
  # Check RTs are still finite and >= ndt when one drift is zero
  expect_true(all(sim_24$rt >= ndt_val), label = "All RTs >= NDT (vzero=0)")
  expect_true(all(is.finite(sim_24$rt)), label = "All RTs are finite (vzero=0)")

})




test_that("drdm = simulation density & log=TRUE works", {
  set.seed(2025)
  n_sim <- 5000
  pars <- list(vzero=0.9, vone=0.7, bs=0.4, bias=0.15, ndt=0.2)

  # simulate RTs
  dat <- do.call(rrdm, c(list(n=n_sim), pars))

  # empirical density (just above ndt)
  emp <- density(dat$rt, bw="SJ", n=64, from=pars$ndt + 1e-8)
  keep <- emp$x > pars$ndt + 1e-6
  x <- emp$x[keep]; y_emp <- emp$y[keep]

  # theoretical density
  y_th <- do.call(drdm, c(list(x=x), pars))
  expect_gt(cor(y_emp, y_th), 0.9)

  # log argument
  pts <- c(0.25, 0.6, 1.0)
  log_d   <- do.call(drdm, c(list(x=pts, log=TRUE), pars))
  raw_d   <- do.call(drdm, c(list(x=pts, log=FALSE), pars))
  nz      <- raw_d > 0
  expect_equal(log_d[nz], log(raw_d[nz]), tol=1e-8)
  expect_true(all(log_d[!nz] == -Inf))
})


test_that("drdm edge-cases x<=ndt give zero density", {
  p <- list(vzero=0.8, vone=0.6, bs=0.5, bias=0.2, ndt=0.15)
  xs <- c(p$ndt - 0.01, p$ndt, p$ndt + 1e-8)

  dens <- drdm(xs,
               vzero = p$vzero,
               vone  = p$vone,
               bs     = p$bs,
               bias     = p$bias,
               ndt   = p$ndt)

  # first two at/below ndt => 0, third slightly above => >0
  expect_equal(dens[1:2], c(0, 0), label = "density = 0 for x <= ndt")
  expect_gt(dens[3], 0, label = "density > 0 for x > ndt")

  # drdm one‐drift‐zero reduces to dwald()
  p <- list(bs=0.5, bias=0.2, ndt=0.15)
  x0 <- 0.6

  d1 <- drdm(x0,
    vzero = 0.8, vone = 0,
    bs = p$bs, bias = p$bias, ndt = p$ndt
  )
  expect_equal(d1,
    cogmod:::.dwald(x0, drift = 0.8, bs = p$bs, bias = p$bias, ndt = p$ndt),
    tol = 1e-8, label = "drdm with vone=0 reduces to dwald()"
  )

  d2 <- drdm(x0,
    vzero = 0, vone = 0.6,
    bs = p$bs, bias = p$bias, ndt = p$ndt
  )
  expect_equal(d2,
    cogmod:::.dwald(x0, drift = 0.6, bs = p$bs, bias = p$bias, ndt = p$ndt),
    tol = 1e-8, label = "drdm with vzero=0 reduces to dwald()"
  )


  # drdm errors on invalid inputs
  # negative drifts
  expect_error(drdm(0.5, vzero=-0.1, vone=0.5, bs=0.5, bias=0.2, ndt=0.1))
  expect_error(drdm(0.5, vzero=0.5,  vone=-0.1, bs=0.5, bias=0.2, ndt=0.1))
  expect_error(drdm(0.5, vzero=0,     vone=0,    bs=0.5, bias=0.2, ndt=0.1))
  # non-positive bs or bias
  expect_error(drdm(0.5, vzero=0.5, vone=0.5, bs=0,   bias=0.2, ndt=0.1))
  expect_error(drdm(0.5, vzero=0.5, vone=0.5, bs=-1,  bias=0.2, ndt=0.1))
  expect_error(drdm(0.5, vzero=0.5, vone=0.5, bs=0.5, bias=0,   ndt=0.1))
  expect_error(drdm(0.5, vzero=0.5, vone=0.5, bs=0.5, bias=-1,  ndt=0.1))
  # negative ndt
  expect_error(drdm(0.5, vzero=0.5, vone=0.5, bs=0.5, bias=0.2, ndt=-0.1))
})
