context("cogmod_warmstart")

# cogmod_warmstart() moves a fit's adaptation onto another model by labelling
# every unconstrained parameter with its Stan name. Most of what is worth
# testing needs no fitted model: the labels, the join, and the starting
# values can all be checked against a made-up source table that has the
# layout a real fit would have. The one end-to-end check that a real brmsfit
# yields the same layout is behind the slow gate.

set.seed(1)
make_data <- function(ids, n = 30) {
  do.call(rbind, lapply(ids, function(id) {
    d <- rcogmod_lnr(n, nuzero = 0.4, nuone = 0.8, sigmazero = 0.6,
                     sigmaone = 0.6, ndt = 0.2, poutlier = 0.01)
    data.frame(id = id, x = rep(c(0, 1), n / 2), RT = d$rt, choice = d$response)
  }))
}
d_pilot <- make_data(c("p3", "p1"))            # deliberately out of order
d_full <- make_data(c("a1", "p1", "p2", "p3"))  # two new participants, one at each end
f <- brms::bf(RT | dec(choice) ~ x + (1 + x | id), nuone ~ 1 + (1 | id),
              ndt ~ 1, poutlier ~ 1, sigmazero ~ 1, sigmaone ~ 1,
              family = cogmod_lnr())

# A stand-in for a fitted pilot: the labels its program declares, with
# made-up variances and means that are easy to recognise afterwards.
fake_source <- function(formula, data, step_size = 0.1) {
  tab <- cogmod:::.warmstart_target(formula, data)$table
  tab$inv_metric <- seq_len(nrow(tab)) / 10
  tab$mean <- seq_len(nrow(tab)) / 100
  tab$mean[grepl("^L_", tab$parameter)] <- -0.3   # a valid correlation
  tab$step_size <- step_size
  tab
}
labels_of <- function(decl, sdata = list()) {
  d <- cogmod:::.parse_stan_decl(decl)
  dims <- vapply(d$dims, cogmod:::.eval_stan_expr, numeric(1), sdata = sdata)
  cogmod:::.unconstrained_labels(d, as.integer(dims))
}


# Labels --------------------------------------------------------------------

test_that("labels follow Stan's storage order and count the free entries", {
  expect_equal(labels_of("real x"), "x")
  expect_equal(labels_of("vector[Kc] b", list(Kc = 2)), c("b[1]", "b[2]"))
  # a matrix is column-major ...
  expect_equal(labels_of("matrix[M_1, N_1] z_1", list(M_1 = 2, N_1 = 3)),
               c("z_1[1,1]", "z_1[2,1]", "z_1[1,2]", "z_1[2,2]", "z_1[1,3]", "z_1[2,3]"))
  # ... an array of vectors has the array index outermost
  expect_equal(labels_of("array[M_2] vector[N_2] z_2", list(M_2 = 2, N_2 = 3)),
               c("z_2[1,1]", "z_2[1,2]", "z_2[1,3]", "z_2[2,1]", "z_2[2,2]", "z_2[2,3]"))
  expect_equal(labels_of("array[2, 2] real a"), c("a[1,1]", "a[1,2]", "a[2,1]", "a[2,2]"))
  # constrained types have fewer free entries than elements
  expect_equal(labels_of("cholesky_factor_corr[M_1] L_1", list(M_1 = 3)),
               c("L_1[2,1]", "L_1[3,1]", "L_1[3,2]"))
  expect_length(labels_of("simplex[4] s"), 3)
  expect_length(labels_of("cov_matrix[3] S"), 6)
  expect_equal(labels_of("vector[Kc] b", list(Kc = 0)), character(0))
})


test_that("the target table sizes the program and names the group-level entries", {
  tab <- cogmod:::.warmstart_target(f, d_full)$table
  # b (1) + six intercepts + sd_1 (2) + z_1 (2 x 4) + L_1 (1) + sd_2 (1) + z_2 (4)
  expect_equal(nrow(tab), 1 + 6 + 2 + 8 + 1 + 1 + 4)
  z1 <- tab[grepl("^z_1\\[", tab$parameter), ]
  expect_equal(unique(z1$group), "id")
  expect_equal(z1$coef[z1$parameter %in% c("z_1[1,1]", "z_1[2,1]")], c("Intercept", "x"))
  expect_equal(z1$level[grepl("^z_1\\[1,", z1$parameter)], c("a1", "p1", "p2", "p3"))
  expect_equal(tab$coef[tab$parameter == "L_1[2,1]"], "Intercept__x")
  expect_equal(tab$coef[tab$parameter == "sd_2[1]"], "nuone_Intercept")
  expect_true(all(is.na(tab$group[!grepl("^(sd|z|L)_", tab$parameter)])))
})


# Joining --------------------------------------------------------------------

test_that("the same model reproduces its source", {
  src <- fake_source(f, d_pilot)
  ws <- cogmod_warmstart(src, f, d_pilot, jitter = 0)
  expect_s3_class(ws, "cogmod_warmstart")
  expect_equal(ws$inv_metric, src$inv_metric)
  expect_equal(ws$step_size, 0.1)
  expect_true(all(ws$table$source == "pilot"))

  init <- ws$init(1)
  decl <- vapply(cogmod:::.stan_param_decls(
    suppressWarnings(brms::make_stancode(f, data = d_pilot))
  ), `[[`, character(1), "name")
  expect_setequal(names(init), decl)
  expect_equal(init$Intercept, src$mean[src$parameter == "Intercept"])
  expect_equal(init$sd_1, src$mean[grepl("^sd_1\\[", src$parameter)])
  # the standardized effects come back as an M x N matrix in Stan's layout
  expect_equal(dim(init$z_1), c(2, 2))
  expect_equal(init$z_1[2, 1], src$mean[src$parameter == "z_1[2,1]"])
  expect_equal(dim(init$z_2), c(1, 2))
  # the Cholesky factor is rebuilt from its strict lower triangle
  expect_equal(init$L_1[2, 1], -0.3)
  expect_equal(init$L_1[2, 2], sqrt(1 - 0.3^2))
  expect_equal(init$L_1[1, ], c(1, 0))
})


test_that("participants follow their names into a bigger model", {
  src <- fake_source(f, d_pilot)   # levels p1, p3 (sorted by brms)
  expect_silent(ws <- cogmod_warmstart(src, f, d_full, jitter = 0))
  tab <- ws$table
  expect_equal(length(ws$inv_metric), 23)

  # population-level entries carry over one to one
  pop <- is.na(tab$group)
  expect_equal(tab$inv_metric[pop], src$inv_metric[match(tab$parameter[pop], src$parameter)])
  expect_true(all(tab$source[pop] == "pilot"))

  # p1 was level 1 of the pilot and is level 2 of the full model; p3 moves 2 -> 4
  pick <- function(t, p) t$inv_metric[t$parameter == p]
  expect_equal(pick(tab, "z_1[1,2]"), pick(src, "z_1[1,1]"))
  expect_equal(pick(tab, "z_1[2,2]"), pick(src, "z_1[2,1]"))
  expect_equal(pick(tab, "z_1[1,4]"), pick(src, "z_1[1,2]"))
  expect_equal(pick(tab, "z_2[1,4]"), pick(src, "z_2[1,2]"))
  expect_equal(tab$level[tab$parameter == "z_1[1,4]"], "p3")

  # a1 and p2 are new: the coefficient's average variance, and a start at zero
  new <- tab[tab$source == "new level", ]
  expect_setequal(unique(new$level), c("a1", "p2"))
  expect_equal(pick(tab, "z_1[1,1]"), mean(src$inv_metric[src$parameter %in% c("z_1[1,1]", "z_1[1,2]")]))
  expect_equal(pick(tab, "z_1[2,3]"), mean(src$inv_metric[src$parameter %in% c("z_1[2,1]", "z_1[2,2]")]))
  init <- ws$init(1)
  expect_equal(unname(init$z_1[, c(1, 3)]), matrix(0, 2, 2))
  expect_equal(init$z_1[1, 2], src$mean[src$parameter == "z_1[1,1]"])
  expect_output(print(ws), "new group levels")
})


test_that("a changed formula falls back to defaults, with a note", {
  src <- fake_source(f, d_pilot)
  d2 <- d_full
  d2$y <- rnorm(nrow(d2))
  # an extra predictor, and the random slope dropped
  f2 <- brms::bf(RT | dec(choice) ~ x + y + (1 | id), nuone ~ 1 + (1 | id),
                 ndt ~ 1, poutlier ~ 1, sigmazero ~ 1, sigmaone ~ 1,
                 family = cogmod_lnr())
  expect_message(ws <- cogmod_warmstart(src, f2, d2, jitter = 0), "no counterpart")
  tab <- ws$table
  expect_equal(tab$source[tab$parameter == "b[2]"], "default")
  expect_equal(tab$inv_metric[tab$parameter == "b[2]"], 1)
  expect_equal(tab$source[tab$parameter == "b[1]"], "pilot")
  # the intercept's random effect is still the same effect, slope or no slope
  expect_equal(tab$source[tab$parameter == "z_1[1,2]"], "pilot")
  expect_equal(tab$inv_metric[tab$parameter == "z_1[1,2]"], src$inv_metric[src$parameter == "z_1[1,1]"])
  expect_equal(tab$source[tab$parameter == "sd_1[1]"], "pilot")
  expect_false(any(grepl("^L_", tab$parameter)))
  expect_output(print(ws), "without a counterpart")
  init <- ws$init(1)
  expect_equal(init$b[2], 0)   # cogmod_inits()'s generic value for a slope
})


test_that("jitter stays inside the bounds and zero jitter is deterministic", {
  src <- fake_source(f, d_pilot)
  ws0 <- cogmod_warmstart(src, f, d_full, jitter = 0)
  expect_identical(ws0$init(1), ws0$init(2))
  ws <- cogmod_warmstart(src, f, d_full, jitter = 0.5)
  a <- ws$init(1)
  b <- ws$init(2)
  expect_false(identical(a$Intercept, b$Intercept))
  expect_true(all(a$sd_1 > 0) && all(a$sd_2 > 0))
  expect_equal(a$L_1, ws0$init(1)$L_1)   # structured values are not jittered
  expect_setequal(names(a), names(cogmod_inits(f, d_full)(1)))
})


test_that("the table round-trips through a data frame and a CSV file", {
  src <- fake_source(f, d_pilot)
  ws <- cogmod_warmstart(src, f, d_full, jitter = 0)
  df <- as.data.frame(ws)
  expect_named(df, c("parameter", "group", "coef", "level", "inv_metric", "mean", "step_size"))
  again <- cogmod_warmstart(df, f, d_full, jitter = 0)
  expect_equal(again$inv_metric, ws$inv_metric)
  expect_equal(again$table$mean, ws$table$mean)

  tmp <- tempfile(fileext = ".csv")
  utils::write.csv(df, tmp, row.names = FALSE)
  from_file <- cogmod_warmstart(tmp, f, d_full, jitter = 0)
  expect_equal(from_file$inv_metric, ws$inv_metric)
  expect_equal(from_file$init(1), ws$init(1))
  expect_equal(from_file$step_size, 0.1)
})


test_that("it refuses what it cannot use", {
  src <- fake_source(f, d_pilot)
  # a table carries neither the formula nor the data, so both are needed
  expect_error(cogmod_warmstart(src, f), "Give `data`")
  expect_error(cogmod_warmstart(src, data = d_full), "Give `formula`")
  expect_error(cogmod_warmstart(src), "Give `formula` and `data`")
  expect_error(cogmod_warmstart(src[, c("parameter", "mean")], f, d_full), "needs the columns")
  expect_error(cogmod_warmstart(42, f, d_full), "must be a brmsfit")
  expect_error(cogmod_warmstart("no-such-file.csv", f, d_full), "No such file")

  # a fit without adaptation (rstan backend, or an approximation)
  bare <- structure(list(fit = structure(list(), class = "stanfit")), class = "brmsfit")
  expect_error(cogmod_warmstart(bare), "no adaptation")
  # a dense metric
  dense <- bare
  attr(dense$fit, "metadata") <- list(inv_metric = list(diag(3)), step_size = list(0.1))
  expect_error(cogmod_warmstart(dense), "dense")
})


# End to end ------------------------------------------------------------------

test_that("a fitted pilot warm-starts the full model", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")
  skip_if_not_slow()

  pilot <- suppressMessages(brms::brm(
    f, data = d_pilot, prior = cogmod_priors(f, d_pilot), init = cogmod_inits(f, d_pilot),
    stanvars = cogmod_stanvars(f), backend = "cmdstanr", chains = 2, iter = 300,
    warmup = 200, refresh = 0, silent = 2, seed = 1
  ))
  md <- attr(pilot$fit, "metadata")

  # the same model: its own adaptation back, and the posterior means as starts
  ws0 <- cogmod_warmstart(pilot, jitter = 0)
  expect_equal(ws0$inv_metric, unname(Reduce(`+`, md$inv_metric) / 2))
  expect_equal(ws0$step_size, mean(unlist(md$step_size)))
  expect_true(all(ws0$table$source == "pilot"))
  init <- ws0$init(1)
  expect_equal(init$Intercept, mean(brms::as_draws_matrix(pilot)[, "Intercept"]))
  # r = sd * L z reproduces the saved group-level effect
  r <- (diag(init$sd_1) %*% init$L_1 %*% init$z_1)[1, 1]
  expect_equal(r, mean(brms::as_draws_matrix(pilot)[, "r_id[p1,Intercept]"]), tolerance = 0.05)

  # a missing formula means the pilot's own; a missing data set its own
  ws <- cogmod_warmstart(pilot, data = d_full)
  expect_equal(ws$table, cogmod_warmstart(pilot, f, d_full)$table)
  expect_equal(cogmod_warmstart(pilot, formula = f, jitter = 0)$inv_metric, ws0$inv_metric)
  # the full model takes the mapped metric and runs
  expect_equal(sum(ws$table$source == "new level"), 2 * 2 + 2)
  m <- suppressMessages(brms::brm(
    f, data = d_full, prior = cogmod_priors(f, d_full), stanvars = cogmod_stanvars(f),
    init = ws$init, inv_metric = ws$inv_metric, step_size = ws$step_size,
    backend = "cmdstanr", chains = 2, iter = 200, warmup = 100, refresh = 0, silent = 2,
    seed = 2
  ))
  expect_s3_class(m, "brmsfit")
  expect_length(attr(m$fit, "metadata")$inv_metric[[1]], length(ws$inv_metric))
})
