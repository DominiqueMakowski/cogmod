library(cmdstanr)
dir <- "C:/Users/domma/AppData/Local/Temp/claude/C--Users-domma-Dropbox-Software-cogmod/b3646c2f-5a2d-4129-99b6-b665de42a8c1/scratchpad"
set.seed(1); N <- 200
y <- 0.25 + rgamma(N, shape = 2, rate = 4)
mod <- cmdstan_model(file.path(dir, "wienerbench.stan"))

run <- function(label, branch, sw = c(1e-6, 2e-6), st0 = c(1e-6, 2e-6), prec = 1e-4) {
  d <- list(N = N, y = y, branch = branch, sw_lo = sw[1], sw_hi = sw[2],
            st0_lo = st0[1], st0_hi = st0[2], prec = prec)
  ini <- list(list(v = 1, a = 1.2, w = 0.5, sv = 0.3,
                   sw = mean(sw), st0 = mean(st0)))
  f <- mod$sample(data = d, chains = 1, iter_warmup = 0, iter_sampling = 15,
                  adapt_engaged = FALSE, step_size = 5e-4, max_treedepth = 5,
                  init = ini, seed = 42, refresh = 0, show_messages = FALSE,
                  output_dir = dir)
  us <- 1e6 * f$time()$total / sum(f$sampler_diagnostics()[, , "n_leapfrog__"])
  cat(sprintf("%-44s | %9.1f us/grad | %6.2f us/obs\n", label, us, us / N)); invisible(us)
}

for (rep in 1:2) for (prec in c(1e-4, 2e-4, 3e-4, 5e-4, 1e-3)) {
  run(sprintf("DDM-7: sw=.10 st0=.05 (2-D) prec=%g", prec), 3, c(0.099, 0.101), c(0.049, 0.051), prec = prec)
}
