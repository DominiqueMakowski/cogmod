# Where does stanli lose time on the multilevel models? Per-gradient cost,
# CmdStan vs stanli, of program variants that each remove one suspect.
#   DDM-4: orig | noroute (the router call dropped from the branch condition)
#          | direct (straight to the 4-arg wiener_lpdf, no branches at all)
#   LNR:   ml (multilevel, N = 1000) | pop (same data, no random effects)
#          | ml_big (multilevel, 10 x 2 x 500 = 10000 trials)
# 40 points around cogmod_inits() (sd 0.1, unconstrained), 11 alternating
# blocks of ~0.5 s. Needs multilevel.R's ddm4 and lnr inputs, and the big LNR:
#
#   Rscript benchmarks/stanli/multilevel.R emit --model ddm4
#   Rscript benchmarks/stanli/multilevel.R emit --model lnr
#   Rscript benchmarks/stanli/multilevel.R emit --model lnr --trials 500 --out $TEMP/stanli_probe/big
#   Rscript benchmarks/stanli/probe_bottlenecks.R
#
# Windows, 2026-10-12 (results/probe_bottlenecks_win.csv), us per gradient
# CmdStan / stanli: ddm4 orig 2000 / 2071, noroute 1556 / 1750, direct
# 1528 / 1694; lnr ml 882 / 441, pop 824 / 309, ml_big 8250 / 6625.
source("benchmarks/gradient_programs.R")
gp_load(".")
suppressPackageStartupMessages(library(stanli))
out <- file.path(Sys.getenv("TEMP"), "stanli_probe")
dir.create(out, showWarnings = FALSE)
subs1 <- function(code, old, new) {
  stopifnot(lengths(regmatches(code, gregexpr(old, code, fixed = TRUE))) == 1)
  sub(old, new, code, fixed = TRUE)
}
prog <- function(f, df) {
  prior <- suppressMessages(cogmod_priors(f, df)); sv <- cogmod_stanvars(f)
  list(code = as.character(brms::make_stancode(f, data = df, prior = prior, stanvars = sv)),
       sdat = lapply(unclass(brms::make_standata(f, data = df, prior = prior, stanvars = sv)), identity),
       init = cogmod_inits(f, df, jitter = 0)(1))
}
d4 <- readRDS("benchmarks/results/stanli_ml/ddm4_inputs.rds")
p4 <- prog(d4$f, d4$df)
ln <- readRDS("benchmarks/results/stanli_ml/lnr_inputs.rds")
pml <- prog(ln$f, ln$df)
fpop <- brms::bf(RT | dec(Error) ~ Condition, nuone ~ Condition, sigmazero ~ 1, sigmaone ~ 1,
                 sigmabias = 0, ndt ~ 1, family = cogmod_lnr())
ppop <- prog(fpop, ln$df)
big <- readRDS(file.path(out, "big", "lnr_inputs.rds"))
pbig <- prog(big$f, big$df)

route_old <- "if (sigmadrift == 0 && t > 0.01 * square(boundary)\n        && cogmod_ddm_log_density_scale(t, v, boundary, w, 0, 0, 0) > -600) {"
variants <- list(
  ddm4_orig = p4,
  ddm4_noroute = modifyList(p4, list(code = subs1(p4$code, route_old,
    "if (sigmadrift == 0 && t > 0.01 * square(boundary)) {"))),
  ddm4_direct = modifyList(p4, list(code = subs1(p4$code, "  real y = t + tau0;\n",
    "  real y = t + tau0;\n  return wiener_lpdf(y | boundary, tau0, w, v);\n"))),
  lnr_ml = pml, lnr_pop = ppop, lnr_ml_big = pbig)

res <- list()
for (v in names(variants)) {
  x <- variants[[v]]
  sf <- file.path(out, paste0(v, ".stan")); writeLines(x$code, sf)
  dj <- file.path(out, paste0(v, ".json")); cmdstanr::write_stan_json(x$sdat, dj)
  mod <- gp_compile(sf)
  fc <- mod$sample(data = dj, init = list(x$init), chains = 1, iter_warmup = 1, iter_sampling = 1,
                   fixed_param = TRUE, refresh = 0, show_messages = FALSE)
  fc$init_model_methods(verbose = FALSE)
  q0 <- as.numeric(fc$unconstrain_variables(x$init))
  sm <- stanli_model(code = x$code, data = x$sdat, threads_per_chain = 1)
  set.seed(3)
  P <- lapply(1:40, function(i) q0 + rnorm(length(q0), 0, 0.1))
  a <- fc$grad_log_prob(P[[1]]); b <- log_prob_grad(sm, P[[1]])
  cat(sprintf("%-13s lp cmdstan %.6f stanli %.6f\n", v, attr(a, "log_prob"), b$lp))
  timers <- list(cmdstan = function(k) for (i in seq_len(k)) for (q in P) fc$grad_log_prob(q),
                 stanli = function(k) for (i in seq_len(k)) for (q in P) log_prob_grad(sm, q))
  t1 <- system.time(timers$cmdstan(1))[["elapsed"]]; timers$stanli(1)
  passes <- max(1L, as.integer(ceiling(0.5 / max(t1, 0.005))))
  for (r in 1:11) for (e in sample(names(timers))) {
    st <- system.time(timers[[e]](passes))[["elapsed"]]
    res[[length(res) + 1]] <- data.frame(variant = v, rep = r, engine = e,
                                         us = st / (passes * length(P)) * 1e6)
  }
}
res <- do.call(rbind, res)
utils::write.csv(res, file.path(out, "probe.csv"), row.names = FALSE)
w <- stats::reshape(res, idvar = c("variant", "rep"), timevar = "engine", direction = "wide")
tab <- do.call(rbind, lapply(split(w, w$variant), function(x) data.frame(
  variant = x$variant[1], cmdstan_us = round(stats::median(x$us.cmdstan)),
  stanli_us = round(stats::median(x$us.stanli)),
  ratio = round(stats::median(x$us.stanli / x$us.cmdstan), 2))))
print(tab[match(names(variants), tab$variant), ], row.names = FALSE)
