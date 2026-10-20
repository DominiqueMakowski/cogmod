# Emit Stan program + data for models whose intercept-only dpars the patch
# turns into reals, from either the installed brms ("base") or a patched
# checkout ("pr"), into benchmarks/results/stanli/brms_scalar_dpar/<which>/.
# compile_check.R then compiles and compares the pairs.
#
#   Rscript benchmarks/stanli/brms_scalar_dpar/emit.R base
#   BRMS_CLONE=<patched checkout> Rscript benchmarks/stanli/brms_scalar_dpar/emit.R pr
which <- match.arg(commandArgs(TRUE)[1], c("base", "pr"))
if (which == "pr") {
  clone <- Sys.getenv("BRMS_CLONE")
  stopifnot(nzchar(clone), dir.exists(clone))
  suppressPackageStartupMessages(pkgload::load_all(clone, quiet = TRUE))
} else {
  suppressPackageStartupMessages(library(brms))
}
out <- file.path("benchmarks/results/stanli/brms_scalar_dpar", which)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

set.seed(1)
N <- 60
d <- data.frame(x = rnorm(N), se = runif(N, 0.1, 0.5), cens = sample(c(-1, 0, 1), N, TRUE),
                dec = sample(0:1, N, TRUE))
d$y <- 1 + 0.5 * d$x + rnorm(N, 0, 0.8)
d$yp <- rpois(N, 3)
d$yo <- factor(sample(1:3, N, TRUE), ordered = TRUE)
ys <- matrix(runif(3 * N), N, 3); d$ys <- ys / rowSums(ys)
d$yb <- d$ys[, 1]
d$rt <- abs(d$y) + 0.2
d$yn <- rnbinom(N, mu = 3, size = 2)

toy <- custom_family("my_lnr", dpars = c("mu", "sigma"), links = c("identity", "log"),
                     lb = c(NA, 0), vars = "vint1[n]", loop = TRUE)
toy_sv <- stanvar(scode = "
  real my_lnr_lpdf(real y, real mu, real sigma, int vint1) {
    return lognormal_lpdf(y | mu + 0.5 * vint1, sigma);
  }", block = "functions")

models <- list(
  gauss_threads_cens = list(bf(y | cens(cens) ~ x, sigma ~ 1), threads = threading(2)),
  gauss_se           = list(bf(y | se(se, sigma = TRUE) ~ x, sigma ~ 1)),
  cumulative_disc    = list(bf(yo ~ x, disc ~ 1), family = cumulative()),
  dirichlet_phi      = list(bf(ys ~ x, phi ~ 1), family = dirichlet()),
  beta_phi           = list(bf(yb ~ x, phi ~ 1), family = Beta()),
  hurdle_lognormal   = list(bf(rt ~ x, hu ~ 1, sigma ~ 1), family = hurdle_lognormal()),
  shifted_lognormal  = list(bf(rt ~ x, ndt ~ 1), family = shifted_lognormal()),
  negbinomial_shape  = list(bf(yn ~ x, shape ~ 1), family = negbinomial()),
  custom_loop        = list(bf(rt | vint(dec) ~ x, sigma ~ 1), family = toy, stanvars = toy_sv),
  gamma_shape_cens   = list(bf(rt | cens(cens) ~ x, shape ~ 1), family = Gamma()),
  negbinomial_rate   = list(bf(yn | rate(yn + 1) ~ x, shape ~ 1), family = negbinomial())
)
for (nm in names(models)) {
  m <- models[[nm]]
  code <- do.call(make_stancode, c(list(m[[1]], data = d), m[-1]))
  writeLines(as.character(code), file.path(out, paste0(nm, ".stan")))
  sd <- do.call(make_standata, c(list(m[[1]], data = d), m[-1]))
  cmdstanr::write_stan_json(lapply(unclass(sd), identity), file.path(out, paste0(nm, ".json")))
  # the count includes lprior (and ptarget under threads); the pr side has
  # one more per intercept-only dpar
  cat(nm, ": ", sum(grepl("real [a-z_]+ = 0;", strsplit(as.character(code), "\n")[[1]])), " 'real x = 0;'\n", sep = "")
}
