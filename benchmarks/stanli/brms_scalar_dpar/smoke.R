# Smoke test for the scalar intercept-only dpar patch (brms_scalar_dpar.patch):
# load a brms checkout with the patch applied and, for each formula, print the
# lines of the generated model / partial_log_lik block that declare or use the
# dpars. "real X = 0;" is the scalar form, "vector[N] X = rep_vector" brms's.
#
#   BRMS_CLONE=<path to the patched brms checkout> Rscript benchmarks/stanli/brms_scalar_dpar/smoke.R
clone <- Sys.getenv("BRMS_CLONE")
stopifnot(nzchar(clone), dir.exists(clone))
suppressPackageStartupMessages(pkgload::load_all(clone, quiet = TRUE))
set.seed(1)
N <- 60
d <- data.frame(x = rnorm(N), g = rep(1:6, 10), se = runif(N, 0.1, 0.5),
                cens = sample(0:1, N, TRUE), w = runif(N, 0.5, 2),
                t = rep(1:10, each = 6), tr = sample(5:10, N, TRUE),
                dec = sample(0:1, N, TRUE))
d$y <- 1 + 0.5 * d$x + rnorm(N, 0, 0.8)
d$yp <- rpois(N, 3)
d$yo <- factor(sample(1:3, N, TRUE), ordered = TRUE)
d$ys <- as.matrix(data.frame(a = runif(N), b = runif(N), c = runif(N)))
d$ys <- d$ys / rowSums(d$ys)
d$y2 <- 2 + d$y + rnorm(N)
d$rt <- abs(d$y) + 0.2
d$yn <- rnbinom(N, mu = 3, size = 2)

show <- function(label, f, ...) {
  code <- tryCatch(as.character(make_stancode(f, data = d, ...)), error = function(e) {
    cat("\n#### ", label, ": ERROR ", conditionMessage(e), "\n"); return(NULL) })
  if (is.null(code)) return(invisible())
  lines <- strsplit(code, "\n")[[1]]
  keep <- grepl("real [a-z_0-9]+ = 0;|rep_vector\\(0\\.0|target \\+=|partial_log_lik|\\[n\\] = |\\[n\\] \\+=", lines) &
    !grepl("lprior|std_normal|reduce_sum\\(", lines)
  cat("\n#### ", label, "\n", paste(trimws(lines[keep]), collapse = "\n"), "\n")
}

cat("\n==== expected SCALAR ====\n")
show("gaussian sigma~1", bf(y ~ x, sigma ~ 1))
show("gaussian sigma~1, threads", bf(y ~ x, sigma ~ 1), threads = threading(2))
show("gaussian sigma~1, cens", bf(y | cens(cens) ~ x, sigma ~ 1))
show("gaussian sigma~1, cens, threads", bf(y | cens(cens) ~ x, sigma ~ 1), threads = threading(2))
show("gaussian sigma~1, weights", bf(y | weights(w) ~ x, sigma ~ 1))
show("gaussian sigma~1, trunc", bf(y | trunc(lb = -10) ~ x, sigma ~ 1))
show("gaussian sigma~1, se", bf(y | se(se, sigma = TRUE) ~ x, sigma ~ 1))
show("student nu~1 sigma~1", bf(y ~ x, sigma ~ 1, nu ~ 1), family = student())
show("student sigma~1 (nu scalar par)", bf(y ~ x, sigma ~ 1), family = student())
show("cumulative disc~1", bf(yo ~ x, disc ~ 1), family = cumulative())
show("dirichlet phi~1", bf(ys ~ x, phi ~ 1), family = dirichlet())
show("hurdle_lognormal hu~1 sigma~1", bf(rt ~ x, hu ~ 1, sigma ~ 1), family = hurdle_lognormal())
show("shifted_lognormal ndt~1", bf(rt ~ x, ndt ~ 1), family = shifted_lognormal())
show("zero_inflated_poisson zi~1", bf(yp ~ x, zi ~ 1), family = zero_inflated_poisson())
show("negbinomial shape~1", bf(yn ~ x, shape ~ 1), family = negbinomial())
show("negbinomial shape~1 + rate", bf(yn | rate(tr) ~ x, shape ~ 1), family = negbinomial())
show("inverse.gaussian shape~1 (stays vectorised)", bf(rt ~ x, shape ~ 1), family = inverse.gaussian())
show("beta phi~1 (y in (0,1))", bf(ys[,1] ~ x, phi ~ 1), family = Beta())
lnr <- custom_family("my_lnr", dpars = c("mu", "sigma"), links = c("identity", "log"),
                     lb = c(NA, 0), vars = "vint1[n]", loop = TRUE)
show("custom loop=TRUE sigma~1", bf(rt | vint(dec) ~ x, sigma ~ 1), family = lnr)
show("gaussian sigma~1 + mu~1 (mu stays vector)", bf(y ~ 1, sigma ~ 1))
show("mv no rescor sigma~1", bf(y ~ x, sigma ~ 1) + bf(y2 ~ x, sigma ~ 1) + set_rescor(FALSE))

cat("\n==== expected VECTOR (gated out) ====\n")
show("sigma~1+(1|g)", bf(y ~ x, sigma ~ 1 + (1 | g)))
show("sigma~0+Intercept", bf(y ~ x, sigma ~ 0 + Intercept))
show("skew_normal sigma~1", bf(y ~ x, sigma ~ 1), family = skew_normal())
show("mixture sigma1~1", bf(y ~ x, sigma1 ~ 1), family = mixture(gaussian, gaussian))
lnr_vec <- custom_family("my_lnr_vec", dpars = c("mu", "sigma"), links = c("identity", "log"),
                         lb = c(NA, 0), vars = "vint1", loop = FALSE)
show("custom loop=FALSE sigma~1", bf(rt | vint(dec) ~ x, sigma ~ 1), family = lnr_vec)
show("mv rescor sigma~1", bf(y ~ x, sigma ~ 1) + bf(y2 ~ x, sigma ~ 1) + set_rescor(TRUE))
show("ar() sigma~1", bf(y ~ x + ar(t, g), sigma ~ 1))
show("nl a+b~1 (not centered)", bf(y ~ a * exp(-b * x), a + b ~ 1, nl = TRUE),
     prior = prior(normal(0, 1), nlpar = "a") + prior(normal(0, 1), nlpar = "b"))
