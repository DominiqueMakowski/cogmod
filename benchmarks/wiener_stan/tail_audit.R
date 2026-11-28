# At the program point where the new density and the base differ most (by
# program_ab.R's check), which side is right? Prints every DDM-4/5 call's
# inputs and value from both programs (stanli, no compile) and scores each
# trial against grid.R's reference.
#
#   Rscript benchmarks/wiener_stan/tail_audit.R --dir <program_ab OUT>/ddm5 --model ddm5 --variant unrolled
source("benchmarks/gradient_programs.R")  # gp_args()
opt <- gp_args(list(dir = "", model = "ddm5", variant = "unrolled"))
suppressPackageStartupMessages(library(stanli))
chk <- utils::read.csv(file.path(opt$dir, "program_check.csv"))
chk <- chk[chk$engine == "cmdstan" & chk$variant == opt$variant & is.finite(chk$lp_absdiff), ]
worst <- chk$point[which.max(chk$lp_absdiff)]
cat("point", worst, "with |dlp|", max(chk$lp_absdiff), "\n")

sdat <- jsonlite::read_json(file.path(opt$dir, "data.json"), simplifyVector = TRUE)
U <- readRDS(file.path(opt$dir, "draws_u.rds"))
D <- lapply(round(seq(1, nrow(U), length.out = 40)), function(i) as.numeric(U[i, ]))
set.seed(7)  # as program_ab.R
tails <- c(lapply(D, function(q) q + rnorm(length(q), 0, 1.5)), lapply(D, function(q) q + rnorm(length(q), 0, 3)))
pts <- c(setNames(D, paste0("draw", seq_along(D))), setNames(tails, paste0("tail", seq_along(tails))))
q <- pts[[worst]]

# Each program gets a print() of the call's inputs and value. stanli prints
# from native code, past R's connections, so each runs in a child R process
# whose output is read back.
child <- function(v, old) {
  f <- gsub("\\", "/", tempfile(fileext = ".R"), fixed = TRUE)
  saveRDS(list(dir = opt$dir, model = opt$model, v = v, old = old, q = q, sdat = sdat), sub(".R$", ".rds", f))
  writeLines(c(
    "suppressPackageStartupMessages(library(stanli))",
    sprintf("a <- readRDS('%s')", sub(".R$", ".rds", f)),
    "code <- paste(readLines(file.path(a$dir, paste0(a$model, '_', a$v, '.stan'))), collapse = '\\n')",
    "stopifnot(grepl(a$old, code, fixed = TRUE))",
    "ret <- sub(';$', '', sub('^ *return ', '', a$old))",
    "code <- sub(a$old, paste0('    real rtrace = ', ret, ';\\n    print(\"TRIAL \", t, \" \", v, \" \", boundary, \" \", w, \" \", sigmadrift, \" \", rtrace);\\n    return rtrace;'), code, fixed = TRUE)",
    "sm <- stanli_model(code = code, data = a$sdat, threads_per_chain = 1)",
    "invisible(log_prob_grad(sm, a$q))"), f)
  out <- system2(file.path(R.home("bin"), "Rscript"), f, stdout = TRUE, stderr = TRUE)
  x <- grep("^TRIAL", out, value = TRUE)
  if (!length(x)) stop("no trace from ", v, ":\n", paste(utils::tail(out, 20), collapse = "\n"))
  m <- do.call(rbind, lapply(strsplit(sub("^TRIAL ", "", x), " "), as.numeric))
  stats::setNames(as.data.frame(m), c("t", "v", "a", "w", "sv", "ll"))
}
base_call <- "    return wiener_lpdf(y | boundary, tau0, w, v, sigmadrift);"
new_call <- "    return cogmod_wfpt_lpdf(t | boundary, w, v, sigmadrift);"
b <- child("base", base_call)
n <- child(opt$variant, new_call)
stopifnot(nrow(b) == nrow(n))

# grid.R's reference (Stan's convention: upper boundary, w from the lower one)
ref_small <- function(u, wl, w, m = 30) {
  r <- 2 / u
  s <- -expm1(-w * r) - w * (1 + exp(-w * r))
  for (k in 1:m) s <- s + (wl + 2 * k) * exp(-k * (k + wl) * r)
  for (j in 2:m) s <- s + (wl - 2 * j) * exp(-j * (j - wl) * r)
  log(s) - wl^2 / (2 * u) - 0.5 * log(2 * pi) - 1.5 * log(u)
}
ref_large <- function(u, wl, w, K = 200) {
  flip <- wl > 0.5; th <- pi * ifelse(flip, w, wl); s <- 0
  for (k in 1:K) s <- s + ifelse(flip & k %% 2 == 0, -1, 1) * k * exp(-(k^2 - 1) * pi^2 * u / 2) * sin(k * th)
  log(pi) + log(s) - pi^2 * u / 2
}
ref_ll <- function(t, v, a, w, sv) {
  wl <- 1 - w; u <- t / a^2; out <- numeric(length(u)); sm <- u < 0.05
  if (any(sm)) out[sm] <- ref_small(u[sm], wl[sm], w[sm])
  if (any(!sm)) out[!sm] <- ref_large(u[!sm], wl[!sm], w[!sm])
  out - 2 * log(a) + (sv^2 * a^2 * wl^2 + 2 * v * a * wl - v^2 * t) / (2 * (1 + sv^2 * t)) - 0.5 * log1p(sv^2 * t)
}
# The printed inputs carry 6 significant digits, so the reference is taken at
# the printed point and each side's value is compared with it only to the
# precision that rounding allows (reported as the reference's own spread).
r <- ref_ll(b$t, b$v, b$a, b$w, b$sv)
d <- data.frame(b[, 1:5], ref = r, base = b$ll, new = n$ll, err_base = b$ll - r, err_new = n$ll - r)
cat(sprintf("%d trials; |base - ref| max %.3g (sum %.3g), |new - ref| max %.3g (sum %.3g)\n",
            nrow(d), max(abs(d$err_base)), sum(abs(d$err_base)), max(abs(d$err_new)), sum(abs(d$err_new))))
print(utils::head(d[order(-abs(d$err_base)), ], 8), row.names = FALSE, digits = 6)
