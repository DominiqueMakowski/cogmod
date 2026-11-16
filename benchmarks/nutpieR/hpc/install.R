# Install brms 2.23.1, the pushed stanli tarball and nutpieR into the
# experiment's own library, ahead of production on R_LIBS. Run on a compute
# node by install.sh, which has unpacked stanli's runtime and set
# STANLI_RUNTIME (stanli_install() is never called; run.sh, cmd_push, says
# why). Then check that every piece bench.R loads is the one it should be,
# and that the engines work: a toy model through stanli's NUTS and through
# nutpie, and the pushed tree through load_all().
lib <- Sys.getenv("NP_LIB")
stopifnot(nzchar(lib), file.exists(Sys.getenv("STANLI_RUNTIME")))
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))
cat("library:", lib, "\nnode   :", Sys.info()[["nodename"]], "\n")

# brms at the laptop's commit, so both platforms sample the same program
# (CRAN has 2.23.0, the R module 2.21.0). It wants loo >= 2.8.0 and
# posterior >= 1.6.0, newer than the module's.
have <- function(p) tryCatch(as.character(packageVersion(p, lib.loc = lib)), error = function(e) "")
if (have("brms") != "2.23.1") {
  for (p in c("loo", "posterior")) if (!nzchar(have(p))) install.packages(p, lib = lib)
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", lib = lib)
  remotes::install_github("paul-buerkner/brms@e71e9d743e382db0bb4d8f741a8c837da992946c",
                          upgrade = "never", lib = lib)
}
tarball <- Sys.glob("stanli_*.tar.gz")
stopifnot(length(tarball) == 1)
# The R bridge compiles with the fixed -march in R_MAKEVARS_USER (run.sh,
# cmd_install); the runtime is the release's prebuilt baseline x86-64.
cat("installing", tarball, "with R_MAKEVARS_USER =", Sys.getenv("R_MAKEVARS_USER", "(unset)"), "\n")
install.packages(tarball, repos = NULL, type = "source", lib = lib)
# From source here (r-universe has binaries for R >= 4.5 only): cargo
# fetches the crates and builds nuts-rs and BridgeStan's client, minutes.
if (!nzchar(have("nutpieR"))) {
  install.packages("nutpieR", lib = lib,
                   repos = c("https://andytimm.r-universe.dev", "https://cloud.r-project.org"))
}

cat("=== CHECK ===\n")
for (p in c("brms", "stanli", "nutpieR", "cmdstanr", "rstan", "posterior", "loo", "pkgload")) {
  cat(sprintf("%-10s %-10s %s\n", p,
              tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING"),
              tryCatch(find.package(p), error = function(e) "")))
}
stopifnot(packageVersion("stanli") >= "0.26.0", packageVersion("nutpieR") >= "1.9.0",
          as.character(packageVersion("brms")) == "2.23.1")
cat("cmdstan:", tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) "MISSING"), "\n")
cat("runtime:", stanli::stanli_runtime_path(), "\n")
stopifnot(normalizePath(stanli::stanli_runtime_path()) == normalizePath(Sys.getenv("STANLI_RUNTIME")))
m <- stanli::stanli_model(code = "parameters { real x; } model { x ~ normal(0, 1); }")
f1 <- stanli::sample_model(m, chains = 2, warmup = 200, samples = 200, refresh = 0)
f2 <- stanli::sample_nutpie(m, chains = 2, warmup = 200, samples = 200)
sd_x <- function(fit) { s <- summary(fit); round(s$sd[s$variable == "x"], 2) }
cat("toy model, sd of x (1): stanli", sd_x(f1), "- nutpie", sd_x(f2), "\n")
tree <- file.path(Sys.getenv("NP_DIR"), "tree")
suppressPackageStartupMessages(pkgload::load_all(tree, quiet = TRUE, export_all = TRUE))
cat("cogmod", as.character(packageVersion("cogmod")), "from", tree, "- brms in use",
    as.character(packageVersion("brms")), "\n")
cat("INSTALL OK\n")
