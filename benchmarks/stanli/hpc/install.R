# Install brms 2.23.1 (the local version, from GitHub at e71e9d7: CRAN has
# 2.23.0, the cluster 2.21.0) and the
# pushed stanli tarball into the experiment's own library, ahead of the
# production one on R_LIBS, then stanli's runtime (into ~/.cache, shared by
# the nodes). Run on a compute node by `run.sh install`. Then check that every
# piece the benchmark loads is the one it should be, and that both backends
# work: a toy model in stanli, and the pushed tree through load_all().
lib <- Sys.getenv("STANLI_LIB")
stopifnot(nzchar(lib))
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))
cat("library:", lib, "\nnode   :", Sys.info()[["nodename"]], "\n")

have <- function(p) tryCatch(as.character(packageVersion(p, lib.loc = lib)), error = function(e) "")
# brms at e71e9d7 wants loo >= 2.8.0 and posterior >= 1.6.0; the production
# library has 2.7.0 and 1.5.0. Every other floor in its DESCRIPTION is met.
for (p in c("loo", "posterior")) if (!nzchar(have(p))) install.packages(p, lib = lib)
if (have("brms") != "2.23.1") {
  # Installs into .libPaths()[1]; upgrade = "never" leaves the production
  # dependencies alone unless brms needs one it does not have.
  remotes::install_github("paul-buerkner/brms@e71e9d743e382db0bb4d8f741a8c837da992946c",
                          upgrade = "never")
}
tarball <- Sys.glob("stanli_*.tar.gz")
stopifnot(length(tarball) == 1)
cat("installing", tarball, "\n")
install.packages(tarball, repos = NULL, type = "source", lib = lib)
stanli::stanli_install()

cat("=== CHECK ===\n")
for (p in c("brms", "stanli", "cmdstanr", "rstan", "posterior", "loo", "pkgload", "rtdists")) {
  cat(sprintf("%-10s %-10s %s\n", p,
              tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING"),
              tryCatch(find.package(p), error = function(e) "")))
}
cat("cmdstan:", tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) "MISSING"), "\n")
m <- stanli::stanli_model(code = "parameters { real x; } model { x ~ normal(0, 1); }")
g <- stanli::log_prob_grad(m, 0.5)
cat("stanli toy model: lp", g$lp, "grad", g$grad, "\n")
tree <- file.path(Sys.getenv("STANLI_DIR"), "tree")
suppressPackageStartupMessages(pkgload::load_all(tree, quiet = TRUE, export_all = TRUE))
cat("cogmod", as.character(packageVersion("cogmod")), "from", tree, "- brms in use",
    as.character(packageVersion("brms")), "\n")
