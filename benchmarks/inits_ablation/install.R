# Install this tree's cogmod (the source tarball run.sh pushed) into the
# ablation's own library, ahead of the production one on R_LIBS. Run on a
# compute node by `run.sh install`.
lib <- Sys.getenv("ABL_LIB")
stopifnot(nzchar(lib))
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
cat("library:", lib, "\nnode   :", Sys.info()[["nodename"]], "\n")

tarball <- Sys.glob("cogmod_*.tar.gz")
stopifnot(length(tarball) == 1)
cat("installing", tarball, "\n")
install.packages(tarball, repos = NULL, type = "source", lib = lib)

cat("=== CHECK ===\n")
for (p in c("brms", "cmdstanr", "posterior", "mgcv", "cogmod")) {
  cat(sprintf("%-10s %s  %s\n", p,
              tryCatch(as.character(packageVersion(p)), error = function(e) "MISSING"),
              tryCatch(find.package(p), error = function(e) "")))
}
ns <- asNamespace("cogmod")
cat("data-aware layer present:", exists(".data_start", envir = ns), "\n")
if (!exists(".data_start", envir = ns)) stop("wrong cogmod installed")
cat("cmdstan:", tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) "MISSING"), "\n")
