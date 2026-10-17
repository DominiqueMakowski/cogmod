# Install this tree's cogmod (the source tarball run.sh pushed) into the
# screen's own library, ahead of the production one on R_LIBS. Run on a
# compute node by `run.sh install`. Then check that each arm will find the
# cogmod it asks for: this tree's with the library on the path, the
# production one (the 0.3.3 defaults) without it.
lib <- Sys.getenv("SCREEN_LIB")
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
probe <- function(paths) {
  cmd <- sprintf(
    ".libPaths(c(%s)); ns <- asNamespace('cogmod'); cat(as.character(packageVersion('cogmod')), exists('.zs_scales', envir = ns), find.package('cogmod'))",
    paste0("'", paths, "'", collapse = ", ")
  )
  system2(file.path(R.home("bin"), "Rscript"), c("-e", shQuote(cmd)), stdout = TRUE)
}
prod <- setdiff(.libPaths(), normalizePath(lib))
cat("0.3.4 arms see:", probe(.libPaths()), "\n")
cat("0.3.3 arm sees:", probe(prod), "\n")
cat("cmdstan:", tryCatch(as.character(cmdstanr::cmdstan_version()), error = function(e) "MISSING"), "\n")
