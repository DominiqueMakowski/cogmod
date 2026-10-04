# Per-node probe behind the SIGILL diagnosis of 2026-10-04 (README, "The
# cluster run"): which CPU, which stanli runtime, which R bridge, and does
# a toy model and then the brms LNR program build and sample. A SIGILL kills
# the process, so how far the output gets says which step faulted, and R's
# own handler prints the fault address (its page offset identifies the
# instruction in the bridge: `objdump -d stanli.so`).
#
# Run on one node under the harness's environment, for example
#
#   srun --partition=short --nodelist=artemis-general-02 --cpus-per-task=2 \
#        --mem=6G --time=00:12:00 bash -c '. /etc/profile.d/lmod.sh;
#        module use /mnt/shared/easybuild/modules/all;
#        module load CmdStanR/0.7.1-foss-2023a-R-4.3.2;
#        export R_LIBS=<lib with the bridge to test>:<experiment lib>:<production lib>;
#        STANLI_RUNTIME=<a libstanli.so> Rscript sigill_probe.R'
#
# STANLI_RUNTIME (stanli's own variable) picks the runtime .so; leave it
# unset for the one the package pins. STANLI_PROBE_DIR (default
# ~/stanli_probe) holds lnr_orig.stan and lnr.data.json, as bench.R emit
# writes them to benchmarks/results/stanli/; without them the LNR step is
# skipped. To test a bridge built on a given node, install the package
# tarball there into its own library and put that library first on R_LIBS:
#
#   R CMD INSTALL --no-test-load -l <lib> stanli_0.19.1.tar.gz
#   objdump -d --no-show-raw-insn <lib>/stanli/libs/stanli.so | grep -c zmm
cat("node    :", Sys.info()[["nodename"]], "|",
    sub(".*: ", "", grep("model name", readLines("/proc/cpuinfo"), value = TRUE)[1]), "\n")
cat("flags   :", sub(".*: ", "", grep("^flags", readLines("/proc/cpuinfo"), value = TRUE)[1]), "\n")
cat("glibc   :", system("ldd --version | head -1", intern = TRUE), "\n")
cat("STANLI_RUNTIME =", Sys.getenv("STANLI_RUNTIME", "(unset)"), "\n")
cat("package :", as.character(packageVersion("stanli")), "from", find.package("stanli"), "\n")
ok <- tryCatch(stanli:::load_runtime(), error = function(e) paste("load error:", conditionMessage(e)))
cat("load_runtime():", format(ok), "\n")
cat("runtime :", stanli::stanli_runtime_path(), "\n")
flush(stdout())

cat("[toy] building...\n"); flush(stdout())
m <- stanli::stanli_model(code = "parameters { real mu; } model { mu ~ normal(0, 1); }", data = list())
cat("[toy] built; sampling...\n"); flush(stdout())
f <- stanli::sample_model(m, chains = 2, warmup = 50, samples = 50, seed = 1, refresh = 0)
cat("[toy] OK\n"); flush(stdout())

dir <- Sys.getenv("STANLI_PROBE_DIR", "~/stanli_probe")
stan <- file.path(dir, "lnr_orig.stan"); json <- file.path(dir, "lnr.data.json")
if (file.exists(stan) && file.exists(json)) {
  code <- paste(readLines(stan), collapse = "\n")
  data <- jsonlite::fromJSON(json)
  cat("[lnr] N =", data$N, "- building...\n"); flush(stdout())
  m <- stanli::stanli_model(code = code, data = data, threads_per_chain = 1)
  cat("[lnr] built; sampling 2 chains x 20 + 20...\n"); flush(stdout())
  f <- stanli::sample_model(m, chains = 2, warmup = 20, samples = 20, seed = 1, refresh = 0)
  cat("[lnr] OK\n")
} else cat("[lnr] program files not found in", dir, "- skipped\n")
