# Reruns only the RDM cells of the warm-start ablation benchmark, resuming from
# the results CSV (strip the family's rows from warmstart_ablation.csv and
# warmstart_ablation_postsd.csv first, or nothing happens). Used on 2026-09-12
# after cogmod_inits() changed the RDM's error-drift start and cogmod_priors()
# widened the default ndt prior; the previous rows are in
# results/archive-2026-09-12-rdm-old-init-and-ndt-prior/, and a second pass
# with only the driftone start changed is in
# results/archive-2026-09-13-rdm-driftone-start-only/. The summary and the
# figure it writes at the end cover the RDM alone; run the script once more
# with its defaults afterwards (every cell is then "already complete") to
# redraw them for all four families.
#
#   "C:\Program Files\R\R-4.5.3\bin\Rscript.exe" benchmarks/run_warmstart_ablation_rdm.R
bench_overrides <- list(families = "rdm")
source("benchmarks/warmstart_ablation.R")
