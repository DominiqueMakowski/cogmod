# Inits ablation on Artemis

**Status: the data-aware layer this measured was removed from the tree on
2026-09-22 (see `../data_inits_ablation.md`, "Decision").** The
implementation is kept here as `data_aware_inits.patch`; to rerun the
ablation, apply it first (`git apply benchmarks/inits_ablation/data_aware_inits.patch`
from the repository root), then `./run.sh build` from the patched tree.
`fit_cell.R` refuses to run against a cogmod without the layer, since the two
arms would then be identical. The result rows of the 2026-09-22 run are not
tracked (`results/` is gitignored); `../data_inits_ablation.md` has the tables.

Does the data-aware layer in `cogmod_inits()` (`.data_start()`, the
`data_init` rules on the registry entries) buy a fit anything a constant start
does not? The laptop runs in `../data_inits_ablation.md` said no on 1000-row
models with two seeds; this is the same question with enough seeds to resolve
a 10% effect, harder data, and the metric that a better start could actually
move: the cold-start transient AGENT.md of the Illusion Game project measured
at 100-150 iterations of maximum treedepth before the first metric window.

`cells.R` is the design, `fit_cell.R` fits one cell and writes one CSV row,
`fit.slurm` is the array task, `run.sh` drives it over SSH, `summarise.R`
compares the arms seed by seed. Everything runs on the `general` and `short`
partitions, whose quotas are separate from the `long` / `sussexneuro` ones the
production fits use, into its own directory and its own R library, so nothing
here can touch a production run.

```bash
cd benchmarks/inits_ablation
./run.sh build              # speed_acc.csv + cogmod_<version>.tar.gz from this tree
./run.sh push
./run.sh install            # this tree's cogmod into the ablation's own library
./run.sh submit A --array=1-2 --partition=short   # smoke test
./run.sh submit A
./run.sh submit B
./run.sh queue ; ./run.sh progress ; ./run.sh log
./run.sh pull
Rscript summarise.R results
```

The "constant" arm is `cogmod_inits()` as it was before the data-aware layer,
rebuilt in `fit_cell.R` from the same internals, so one installed cogmod
serves both arms and nothing else differs between them.
