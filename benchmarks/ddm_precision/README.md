# The tolerance passed to Stan's 7-parameter `wiener_lpdf()`

What `.DDM_WIENER_PRECISION` (R/model_ddm.R) cost and bought, measured
2026-09-18 on CmdStan 2.38.0, Windows, 16 logical cores. Results in
`benchmarks/results/ddm_precision/`. These scripts are records of that
measurement, not a CI job: a 7-parameter DDM fit takes hours even on 100
trials, and every script here takes `--base` and `--pr` directories emitted
from two working trees that differ only in that constant.

| file | what |
| --- | --- |
| `wienerbench.stan`, `wienerbench2.R` | The report's standalone benchmark of `wiener_lpdf()`'s forms, with the tolerance as data. Gives the per-observation microseconds by branch and tolerance at one benign operating point |
| `emit.R` | `emit --pkg <tree> --out <dir> [--n]`: the intercept-only brms program, data and start for `cogmod_ddm`, through `benchmarks/gradient_programs.R`. Its `run` mode is a cold-start fit; it was abandoned after 90 minutes without leaving warmup and is kept only because `emit` lives in the same file |
| `ess_warm.R` | Stage A, end to end. Reads a posterior mode off an optimiser CSV (`--mode`), starts every chain there under one hand-set diagonal metric, runs the same short warmup and sampling on both trees and reports leapfrogs, step size, acceptance, divergences, ESS per second and per gradient. With 100 warmup iterations CmdStan adapts the step size only, so the metric stays as given - identical for both trees |
| `ess_fixed.R` | Stage B, adaptation off. Takes the step size, metric and last draws of a finished run (`--from`) and samples both trees from that one state, so every difference is the density and its gradient |

The runs behind the numbers in NEWS and in the comment on the constant:

```
Rscript benchmarks/ddm_precision/emit.R emit --pkg <base tree> --out <b> --n 100
Rscript benchmarks/ddm_precision/emit.R emit --pkg .           --out <p> --n 100
Rscript benchmarks/ddm_precision/ess_warm.R  --base <b> --pr <p> --out <A> --mode <optimiser csv> --warmup 100 --sampling 200 --treedepth 7
Rscript benchmarks/ddm_precision/ess_fixed.R --from <A>/warm1_base --base <b> --pr <p> --out <B> --iter 200 --treedepth 7 --order base,pr
```

Two things learned the hard way, for the next fit-level benchmark of this
family: run with `save_warmup = TRUE` and a `refresh` from the start, or a
slow warmup is indistinguishable from a hung one; and do not take a metric
from Laplace draws - they land in the tails, where the 7-parameter density
costs seconds per evaluation (cubature's 6000-evaluation cap, eight times
over), and 2000 of them did not finish in half an hour.
