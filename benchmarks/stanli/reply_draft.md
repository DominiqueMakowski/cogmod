**0.21.0 on the #455 models: faster end to end on four of six**

Thanks for turning all of this around so fast. We reran the #455 benchmark on 0.21.0 the day it came out, with the same programs and design (10 participants x 2 conditions x 50 trials, 4 chains x 500 + 500, 8 seeds per model, 4 whole cores per task). Ratios are stanli / cmdstanr within seed, median [min, max], with 0.20.0's in parentheses:

| model | per gradient | per leapfrog, in the fit | fit wall time | ESS per wall second |
| --- | ---: | ---: | ---: | ---: |
| LNR | 0.65 [0.64, 0.66] (0.86) | 0.70 [0.69, 0.75] (1.04) | **0.68** [0.63, 0.86] (1.06) | **1.45** [1.11, 1.89] (1.04) |
| LNR-6 | 0.69 [0.69, 0.70] (0.81) | 0.71 [0.68, 0.75] (0.97) | **0.70** [0.68, 0.83] (1.01) | **1.38** [1.06, 1.52] (1.17) |
| DDM-4 | 1.22 [0.98, 1.27] (1.79) | 1.15 [0.97, 1.22] (1.56) | 1.11 [0.87, 1.79] (1.57) | 0.93 [0.42, 1.20] (0.65) |
| DDM-5 | 0.91 [0.89, 0.92] (1.04) | 0.95 [0.75, 1.28] (1.10) | 0.97 [0.67, 1.62] (1.17) | 1.05 [0.69, 1.62] (0.73) |
| DDM-7 | 0.97 [0.95, 0.99] | | | |
| LBA | 0.72 [0.70, 0.73] (0.86) | 0.74 [0.59, 0.76] (1.16) | **0.75** [0.50, 0.85] (1.32) | **1.65** [1.30, 2.08] (0.97) |
| RDM | 0.56 [0.52, 0.57] (0.65) | 0.57 [0.47, 0.71] (0.85) | **0.53** [0.46, 0.81] (0.79) | **1.43** [0.99, 2.00] (1.06) |

The wall times leave out cmdstanr's 37-46 s of compiling, against stanli's 1 s or so. The DDM-7 row is per gradient only (2 seeds, at 40 starting points): one cmdstanr fit of it ran 25 h on 0.20.0 without finishing.

- **stanli is now the faster backend end to end on the four race models**: 0.53-0.75x cmdstanr's wall time and 1.38-1.65x its ESS per wall second. The DDM-5 is even. The DDM-4 is still slower, but less so than before (see below).
- **The multi-chain slowdown is gone**, as your tests predicted. In the 4-chain fits a leapfrog costs stanli 1.02-1.05x its single-threaded gradient (1.28x on the DDM-5), against 1.04-1.51x on 0.20.0. CmdStan's is 0.94-1.19x. In `chains_probe.R`, four chains cost stanli 1.01-1.17x one chain alone, and CmdStan 1.01-1.09x. On 0.20.0, on the same CPU, stanli's was 1.9x on the RDM and 2.5x on the LBA.
- **Still exact.** The log density is within 27 ULP of CmdStan's at the 1920 posterior draws (the RDM's 53000 ULP are at a log density of 0.02, so 1e-13 absolute), and 7 ULP at the DDM-7's 80 starting points. Gradients agree to 5e-12. The posteriors agree (largest difference in a mean: 3.5 joint MCSEs), and the total leapfrogs per fit are 0.96-1.03x.
- **The build confound we raised is small.** We rebuilt CmdStan's model methods at CmdStan's own flags (`-O3`, no `-march=native`, through a user Makevars; checked on the compiler lines). CmdStan's gradient moves by -4% to +7%. The ratios move by at most 0.02, except the DDM-5's and DDM-7's (0.06). The table uses the matched build. So the build does not explain why stanli's lead on the population-level LNR is larger on Windows (0.42x) than on Linux (0.75x), and our guess in the last comment was wrong.
- **Not the node mix.** 45 of the 48 tasks ran on EPYC 7513s. On the 7513s alone, the wall ratios go from 1.14 / 0.99 / 1.76 / 1.23 / 1.22 / 0.92 on 0.20.0 to 0.68 / 0.69 / 1.11 / 0.97 / 0.74 / 0.53 (LNR, LNR-6, DDM-4, DDM-5, LBA, RDM).

**The DDM-4.** Part of what is left was on our side. About a fifth of its gradient was a routing estimate in cogmod's Stan code: its value only picks a branch, but it was recorded on the tape (the "compute values that only steer a branch without the tape" item of the checklist). We now run a cheap test first that proves the estimate's answer for nearly every trial, so the estimate rarely runs, and every log density and gradient is the same to the last bit ([details](https://github.com/DominiqueMakowski/cogmod/blob/dev/benchmarks/ddm_router/README.md)). On a Windows laptop (CmdStan 2.38.0, stanli 0.21.0, this DDM-4 program at 40 posterior draws), that takes a gradient to 0.77x in CmdStan and 0.83x in stanli. So it helps CmdStan slightly more, and it will not close the 1.22x on its own. The table above is from before this change, and we have not rerun it on the cluster yet. On that laptop, stanli is already the faster engine on this program: 0.76-0.79x CmdStan's gradient before the change, 0.82-0.85x after. So the DDM-4 gap looks Linux-specific, like the population-level LNR's lead (0.42x on Windows, 0.75x on Linux). The build flags are not the reason (see above), and we do not know what is.

**One thing we are watching.** On the RDM, the lowest-ESS parameter (always a participant SD) gets 0.85x cmdstanr's ESS per leapfrog in stanli. It is below 1 in 12 of 16 seeds over both releases. Over all population-level parameters it is 0.94x, and the LBA goes the other way (1.22x). We think 16 seeds are too few to call it, and mention it in case it rings a bell.

**The checklist, as we see it now:**

- [x] What a leapfrog costs inside `sample_model()` beyond `log_prob_grad()`: gone (above).
- [x] Batched gathers for random effects: #456; the ratios above show it.
- [x] A dense metric and a user-supplied `inv_metric`: #458. Not tried yet; it makes our `cogmod_warmstart()` usable with stanli, which is next.
- [x] Values that only steer a branch: fixed on our side for the DDM. It may still matter for other programs that branch on a computed value.
- [ ] Within-chain parallelism: #464, not tested by us. Four chains on four cores leave no core for threads within a chain.
- [ ] Memory a gradient touches: your loop fix may have covered the growth with N, but we have not rerun the N = 10000 probe.
- [ ] How the kernels are built: the 4-argument `wiener_lpdf()` is still 1.22x on Linux, and stanli is faster than CmdStan on the same program on Windows.
- [ ] Hoisting loop-invariant work: brms#1946 is still under review on brms's side.
- [ ] Untested: an x86-64-v3 runtime, SIMD libm under `fast_math`, and the AVX-512 nodes (3 tasks here).

Everything (numbers, harness, commands) is in the [README section "0.21.0"](https://github.com/DominiqueMakowski/cogmod/blob/dev/benchmarks/stanli/README.md#0210-multilevel-models-faster-end-to-end).

From our side, what #455 reported is answered. Close it whenever suits you, or keep it open for the Linux DDM-4 gap.
