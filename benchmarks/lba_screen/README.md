# LBA screen on Artemis

Does the chain split of the Illusion Game's full-data `gam_lba` fit come from
the LBA's |v|/s² ray, or from the start-point range that every split fit
shares? Set up 2026-10-01, after the cross-check of the four full-data
MullerLyer fits: `gam_lba` (5 chains against 3), `gam_rdm5` (one chain 1,671
lp__ units above the other seven) and `gam_lnr6` (one chain frozen, one 1,030
below) all smooth `sigmabias`; `gam_ddm5`, with `sigmabias = 0`, is clean.
The review that motivated it: https://claude.ai/artifact/15vwAHwEJTKZ5Y47C3f6DS

Each arm is 8 cold-start chains on the first 200 participants (about 29,000
rows), warmup 1000 + 500, the production formula, priors and settings. One
array task per chain.

| arm | model | cogmod | question |
| --- | --- | --- | --- |
| `lba_old` | `gam_lba` | production (d04c7f8) | does the split exist at this size? |
| `lba_new` | `gam_lba` | this tree (0.3.4) | do the new defaults change anything? |
| `lba_a` | `sigmaone ~ 1 + (1 \| P)` | 0.3.4 | is it the ray? |
| `lba_b` | `sigmabias ~ 1 + (1 \| P)` | 0.3.4 | is it the start-point range? |
| `rdm5_new` | `gam_rdm5` | 0.3.4 | does the RDM split at this size? |
| `rdm5_b` | `gam_rdm5`, `sigmabias ~ 1 + (1 \| P)` | 0.3.4 | the same hypothesis, with no ray |

```bash
cd benchmarks/lba_screen
./run.sh build        # tarball of this tree
./run.sh push
./run.sh install      # this tree's cogmod into the screen's own library; checks both arms' cogmod
./run.sh submit lba --array=1,9 --partition=short   # smoke test: chain 1 of lba_old and lba_new
./run.sh submit lba   # 32 chains, at most 32 at once (SCREEN_MAX)
./run.sh submit rdm   # 16 chains, optional
./run.sh queue ; ./run.sh progress ; ./run.sh log
./run.sh pull
Rscript summarise.R results
```

`fit_cell.R` checks that each arm loaded the cogmod it asked for, and stops
otherwise. A local smoke test runs from this directory with
`SCREEN_LOAD_ALL=../.. SCREEN_PARTICIPANTS=10 SCREEN_WARMUP=30
SCREEN_SAMPLES=20 SLURM_ARRAY_TASK_ID=9 Rscript fit_cell.R` (the 0.3.3 arm
uses whatever cogmod is installed, which must be 0.3.3).

## Reading the result

A chain more than 100 lp__ units from its arm's median is in another place;
an arm "agrees" when all eight chains are within that, none is frozen (step
size below a hundredth of the median) and the pooled R-hat over the
population-level parameters is below 1.05. `summarise.R` prints the tables and
applies the rules:

- `lba_old` agrees: 200 participants cannot see the problem. Step up to 480
  before spending full-data time.
- `lba_old` and `lba_new` both split: the defaults were not the cause (the
  expected result). Then (b) agrees and (a) splits means the start-point
  range, so refit full data as (b) and screen `gam_rdm5` / `gam_lnr6` the same
  way; (a) agrees and (b) splits means the ray, option (a); both agree, choose
  on elpd at 200; neither, try (a) and (b) together, then `sigmaone = 1`, then
  drop the LBA.
- `lba_new` agrees where `lba_old` splits: the defaults mattered after all.
  Re-run full data unchanged with 0.3.4.

Agreement at 200 participants is necessary, not sufficient, for full data.
