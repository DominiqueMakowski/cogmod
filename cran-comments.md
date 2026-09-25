# cogmod 0.3.3

This is an update of cogmod 0.3.0 (published 2026-09-12). Versions 0.3.1 and
0.3.2 were development versions and were never submitted; their changes are
included here and listed in `NEWS.md`. The main additions are warm starts for
repeated fits (`cogmod_warmstart()`), default priors and initial values for
the bounded-scale rating families, and a numerically exact gradient for the
racing diffusion model.

## Test environments

* Local: Windows 11, R 4.5.3, `R CMD check --as-cran` with every suggested
  package and CmdStan 2.38.0 installed, so the `\donttest{}` examples that fit
  a model were run as well.
* GitHub Actions: ubuntu-latest (R devel, release and oldrel-1), macOS-latest
  (release) and windows-latest (release). The ubuntu release job also installs
  CmdStan and runs the examples with `--run-donttest`.

## R CMD check results

0 errors | 0 warnings | 0 notes

The CRAN incoming check reports, as for 0.3.0:

    Suggests or Enhances not in mainstream repositories: cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://mc-stan.org/r-packages/

## Reverse dependencies

There are no reverse dependencies on CRAN.

## Notes for the reviewer

* `cmdstanr`, in Suggests, is not on CRAN; it is distributed from the Stan
  developers' repository, declared in `Additional_repositories`. As in 0.3.0,
  every use of it in code, examples and tests is conditional, and the package
  installs and checks normally without it.

* **`\dontrun{}` in six examples.** In response to the review of 0.3.0, every
  `\dontrun{}` was replaced with `\donttest{}`, and all examples that fit a
  model still use `\donttest{}`. Six examples are back in `\dontrun{}`, each
  with a comment saying why: those of the `*_lpdf_expose()` / `*_lpmf_expose()`
  functions in `rcogmod_betadiscrete`, `rcogmod_betagate`, `rcogmod_ddm`,
  `rcogmod_lba2`, `rcogmod_lnr` and `rcogmod_rdm`. These functions compile a
  Stan function and load it into the R session. R CMD check runs every example
  in one R session, and the `\donttest{}` examples of `cogmod_inits()` and
  `p_outlier()` fit a model there with 'brms', which loads 'rstan'. Loading an
  exposed Stan function after that crashes R with a segmentation fault on
  Linux, so these examples cannot be run under `--run-donttest` wherever
  CmdStan is installed. They need CmdStan in any case, so they would not run
  on CRAN's machines either. The interface is exercised by the package tests
  wherever CmdStan is available.
