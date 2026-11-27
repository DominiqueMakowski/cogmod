# cogmod 0.3.6

This is an update of cogmod 0.3.3 (published 2026-09-25). Versions 0.3.4 and
0.3.5 were development versions; their changes are included here and listed
in `NEWS.md`. The main changes are a faster and more accurate Stan likelihood
for the drift diffusion model, which no longer depends on the accuracy of
Stan's 5-parameter `wiener_lpdf()` in the tails; a fix for a density of the
linear ballistic accumulator that could be NaN, with exact gradients in its
tails; revised default priors (smooth terms, the outlier proportion, and the
DDM's boundary and bias); and warm starts from fits made with the new
`stanli` backend of 'brms'.

## Test environments

* Local: Windows 10, R 4.5.3, `R CMD check --as-cran` with every suggested
  package and CmdStan 2.38.0 installed, so the `\donttest{}` examples that fit
  a model were run as well.
* GitHub Actions: ubuntu-latest (R devel, release and oldrel-1), macOS-latest
  (release) and windows-latest (release). The ubuntu release job also installs
  CmdStan and runs the examples with `--run-donttest`.

## R CMD check results

0 errors | 0 warnings | 0 notes

(Locally, two NOTEs came from the check machine rather than the package: the
current time could not be verified, and `pandoc` was not on the PATH for
checking `NEWS.md`, which converts cleanly with it.)

The website moved from `dominiquemakowski.github.io/cogmod` to
`dominiquemakowski.com/cogmod`; every URL in the package has been updated to
the new address.

The CRAN incoming check reports, as for previous versions:

    Suggests or Enhances not in mainstream repositories: cmdstanr
    Availability using Additional_repositories specification:
      cmdstanr   yes   https://mc-stan.org/r-packages/

## Reverse dependencies

There are no reverse dependencies on CRAN.

## Notes for the reviewer

* `cmdstanr`, in Suggests, is not on CRAN; it is distributed from the Stan
  developers' repository, declared in `Additional_repositories`. As in
  previous versions, every use of it in code, examples and tests is
  conditional, and the package installs and checks normally without it.

* **`\dontrun{}` in six examples**, unchanged since 0.3.3. All examples that
  fit a model use `\donttest{}`. The six in `\dontrun{}` are those of the
  `*_lpdf_expose()` / `*_lpmf_expose()` functions in `rcogmod_betadiscrete`,
  `rcogmod_betagate`, `rcogmod_ddm`, `rcogmod_lba2`, `rcogmod_lnr` and
  `rcogmod_rdm`, each with a comment saying why. These functions compile a
  Stan function and load it into the R session. R CMD check runs every example
  in one R session, and the `\donttest{}` examples of `cogmod_inits()` and
  `p_outlier()` fit a model there with 'brms', which loads 'rstan'. Loading an
  exposed Stan function after that crashes R with a segmentation fault on
  Linux, so these examples cannot be run under `--run-donttest` wherever
  CmdStan is installed. They need CmdStan in any case, so they would not run
  on CRAN's machines either. The interface is exercised by the package tests
  wherever CmdStan is available.
