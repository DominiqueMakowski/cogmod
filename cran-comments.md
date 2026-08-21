# cogmod 0.3.0

This is a new submission.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE

  New submission.

  Suggests or Enhances not in mainstream repositories: cmdstanr

  `cmdstanr` provides the Stan backend used to fit the models in this package.
  It is not on CRAN and is distributed from the Stan developers' repository,
  which is declared in `Additional_repositories` and which the check confirms
  resolves it. Every use of it, in code, examples and tests, is conditional on
  `insight::check_if_installed()` or `testthat::skip_if_not_installed()`, so
  the package installs and checks normally without it.

## Notes for the reviewer

* Model fitting requires compiling Stan programs, which is far too slow for
  examples and tests. The functions that do so are exercised in tests that are
  skipped on CRAN and run in continuous integration instead.
* The longer, end-to-end worked analyses are published on the package website
  rather than shipped as vignettes, since they depend on pre-fitted models.
