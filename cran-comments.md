# cogmod 0.3.0

This is a resubmission. Thank you for the review - all four points are
addressed below.

## Response to the previous review

* **References in `DESCRIPTION`.** Added, in the requested
  `authors (year) <doi:...>` form with no space after `doi:`, for the drift
  diffusion model (Ratcliff and McKoon, 2008), the linear ballistic accumulator
  (Brown and Heathcote, 2008), the lognormal race (Rouder et al., 2015), the
  racing diffusion model (Tillman et al., 2020), the ordered beta model
  (Kubinec, 2023) and the discrete beta model (Sciandra et al., 2024).

* **Missing `\value` tags.** Every exported function now has one. Each help page
  documents a whole model family, so the `\value` section says what each of its
  functions returns and what the output means: the class of the family object
  and of the `stanvars`, the length and meaning of the density and
  random-generation vectors, the two-column data frame the choice-and-RT
  generators return, and the shape of the `log_lik()`, `posterior_predict()` and
  `posterior_epred()` output. Where `posterior_epred()` errors by design -
  because the decision time of that family has no finite mean - the `\value`
  section says so.

  `cogmod-deprecated.Rd` is gone rather than fixed. It documented 140 aliases
  kept for names used before a rename in the previous development version;
  since this is the first release there is nothing to be compatible with, so
  the aliases were removed and the package now exports one name per function.

  `badlm.Rd` is a data set and documents its contents in `\format`.

* **Commented-out example code.** Removed everywhere. The plots now run, the
  `brms::bf()` formulas are built, and `p_outlier()`, `with_outliers()` and
  `cogmod_inits()` - whose examples were entirely commented out - have real,
  executable ones.

* **`\dontrun{}`.** Replaced with `\donttest{}` throughout; there is no
  `\dontrun{}` left in the package. The examples inside it need `cmdstanr` and a
  CmdStan installation, so each is additionally guarded by
  `requireNamespace("cmdstanr")` and `cmdstan_version(error_on_NA = FALSE)`, and
  is skipped rather than failing where the toolchain is absent.

## R CMD check results

0 errors | 0 warnings | 1 note

* checking CRAN incoming feasibility ... NOTE

  New submission.

  Suggests or Enhances not in mainstream repositories: cmdstanr

  `cmdstanr` provides the Stan backend used to fit the models in this package.
  It is not on CRAN and is distributed from the Stan developers' repository,
  which is declared in `Additional_repositories` and which the check confirms
  resolves it. Every use of it, in code, examples and tests, is conditional on
  `requireNamespace()`, `insight::check_if_installed()` or
  `testthat::skip_if_not_installed()`, so the package installs and checks
  normally without it.

## Notes for the reviewer

* Model fitting requires compiling Stan programs, which is far too slow for an
  unconditional example. The three examples that fit a model do so on 200
  simulated observations with one short chain, inside `\donttest{}` and behind
  the `cmdstanr` guard described above.
* The longer, end-to-end worked analyses are published on the package website
  rather than shipped as vignettes, since they depend on pre-fitted models.
