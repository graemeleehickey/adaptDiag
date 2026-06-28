# adaptDiag 0.1.2.9000

* Fixed a bug in `multi_trial()` where the Windows `foreach` parallel branch called `single_trial_wrapper()` without passing the iteration argument, causing silent incorrect behaviour on Windows.
* Fixed a bug in `multi_trial()` where specifying `ncores < 1` triggered a warning but did not reset the value to 1, causing a downstream error.
* Fixed a bug in `analysis()` where a simulated dataset with zero observations in any cell of the 2x2 contingency table (possible at early looks with low prevalence) would cause an error. The table is now constructed using `factor()` with explicit levels to guarantee all cells are always present.
* Fixed scalar logical checks throughout `multi_trial()` that used `|` and `&` (vectorised) instead of `||` and `&&` (short-circuit).
* Fixed `1:n` idioms that are unsafe for zero-length inputs; replaced with `seq_len()` and `seq_along()` throughout.
* Fixed all roxygen2 comment prefixes in `R/binom_test_power.R` from `##'` to `#'`; documentation for `binom_sample_size()` was previously never generated from source.
* Fixed a broken `\link{\code{}}` tag in internal documentation for `analysis()`.
* Fixed several typos in `multi_trial()` documentation: `\code{code = "both"}` corrected to `\code{endpoint = "both"}`, and doubled words ("the the", "for for") removed.
* Fixed a duplicate `\item` bullet in `evaluate_trial()` internal documentation; the "insufficient positive cases" outcome now has the correct label.
* Fixed the test description in `test-both.R` which incorrectly described a sensitivity-only test as "both test works".
* Removed dead code in `simulate_data()` that handled a `NULL` `n_at_looks` argument that could never occur.
* Added a pkgdown site with Bootstrap 5, structured reference index, articles, and a GitHub Actions deploy workflow targeting GitHub Pages (#).
* Updated GitHub Actions workflows to latest standards: `actions/checkout@v6`, `codecov/codecov-action@v7` (SHA-pinned), `actions/upload-artifact@v7`, and removed deprecated `use-public-rspm` option.
* Added `library(adaptDiag)` setup chunks to both vignettes, which were missing and caused failures when vignettes were built in a clean subprocess.
* Fixed `%\VignetteIndexEntry{}` titles in both vignettes to match their YAML titles.
* Improved test coverage from 80% to 95%, adding tests for `binom_sample_size()`, `multi_trial()` input validation, and `summarise_trials()` decision branches.

# adaptDiag 0.1.1

* Updated GitHub actions workflows
* Updated README badges

# adaptDiag 0.1.0

* First package release to CRAN.
