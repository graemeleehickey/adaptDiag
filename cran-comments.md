## Test environments

* local macOS (Tahoe 26.5.1, aarch64), R 4.5.3
* ubuntu-latest (on GitHub Actions), R devel
* ubuntu-latest (on GitHub Actions), R release
* ubuntu-latest (on GitHub Actions), R oldrel-1
* windows-latest (on GitHub Actions), R release
* macos-latest (on GitHub Actions), R release
* win-builder (devel)
* win-builder (release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Changes since last submission (0.1.0)

* `multi_trial()`: added a `seed` argument for reproducible simulations.
  Reproducibility is handled via `doRNG::%dorng%`, which correctly manages
  RNG stream splitting across all parallel backends.
* `multi_trial()`: replaced platform-specific parallel dispatch
  (`pbmcapply` on Unix, `foreach` on Windows) with a single unified
  `foreach` + `doRNG::%dorng%` path. `pbmcapply` has been removed as a
  dependency; `doRNG` added.
* `multi_trial()`: fixed an error at early interim looks with low
  prevalence, when zero observations in a contingency table cell caused
  `table()` to return a non-2x2 result.
* `binom_sample_size()`: corrected the normal approximation formula to
  match the cited Chow et al. (2017) reference. The exact search result
  is unchanged.
* Improved prior parameter documentation to explain the implications of
  the Beta(0.1, 0.1) default.
