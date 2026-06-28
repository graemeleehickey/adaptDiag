# Changelog

## adaptDiag 0.1.2.9000

- [`multi_trial()`](https://graemeleehickey.github.io/adaptDiag/reference/multi_trial.md)
  gains a `seed` argument for reproducible simulations, addressing a
  user request in
  [\#5](https://github.com/graemeleehickey/adaptDiag/issues/5).
  Reproducibility is implemented via `doRNG::%dorng%`, which correctly
  handles RNG stream splitting across all parallel backends.
- [`multi_trial()`](https://graemeleehickey.github.io/adaptDiag/reference/multi_trial.md)
  now uses a single `foreach` + `doRNG::%dorng%` dispatch path for all
  platforms, replacing the previous platform-specific branches
  ([`pbmcapply::pbmclapply`](https://rdrr.io/pkg/pbmcapply/man/pbmclapply.html)
  on Unix, `foreach::%dopar%` on Windows). This fixes a latent Windows
  bug and simplifies the implementation. `pbmcapply` has been removed as
  a dependency and `doRNG` added.
- Fixed the contingency table error reported in
  [\#5](https://github.com/graemeleehickey/adaptDiag/issues/5) that
  occurred when a simulated dataset had zero observations in a cell at
  early interim looks.

## adaptDiag 0.1.1

- Updated GitHub actions workflows
- Updated README badges

## adaptDiag 0.1.0

CRAN release: 2021-08-17

- First package release to CRAN.
