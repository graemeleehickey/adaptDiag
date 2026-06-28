# adaptDiag 0.1.2.9000

* `binom_sample_size()`: fixed the normal approximation used to initialise the sample size search. The previous formula used only \eqn{p_1} in the variance term; the correct Chow et al. (2017) formula uses \eqn{z_\alpha\sqrt{p_0(1-p_0)} + z_\beta\sqrt{p_1(1-p_1)}}. The exact discrete search result is unchanged, but the approximation now matches the cited reference.
* `multi_trial()`: improved documentation for `prior_sens`, `prior_spec`, and `prior_prev` to explain that `c(0.1, 0.1)` is a bimodal near-Jeffreys prior that places most mass near 0 and 1, and that `c(1, 1)` (uniform) or an informative prior from pilot data may be more appropriate.
* `multi_trial()`: added a `seed` argument for reproducible simulations, addressing a user request in #5. Reproducibility is implemented via `doRNG::%dorng%`, which correctly handles RNG stream splitting across all parallel backends.
* `multi_trial()` now uses a single `foreach` + `doRNG::%dorng%` dispatch path for all platforms, replacing the previous platform-specific branches (`pbmcapply::pbmclapply` on Unix, `foreach::%dopar%` on Windows). This fixes a latent Windows bug and simplifies the implementation. `pbmcapply` has been removed as a dependency and `doRNG` added.
* Fixed the contingency table error reported in #5 that occurred when a simulated dataset had zero observations in a cell at early interim looks.



# adaptDiag 0.1.1

* Updated GitHub actions workflows
* Updated README badges

# adaptDiag 0.1.0

* First package release to CRAN.
