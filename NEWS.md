# adaptDiag 0.1.1

* `binom_sample_size()`: fixed the normal approximation used to initialise the sample size search. The previous formula used only $p_1$ in the variance term; the correct Chow et al. (2017) formula uses $z_\alpha\sqrt{p_0(1-p_0)} + z_\beta\sqrt{p_1(1-p_1)}$. The exact discrete search result is unchanged, but the approximation now matches the cited reference.
* `multi_trial()`: improved documentation for `prior_sens`, `prior_spec`, and `prior_prev` to explain that `c(0.1, 0.1)` is a bimodal near-Jeffreys prior that places most mass near 0 and 1, and that `c(1, 1)` (uniform) or an informative prior from pilot data may be more appropriate.
* `multi_trial()`: added a `seed` argument for reproducible simulations (#5). Reproducibility is implemented via `doRNG::%dorng%`, which correctly handles RNG stream splitting across all parallel backends.
* `multi_trial()` now uses a single `foreach` + `doRNG::%dorng%` dispatch path for all platforms, replacing the previous platform-specific branches (`pbmcapply::pbmclapply` on Unix, `foreach::%dopar%` on Windows). This fixes a latent Windows bug and simplifies the implementation. `pbmcapply` has been removed as a dependency and `doRNG` added.
* Fixed an error in `multi_trial()` that occurred at early interim looks with low prevalence, when zero observations fell into a cell of the 2×2 contingency table (#5).

# adaptDiag 0.1.0

* First package release to CRAN.
