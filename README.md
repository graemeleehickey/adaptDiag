
<!-- README.md is generated from README.Rmd. Please edit that file -->

# adaptDiag <img src="man/figures/hex.png" alt="adaptDiag hex logo" width = "175" height = "200" align="right" />

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/adaptDiag)](https://CRAN.R-project.org/package=adaptDiag)
[![CRAN
downloads](https://cranlogs.r-pkg.org/badges/adaptDiag)](https://CRAN.R-project.org/package=adaptDiag)
[![Codecov test
coverage](https://codecov.io/gh/graemeleehickey/adaptDiag/graph/badge.svg)](https://app.codecov.io/gh/graemeleehickey/adaptDiag)
[![R-CMD-check](https://github.com/graemeleehickey/adaptDiag/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/graemeleehickey/adaptDiag/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/graemeleehickey/adaptDiag/actions/workflows/pkgdown.yaml/badge.svg)](https://graemeleehickey.github.io/adaptDiag/)
[![License: GPL
v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

The goal of `adaptDiag` is to simplify the process of designing adaptive
trials for diagnostic test studies. With accumulating data in a clinical
trial of a new diagnostic test compared to a gold-standard reference,
decisions can be made at interim analyses to either stop the trial for
early success, stop the trial for expected futility, or continue to the
next sample size look. Designs can be focused around test sensitivity,
specificity, or both. The package is heavily influenced by the seminal
article by Broglio et al. (2014).

## References

Broglio KR, Connor JT, Berry SM. Not too big, not too small: a
Goldilocks approach to sample size selection. *Journal of
Biopharmaceutical Statistics*, 2014; **24(3)**: 685–705.

## Installation

You can install the development version of `adaptDiag`
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("graemeleehickey/adaptDiag")
```
