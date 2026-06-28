# Calculate the minimum number of samples required for a one-sided exact binomial test

Calculate the minimum number of samples required for a one-sided exact
binomial test to distinguish between two success probabilities with
specified alpha and power.

## Usage

``` r
binom_sample_size(alpha = 0.05, power = 0.9, p0 = 0.9, p1 = 0.95)
```

## Arguments

- alpha:

  scalar. The desired false positive rate (probability of incorrectly
  rejecting the null). Must be be between 0 and 1. Default value is
  `alpha = 0.05`.

- power:

  scalar. The minimum probability of correctly rejecting the null when
  the alternate is true.

- p0:

  scalar. The expected proportion of successes under the null.

- p1:

  scalar. The proportion of successes under the alternate hypothesis.

## Value

A list containing the required sample size and the number of successful
trials required.

## Details

This is a one-sided function, such that \\p_0 \< p_1\\. It determines
the minimum sample size to evaluate the hypothesis test:

\$\$H_0: \\ p_1 \le p_0, \\ vs.\$\$ \$\$H_1: \\ p_1 \> p_0\$\$

## References

Chow S-C, Shao J, Wang H, Lokhnygina Y. (2017) *Sample Size Calculations
in Clinical Research*, Boca Raton, FL: CRC Press.

## Examples

``` r

# The minimum number of reference positive cases required to demonstrate
# the true sensitivity is >0.7, assuming that the true value is 0.824, with
# 90% power is

binom_sample_size(alpha = 0.05, power = 0.9, p0 = 0.7, p1 = 0.824)
#> $successes
#> [1] 81
#> 
#> $N
#> [1] 104
#> 

# With a sample size of n = 104, if the true prevalence is 0.2, we would
# require a sample size of at least n = 520 randomly sampled subjects to
# have adequate power to demonstrate the sensitivity of the new test.

# The minimum number of reference negative cases required to demonstrate
# the true specificity is >0.9, assuming that the true value is 0.963, with
# 90% power is

binom_sample_size(alpha = 0.05, power = 0.9, p0 = 0.9, p1 = 0.963)
#> $successes
#> [1] 134
#> 
#> $N
#> [1] 142
#> 

# The proposed total sample size of n = 520 would be sufficient to
# demonstrate both endpoint goals are met.
```
