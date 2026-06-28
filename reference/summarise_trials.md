# Summarise results of multiple simulated trials to give the operating characteristics

Summarise results of multiple simulated trials to give the operating
characteristics

## Usage

``` r
summarise_trials(data, min_pos = 1, fut = 0)
```

## Arguments

- data:

  list. Output from the
  [`multi_trial`](https://graemeleehickey.github.io/adaptDiag/reference/multi_trial.md)
  function.

- min_pos:

  integer. The minimum number of reference positive cases before
  stopping is allowed. Default is `min_pos = 1`.

- fut:

  scalar. A probability threshold at which the posterior predictive
  probability of eventual success is compared to. If the probability is
  less than `fut`, the trial stops for binding futility. Default is
  `fut = 0`, which corresponds to no stopping for futility.

## Value

A data frame of row length 1, with the following columns:

- `power`: Power is defined as the proportion of trials that result in
  success, irrespective of whether it is an early stop for success or
  not. Trials that stop for futility, but which subsequently go on to be
  successful, are not considered as a success. In other words, the
  futility decision is binding, and in practice, if a trial triggered a
  futility rule, the sponsor would not see the eventual outcome if the
  trial were to continue enrolling. When the performance goals are set
  equal to the respective true values, the power returned is the type I
  error.

- `stop_futility`: The proportion of trials that stopped early for
  expected futility.

- `n_avg`: The average sample size for trials at the stage they stopped.

- `sens`: The average sensitivity for trials at the stage they stopped.

- `spec`: The average specificity for trials at the stage they stopped.

- `mean_pos`: The average number of reference positive cases for trials
  at the stage they stopped.

## Examples

``` r
data <- multi_trial(
    sens_true = 0.9,
    spec_true = 0.95,
    prev_true = 0.1,
    endpoint = "both",
    sens_pg = 0.8,
    spec_pg = 0.8,
    prior_sens = c(1, 1),
    prior_spec = c(1, 1),
    prior_prev = c(1, 1),
    succ_sens = 0.95,
    succ_spec = 0.95,
    n_at_looks = c(200, 400, 600, 800, 1000),
    n_mc = 10000,
    n_trials = 20,
    ncores = 1
    )
#> 

summarise_trials(data, fut = 0.05, min_pos = 10)
#>                    n
#> decision            200 400 600 800 1000
#>   early win           5   6   3   2    0
#>   late win            0   0   0   0    2
#>   stop for futility   1   0   0   1    0
#> 
#>   power stop_futility n_avg      sens      spec mean_pos
#> 1   0.9           0.1   490 0.9024721 0.9525038    49.25
```
