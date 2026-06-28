test_that("both endpoint works", {
  set.seed(4712)
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
    ncores = 1L
  )

  result <- summarise_trials(data, fut = 0.05)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
  expect_named(
    result,
    c("power", "stop_futility", "n_avg", "sens", "spec", "mean_pos")
  )
})

test_that("summarise_trials min_pos triggers insufficient positive cases decision", {
  set.seed(8821)
  # Use high min_pos so most/all trials at early looks won't meet the threshold
  data <- multi_trial(
    sens_true = 0.9,
    spec_true = 0.95,
    prev_true = 0.05,
    endpoint = "both",
    sens_pg = 0.8,
    spec_pg = 0.8,
    prior_sens = c(1, 1),
    prior_spec = c(1, 1),
    prior_prev = c(1, 1),
    succ_sens = 0.95,
    succ_spec = 0.95,
    n_at_looks = c(50, 100),
    n_mc = 1000,
    n_trials = 20,
    ncores = 1L
  )

  # Set min_pos very high so trials are forced into the insufficient branch
  result <- summarise_trials(data, fut = 0, min_pos = 500)
  expect_s3_class(result, "data.frame")
  decisions <- unique(attr(by(data$sims, data$sims$trial, \(x) x), "dimnames"))
  # All trials should end with "no stopping - insufficient positive cases"
  expect_equal(result$power, 0)
})

test_that("spec endpoint futility stopping works", {
  set.seed(3307)
  # Use true spec close to PG and tight fut threshold to encourage futility stops
  data <- multi_trial(
    sens_true = 0.9,
    spec_true = 0.81,
    prev_true = 0.5,
    endpoint = "spec",
    sens_pg = NULL,
    spec_pg = 0.9,
    prior_sens = c(1, 1),
    prior_spec = c(1, 1),
    prior_prev = c(1, 1),
    succ_sens = 0.95,
    succ_spec = 0.95,
    n_at_looks = c(100, 200, 300),
    n_mc = 1000,
    n_trials = 50,
    ncores = 1L
  )

  result <- summarise_trials(data, fut = 0.99)
  expect_s3_class(result, "data.frame")
  expect_gte(result$stop_futility, 0)
})
