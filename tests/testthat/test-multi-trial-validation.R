test_that("multi_trial errors on invalid endpoint", {
  expect_snapshot(
    error = TRUE,
    multi_trial(
      sens_true = 0.9,
      spec_true = 0.95,
      prev_true = 0.1,
      endpoint = "invalid",
      n_at_looks = c(100, 200),
      n_trials = 2,
      ncores = 1
    )
  )
})

test_that("multi_trial errors when sens_pg missing for 'both' endpoint", {
  expect_snapshot(
    error = TRUE,
    multi_trial(
      sens_true = 0.9,
      spec_true = 0.95,
      prev_true = 0.1,
      endpoint = "both",
      sens_pg = NULL,
      spec_pg = 0.8,
      succ_sens = 0.95,
      succ_spec = 0.95,
      n_at_looks = c(100, 200),
      n_trials = 2,
      ncores = 1
    )
  )
})

test_that("multi_trial errors when succ thresholds missing for 'both' endpoint", {
  expect_snapshot(
    error = TRUE,
    multi_trial(
      sens_true = 0.9,
      spec_true = 0.95,
      prev_true = 0.1,
      endpoint = "both",
      sens_pg = 0.8,
      spec_pg = 0.8,
      succ_sens = NULL,
      succ_spec = 0.95,
      n_at_looks = c(100, 200),
      n_trials = 2,
      ncores = 1
    )
  )
})

test_that("multi_trial errors when sens_pg missing for 'sens' endpoint", {
  expect_snapshot(
    error = TRUE,
    multi_trial(
      sens_true = 0.9,
      spec_true = 0.95,
      prev_true = 0.1,
      endpoint = "sens",
      sens_pg = NULL,
      n_at_looks = c(100, 200),
      n_trials = 2,
      ncores = 1
    )
  )
})

test_that("multi_trial errors when spec_pg missing for 'spec' endpoint", {
  expect_snapshot(
    error = TRUE,
    multi_trial(
      sens_true = 0.9,
      spec_true = 0.95,
      prev_true = 0.1,
      endpoint = "spec",
      spec_pg = NULL,
      n_at_looks = c(100, 200),
      n_trials = 2,
      ncores = 1
    )
  )
})

test_that("multi_trial errors when prior arguments are NULL", {
  expect_snapshot(
    error = TRUE,
    multi_trial(
      sens_true = 0.9,
      spec_true = 0.95,
      prev_true = 0.1,
      endpoint = "both",
      sens_pg = 0.8,
      spec_pg = 0.8,
      succ_sens = 0.95,
      succ_spec = 0.95,
      prior_sens = NULL,
      prior_spec = c(1, 1),
      prior_prev = c(1, 1),
      n_at_looks = c(100, 200),
      n_trials = 2,
      ncores = 1
    )
  )
})

test_that("multi_trial warns and corrects when ncores < 1", {
  expect_warning(
    multi_trial(
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
      n_at_looks = c(500),
      n_mc = 100,
      n_trials = 2,
      ncores = 0
    ),
    "Must use at least 1 core"
  )
})

test_that("multi_trial is reproducible with a seed", {
  args <- list(
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
    n_at_looks = c(200, 400),
    n_mc = 500,
    n_trials = 5,
    seed = 42L,
    ncores = 1
  )
  result1 <- do.call(multi_trial, args)
  result2 <- do.call(multi_trial, args)
  expect_equal(result1$sims, result2$sims)
})

test_that("multi_trial uses ncores default when missing", {
  result <- multi_trial(
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
    n_at_looks = c(500),
    n_mc = 100,
    n_trials = 2,
    ncores = 1
  )
  expect_type(result, "list")
  expect_named(result, c("sims", "call", "args"))
})
