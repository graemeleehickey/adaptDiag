test_that("sens only test works", {
  data <- multi_trial(
    sens_true = 0.9,
    spec_true = 0.95,
    prev_true = 0.1,
    endpoint = "sens",
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

  result <- summarise_trials(data, fut = 0.05)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
})

test_that("spec only test works", {
  data <- multi_trial(
    sens_true = 0.9,
    spec_true = 0.95,
    prev_true = 0.1,
    endpoint = "spec",
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

  result <- summarise_trials(data, fut = 0.05)
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)
})
