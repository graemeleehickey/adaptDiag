test_that("binom_sample_size returns correct structure", {
  result <- binom_sample_size(alpha = 0.05, power = 0.9, p0 = 0.7, p1 = 0.824)
  expect_type(result, "list")
  expect_named(result, c("successes", "N"))
  expect_true(is.numeric(result$N))
  expect_true(is.numeric(result$successes))
})

test_that("binom_sample_size N is at least as large as successes", {
  result <- binom_sample_size(alpha = 0.05, power = 0.9, p0 = 0.7, p1 = 0.824)
  expect_gte(result$N, result$successes)
})

test_that("binom_sample_size returns known reference values", {
  # From the package documentation example
  result <- binom_sample_size(alpha = 0.05, power = 0.9, p0 = 0.7, p1 = 0.824)
  expect_equal(result$N, 104L)
})

test_that("binom_sample_size errors when p0 >= p1", {
  expect_snapshot(error = TRUE, binom_sample_size(p0 = 0.9, p1 = 0.8))
  expect_snapshot(error = TRUE, binom_sample_size(p0 = 0.8, p1 = 0.8))
})

test_that("binom_sample_size increases N with stricter power", {
  low <- binom_sample_size(alpha = 0.05, power = 0.8, p0 = 0.7, p1 = 0.85)
  high <- binom_sample_size(alpha = 0.05, power = 0.95, p0 = 0.7, p1 = 0.85)
  expect_lt(low$N, high$N)
})

test_that("binom_sample_size increases N with smaller effect size", {
  large_effect <- binom_sample_size(
    alpha = 0.05,
    power = 0.9,
    p0 = 0.7,
    p1 = 0.95
  )
  small_effect <- binom_sample_size(
    alpha = 0.05,
    power = 0.9,
    p0 = 0.7,
    p1 = 0.75
  )
  expect_lt(large_effect$N, small_effect$N)
})
