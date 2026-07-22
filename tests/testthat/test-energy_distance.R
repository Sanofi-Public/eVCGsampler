test_that("energy_distance returns numeric in [0, 1] for standardized=TRUE", {
  set.seed(1)
  dat <- data.frame(
    treated = rep(0:1, c(40, 20)),
    x = c(rnorm(40, 10, 2), rnorm(20, 12, 2))
  )
  result <- energy_distance(treated ~ x, data = dat)
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("energy_distance ~0 for identical distributions", {
  set.seed(42)
  x <- rnorm(60, 10, 2)
  dat <- data.frame(treated = rep(0:1, c(30, 30)), x = x)
  result <- energy_distance(treated ~ x, data = dat)
  expect_true(result < 0.15)
})

test_that("energy_distance higher for clearly separated groups", {
  set.seed(42)
  dat_diff <- data.frame(
    treated = rep(0:1, c(40, 20)),
    x = c(rnorm(40, 10, 1), rnorm(20, 20, 1))
  )
  dat_same <- data.frame(
    treated = rep(0:1, c(40, 20)),
    x = c(rnorm(40, 10, 1), rnorm(20, 10, 1))
  )
  d_diff <- energy_distance(treated ~ x, data = dat_diff)
  d_same <- energy_distance(treated ~ x, data = dat_same)
  expect_true(d_diff > d_same)
})

test_that("energy_distance non-standardized returns numeric", {
  set.seed(42)
  dat <- data.frame(
    treated = rep(0:1, c(30, 15)),
    x = c(rnorm(30, 10, 1), rnorm(15, 20, 1))
  )
  result <- energy_distance(treated ~ x, data = dat, standardized = FALSE)
  expect_true(is.numeric(result))
})

test_that("energy_distance multi-covariate works", {
  set.seed(42)
  dat <- data.frame(
    treated = rep(0:1, c(30, 15)),
    x = c(rnorm(30, 10, 2), rnorm(15, 12, 2)),
    y = c(rnorm(30, 5, 1),  rnorm(15, 6, 1))
  )
  result <- energy_distance(treated ~ x + y, data = dat)
  expect_true(is.numeric(result))
  expect_true(result >= 0 && result <= 1)
})

test_that("energy_distance errors on stratum formula", {
  dat <- data.frame(treated = rep(0:1, 10), x = rnorm(20), s = rep(1:2, 10))
  expect_error(energy_distance(treated ~ x | s, data = dat))
})

test_that("energy_distance warns and handles NA values", {
  dat <- data.frame(
    treated = c(rep(0, 20), rep(1, 10)),
    x = c(NA, rnorm(29))
  )
  expect_warning(energy_distance(treated ~ x, data = dat))
})
