test_that("combine_variables returns numeric vector of correct length", {
  set.seed(42)
  dat <- data.frame(x = rnorm(50), y = rnorm(50))
  result <- combine_variables(1 ~ x + y, data = dat)
  expect_true(is.numeric(result))
  expect_equal(length(result), 50)
})

test_that("combine_variables NULL weights equals all-ones weights", {
  set.seed(42)
  dat <- data.frame(x = rnorm(50), y = rnorm(50))
  r1 <- combine_variables(1 ~ x + y, data = dat)
  r2 <- combine_variables(1 ~ x + y, data = dat, weights = c(1, 1))
  expect_equal(r1, r2)
})

test_that("combine_variables higher weight increases variable contribution", {
  set.seed(42)
  dat <- data.frame(x = rnorm(50), y = rnorm(50))
  r_eq  <- combine_variables(1 ~ x + y, data = dat, weights = c(1, 1))
  r_wt  <- combine_variables(1 ~ x + y, data = dat, weights = c(5, 1))
  x_scaled <- as.numeric(scale(dat$x))
  expect_true(cor(r_wt, x_scaled) > cor(r_eq, x_scaled))
})

test_that("combine_variables errors on wrong weight length", {
  dat <- data.frame(x = rnorm(50), y = rnorm(50))
  expect_error(combine_variables(1 ~ x + y, data = dat, weights = c(1)))
})

test_that("combine_variables errors on non-positive weights", {
  dat <- data.frame(x = rnorm(50), y = rnorm(50))
  expect_error(combine_variables(1 ~ x + y, data = dat, weights = c(1, -1)))
  expect_error(combine_variables(1 ~ x + y, data = dat, weights = c(1, 0)))
})

test_that("combine_variables errors on non-finite weights", {
  dat <- data.frame(x = rnorm(50), y = rnorm(50))
  expect_error(combine_variables(1 ~ x + y, data = dat, weights = c(1, Inf)))
})

test_that("combine_variables errors when variable not in data", {
  dat <- data.frame(x = rnorm(50))
  expect_error(combine_variables(1 ~ x + z, data = dat))
})

test_that("combine_variables errors on stratum syntax", {
  dat <- data.frame(x = rnorm(50), y = rnorm(50), s = rep(1:2, 25))
  expect_error(combine_variables(1 ~ x + y | s, data = dat))
})

test_that("combine_variables warns on NA in data", {
  dat <- data.frame(x = c(NA, rnorm(49)), y = rnorm(50))
  expect_warning(combine_variables(1 ~ x + y, data = dat))
})
