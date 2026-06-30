test_that("energy_test returns htest list (plot=FALSE)", {
  set.seed(42)
  dat <- data.frame(
    treated = rep(0:1, c(30, 15)),
    x = c(rnorm(30, 10, 2), rnorm(15, 12, 2))
  )
  result <- energy_test(treated ~ x, data = dat, R = 50, plot = FALSE)
  expect_s3_class(result, "htest")
  expect_true("p.value" %in% names(result))
  expect_true("estimate" %in% names(result))
  expect_true("critical.value" %in% names(result))
  expect_true("n.permutations" %in% names(result))
  expect_equal(result$n.permutations, 50)
})

test_that("energy_test p-value is between 0 and 1", {
  set.seed(1)
  dat <- data.frame(
    treated = rep(0:1, c(30, 15)),
    x = c(rnorm(30, 10, 2), rnorm(15, 10, 2))
  )
  result <- energy_test(treated ~ x, data = dat, R = 50, plot = FALSE)
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})

test_that("energy_test permutation vector has length R", {
  set.seed(42)
  dat <- data.frame(
    treated = rep(0:1, c(30, 15)),
    x = c(rnorm(30, 10, 2), rnorm(15, 12, 2))
  )
  result <- energy_test(treated ~ x, data = dat, R = 80, plot = FALSE)
  expect_equal(length(result$permutations), 80)
})

test_that("energy_test with plot=TRUE returns list with ggplot", {
  set.seed(42)
  dat <- data.frame(
    treated = rep(0:1, c(30, 15)),
    x = c(rnorm(30, 10, 2), rnorm(15, 12, 2))
  )
  result <- energy_test(treated ~ x, data = dat, R = 30, plot = TRUE)
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_s3_class(result[[2]], "gg")
})

test_that("energy_test with stratum returns table/matrix", {
  set.seed(42)
  # Each stratum must contain both treated (0) and TG (1) observations
  dat <- data.frame(
    treated = c(rep(0, 15), rep(1, 8), rep(0, 15), rep(1, 7)),
    x       = c(rnorm(15, 10, 2), rnorm(8, 12, 2), rnorm(15, 10, 2), rnorm(7, 12, 2)),
    stratum = c(rep("A", 23), rep("B", 22))
  )
  result <- suppressWarnings(
    energy_test(treated ~ x | stratum, data = dat, R = 30, plot = FALSE)
  )
  expect_true(is.matrix(result) || is.data.frame(result))
  expect_true("Stratum" %in% colnames(result))
})

test_that("energy_test errors on interaction formula", {
  dat <- data.frame(treated = rep(0:1, 10), x = rnorm(20), y = rnorm(20))
  expect_error(energy_test(treated ~ x * y, data = dat, R = 10, plot = FALSE))
})
