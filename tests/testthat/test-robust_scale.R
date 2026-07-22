test_that("robust_scale numeric: centered on TG group (level 2) median", {
  set.seed(42)
  x <- c(rnorm(60, 10, 3), rnorm(20, 15, 2))
  group <- rep(c(0, 1), c(60, 20))
  res <- robust_scale(x, group)
  tg_idx <- which(group == 1)
  expect_true(abs(median(res[tg_idx])) < 0.05)
  expect_true(abs(mad(res[tg_idx]) - 1) < 0.05)
})

test_that("robust_scale numeric: length preserved", {
  set.seed(1)
  x <- rnorm(50, 10, 2)
  group <- rep(c(0, 1), 25)
  res <- robust_scale(x, group)
  expect_equal(length(res), 50)
})

test_that("robust_scale binary factor (nn=2): correct transformation", {
  x <- c(1, 2, 1, 2, 1, 2)
  group <- c(0, 0, 0, 1, 1, 1)
  res <- robust_scale(x, group)
  expected_vals <- sort(unique(((c(1, 2) - 1.5) * 2) * 0.675))
  expect_equal(sort(unique(res)), expected_vals)
})

test_that("robust_scale 3-level factor (nn=3): correct transformation", {
  x <- c(1, 2, 3, 1, 2, 3)
  group <- rep(c(0, 1), 3)
  res <- robust_scale(x, group)
  expected <- (c(1, 2, 3, 1, 2, 3) - 2) * 0.675
  expect_equal(res, expected)
})

test_that("robust_scale 4-level factor (nn=4): correct transformation", {
  x <- c(1, 2, 3, 4, 1, 2, 3, 4)
  group <- rep(c(0, 1), 4)
  res <- robust_scale(x, group)
  expected <- (c(1, 2, 3, 4, 1, 2, 3, 4) - 2.5) / 1.5
  expect_equal(res, expected)
})

test_that("robust_scale data frame: returns data frame with same dimensions", {
  set.seed(42)
  df <- data.frame(x = rnorm(50, 10, 2), y = rnorm(50, 5, 1))
  group <- rep(c(0, 1), c(30, 20))
  res <- robust_scale(df, group)
  expect_true(is.data.frame(res))
  expect_equal(ncol(res), 2)
  expect_equal(nrow(res), 50)
})

test_that("robust_scale data frame: each column is independently scaled", {
  set.seed(42)
  df <- data.frame(x = c(rnorm(40, 10, 2), rnorm(20, 15, 2)),
                   y = c(rnorm(40, 5, 1),  rnorm(20, 8, 1)))
  group <- rep(c(0, 1), c(40, 20))
  res <- robust_scale(df, group)
  tg_idx <- which(group == 1)
  expect_true(abs(median(res$x[tg_idx])) < 0.05)
  expect_true(abs(median(res$y[tg_idx])) < 0.05)
})
