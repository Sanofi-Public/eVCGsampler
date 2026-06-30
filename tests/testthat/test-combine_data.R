test_that("combine_data returns data frame with correct indicator column", {
  pool <- data.frame(cov1 = rnorm(50), cov2 = rnorm(50))
  tg   <- data.frame(cov1 = rnorm(20), cov2 = rnorm(20))
  result <- combine_data(pool, tg)
  expect_true(is.data.frame(result))
  expect_true("treated" %in% colnames(result))
  expect_equal(nrow(result), 70)
})

test_that("combine_data indicator is a factor with POOL and TG levels", {
  pool <- data.frame(cov1 = rnorm(50), cov2 = rnorm(50))
  tg   <- data.frame(cov1 = rnorm(20), cov2 = rnorm(20))
  result <- combine_data(pool, tg)
  expect_true(is.factor(result$treated))
  expect_equal(sort(levels(result$treated)), c("POOL", "TG"))
  expect_equal(sum(result$treated == "POOL"), 50)
  expect_equal(sum(result$treated == "TG"), 20)
})

test_that("combine_data keeps only overlapping columns", {
  pool <- data.frame(cov1 = rnorm(50), cov2 = rnorm(50), cov3 = rnorm(50))
  tg   <- data.frame(cov1 = rnorm(20), cov2 = rnorm(20))
  result <- combine_data(pool, tg)
  expect_true("cov1" %in% colnames(result))
  expect_true("cov2" %in% colnames(result))
  expect_false("cov3" %in% colnames(result))
})

test_that("combine_data custom indicator_name is used", {
  pool <- data.frame(cov1 = rnorm(30), cov2 = rnorm(30))
  tg   <- data.frame(cov1 = rnorm(10), cov2 = rnorm(10))
  result <- combine_data(pool, tg, indicator_name = "group")
  expect_true("group" %in% colnames(result))
  expect_false("treated" %in% colnames(result))
})

test_that("combine_data warns when TG is larger than POOL", {
  pool <- data.frame(cov1 = rnorm(10), cov2 = rnorm(10))
  tg   <- data.frame(cov1 = rnorm(30), cov2 = rnorm(30))
  expect_warning(combine_data(pool, tg))
})

test_that("combine_data errors when no overlapping columns exist", {
  pool <- data.frame(cov1 = rnorm(30))
  tg   <- data.frame(cov2 = rnorm(10))
  expect_error(combine_data(pool, tg))
})
