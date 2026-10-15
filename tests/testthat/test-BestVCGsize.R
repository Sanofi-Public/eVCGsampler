test_that("BestVCGsize returns a value (plot=FALSE)", {
  set.seed(42)
  dat <- data.frame(
    treat = rep(c(0, 1), c(15, 10)),
    cov1  = c(rnorm(15, 10, 2), rnorm(10, 11, 2))
  )
  result <- BestVCGsize(treat ~ cov1, data = dat, plot = FALSE)
  expect_true(is.numeric(result) || is.na(result))
})

test_that("BestVCGsize with plot=TRUE returns list with ggplot", {
  set.seed(42)
  dat <- data.frame(
    treat = rep(c(0, 1), c(15, 10)),
    cov1  = c(rnorm(15, 10, 2), rnorm(10, 11, 2))
  )
  result <- BestVCGsize(treat ~ cov1, data = dat, plot = TRUE)
  expect_true(is.list(result))
  expect_s3_class(result[[2]], "gg")
})

test_that("BestVCGsize with stratum returns table (plot=FALSE)", {
  set.seed(42)
  dat <- data.frame(
    treat = rep(c(0, 1), c(30, 15)),
    cov1  = c(rnorm(30, 10, 2), rnorm(15, 11, 2)),
    sex   = rep(c("M", "F"), 22)[seq_len(45)]
  )
  result <- BestVCGsize(treat ~ cov1 | sex, data = dat, plot = FALSE)
  expect_true(is.matrix(result) || is.data.frame(result))
  expect_true("Stratum" %in% colnames(result))
})

test_that("BestVCGsize errors on interaction formula", {
  dat <- data.frame(treat = rep(0:1, 10), x = rnorm(20), y = rnorm(20))
  expect_error(BestVCGsize(treat ~ x * y, data = dat, plot = FALSE))
})

test_that("BestVCGsize handles small pools and returns a single NA without a plateau", {
  set.seed(3)
  dat <- data.frame(treat = rep(0:1, c(4, 4)), cov1 = rnorm(8))
  expect_no_warning(res <- BestVCGsize(treat ~ cov1, data = dat, plot = FALSE))
  expect_length(res, 1)
})

test_that("BestVCGsize only scans POOL sizes with complete data", {
  set.seed(4)
  dat <- data.frame(treat = rep(0:1, c(40, 10)), cov1 = c(rnorm(40), rnorm(10, 2)))
  dat$cov1[1:5] <- NA
  expect_no_error(res <- BestVCGsize(treat ~ cov1, data = dat, plot = FALSE))
  expect_length(res, 1)
})
