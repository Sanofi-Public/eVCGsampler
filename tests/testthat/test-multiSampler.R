test_that("multiSampler returns data frame with VCG columns (plot=FALSE)", {
  set.seed(42)
  dat <- data.frame(
    treat = rep(c(0, 1), c(30, 15)),
    cov1  = c(rnorm(30, 10, 2), rnorm(15, 11, 2))
  )
  result <- multiSampler(treat ~ cov1, data = dat, n = 5, Nsamples = 3, plot = FALSE)
  expect_true(is.data.frame(result))
  expect_true("VCG_1" %in% colnames(result))
  expect_true("VCG_3" %in% colnames(result))
})

test_that("multiSampler VCG columns contain only 0/1/NA", {
  set.seed(42)
  dat <- data.frame(
    treat = rep(c(0, 1), c(30, 15)),
    cov1  = c(rnorm(30, 10, 2), rnorm(15, 11, 2))
  )
  result <- multiSampler(treat ~ cov1, data = dat, n = 5, Nsamples = 3, plot = FALSE)
  vcg_vals <- unlist(result[, paste0("VCG_", 1:3)])
  expect_true(all(vcg_vals %in% c(0, 1, NA)))
})

test_that("multiSampler with plot=TRUE returns list with data and ggplot", {
  set.seed(42)
  dat <- data.frame(
    treat = rep(c(0, 1), c(30, 15)),
    cov1  = c(rnorm(30, 10, 2), rnorm(15, 11, 2))
  )
  result <- multiSampler(treat ~ cov1, data = dat, n = 5, Nsamples = 3, plot = TRUE)
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_true(is.data.frame(result[[1]]))
  expect_s3_class(result[[2]], "gg")
})

test_that("multiSampler errors on interaction formula", {
  dat <- data.frame(treat = rep(0:1, 10), x = rnorm(20), y = rnorm(20))
  expect_error(multiSampler(treat ~ x * y, data = dat, n = 3, Nsamples = 2, plot = FALSE))
})

test_that("multiSampler with stratum returns data frame", {
  set.seed(42)
  dat <- data.frame(
    treat = rep(c(0, 1), c(40, 20)),
    cov1  = c(rnorm(40, 10, 2), rnorm(20, 11, 2)),
    sex   = rep(c("M", "F"), 30)
  )
  result <- multiSampler(treat ~ cov1 | sex, data = dat, n = 4, Nsamples = 3, plot = FALSE)
  expect_true(is.data.frame(result))
})
