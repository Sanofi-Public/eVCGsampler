test_that("VCG_sampler returns data frame with expected columns (plot=FALSE)", {
  set.seed(42)
  dat <- data.frame(
    treated = rep(c(0, 1), c(30, 15)),
    cov1 = c(rnorm(30, 10, 2), rnorm(15, 11, 2)),
    cov2 = c(rnorm(30, 5, 1),  rnorm(15, 6, 1))
  )
  result <- VCG_sampler(treated ~ cov1 + cov2, data = dat, n = 10, plot = FALSE)
  expect_true(is.data.frame(result))
  expect_true("VCG" %in% colnames(result))
  expect_true("e_weights" %in% colnames(result))
  expect_true("treated_balanced" %in% colnames(result))
})

test_that("VCG_sampler selects exactly n pool observations (deterministic)", {
  set.seed(42)
  dat <- data.frame(
    treated = rep(c(0, 1), c(30, 15)),
    cov1 = c(rnorm(30, 10, 2), rnorm(15, 11, 2))
  )
  result <- VCG_sampler(treated ~ cov1, data = dat, n = 10, plot = FALSE)
  expect_equal(sum(result$VCG == 1, na.rm = TRUE), 10)
})

test_that("VCG_sampler deterministic mode is reproducible", {
  dat <- data.frame(
    treated = rep(c(0, 1), c(30, 15)),
    cov1 = c(rnorm(30, 10, 2), rnorm(15, 11, 2))
  )
  r1 <- VCG_sampler(treated ~ cov1, data = dat, n = 8, random = FALSE, plot = FALSE)
  r2 <- VCG_sampler(treated ~ cov1, data = dat, n = 8, random = FALSE, plot = FALSE)
  expect_equal(r1$VCG, r2$VCG)
})

test_that("VCG_sampler random mode selects exactly n observations", {
  set.seed(42)
  dat <- data.frame(
    treated = rep(c(0, 1), c(30, 15)),
    cov1 = c(rnorm(30, 10, 2), rnorm(15, 11, 2))
  )
  result <- VCG_sampler(treated ~ cov1, data = dat, n = 8, random = TRUE, plot = FALSE)
  expect_equal(sum(result$VCG == 1, na.rm = TRUE), 8)
})

test_that("VCG_sampler TG rows have NA in VCG column", {
  set.seed(42)
  dat <- data.frame(
    treated = rep(c(0, 1), c(30, 15)),
    cov1 = c(rnorm(30, 10, 2), rnorm(15, 11, 2))
  )
  result <- VCG_sampler(treated ~ cov1, data = dat, n = 8, plot = FALSE)
  tg_rows <- which(dat$treated == 1)
  expect_true(all(is.na(result$VCG[tg_rows])))
})

test_that("VCG_sampler with plot=TRUE returns list of length 2 with ggplot", {
  set.seed(42)
  dat <- data.frame(
    treated = rep(c(0, 1), c(30, 15)),
    cov1 = c(rnorm(30, 10, 2), rnorm(15, 11, 2))
  )
  result <- VCG_sampler(treated ~ cov1, data = dat, n = 8, plot = TRUE)
  expect_true(is.list(result))
  expect_equal(length(result), 2)
  expect_true(is.data.frame(result[[1]]))
  expect_s3_class(result[[2]], "gg")
})

test_that("VCG_sampler with covariate weights (c_w) works", {
  set.seed(42)
  dat <- data.frame(
    treated = rep(c(0, 1), c(30, 15)),
    cov1 = c(rnorm(30, 10, 2), rnorm(15, 11, 2)),
    cov2 = c(rnorm(30, 5, 1),  rnorm(15, 6, 1))
  )
  result <- VCG_sampler(treated ~ cov1 + cov2, data = dat, n = 8,
                        c_w = c(2, 1), plot = FALSE)
  expect_equal(sum(result$VCG == 1, na.rm = TRUE), 8)
})

test_that("VCG_sampler errors on interaction formula", {
  dat <- data.frame(treated = rep(0:1, 10), x = rnorm(20), y = rnorm(20))
  expect_error(VCG_sampler(treated ~ x * y, data = dat, n = 5, plot = FALSE))
})

test_that("VCG_sampler with stratum: returns data frame", {
  set.seed(42)
  dat <- data.frame(
    treated = rep(c(0, 1), c(40, 20)),
    cov1 = c(rnorm(40, 10, 2), rnorm(20, 11, 2)),
    sex = rep(c("M", "F"), 30)
  )
  result <- VCG_sampler(treated ~ cov1 | sex, data = dat, n = 5, plot = FALSE)
  expect_true(is.data.frame(result))
  expect_true("VCG" %in% colnames(result))
})

test_that("VCG_sampler size warnings count the POOL correctly for combine_data() output", {
  set.seed(1)
  pool <- data.frame(cov1 = rnorm(60, 10), cov2 = rnorm(60, 5))
  tg   <- data.frame(cov1 = rnorm(20, 9.5), cov2 = rnorm(20, 4.8))
  dat  <- combine_data(pool, tg)
  expect_no_warning(VCG_sampler(treated ~ cov1 + cov2, data = dat, n = 20, plot = FALSE))
})

test_that("VCG_sampler stops when n is larger than the POOL", {
  dat <- data.frame(treated = rep(0:1, c(10, 10)), cov1 = rnorm(20))
  expect_error(suppressWarnings(VCG_sampler(treated ~ cov1, data = dat, n = 15, plot = FALSE)), "larger than")
  expect_error(VCG_sampler(treated ~ cov1, data = dat, n = 0, plot = FALSE), "positive whole number")
})

test_that("VCG_sampler with stratum keeps the input row order and rows with missing stratum", {
  set.seed(2)
  dat <- data.frame(
    treated = rep(0:1, c(60, 20)),
    cov1    = rnorm(80),
    sex     = rep(c("m", "f"), 40)
  )
  dat$sex[1:3] <- NA
  res <- VCG_sampler(treated ~ cov1 | sex, data = dat, n = 8, plot = FALSE)
  expect_equal(nrow(res), nrow(dat))
  expect_equal(res$cov1, dat$cov1)
  expect_true(all(is.na(res$VCG[1:3])))
  expect_equal(sum(res$VCG == 1, na.rm = TRUE), 16)
})
