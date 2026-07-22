test_that("plot_var returns ggplot for continuous variable", {
  set.seed(42)
  dat <- data.frame(
    treated = rep(c(0, 1), c(30, 15)),
    cov1 = c(rnorm(30, 10, 2), rnorm(15, 11, 2))
  )
  out <- VCG_sampler(treated ~ cov1, data = dat, n = 8, plot = FALSE)
  p <- plot_var(out, what = "cov1", group = "VCG")
  expect_s3_class(p, "gg")
})

test_that("plot_var returns ggplot for categorical variable (2 levels)", {
  set.seed(42)
  dat <- data.frame(
    treated = rep(c(0, 1), c(30, 15)),
    cov1 = c(rnorm(30, 10, 2), rnorm(15, 11, 2)),
    sex  = rep(0:1, 22)[seq_len(45)]
  )
  out <- VCG_sampler(treated ~ cov1, data = dat, n = 8, plot = FALSE)
  p <- plot_var(out, what = "sex", group = "VCG")
  expect_s3_class(p, "gg")
})

test_that("plot_var accepts optional title", {
  set.seed(42)
  dat <- data.frame(
    treated = rep(c(0, 1), c(30, 15)),
    cov1 = c(rnorm(30, 10, 2), rnorm(15, 11, 2))
  )
  out <- VCG_sampler(treated ~ cov1, data = dat, n = 8, plot = FALSE)
  p <- plot_var(out, what = "cov1", group = "VCG", title = "My Title")
  expect_s3_class(p, "gg")
})
