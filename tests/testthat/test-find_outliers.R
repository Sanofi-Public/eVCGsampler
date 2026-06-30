test_that("find_outliers returns list with expected elements (plot=FALSE)", {
  set.seed(42)
  dat <- data.frame(
    study = factor(rep(paste0("S", 1:5), each = 15)),
    x = c(rnorm(15, 10, 1), rnorm(15, 10, 1), rnorm(15, 10, 1),
          rnorm(15, 10, 1), rnorm(15, 20, 1))
  )
  result <- find_outliers(study ~ x, data = dat, R = 100, plot = FALSE)
  expect_true(is.list(result))
  expect_true("cutoff_value" %in% names(result))
  expect_true("summary" %in% names(result))
})

test_that("find_outliers summary has required columns", {
  set.seed(42)
  dat <- data.frame(
    study = factor(rep(paste0("S", 1:4), each = 10)),
    x = rnorm(40, 10, 2)
  )
  result <- find_outliers(study ~ x, data = dat, R = 100, plot = FALSE)
  expect_true(all(c("group", "median_distance", "is_outlier") %in% colnames(result$summary)))
  expect_equal(nrow(result$summary), 4)
})

test_that("find_outliers detects a clear outlier group", {
  set.seed(42)
  dat <- data.frame(
    study = factor(rep(paste0("S", 1:5), each = 20)),
    x = c(rnorm(20, 10, 1), rnorm(20, 10, 1), rnorm(20, 10, 1),
          rnorm(20, 10, 1), rnorm(20, 40, 1))
  )
  result <- find_outliers(study ~ x, data = dat, R = 100, plot = FALSE)
  s5_row <- result$summary[result$summary$group == "S5", ]
  expect_true(s5_row$is_outlier)
})

test_that("find_outliers with plot=TRUE returns heatmap and barplot", {
  set.seed(42)
  dat <- data.frame(
    study = factor(rep(paste0("S", 1:4), each = 15)),
    x = c(rnorm(15, 10, 1), rnorm(15, 10, 1), rnorm(15, 10, 1), rnorm(15, 20, 1))
  )
  result <- find_outliers(study ~ x, data = dat, R = 100, plot = TRUE)
  expect_s3_class(result$heatmap, "gg")
  expect_s3_class(result$barplot, "gg")
})

test_that("find_outliers cutoff_value is numeric", {
  set.seed(42)
  dat <- data.frame(
    study = factor(rep(paste0("S", 1:4), each = 10)),
    x = rnorm(40, 10, 2)
  )
  result <- find_outliers(study ~ x, data = dat, R = 100, plot = FALSE)
  expect_true(is.numeric(result$cutoff_value))
  expect_equal(length(result$cutoff_value), 1)
})

test_that("find_outliers errors on stratum formula", {
  dat <- data.frame(study = rep(letters[1:3], 10), x = rnorm(30), s = rep(1:2, 15))
  expect_error(find_outliers(study ~ x | s, data = dat))
})

test_that("find_outliers errors with fewer than 2 groups", {
  dat <- data.frame(study = rep("A", 20), x = rnorm(20))
  expect_error(find_outliers(study ~ x, data = dat, R = 10, plot = FALSE))
})
