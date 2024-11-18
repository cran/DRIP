print("This is the beginning of cv.jpex")

test_that("check type argument works", {
  foo <- list(LLK = matrix(0.0, nrow = 2, ncol = 2), sigma = 0.0,
              cv_scores = c(1.2, 1.1), input = as.integer(3:4),
              bandwidth = as.integer(3))
  class(foo) <- "JPEX_Parameters"
  expect_error(print.JPEX_Parameters(x = foo, type = "numeric"),
               "Wrong type value")
  expect_no_error(print.JPEX_Parameters(x = foo, type = "cv_scores"))
  expect_no_error(print.JPEX_Parameters(x = foo, type = "bandwidth"))
  expect_no_error(print.JPEX_Parameters(x = foo, type = "sigma"))
  expect_no_error(print.JPEX_Parameters(x = foo, type = "all"))
  expect_no_error(print.JPEX_Parameters(x = foo))
})

test_that("accepts S3 JPEX_Parameters only", {
  foo <- list(LLK = matrix(0.0, nrow = 2, ncol = 2), sigma = 0.0,
              cv_scores = c(1.2, 1.1), input = as.integer(3:4),
              bandwidth = as.integer(3))
  expect_error(summary.JPEX_Parameters(object = foo))
  class(foo) <- "JPEX_Parameters"
  expect_no_error(summary.JPEX_Parameters(object = foo))
})

test_that("accepts S3 JPEX_Parameters only", {
  foo <- list(LLK = matrix(0.0, nrow = 2, ncol = 2), sigma = 0.0,
              cv_scores = c(1.2, 1.1), input = as.integer(3:4),
              bandwidth = as.integer(3))
  expect_error(plot.JPEX_Parameters(x = foo))
  class(foo) <- "JPEX_Parameters"
  expect_no_error(plot.JPEX_Parameters(x = foo))
})

test_that("only accept proper input type", {
  expect_error(cv.jpex(1:100, c(2, 3)), "image data must be a matrix")
  expect_error(cv.jpex(matrix(1:6, 2, 3), c(2, 3)),
               "image data must be a square matrix")
  expect_error(cv.jpex(matrix(1:4, 2, 2), letters[1:3]),
               "bandwidths must be numeric")
})

test_that("returns S3 class JPEX_Parameters", {
  set.seed(100)
  img <- matrix(rnorm(100), 10, 10)
  out <- cv.jpex(img, c(2,3))
  expect_s3_class(out, class = "JPEX_Parameters", exact = TRUE)
})

print("This is the end of cv.jpex")
