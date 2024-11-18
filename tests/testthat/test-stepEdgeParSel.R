print("This is the beginning of test-stepEdgeParSel")

test_that("check type argument works", {
  foo <- list(dKQ = matrix(0.0, nrow = 2, ncol = 2), bandwidth = as.integer(3),
              threshold = 1.2)
  class(foo) <- "Edge_Parameters"
  expect_error(print.Edge_Parameters(x = foo, type = "bandwidth"),
               "Wrong type value")
  expect_no_error(print.Edge_Parameters(x = foo, type = "matrix"))
  expect_no_error(print.Edge_Parameters(x = foo, type = "parameters"))
  expect_no_error(print.Edge_Parameters(x = foo, type = "all"))
})

test_that("accepts S3 Edge_Parameters only", {
  foo <- list(dKQ = matrix(0.0, nrow = 2, ncol = 2), bandwidth = as.integer(3),
              threshold = 1.2)
  expect_error(summary.Edge_Parameters(object = foo))
  class(foo) <- "Edge_Parameters"
  expect_no_error(summary.Edge_Parameters(object = foo))
})

test_that("returns S3 class Edge_Parameters", {
  set.seed(100)
  img <- matrix(rnorm(100), 10, 10)
  parSel <- stepEdgeParSel(image = img, bandwidth = 2, thresh = c(0.1, 0.2),
                           nboot = 1)
  expect_s3_class(parSel, class = "Edge_Parameters")
  expect_no_error(parSel <- stepEdgeParSel(image = img, bandwidth = 2,
                                           thresh = c(0.1, 0.2), nboot = 1,
                                           degree = 0))
  expect_no_error(parSel <- stepEdgeParSel(image = img, bandwidth = 2,
                                           thresh = c(0.1, 0.2), nboot = 1,
                                           degree = 1, blur = TRUE))
  expect_no_error(parSel <- stepEdgeParSel(image = img, bandwidth = 2,
                                           thresh = c(0.1, 0.2), nboot = 1,
                                           degree = 0, blur = TRUE))
})

print("This is the end of test-stepEdgeParSel")
