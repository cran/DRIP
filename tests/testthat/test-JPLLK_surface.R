print("This is the beginning of test-JPLLK_surface")

test_that("only accept S3 class JPLLK_Parameters", {
  foo <- list(fitted = matrix(0.0, nrow = 2, ncol = 0), bandwidth = 3.0, resid
              = matrix(0.0, nrow = 2, ncol = 0), sigma = 0.0)
  expect_error(print.JPLLK_Parameters(x = foo, type = "sigma"))
})

test_that("only accept correct type values", {
  foo <- list(fitted = matrix(0.0, nrow = 2, ncol = 0), bandwidth = 3.0, resid
              = matrix(0.0, nrow = 2, ncol = 0), sigma = 0.0)
  class(foo) <- "JPLLK_Parameters"
  expect_error(print.JPLLK_Parameters(x = foo, type = "thresh"),
               "Wrong values for type")
  expect_no_error(print.JPLLK_Parameters(x = foo, type = "bandwidth"))
  expect_no_error(print.JPLLK_Parameters(x = foo, type = "sigma"))
  expect_no_error(print.JPLLK_Parameters(x = foo, type = "all"))
})

test_that("only accept S3 class JPLLK_Parameters", {
  foo <- list(fitted = matrix(0.0, nrow = 2, ncol = 0), bandwidth = 3.0, resid
              = matrix(0.0, nrow = 2, ncol = 0), sigma = 0.0)
  expect_error(summary.JPLLK_Parameters(object = foo))
  class(foo) <- "JPLLK_Parameters"
  expect_no_error(summary.JPLLK_Parameters(object = foo))
})

test_that("only accept S3 class JPLLK_Parameters", {
  foo <- list(fitted = matrix(0.0, nrow = 2, ncol = 2), bandwidth = 3.0, resid
              = matrix(0.0, nrow = 2, ncol = 2), sigma = 0.0)
  expect_error(plot.JPLLK_Parameters(x = foo))
  class(foo) <- "JPLLK_Parameters"
  expect_no_error(plot.JPLLK_Parameters(x = foo))
})

test_that("An S3 class JPLLK_Parameters is returned", {
  fit <- JPLLK_surface(image = sar, bandwidth = c(3, 4), plot = FALSE)
  expect_identical(class(fit), "JPLLK_Parameters")
  expect_no_error(parsel <- JPLLK_surface(image = matrix(0, 3, 3),
                                          bandwidth = 2, plot = TRUE))
})

test_that("only accept proper input", {
  expect_error(JPLLK_surface(image = 1:4, bandwidth = 2),
               "image data must be a matrix")
  expect_error(JPLLK_surface(image = matrix(0, 3, 4), bandwidth = 2),
               "image data must be a square matrix")
  expect_error(JPLLK_surface(image = matrix(1, 3, 3), bandwidth = "A"),
               "bandwidth must be numeric")
})

print("This is the end of test-JPLLK_surface")
