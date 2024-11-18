print("This is the beginning of test-stepEdgeParSelLCK")

test_that("Only accept proper input", {
  expect_error(stepEdgeParSelLCK(image = 1:4, bandwidth = 2, thresh = 0.1,
                                 nboot = 2), "image data must be a matrix")
  expect_error(stepEdgeParSelLCK(image = matrix(0, 3, 4), bandwidth = 2,
                                 thresh = 0.1, nboot = 2),
               "image data must be a square matrix")
  expect_error(stepEdgeParSelLCK(image = matrix(0, 4, 4), bandwidth = "A",
                                 thresh = 0.1, nboot = 2),
               "bandwidth must be numeric")
  expect_error(stepEdgeParSelLCK(image = matrix(0, 4, 4), bandwidth = c(-1, 2),
                                 thresh = 0.1, nboot = 2),
               "bandwidth must be a positive integer")
  expect_error(stepEdgeParSelLCK(image = matrix(0, 4, 4), bandwidth = 2,
                                 thresh = "A", nboot = 2),
               "threshold must be numeric")
  expect_error(stepEdgeParSelLCK(image = matrix(0, 4, 4), bandwidth = 2,
                                 thresh = 0.1, nboot = 2:3),
               "nboot must be a positive integer")
  expect_no_error(parsel <- stepEdgeParSelLCK(image = sar, bandwidth = 10,
                                               thresh = 17, nboot = 2))
})

print("This is the end of test-stepEdgeParSelLCK")
