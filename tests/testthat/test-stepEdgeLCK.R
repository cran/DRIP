print("This is the beginning of test-stepEdgeLCK")

test_that("only accept proper input", {
  expect_error(stepEdgeLCK(image = 1:16, bandwidth = 2, thresh = 0.1),
               "image data must be a matrix")
  expect_error(stepEdgeLCK(image = matrix(10, 3, 4), bandwidth = 2,
                           thresh = 0.1), "image data must be a square matrix")
  expect_error(stepEdgeLCK(image = matrix(0, 3, 3), bandwidth = "A",
                           thresh = 0.1), "bandwidth must be numeric")
  expect_error(stepEdgeLCK(image = matrix(0, 3, 3), bandwidth = -1,
                           thresh = 0.1),
               "bandwidth must be a positive integer")
  expect_error(stepEdgeLCK(image = matrix(0, 3, 3), bandwidth = 2,
                           thresh = "A"), "threshold  must be numeric")
  expect_no_error(stepEdgeLCK(image = matrix(0, 3, 3), bandwidth = 2,
                              thresh = 10, plot = TRUE))
})

print("This is the end of test-stepEdgeLCK")
