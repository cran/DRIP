print("This is the beginning of test-diffLCK")

test_that("only accept proper input", {
  expect_error(diffLCK(image = 1:4, bandwidth = 2),
               "image data must be a matrix")
  expect_error(diffLCK(image = matrix(0, 2, 3), bandwidth = 2),
               "image data must be a square matrix")
  expect_error(diffLCK(image = matrix(0, 2, 2), bandwidth = "A"),
               "bandwidth must be numeric")
  expect_error(diffLCK(image = matrix(0, 2, 2), bandwidth = 0.2),
               "bandwidth must be a positive integer")
  expect_error(diffLCK(image = matrix(0, 2, 2), bandwidth = c(2, 3)),
               "bandwidth must be a positive integer")
  expect_error(diffLCK(image = matrix(0, 2, 2), bandwidth = 2, plot = 1),
               "plot must be TRUE or FALSE")
  expect_no_error(diffLCK(image = matrix(0, 3, 3), bandwidth = 1))
  expect_no_error(diffLCK(image = matrix(0, 3, 3), bandwidth = 1, plot = TRUE))
})

print("This is the end of test-diffLCK")
