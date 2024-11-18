print("This is the beginning of test-stepEdge")

test_that("degree either 0 or 1", {
  expect_error(stepEdge(image = sar, bandwidth = 4, thresh = 17, degree = 2))
  expect_no_error(stepEdge(image = matrix(0, 10, 10), bandwidth = 4,
                           thresh = 17, degree = 0))
  expect_no_error(stepEdge(image = matrix(0, 10, 10), bandwidth = 4,
                           thresh = 17, degree = 0, blur = TRUE))
  expect_no_error(stepEdge(image = matrix(0, 10, 10), bandwidth = 4,
                           thresh = 17, degree = 1))
  expect_no_error(stepEdge(image = matrix(0, 10, 10), bandwidth = 4,
                           thresh = 17, degree = 1, blur = TRUE))
})

print("This is the end of test-stepEdge")
