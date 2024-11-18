print("This is the beginning of test-stepDiff")

test_that("degree either 0 or 1", {
  expect_error(stepDiff(image = sar, bandwidth = 4, degree = 2.7))
  expect_no_error(stepDiff(image = matrix(0, 10, 10), bandwidth = 3,
                           degree = 0))
  expect_no_error(stepDiff(image = matrix(0, 10, 10), bandwidth = 3, degree = 0,
                           blur = TRUE))
  expect_no_error(stepDiff(image = matrix(0, 10, 10), bandwidth = 3,
                           degree = 1))
  expect_no_error(stepDiff(image = matrix(0, 10, 10), bandwidth = 3, degree = 1,
                           blur = TRUE))
})

print("This is the end of test-stepDiff")
