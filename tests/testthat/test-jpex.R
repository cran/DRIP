print("This is the beginning of test-jpex")

test_that("only accept matrix as image input", {
  expect_error(jpex(image = 1:100, bandwidth = as.integer(2), sigma =
                      0.00623, alpha = 0.001), "image data must be a matrix")
})

test_that("only accept proper input", {
  expect_error(jpex(image = matrix(1, 2, 3), bandwidth = 2, sigma = 0.1,
                    alpha = 0.1), "image data must be a square matrix")
  expect_error(jpex(image = matrix(0, 2, 2), bandwidth = "A", sigma = 0.1,
                    alpha = 0.1), "bandwidth must be numeric")
  expect_error(jpex(image = matrix(0, 2, 2), bandwidth = 2:3, sigma = 0.1,
                    alpha = 0.1), "bandwidth must be of length 1")
  expect_error(jpex(image = matrix(0, 2, 2), bandwidth = 0.8, sigma = 0.1,
                    alpha = 0.1), "bandwidth must be a positive integer")
  expect_error(jpex(image = matrix(0, 2, 2), bandwidth = 2, sigma = 0.1,
                    alpha = "A"), "alpha must be numeric")
  expect_error(jpex(image = matrix(0, 2, 2), bandwidth = 2, sigma = 0.1,
                    alpha = -1.2), "alpha must be a number bewteen 0 and 1")
  expect_error(jpex(image = matrix(0, 2, 2), bandwidth = 2, sigma = "A",
                    alpha = 0.1), "sigma must be numeric")
  expect_error(jpex(image = matrix(0, 2, 2), bandwidth = 2, sigma = -1.2,
                    alpha = 0.1), "sigma must be positive")
  img <- matrix(0, 9, 9)
  img[1:4, 1:4] <- 1
  set.seed(100)
  img <- img + matrix(rnorm(9 * 9, sd = 0.1), 9, 9)
  expect_no_error(jpex(image = img, bandwidth = 2, sigma = 0.1, alpha = 0.1))
})

print("This is the end of test-jpex")
