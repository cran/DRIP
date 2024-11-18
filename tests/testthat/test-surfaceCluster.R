print("This is the beginning of test-surfaceCluster")

test_that("only accept proper input", {
  expect_error(surfaceCluster(image = 1:4, bandwidth = 4, sig.level = .9995,
                              cw = 3, blur = FALSE),
               "image data must be a matrix")
  expect_error(surfaceCluster(image = matrix(0, 2, 3), bandwidth = 2,
                              sig.level = 0.7, sigma = 0.1, phi0 = 0.2,
                              mean_std_abs = 1.2),
               "image data must be a square matrix")
  expect_error(surfaceCluster(image = matrix(0, 2, 2), bandwidth = "A",
                              sig.level = 0.7, sigma = 0.1, phi0 = 0.2,
                              mean_std_abs = 1.2),
               "bandwidth must be numeric")
  expect_error(surfaceCluster(image = matrix(0, 2, 2), bandwidth = 2:3,
                              sig.level = 0.7, sigma = 0.1, phi0 = 0.2,
                              mean_std_abs = 1.2),
               "bandwidth must be a positive integer")
  expect_error(surfaceCluster(image = matrix(0, 2, 2), bandwidth = 2,
                              sig.level = -0.7, sigma = 0.1, phi0 = 0.2,
                              mean_std_abs = 1.2),
               "sig.level must be a number between 0 and 1")
  img <- matrix(0, 100, 100)
  img[45:55, 50:70] <- 3
  set.seed(100)
  img <- img + matrix(rnorm(10^4), 100, 100)
  expect_no_error(surfaceCluster(image = img, bandwidth = 3, sig.level = 0.7,
                                 blur = FALSE))
  expect_no_error(surfaceCluster(image = img, bandwidth = 3, sig.level = 0.7,
                                 blur = TRUE))
  expect_no_error(surfaceCluster(image = img, bandwidth = 3, sig.level = 0.7,
                                 blur = FALSE, plot = TRUE))
  expect_no_error(surfaceCluster(image = img, bandwidth = 3, sig.level = 0.7,
                                 blur = TRUE, plot = TRUE))
})

print("This is the end of test-surfaceCluster")
