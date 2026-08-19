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

test_that("output test -- noise removal as expected", {
  true_img <- matrix(0, 100, 100)
  true_img[50:100, ] <- 1
  image(true_img, col = gray.colors(256))
  set.seed(2026)
  obs_img <- true_img + matrix(rnorm(100 * 100, sd = 0.1), 100, 100)
  image(obs_img, col = gray.colors(256))
  res <- surfaceCluster(image = obs_img, bandwidth = 3, sig.level = 0.95, sigma = 0.1, plot = TRUE)
  expect_lte(mean((res$estImg - true_img)^2), mean((obs_img - true_img)^2)) # reduces the noise in observed image
})

print("This is the end of test-surfaceCluster")
