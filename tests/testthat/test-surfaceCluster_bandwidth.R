print("This is the beginning of test-surfaceCluster_bandwidth")

test_that("only accept S3 class Surface_Cluster_Parameters", {
  foo <- list(cv_dataframe = data.frame(x = 1:2, y = 3:4),
              bandwidth = as.integer(3), sigma = 0.0, phi0 = 1.0,
              mean_std_abs = 1.0)
  expect_error(print.Surface_Cluster_Parameters(x = foo, type = "bandwidth"))
})

test_that("check type value", {
  foo <- list(cv_dataframe = data.frame(x = 1:2, y = 3:4),
              bandwidth = as.integer(3), sigma = 0.0, phi0 = 1.0,
              mean_std_abs = 1.0)
  class(foo) <- "Surface_Cluster_Parameters"
  expect_error(print.Surface_Cluster_Parameters(x = foo, type = "thresh"),
               "Wrong type value")
  expect_no_error(print.Surface_Cluster_Parameters(x = foo,
                                                   type = "cv_scores"))
  expect_no_error(print.Surface_Cluster_Parameters(x = foo, type = "bandwidth"))
  expect_no_error(print.Surface_Cluster_Parameters(x = foo, type = "sigma"))
  expect_no_error(print.Surface_Cluster_Parameters(x = foo, type = "phi0"))
  expect_no_error(print.Surface_Cluster_Parameters(x = foo,
                                                   type = "mean_std_abs"))
  expect_no_error(print.Surface_Cluster_Parameters(x = foo, type = "all"))
})

test_that("only accept S3 class Surface_Cluster_Parameters", {
  foo <- list(cv_dataframe = data.frame(x = 1:2, y = 3:4),
              bandwidth = as.integer(3), sigma = 0.0, phi0 = 1.0,
              mean_std_abs = 1.0)
  expect_error(summary.Surface_Cluster_Parameters(object = foo))
  class(foo) <- "Surface_Cluster_Parameters"
  expect_no_error(summary.Surface_Cluster_Parameters(object = foo))
})

test_that("only accept S3 class Surface_Cluster_Parameters", {
  foo <- list(cv_dataframe = data.frame(x = 1:2, y = 3:4),
              bandwidth = as.integer(3), sigma = 0.0, phi0 = 1.0,
              mean_std_abs = 1.0)
  expect_error(plot.Surface_Cluster_Parameters(x = foo))
  class(foo) <- "Surface_Cluster_Parameters"
  expect_no_error(plot.Surface_Cluster_Parameters(x = foo))
})

test_that("only accept proper input", {
  expect_error(surfaceCluster_bandwidth(image = 1:4, bandwidths = 2,
                                        sig.level = 0.7, sigma = 0.1,
                                        phi0 = 0.2, mean_std_abs = 1.2),
               "image data must be a matrix")
  expect_error(surfaceCluster_bandwidth(image = matrix(0, 2, 3), bandwidths = 2,
                                        sig.level = 0.7, sigma = 0.1,
                                        phi0 = 0.2, mean_std_abs = 1.2),
               "image data must be a square matrix")
  expect_error(surfaceCluster_bandwidth(image = matrix(0, 2, 2),
                                        bandwidths = "A", sig.level = 0.7,
                                        sigma = 0.1, phi0 = 0.2,
                                        mean_std_abs = 1.2),
               "bandwidth must be numeric")
  expect_error(surfaceCluster_bandwidth(image = matrix(0, 2, 2),
                                        bandwidths = -1, sig.level = 0.7,
                                        sigma = 0.1, phi0 = 0.2,
                                        mean_std_abs = 1.2),
               "All bandwidths must be positive integers")
  expect_error(surfaceCluster_bandwidth(image = matrix(0, 2, 2), bandwidths = 2,
                                        sig.level = -0.7, sigma = 0.1,
                                        phi0 = 0.2, mean_std_abs = 1.2),
               "sig.level must be a number between 0 and 1")
  expect_error(surfaceCluster_bandwidth(image = matrix(0, 2, 2), bandwidths = 2,
                                        sig.level = 0.7, sigma = 0.1,
                                        phi0 = 0.2, mean_std_abs = 1.2,
                                        relwt = 1.2),
               "relwt must be a number between 0 and 1")
  set.seed(100)
  img <- matrix(rnorm(10^4), 100, 100)
  expect_no_error(parsel <- surfaceCluster_bandwidth(image = img,
                                                     bandwidths = 3,
                                                     sig.level = 0.7,
                                                     blur = FALSE))
  expect_no_error(parsel <- surfaceCluster_bandwidth(image = img,
                                                     bandwidths = 3,
                                                     sig.level = 0.7,
                                                     blur = TRUE))
})

print("This is the end of test-surfaceCluster_bandwidth")
