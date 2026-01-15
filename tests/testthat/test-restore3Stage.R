print("This is the beginning of test-restore3Stage")

test_that("only accept proper input", {
  expect_error(restore3Stage(image = 1:4, bandwidth = 3, edge1 = matrix(0, 2, 2),
                          edge2 = matrix(0, 2, 2)),
               "image, step_edge and roof_edge must be a matrix")
  expect_error(restore3Stage(image = matrix(0, 2, 3), bandwidth = 3,
                          edge1 = matrix(0, 2, 2), edge2 = matrix(0, 2, 2)),
               "image data must be a square matrix")
  expect_error(restore3Stage(image = matrix(0, 2, 2), bandwidth = "A",
                          edge1 = matrix(0, 2, 2), edge2 = matrix(0, 2, 2)),
               "bandwidth must be numeric")
  expect_error(restore3Stage(image = matrix(0, 2, 2), bandwidth = 2:3,
                          edge1 = matrix(0, 2, 2), edge2 = matrix(0, 2, 2)),
               "bandwidth must be a positive integer")
  expect_error(restore3Stage(image = matrix(0, 2, 2), bandwidth = 2,
                          edge1 = matrix(1.2, 2, 2), edge2 = matrix(0, 2, 2)),
               "step_edge must be either 0 or 1")
  expect_error(restore3Stage(image = matrix(0, 2, 2), bandwidth = 2,
                          edge1 = matrix(0, 2, 2), edge2 = matrix(1.2, 2, 2)),
               "roof_edge must be either 0 or 1")
  edge1 <- matrix(0, 10, 10)
  edge1[, 5] <- 1
  edge1[5, ] <- 1
  edge1[1, 1] <- 1
  edge1[2, 2] <- 1
  edge1[3, 3] <- 1
  edge2 <- matrix(0, 10, 10)
  edge2[2, ] <- 1
  edge2[, 2] <- 1
  edge2[7, 7] <- 1
  edge2[8, 8] <- 1
  edge2[9, 9] <- 1
  set.seed(100)
  img <- matrix(rnorm(100), 10, 10)
  expect_no_error(restore3Stage(image = img, bandwidth = 2, edge1 = edge1,
                             edge2 = edge2))
  expect_no_error(restore3Stage(image = img, bandwidth = 2, edge1 = edge1,
                             edge2 = edge2, blur = TRUE))
  expect_no_error(restore3Stage(image = img, bandwidth = 2, edge1 = edge1,
                             edge2 = edge2, blur = TRUE, plot = TRUE))
  expect_no_error(restore3Stage(image = img, bandwidth = 2, edge1 = edge1,
                             edge2 = edge2, blur = FALSE, plot = TRUE))
})

test_that("edges and image are of the same size", {
  expect_error(restore3Stage(image = matrix(1:4, 2, 2), bandwidth = 3,
                          edge1 = matrix(0, 3, 3), edge2 = matrix(0, 2, 2)),
               "different size in step_edge and image")
  expect_error(restore3Stage(image = matrix(1:4, 2, 2), bandwidth = 3,
                          edge1 = matrix(0, 2, 2), edge2 = matrix(0, 3, 3)),
               "different size in roof_edge and image")
})

print("This is the end of test-restore3Stage")
