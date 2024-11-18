print("This is the beginning of test-roofEdge")

test_that("only accept proper input", {
  expect_error(roofEdge(image = matrix(0, 3, 4), bandwidth = 3, thresh = 17,
                        edge1 = matrix(0, 2, 2)),
               "image data must be a square matrix")
  expect_error(roofEdge(image = 1:4, bandwidth = 3, thresh = 17,
                        edge1 = matrix(0, 2, 2)), "image data must be a matrix")
  expect_error(roofEdge(image = matrix(0, 3, 3), bandwidth = -1, thresh = 17,
                        edge1 = matrix(0, 2, 2)),
               "bandwidth must be a positive integer.")
  expect_error(roofEdge(image = matrix(0, 3, 3), bandwidth = 2, thresh = -1,
                        edge1 = matrix(0, 2, 2)),
               "threshold  must be a positive number")
  expect_error(roofEdge(image = matrix(0, 3, 3), bandwidth = 2, thresh = 17,
                        edge1 = 1:4), "edge1 must be a square matrix")
  expect_error(roofEdge(image = matrix(0, 3, 3), bandwidth = 2, thresh = 17,
                        edge1 = matrix(1.2, 3, 3)),
               "edge1's must be either 0 or 1")
  expect_error(roofEdge(image = matrix(0, 3, 3), bandwidth = 2, thresh = 17,
                        edge1 = matrix(0, 2, 2)),
               "edge1 and image are not of the same size")
  img <- matrix(0, 9, 9)
  img[1:4, 1:4] <- 1
  expect_no_error(roofEdge(image = img, bandwidth = 2, thresh = 0.1,
                           edge1 = matrix(1, 9, 9)))
  expect_no_error(roofEdge(image = img, bandwidth = 2, thresh = 0.1,
                           edge1 = matrix(1, 9, 9), plot = TRUE))
  expect_no_error(roofEdge(image = img, bandwidth = 2, thresh = 0.1,
                           edge1 = matrix(1, 9, 9), blur = TRUE))
  expect_no_error(roofEdge(image = img, bandwidth = 2, thresh = 0.1,
                           edge1 = matrix(0, 9, 9)))
  expect_no_error(roofEdge(image = img, bandwidth = 2, thresh = 0.1,
                           edge1 = matrix(0, 9, 9), plot = TRUE))
  expect_no_error(roofEdge(image = img, bandwidth = 2, thresh = 0.1,
                           edge1 = matrix(0, 9, 9), blur = TRUE))
})

print("This is the end of test-roofEdge")
