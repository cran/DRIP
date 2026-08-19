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
                        edge1 = 1:4), "step_edge must be a square matrix")
  expect_error(roofEdge(image = matrix(0, 3, 3), bandwidth = 2, thresh = 17,
                        edge1 = matrix(1.2, 3, 3)),
               "step_edge's must be either 0 or 1")
  expect_error(roofEdge(image = matrix(0, 3, 3), bandwidth = 2, thresh = 17,
                        edge1 = matrix(0, 2, 2)),
               "step_edge and image are not of the same size")
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

test_that("output test -- detect roof edges as expected", {
  # 1. Define the grid on [0, 1] x [0, 1]
  n <- 100
  x <- seq(0, 1, length.out = n)
  y <- seq(0, 1, length.out = n)

  # 2. Define 1D intensity profile across x
  profile_x <- function(x) {
    if (x < 0.3) {
      return(0.2)
    } else if (x < 0.5) {
      return(0.2 + 1.0 * (x - 0.3))  # slope changes from 0 to 1.0 (roof edge at x = 0.3)
    } else {
      if (x < 0.7) {
        return(0.4 - 1.0 * (x - 0.5)) # slope changes from 1.0 to -1.0 (roof edge at x = 0.5)
      } else {
        return(0.2)  # slope changes from -1.0 to 0.0 (roof edge at x = 0.7)
      }
    }
  }

  # Vectorize the profile function
  profile_vec <- Vectorize(profile_x)

  # 3. Create the 100x100 matrix (extruding profile along the y-axis)
  intensity_1d <- profile_vec(x)
  img_matrix <- matrix(rep(intensity_1d, times = n), nrow = n, ncol = n, byrow = FALSE)

  # 4. Plot both the 2D Image and the 1D Cross-Section
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

  # Plot 1: 2D Grayscale Image
  image(x, y, img_matrix, col = gray.colors(256),
        main = "2D Image with Edges", xlab = "x", ylab = "y", asp = 1)

  # Plot 2: 1D Intensity Cross-Section along x
  plot(x, intensity_1d, type = "l", lwd = 2, col = "black",
       main = "1D Cross-Section Profile", xlab = "x", ylab = "Intensity f(x)",
       ylim = c(0, 0.6))

  # Detect roof edges
  edge1 <- matrix(0, 100, 100)
  true_roof_edge <- edge1
  true_roof_edge[c(30, 50, 70), ] <- 1
  image(true_roof_edge, col = gray.colors(256))
  res <- roofEdge(image = img_matrix, bandwidth = 5, thresh = 0.5, step_edge = edge1, plot = TRUE) # the min jump in slopes is 1.0

  expect_lte(dKQ(res, true_roof_edge), 5/nrow(img_matrix)) # Hausdorff distance bounded by the bandwidth
})

print("This is the end of test-roofEdge")
