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

test_that("output test -- remove blur", {
  # 1. Create a 100x100 true underlying image (a sharp square in the center)
  true_img <- matrix(0, nrow = 100, ncol = 100)
  true_img[30:70, 30:70] <- 1

  # 2. Define the 3x3 uniform blurring kernel (mean filter)
  kernel <- matrix(1/9, nrow = 3, ncol = 3)

  # 3. 2D Convolution Function
  convolve2d <- function(img, kernel) {
    nr <- nrow(img)
    nc <- ncol(img)
    knr <- nrow(kernel)
    knc <- ncol(kernel)

    pad_r <- floor(knr / 2)
    pad_c <- floor(knc / 2)

    # Pad the image edges with zero (zero-padding)
    padded <- matrix(0, nrow = nr + 2 * pad_r, ncol = nc + 2 * pad_c)
    padded[(pad_r + 1):(pad_r + nr), (pad_c + 1):(pad_c + nc)] <- img

    out <- matrix(0, nrow = nr, ncol = nc)

    # Perform convolution
    for (i in 1:nr) {
      for (j in 1:nc) {
        sub_mat <- padded[i:(i + knr - 1), j:(j + knc - 1)]
        out[i, j] <- sum(sub_mat * kernel)
      }
    }
    return(out)
  }

  # 4. Generate the blurred image
  blurred_img <- convolve2d(true_img, kernel)

  # 5. Plot the original vs. blurred images side by side
  par(mfrow = c(1, 2), mar = c(2, 2, 2, 2))

  image(true_img, col = gray.colors(256), main = "Original True Image", axes = FALSE)
  image(blurred_img, col = gray.colors(256), main = "Blurred Image (3x3 Uniform)", axes = FALSE)

  set.seed(2026)
  obs_img <- blurred_img + matrix(rnorm(100 * 100, sd = 0.1), 100, 100)
  image(obs_img, col = gray.colors(256))
  res <- jpex(image = obs_img, bandwidth = 3, alpha = 0.05, sigma = 0.1)$deblurred
  image(res, col = gray.colors(256))

  expect_lte(mean((res - true_img)^2), mean((obs_img - true_img)^2)) # expect to improve upon the observed image in terms of MSE.

})

print("This is the end of test-jpex")
