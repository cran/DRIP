print("This is the beginning of test-roofDiff")

test_that("only accept matrix as image input", {
  expect_error(roofDiff(image = 1:4, bandwidth = 3),
               "image data must be a matrix")
  expect_error(roofDiff(image = matrix(0, 2, 3), bandwidth = 2),
               "image data must be a square matrix")
  expect_error(roofDiff(image = matrix(0, 3, 3), bandwidth = "A"),
               "bandwidth must be numeric")
  expect_error(roofDiff(image = matrix(0, 3, 3), bandwidth = -1),
               "bandwidth must be a positive integer")
  img <- matrix(0, 9, 9)
  img[1:4, 1:4] <- 1
  expect_no_error(roofDiff(image = img, bandwidth = 2, blur = FALSE))
  expect_no_error(roofDiff(image = img, bandwidth = 2, blur = TRUE))
})

print("This is the end of test-roofDiff")
