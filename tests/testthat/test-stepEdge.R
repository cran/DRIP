print("This is the beginning of test-stepEdge")

test_that("degree either 0 or 1", {
  expect_error(stepEdge(image = sar, bandwidth = 4, thresh = 17, degree = 2))
  expect_no_error(stepEdge(image = matrix(0, 10, 10), bandwidth = 4,
                           thresh = 17, degree = 0))
  expect_no_error(stepEdge(image = matrix(0, 10, 10), bandwidth = 4,
                           thresh = 17, degree = 0, blur = TRUE))
  expect_no_error(stepEdge(image = matrix(0, 10, 10), bandwidth = 4,
                           thresh = 17, degree = 1))
  expect_no_error(stepEdge(image = matrix(0, 10, 10), bandwidth = 4,
                           thresh = 17, degree = 1, blur = TRUE))
})

test_that("output test -- detect step edges as expected", {
  testImg <- matrix(0, 100, 100)
  testImg[50:100, ] <- 1
  edgehat <- stepEdge(image = testImg, bandwidth = 2, thresh = 0.3, degree = 0, plot = FALSE)
  edge <- matrix(0, 100, 100)
  edge[50, ] <- 1
  res <- dKQ(edge, edgehat)
  expect_lte(abs(res), 2.0/nrow(testImg)) # The Hausdorff distance is bounded by bandwidth
  edgehat1 <- stepEdge(image = testImg, bandwidth = 3, thresh = 0.3, degree = 1, plot = TRUE)
  res1 <- dKQ(edge, edgehat1)
  expect_lte(abs(res1), 3.0/nrow(testImg)) # The Hausdorff distance is bounded by bandwidth
})

print("This is the end of test-stepEdge")
