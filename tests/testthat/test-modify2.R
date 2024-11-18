print("This is the beginning of test-modify2")

test_that("only accept proper input", {
  expect_error(modify2(bandwidth = 2, edge = 1:4),
               "edge must be a square matrix")
  expect_error(modify2(bandwidth = 2, edge = matrix(1.2, 4, 4), plot = FALSE),
               "edge can only have entry equal to 0 or 1")
  expect_error(modify2(bandwidth = -1, edge = matrix(1, 4, 4)),
               "bandwidth must be a positive integer")
  stepedge <- stepEdge(image = sar, bandwidth = 10, thresh = 17, degree = 1)
  expect_no_error(modify2(bandwidth = 10, edge = stepedge, plot = FALSE))
  expect_no_error(modify2(bandwidth = 10, edge = stepedge, plot = TRUE))
})

print("This is the end of test-modify2")
