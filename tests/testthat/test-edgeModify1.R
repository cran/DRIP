print("This is the beginning of test-edgeModify1")

test_that("only accept proper input", {
  expect_error(edgeModify1(bandwidth = 2, image = 1:4, edge = matrix(0, 2, 2)),
               "obsImg must be a square matrix")
  expect_error(edgeModify1(bandwidth = 2, image = matrix(0, 2, 2), edge = 1:4),
               "edge must be a square matrix")
  expect_error(edgeModify1(bandwidth = 2, image = matrix(0, 4, 4),
                       edge = matrix(1.2, 4, 4), plot = FALSE),
               "edge can only have entry equal to 0 or 1")
  expect_error(edgeModify1(bandwidth = 2, image = matrix(0, 3, 3),
                       edge = matrix(0, 2, 2)),
               "obsImg and edge must have the same size")
  expect_error(edgeModify1(bandwidth = -1, image = matrix(0, 2, 2),
                       edge = matrix(0, 2, 2)),
               "bandwidth must be a positive integer")
  expect_no_error(edgeModify1(bandwidth = 2, image = matrix(0, 10, 10),
                          edge = matrix(0, 10, 10), plot = FALSE))
  expect_no_error(edgeModify1(bandwidth = 2, image = matrix(0, 10, 10),
                          edge = matrix(0, 10, 10), plot = TRUE))
  img <- matrix(0, 10, 10)
  img[1:5, 1:5] <- 1
  edge <- matrix(0, 10, 10)
  edge[4:6, 1:5] <- 1
  edge[1:5, 4:6] <- 1
  expect_no_error(edgeModify1(bandwidth = 3, image = img, edge = edge, plot = TRUE))
  img <- matrix(0, 10, 10)
  img[upper.tri(img, diag = TRUE)] <- 1
  edge <- matrix(0, 10, 10)
  for (i in 1:10) {
    for (j in 1:10) {
     if (abs(i - j) < 2) {
       edge[i, j] <- 1
     }
    }
  }
  expect_no_error(edgeModify1(bandwidth = 2, image = img, edge = edge, plot = TRUE))
  edge <- stepEdge(sar, bandwidth = 4, thresh = 20, degree = 0)
  expect_no_error(edgeModify1(4, sar, edge, plot = TRUE))
})

print("This is the end of test-edgeModify1")
