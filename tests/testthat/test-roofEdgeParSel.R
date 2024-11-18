print("This is the beginning of test-roofEdgeParSel")

test_that("only accept proper input", {
  expect_error(roofEdgeParSel(image = 1:4, bandwidth = 3, thresh = 17,
                              nboot = 2, edge1 = matrix(0, 2, 2)),
               "image data must be a matrix")
  expect_error(roofEdgeParSel(image = matrix(0, 3, 4), bandwidth = 3,
                              thresh = 17, nboot = 2, edge1 = matrix(0, 2, 2)),
               "image data must be a square matrix")
  expect_error(roofEdgeParSel(image = matrix(0, 3, 3), bandwidth = "A",
                              thresh = 17, nboot = 2, edge1 = matrix(0, 3, 3)),
               "bandwidth must be numeric")
  expect_error(roofEdgeParSel(image = matrix(0, 3, 3), bandwidth = c(-1, 3),
                              thresh = 17, nboot = 2, edge1 = matrix(0, 3, 3)),
               "every candidate bandwidth must be a positive integer")
  expect_error(roofEdgeParSel(image = matrix(0, 3, 3), bandwidth = 3,
                              thresh = "A", nboot = 2, edge1 = matrix(0, 3, 3)),
               "threshold must be numeric")
  expect_error(roofEdgeParSel(image = matrix(0, 3, 3), bandwidth = 3,
                              thresh = 17, nboot = 2:3,
                              edge1 = matrix(0, 3, 3)),
               "nboot must be a positive integer.")
  expect_error(roofEdgeParSel(image = matrix(0, 3, 3), bandwidth = 3,
                              thresh = 17, nboot = 2,
                              edge1 = matrix(1.2, 3, 3)),
               "edge1's must be either 0 or 1")
  expect_error(roofEdgeParSel(image = matrix(0, 3, 3), bandwidth = 3,
                              thresh = 17, nboot = 2, edge1 = matrix(0, 2, 3)),
               "edge1 must be a square matrix")
  expect_error(roofEdgeParSel(image = matrix(0, 3, 3), bandwidth = 3,
                              thresh = 17, nboot = 2, edge1 = matrix(0, 2, 2)),
               "edge1 and image are not of the same size")
  edge1 <- stepEdge(image = sar, bandwidth = 10, thresh = 17)
  set.seed(100)
  expect_no_error(parSel <- roofEdgeParSel(image = sar,
                                           bandwidth = 10, thresh = 800,
                                           nboot = 2, edge1 = edge1,
                                           blur = TRUE))
})

test_that("return an S3 class Edge_Parameters", {
  edge1 <- stepEdge(image = sar, bandwidth = 10, thresh = 17)
  set.seed(100)
  parSel <- roofEdgeParSel(image = sar, bandwidth = 10,
                           thresh = 800, nboot = 2, edge1 = edge1,
                           blur = FALSE)
  expect_s3_class(parSel, "Edge_Parameters", exact = TRUE)
})

print("This is the end of test-roofEdgeParSel")
