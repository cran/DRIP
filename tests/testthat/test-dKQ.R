print("This is the beginning of test-dKQ")

test_that("only accept matrix of 0 or 1", {
  expect_error(dKQ(matrix(0, nrow = 2, ncol = 2), matrix(0.5, nrow = 2,
                                                         ncol = 2)),
               "can only have entries of 0 or 1")
})

test_that("only accept proper input", {
  expect_error(dKQ(rep(1, 4), matrix(0, 2, 2)),
               "Both arguments must be matrices")
  expect_error(dKQ(matrix(0, 2, 2), matrix(1, 3, 3)),
               "Two matrices must be of the same size.")
  expect_error(dKQ(matrix(0, 2, 3), matrix(1, 2, 3)),
               "Argument must be a square matrix")
})

test_that("return a numeric", {
  expect_true(is.double(dKQ(matrix(c(1, 0, 0, 1), 2, 2),
                            matrix(c(1, 1, 0, 0), 2, 2))))
})

print("This is the end of test-dKQ")
