print("This is the beginning of test-edgeParSelPilot")

test_that("only accept proper input", {
  expect_error(edgeParSelPilot(sar, edgeType = "step.edges", degree = 0, blur = FALSE, bandwidth = c(3, 5, 7), probs = c(0.8, 0.9, 0.95)))
  expect_error(edgeParSelPilot(sar, edgeType = c("step", "roof"), degree = 0, blur = FALSE, bandwidth = c(3, 5, 7), probs = c(0.8, 0.9, 0.95)))
  expect_error(edgeParSelPilot(sar, edgeType = 1, degree = 0, blur = FALSE, bandwidth = c(3, 5, 7), probs = c(0.8, 0.9, 0.95)))
})

test_that("the step edge pilot works", {
  expect_equal(nrow(edgeParSelPilot(sar, edgeType = "step", degree = 0, blur = FALSE, bandwidth = c(3, 5, 7), probs = c(0.8, 0.9, 0.95))), 3)
  expect_equal(ncol(edgeParSelPilot(sar, edgeType = "step", degree = 0, blur = FALSE, bandwidth = c(3, 5, 7), probs = c(0.8, 0.9))), 2)
  expect_true(is.matrix(edgeParSelPilot(sar, edgeType = "step", degree = 0, blur = FALSE, bandwidth = c(3, 5, 7), probs = c(0.8, 0.9))))
})

test_that("the roof edge pilot works", {
  expect_equal(nrow(edgeParSelPilot(peppers, edgeType = "roof", blur = FALSE, bandwidth = c(5, 7), probs = c(0.8, 0.9, 0.95))), 2)
  expect_equal(ncol(edgeParSelPilot(peppers, edgeType = "roof", blur = FALSE, bandwidth = c(5, 7), probs = c(0.8, 0.9))), 2)
  expect_true(is.matrix(edgeParSelPilot(peppers, edgeType = "roof", blur = FALSE, bandwidth = c(5, 7), probs = c(0.8, 0.9))))
})


print("This is the end of test-edgeParSelPilot")
