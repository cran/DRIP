print("This is the beginning of test-restore3StageParSel")

test_that("only accept S3 class Three_Stage_Parameters", {
  foo <- list(cv_scores = matrix(c(1.2, 1.3), nrow = 1),
              input = as.integer(3:4), bandwidth = as.integer(4))
  expect_error(print.Three_Stage_Parameters(x = foo, type = "all"))
})

test_that("check the type value", {
  foo <- list(cv_scores = matrix(c(1.2, 1.3), nrow = 1),
              input = as.integer(3:4), bandwidth = as.integer(4))
  class(foo) <- "Three_Stage_Parameters"
  expect_error(print.Three_Stage_Parameters(foo, type = "thresh"),
               "Wrong type value")
  expect_no_error(print.Three_Stage_Parameters(foo, type = "cv_scores"))
  expect_no_error(print.Three_Stage_Parameters(foo, type = "bandwidth"))
  expect_no_error(print.Three_Stage_Parameters(foo, type = "all"))
  expect_no_error(print.Three_Stage_Parameters(foo))
})

test_that("only accept S3 class Three_Stage_Parameters", {
  foo <- list(cv_scores = matrix(c(1.2, 1.3), nrow = 1),
              input = as.integer(3:4), bandwidth = as.integer(4))
  expect_error(summary.Three_Stage_Parameters(object = foo))
  class(foo) <- "Three_Stage_Parameters"
  expect_no_error(summary.Three_Stage_Parameters(object = foo))
})

test_that("only accept S3 class Three_Stage_Parameters", {
  foo <- list(cv_scores = matrix(c(1.2, 1.3), nrow = 1),
              input = as.integer(3:4), bandwidth = as.integer(4))
  expect_error(plot.Three_Stage_Parameters(x = foo))
  class(foo) <- "Three_Stage_Parameters"
  expect_no_error(plot.Three_Stage_Parameters(x = foo))
})

test_that("return an S3 class Three_Stage_Parameters", {
  set.seed(100)
  parSel <- restore3StageParSel(image = matrix(rnorm(100), 10, 10),
                             edge1 = matrix(0, 10, 10),
                             edge2 = matrix(0, 10, 10), bandwidth = 3:4,
                             nboot = 1, blur = FALSE)
  expect_s3_class(parSel, "Three_Stage_Parameters", exact = TRUE)
})

test_that("only accept proper input", {
  expect_error(restore3StageParSel(image = 1:4, bandwidth = 3,
                                edge1 = matrix(0, 2, 2),
                                edge2 = matrix(0, 2, 2), nboot = 1),
               "image, step_edge and roof_edge must be a matrix")
  expect_error(restore3StageParSel(image = matrix(0, 2, 3), bandwidth = 3,
                                edge1 = matrix(0, 2, 2),
                                edge2 = matrix(0, 2, 2), nboot = 1),
               "image data must be a square matrix")
  expect_error(restore3StageParSel(image = matrix(0, 2, 2), bandwidth = "A",
                                edge1 = matrix(0, 2, 2),
                                edge2 = matrix(0, 2, 2), nboot = 1),
               "bandwidth must be numeric")
  expect_error(restore3StageParSel(image = matrix(0, 2, 2), bandwidth = 2,
                                edge1 = matrix(1.2, 2, 2),
                                edge2 = matrix(0, 2, 2)),
               "step_edge must be either 0 or 1")
  expect_error(restore3StageParSel(image = matrix(0, 2, 2), bandwidth = 2,
                                edge1 = matrix(0, 2, 2),
                                edge2 = matrix(1.2, 2, 2)),
               "roof_edge must be either 0 or 1")
  expect_error(restore3StageParSel(image = matrix(0, 2, 2), bandwidth = 2,
                                edge1 = matrix(0, 2, 2),
                                edge2 = matrix(0, 2, 2), blur = TRUE),
               "nboot must be specified when blur is TRUE")
  expect_error(restore3StageParSel(image = matrix(0, 2, 2), bandwidth = 2,
                                edge1 = matrix(0, 2, 2),
                                edge2 = matrix(0, 2, 2), nboot = 1:2,
                                blur = TRUE),
               "nboot must be an integer number")
  set.seed(100)
  img <- matrix(rnorm(100), 10, 10)
  expect_no_error(parsel <- restore3StageParSel(image = img, bandwidth = 2:3,
                                             edge1 = matrix(0, 10, 10),
                                             edge2 = matrix(0, 10, 10),
                                             blur = FALSE))
  edge1 <- stepEdge(image = sar, bandwidth = 10, thresh = 17, degree = 1,
                    blur = FALSE, plot = TRUE)
  edge2 <- roofEdge(image = sar, bandwidth = 10, thresh = 800, edge1 = edge1,
                    blur = FALSE, plot = TRUE)
  set.seed(100)
  expect_no_error(parsel <- restore3StageParSel(image = sar, bandwidth = c(5, 10),
                                             edge1 = edge1, nboot = 2,
                                             edge2 = edge2, blur = TRUE))
  expect_no_error(parsel <- restore3StageParSel(image = sar, bandwidth = c(5, 10),
                                             edge1 = edge1, nboot = 1,
                                             edge2 = edge2, blur = FALSE))
})

test_that("edges and image are of the same size", {
  expect_error(restore3StageParSel(image = matrix(1:4, 2, 2), bandwidth = 3,
                                edge1 = matrix(0, 3, 3),
                                edge2 = matrix(0, 2, 2), nboot = 1),
               "different size in step_edge and image")
  expect_error(restore3StageParSel(image = matrix(1:4, 2, 2), bandwidth = 3,
                                edge1 = matrix(0, 2, 2),
                                edge2 = matrix(0, 3, 3), nboot = 1),
               "different size in roof_edge and image")
})

print("This is the end of test-restore3StageParSel")
