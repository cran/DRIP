# The print method for S3 class
foo <- list(LLK = matrix(0.0, nrow = 2, ncol = 2), sigma = 0.0,
            cv_scores = c(1.2, 1.1), input = as.integer(3:4),
            bandwidth = as.integer(3))
class(foo) <- "JPEX_Parameters"

print.JPEX_Parameters <- function(x, type = "all", ...) {
  stopifnot(class(x) == "JPEX_Parameters")
  if (!(type %in% c("cv_scores", "bandwidth", "sigma", "all"))) {
    stop("Wrong type value.")
  } else {
    if (type == "cv_scores") {
      cat("The cross validation scores:\n")
      print.table(x$cv_scores)
    }
    if (type == "bandwidth") {
      cat("The selected bandwidth: ", x$bandwidth, "\n")
    }
    if (type == "sigma") {
      cat("The estimated sigma: ", x$sigma, "\n")
    }
    if (type == "all") {
      cat("The cross validation scores:\n")
      print.table(x$cv_scores)
      cat("The selected bandwidth: ", x$bandwidth, "\n")
      cat("The estimated sigma: ", x$sigma, "\n")
    }
  }
}

# The summary method for S3 class

summary.JPEX_Parameters <- function(object, ...) {
  stopifnot(class(object) == "JPEX_Parameters")
  cat("The selected bandwidth: ", object$bandwidth, "\n")
  cat("The estimated sigma: ", object$sigma, "\n")
}

# The plot method for S3 class

plot.JPEX_Parameters <- function(x, ...) {
  stopifnot(class(x) == "JPEX_Parameters")
  plot.default(x = x$input, y = x$cv_scores, type = "b", xlab = "Bandwidth",
               ylab = "CV Scores")
}

onecv.jpex <- function(image, bandwidth) {
  n1 <- nrow(image)
  z <- as.double(c(t(image)))
  LLK <- double(n1 * n1)
  k <- as.integer(bandwidth)
  out <- .C(C_LOOCV, Zin = z, nin = n1, kin = k, LLK = LLK, cv = as.double(0))
  return(out$cv)
}

cv.jpex <- function(image, bandwidths, ncpus = 1) {
  ncores <- detectCores()
  ncores <- min(c(ncores, ncpus))

  if (!is.matrix(image)) {
    stop("image data must be a matrix")
  }
  else {
    n1 <- dim(image)[1]
    n2 <- dim(image)[2]
  }
  if (n1 != n2)
    stop("image data must be a square matrix")
  if (!is.numeric(bandwidths))
    stop("bandwidths must be numeric")
  cvs <- mcmapply(onecv.jpex, bandwidth = bandwidths,
                  MoreArgs = list(image = image), mc.cores = ncores)
  n1 <- nrow(image)
  z <- as.double(c(t(image)))
  LLK <- double(n1 * n1)
  band.min <- as.integer(bandwidths[which.min(cvs)])
  out <- .C(C_LOOCV, Zin = z, nin = n1, kin = band.min, LLK = LLK,
            cv = as.double(0))
  LLK <- out$LLK
  sigma <- sqrt(mean((z - LLK)^2))
  LLK <- matrix(LLK, nrow = n1, byrow = TRUE)
  out1 <- list(LLK = LLK, sigma = sigma, cv_scores = cvs, input = bandwidths,
               bandwidth = band.min)
  class(out1) <- "JPEX_Parameters"
  return(out1)
}

