# This is R source code for function 'threeStageParSel', in the
# R package "image".
# Date: Aug 30, 2015
# Creator: Yicheng Kang

# S3 Class
foo <- list(cv_scores = matrix(c(1.2, 1.3), nrow = 1), input = as.integer(3:4),
            bandwidth = as.integer(4))
class(foo) <- "Three_Stage_Parameters"

# The print method for the S3 class

print.Three_Stage_Parameters <- function(x, type = "all", ...) {
  stopifnot(class(x) == "Three_Stage_Parameters")
  if (!(type %in% c("cv_scores", "bandwidth", "all"))) {
    stop("Wrong type value")
  } else {
    if (type == "cv_scores") {
      cat("The cross validation scores:\n")
      print.table(x$cv_scores)
    }
    if (type == "bandwidth") {
      cat("The selected bandwidth: ", x$bandwidth, "\n")
    }
    if (type == "all") {
      cat("The cross validation scores:\n")
      print.table(x$cv_scores)
      cat("The selected bandwidth: ", x$bandwidth, "\n")
    }
  }
}

# The summary method for the S3 class

summary.Three_Stage_Parameters <- function(object, ...) {
  stopifnot(class(object) == "Three_Stage_Parameters")
  cat("The cross validation scores:\n")
  print.table(object$cv_scores)
  cat("The selected bandwidth: ", object$bandwidth, "\n")
}

# The plot method for the S3 class

plot.Three_Stage_Parameters <- function(x, ...) {
  stopifnot(class(x) == "Three_Stage_Parameters")
  plot.default(x = x$input, y = x$cv_scores, type = "b", xlab = "Bandwidth",
               ylab = "(MSE)CV Scores")
}

threeStageParSel <- function(image, bandwidth, edge1, edge2, nboot,
                             blur = FALSE){
  if (!is.matrix(image) || !is.matrix(edge1) || !is.matrix(edge2)) {
    stop("image, edge1 and edge2 must be a matrix")
  } else {
    n1 <- dim(image)[1]
    n2 <- dim(image)[2]
  }
  if (n1 != n2)
    stop("image data must be a square matrix")
  if (!is.numeric(bandwidth))
    stop("bandwidth must be numeric")
#   if (n1 + 2 * max(bandwidth) + 2 > 600)
#     stop("some choice of bandwidth or the resolution of the
# image is too large")
  n1 <- dim(image)[1]
  z <- matrix(as.double(image), ncol = n1)
  if (nrow(edge1) != n1 || ncol(edge1) != n1) {
    stop("different size in edge1 and image")
  }
  if (nrow(edge2) != n1 || ncol(edge2) != n1) {
    stop("different size in edge2 and image")
  }
  if(!all(edge1 == 0 | edge1 == 1))
    stop("edge1 must be either 0 or 1.")
  if(!all(edge2 == 0 | edge2 == 1))
    stop("edge2 must be either 0 or 1.")
  edge1 <- matrix(as.integer(edge1), ncol = n1)
  edge2 <- matrix(as.integer(edge2), ncol = n1)
  n_band <- length(bandwidth)
  if (blur == FALSE) {
    out <- .Fortran(C_denoise_3stage_bandwidth, n = as.integer(n1 - 1),
                    obsImg = z, nband = n_band,
                    bandwidth = as.integer(bandwidth), edge1 = edge1,
                    edge2 = edge2, cv = rep(as.double(0), n_band))
    k.cv <- out$cv
    band_sel <- mean(bandwidth[k.cv ==  min(k.cv)])
    out.mat <- matrix(0, ncol = n_band, nrow = 1)
    out.mat[1, ] <- k.cv
    colnames(out.mat) <- paste('bandwidth=', bandwidth, sep='')
    rownames(out.mat) <- "CV-score"

  }
  else {
    if(missing(nboot))
      stop("nboot must be specified when blur is TRUE.")
    if (length(nboot) > 1)
      stop("nboot must be an integer number.")
    n_boot <- as.integer(nboot)
    out <- .Fortran(C_deblur_3stage_bandwidth, n = as.integer(n1 - 1),
                    obsImg = z, nband = n_band,
                    bandwidth = as.integer(bandwidth), edge1 = edge1,
                    edge2 = edge2, nboot = n_boot,
                    msecv = rep(as.double(0), n_band))
    k.msecv <- out$msecv
    band_sel <- mean(bandwidth[k.msecv ==  min(k.msecv)])
    out.mat <- matrix(0, ncol = n_band, nrow = 1)
    out.mat[1, ] <- k.msecv
    colnames(out.mat) <- paste('bandwidth=', bandwidth, sep='')
    rownames(out.mat) <- "MSECV-score"
  }

  # message(paste('The selected bandwidth is', band_sel))
  out1 <- list(cv_scores = out.mat, input = bandwidth, bandwidth = band_sel)
  class(out1) <- "Three_Stage_Parameters"
  return(out1)
}

