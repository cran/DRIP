# This is R source code for function 'restore3StageParSel', in the
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

restore3StageParSel <- function(image, bandwidth, step_edge = NULL,
                                roof_edge = NULL, edge1 = NULL, edge2 = NULL,
                                nboot, blur = FALSE){
  if (!is.null(edge1)) {
    step_edge <- edge1
  }
  if (!is.null(edge2)) {
    roof_edge <- edge2
  }
  if (!is.matrix(image) || !is.matrix(step_edge) || !is.matrix(roof_edge)) {
    stop("image, step_edge and roof_edge must be a matrix")
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
  if (nrow(step_edge) != n1 || ncol(step_edge) != n1) {
    stop("different size in step_edge and image")
  }
  if (nrow(roof_edge) != n1 || ncol(roof_edge) != n1) {
    stop("different size in roof_edge and image")
  }
  if(!all(step_edge == 0 | step_edge == 1))
    stop("step_edge must be either 0 or 1.")
  if(!all(roof_edge == 0 | roof_edge == 1))
    stop("roof_edge must be either 0 or 1.")
  step_edge <- matrix(as.integer(step_edge), ncol = n1)
  roof_edge <- matrix(as.integer(roof_edge), ncol = n1)
  n_band <- length(bandwidth)
  if (blur == FALSE) {
    out <- .Fortran(C_denoise_3stage_bandwidth, n = as.integer(n1 - 1),
                    obsImg = z, nband = n_band,
                    bandwidth = as.integer(bandwidth), edge1 = step_edge,
                    edge2 = roof_edge, cv = rep(as.double(0), n_band))
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
                    bandwidth = as.integer(bandwidth), edge1 = step_edge,
                    edge2 = roof_edge, nboot = n_boot,
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

