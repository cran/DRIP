# This is R source code for function 'JPLLK_surface', in the
# R package "image".
# Date: April 25, 2013
# Creator: Yicheng Kang

# S3 Class
foo <- list(fitted = matrix(0.0, nrow = 2, ncol = 0), bandwidth = 3.0, resid
             = matrix(0.0, nrow = 2, ncol = 0), sigma = 0.0)
class(foo) <- "JPLLK_Parameters"

# The print method for S3 class

print.JPLLK_Parameters <- function(x, type = "all", ...) {
  stopifnot(class(x) == "JPLLK_Parameters")
  if (type == "bandwidth") {
    cat("The selected bandwidth: ", x$bandwidth, "\n")
  } else {
    if (type == "sigma") {
      cat("The estimated noise level: ", x$sigma, "\n")
    } else {
      if (type == "all") {
        cat("The selected bandwidth: ", x$bandwidth, "\n")
        cat("The estimated noise level: ", x$sigma, "\n")
      } else {
        stop("Wrong values for type.")
      }
    }
  }
}

# The summary method for S3 class

summary.JPLLK_Parameters <- function(object, ...) {
  stopifnot(class(object) == "JPLLK_Parameters")
  cat("The selected bandwidth: ", object$bandwidth, "\n")
  cat("The estimated noise level: ", object$sigma, "\n")
}

# The plot method for S3 class

plot.JPLLK_Parameters <- function(x, ...) {
  stopifnot(class(x) == "JPLLK_Parameters")
  image(x$resid, col = gray(c(0:255)/255), xlab = "", ylab = "")
}

JPLLK_surface <- function(image, bandwidth, plot = FALSE){
  if (!is.matrix(image)) {
    stop("image data must be a matrix")
  } else {
    n1 <- dim(image)[1]
    n2 <- dim(image)[2]
    if (n1 != n2)
      stop("image data must be a square matrix")
    if (!is.numeric(bandwidth))
      stop("bandwidth must be numeric")
#     if (n1 + 2 * max(bandwidth) + 2 > 600)
#       stop("some choice of bandwidth or the resolution of the
# image is too large")
    n1 <- dim(image)[1]
    z <- matrix(as.double(image), ncol = n1)
    n_band <- length(bandwidth)
    out <- .Fortran(C_jp_llk_cv, n = as.integer(n1 - 1), obsImg = z,
                    nband = n_band, bandwidth = as.integer(bandwidth),
                    cv = rep(as.double(0), n_band))
    k.cv <- out$cv
    #cv.band <- mean(bandwidth[k.cv ==  min(k.cv)])
    cv.band <- bandwidth[which.min(k.cv)]
    jp.llk <- .Fortran(C_jp_llk_fit, n = as.integer(n1 - 1),
                      obsImg = z, bandwidth = as.integer(cv.band), fitted = z,
                      resid = z, sigma = as.double(0))
    out1 <- list(fitted = jp.llk$fitted, bandwidth = cv.band, resid
                 = jp.llk$resid, sigma = jp.llk$sigma)
    class(out1) <- "JPLLK_Parameters"
    if (plot == FALSE) {
      return(out1)
    }
    else { image(jp.llk$fitted, col = gray((0:255)/255))
      return(out1)
    }
  }
}

