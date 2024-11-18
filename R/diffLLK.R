# This is R source code for function 'diffLLK' in the R
# package 'image'.
# Creator: Yicheng Kang
# Date: April 29, 2013

diffLLK <- function(image, bandwidth, plot = FALSE){
  if (!is.matrix(image))
    stop("image data must be a matrix")
  n1 <- as.integer(dim(image)[1])
  n2 <- as.integer(dim(image)[2])
  if (n1 != n2)
    stop("image data must be a square matrix")
  if (!is.numeric(bandwidth))
    stop("bandwidth must be numeric")
  if (length(bandwidth) != 1)
    stop("bandwidth must be a positive integer")
  if (bandwidth < 1)
    stop("bandwidth must be a positive integer")
  if (is.logical(plot) == FALSE)
    stop("plot must be TRUE or FALSE")
  z <- matrix(as.double(image), ncol = n1)
  k <- as.integer(bandwidth)
  out <- .Fortran(C_llk_diff, n = as.integer(n1 - 1), obsImg = z,
                  bandwidth = as.integer(k), diff = z)
  if (plot == FALSE){
    return(out$diff)
  } else {
    image(out$diff, col = gray((0:255)/255))
    return(out$diff)
  }
}
