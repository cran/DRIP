# This is R source code for function 'roofDiff' in the R
# package 'image'.
# Creator: Yicheng Kang
# Date: May 6, 2013

roofDiff <- function(image, bandwidth, blur = FALSE){
  if (!is.matrix(image))
    stop("image data must be a matrix")
  n1 <- as.integer(dim(image)[1])
  n2 <- as.integer(dim(image)[2])
  if (n1 != n2)
    stop("image data must be a square matrix")
  if (!is.numeric(bandwidth))
    stop("bandwidth must be numeric")
  if (length(bandwidth) != 1 || as.integer(bandwidth) < 1)
    stop("bandwidth must be a positive integer")
  n1 <- dim(image)[1]
  z <- matrix(as.double(image), ncol = n1)
  k <- as.integer(bandwidth)
  if (blur == FALSE) {
    out <- .Fortran(C_roofdiff_denoise, n = as.integer(n1 - 1),
                    obsImg = z, bandwidth = as.integer(k), diff = z)
  }
  else {
    out <- .Fortran(C_roofdiff_deblur, n = as.integer(n1 - 1),
                    obsImg = z, bandwidth = as.integer(k), diff = z)
  }
  return(out$diff)
}
