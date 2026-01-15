# This is R source code for function 'roofEdge', in the R
# package "image".
# Date: May 6, 2013
# Creator: Yicheng Kang

roofEdge <- function(image, bandwidth, thresh, step_edge = NULL, edge1 = NULL,
                     blur = FALSE, plot = FALSE){
  if (!is.null(edge1)) {
    step_edge <- edge1
  }
  if (!is.matrix(image))
    stop("image data must be a matrix.")
  n1 <- as.integer(dim(image)[1])
  n2 <- as.integer(dim(image)[2])
  if (n1 != n2)
    stop("image data must be a square matrix.")
  if (!is.numeric(bandwidth) || length(bandwidth) > 1 ||
      as.integer(bandwidth) < 1)
    stop("bandwidth must be a positive integer.")
  if (!is.numeric(thresh) || length(thresh) > 1 || thresh < 0)
    stop("threshold  must be a positive number.")
  if(!is.matrix(step_edge) || ncol(step_edge) != nrow(step_edge))
    stop("step_edge must be a square matrix")
  if(!all(step_edge == 0 | step_edge == 1))
    stop("step_edge's must be either 0 or 1.")
  if(ncol(step_edge) != n1)
    stop("step_edge and image are not of the same size.")
  n1 <- dim(image)[1]
  z <- matrix(as.double(image), ncol = n1)
  step_edge <- matrix(as.integer(step_edge), ncol=n1)
  edge2 <- array(as.integer(0), c(n1, n1))
  k <- as.integer(bandwidth)
  u <- as.double(thresh)
  if (blur == FALSE) {
    out <- .Fortran(C_roofdetect_denoise, n = as.integer(n1 - 1),
                    obsImg = z, bandwidth = as.integer(k), thresh = u,
                    edge1 = step_edge, edge2 = edge2)
  }
  else {
    out <- .Fortran(C_roofdetect_deblur, n = as.integer(n1 - 1),
                    obsImg = z, bandwidth = as.integer(k), thresh = u,
                    edge1 = step_edge, edge2 = edge2)
  }
  edge <- out$edge2
  if (plot == FALSE) { return(edge) }
  if (plot == TRUE) {
    x <- seq(0, 1, length = n1); y <- x
    image(x, y, 1 - edge, col = gray(0:1))
    return(edge)
  }
}
