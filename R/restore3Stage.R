# This is R source code for function 'restore3Stage', in the
# R package "image".
# Date: May 7, 2013
# Creator: Yicheng Kang

restore3Stage <- function(image, bandwidth, step_edge = NULL, roof_edge = NULL, edge1 = NULL, edge2 = NULL,
                       blur = FALSE, plot = FALSE){
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
  if (length(bandwidth) > 1) {
    stop('bandwidth must be a positive integer')
  }
#   if (n1 + 2 * bandwidth + 2 > 600)
#     stop("some choice of bandwidth or the resolution of the
# image is too large")
  n1 <- dim(image)[1]
  z <- matrix(as.double(image), ncol = n1)
  k <- as.integer(bandwidth)
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
  roof_edge <- matrix(as.integer(roof_edge), ncol = n2)
  if (blur == TRUE) {
    out <- .Fortran(C_deblur_3stage, n = as.integer(n1 - 1),
                    obsImg = z, bandwidth = k, edge1 = step_edge, edge2 = roof_edge,
                    estImg = z)
  }
  else {
    out <- .Fortran(C_denoise_3stage, n = as.integer(n1 - 1),
                    obsImg = z, bandwidth = k, edge1 = step_edge, edge2 = roof_edge,
                    estImg = z)
  }

  if (plot == FALSE) {
    return(out$estImg)
  }
  else { x <- seq(0, 1, length = n1); y <- x
  image(x, y, out$estImg, zlim = c(0, 255),
        col = gray((0:255)/255))
  return(out$estImg)
  }
}
