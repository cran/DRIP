jpex <- function(image, bandwidth, alpha, sigma) {
    if (!is.matrix(image)) {
        stop("image data must be a matrix")
    }
    n1 <- dim(image)[1]
    n2 <- dim(image)[2]
    if (n1 != n2)
        stop("image data must be a square matrix")
    if (!is.numeric(bandwidth))
        stop("bandwidth must be numeric")
    if (length(bandwidth) != 1)
        stop("bandwidth must be of length 1")
    if (bandwidth < 1)
        stop("bandwidth must be a positive integer")
    if (!is.numeric(alpha))
        stop("alpha must be numeric")
    if (alpha <= 0 || alpha >= 1)
        stop("alpha must be a number bewteen 0 and 1")
    if (!is.numeric(sigma))
        stop("sigma must be numeric")
    if (sigma <= 0)
        stop("sigma must be positive")
    z <- as.double(t(image))
    edge <- double(n1 * n1)
    fhat <- double(n1 * n1)
    out <- .C(C_JPEX0, Z = z, nin = n1, kin = as.integer(bandwidth),
              alphain = as.double(alpha), sigmain = as.double(sigma),
              EDGE = edge, fhat = fhat)
    fhat <- matrix(out$fhat, nrow = n1, byrow = TRUE)
    edge <- matrix(out$EDGE, nrow = n1, byrow = TRUE)
    return(list(deblurred = fhat, edge = edge))
}
