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
    n1 = dim(image)[1]
    n2 = dim(image)[2]
  }
  if (n1 != n2)
    stop("image data must be a square matrix")
  if (!is.numeric(bandwidths))
    stop("bandwidths must be numeric")
  cvs <- mcmapply(onecv.jpex, bandwidth = bandwidths, MoreArgs = list(image = image), mc.cores = ncores)
  n1 <- nrow(image)
  z <- as.double(c(t(image)))
  LLK <- double(n1 * n1)
  band.min <- as.integer(bandwidths[which.min(cvs)])
  out <- .C(C_LOOCV, Zin = z, nin = n1, kin = band.min, LLK = LLK, cv = as.double(0))
  LLK <- out$LLK
  sigma <- sqrt(mean((z - LLK)^2))
  LLK <- matrix(LLK, nrow = n1, byrow = TRUE)
  return(list(LLK = LLK, sigma = sigma, cv = cvs, bandwidths = bandwidths, band.min = band.min))
}
