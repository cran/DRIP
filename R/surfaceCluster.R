# This is R source code for function 'surfaceCluster', in the
# R package "DRIP".
# Date: September 05, 2017
# Creator: Yicheng Kang

surfaceCluster <- function(image, bandwidth, sig.level, sigma, phi0,
                           mean_std_abs, cw = 3, blur = FALSE, plot = FALSE){
  if (!is.matrix(image)) {
    stop("image data must be a matrix")
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
  if (!is.numeric(sig.level) || abs(sig.level - 0.5) > 0.5)
  { stop('sig.level must be a number between 0 and 1') }
  n1 <- dim(image)[1]
  z <- matrix(as.double(image), ncol = n1)
  zq <- as.double(qnorm(sig.level))
  if (missing(sigma) || missing(phi0) || missing(mean_std_abs)) {
    jp.llk <- JPLLK_surface(z, 2:7)
    fitted <- jp.llk$fitted
    resid <- jp.llk$resid
    sigma <- as.double(jp.llk$sigma)
    std_resid <- resid / sigma
    phi0 <- as.double(density(x = std_resid, bw = 1.06*n1^(-2/5),
                              kernel = "gaussian", n = 4, from = -1,
                              to = 2)$y[2])
    mean_std_abs <- as.double(mean(abs(std_resid)))
  }
  k <- as.integer(bandwidth)
  cw <- as.integer(cw)
  if (blur == FALSE) {
    out <- .Fortran(C_cluster_cwm_denoise, n = as.integer(n1 - 1),
                    obsImg = z, k = k, zq = zq, sigma = sigma, phi0 = phi0,
                    mean_std_abs = mean_std_abs, cw = cw, estImg = z)
  }
  else {
    out <- .Fortran(C_cluster_cwm_deblur, n = as.integer(n1 - 1), obsImg = z,
                    k = k, zq = zq, sigma = sigma, phi0 = phi0,
                    mean_std_abs = mean_std_abs, cw = cw, estImg = z)
  }
  if (plot == FALSE) {
    return(list(estImg = out$estImg, sigma = sigma, phi0 = phi0,
                mean_std_abs = mean_std_abs))
  }
  else { image(out$estImg, col = gray((0:255)/255))
    return(list(estImg = out$estImg, sigma = sigma, phi0 = phi0,
                mean_std_abs = mean_std_abs))
  }
}
