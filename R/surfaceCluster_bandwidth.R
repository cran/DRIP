# This is R source code for function 'surfaceCluster_bandwidth', in the
# R package "DRIP".
# Date: September 07, 2017
# Creator: Yicheng Kang

# The print method for S3 class
foo <- list(cv_dataframe = data.frame(x = 1:2, y = 3:4),
            bandwidth = as.integer(3), sigma = 0.0, phi0 = 1.0,
            mean_std_abs = 1.0)
class(foo) <- "Surface_Cluster_Parameters"

print.Surface_Cluster_Parameters <- function(x, type = "all", ...) {
  stopifnot(class(x) == "Surface_Cluster_Parameters")
  if (!(type %in% c("cv_scores", "bandwidth", "sigma", "phi0", "mean_std_abs",
                    "all"))) {
    stop("Wrong type value.")
  } else {
    if (type == "cv_scores") {
      cat("The cross validation scores:\n")
      print.data.frame(x$cv_dataframe)
    }
    if (type == "bandwidth") {
      cat("The selected bandwidth: ", x$bandwidth, "\n")
    }
    if (type == "sigma") {
      cat("The estimated sigma: ", x$sigma, "\n")
    }
    if (type == "phi0") {
      cat("The estimated value of the density at 0: ", x$phi0, "\n")
    }
    if (type == "mean_std_abs") {
      cat("The estimated mean of absolute error: ", x$mean_std_abs, "\n")
    }
    if (type == "all") {
      cat("The cross validation scores:\n")
      print.data.frame(x$cv_dataframe)
      cat("The selected bandwidth: ", x$bandwidth, "\n")
      cat("The estimated sigma: ", x$sigma, "\n")
      cat("The estimated value of the density at 0: ", x$phi0, "\n")
      cat("The estimated mean of absolute error: ", x$mean_std_abs, "\n")
    }
  }
}

# The summary method for the S3 class

summary.Surface_Cluster_Parameters <- function(object, ...) {
  stopifnot(class(object) == "Surface_Cluster_Parameters")
  cat("The selected bandwidth: ", object$bandwidth, "\n")
  cat("The estimated sigma: ", object$sigma, "\n")
  cat("The estimated value of the density at 0: ", object$phi0, "\n")
  cat("The estimated mean of absolute error: ", object$mean_std_abs, "\n")
}

# The plot method for the S3 class

plot.Surface_Cluster_Parameters <- function(x, ...) {
  stopifnot(class(x) == "Surface_Cluster_Parameters")
  plot.default(x = x$cv_dataframe[, 1], y = x$cv_dataframe[, 2], type = "b",
               xlab = "Bandwidth", ylab = "(Modified) CV")
}

surfaceCluster_bandwidth <- function(image, bandwidths, sig.level, sigma, phi0,
                                     mean_std_abs, relwt = 0.5, cw = 3,
                                     blur = FALSE){
  if (!is.matrix(image)) {
    stop("image data must be a matrix")
  } else {
    n1 <- dim(image)[1]
    n2 <- dim(image)[2]
  }
  if (n1 != n2)
    stop("image data must be a square matrix")
  if (!is.numeric(bandwidths))
    stop("bandwidth must be numeric")
  if (any(bandwidths < 1))
    stop("All bandwidths must be positive integers.")
  if (!is.numeric(sig.level) || abs(sig.level - 0.5) > 0.5)
    stop('sig.level must be a number between 0 and 1')
  if (!is.numeric(relwt) || abs(relwt - 0.5) > 0.5)
    stop("relwt must be a number between 0 and 1.")
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
  nband <- length(bandwidths)
  bandwidths <- as.integer(bandwidths)
  cw <- as.integer(cw)
  bandwidth_hat <- as.integer(0)
  cv <- as.double(rep(0, nband))
  cv_cty <- cv
  cv_jump <- cv
  if (blur == FALSE) {
    out <- .Fortran(C_cluster_cwm_denoise_bandwidth, n = as.integer(n1 - 1),
                    obsImg = z, nband = nband, bandwidths = bandwidths, zq = zq,
                    sigma = sigma, phi0 = phi0, mean_std_abs = mean_std_abs,
                    cw = cw, bandwidth_hat = bandwidth_hat, cv = cv)
  }
  else {
    out <- .Fortran(C_cluster_cwm_deblur_bandwidth, n = as.integer(n1 - 1),
                    obsImg = z, nband = nband, bandwidths = bandwidths, zq = zq,
                    sigma = sigma, phi0 = phi0, mean_std_abs = mean_std_abs,
                    cw = cw, relwt = as.double(relwt),
                    bandwidth_hat = bandwidth_hat, cv = cv, cv_jump = cv_jump,
                    cv_cty = cv_cty)
  }

  if (blur == FALSE){
    cv_dataframe <- data.frame(bandwidths = bandwidths, cv = out$cv)
  } else {
    cv_dataframe <- data.frame(bandwidths = bandwidths, mcv = out$cv,
                               cv_jump = out$cv_jump, cv_cty = out$cv_cty)
  }

  out1 <- list(cv_dataframe = cv_dataframe, bandwidth = out$bandwidth_hat,
               sigma = sigma, phi0 = phi0, mean_std_abs = mean_std_abs)
  class(out1) <- "Surface_Cluster_Parameters"

  return(out1)

}
